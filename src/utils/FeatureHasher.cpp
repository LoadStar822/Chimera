#include "utils/FeatureHasher.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <span>
#include <stdexcept>
#include <string>

#include <xxhash.h>

#ifdef CHIMERA_HAS_STROBEMERS
#include <strobemers/indexparameters.hpp>
#include <strobemers/randstrobes.hpp>
#endif

#include "utils/Syncmer.hpp"

namespace chimera::feature
{

namespace
{

#ifdef CHIMERA_HAS_STROBEMERS
constexpr uint64_t kDefaultRandstrobeQ = (1ULL << 8) - 1; // matches c=8 in strobealign
constexpr unsigned kDefaultAuxLen = 17;

inline uint64_t mix64(const uint64_t * data, size_t count, uint64_t seed)
{
    return XXH3_64bits_withSeed(data, static_cast<size_t>(count) * sizeof(uint64_t), seed);
}

inline uint64_t syncmer_packed_hash(uint64_t input) noexcept
{
    constexpr uint64_t prime1 = 0x9E3779B185EBCA87ULL;
    constexpr uint64_t prime2 = 0xC2B2AE3D27D4EB4FULL;
    constexpr uint64_t prime3 = 0x165667B19E3779F9ULL;
    constexpr uint64_t prime4 = 0x85EBCA77C2B2AE63ULL;
    constexpr uint64_t prime5 = 0x27D4EB2F165667C5ULL;

    uint64_t result = prime5 + 8;
    input *= prime2;
    input = std::rotl(input, 31);
    result ^= input * prime1;
    result = std::rotl(result, 27);
    result = result * prime1 + prime4;
    result ^= result >> 33;
    result *= prime2;
    result ^= result >> 29;
    result *= prime3;
    result ^= result >> 32;
    return result;
}

inline int derive_syncmer_s(int k)
{
    // 目标：k - s ≈ 16 且为偶数，避免过密
    const int span = 16;
    int s = (k > span) ? (k - span) : (k - 2);
    if (((k - s) & 1) != 0)
    {
        --s;
    }
    if (s < 4)
    {
        s = 4;
    }
    if (s > k - 2)
    {
        s = k - 2;
    }
    if (((k - s) & 1) != 0)
    {
        if (s > 4)
            --s;
        else if (s < k - 2)
            ++s;
    }
    return std::clamp(s, 2, k);
}

inline int derive_max_dist(int k, unsigned w_max)
{
    // Allow the second strobe to extend roughly a read length (k + w_max) but clamp to hardware limit.
    int suggested = static_cast<int>(k) + static_cast<int>(w_max);
    suggested = std::max(suggested, k);
    return std::clamp(suggested, k, 255);
}

inline uint64_t make_strobemer_hash(const Randstrobe & strobe,
                                    const StrobemerParams & params)
{
    const uint64_t content = params.canonical
                                 ? std::min(strobe.hash, strobe.hash_revcomp)
                                 : strobe.hash;
    return mix64(&content, 1, params.seed);
}

void collect_syncmers_encoded(std::span<const uint8_t> seq,
                              SyncmerParameters parameters,
                              std::vector<Syncmer> & out)
{
    out.clear();
    constexpr size_t kQueueCapacity = 64;
    constexpr size_t kQueueMask = kQueueCapacity - 1;
    const uint8_t * seq_ptr = seq.data();
    const size_t seq_size = seq.size();

    const size_t window_size = static_cast<size_t>(parameters.k - parameters.s + 1);
    const size_t target_index = static_cast<size_t>(parameters.t_syncmer - 1);
    if (window_size >= kQueueCapacity)
    {
        throw std::runtime_error("syncmer window exceeds encoded collector capacity");
    }

    const uint64_t kmask = (1ULL << (2 * parameters.k)) - 1;
    const uint64_t smask = (1ULL << (2 * parameters.s)) - 1;
    const uint64_t kshift = static_cast<uint64_t>(parameters.k - 1) * 2;
    const uint64_t sshift = static_cast<uint64_t>(parameters.s - 1) * 2;

    uint64_t qs[kQueueCapacity] = {};
    uint64_t qs_min_candidates[kQueueCapacity] = {};
    size_t qs_head = 0;
    size_t qs_size = 0;
    size_t qs_min_head = 0;
    size_t qs_min_size = 0;
    size_t l = 0;
    uint64_t xk[2] = {0, 0};
    uint64_t xs[2] = {0, 0};

    for (size_t i = 0; i < seq_size; ++i)
    {
        const uint8_t base = seq_ptr[i];
        if (base >= 4)
        {
            qs_head = 0;
            qs_size = 0;
            qs_min_head = 0;
            qs_min_size = 0;
            l = 0;
            xk[0] = xk[1] = xs[0] = xs[1] = 0;
            continue;
        }

        const uint64_t c = static_cast<uint64_t>(base);
        const uint64_t rc = 3ULL - c;
        xk[0] = (xk[0] << 2 | c) & kmask;
        xk[1] = xk[1] >> 2 | rc << kshift;
        xs[0] = (xs[0] << 2 | c) & smask;
        xs[1] = xs[1] >> 2 | rc << sshift;
        const bool has_smer = (++l >= static_cast<size_t>(parameters.s));
        if (!has_smer)
        {
            continue;
        }

        const uint64_t ys = xs[0] < xs[1] ? xs[0] : xs[1];
        const uint64_t hash_s = syncmer_packed_hash(ys);

        const size_t qs_tail = (qs_head + qs_size) & kQueueMask;
        qs[qs_tail] = hash_s;
        ++qs_size;

        while (qs_min_size != 0)
        {
            const size_t min_tail = (qs_min_head + qs_min_size - 1) & kQueueMask;
            if (qs_min_candidates[min_tail] <= hash_s)
            {
                break;
            }
            --qs_min_size;
        }
        const size_t qs_min_tail = (qs_min_head + qs_min_size) & kQueueMask;
        qs_min_candidates[qs_min_tail] = hash_s;
        ++qs_min_size;

        if (qs_size < window_size)
        {
            continue;
        }
        if (qs_size > window_size)
        {
            const uint64_t front = qs[qs_head];
            qs_head = (qs_head + 1) & kQueueMask;
            --qs_size;
            if (qs_min_size != 0 && front == qs_min_candidates[qs_min_head])
            {
                qs_min_head = (qs_min_head + 1) & kQueueMask;
                --qs_min_size;
            }
        }

        const size_t target_slot = (qs_head + target_index) & kQueueMask;
        if (qs[target_slot] == qs_min_candidates[qs_min_head])
        {
            const uint64_t yk = xk[0] < xk[1] ? xk[0] : xk[1];
            out.push_back(
                Syncmer{syncmer_packed_hash(yk), i - static_cast<size_t>(parameters.k) + 1});
        }
    }
}

bool append_randstrobes_from_syncmers(const std::vector<Syncmer> & syncmers,
                                      RandstrobeParameters parameters,
                                      const StrobemerParams & strobe_params,
                                      std::vector<uint64_t> & out)
{
    if (syncmers.empty())
        return false;

    const size_t initial_size = out.size();
    out.reserve(out.size() + syncmers.size());

    RandstrobeIterator iterator(syncmers, parameters);
    while (iterator.has_next())
    {
        out.push_back(make_strobemer_hash(iterator.next(), strobe_params));
    }

    return out.size() > initial_size;
}
#endif

std::vector<uint64_t> compute_hashes_syncmer_(const std::vector<seqan3::dna4> & seq,
                                              const SyncmerParams & sp)
{
    const std::array<size_t, 1> pos{static_cast<size_t>(sp.pos)};
    const std::span<const size_t> span_pos{pos};
    return chimera::syncmer::compute_hashes(seq, sp.s, sp.k, span_pos, sp.seed, sp.canonical);
}

void compute_hashes_syncmer_append(const std::vector<seqan3::dna4> &seq,
                                   const SyncmerParams &sp,
                                   std::vector<uint64_t> &out)
{
    const std::array<size_t, 1> pos{static_cast<size_t>(sp.pos)};
    const std::span<const size_t> span_pos{pos};
    chimera::syncmer::compute_hashes_append(seq, sp.s, sp.k, span_pos, sp.seed, sp.canonical, out);
}

bool compute_hashes_strobemer_append_(const std::vector<seqan3::dna4> & seq,
                                      const StrobemerParams & params,
                                      std::vector<uint64_t> & out)
{
#ifdef CHIMERA_HAS_STROBEMERS
    if (params.order != 2)
        return false;

    const size_t sequence_length = seq.size();
    if (sequence_length < params.k || params.k < 8)
        return false;

    const std::span<const uint8_t> seq_span{
        reinterpret_cast<const uint8_t *>(seq.data()), seq.size()};
    const size_t initial_size = out.size();
    try
    {
        const int k = static_cast<int>(params.k);
        const int s_value = derive_syncmer_s(k);

        const unsigned w_min = std::max<unsigned>(params.w_min, 1u);
        const unsigned w_max = std::max<unsigned>(params.w_max, w_min);
        const int max_dist = derive_max_dist(k, w_max);
        const uint64_t main_hash_mask = ~0ULL << (9 + kDefaultAuxLen);

        SyncmerParameters syncmer_params{k, s_value};
        RandstrobeParameters rand_params{kDefaultRandstrobeQ,
                                         max_dist,
                                         w_min,
                                         w_max,
                                         main_hash_mask};

        std::vector<Syncmer> syncmers;
        const size_t syncmer_reserve =
            sequence_length / std::max<size_t>(1, static_cast<size_t>(k - s_value + 1));
        if (syncmer_reserve > 0)
            syncmers.reserve(syncmer_reserve);

        collect_syncmers_encoded(seq_span, syncmer_params, syncmers);

        if (syncmers.empty())
            return false;

        const bool appended =
            append_randstrobes_from_syncmers(syncmers, rand_params, params, out);

        return appended;
    }
    catch (const std::exception &)
    {
        out.resize(initial_size);
        return false;
    }
#else
    (void)seq;
    (void)params;
    (void)out;
    return false;
#endif
}

std::vector<uint64_t> compute_hashes_strobemer_(const std::vector<seqan3::dna4> & seq,
                                                const StrobemerParams & params)
{
    std::vector<uint64_t> hashes;
    compute_hashes_strobemer_append_(seq, params, hashes);
    return hashes;
}

} // namespace

std::vector<uint64_t> compute_hashes(const std::vector<seqan3::dna4> & seq, const Params & p)
{
    switch (p.method)
    {
    case Method::Syncmer:
        return compute_hashes_syncmer_(seq, p.sync);
    case Method::Strobemer:
        return compute_hashes_strobemer_(seq, p.strobe);
    case Method::Auto: {
        auto preferred = p;
        preferred.method = Method::Strobemer;
        auto hashes = compute_hashes_strobemer_(seq, preferred.strobe);
        if (!hashes.empty())
            return hashes;
        preferred.method = Method::Syncmer;
        return compute_hashes_syncmer_(seq, preferred.sync);
    }
    }
    return {};
}

void compute_hashes_append(const std::vector<seqan3::dna4> &seq,
                           const Params &p,
                           std::vector<uint64_t> &out)
{
    switch (p.method)
    {
    case Method::Syncmer:
        compute_hashes_syncmer_append(seq, p.sync, out);
        break;
    case Method::Strobemer:
        compute_hashes_strobemer_append_(seq, p.strobe, out);
        break;
    case Method::Auto: {
        auto preferred = p;
        preferred.method = Method::Strobemer;
        const size_t initial_size = out.size();
        if (compute_hashes_strobemer_append_(seq, preferred.strobe, out) && out.size() > initial_size)
            return;
        preferred.method = Method::Syncmer;
        compute_hashes_syncmer_append(seq, preferred.sync, out);
        break;
    }
    }
}

Params auto_params_from_readlen(size_t read_len)
{
    Params params{};
    params.method = Method::Strobemer;

    if (read_len < 80)
    {
        params.method = Method::Syncmer;
    }
    else if (read_len < 500)
    {
        params.strobe = {.k = 28, .order = 2, .w_min = 12, .w_max = 32, .seed = 0, .canonical = true};
    }
    else
    {
        params.strobe = {.k = 24, .order = 2, .w_min = 64, .w_max = 180, .seed = 0, .canonical = true};
    }

    if (params.method == Method::Strobemer && !strobemer_available())
    {
        params.method = Method::Syncmer;
    }

    return params;
}

bool strobemer_available() noexcept
{
#ifdef CHIMERA_HAS_STROBEMERS
    return true;
#else
    return false;
#endif
}

size_t min_required_length(const Params &p)
{
    if (p.method == Method::Strobemer)
    {
        const size_t k = p.strobe.k;
        const size_t w_min = p.strobe.w_min;
        return (2 * k) + w_min - 1;
    }
    return p.sync.k;
}

} // namespace chimera::feature
