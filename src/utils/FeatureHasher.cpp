#include "utils/FeatureHasher.hpp"

#include <algorithm>
#include <array>
#include <span>
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

    const int k = static_cast<int>(params.k);
    const int s_value = derive_syncmer_s(k);

    const unsigned w_min = std::max<unsigned>(params.w_min, 1u);
    const unsigned w_max = std::max<unsigned>(params.w_max, w_min);
    const int max_dist = derive_max_dist(k, w_max);
    const uint64_t main_hash_mask = ~0ULL << (9 + kDefaultAuxLen);

    const size_t initial_size = out.size();
    try
    {
        SyncmerParameters syncmer_params{k, s_value};
        RandstrobeParameters rand_params{kDefaultRandstrobeQ,
                                         max_dist,
                                         w_min,
                                         w_max,
                                         main_hash_mask};

        RandstrobeGenerator generator(seq_span, syncmer_params, rand_params);
        const Randstrobe sentinel = generator.end();

        const size_t expected = sequence_length / std::max<int>(1, k);
        if (expected > 0)
            out.reserve(out.size() + expected);

        for (auto strobe = generator.next(); strobe != sentinel; strobe = generator.next())
        {
            if (strobe.hash == 0 && strobe.hash_revcomp == 0 &&
                strobe.strobe1_pos == 0 && strobe.strobe2_pos == 0)
            {
                break;
            }
            out.push_back(make_strobemer_hash(strobe, params));
        }

        return out.size() > initial_size;
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
