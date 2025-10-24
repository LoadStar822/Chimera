#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <ranges>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "strobemers/indexparameters.hpp"
#include "strobemers/randstrobes.hpp"
#include "xxhash.h"

namespace chimera::strobemer {

struct StrobeParams {
    uint16_t k;
    uint16_t s;
    uint16_t w_min;
    uint16_t w_max;
    uint16_t q;
    uint32_t max_dist;
    uint16_t aux_len;
    uint64_t seed;

    bool operator==(const StrobeParams& other) const {
        return k == other.k && s == other.s &&
               w_min == other.w_min && w_max == other.w_max &&
               q == other.q && max_dist == other.max_dist &&
               aux_len == other.aux_len && seed == other.seed;
    }
};

struct StrobeProfile {
    StrobeParams params;
    double weight;
};

inline uint16_t normalize_syncmer_span(uint16_t k, uint16_t candidate) {
    uint16_t maxAllowed = static_cast<uint16_t>(k > 1 ? k - 1 : 1);
    uint16_t minAllowed = static_cast<uint16_t>(std::min<uint16_t>(maxAllowed, 4));
    if (minAllowed > maxAllowed) {
        minAllowed = maxAllowed;
    }
    candidate = std::clamp(candidate, minAllowed, maxAllowed);
    auto parity_ok = [k](uint16_t s) -> bool {
        return ((k - s) & 1u) == 0;
    };
    if (parity_ok(candidate)) {
        return candidate;
    }
    int bestDiff = std::numeric_limits<int>::max();
    uint16_t best = candidate;
    for (uint16_t s = minAllowed; s <= maxAllowed; ++s) {
        if (!parity_ok(s)) {
            continue;
        }
        int diff = static_cast<int>(s) - static_cast<int>(candidate);
        diff = diff < 0 ? -diff : diff;
        if (diff < bestDiff) {
            bestDiff = diff;
            best = s;
            if (bestDiff == 0) {
                break;
            }
        }
    }
    return best;
}

inline void sanitize_params(StrobeParams& params) {
    params.s = normalize_syncmer_span(params.k, params.s);
    if (params.w_min == 0) {
        params.w_min = 1;
    }
    if (params.w_max == 0) {
        params.w_max = params.w_min;
    }
    if (params.w_max < params.w_min) {
        params.w_max = params.w_min;
    }
}

inline constexpr uint64_t encode_key(bool is_strobemer, uint64_t hash_value) {
    constexpr uint64_t type_bit = 1ull << 63;
    hash_value &= ~type_bit;
    return is_strobemer ? (hash_value | type_bit) : hash_value;
}

class StrobemerRunner {
public:
    explicit StrobemerRunner(StrobeParams params)
        : params_(params)
        , syncmer_params_(static_cast<int>(params.k), static_cast<int>(params.s))
        , rand_params_(params.q,
                       static_cast<int>(params.max_dist),
                       params.w_min,
                       params.w_max,
                       (~0ull) << (9 + std::min<uint16_t>(params.aux_len, static_cast<uint16_t>(27))))
    {}

    template <std::ranges::forward_range rng_t, typename Callback>
    void for_each(rng_t const& sequence, Callback&& callback) const {
        using alphabet_t = std::remove_cvref_t<std::ranges::range_value_t<rng_t>>;
        static_assert(seqan3::semialphabet<alphabet_t>, "sequence must contain SeqAn nucleotides");

        if (std::ranges::size(sequence) < params_.k) {
            return;
        }

        buffer_.assign(std::ranges::size(sequence), 'A');
        size_t idx = 0;
        for (auto const& nuc : sequence) {
            buffer_[idx++] = seqan3::to_char(nuc);
        }

        RandstrobeGenerator generator(buffer_, syncmer_params_, rand_params_);
        const Randstrobe end = generator.end();
        for (auto stro = generator.next(); stro != end; stro = generator.next()) {
            uint64_t canonical = stro.hash < stro.hash_revcomp ? stro.hash : stro.hash_revcomp;
            uint64_t mixed = XXH3_64bits_withSeed(&canonical, sizeof(canonical), params_.seed);
            std::forward<Callback>(callback)(mixed);
        }
    }

private:
    mutable std::string buffer_;
    StrobeParams params_;
    SyncmerParameters syncmer_params_;
    RandstrobeParameters rand_params_;
};

inline StrobeProfile auto_profile_from_read(size_t read_length, uint16_t kmer_size, uint64_t base_seed) {
    StrobeProfile profile{};
    if (read_length <= 250) {
        profile.params = StrobeParams{
            static_cast<uint16_t>(std::min<size_t>(kmer_size, 30)),
            static_cast<uint16_t>(std::min<size_t>(kmer_size > 2 ? kmer_size - 2 : kmer_size, 16)),
            25, 55, 8, 120, 15, base_seed
        };
        profile.weight = 1.8;
    } else if (read_length <= 1200) {
        profile.params = StrobeParams{
            static_cast<uint16_t>(std::min<size_t>(kmer_size, 32)),
            static_cast<uint16_t>(std::min<size_t>(kmer_size > 2 ? kmer_size - 2 : kmer_size, 18)),
            32, 72, 9, 180, 15, base_seed
        };
        profile.weight = 1.6;
    } else {
        profile.params = StrobeParams{
            static_cast<uint16_t>(std::min<size_t>(kmer_size, 40)),
            static_cast<uint16_t>(std::min<size_t>(kmer_size > 4 ? kmer_size - 4 : kmer_size, 20)),
            48, 128, 11, 500, 15, base_seed
        };
        profile.weight = 1.6;
    }
    if (profile.params.k < kmer_size) {
        profile.params.k = kmer_size;
    }
    sanitize_params(profile.params);
    return profile;
}

inline StrobeParams override_params(const StrobeParams& base, const StrobeParams& override) {
    StrobeParams out = base;
    if (override.w_min > 0) out.w_min = override.w_min;
    if (override.w_max > 0) out.w_max = override.w_max;
    if (override.q > 0) out.q = override.q;
    if (override.max_dist > 0) out.max_dist = override.max_dist;
    if (override.aux_len > 0) out.aux_len = override.aux_len;
    if (override.seed != 0) out.seed = override.seed;
    sanitize_params(out);
    return out;
}

} // namespace chimera::strobemer
