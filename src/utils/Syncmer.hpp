#pragma once

#include <cstdint>
#include <ranges>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "syncmer_hash.hpp"

namespace chimera::syncmer
{

struct SyncmerOptions
{
    size_t kmer_size{};
    size_t smer_size{};
    std::vector<size_t> positions{};
    uint64_t seed{0x8F3F73B5CF1C9ADEULL};
    bool canonical{true};
};

class SyncmerRunner
{
public:
    SyncmerRunner() = default;

    SyncmerRunner(size_t const smer_size,
                  size_t const kmer_size,
                  std::span<const size_t> positions,
                  uint64_t const seed,
                  bool const canonical) :
        kmer_size_{kmer_size},
        smer_size_{smer_size},
        seed_{seed},
        canonical_{canonical}
    {
        if (positions.empty())
            throw std::invalid_argument{"syncmer positions must not be empty"};
        positions_.reserve(positions.size());
        for (size_t const pos : positions)
        {
            if (pos >= kmer_size - smer_size + 1)
                throw std::invalid_argument{"syncmer position out of window range"};
            positions_.push_back(static_cast<int>(pos));
        }
    }

    SyncmerRunner(SyncmerOptions const & options,
                  std::span<const size_t> positions_override = {}) :
        SyncmerRunner(options.smer_size,
                      options.kmer_size,
                      positions_override.empty()
                          ? std::span<const size_t>(options.positions.data(), options.positions.size())
                          : positions_override,
                      options.seed,
                      options.canonical)
    {}

    template <std::ranges::forward_range rng_t, typename callback_t>
    void for_each(rng_t const & sequence, callback_t && callback) const
    {
        if (canonical_)
        {
            auto view = seqan3::detail::syncmer_hash_fn{}(sequence,
                                                          smer_size_,
                                                          kmer_size_,
                                                          positions_,
                                                          seed_);
            for (auto const value : view)
                std::forward<callback_t>(callback)(static_cast<uint64_t>(value));
        }
        else
        {
            auto view = seqan3::detail::syncmer_hash_no_reverse_fn{}(sequence,
                                                                     smer_size_,
                                                                     kmer_size_,
                                                                     positions_);
            for (auto const value : view)
                std::forward<callback_t>(callback)(static_cast<uint64_t>(value));
        }
    }

private:
    size_t kmer_size_{0};
    size_t smer_size_{0};
    std::vector<int> positions_{};
    seqan3::seed seed_{0x8F3F73B5CF1C9ADEULL};
    bool canonical_{true};
};

template <std::ranges::forward_range rng_t, typename callback_t>
inline void stream_hashes(rng_t const & sequence,
                          size_t const smer_size,
                          size_t const kmer_size,
                          std::span<const size_t> positions,
                          uint64_t const seed,
                          bool const canonical,
                          callback_t && callback)
{
    using alphabet_t = std::remove_cvref_t<std::ranges::range_value_t<rng_t>>;
    static_assert(seqan3::semialphabet<alphabet_t>, "Sequence must contain SeqAn3 nucleotides");

    if (smer_size == 0)
        throw std::invalid_argument{"syncmer smer_size must be greater than 0"};
    if (kmer_size <= smer_size)
        throw std::invalid_argument{"syncmer kmer_size must be greater than smer_size"};

    SyncmerRunner runner{smer_size, kmer_size, positions, seed, canonical};
    runner.for_each(sequence, std::forward<callback_t>(callback));
}

template <std::ranges::forward_range rng_t>
inline std::vector<uint64_t> compute_hashes(rng_t const & sequence,
                                            size_t const smer_size,
                                            size_t const kmer_size,
                                            std::span<const size_t> positions,
                                            uint64_t const seed = 0x8F3F73B5CF1C9ADEULL,
                                            bool const canonical = true)
{
    std::vector<uint64_t> values;
    stream_hashes(sequence,
                  smer_size,
                  kmer_size,
                  positions,
                  seed,
                  canonical,
                  [&](uint64_t hash)
                  {
                      values.push_back(hash);
                  });
    return values;
}

inline std::vector<uint64_t> compute_hashes(seqan3::dna5_vector const & sequence,
                                            SyncmerOptions const & options)
{
    return compute_hashes(sequence,
                          options.smer_size,
                          options.kmer_size,
                          options.positions,
                          options.seed,
                          options.canonical);
}

} // namespace chimera::syncmer
