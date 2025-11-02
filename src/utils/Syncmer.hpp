#pragma once

#include <cstdint>
#include <ranges>
#include <span>
#include <stdexcept>
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

template <std::ranges::forward_range rng_t>
inline std::vector<uint64_t> compute_hashes(rng_t const & sequence,
                                            size_t const smer_size,
                                            size_t const kmer_size,
                                            std::span<const size_t> positions,
                                            uint64_t const seed = 0x8F3F73B5CF1C9ADEULL,
                                            bool const canonical = true)
{
    using alphabet_t = std::remove_cvref_t<std::ranges::range_value_t<rng_t>>;
    static_assert(seqan3::semialphabet<alphabet_t>, "Sequence must contain SeqAn3 nucleotides");

    if (smer_size == 0)
        throw std::invalid_argument{"syncmer smer_size must be greater than 0"};
    if (kmer_size <= smer_size)
        throw std::invalid_argument{"syncmer kmer_size must be greater than smer_size"};

    if (positions.empty())
        throw std::invalid_argument{"syncmer positions must not be empty"};

    std::vector<int> converted_positions;
    converted_positions.reserve(positions.size());
    for (size_t const pos : positions)
    {
        if (pos >= kmer_size - smer_size + 1)
            throw std::invalid_argument{"syncmer position out of window range"};
        converted_positions.push_back(static_cast<int>(pos));
    }

    if (canonical)
    {
        auto view = seqan3::detail::syncmer_hash_fn{}(sequence,
                                                      smer_size,
                                                      kmer_size,
                                                      converted_positions,
                                                      seqan3::seed{seed});
        std::vector<uint64_t> values;
        for (auto const value : view)
            values.push_back(static_cast<uint64_t>(value));
        return values;
    }
    else
    {
        auto view = seqan3::detail::syncmer_hash_no_reverse_fn{}(sequence,
                                                                 smer_size,
                                                                 kmer_size,
                                                                 converted_positions);
        std::vector<uint64_t> values;
        for (auto const value : view)
            values.push_back(static_cast<uint64_t>(value));
        return values;
    }
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

template <std::ranges::forward_range rng_t>
inline void compute_hashes_append(rng_t const &sequence,
                                  size_t const smer_size,
                                  size_t const kmer_size,
                                  std::span<const size_t> positions,
                                  uint64_t const seed,
                                  bool const canonical,
                                  std::vector<uint64_t> &out)
{
    using alphabet_t = std::remove_cvref_t<std::ranges::range_value_t<rng_t>>;
    static_assert(seqan3::semialphabet<alphabet_t>, "Sequence must contain SeqAn3 nucleotides");

    if (smer_size == 0)
        throw std::invalid_argument{"syncmer smer_size must be greater than 0"};
    if (kmer_size <= smer_size)
        throw std::invalid_argument{"syncmer kmer_size must be greater than smer_size"};

    if (positions.empty())
        throw std::invalid_argument{"syncmer positions must not be empty"};

    std::vector<int> converted_positions;
    converted_positions.reserve(positions.size());
    for (size_t const pos : positions)
    {
        if (pos >= kmer_size - smer_size + 1)
            throw std::invalid_argument{"syncmer position out of window range"};
        converted_positions.push_back(static_cast<int>(pos));
    }

    if (canonical)
    {
        auto view = seqan3::detail::syncmer_hash_fn{}(sequence,
                                                      smer_size,
                                                      kmer_size,
                                                      converted_positions,
                                                      seqan3::seed{seed});
        for (auto const value : view)
            out.push_back(static_cast<uint64_t>(value));
    }
    else
    {
        auto view = seqan3::detail::syncmer_hash_no_reverse_fn{}(sequence,
                                                                 smer_size,
                                                                 kmer_size,
                                                                 converted_positions);
        for (auto const value : view)
            out.push_back(static_cast<uint64_t>(value));
    }
}

} // namespace chimera::syncmer
