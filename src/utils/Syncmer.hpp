#pragma once

#include <cstdint>
#include <ranges>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "strobemers/chimera_syncmer.hpp"

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
        config_{SyncmerConfigFrom(kmer_size, smer_size, positions, seed, canonical)}
    {
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
        using alphabet_t = std::remove_cvref_t<std::ranges::range_value_t<rng_t>>;
        static_assert(seqan3::semialphabet<alphabet_t>, "Sequence must contain SeqAn3 nucleotides");

        if (config_.k == 0)
            throw std::logic_error{"SyncmerRunner is not configured"};

        strobemers::chimera::SyncmerHasher hasher{config_};
        auto rank_of = [](auto const & symbol) -> uint8_t {
            return static_cast<uint8_t>(seqan3::to_rank(symbol));
        };

        hasher.process(std::ranges::begin(sequence),
                       std::ranges::end(sequence),
                       rank_of,
                       std::forward<callback_t>(callback));
    }

private:
    static strobemers::chimera::SyncmerConfig SyncmerConfigFrom(size_t kmer_size,
                                                                size_t smer_size,
                                                                std::span<const size_t> positions,
                                                                uint64_t seed,
                                                                bool canonical)
    {
        if (smer_size == 0)
            throw std::invalid_argument{"syncmer smer_size must be greater than 0"};
        if (kmer_size <= smer_size)
            throw std::invalid_argument{"syncmer kmer_size must be greater than smer_size"};
        if (positions.empty())
            throw std::invalid_argument{"syncmer positions must not be empty"};

        size_t const window_size = kmer_size - smer_size + 1;
        for (size_t pos : positions)
            if (pos >= window_size)
                throw std::invalid_argument{"syncmer position out of window range"};

        strobemers::chimera::SyncmerConfig cfg{};
        cfg.k = kmer_size;
        cfg.s = smer_size;
        cfg.seed = seed;
        cfg.canonical = canonical;
        cfg.positions.assign(positions.begin(), positions.end());
        return cfg;
    }

    strobemers::chimera::SyncmerConfig config_{};
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
