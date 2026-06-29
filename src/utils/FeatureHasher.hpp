#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

#ifdef CHIMERA_HAS_STROBEMERS
#include <strobemers/randstrobes.hpp>
#endif

namespace chimera::feature
{

struct StrobemerParams
{
    uint8_t k{28};
    uint8_t order{2};
    uint16_t w_min{12};
    uint16_t w_max{32};
    uint64_t seed{0};
    bool canonical{true};
};

struct Params
{
    StrobemerParams strobe{};
};

struct FeatureHashScratch
{
#ifdef CHIMERA_HAS_STROBEMERS
    std::vector<Syncmer> syncmers;
#endif
};

std::vector<uint64_t> compute_hashes(const std::vector<seqan3::dna4> & seq, const Params & p);
void compute_hashes_append(const std::vector<seqan3::dna4> &seq,
	                           const Params &p,
	                           std::vector<uint64_t> &out);
void compute_hashes_append(const std::vector<seqan3::dna4> &seq,
                           const Params &p,
                           std::vector<uint64_t> &out,
                           FeatureHashScratch &scratch);

bool strobemer_available() noexcept;

size_t min_required_length(const Params &p);

} // namespace chimera::feature
