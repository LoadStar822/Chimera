#pragma once

#include <cstdint>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

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

std::vector<uint64_t> compute_hashes(const std::vector<seqan3::dna4> & seq, const Params & p);
void compute_hashes_append(const std::vector<seqan3::dna4> &seq,
                           const Params &p,
                           std::vector<uint64_t> &out);

bool strobemer_available() noexcept;

size_t min_required_length(const Params &p);

} // namespace chimera::feature
