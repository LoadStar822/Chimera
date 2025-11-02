#pragma once

#include <cstdint>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

namespace chimera::feature
{

enum class Method : uint8_t
{
    Syncmer = 0,
    Strobemer = 1,
    Auto = 2
};

struct SyncmerParams
{
    uint8_t k{31};
    uint16_t s{16};
    uint16_t pos{7};
    uint64_t seed{0};
    bool canonical{true};
};

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
    Method method{Method::Syncmer};
    SyncmerParams sync{};
    StrobemerParams strobe{};
};

std::vector<uint64_t> compute_hashes(const std::vector<seqan3::dna4> & seq, const Params & p);
void compute_hashes_append(const std::vector<seqan3::dna4> &seq,
                           const Params &p,
                           std::vector<uint64_t> &out);

Params auto_params_from_readlen(size_t read_len);

bool strobemer_available() noexcept;

size_t min_required_length(const Params &p);

} // namespace chimera::feature
