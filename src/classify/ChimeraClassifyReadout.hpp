#pragma once

#include "ChimeraClassify.hpp"

namespace ChimeraClassify {
struct NcbiTaxdump;

namespace readout {

// Selective readout contract:
// 1. postEmDecision may abstain when posterior/sample weight is unsafe.
// 2. Only posterior_weight abstentions may be released by sample-level evidence.
// 3. Accepted/released reads may switch only within read-local candidates.
// 4. Sample-mixture posterior is a calibration surface, not a majority label.
// 5. Readout must not create taxa absent from sampleMixturePosteriors.
bool apply_selective_readout(classifyResult &result,
                             const NcbiTaxdump *ncbiTaxdump);

bool apply_deferred_abstention_release(classifyResult &result,
                                       const NcbiTaxdump *ncbiTaxdump);

bool apply_local_correction(classifyResult &result,
                            const NcbiTaxdump *ncbiTaxdump);

} // namespace readout
} // namespace ChimeraClassify
