#include <iostream>
#include <string>

#include "classify/ChimeraClassifyCommon.hpp"

namespace {

bool expect_bool(const std::string &name, bool got, bool expected,
                 std::string &message) {
  if (got == expected) {
    return true;
  }
  message = name + " 期望 " + (expected ? "true" : "false") + ", 实际 " +
            (got ? "true" : "false");
  return false;
}

} // namespace

int main() {
  int failures = 0;
  std::string message;

  if (!expect_bool("xor_00",
                   ChimeraClassify::local_contrast_xor_gate(false, false),
                   false, message)) {
    ++failures;
    std::cerr << message << std::endl;
  }
  if (!expect_bool("xor_11",
                   ChimeraClassify::local_contrast_xor_gate(true, true), false,
                   message)) {
    ++failures;
    std::cerr << message << std::endl;
  }
  if (!expect_bool("xor_10",
                   ChimeraClassify::local_contrast_xor_gate(true, false), true,
                   message)) {
    ++failures;
    std::cerr << message << std::endl;
  }
  if (!expect_bool("xor_01",
                   ChimeraClassify::local_contrast_xor_gate(false, true), true,
                   message)) {
    ++failures;
    std::cerr << message << std::endl;
  }

  if (failures == 0) {
    std::cout << "All tests passed." << std::endl;
    return 0;
  }
  std::cerr << "Failures: " << failures << std::endl;
  return 1;
}

