#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "classify/post_em_debug_dump.hpp"

namespace {

bool expect_true(const std::string &name, bool got, std::string &message) {
  if (got) {
    return true;
  }
  message = name + " 期望 true, 实际 false";
  return false;
}

bool expect_false(const std::string &name, bool got, std::string &message) {
  if (!got) {
    return true;
  }
  message = name + " 期望 false, 实际 true";
  return false;
}

bool expect_eq_size_t(const std::string &name, std::size_t got, std::size_t want,
                      std::string &message) {
  if (got == want) {
    return true;
  }
  message = name + " 期望 " + std::to_string(want) + ", 实际 " +
            std::to_string(got);
  return false;
}

bool expect_near(const std::string &name, double got, double want, double eps,
                 std::string &message) {
  if (std::abs(got - want) <= eps) {
    return true;
  }
  message = name + " 期望约 " + std::to_string(want) + ", 实际 " +
            std::to_string(got);
  return false;
}

} // namespace

int main() {
  int failures = 0;
  std::vector<std::string> failure_messages;

  {
    std::string message;
    std::vector<std::pair<std::string, double>> posterior{
        {"A", 0.6}, {"B", 0.3}, {"C", 0.1}};
    auto hit = ChimeraClassify::lookup_taxid_in_topk(posterior, "B", 2);
    bool ok = expect_true("lookup_found", hit.found, message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(message);
    }
    message.clear();
    ok = expect_eq_size_t("lookup_rank", hit.rank, 1, message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(message);
    }
    message.clear();
    ok = expect_near("lookup_prob", hit.prob, 0.3, 1e-12, message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(message);
    }
  }

  {
    std::string message;
    std::vector<std::pair<std::string, double>> posterior{
        {"A", 0.6}, {"B", 0.3}, {"C", 0.1}};
    auto hit = ChimeraClassify::lookup_taxid_in_topk(posterior, "C", 2);
    bool ok = expect_false("lookup_not_found_when_outside_topk", hit.found, message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(message);
    }
  }

  if (failures == 0) {
    std::cout << "All tests passed." << std::endl;
    return 0;
  }

  std::cerr << "Failures: " << failures << std::endl;
  for (const auto &msg : failure_messages) {
    std::cerr << "- " << msg << std::endl;
  }
  return 1;
}

