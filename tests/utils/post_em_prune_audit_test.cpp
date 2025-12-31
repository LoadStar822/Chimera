#include <cmath>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "classify/post_em_prune_audit.hpp"

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

bool expect_eq_str(const std::string &name, const std::string &got,
                   const std::string &want, std::string &message) {
  if (got == want) {
    return true;
  }
  message = name + " 期望 " + want + ", 实际 " + got;
  return false;
}

bool expect_eq_size(const std::string &name, size_t got, size_t want,
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
  if (std::fabs(got - want) <= eps) {
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

  using ChimeraClassify::PostEmPresenceLevel;

  auto always_unknown = [](const std::string &) {
    return PostEmPresenceLevel::kUnknown;
  };

  {
    std::string message;
    std::vector<std::pair<std::string, double>> full{{"A", 0.6}, {"B", 0.4}};
    std::unordered_map<std::string, double> weights{{"A", 1e-3}, {"B", 1e-6}};
    auto res = ChimeraClassify::prune_post_em_posterior(
        full, /*min_class_weight=*/0.0, weights, always_unknown,
        /*presence_pi_floor=*/1e-6, /*reject_factor=*/2.0);
    bool ok = expect_false("min_weight_0_disables_pruning", res.audit.pruning_active,
                           message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
    message.clear();
    ok = expect_eq_size("min_weight_0_keeps_full_size", res.posterior.size(), 2,
                        message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    std::vector<std::pair<std::string, double>> full{{"A", 0.6}, {"B", 0.4}};
    std::unordered_map<std::string, double> weights{{"A", 1e-3}, {"B", 1e-6}};
    auto res = ChimeraClassify::prune_post_em_posterior(
        full, /*min_class_weight=*/1e-4, weights, always_unknown,
        /*presence_pi_floor=*/1e-6, /*reject_factor=*/2.0);
    bool ok = expect_true("pruning_active_basic", res.audit.pruning_active, message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
    message.clear();
    ok = expect_true("pruned_any_when_tail_removed", res.audit.pruned_any, message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
    message.clear();
    ok = expect_false("top1_not_removed_when_top_survives", res.audit.top1_removed,
                      message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
    message.clear();
    ok = expect_eq_size("posterior_pruned_to_1", res.posterior.size(), 1, message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
    message.clear();
    ok = expect_eq_str("posterior_keeps_A", res.posterior[0].first, "A", message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
    message.clear();
    ok = expect_near("posterior_renorm_to_1", res.posterior[0].second, 1.0, 1e-12,
                     message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
    message.clear();
    ok = expect_near("dropped_mass_0p4", res.audit.dropped_mass, 0.4, 1e-12, message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    std::vector<std::pair<std::string, double>> full{{"A", 0.6}, {"B", 0.4}};
    std::unordered_map<std::string, double> weights{{"A", 1e-6}, {"B", 1e-3}};
    auto res = ChimeraClassify::prune_post_em_posterior(
        full, /*min_class_weight=*/1e-4, weights, always_unknown,
        /*presence_pi_floor=*/1e-6, /*reject_factor=*/2.0);
    bool ok = expect_true("top1_removed_when_A_pruned", res.audit.top1_removed,
                          message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
    message.clear();
    ok = expect_eq_str("posterior_switches_to_B", res.posterior[0].first, "B", message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
    message.clear();
    ok = expect_near("dropped_mass_0p6", res.audit.dropped_mass, 0.6, 1e-12, message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    std::vector<std::pair<std::string, double>> full{{"A", 0.6}, {"B", 0.4}};
    std::unordered_map<std::string, double> weights{{"A", 1e-6}, {"B", 1e-6}};
    auto res = ChimeraClassify::prune_post_em_posterior(
        full, /*min_class_weight=*/1e-4, weights, always_unknown,
        /*presence_pi_floor=*/1e-6, /*reject_factor=*/2.0);
    bool ok = expect_true("fallback_full_when_all_pruned", res.audit.fallback_full,
                          message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
    message.clear();
    ok = expect_eq_size("fallback_full_keeps_size", res.posterior.size(), 2, message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    std::vector<std::pair<std::string, double>> full{{"A", 0.6}, {"B", 0.4}};
    std::unordered_map<std::string, double> weights{{"A", 1e-3}, {"B", 5e-6}};
    auto presence = [](const std::string &tid) {
      return (tid == "B") ? PostEmPresenceLevel::kAccepted
                           : PostEmPresenceLevel::kUnknown;
    };
    auto res = ChimeraClassify::prune_post_em_posterior(
        full, /*min_class_weight=*/1e-4, weights, presence,
        /*presence_pi_floor=*/1e-6, /*reject_factor=*/2.0);
    bool ok = expect_false("accepted_taxon_survives_pruning", res.audit.pruned_any,
                           message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
    message.clear();
    ok = expect_eq_size("accepted_keeps_both", res.posterior.size(), 2, message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
  }

  {
    std::string message;
    std::vector<std::pair<std::string, double>> full{{"A", 0.6}, {"B", 0.4}};
    std::unordered_map<std::string, double> weights{{"A", 1.5e-4}, {"B", 1e-3}};
    auto presence = [](const std::string &tid) {
      return (tid == "A") ? PostEmPresenceLevel::kRejected
                           : PostEmPresenceLevel::kUnknown;
    };
    auto res = ChimeraClassify::prune_post_em_posterior(
        full, /*min_class_weight=*/1e-4, weights, presence,
        /*presence_pi_floor=*/1e-6, /*reject_factor=*/2.0);
    bool ok = expect_true("rejected_taxon_pruned_more_strictly",
                          res.audit.top1_removed, message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
    }
    message.clear();
    ok = expect_eq_str("rejected_top1_switches", res.posterior[0].first, "B", message);
    if (!ok) {
      ++failures;
      failure_messages.push_back(std::move(message));
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

