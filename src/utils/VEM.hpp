/*
 * -----------------------------------------------------------------------------
 * Filename:      VEM.hpp
 * -----------------------------------------------------------------------------
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <classifyConfig.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ChimeraClassify {

struct VEMOptions {
	double eps = 1e-9;
	double alpha = 1e-6;
	double temp = 1.0;
	double prior_strength = 0.0;
	double coexist_penalty = 0.0;
};

namespace vem_detail {
	constexpr const char* kUnclassifiedTaxid = "unclassified";

	inline bool is_unclassified(const std::string& taxid) {
		return taxid == kUnclassifiedTaxid;
	}

	inline double log_sum_exp(const std::vector<double>& values) {
		double max_val = -std::numeric_limits<double>::infinity();
		for (double v : values) {
			max_val = std::max(max_val, v);
		}
		if (!std::isfinite(max_val)) {
			return max_val;
		}
		double sum = 0.0;
		for (double v : values) {
			sum += std::exp(v - max_val);
		}
		return max_val + std::log(sum);
	}
}

inline std::pair<std::vector<classifyResult>, std::unordered_map<std::string, double>>
VEMAlgorithm(const std::vector<classifyResult>& input,
	size_t max_iter,
	double tol,
	const VEMOptions& options,
	const std::unordered_map<std::string, double>* prior_scale = nullptr) {
	std::vector<classifyResult> results = input;
	(void)tol;
	std::unordered_set<std::string> taxid_set;
	for (const auto& record : results) {
		for (const auto& [taxid, count] : record.taxidCount) {
			if (vem_detail::is_unclassified(taxid)) {
				continue;
			}
			taxid_set.insert(taxid);
		}
	}

	std::vector<std::string> taxid_list;
	taxid_list.reserve(taxid_set.size());
	for (const auto& taxid : taxid_set) {
		taxid_list.emplace_back(taxid);
	}

	std::unordered_map<std::string, double> abundance_prior;
	abundance_prior.reserve(taxid_list.size());
	for (const auto& taxid : taxid_list) {
		abundance_prior[taxid] = 0.0;
	}
	double abundance_sum = 0.0;
	auto resolve_scale = [&](const std::string& taxid) -> double {
		if (!prior_scale) {
			return 1.0;
		}
		auto it = prior_scale->find(taxid);
		if (it == prior_scale->end() || !(it->second > 0.0)) {
			return 1.0;
		}
		return it->second;
	};
	for (const auto& record : results) {
		for (const auto& [taxid, count] : record.taxidCount) {
			if (vem_detail::is_unclassified(taxid)) {
				continue;
			}
			double scale = resolve_scale(taxid);
			double contrib = count / scale;
			abundance_prior[taxid] += contrib;
			abundance_sum += contrib;
		}
	}
	if (abundance_sum <= 0.0) {
		double uniform = 1.0 / static_cast<double>(std::max<size_t>(1, taxid_list.size()));
		for (auto& entry : abundance_prior) {
			entry.second = uniform;
		}
	} else {
		constexpr double prior_exponent = 0.75;
		double norm = 0.0;
		for (auto& entry : abundance_prior) {
			double value = std::max(entry.second, options.eps);
			double smooth = std::pow(value, prior_exponent);
			entry.second = smooth;
			norm += smooth;
		}
		norm = std::max(norm, options.eps * static_cast<double>(taxid_list.size()));
		for (auto& entry : abundance_prior) {
			entry.second /= norm;
		}
	}

	std::unordered_map<std::string, double> pi;
	if (taxid_list.empty()) {
		return { results, pi };
	}

	double uniform = 1.0 / static_cast<double>(taxid_list.size());
	for (const auto& taxid : taxid_list) {
		double prior = abundance_prior[taxid];
		if (!(prior > 0.0)) {
			prior = uniform;
		}
		pi[taxid] = prior;
	}

	size_t iteration = 0;
	while (iteration < max_iter) {
		std::unordered_map<std::string, double> expected_log_pi;
		for (const auto& taxid : taxid_list) {
			double prior_mix = std::max(pi[taxid], options.eps);
			if (options.prior_strength > 0.0) {
				double blended = (pi[taxid] + options.prior_strength * abundance_prior[taxid]) /
				                (1.0 + options.prior_strength);
				prior_mix = std::max(blended, options.eps);
			}
			expected_log_pi[taxid] = std::log(prior_mix);
		}

#ifdef _OPENMP
		int num_threads = omp_get_max_threads();
#else
		int num_threads = 1;
#endif
		std::vector<std::unordered_map<std::string, double>> thread_expected(num_threads);

#ifdef _OPENMP
#pragma omp parallel
#endif
		{
#ifdef _OPENMP
			int tid = omp_get_thread_num();
#else
			int tid = 0;
#endif
			auto& local_expected = thread_expected[tid];

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
			for (size_t i = 0; i < results.size(); ++i) {
				const auto& source = input[i];
				if (source.evaluated <= 0.0) {
					continue;
				}
				auto& destination = results[i];
				destination.posteriors.clear();

				std::vector<std::pair<std::string, double>> candidates;
				double max_count = 0.0;
				double second_count = 0.0;
				std::string best_taxid;
				for (const auto& [taxid, count] : source.taxidCount) {
					if (vem_detail::is_unclassified(taxid)) {
						continue;
					}
					double c = count / source.evaluated;
					candidates.emplace_back(taxid, c);
					if (c > max_count) {
						second_count = max_count;
						max_count = c;
						best_taxid = taxid;
					} else if (c > second_count) {
						second_count = c;
					}
				}

				if (candidates.empty()) {
					continue;
				}

				double denom = max_count + options.eps * static_cast<double>(candidates.size());
				if (denom <= 0.0) {
					denom = options.eps * static_cast<double>(candidates.size());
				}
				double penalty_scale = 0.0;
				if (options.coexist_penalty > 0.0 && second_count > 0.0) {
					double ratio_est = max_count /
					                 std::max(second_count, options.eps);
					double tightness = std::clamp((1.5 - ratio_est) / 0.5, 0.0, 1.0);
					penalty_scale = options.coexist_penalty * tightness;
				}

				std::vector<std::pair<std::string, double>> log_components;
				log_components.reserve(candidates.size());
				for (const auto& [taxid, c] : candidates) {
					auto it_pi = expected_log_pi.find(taxid);
					if (it_pi == expected_log_pi.end()) {
						continue;
					}
					double likelihood = (c + options.eps) / denom;
					double log_likelihood = options.temp * std::log(std::max(likelihood, options.eps));
					double coexist = 0.0;
					if (!best_taxid.empty() && taxid != best_taxid) {
						coexist = penalty_scale;
					}
					log_components.emplace_back(taxid, it_pi->second + log_likelihood - coexist);
				}

				if (log_components.empty()) {
					continue;
				}

				std::vector<double> log_values;
				log_values.reserve(log_components.size());
				for (const auto& entry : log_components) {
					log_values.emplace_back(entry.second);
				}

				double normalizer = vem_detail::log_sum_exp(log_values);
				std::vector<std::pair<std::string, double>> posterior;
				posterior.reserve(log_components.size());

				for (const auto& entry : log_components) {
					double q = std::exp(entry.second - normalizer);
					posterior.emplace_back(entry.first, q);
					local_expected[entry.first] += q;
				}

				destination.posteriors = std::move(posterior);
			}
		}

		std::unordered_map<std::string, double> expected_counts;
		double sum_expected = 0.0;
		for (auto& local : thread_expected) {
			for (const auto& [taxid, value] : local) {
				double& slot = expected_counts[taxid];
				slot += value;
				sum_expected += value;
			}
		}

		double prior_mass = std::max(0.0, options.prior_strength);
		double uniform_mass = options.alpha * static_cast<double>(taxid_list.size());
		double pseudo_mass = (prior_mass > 0.0) ? prior_mass : uniform_mass;
		double denominator = pseudo_mass + sum_expected;
		if (denominator <= 0.0) {
			denominator = 1.0;
		}

		for (const auto& taxid : taxid_list) {
			double expected = 0.0;
			auto it = expected_counts.find(taxid);
			if (it != expected_counts.end()) {
				expected = it->second;
			}
			double prior_component = (prior_mass > 0.0)
			                             ? prior_mass * abundance_prior[taxid]
			                             : options.alpha;
			double updated = (expected + prior_component) / denominator;
			pi[taxid] = updated;
		}

		++iteration;
	}

	return { results, pi };
}

} // namespace ChimeraClassify
