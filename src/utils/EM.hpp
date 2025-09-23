/*
 * -----------------------------------------------------------------------------
 * Filename:      EM.hpp
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

struct EMOptions {
	double eps = 1e-9;
	double alpha = 1e-6;
	double temp = 1.0;
};

namespace detail {
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
EMAlgorithm(const std::vector<classifyResult>& input,
	size_t max_iter,
	double tol,
	const EMOptions& options) {
	std::vector<classifyResult> results = input;
	std::unordered_set<std::string> taxid_set;
	for (const auto& record : results) {
		for (const auto& [taxid, count] : record.taxidCount) {
			if (detail::is_unclassified(taxid)) {
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

	std::unordered_map<std::string, double> pi;
	if (taxid_list.empty()) {
		return { results, pi };
	}

	double uniform = 1.0 / static_cast<double>(taxid_list.size());
	for (const auto& taxid : taxid_list) {
		pi[taxid] = uniform;
	}

	size_t iteration = 0;
	while (iteration < max_iter) {
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
				for (const auto& [taxid, count] : source.taxidCount) {
					if (detail::is_unclassified(taxid)) {
						continue;
					}
					double c = static_cast<double>(count) / source.evaluated;
					candidates.emplace_back(taxid, c);
					max_count = std::max(max_count, c);
				}

				if (candidates.empty()) {
					continue;
				}

				double denom = max_count + options.eps * static_cast<double>(candidates.size());
				if (denom <= 0.0) {
					denom = options.eps * static_cast<double>(candidates.size());
				}

				std::vector<std::pair<std::string, double>> log_components;
				log_components.reserve(candidates.size());
				for (const auto& [taxid, c] : candidates) {
					auto it = pi.find(taxid);
					if (it == pi.end()) {
						continue;
					}
					double likelihood = (c + options.eps) / denom;
					double log_likelihood = options.temp * std::log(std::max(likelihood, options.eps));
					double prior = std::max(it->second, options.eps);
					double log_prior = std::log(prior);
					log_components.emplace_back(taxid, log_prior + log_likelihood);
				}

				if (log_components.empty()) {
					continue;
				}

				std::vector<double> log_values;
				log_values.reserve(log_components.size());
				for (const auto& entry : log_components) {
					log_values.emplace_back(entry.second);
				}

				double normalizer = detail::log_sum_exp(log_values);
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

		double denominator = options.alpha * static_cast<double>(taxid_list.size()) + sum_expected;
		if (denominator <= 0.0) {
			denominator = 1.0;
		}

		double diff = 0.0;
		for (const auto& taxid : taxid_list) {
			double expected = 0.0;
			auto it = expected_counts.find(taxid);
			if (it != expected_counts.end()) {
				expected = it->second;
			}
			double updated = (expected + options.alpha) / denominator;
			diff += std::abs(pi[taxid] - updated);
			pi[taxid] = updated;
		}

		if (diff < tol) {
			break;
		}

		++iteration;
	}

	return { results, pi };
}

} // namespace ChimeraClassify
