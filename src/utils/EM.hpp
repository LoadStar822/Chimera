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
	double prior_strength = 0.0;
	double prune_ratio = 1e-3; // relative threshold to max_expected for sparsity
	double conf_power = 0.0;   // confidence weighting exponent; 0 disables
};

namespace detail {
	constexpr const char* kUnclassifiedTaxid = "unclassified";

	inline bool is_unclassified(const std::string& taxid) {
		return taxid == kUnclassifiedTaxid;
	}
}

inline std::pair<std::vector<classifyResult>, std::unordered_map<std::string, double>>
EMAlgorithm(std::vector<classifyResult> input,
	size_t max_iter,
	double tol,
	const EMOptions& options,
	const std::unordered_map<std::string, double>* prior_scale = nullptr) {
	std::vector<classifyResult> results = std::move(input);
	(void)tol;
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

	if (taxid_list.empty()) {
		std::unordered_map<std::string, double> pi;
		return { results, pi };
	}

	const size_t taxid_count = taxid_list.size();
	std::unordered_map<std::string, size_t> taxid_to_idx;
	taxid_to_idx.reserve(taxid_count);
	for (size_t i = 0; i < taxid_count; ++i) {
		taxid_to_idx.emplace(taxid_list[i], i);
	}

	struct PreparedCandidate {
		size_t idx{0};
		double log_likelihood{0.0};
	};
	struct PreparedRead {
		double base_weight{1.0};
		bool valid{false};
		std::vector<PreparedCandidate> candidates;
	};

	auto normalize_distribution = [&](std::vector<double>& dist) {
		double sum = 0.0;
		for (double v : dist) {
			sum += v;
		}
		if (!(sum > 0.0)) {
			double uniform = 1.0 / static_cast<double>(taxid_count);
			for (double& v : dist) {
				v = uniform;
			}
			return;
		}
		double inv_sum = 1.0 / sum;
		for (double& v : dist) {
			v *= inv_sum;
		}
	};

	// 预处理每条 read 的候选集合（taxid->idx 映射、常量似然项）；
	// 后续迭代只更新先验项，避免重复字符串哈希与重复分母计算。
	std::vector<PreparedRead> prepared(results.size());
	std::vector<double> weighted_evidence(taxid_count, 0.0);
	for (size_t i = 0; i < results.size(); ++i) {
		const auto& source = results[i];
		auto& prep = prepared[i];
		if (source.evaluated <= 0.0 || source.taxidCount.empty()) {
			continue;
		}

		struct TempCandidate {
			size_t idx{0};
			double c{0.0}; // count / evaluated
			double raw{0.0};
		};
		std::vector<TempCandidate> tmp;
		tmp.reserve(source.taxidCount.size());

		double total_count = 0.0;
		double max_count = 0.0;
		for (const auto& kv : source.taxidCount) {
			if (detail::is_unclassified(kv.first)) {
				continue;
			}
			auto it_idx = taxid_to_idx.find(kv.first);
			if (it_idx == taxid_to_idx.end()) {
				continue;
			}
			double c = kv.second / source.evaluated;
			tmp.push_back(TempCandidate{it_idx->second, c, kv.second});
			total_count += kv.second;
			if (c > max_count) {
				max_count = c;
			}
		}
		if (tmp.empty()) {
			continue;
		}

		if (total_count > 0.0) {
			for (const auto& cand : tmp) {
				weighted_evidence[cand.idx] += cand.raw / total_count;
			}
		}

		double denom = max_count + options.eps * static_cast<double>(tmp.size());
		if (denom <= 0.0) {
			denom = options.eps * static_cast<double>(tmp.size());
		}
		prep.candidates.reserve(tmp.size());
		for (const auto& cand : tmp) {
			double likelihood = (cand.c + options.eps) / denom;
			double log_likelihood =
				options.temp * std::log(std::max(likelihood, options.eps));
			prep.candidates.push_back(PreparedCandidate{cand.idx, log_likelihood});
		}

		prep.base_weight = source.evaluated;
		if (!(prep.base_weight > 0.0)) {
			prep.base_weight = 1.0;
		}
		prep.valid = !prep.candidates.empty();
	}

	auto get_scale = [&](size_t idx) -> double {
		if (!prior_scale) {
			return 1.0;
		}
		auto it = prior_scale->find(taxid_list[idx]);
		if (it == prior_scale->end() || !(it->second > 0.0)) {
			return 1.0;
		}
		return std::sqrt(it->second);
	};

	std::vector<double> pi_vec(taxid_count, 0.0);
	double background_prior = 1e-9;
	for (size_t idx = 0; idx < taxid_count; ++idx) {
		double scale = get_scale(idx);
		double init_val = weighted_evidence[idx] * scale;
		if (init_val <= 0.0 && prior_scale && scale > 0.0) {
			init_val = 1e-5;
		}
		if (init_val <= 0.0) {
			init_val = background_prior;
		}
		pi_vec[idx] = init_val;
	}

	normalize_distribution(pi_vec);
	std::vector<double> abundance_prior = pi_vec;

#ifdef _OPENMP
	const int num_threads = omp_get_max_threads();
#else
	const int num_threads = 1;
#endif
	std::vector<std::vector<double>> thread_expected(
		num_threads, std::vector<double>(taxid_count, 0.0));
	std::vector<double> expected_counts(taxid_count, 0.0);

	size_t iteration = 0;
	while (iteration < max_iter) {
		const bool final_iter = (iteration + 1 == max_iter);
		for (auto& local : thread_expected) {
			std::fill(local.begin(), local.end(), 0.0);
		}

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
			std::vector<std::pair<size_t, double>> log_components;
			std::vector<double> q_values;

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
			for (size_t i = 0; i < results.size(); ++i) {
				auto& destination = results[i];
				if (final_iter) {
					destination.posteriors.clear();
				}

				const auto& prep = prepared[i];
				if (!prep.valid) {
					continue;
				}

				log_components.clear();
				log_components.reserve(prep.candidates.size());
				q_values.clear();
				q_values.reserve(prep.candidates.size());
				double max_log = -std::numeric_limits<double>::infinity();
				for (const auto& cand : prep.candidates) {
					double prior = std::max(pi_vec[cand.idx], options.eps);
					double score = std::log(prior) + cand.log_likelihood;
					log_components.emplace_back(cand.idx, score);
					if (score > max_log) {
						max_log = score;
					}
				}
				if (log_components.empty()) {
					continue;
				}
				if (!std::isfinite(max_log)) {
					continue;
				}

				double sum_exp = 0.0;
				for (const auto& entry : log_components) {
					sum_exp += std::exp(entry.second - max_log);
				}
				double normalizer = max_log + std::log(sum_exp);
				double max_post = 0.0;
				for (const auto& entry : log_components) {
					double q = std::exp(entry.second - normalizer);
					q_values.emplace_back(q);
					if (q > max_post) {
						max_post = q;
					}
				}

				double conf_w = 1.0;
				if (options.conf_power > 0.0 && max_post > 0.0) {
					conf_w = std::pow(max_post, options.conf_power);
				}
				double weight = prep.base_weight * conf_w;

				if (final_iter) {
					auto& posterior = destination.posteriors;
					posterior.reserve(log_components.size());
					for (size_t k = 0; k < log_components.size(); ++k) {
						size_t idx = log_components[k].first;
						double q = q_values[k];
						posterior.emplace_back(taxid_list[idx], q);
						local_expected[idx] += weight * q;
					}
				} else {
					for (size_t k = 0; k < log_components.size(); ++k) {
						size_t idx = log_components[k].first;
						double q = q_values[k];
						local_expected[idx] += weight * q;
					}
				}
			}
		}

		std::fill(expected_counts.begin(), expected_counts.end(), 0.0);
		double sum_expected = 0.0;
		for (const auto& local : thread_expected) {
			for (size_t idx = 0; idx < taxid_count; ++idx) {
				expected_counts[idx] += local[idx];
			}
		}
		for (double v : expected_counts) {
			sum_expected += v;
		}

		double prior_mass = std::max(0.0, options.prior_strength);
		double uniform_mass = options.alpha * static_cast<double>(taxid_count);
		double pseudo_mass = (prior_mass > 0.0) ? prior_mass : uniform_mass;
		double denominator = pseudo_mass + sum_expected;
		if (denominator <= 0.0) {
			denominator = 1.0;
		}

		double max_expected = 0.0;
		for (double v : expected_counts) {
			if (v > max_expected) {
				max_expected = v;
			}
		}
		double prune_thres = max_expected * options.prune_ratio;
		for (size_t idx = 0; idx < taxid_count; ++idx) {
			double expected = expected_counts[idx];
			if (expected < prune_thres) {
				expected = 0.0;
			}
			double prior_component = (prior_mass > 0.0)
			                             ? prior_mass * abundance_prior[idx]
			                             : options.alpha;
			pi_vec[idx] = (expected + prior_component) / denominator;
		}

		++iteration;
	}

	std::unordered_map<std::string, double> pi;
	pi.reserve(taxid_count);
	for (size_t idx = 0; idx < taxid_count; ++idx) {
		pi.emplace(taxid_list[idx], pi_vec[idx]);
	}
	return { results, pi };
}

} // namespace ChimeraClassify
