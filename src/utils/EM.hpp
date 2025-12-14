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
	double coexist_penalty = 0.0;
	double prune_ratio = 1e-3; // relative threshold to max_expected for sparsity
	double conf_power = 0.0;   // confidence weighting exponent; 0 disables
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
	const EMOptions& options,
	const std::unordered_map<std::string, double>* prior_scale = nullptr) {
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

	if (taxid_list.empty()) {
		std::unordered_map<std::string, double> pi;
		return { results, pi };
	}

	auto normalize_distribution = [&](std::unordered_map<std::string, double>& dist) {
		double sum = 0.0;
		for (const auto& kv : dist) {
			sum += kv.second;
		}
		if (!(sum > 0.0)) {
			double uniform = 1.0 / static_cast<double>(taxid_list.size());
			for (auto& kv : dist) {
				kv.second = uniform;
			}
			return;
		}
		double inv_sum = 1.0 / sum;
		for (auto& kv : dist) {
			kv.second *= inv_sum;
		}
	};

	// 使用加权证据初始化，让共享 reads 也能提供初始贡献
	std::unordered_map<std::string, double> weighted_evidence;
	weighted_evidence.reserve(taxid_list.size());
	for (const auto& record : results) {
		if (record.taxidCount.empty()) continue;

		double total_count = 0.0;
		for (const auto& kv : record.taxidCount) {
			if (!detail::is_unclassified(kv.first)) {
				total_count += kv.second;
			}
		}

		if (total_count <= 0.0) continue;

		// Soft assignment：按比例分配该 read 的权重
		for (const auto& kv : record.taxidCount) {
			if (detail::is_unclassified(kv.first)) continue;
			weighted_evidence[kv.first] +=
				kv.second / total_count;
		}
	}

	auto get_scale = [&](const std::string& taxid) -> double {
		if (!prior_scale) return 1.0;
		auto it = prior_scale->find(taxid);
		if (it == prior_scale->end() || !(it->second > 0.0)) return 1.0;
		return std::sqrt(it->second);
	};

	std::unordered_map<std::string, double> pi;
	pi.reserve(taxid_list.size());
	double background_prior = 1e-9;
	for (const auto& taxid : taxid_list) {
		double u = 0.0;
		auto it_w = weighted_evidence.find(taxid);
		if (it_w != weighted_evidence.end()) {
			u = it_w->second;
		}

		double s = get_scale(taxid);
		double init_val = u * s;

		// 若没有 evidence 但 prior_scale>0，则给极小值兜底
		if (init_val <= 0.0 && prior_scale && s > 0.0) {
			init_val = 1e-5;
		}
		if (init_val <= 0.0) {
			init_val = background_prior;
		}

		pi[taxid] = init_val;
	}

	normalize_distribution(pi);

	std::unordered_map<std::string, double> abundance_prior = pi;

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
				double second_count = 0.0;
				std::string best_taxid;
				for (const auto& [taxid, count] : source.taxidCount) {
					if (detail::is_unclassified(taxid)) {
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
					auto it = pi.find(taxid);
					if (it == pi.end()) {
						continue;
					}
					double likelihood = (c + options.eps) / denom;
					double log_likelihood = options.temp * std::log(std::max(likelihood, options.eps));
					double prior = std::max(it->second, options.eps);
					double log_prior = std::log(prior);
					double coexist = 0.0;
					if (!best_taxid.empty() && taxid != best_taxid) {
						coexist = penalty_scale;
					}
					log_components.emplace_back(taxid, log_prior + log_likelihood - coexist);
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
				}

				// confidence-weighted expected counts
				double max_post = 0.0;
				for (const auto& p : posterior) {
					if (p.second > max_post) max_post = p.second;
				}
				double conf_w = 1.0;
				if (options.conf_power > 0.0 && max_post > 0.0) {
					conf_w = std::pow(max_post, options.conf_power);
				}
					double base_w = (source.sample_weight > 0.0)
					                     ? source.sample_weight
					                     : ((source.evaluated > 0.0) ? source.evaluated : 1.0);
					double weight = base_w * conf_w;

				for (const auto& p : posterior) {
					local_expected[p.first] += weight * p.second;
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

		// Sparsity: 相对剪枝，抑制低丰度噪声
		double max_expected = 0.0;
		for (const auto& kv : expected_counts) {
			if (kv.second > max_expected) {
				max_expected = kv.second;
			}
		}
		double prune_thres = max_expected * options.prune_ratio; // 低于主导物种一定比例的视为噪声

		double diff = 0.0;
		for (const auto& taxid : taxid_list) {
			double expected = 0.0;
			auto it = expected_counts.find(taxid);
			if (it != expected_counts.end()) {
				expected = it->second;
			}
			if (expected < prune_thres) {
				expected = 0.0;
			}
			double prior_component = (prior_mass > 0.0)
			                             ? prior_mass * abundance_prior[taxid]
			                             : options.alpha;
			double updated = (expected + prior_component) / denominator;
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
