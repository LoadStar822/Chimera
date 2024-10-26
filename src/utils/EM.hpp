/*
 * -----------------------------------------------------------------------------
 * Filename:      EM.cpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-09-19
 *
 * Last Modified: 2024-10-20
 *
 * Description:
 * Expectation-Maximization algorithm for Chimera (Parallelized Version with OpenMP)
 *
 * Version:
 *  1.2
 * -----------------------------------------------------------------------------
 */
#include <iostream>
#include <vector>
#include <robin_hood.h>
#include <string>
#include <cmath>
#include <limits>
#include <classifyConfig.hpp>
#include <omp.h>  
using namespace ChimeraClassify;

/*
* @brief Get the top match taxid based on probabilities
*
* @param matches A vector of taxid-count pairs
* @param prob A map of taxid-probability pairs
*/
inline std::string getTopMatch(const std::vector<std::pair<std::string, size_t>>& matches,
	const robin_hood::unordered_map<std::string, double>& prob) {
	std::string top_taxid = matches[0].first;
	double max_prob = 0.0;

	for (const auto& m : matches) {
		const std::string& taxid = m.first;
		auto it = prob.find(taxid);
		if (it != prob.end() && it->second > max_prob) {
			max_prob = it->second;
			top_taxid = taxid;
		}
	}
	return top_taxid;
}

/*
* @brief Initialize the probabilities uniformly for all taxids
*
* @param classifyResults A vector of classifyResult objects
* @param prob A map of taxid-probability pairs
*/
inline void initializeProbabilities(const std::vector<classifyResult>& classifyResults,
	robin_hood::unordered_map<std::string, double>& prob) {
	robin_hood::unordered_set<std::string> taxid_set;
	for (const auto& result : classifyResults) {
		for (const auto& pair : result.taxidCount) {
			taxid_set.insert(pair.first);
		}
	}
	double initial_prob = 1.0 / taxid_set.size();
	for (const auto& taxid : taxid_set) {
		prob[taxid] = initial_prob;
	}
}

/*
* @brief Expectation-Maximization algorithm for Chimera with OpenMP parallelization
*
* @param classifyResults A vector of classifyResult objects
* @param max_iter Maximum number of iterations
* @param threshold Convergence threshold
*/
inline std::vector<classifyResult> EMAlgorithm(std::vector<classifyResult>& classifyResults,
	size_t max_iter = 100, double threshold = 0.001) {
	robin_hood::unordered_map<std::string, double> prob;  // Store the probability of each taxid
	size_t total_weight = classifyResults.size();          // Total number of sequences

	initializeProbabilities(classifyResults, prob);

	size_t em_iter_cnt = 0;
	while (em_iter_cnt < max_iter) {
		robin_hood::unordered_map<std::string, double> expected_counts;

		int num_threads = omp_get_max_threads();
		std::vector<robin_hood::unordered_map<std::string, double>> thread_expected_counts(num_threads);

#pragma omp parallel
		{
			int thread_id = omp_get_thread_num();
			auto& local_counts = thread_expected_counts[thread_id];

#pragma omp for nowait
			for (size_t i = 0; i < classifyResults.size(); ++i) {
				const auto& result = classifyResults[i];
				double total_prob = 0.0;
				for (const auto& pair : result.taxidCount) {
					total_prob += prob[pair.first];
				}
				for (const auto& pair : result.taxidCount) {
					const std::string& taxid = pair.first;
					double weight = prob[taxid] / total_prob;
					local_counts[taxid] += weight;
				}
			}
		}

		for (const auto& local_counts : thread_expected_counts) {
			for (const auto& entry : local_counts) {
				expected_counts[entry.first] += entry.second;
			}
		}

		// Maximization step: Update probabilities
		double sum_counts = 0.0;
		for (const auto& entry : expected_counts) {
			sum_counts += entry.second;
		}

		double diff = 0.0;

		std::vector<std::string> taxid_list;
		for (const auto& entry : prob) {
			taxid_list.push_back(entry.first);
		}

#pragma omp parallel for reduction(+:diff)
		for (size_t i = 0; i < taxid_list.size(); ++i) {
			const std::string& taxid = taxid_list[i];
			double new_prob = 0.0;
			auto it = expected_counts.find(taxid);
			if (it != expected_counts.end()) {
				new_prob = it->second / sum_counts;
			}
			diff += std::abs(prob[taxid] - new_prob);
			prob[taxid] = new_prob;
		}

		// Check convergence condition
		if (diff <= threshold) {
			break;
		}

		em_iter_cnt++;
	}

	std::vector<classifyResult> updatedClassifyResults(classifyResults.size());

#pragma omp parallel for
	for (size_t i = 0; i < classifyResults.size(); ++i) {
		const auto& result = classifyResults[i];
		std::string top_taxid = getTopMatch(result.taxidCount, prob);
		updatedClassifyResults[i] = { result.id, {{top_taxid, 1}} };
	}

	return updatedClassifyResults;
}
