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
 * Last Modified: 2024-09-19
 *
 * Description:
 * Expectation-Maximization algorithm for Chimera
 *
 * Version:
 *  1.0
 * -----------------------------------------------------------------------------
 */
#include <iostream>
#include <vector>
#include <robin_hood.h>
#include <string>
#include <cmath>
#include <limits>
#include <classifyConfig.hpp>
using namespace ChimeraClassify;
/*
* @brief Get the top match taxid
* 
* This function returns the top match taxid based on the probability.
* 
* @param matches A vector of taxid-count pairs
* @param prob A map of taxid-probability pairs
*/
inline std::string getTopMatch(const std::vector<std::pair<std::string, size_t>>& matches,
	const robin_hood::unordered_map<std::string, double>& prob) {
	std::string top_taxid = matches[0].first;
	size_t top_kcount = matches[0].second;
	double max_prob = 0.0;

	for (const auto& m : matches) {
		const std::string& taxid = m.first;
		size_t kcount = m.second;
		if (prob.at(taxid) > max_prob) {
			max_prob = prob.at(taxid);
			top_taxid = taxid;
			top_kcount = kcount;
		}
	}
	return top_taxid;
}

/*
* @brief Initialize the probabilities
* 
* This function initializes the probabilities for each taxid.
* 
* @param classifyResults A vector of classifyResult objects
* @param prob A map of taxid-probability pairs
*/
inline void initializeProbabilities(const std::vector<classifyResult>& classifyResults,
	robin_hood::unordered_map<std::string, double>& prob) {
	for (const auto& result : classifyResults) {
		for (const auto& pair : result.taxidCount) {
			const std::string& taxid = pair.first;
			if (prob.find(taxid) == prob.end()) {
				prob[taxid] = 0.0;  
			}
		}
	}
}

/*
* @brief Expectation-Maximization algorithm for Chimera
* 
* This function implements the Expectation-Maximization algorithm for Chimera.
* 
* @param classifyResults A vector of classifyResult objects
* @param max_iter Maximum number of iterations
* @param threshold Convergence threshold
*/
inline std::vector<classifyResult> EMAlgorithm(std::vector<classifyResult>& classifyResults,
	size_t max_iter = 100, double threshold = 0.001) {
	robin_hood::unordered_map<std::string, double> prob;  // Store the probability of each taxid
	robin_hood::unordered_map<std::string, size_t> initial_weight;  // Initial weights
	size_t total_weight = classifyResults.size();  // Total number of sequences

	initializeProbabilities(classifyResults, prob);

	// Initialize the initial weights for each taxid
	for (const auto& result : classifyResults) {
		if (result.taxidCount.size() == 1) {
			const std::string& taxid = result.taxidCount[0].first;
			initial_weight[taxid] += 1;
		}
	}

	// Calculate the initial probabilities
	for (const auto& entry : initial_weight) {
		prob[entry.first] = static_cast<double>(entry.second) / total_weight;
	}

	size_t em_iter_cnt = 0;
	while (em_iter_cnt < max_iter) {
		robin_hood::unordered_map<std::string, size_t> reassigned_matches = initial_weight;
		double diff = 0.0;

		// Expectation step: Reassign sequences to the most likely taxid
		for (const auto& result : classifyResults) {
			if (result.taxidCount.size() > 1) {
				std::string top_taxid = getTopMatch(result.taxidCount, prob);
				reassigned_matches[top_taxid] += 1;
			}
		}

		// Maximization step: Calculate new probabilities
		for (auto& entry : reassigned_matches) {
			const std::string& taxid = entry.first;
			size_t count = entry.second;
			double new_prob = static_cast<double>(count) / total_weight;
			diff += std::abs(prob[taxid] - new_prob);
			prob[taxid] = new_prob;
		}

		// Check convergence condition
		if (diff <= threshold) {
			break;
		}

		em_iter_cnt++;
	}

	std::vector<classifyResult> updatedClassifyResults;
	for (const auto& result : classifyResults) {
		std::string top_taxid = getTopMatch(result.taxidCount, prob);
		updatedClassifyResults.push_back({ result.id, {{top_taxid, 1}} });  // Update the taxidCount
	}

	return updatedClassifyResults;
}