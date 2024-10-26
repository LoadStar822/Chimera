/*
 * -----------------------------------------------------------------------------
 * Filename:      VEM.hpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-10-22
 *
 * Last Modified: 2024-10-22
 *
 * Description:
 * Variational Expectation-Maximization algorithm for Chimera
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
#include <omp.h>  
using namespace ChimeraClassify;

/*
* @brief Get the top match taxid based on variational probabilities
*
* @param variational_params A vector of taxid-probability pairs for a sample
*/
inline std::string getTopMatch(const robin_hood::unordered_map<std::string, double>& variational_params) {
	std::string top_taxid;
	double max_prob = -1.0;

	for (const auto& m : variational_params) {
		const std::string& taxid = m.first;
		double prob = m.second;
		if (prob > max_prob) {
			max_prob = prob;
			top_taxid = taxid;
		}
	}
	return top_taxid;
}

/*
* @brief Initialize the model probabilities and variational parameters uniformly
*
* @param classifyResults A vector of classifyResult objects
* @param model_probs A map of taxid-probability pairs (model parameters)
* @param variational_params A vector of maps for each sample, storing taxid-probability pairs (variational parameters)
*/
inline void initializeVariationalParameters(const std::vector<classifyResult>& classifyResults,
	robin_hood::unordered_map<std::string, double>& model_probs,
	std::vector<robin_hood::unordered_map<std::string, double>>& variational_params) {

	robin_hood::unordered_set<std::string> taxid_set;
	for (const auto& result : classifyResults) {
		for (const auto& pair : result.taxidCount) {
			taxid_set.insert(pair.first);
		}
	}

	double initial_model_prob = 1.0 / taxid_set.size();
	for (const auto& taxid : taxid_set) {
		model_probs[taxid] = initial_model_prob;
	}

	variational_params.resize(classifyResults.size());
#pragma omp parallel for
	for (size_t i = 0; i < classifyResults.size(); ++i) {
		const auto& result = classifyResults[i];
		robin_hood::unordered_map<std::string, double> local_params;
		double initial_var_prob = 1.0 / result.taxidCount.size();
		for (const auto& pair : result.taxidCount) {
			local_params[pair.first] = initial_var_prob;
		}
		variational_params[i] = std::move(local_params);
	}
}

/*
* @brief Variational Expectation-Maximization algorithm for Chimera
*
* @param classifyResults A vector of classifyResult objects
* @param max_iter Maximum number of iterations
* @param threshold Convergence threshold
*/
inline std::vector<classifyResult> VEMAlgorithm(std::vector<classifyResult>& classifyResults,
	size_t max_iter = 100, double threshold = 0.001) {

	robin_hood::unordered_map<std::string, double> model_probs;  // Model parameters
	std::vector<robin_hood::unordered_map<std::string, double>> variational_params;  // Variational parameters

	initializeVariationalParameters(classifyResults, model_probs, variational_params);

	double prev_elbo = -std::numeric_limits<double>::infinity();
	size_t em_iter_cnt = 0;

	while (em_iter_cnt < max_iter) {
		// Variational E-step: Update variational parameters
		double elbo = 0.0;

#pragma omp parallel for reduction(+:elbo)
		for (size_t i = 0; i < classifyResults.size(); ++i) {
			const auto& result = classifyResults[i];
			auto& q = variational_params[i];
			double sum_q = 0.0;

			// Update q(taxid) ¡Ø exp(E[log P(taxid | data)])
			for (const auto& pair : result.taxidCount) {
				const std::string& taxid = pair.first;
				double log_prob = std::log(model_probs[taxid]);
				q[taxid] = std::exp(log_prob);
				sum_q += q[taxid];
			}

			// Normalize q(taxid)
			for (auto& pair : q) {
				pair.second /= sum_q;
			}

			// Compute ELBO contribution for this sample
			for (const auto& pair : q) {
				double q_val = pair.second;
				if (q_val > 0) {
					elbo += q_val * (std::log(model_probs[pair.first]) - std::log(q_val));
				}
			}
		}

		// M-step: Update model parameters
		robin_hood::unordered_map<std::string, double> expected_counts;
		double sum_counts = 0.0;

#pragma omp parallel
		{
			robin_hood::unordered_map<std::string, double> local_counts;

#pragma omp for nowait
			for (size_t i = 0; i < variational_params.size(); ++i) {
				const auto& q = variational_params[i];
				for (const auto& pair : q) {
#pragma omp atomic
					local_counts[pair.first] += pair.second;
				}
			}

#pragma omp critical
			{
				for (const auto& pair : local_counts) {
					expected_counts[pair.first] += pair.second;
				}
			}
		}

		for (const auto& pair : expected_counts) {
			sum_counts += pair.second;
		}

		double diff = 0.0;

#pragma omp parallel for reduction(+:diff)
		for (size_t i = 0; i < expected_counts.size(); ++i) {
			auto it = expected_counts.begin();
			std::advance(it, i);
			const std::string& taxid = it->first;
			double new_prob = it->second / sum_counts;
			diff += std::abs(model_probs[taxid] - new_prob);
			model_probs[taxid] = new_prob;
		}

		// Check convergence based on ELBO
		if (std::abs(elbo - prev_elbo) <= threshold) {
			break;
		}
		prev_elbo = elbo;

		em_iter_cnt++;
	}

	// Update classification results
	std::vector<classifyResult> updatedClassifyResults(classifyResults.size());

#pragma omp parallel for
	for (size_t i = 0; i < classifyResults.size(); ++i) {
		const auto& result = classifyResults[i];
		const auto& q = variational_params[i];
		std::string top_taxid = getTopMatch(q);
		updatedClassifyResults[i] = { result.id, {{top_taxid, 1}} };
	}

	return updatedClassifyResults;
}
