// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm>  // for sort
#include <cassert>    // for assert
#include <cinttypes>  // for uint64_t, int64_t
#include <cmath>      // for floor, sqrt
#include <cstddef>    // for size_t
#include <functional> // for equal_to
#include <limits>     // for numeric_limits
#include <omp.h>      // for omp_get_thread_num
#include <random>     // for uniform_int_distribution, mt19937_64
#include <utility>    // for swap, move
#include <vector>     // for vector

#include <robin_hood.h> // for hash, unordered_flat_map
#include <HyperLogLog.hpp> // for HyperLogLog
#include <toolbox.hpp>     // for clustering_node, entry, neighbor, prio_queue, distance_matrix, cluste...

namespace toolbox
{
	void sort_by_cardinalities(std::vector<size_t> const& kmer_counts, std::vector<size_t>& positions)
	{
		assert(positions.size() <= kmer_counts.size());

		auto cardinality_compare = [&](size_t const index1, size_t const index2)
			{
				return kmer_counts[index1] > kmer_counts[index2];
			};

		std::sort(positions.begin(), positions.end(), cardinality_compare);
	}

	void precompute_union_estimates_for(std::vector<uint64_t>& estimates,
		std::vector<HyperLogLog> const& sketches,
		std::vector<size_t> const& counts,
		std::vector<size_t> const& positions,
		int64_t const j)
	{
		assert(counts.size() == sketches.size());
		assert(positions.size() <= counts.size());
		assert(estimates.size() == sketches.size()); // Resize happens in precompute_init_interval_union_estimations
		assert(estimates.size() > static_cast<size_t>(j));

		HyperLogLog temp_hll = sketches[positions[j]];
		estimates[j] = counts[positions[j]];

		for (int64_t j_prime = j - 1; j_prime >= 0; --j_prime)
		{
			estimates[j_prime] = static_cast<uint64_t>(temp_hll.merge_and_estimate(sketches[positions[j_prime]]));
		}
	}

	void precompute_initial_union_estimates(std::vector<uint64_t>& estimates,
		std::vector<HyperLogLog> const& sketches,
		std::vector<size_t> const& counts,
		std::vector<size_t> const& positions)
	{
		assert(counts.size() == sketches.size());
		assert(positions.size() <= counts.size());
		assert(sketches.size() > 0u);

		estimates.resize(sketches.size());

		HyperLogLog temp_hll = sketches[positions[0]];
		estimates[0] = counts[positions[0]];

		for (size_t j = 1; j < positions.size(); ++j)
			estimates[j] = static_cast<uint64_t>(temp_hll.merge_and_estimate(sketches[positions[j]]));
	}

#if 0 // Currently unused
	uint64_t estimate_interval(std::vector<HyperLogLog> const& sketches, std::vector<size_t> const& positions)
	{
		assert(positions.size() <= sketches.size());
		assert(!positions.empty());

		HyperLogLog temp_hll = sketches[positions[0]];

		for (size_t i = 1; i < positions.size(); ++i)
			temp_hll.merge(sketches[positions[i]]);

		return temp_hll.estimate();
	}
#endif

	void rearrange_bins(std::vector<HyperLogLog> const& sketches,
		std::vector<size_t> const& kmer_counts,
		std::vector<size_t>& positions,
		double const max_ratio,
		size_t const num_threads)
	{
		std::vector<size_t> permutation;

		size_t first = 0;
		size_t last = 1;

		while (first < positions.size())
		{
			// size difference is too large or sequence is over -> do the clustering
			if (last == positions.size() || kmer_counts[positions[first]] * max_ratio > kmer_counts[positions[last]])
			{
				// if this is not the first group, we want one bin overlap
				cluster_bins(sketches, positions, permutation, first, last, num_threads);
				first = last;
			}
			++last;
		}

		for (size_t i{ 0 }; i < permutation.size(); ++i)
		{
			size_t swap_index = permutation[i];
			while (swap_index < i)
				swap_index = permutation[swap_index];

			std::swap(positions[i], positions[swap_index]);
		}
	}

	void cluster_bins(std::vector<HyperLogLog> const& sketches,
		std::vector<size_t>& positions,
		std::vector<size_t>& permutation,
		size_t const first,
		size_t const last,
		size_t const num_threads)
	{
		assert(num_threads >= 1);
		assert(positions.size() <= sketches.size());
		assert((first == 0) == permutation.empty());

		size_t const n = sketches.size();
		size_t const chunk_size = std::floor(std::sqrt(n));

		size_t const prune_steps = chunk_size;
		size_t steps_without_prune = 0;

		size_t const none = std::numeric_limits<size_t>::max();
		/* internal map that stores the distances
			*
			* The first layer is a hash map with the ids of active clusters as keys.
			* The values (second layer) are priority queues with neighbors of the cluster
			* with the respective key in the first layer.
			* These neighbors are themselves clusters with an id and store a distance to the
			* cluster of the first layer.
			*/
		distance_matrix dist;
		dist.reserve(n + 1);

		// map that indicates which ids of clusters are still in the distance matrix
		// the values are the indices where the priority queue for the given id as key can be found in dist
		robin_hood::unordered_flat_map<size_t, size_t> remaining_ids;

		// clustering tree stored implicitly in a vector
		std::vector<clustering_node> clustering;
		clustering.reserve(2 * n);

		// cache for hll cardinality estimates
		std::vector<double> estimates;
		estimates.reserve(2 * n);

		// every thread will write its observed id with minimal distance to some other here
		// id == none means that the thread observed only empty or no priority queues
		std::vector<size_t> min_ids(num_threads, none);

		// these will be the new ids for new clusters
		// the first one is invalid, but it will be incremented before it is used for the first time
		size_t new_id = last - 1;

		// initialize clustering and estimates
		for (size_t id = first; id < last; ++id)
		{
			// id i is at the index i - first
			clustering.push_back({ none, none, sketches[positions[id]] });
			estimates.emplace_back(sketches[positions[id]].estimate());
		}

		// if this is not the first group, we want to have one overlapping bin
		size_t previous_rightmost = none;
		if (first != 0)
		{
			// For all other clusters, their id is also their index in filesnames, sketches etc. .
			// This is important, because their id is then inserted into the clustering.
			// This does not work for previous rightmost, because its index does not necessarily lie on
			// the continuous spectrum from first to last. We run into a problem, because the entries are
			// stored in vectors. Therefore we give previous_rightmost a different id (==last). This is
			// fine, because we only need the HLL sketch of the actual index. previous_rightmost will be ignored
			// in the traceback anyway and won't be added to the permutation in this step.
			size_t const actual_previous_rightmost = permutation.back();
			++new_id;
			previous_rightmost = new_id;

			clustering.push_back({ none, none, sketches[positions[actual_previous_rightmost]] });
			estimates.emplace_back(sketches[positions[actual_previous_rightmost]].estimate());
		}

		// initialize priority queues in the distance matrix (sequentially)
		for (size_t id = first; id < first + clustering.size(); ++id)
		{
			// empty priority queue for every item in clustering
			dist.push_back({ id, prio_queue{} });
			remaining_ids[id] = id - first;
		}

#pragma omp parallel num_threads(num_threads)
		{
			double min_dist = std::numeric_limits<double>::max();
			// minimum distance exclusively for this thread

			// initialize all the priority queues of the distance matrix
			// while doing that, compute the first min_id
#pragma omp for schedule(nonmonotonic : dynamic, chunk_size)
			for (size_t i = 0; i < clustering.size(); ++i)
			{
				for (size_t j = 0; j < clustering.size(); ++j)
				{
					// we only want one diagonal of the distance matrix
					if (i < j)
					{
						// this must be a copy, because merging changes the hll sketch
						HyperLogLog temp_hll = clustering[i].hll;
						double const estimate_ij = temp_hll.merge_and_estimate(clustering[j].hll);
						// Jaccard distance estimate
						double const distance = 2 - (estimates[i] + estimates[j]) / estimate_ij;
						dist[i].pq.push({ j + first, distance });
					}
				}
				if (dist[i].pq.empty())
					continue;

				// check if the just initialized priority queue contains the minimum value for this thread
				neighbor const& curr = dist[i].pq.top();
				if (curr.dist < min_dist)
				{
					min_dist = curr.dist;
					min_ids[omp_get_thread_num()] = dist[i].id;
				}
			} // implicit barrier

	// a single thread shuffles dist to approximately balance loads in static scheduling
#pragma omp single
			random_shuffle(dist, remaining_ids);

			// main loop of the clustering
			// keep merging nodes until we have a complete tree
			while (remaining_ids.size() > 1)
			{
				// Wait for all threads to have evaluated remaining_ids.size() as remaining_ids
				// may be modified by the following pragma omp single.
#pragma omp barrier

#pragma omp single
				{
					// perform critical update
					// increment id for the new cluster (must be done at the beginning)
					++new_id;

					// compute the final min_id from the min_ids of the worker threads
					size_t min_id = min_ids[0];
					double min_dist = std::numeric_limits<double>::max();
					for (auto candidate_id : min_ids)
					{
						// check if the thread saw any id
						if (candidate_id == none)
							continue;

						size_t const dist_index = remaining_ids.at(candidate_id);
						neighbor const& curr = dist[dist_index].pq.top();
						if (curr.dist < min_dist)
						{
							min_dist = curr.dist;
							min_id = candidate_id;
						}
					}

					size_t const min_index = remaining_ids.at(min_id); // how can min_id be none?
					size_t const neighbor_id = dist[min_index].pq.top().id;

					// merge the two nodes with minimal distance together insert the new node into the clustering
					clustering.push_back({ min_id, neighbor_id, std::move(clustering[min_id - first].hll) });
					estimates.emplace_back(clustering.back().hll.merge_and_estimate(clustering[neighbor_id - first].hll));

					// remove old ids
					remaining_ids.erase(min_id);
					remaining_ids.erase(neighbor_id);

					// overwrite one of the removed entries with the new one
					remaining_ids[new_id] = min_index;
					dist[min_index] = { new_id, prio_queue{} };

					// prune the distance matrix to reduce overhead due to inactive entries
					++steps_without_prune;
					if (steps_without_prune > prune_steps)
					{
						prune(dist, remaining_ids);
						steps_without_prune = 0;
					}
				} // implicit barrier

				// reset values for the computation of the new minimum
				min_ids[omp_get_thread_num()] = none;
				min_dist = std::numeric_limits<double>::max();

				HyperLogLog const new_hll = clustering.back().hll;

				// update distances in dist
				// while doing that, compute the new min_id
#pragma omp for schedule(static)
				for (size_t i = 0; i < dist.size(); ++i)
				{
					size_t other_id = dist[i].id;
					if (other_id == new_id || !remaining_ids.contains(other_id))
						continue;

					// this must be a copy, because merge_and_estimate() changes the hll
					HyperLogLog temp_hll = new_hll;
					double const estimate_ij = temp_hll.merge_and_estimate(clustering[other_id - first].hll);
					// Jaccard distance estimate
					double const distance = 2 - (estimates[other_id - first] + estimates.back()) / estimate_ij;
					dist[i].pq.push({ new_id, distance });

					// make sure the closest neighbor is not yet deleted (this is a lazy update)
					while (!remaining_ids.contains(dist[i].pq.top().id))
					{
						dist[i].pq.pop();
					}

					// check if the just updated priority queue contains the minimum value for this thread
					neighbor const& curr = dist[i].pq.top();
					if (curr.dist < min_dist)
					{
						min_dist = curr.dist;
						min_ids[omp_get_thread_num()] = other_id;
					}
				} // implicit barrier
			}
		} // end of the parallel region

		size_t final_root_index = remaining_ids.begin()->second;
		size_t final_root_id = dist[final_root_index].id;

		// rotate the previous rightmost to the left so that it has the correct place in the permutation
		if (first != 0)
		{
			rotate(clustering, previous_rightmost, first, final_root_id);
		}

		// traceback into permutation and ignore the previous rightmost
		trace(clustering, permutation, previous_rightmost, first, final_root_id);
	}

	void random_shuffle(distance_matrix& dist, robin_hood::unordered_flat_map<size_t, size_t>& remaining_ids)
	{
		size_t const n = dist.size();

		std::mt19937_64 gen(0x7E1E5665D46800E5ULL);

		for (size_t i = 0; i < n - 1; ++i)
		{
			std::uniform_int_distribution<size_t> distrib(i, n - 1);
			size_t const swap_i = distrib(gen);

			size_t const id = dist[i].id;
			size_t const swap_id = dist[swap_i].id;

			// swap entries and update the reming ids, because the indices in dist changed
			std::swap(dist[i], dist[swap_i]);
			std::swap(remaining_ids[id], remaining_ids[swap_id]);
		}
	}

	void prune(distance_matrix& dist, robin_hood::unordered_flat_map<size_t, size_t>& remaining_ids)
	{
		if (dist.empty())
			return; // GCOVR_EXCL_LINE

		// index of the first entry after the valid range
		size_t valid_range_end = 0;
		// index of the first entry before the invalid range
		size_t invalid_range_start = dist.size() - 1;

		while (valid_range_end != invalid_range_start)
		{
			if (remaining_ids.contains(dist[valid_range_end].id))
			{
				++valid_range_end;
			}
			else if (!remaining_ids.contains(dist[invalid_range_start].id))
			{
				--invalid_range_start;
			}
			else
			{
				// If we arrive here, then valid_range_end has an invalid id
				// and invalid_range_start has a valid id. The correspoding entries should be swapped
				std::swap(dist[valid_range_end], dist[invalid_range_start]);

				// update the index of the valid entry
				remaining_ids.at(dist[valid_range_end].id) = valid_range_end;
			}
		}

		// check the last element between the valid and invalid range
		if (remaining_ids.contains(dist[valid_range_end].id))
		{
			++valid_range_end;
		}

		// cut off invalid values
		dist.resize(valid_range_end);
	}

	bool rotate(std::vector<clustering_node>& clustering,
		size_t const previous_rightmost,
		size_t const first,
		size_t const id)
	{
		if (id == previous_rightmost) // we are at the leaf that is previous_rightmost (Anchor of the recursion)
			return true;

		clustering_node& curr = clustering[id - first];

		if (curr.left == std::numeric_limits<size_t>::max()) // we are at a leaf that is not previous_rightmost
		{
			return false;
		}
		// nothing to do if previous_rightmost is in the left subtree
		else if (rotate(clustering, previous_rightmost, first, curr.left))
		{
			return true;
		}
		// rotate if previous_rightmost is in the right subtree
		else if (rotate(clustering, previous_rightmost, first, curr.right))
		{
			std::swap(curr.left, curr.right);
			return true;
		}

		// else: previous_rightmost is not in this subtree
		return false;
	}

	void trace(std::vector<clustering_node> const& clustering,
		std::vector<size_t>& permutation,
		size_t const previous_rightmost,
		size_t const first,
		size_t const id)
	{
		clustering_node const& curr = clustering[id - first];

		if (curr.left == std::numeric_limits<size_t>::max()) // I am at a leaf
		{
			if (id != previous_rightmost)
				permutation.push_back(id);
			return;
		}

		trace(clustering, permutation, previous_rightmost, first, curr.left);
		trace(clustering, permutation, previous_rightmost, first, curr.right);
	}
} // namespace seqan::hibf::sketch::toolbox