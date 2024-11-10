// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cinttypes>  // for uint64_t, int64_t
#include <cstddef>    // for size_t
#include <functional> // for greater
#include <queue>      // for priority_queue
#include <vector>     // for vector

#include <robin_hood.h> // for unordered_flat_map
#include <HyperLogLog.hpp> // for HyperLogLog

namespace toolbox
{

/*!\brief type for a node in the clustering tree when for the rearrangement
 * \ingroup hibf_sketch_toolbox
 */
struct clustering_node
{
    // children in the tree
    size_t left;
    size_t right;
    // hll sketch of the union if the node is still a root
    HyperLogLog hll;
};

/*!\brief element of the second priority queue layer of the distance matrix
 * \ingroup hibf_sketch_toolbox
 */
struct neighbor
{
    size_t id;
    double dist;

    bool operator>(neighbor const & other) const
    {
        return dist > other.dist;
    }
};

/*!\brief type of a min heap based priority queue
 * \ingroup hibf_sketch_toolbox
 */
using prio_queue = std::priority_queue<neighbor, std::vector<neighbor>, std::greater<neighbor>>;

/*!\brief entry of the distance matrix that has the id of a cluster with its neighbors in a prio queue
 * \ingroup hibf_sketch_toolbox
 */
struct entry
{
    size_t id;
    prio_queue pq;
};

/*!\brief type of the distance matrix for the clustering for the rearrangement
 * \ingroup hibf_sketch_toolbox
 */
using distance_matrix = std::vector<entry>;

/*!\brief Sorts filenames and cardinalities by looking only at the cardinalities.
 * \ingroup hibf_sketch_toolbox
 */
void sort_by_cardinalities(std::vector<size_t> const & counts, std::vector<size_t> & positions);

/*!\brief Estimate the cardinality of the union for a single user bin j with all prior ones j' < j.
 * \ingroup hibf_sketch_toolbox
 * \param[out] estimates output row
 * \param[in] sketches The HyperLogLog sketches of the respective user bins.
 * \param[in] counts The counts/sketch.estimates() of the respective user bins.
 * \param[in] positions The realtive positions of the input information to correctly access sketches and counts.
 * \param[in] j The current user bin (column in the DP matrix)
 *
 * estimates[j_prime] will be the union cardinality estimate of the interval {j_prime, ..., j}.
 */
void precompute_union_estimates_for(std::vector<uint64_t> & estimates,
                                    std::vector<HyperLogLog> const & sketches,
                                    std::vector<size_t> const & counts,
                                    std::vector<size_t> const & positions,
                                    int64_t const j);

/*!\brief Estimate the cardinality of the union for each interval [0, j] for all user bins j.
 * \ingroup hibf_sketch_toolbox
 * \param[out] estimates output row
 * \param[in] sketches The HyperLogLog sketches of the respective user bins.
 * \param[in] counts The counts/sketch.estimates() of the respective user bins.
 * \param[in] positions The realtive positions of the input information to correctly access sketches and counts.
 *
 * estimates[j] will be the union cardinality estimate of the interval {0, ..., j}.
 */
void precompute_initial_union_estimates(std::vector<uint64_t> & estimates,
                                        std::vector<HyperLogLog> const & sketches,
                                        std::vector<size_t> const & counts,
                                        std::vector<size_t> const & positions);

#if 0 // Currently unused
/*!\brief Estimate the cardinality of the union for a single interval.
 * \ingroup hibf_sketch_toolbox
 * \param[in] sketches The HyperLogLog sketches to be used for estimation.
 * \param[in] positions The realtive positions of the input information to correctly access sketches and counts.
 * \returns The the cardinality of the union for the interval [start, end).
 */
uint64_t estimate_interval(std::vector<HyperLogLog> const & sketches, std::vector<size_t> const & positions);
#endif

/*!\brief Rearrange filenames, sketches and counts such that similar bins are close to each other
 * \ingroup hibf_sketch_toolbox
 * \param[in] sketches The HyperLogLog sketches of the respective user bins.
 * \param[in] counts The counts/sketch.estimates() of the respective user bins.
 * \param[in,out] positions The realtive positions of the input information to correctly access sketches and counts.
 * \param[in] max_ratio the maximal cardinality ratio in the clustering intervals (must be <= 1 and >= 0)
 * \param[in] num_threads the number of threads to use
 */
void rearrange_bins(std::vector<HyperLogLog> const & sketches,
                    std::vector<size_t> const & counts,
                    std::vector<size_t> & positions,
                    double const max_ratio,
                    size_t const num_threads);

/*!\brief Perform an agglomerative clustering variant on the index range [first:last)
 * \ingroup hibf_sketch_toolbox
 * \param[in] sketches The HyperLogLog sketches of the respective user bins.
 * \param[in,out] positions The realtive positions of the input information to correctly access sketches and counts.
 * \param[in,out] permutation append the new order to this
 * \param[in] first id of the first cluster of the interval
 * \param[in] last id of the last cluster of the interval plus one
 * \param[in] num_threads the number of threads to use
 */

void cluster_bins(std::vector<HyperLogLog> const & sketches,
                  std::vector<size_t> & positions,
                  std::vector<size_t> & permutation,
                  size_t const first,
                  size_t const last,
                  size_t const num_threads);

/*!\brief Randomly swap entries in dist while keeping track of the changes of indices.
 * \ingroup hibf_sketch_toolbox
 * \param[in] dist the distance matrix (vector of priority queues) to shuffle
 * \param[in] remaining_ids the map with information about which ids remain at which index
 */
void random_shuffle(distance_matrix & dist, robin_hood::unordered_flat_map<size_t, size_t> & remaining_ids);

/*!\brief Delete inactive entries out of dist and shrink to fit its size while keeping track of the changes of indices
 * \ingroup hibf_sketch_toolbox
 * \param[in] dist the distance matrix (vector of priority queues) to prune
 * \param[in] remaining_ids the map with information about which ids remain at which index
 */
void prune(distance_matrix & dist, robin_hood::unordered_flat_map<size_t, size_t> & remaining_ids);

/*!\brief Rotate the previous rightmost bin to the left of the clustering tree
 * \ingroup hibf_sketch_toolbox
 * \param[in, out] clustering the tree to do the rotation on
 * \param[in] previous_rightmost the id of the node to be rotated to the left
 * \param[in] first the id of the first node in the interval to shift the index
 * \param[in] id the id of the current node
 *
 * If called with the root of the tree, this function recursively calls itself while rotating
 * several subtrees until previous_rightmost is at the very left end of the whole clustering tree.
 *
 * \return whether previous rightmost was in the subtree rooted at id
 */
bool rotate(std::vector<clustering_node> & clustering,
            size_t const previous_rightmost,
            size_t const first,
            size_t const id);

/*!\brief Do a recursive traceback to find the order of leaves in the clustering tree
 * \ingroup hibf_sketch_toolbox
 * \param[in] clustering the tree to do the traceback on
 * \param[out] permutation append the new order to this
 * \param[in] previous_rightmost the id of the node on the left which should be ignored
 * \param[in] first the id of the first node in the interval to shift the index
 * \param[in] id the id of the current node
 *
 * This function traverses the tree in depth-first-search accessing the leaves from left to right.
 * 'Left to right' refers to the order of nodes in `clustering`.
 */
void trace(std::vector<clustering_node> const & clustering,
           std::vector<size_t> & permutation,
           size_t const previous_rightmost,
           size_t const first,
           size_t const id);

} // namespace seqan::hibf::sketch::toolbox
