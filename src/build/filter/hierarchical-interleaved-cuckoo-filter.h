/*
 * -----------------------------------------------------------------------------
 * Filename:      hierarchical-interleaved-cuckoo-filter.h
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-10-27
 *
 * Last Modified: 2024-10-
 *
 * Description:
 *  This file defines the Hierarchical Interleaved Cuckoo Filter
 *
 * Version:
 *  1.0
 * -----------------------------------------------------------------------------
 */

#pragma once

#include <HyperLogLog.hpp>
#include <buildConfig.hpp>
#include <optional>
#include <robin_hood.h>
#include <atomic>
#include <kvec.h>

namespace chimera::hicf {
	struct Layout
	{
		struct maxBin
		{
			std::vector<size_t> previousTechnical{};
			size_t id{};
			friend auto operator<=>(maxBin const&, maxBin const&) = default;
		};
		struct userBin
		{
			std::vector<size_t> previousTechnical{};
			size_t technicalBinId{};
			size_t technicalBinNum{};
			size_t index{};
			friend auto operator<=>(userBin const&, userBin const&) = default;
		};
		size_t topMaxBinId{};
		std::vector<maxBin> maxBins{};
		std::vector<userBin> userBins{};
		bool operator==(Layout const& other) const = default;
	};

	struct DataStore
	{
		struct previousLevel
		{
			std::vector<size_t> binIndices{};
			std::string binsNum;
			bool empty() const
			{
				return binIndices.empty();
			};
		};
		Layout* layout{};
		std::vector<size_t>* estimateKmerCounts{};
		std::vector<HyperLogLog>* hllVec{};

		previousLevel previous{};
		std::vector<uint64_t> estimates;
		long long totalTime = 0;
		std::vector<size_t> location;
		double loadFactor;
		double relaxedLoadFactor;

		void createLocation()
		{
			location.resize(estimateKmerCounts->size());
			std::iota(location.begin(), location.end(), 0);
		}

	};

	struct graph {
		struct node {
			std::vector<node> children{};
			size_t parentBinId{};
			size_t maxBinId{};
			size_t technicalBinNum{};
			std::optional<size_t> favouriteChildId{ std::nullopt };
			std::vector<Layout::userBin> remainingUserBins{};

			bool maxBinIsMerged() const
			{
				return favouriteChildId.has_value();
			}
		};
		node root;

		void updateHeaderNodeData(std::vector<Layout::maxBin>& headerMaxBins, graph& icfGraph)
		{
			for (auto const& [previousTechnical, id] : headerMaxBins)
			{
				node* parent = &icfGraph.root;
				auto binId = std::ranges::prev(previousTechnical.end());
				for (auto it = previousTechnical.begin(); it != binId; ++it)
				{
					for (node& child : parent->children)
					{
						if (child.parentBinId == *it)
						{
							parent = &child;
							break;
						}
					}
				}

				node newChild{ {}, previousTechnical.back(), id, 0u, std::nullopt, {} };
				parent->children.push_back(newChild);
				if (parent->maxBinId == previousTechnical.back())
				{
					parent->favouriteChildId = parent->children.size() - 1;
				}

			}
		}

		void updateContentNodeData(std::vector<Layout::userBin> const& layoutUserBins, graph& icfGraph)
		{
			for (size_t userBinIndex = 0; userBinIndex < layoutUserBins.size(); userBinIndex++)
			{
				auto const& record = layoutUserBins[userBinIndex];
				node* currentNode = &icfGraph.root;
				for (size_t technicalBinIndex = 0; technicalBinIndex < record.previousTechnical.size(); technicalBinIndex++)
				{
					size_t const bin = record.previousTechnical[technicalBinIndex];
					size_t const technicalBinNum = 1;
					currentNode->technicalBinNum = std::max(bin + technicalBinNum, currentNode->technicalBinNum);
					for (node& child : currentNode->children)
					{
						if (child.parentBinId == bin)
						{
							currentNode = &child;
							break;
						}
					}
				}
				size_t const bin = record.technicalBinId;
				size_t const technicalBinNum = record.technicalBinNum;
				currentNode->technicalBinNum = std::max(bin + technicalBinNum, currentNode->technicalBinNum);
				if (record.technicalBinId == currentNode->maxBinId)
				{
					currentNode->remainingUserBins.insert(currentNode->remainingUserBins.begin(), record);
				}
				else
				{
					currentNode->remainingUserBins.emplace_back(record);
				}
			}
		}


		graph(Layout& layout)
		{
			root.children = {};
			root.parentBinId = 0;
			root.maxBinId = layout.topMaxBinId;
			root.technicalBinNum = 0;
			root.favouriteChildId = std::nullopt;
			root.remainingUserBins = {};

			updateHeaderNodeData(layout.maxBins, *this);
			updateContentNodeData(layout.userBins, *this);
		}

		graph() = default;
		

	};

	struct buildData {
		std::atomic<size_t> icfNumber{};
		ChimeraBuild::HICFConfig const& config;
		graph icfGraph{};

		size_t requestICFIndex()
		{
			return std::atomic_fetch_add(&icfNumber, 1u);
		}
	};

	class HierarchicalInterleavedCuckooFilter {
		typedef kvec_t(int) kvector;
		size_t userBinsNum{};
		size_t technicalBinsNum{};
		uint8_t tagNum{ 4 };
		uint8_t tagBits{ 16 };
		uint8_t hyperBits{ 12 };
		int MaxCuckooCount{ 500 };
		double loadFactor{ 0.95 };
		double relaxedLoadFactor{ 0.95 };
		// 存储了多个ICF；
		std::vector<chimera::InterleavedCuckooFilter> icfVec;
		// 存储了每个ICF的每个BIN对应的下一级ICF的ID；
		std::vector<std::vector<size_t>> nextICFId;
		// 存储了每个ICF的每个BIN对应的用户bin的ID；
		std::vector<std::vector<size_t>> icfBinToUserBinId{};

		friend std::ostream& operator<<(std::ostream& os, const HierarchicalInterleavedCuckooFilter& hicf) {
			os << "Hierarchical Interleaved Cuckoo Filter" << std::endl;
			os << "User bins: " << hicf.userBinsNum << std::endl;
			os << "Technical bins: " << hicf.technicalBinsNum << std::endl;
			os << "Tag bits: " << (int)hicf.tagBits << std::endl;
			os << "Hyper bits: " << (int)hicf.hyperBits << std::endl;
			os << "Max cuckoo count: " << hicf.MaxCuckooCount << std::endl;
			os << "Load factor: " << hicf.loadFactor << std::endl;
			os << "Relaxed load factor: " << hicf.relaxedLoadFactor << std::endl;
			return os;

		}

		void insertICF(robin_hood::unordered_flat_set<uint64_t> const& kmers,
			size_t const binNum,
			size_t const binIndex,
			chimera::InterleavedCuckooFilter& icf)
		{
			size_t const kmersPerBin = std::ceil(kmers.size() / binNum);
			size_t chunkNum{};
			for (auto chunk : kmers | seqan::std::views::chunk(kmersPerBin))
			{
				size_t Index = binIndex + chunkNum;
				chunkNum++;
				for (auto kmer : chunk)
				{
					icf.insertTag(Index, kmer);
				}
			}
		}

		size_t hierarchicalBuild(	HierarchicalInterleavedCuckooFilter& hicf,
									robin_hood::unordered_flat_set<uint64_t>& parentKmers,
									graph::node const& currentNode,
									buildData& data,
									bool isRoot,
									std::vector<std::string>& indexToTaxid)
		{
			size_t const icfIndex = data.requestICFIndex();
			auto& technicalBinToICFID = hicf.nextICFId[icfIndex];
			technicalBinToICFID.resize(currentNode.technicalBinNum, icfIndex);
			auto& technicalBinToUserBinID = hicf.icfBinToUserBinId[icfIndex];
			// use std::numeric_limits<size_t>::max() to indicate that the bin is merged
			technicalBinToUserBinID.resize(currentNode.technicalBinNum, std::numeric_limits<size_t>::max());

			auto& icf = hicf.icfVec[icfIndex];
			robin_hood::unordered_flat_set<uint64_t> kmers{};
			
			auto initialiseMaxBinKers = [&]() -> size_t {
				if (currentNode.maxBinIsMerged())
				{
					technicalBinToICFID[currentNode.maxBinId] = hierarchicalBuild(hicf, kmers, currentNode.children[currentNode.favouriteChildId.value()], data, false, indexToTaxid);
					return 1;
				}
				else
				{
					auto const& record = currentNode.remainingUserBins[0];
					std::string miniFilePath = "tmp/" + indexToTaxid[record.index] + ".mini";
					std::ifstream ifile(miniFilePath, std::ios::binary);
					if (ifile.fail()) {
						std::cerr << "Failed to open minimiser file: " << indexToTaxid[record.index] << ".mini" << std::endl;
					}
					uint64_t hash;
					while (ifile.read(reinterpret_cast<char*>(&hash), sizeof(hash))) {
						kmers.insert(hash);
					}

					std::fill_n(technicalBinToUserBinID.begin() + record.technicalBinId, record.technicalBinNum, record.index);
					ifile.close();

					return record.technicalBinNum;
				}
			};


			// Construct the ICF
			size_t const maxBinKmers = initialiseMaxBinKers();
			bool const maxBinIsMerged = currentNode.maxBinIsMerged();
			size_t const kmersPerBin = std::ceil(kmers.size() / maxBinKmers);
			size_t const binSize = (maxBinIsMerged) ? kmersPerBin / relaxedLoadFactor : kmersPerBin / loadFactor;
			size_t const binNum = currentNode.technicalBinNum;
			icf = chimera::InterleavedCuckooFilter{ binNum, binSize, 16 };
			//std::cout << icf << std::endl;
			insertICF(kmers, maxBinKmers, currentNode.maxBinId, icf);
			if (!isRoot)
			{
				parentKmers.insert(kmers.begin(), kmers.end());
			}
			kmers.clear();


			std::mt19937 rng(std::random_device{}());
			auto loopOverChildren = [&]() {
				if (currentNode.children.empty())
				{
					return;
				}
				std::vector<graph::node> children = currentNode.children;
				std::vector<size_t> childBinIndices(children.size());
				std::iota(childBinIndices.begin(), childBinIndices.end(), size_t{});

				if (currentNode.maxBinIsMerged())
				{
					std::erase(childBinIndices, currentNode.favouriteChildId.value());
				}

				if (isRoot)
				{
					std::shuffle(childBinIndices.begin(), childBinIndices.end(), rng);
				}
				size_t const mutexNum = currentNode.technicalBinNum;
				std::vector<std::mutex> mutexes(mutexNum);

#pragma omp parallel for schedule(dynamic)
				for (auto&& index : childBinIndices)
				{
					auto& child = children[index];
					robin_hood::unordered_flat_set<uint64_t> kmers{};
					size_t const icfIndex = hierarchicalBuild(hicf, kmers, child, data, false, indexToTaxid);
					auto parentBinIndex = child.parentBinId;
					{
						size_t const mutexID{ parentBinIndex };
						std::lock_guard<std::mutex> lock(mutexes[mutexID]);
						technicalBinToICFID[parentBinIndex] = icfIndex;
						insertICF(kmers, 1, parentBinIndex, icf);
						if (!isRoot)
						{
							parentKmers.insert(kmers.begin(), kmers.end());
						}
					}
				}
			};

			loopOverChildren();

			size_t const start = currentNode.maxBinIsMerged() ? 0u : 1u;
			for (size_t i = start; i < currentNode.remainingUserBins.size(); i++)
			{
				auto const& record = currentNode.remainingUserBins[i];
				if (isRoot && record.technicalBinNum == 1)
				{
					std::string miniFilePath = "tmp/" + indexToTaxid[record.index] + ".mini";
					std::ifstream ifile(miniFilePath, std::ios::binary);
					if (ifile.fail()) {
						std::cerr << "Failed to open minimiser file: " << indexToTaxid[record.index] << ".mini" << std::endl;
					}
					uint64_t hash;
					while (ifile.read(reinterpret_cast<char*>(&hash), sizeof(hash))) {
						icf.insertTag(record.technicalBinId, hash);
					}
					ifile.close();
				}
				else
				{
					std::string miniFilePath = "tmp/" + indexToTaxid[record.index] + ".mini";
					std::ifstream ifile(miniFilePath, std::ios::binary);
					if (ifile.fail()) {
						std::cerr << "Failed to open minimiser file: " << indexToTaxid[record.index] << ".mini" << std::endl;
					}
					uint64_t hash;
					while (ifile.read(reinterpret_cast<char*>(&hash), sizeof(hash))) {
						kmers.insert(hash);
					}
					insertICF(kmers, record.technicalBinNum, record.technicalBinId, icf);
					ifile.close();

					if (!isRoot)
					{
						parentKmers.insert(kmers.begin(), kmers.end());
					}
				}
				std::fill_n(technicalBinToUserBinID.begin() + record.technicalBinId, record.technicalBinNum, record.index);
				kmers.clear();
			}

			return icfIndex;
		}

	public:

		HierarchicalInterleavedCuckooFilter() = default;
		HierarchicalInterleavedCuckooFilter(Layout& layout, ChimeraBuild::HICFConfig& config, std::vector<std::string>& indexToTaxid)
		{
			userBinsNum = config.userBinsNum;
			technicalBinsNum = config.technicalBinsNum;
			size_t icfNum = layout.maxBins.size() + 1;
			icfVec.resize(icfNum);
			icfBinToUserBinId.resize(icfNum);
			nextICFId.resize(icfNum);
			loadFactor = config.loadFactor;
			relaxedLoadFactor = config.relaxedLoadFactor;

			graph icfGraph(layout);
			graph::node const& rootNode = icfGraph.root;
			size_t const maxTechnicalBinNum = rootNode.technicalBinNum;
			robin_hood::unordered_flat_set<uint64_t> rootKmers{};
			buildData data{ 0, config, icfGraph };

			hierarchicalBuild(*this, rootKmers, rootNode, data, true, indexToTaxid);
		}

		template <class Archive>
		void serialize(Archive& archive) {
			archive(userBinsNum, technicalBinsNum, tagNum, tagBits, hyperBits, MaxCuckooCount, loadFactor, icfVec, nextICFId, icfBinToUserBinId);
		}


		typedef kvec_t(bool) kvectorBool;
		template <std::ranges::range value_range_t>
		inline void bulkContain(value_range_t&& values, kvector& result, size_t const icfIndex, size_t const threshold)
		{
			const auto& bulkResult = icfVec[icfIndex].bulk_count(values);
			size_t sum{};
			for (size_t bin{}; bin < bulkResult.n; bin++)
			{
				sum += kv_A(bulkResult, bin);
				auto const userBinId = icfBinToUserBinId[icfIndex][bin];
				if (userBinId == std::numeric_limits<size_t>::max())
				{
					if (sum >= threshold)
					{
						bulkContain(values, result, nextICFId[icfIndex][bin], threshold);
					}
					sum = 0u;
				}
				else if (bin + 1u == bulkResult.n || userBinId != icfBinToUserBinId[icfIndex][bin + 1])
				{
					if (sum >= threshold)
					{
						kv_A(result, userBinId) = sum;
					}
					sum = 0u;
				}
			}
		}

		template <std::ranges::range value_range_t>
		inline kvector bulkCount(value_range_t&& values,
									size_t const threshold)
		{
			kvector result;
			kv_init(result);
			kv_resize(int, result, userBinsNum);
			std::memset(result.a, 0, sizeof(int) * userBinsNum);
			kv_size(result) = userBinsNum;
			
			bulkContain(values, result, 0, threshold);

			return result;
		}
	};
}