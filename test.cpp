//#include <seqan3/alphabet/nucleotide/dna4.hpp>
//#include <seqan3/core/debug_stream.hpp>
//#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
//#include <seqan3/search/views/kmer_hash.hpp>
//#include <seqan3/search/views/minimiser_hash.hpp>
//#include <iostream>
//#include <vector>
//#include <algorithm>
//#include <chrono>
//#include <sdsl/int_vector.hpp>
//#include <numeric> // for std::iota
//#include <kvec.h>
//#include <interleaved-cuckoo-filter.h>
//#include <random>
//#include <cereal/archives/binary.hpp>
//#include <cereal/types/memory.hpp>
//#include <CLI11.hpp>
//using namespace seqan3;
//using namespace seqan3::literals;
//// 随机生成DNA序列
//std::vector<dna4> generate_random_dna_sequence(size_t length, std::mt19937& gen) {
//	std::vector<dna4> sequence;
//	sequence.reserve(length);
//	std::uniform_int_distribution<> dis(0, 3);
//
//	for (size_t i = 0; i < length; ++i) {
//		sequence.push_back(dna4{}.assign_rank(dis(gen)));
//	}
//
//	return sequence;
//}
//
//std::unordered_map<size_t, std::array<size_t, 3>> calculate_combined_expected_counts(auto&& sequence1, auto&& sequence2, auto&& sequence3, auto&& hash_adaptor) {
//	std::unordered_map<size_t, std::array<size_t, 3>> combined_counts;
//
//	// 统计每个序列的哈希值计数
//	for (auto&& value : sequence1 | hash_adaptor) {
//		combined_counts[value][0]++;
//	}
//	for (auto&& value : sequence2 | hash_adaptor) {
//		combined_counts[value][1]++;
//	}
//	for (auto&& value : sequence3 | hash_adaptor) {
//		combined_counts[value][2]++;
//	}
//
//	return combined_counts;
//}
//
//// 计算预期结果
//std::vector<std::vector<size_t>> calculate_expected_counts(auto&& sequence1, auto&& sequence2, auto&& sequence3, auto&& hash_adaptor) {
//	std::unordered_map<size_t, std::array<size_t, 3>> combined_counts;
//
//	// 统计每个序列的哈希值计数
//	for (auto&& value : sequence1 | hash_adaptor) {
//		combined_counts[value][0]++;
//	}
//	for (auto&& value : sequence2 | hash_adaptor) {
//		combined_counts[value][1]++;
//	}
//	for (auto&& value : sequence3 | hash_adaptor) {
//		combined_counts[value][2]++;
//	}
//
//	// 初始化结果
//	std::vector<size_t> result1(3, 0); // 序列1的结果
//	std::vector<size_t> result2(3, 0); // 序列2的结果
//	std::vector<size_t> result3(3, 0); // 序列3的结果
//
//	// 填充预期结果
//	for (const auto& [hash_value, count_array] : combined_counts) {
//		// 检查每个序列的哈希值
//		if (count_array[0] > 0) { // 序列1的所有hash都存入bin 0
//			result1[0] += count_array[0];
//			if (count_array[1] > 0) result1[1] += count_array[0]; // 序列1的hash在序列2出现的次数
//			if (count_array[2] > 0) result1[2] += count_array[0]; // 序列1的hash在序列3出现的次数
//		}
//
//		if (count_array[1] > 0) { // 序列2的所有hash都存入bin 4
//			result2[1] += count_array[1];
//			if (count_array[0] > 0) result2[0] += count_array[1]; // 序列2的hash在序列1出现的次数
//			if (count_array[2] > 0) result2[2] += count_array[1]; // 序列2的hash在序列3出现的次数
//		}
//
//		if (count_array[2] > 0) { // 序列3的所有hash都存入bin 7
//			result3[2] += count_array[2];
//			if (count_array[0] > 0) result3[0] += count_array[2]; // 序列3的hash在序列1出现的次数
//			if (count_array[1] > 0) result3[1] += count_array[2]; // 序列3的hash在序列2出现的次数
//		}
//	}
//
//	return { result1, result2, result3 };
//}
//
//// 输出计数结果
//void print_results(const std::vector<size_t>& result, const std::string& label) {
//	std::cout << label << ": ";
//	for (const auto& count : result) {
//		std::cout << count << " ";
//	}
//	std::cout << std::endl;
//}
//
//int main() {
//	// 声明选项变量，不能是 constexpr
//	size_t sequence_length = 5000; // 每个序列的长度
//	constexpr size_t num_bins = 5;         // Bloom 和 Cuckoo 过滤器的 bin 数量
//
//	auto hash_adaptor = views::minimiser_hash(shape{ ungapped{4} }, window_size{ 8 });
//	std::mt19937 gen(19);
//	// 生成随机 DNA 序列
//	auto sequence1 = generate_random_dna_sequence(sequence_length, gen);
//	auto sequence2 = generate_random_dna_sequence(sequence_length, gen);
//	auto sequence3 = generate_random_dna_sequence(sequence_length, gen);
//	// 计算预期结果
//	auto expected_counts = calculate_expected_counts(sequence1, sequence2, sequence3, hash_adaptor);
//
//	// 输出预期计数结果
//	print_results(expected_counts[0], "预期结果1");
//	print_results(expected_counts[1], "预期结果2");
//	print_results(expected_counts[2], "预期结果3");
//
//	// IBF 的构建和测试
//	{
//		// 记录构建时间
//		auto start = std::chrono::high_resolution_clock::now();
//		interleaved_bloom_filter ibf{ bin_count{num_bins},
//									 bin_size{8192u},
//									 hash_function_count{2u} };
//		auto end = std::chrono::high_resolution_clock::now();
//		std::chrono::duration<double> elapsed_seconds = end - start;
//		std::cout << "IBF的构建时间：" << elapsed_seconds.count() << "秒\n";
//
//		// 记录插入时间
//		start = std::chrono::high_resolution_clock::now();
//		for (auto&& value : sequence1 | hash_adaptor)
//			ibf.emplace(value, bin_index{ 0u });
//		for (auto&& value : sequence2 | hash_adaptor)
//			ibf.emplace(value, bin_index{ 1u });
//		for (auto&& value : sequence3 | hash_adaptor)
//			ibf.emplace(value, bin_index{ 2u });
//		end = std::chrono::high_resolution_clock::now();
//		elapsed_seconds = end - start;
//		std::cout << "IBF的插入时间：" << elapsed_seconds.count() << "秒\n";
//
//		// 记录查询时间
//		start = std::chrono::high_resolution_clock::now();
//		auto agent = ibf.counting_agent();
//		auto ibf_result1 = agent.bulk_count(sequence1 | hash_adaptor);
//		auto ibf_result2 = agent.bulk_count(sequence2 | hash_adaptor);
//		auto ibf_result3 = agent.bulk_count(sequence3 | hash_adaptor);
//		end = std::chrono::high_resolution_clock::now();
//		elapsed_seconds = end - start;
//		std::cout << "IBF的查询时间：" << elapsed_seconds.count() << "秒\n";
//
//		// 输出查询结果
//		debug_stream << "IBF查询结果1: " << ibf_result1 << '\n';
//		debug_stream << "IBF查询结果2: " << ibf_result2 << '\n';
//		debug_stream << "IBF查询结果3: " << ibf_result3 << '\n';
//
//		std::ofstream os("ibf.cereal", std::ios::binary);
//		cereal::BinaryOutputArchive archive(os);
//		archive(ibf);
//	}
//
//	// ICF 的构建和测试
//	{
//		// 记录构建时间
//		auto start = std::chrono::high_resolution_clock::now();
//		chimera::InterleavedCuckooFilter icf(num_bins, 8192);
//		auto end = std::chrono::high_resolution_clock::now();
//		std::chrono::duration<double> elapsed_seconds = end - start;
//		std::cout << "ICF的构建时间：" << elapsed_seconds.count() << "秒\n";
//
//		// 记录插入时间
//		start = std::chrono::high_resolution_clock::now();
//		for (auto&& value : sequence1 | hash_adaptor)
//			icf.insertTag(0, value);
//		for (auto&& value : sequence2 | hash_adaptor)
//			icf.insertTag(1, value);
//		for (auto&& value : sequence3 | hash_adaptor)
//			icf.insertTag(2, value);
//		end = std::chrono::high_resolution_clock::now();
//		elapsed_seconds = end - start;
//		std::cout << "ICF的插入时间：" << elapsed_seconds.count() << "秒\n";
//
//		// 记录查询时间
//		start = std::chrono::high_resolution_clock::now();
//		auto icf_result1 = icf.bulk_count(sequence1 | hash_adaptor);
//		auto icf_result2 = icf.bulk_count(sequence2 | hash_adaptor);
//		auto icf_result3 = icf.bulk_count(sequence3 | hash_adaptor);
//		end = std::chrono::high_resolution_clock::now();
//		elapsed_seconds = end - start;
//		std::cout << "ICF的查询时间：" << elapsed_seconds.count() << "秒\n";
//
//		// 输出查询结果
//		std::cout << "ICF查询结果1: ";
//		for (size_t i = 0; i < kv_size(icf_result1); ++i) {
//			std::cout << kv_A(icf_result1, i) << " ";
//		}
//		std::cout << std::endl;
//
//		std::cout << "ICF查询结果2: ";
//		for (size_t i = 0; i < kv_size(icf_result2); ++i) {
//			std::cout << kv_A(icf_result2, i) << " ";
//		}
//		std::cout << std::endl;
//
//		std::cout << "ICF查询结果3: ";
//		for (size_t i = 0; i < kv_size(icf_result3); ++i) {
//			std::cout << kv_A(icf_result3, i) << " ";
//		}
//		std::cout << std::endl;
//
//		// 序列化过程
//		{
//			std::ofstream os("icf.cereal", std::ios::binary);
//			if (!os) {
//				std::cerr << "Error opening file for writing: icf.cereal" << std::endl;
//				return -1;
//			}
//			cereal::BinaryOutputArchive outputArchive(os); // 改为 outputArchive
//			outputArchive(icf);
//		}
//
//		// 反序列化过程
//		chimera::InterleavedCuckooFilter loaded_icf(num_bins, 8192);
//		{
//			std::ifstream is("icf.cereal", std::ios::binary);
//			if (!is) {
//				std::cerr << "Error opening file for reading: icf.cereal" << std::endl;
//				return -1;
//			}
//			cereal::BinaryInputArchive inputArchive(is); // 改为 inputArchive
//			inputArchive(loaded_icf);
//		}
//
//		// 使用loaded_icf再查询一遍
//		auto loaded_icf_result1 = loaded_icf.bulk_count(sequence1 | hash_adaptor);
//		auto loaded_icf_result2 = loaded_icf.bulk_count(sequence2 | hash_adaptor);
//		auto loaded_icf_result3 = loaded_icf.bulk_count(sequence3 | hash_adaptor);
//		// 输出查询结果
//		std::cout << "Loaded ICF查询结果1: ";
//		for (size_t i = 0; i < kv_size(loaded_icf_result1); ++i) {
//			std::cout << kv_A(loaded_icf_result1, i) << " ";
//		}
//		std::cout << std::endl;
//		std::cout << "Loaded ICF查询结果2: ";
//		for (size_t i = 0; i < kv_size(loaded_icf_result2); ++i) {
//			std::cout << kv_A(loaded_icf_result2, i) << " ";
//		}
//		std::cout << std::endl;
//		std::cout << "Loaded ICF查询结果3: ";
//		for (size_t i = 0; i < kv_size(loaded_icf_result3); ++i) {
//			std::cout << kv_A(loaded_icf_result3, i) << " ";
//		}
//		std::cout << std::endl;
//	}
//
//	return 0;
//}
//#include <khashOperation.hpp>
//
//int main() {
//	// 初始化哈希表
//	khash_t(str_vec_str)* h = kh_init(str_vec_str);
//
//	// 插入元素
//	khash_insert_vec(h, "fruits", std::vector<std::string>{ "apple", "banana", "cherry" });
//	khash_insert_vec(h, "colors", std::vector<std::string>{ "red", "green", "blue" });
//
//	// 检索元素
//	std::vector<std::string> values;
//	if (khash_find_vec(h, "fruits", values)) {
//		std::cout << "Found 'fruits' with values: [";
//		for (const auto& val : values) {
//			std::cout << val << " ";
//		}
//		std::cout << "]" << std::endl;
//	}
//	else {
//		std::cout << "'fruits' not found" << std::endl;
//	}
//
//	khash_insert_vec(h, "fruits", std::vector<std::string>{ "orange", "pear" });
//	khash_insert_vec(h, "fruits", std::string("peach"));
//	khash_insert_vec(h, "people", std::string("man"));
//
//	// 遍历所有元素
//	khash_iterate_vec(h);
//
//	// 删除元素
//	khash_delete_vec(h, "fruits");
//
//	// 再次遍历以验证删除效果
//	khash_iterate_vec(h);
//
//	// 删除所有元素以释放内存
//	for (khiter_t k = kh_begin(h); k != kh_end(h); ++k) {
//		if (kh_exist(h, k)) {
//			free((char*)kh_key(h, k)); // 释放每个键的内存
//		}
//	}
//
//	// 销毁哈希表
//	kh_destroy(str_vec_str, h);
//
//	return 0;
//}