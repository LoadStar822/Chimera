//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <vector>
//#include <zlib.h>
//#include "kseq.h"
//#include <seqan3/alphabet/nucleotide/dna4.hpp>
//#include <seqan3/search/views/minimiser_hash.hpp>
//#include "cuckoofilter.h"
//#include <cereal/types/vector.hpp>
//#include <cereal/archives/binary.hpp>
//#include <cereal/types/memory.hpp>
//
//#include <sdsl/bit_vectors.hpp>
//
//#include <seqan3/core/debug_stream.hpp>
//#include <seqan3/search/views/kmer_hash.hpp>
//
//#include <kvec.h>
//#include <interleaved-cuckoo-filter.h>
//#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
//using namespace seqan3::literals;
//KSEQ_INIT(gzFile, gzread)
//
//std::vector<std::pair<std::string, uint32_t>> read_target_info(const std::string& filepath) {
//	std::ifstream file(filepath);
//	std::vector<std::pair<std::string, uint32_t>> targets;
//	std::string line;
//
//	while (std::getline(file, line)) {
//		std::istringstream iss(line);
//		std::string path;
//		uint32_t taxid;
//		if (!(iss >> path >> taxid)) { break; }
//		targets.emplace_back(path, taxid);
//	}
//
//	return targets;
//}
//
//std::vector<seqan3::dna4> read_fna_gz(const std::string& filepath) {
//	using namespace seqan3::literals;
//	gzFile fp = gzopen(filepath.c_str(), "r");
//	kseq_t* seq = kseq_init(fp);
//	std::vector<seqan3::dna4> sequence;
//
//	while (kseq_read(seq) >= 0) {
//		for (size_t i = 0; i < seq->seq.l; ++i) {
//			switch (seq->seq.s[i]) {
//			case 'A': sequence.push_back('A'_dna4); break;
//			case 'C': sequence.push_back('C'_dna4); break;
//			case 'G': sequence.push_back('G'_dna4); break;
//			case 'T': sequence.push_back('T'_dna4); break;
//			default: break;
//			}
//		}
//	}
//
//	kseq_destroy(seq);
//	gzclose(fp);
//
//	return sequence;
//}
//
//int main() {
//	// 记录运行时间
//	auto start = std::chrono::high_resolution_clock::now();
//
//	std::string target_info_file = "/mnt/d/code/src/ganon/test_files/build/target_info.tsv";
//	auto targets = read_target_info(target_info_file);
//
//	//size_t total_items = 100000;
//	//cuckoofilter::CuckooFilter<size_t, 12> filter(total_items);
//
//	seqan3::interleaved_bloom_filter ibf{ seqan3::bin_count{ 20u },
//									  seqan3::bin_size{ 3800681u },
//									  seqan3::hash_function_count{ 2u } };
//	interleaved_cuckoo_filter::InterleavedCuckooFilter icf(20, 3800681);
//	std::cout << "Size of ibf: " << sizeof(ibf) << " bytes" << std::endl;
//	std::cout << "Size of icf: " << sizeof(icf) << " bytes" << std::endl;
//	auto agent = ibf.counting_agent();
//	int count = 0;
//	for (const auto& [path, taxid] : targets) {
//		auto sequence = read_fna_gz(path);
//		auto minimizer_hash = sequence | seqan3::views::minimiser_hash(seqan3::shape{ seqan3::ungapped{4} }, seqan3::window_size{ 8 });
//		for (auto hash : minimizer_hash) {
//			ibf.emplace(hash, seqan3::bin_index{ count });
//			icf.insertTag(count, hash);
//		}
//		//for (auto hash : minimizer_hash) {
//		//	if (!icf.lookupTag(count, hash)) {
//		//		std::cout << "Hash does not exist in the icf.\n";
//		//	}
//		//}
//		//seqan3::debug_stream << agent.bulk_count(minimizer_hash) << '\n';
//		count++;
//	}
//
//	////分别序列化icf和ibf
//	//{
//	//	std::ofstream os("icf.cereal", std::ios::binary);
//	//	cereal::BinaryOutputArchive archive(os);
//	//	archive(icf);
//	//}
//	//{
//	//	std::ofstream os("ibf.cereal", std::ios::binary);
//	//	cereal::BinaryOutputArchive archive(os);
//	//	archive(ibf);
//	//}
//
//	//// 序列化 Cuckoo Filter
//	//{
//	//	std::ofstream os("cuckoofilter.cereal", std::ios::binary);
//	//	cereal::BinaryOutputArchive archive(os);
//	//	archive(filter);
//	//}
//
//	//// 反序列化 Cuckoo Filter
//	//cuckoofilter::CuckooFilter<size_t, 12> loaded_filter(total_items); // 提供必要的参数
//	//{
//	//	std::ifstream is("cuckoofilter.cereal", std::ios::binary);
//	//	cereal::BinaryInputArchive archive(is);
//	//	archive(loaded_filter);
//	//}
//
//	//// 检查第一个 hash 是否存在于反序列化后的过滤器中
//	//if (!targets.empty()) {
//	//	auto first_sequence = read_fna_gz(targets[0].first);
//	//	auto first_minimizer_hash = first_sequence | seqan3::views::minimiser_hash(seqan3::shape{ seqan3::ungapped{4} }, seqan3::window_size{ 8 });
//	//	auto first_hash = *first_minimizer_hash.begin();
//	//	if (loaded_filter.Contain(first_hash) == cuckoofilter::Ok) {
//	//		std::cout << "Hash exists in the loaded filter.\n";
//	//	}
//	//	else {
//	//		std::cout << "Hash does not exist in the loaded filter.\n";
//	//	}
//	//}
//	return 0;
//}

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <sdsl/int_vector.hpp>
#include <numeric> // for std::iota
#include <kvec.h>
#include <interleaved-cuckoo-filter.h>
using namespace seqan3::literals;

//int main()
//{
//	//记录构建时间
//	auto start = std::chrono::high_resolution_clock::now();
//	seqan3::interleaved_bloom_filter ibf{ seqan3::bin_count{ 8u },
//									  seqan3::bin_size{ 8192u },
//									  seqan3::hash_function_count{ 2u } };
//	auto end = std::chrono::high_resolution_clock::now();
//	std::chrono::duration<double> elapsed_seconds = end - start;
//	std::cout << "IBF的构建时间：" << elapsed_seconds.count() << "s\n";
//
//	auto const sequence1 = "ACTGACTGACTGATC"_dna4;
//	auto const sequence2 = "GTGACTGACTGACTCG"_dna4;
//	auto const sequence3 = "AAAAAAACGATCGACA"_dna4;
//	auto       hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{ seqan3::ungapped{4} }, seqan3::window_size{ 8 });
//
//	std::cout << "bin0: " << std::endl;
//	for (auto&& value : sequence1 | hash_adaptor)
//		seqan3::debug_stream << value << '\n';
//	std::cout << "bin4: " << std::endl;
//	for (auto&& value : sequence2 | hash_adaptor)
//		seqan3::debug_stream << value << '\n';
//	std::cout << "bin7: " << std::endl;
//	for (auto&& value : sequence3 | hash_adaptor)
//		seqan3::debug_stream << value << '\n';
//
//	//记录插入的时间
//	start = std::chrono::high_resolution_clock::now();
//
//	// Insert all 5-mers of sequence1 into bin 0
//	for (auto&& value : sequence1 | hash_adaptor)
//		ibf.emplace(value, seqan3::bin_index{ 0u });
//
//	// Insert all 5-mers of sequence2 into bin 4
//	for (auto&& value : sequence2 | hash_adaptor)
//		ibf.emplace(value, seqan3::bin_index{ 4u });
//
//	// Insert all 5-mers of sequence3 into bin 7
//	for (auto&& value : sequence3 | hash_adaptor)
//		ibf.emplace(value, seqan3::bin_index{ 7u });
//
//	//记录插入的时间
//	end = std::chrono::high_resolution_clock::now();
//	elapsed_seconds = end - start;
//	std::cout << "IBF的插入时间：" << elapsed_seconds.count() << "s\n";
//
//	//记录查询时间
//	start = std::chrono::high_resolution_clock::now();
//	auto agent = ibf.counting_agent();
//
//	// Count all 5-mers of sequence1 for all bins
//	seqan3::debug_stream << agent.bulk_count(sequence1 | hash_adaptor) << '\n'; // [11,0,0,0,9,0,0,0]
//	seqan3::debug_stream << agent.bulk_count(sequence2 | hash_adaptor) << '\n'; // [0,0,0,0,0,0,0,0]
//	seqan3::debug_stream << agent.bulk_count(sequence3 | hash_adaptor) << '\n'; // [0,0,0,0,0,0,0,0]
//
//	//记录查询时间
//	end = std::chrono::high_resolution_clock::now();
//	elapsed_seconds = end - start;
//	std::cout << "IBF的查询时间：" << elapsed_seconds.count() << "s\n";
//
//	//记录构建时间
//	start = std::chrono::high_resolution_clock::now();
//	chimera::InterleavedCuckooFilter icf(8, 8192);
//	end = std::chrono::high_resolution_clock::now();
//	elapsed_seconds = end - start;
//	std::cout << "ICF的构建时间：" << elapsed_seconds.count() << "s\n";
//
//	//记录插入的时间
//	start = std::chrono::high_resolution_clock::now();
//	for (auto&& value : sequence1 | hash_adaptor)
//		icf.insertTag(0, value);
//	for (auto&& value : sequence2 | hash_adaptor)
//		icf.insertTag(4, value);
//	for (auto&& value : sequence3 | hash_adaptor)
//		icf.insertTag(7, value);
//	end = std::chrono::high_resolution_clock::now();
//	elapsed_seconds = end - start;
//	std::cout << "ICF的插入时间：" << elapsed_seconds.count() << "s\n";
//
//	//记录查询时间
//	start = std::chrono::high_resolution_clock::now();
//	auto result = icf.bulk_count(sequence1 | hash_adaptor);
//	//输出result的所有内容
//	for (size_t i = 0; i < kv_size(result); ++i) {
//		std::cout << kv_A(result, i) << " ";
//	}
//	std::cout << std::endl;
//	auto result2 = icf.bulk_count(sequence2 | hash_adaptor);
//	//输出result2的所有内容
//	for (size_t i = 0; i < kv_size(result2); ++i) {
//		std::cout << kv_A(result2, i) << " ";
//	}
//	std::cout << std::endl;
//	auto result3 = icf.bulk_count(sequence3 | hash_adaptor);
//	//输出result3的所有内容
//	for (size_t i = 0; i < kv_size(result3); ++i) {
//		std::cout << kv_A(result3, i) << " ";
//	}
//	std::cout << std::endl;
//
//	//记录查询时间
//	end = std::chrono::high_resolution_clock::now();
//	elapsed_seconds = end - start;
//	std::cout << "ICF的查询时间：" << elapsed_seconds.count() << "s\n";
//
//	//序列化保存ibf和icf
//	{
//		std::ofstream os("ibf.cereal", std::ios::binary);
//		cereal::BinaryOutputArchive archive(os);
//		archive(ibf);
//	}
//	{
//		std::ofstream os("icf.cereal", std::ios::binary);
//		cereal::BinaryOutputArchive archive(os);
//		archive(icf);
//	}
//
//	return 0;
//}

using namespace sdsl;

// 使用掩码批量插入八位数到 bit_vector
void batch_insert_to_bit_vector(bit_vector& bv, uint8_t value, size_t position) {
	uint64_t mask = static_cast<uint64_t>(value) << (position % 64);  // 将8位数移到正确位置
	size_t idx = position / 64;
	bv.data()[idx] |= mask;  // 批量写入
	if ((position % 64) > 56) {
		// 处理跨越64位边界的情况
		bv.data()[idx + 1] |= (value >> (64 - (position % 64)));
	}
}

// 插入八位数到 int_vector<8>
void insert_to_int_vector(int_vector<8>& iv, uint8_t value, size_t position) {
	iv[position] = value;
}

// 查询 bit_vector 中的值
uint8_t query_bit_vector(bit_vector& bv, size_t position) {
	uint64_t mask = 0xFFULL << (position % 64);
	size_t idx = position / 64;
	uint64_t chunk = (bv.data()[idx] & mask) >> (position % 64);
	if ((position % 64) > 56) {
		chunk |= (bv.data()[idx + 1] & 0xFFULL) << (64 - (position % 64));
	}
	return static_cast<uint8_t>(chunk);
}

// 查询 int_vector<8> 中的值
uint8_t query_int_vector(int_vector<8>& iv, size_t position) {
	return iv[position];
}

int main() {
	const size_t num_elements = 100000000;  // 测试元素数量
	const uint8_t test_value = 0b10101010; // 测试插入值

	// 构建 bit_vector 和 int_vector<8>
	auto start_bv_construct = std::chrono::high_resolution_clock::now();
	bit_vector bv(num_elements * 8, 0); // 每个元素8位
	auto end_bv_construct = std::chrono::high_resolution_clock::now();

	auto start_iv_construct = std::chrono::high_resolution_clock::now();
	int_vector<8> iv(num_elements, 0);
	auto end_iv_construct = std::chrono::high_resolution_clock::now();

	// 测试 bit_vector 的插入性能
	auto start_bv_insert = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < num_elements; ++i) {
		batch_insert_to_bit_vector(bv, test_value, i * 8);
	}
	auto end_bv_insert = std::chrono::high_resolution_clock::now();

	// 测试 int_vector<8> 的插入性能
	auto start_iv_insert = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < num_elements; ++i) {
		insert_to_int_vector(iv, test_value, i);
	}
	auto end_iv_insert = std::chrono::high_resolution_clock::now();

	// 测试 bit_vector 的查询性能
	auto start_bv_query = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < num_elements; ++i) {
		volatile uint8_t value = query_bit_vector(bv, i * 8); // 使用volatile防止优化
	}
	auto end_bv_query = std::chrono::high_resolution_clock::now();

	// 测试 int_vector<8> 的查询性能
	auto start_iv_query = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < num_elements; ++i) {
		volatile uint8_t value = query_int_vector(iv, i);
	}
	auto end_iv_query = std::chrono::high_resolution_clock::now();

	// 测试 bit_vector 的 get_int 性能
	auto start_bv_get_int = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < num_elements; ++i) {
		size_t index = i * 8;
		volatile uint64_t value = bv.get_int(index, 64); // 使用volatile防止优化
	}
	auto end_bv_get_int = std::chrono::high_resolution_clock::now();

	// 测试 int_vector<8> 的 get_int 性能
	auto start_iv_get_int = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < num_elements; ++i) {
		size_t index = i * 8;
		volatile uint8_t value = iv.get_int(index, 64);
	}
	auto end_iv_get_int = std::chrono::high_resolution_clock::now();

	// 输出测试结果
	std::cout << "bit_vector 构建时间: " << std::chrono::duration<double>(end_bv_construct - start_bv_construct).count() << " 秒" << std::endl;
	std::cout << "int_vector<8> 构建时间: " << std::chrono::duration<double>(end_iv_construct - start_iv_construct).count() << " 秒" << std::endl;

	std::cout << "bit_vector 插入时间: " << std::chrono::duration<double>(end_bv_insert - start_bv_insert).count() << " 秒" << std::endl;
	std::cout << "int_vector<8> 插入时间: " << std::chrono::duration<double>(end_iv_insert - start_iv_insert).count() << " 秒" << std::endl;

	std::cout << "bit_vector 查询时间: " << std::chrono::duration<double>(end_bv_query - start_bv_query).count() << " 秒" << std::endl;
	std::cout << "int_vector<8> 查询时间: " << std::chrono::duration<double>(end_iv_query - start_iv_query).count() << " 秒" << std::endl;

	std::cout << "bit_vector get_int 时间: " << std::chrono::duration<double>(end_bv_get_int - start_bv_get_int).count() << " 秒" << std::endl;
	std::cout << "int_vector<8> get_int 时间: " << std::chrono::duration<double>(end_iv_get_int - start_iv_get_int).count() << " 秒" << std::endl;

	return 0;
}