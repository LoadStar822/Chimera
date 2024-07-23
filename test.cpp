#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <zlib.h>
#include "kseq.h"
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include "cuckoofilter.h"
#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/memory.hpp>

#include <sdsl/bit_vectors.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <kvec.h>
using namespace seqan3::literals;
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
//	size_t total_items = 100000;
//	cuckoofilter::CuckooFilter<size_t, 12> filter(total_items);
//	sdsl::int_vector<8> sdsltest(100000000, 0);
//	for (const auto& [path, taxid] : targets) {
//		auto sequence = read_fna_gz(path);
//		auto minimizer_hash = sequence | seqan3::views::minimiser_hash(seqan3::shape{ seqan3::ungapped{4} }, seqan3::window_size{ 8 });
//		auto i = 0;
//		for (auto hash : minimizer_hash) {
//			filter.Add(hash);
//			sdsltest[i++] = hash;
//		}
//	}
//
//	// 序列化 Cuckoo Filter
//	{
//		std::ofstream os("cuckoofilter.cereal", std::ios::binary);
//		cereal::BinaryOutputArchive archive(os);
//		archive(filter);
//		archive(sdsltest);
//	}
//
//	// 反序列化 Cuckoo Filter
//	cuckoofilter::CuckooFilter<size_t, 12> loaded_filter(total_items); // 提供必要的参数
//	{
//		std::ifstream is("cuckoofilter.cereal", std::ios::binary);
//		cereal::BinaryInputArchive archive(is);
//		archive(loaded_filter);
//	}
//
//	// 检查第一个 hash 是否存在于反序列化后的过滤器中
//	if (!targets.empty()) {
//		auto first_sequence = read_fna_gz(targets[0].first);
//		auto first_minimizer_hash = first_sequence | seqan3::views::minimiser_hash(seqan3::shape{ seqan3::ungapped{4} }, seqan3::window_size{ 8 });
//		auto first_hash = *first_minimizer_hash.begin();
//		if (loaded_filter.Contain(first_hash) == cuckoofilter::Ok) {
//			std::cout << "Hash exists in the loaded filter.\n";
//		}
//		else {
//			std::cout << "Hash does not exist in the loaded filter.\n";
//		}
//	}
//
//	//输出两个过滤器的各种大小属性
//	std::cout << "Original filter: " << filter.SizeInBytes() << " bytes\n";
//	std::cout << "Loaded filter: " << loaded_filter.SizeInBytes() << " bytes\n";
//	std::cout << "Original filter: " << filter.Size() << " items\n";
//	std::cout << "Loaded filter: " << loaded_filter.Size() << " items\n";
//	std::cout << "Original filter: " << filter.Info() << "\n";
//	std::cout << "Loaded filter: " << loaded_filter.Info() << "\n";
//	auto end = std::chrono::high_resolution_clock::now();
//	std::chrono::duration<double> elapsed = end - start;
//	std::cout << "Elapsed time: " << elapsed.count() << " s\n";
//	return 0;
//}
//using cuckoofilter::CuckooFilter;   // 使用cuckoofilter命名空间中的CuckooFilter类
//
//int main(int argc, char** argv) {
//	size_t total_items = 1;  // 要插入过滤器的总项目数
//
//
//	// 创建一个Cuckoo Filter，项目类型为size_t，每个项目使用12位：
//	// CuckooFilter<size_t, 12> filter(total_items);
//	// 要启用半排序，可以将Cuckoo Filter的存储定义为PackedTable，接受size_t类型的键，并为每个键分配13位：
//	// CuckooFilter<size_t, 13, cuckoofilter::PackedTable> filter(total_items);
//	CuckooFilter<size_t, 12> filter(total_items); // 实例化一个Cuckoo Filter对象
//
//	// 向这个Cuckoo Filter中插入项目
//	size_t num_inserted = 0;  // 记录成功插入的项目数量
//	for (size_t i = 0; i < total_items; i++, num_inserted++) {
//		if (filter.Add(i) != cuckoofilter::Ok) {    // 尝试添加项目，如果添加失败则跳出循环
//			break;
//		}
//	}
//
//	// 检查之前插入的项目是否在过滤器中，期望所有项目都存在
//	for (size_t i = 0; i < num_inserted; i++) {
//		assert(filter.Contain(i) == cuckoofilter::Ok);  // 如果检查失败，程序会终止
//	}
//
//	// 检查不存在的项目，预计会有一些误报
//	size_t total_queries = 0; // 总查询次数
//	size_t false_queries = 0; // 误报次数
//	for (size_t i = total_items; i < 2 * total_items; i++) {
//		if (filter.Contain(i) == cuckoofilter::Ok) {
//			false_queries++;
//		}
//		total_queries++;
//	}
//
//	// 输出测量的误报率
//	std::cout << "false positive rate is "
//		<< 100.0 * false_queries / total_queries << "%\n";
//
//	return 0;
//}
#include <iostream>
#include <chrono>
#include <random>

int main() {
	// 模运算测试
	auto start_modulo = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 1000000; ++i) {
		int result = i % 10;
	}
	auto end_modulo = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_modulo = end_modulo - start_modulo;
	std::cout << "Elapsed time for modulo operation: " << elapsed_modulo.count() << " s\n";

	// 随机数测试
	auto start_random = std::chrono::high_resolution_clock::now();
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> dis(0, 9);
	for (int i = 0; i < 1000000; ++i) {
		int result = dis(gen);
	}
	auto end_random = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_random = end_random - start_random;
	std::cout << "Elapsed time for random number generation: " << elapsed_random.count() << " s\n";

	return 0;
}