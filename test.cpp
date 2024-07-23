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
//	// ��¼����ʱ��
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
//	// ���л� Cuckoo Filter
//	{
//		std::ofstream os("cuckoofilter.cereal", std::ios::binary);
//		cereal::BinaryOutputArchive archive(os);
//		archive(filter);
//		archive(sdsltest);
//	}
//
//	// �����л� Cuckoo Filter
//	cuckoofilter::CuckooFilter<size_t, 12> loaded_filter(total_items); // �ṩ��Ҫ�Ĳ���
//	{
//		std::ifstream is("cuckoofilter.cereal", std::ios::binary);
//		cereal::BinaryInputArchive archive(is);
//		archive(loaded_filter);
//	}
//
//	// ����һ�� hash �Ƿ�����ڷ����л���Ĺ�������
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
//	//��������������ĸ��ִ�С����
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
//using cuckoofilter::CuckooFilter;   // ʹ��cuckoofilter�����ռ��е�CuckooFilter��
//
//int main(int argc, char** argv) {
//	size_t total_items = 1;  // Ҫ���������������Ŀ��
//
//
//	// ����һ��Cuckoo Filter����Ŀ����Ϊsize_t��ÿ����Ŀʹ��12λ��
//	// CuckooFilter<size_t, 12> filter(total_items);
//	// Ҫ���ð����򣬿��Խ�Cuckoo Filter�Ĵ洢����ΪPackedTable������size_t���͵ļ�����Ϊÿ��������13λ��
//	// CuckooFilter<size_t, 13, cuckoofilter::PackedTable> filter(total_items);
//	CuckooFilter<size_t, 12> filter(total_items); // ʵ����һ��Cuckoo Filter����
//
//	// �����Cuckoo Filter�в�����Ŀ
//	size_t num_inserted = 0;  // ��¼�ɹ��������Ŀ����
//	for (size_t i = 0; i < total_items; i++, num_inserted++) {
//		if (filter.Add(i) != cuckoofilter::Ok) {    // ���������Ŀ��������ʧ��������ѭ��
//			break;
//		}
//	}
//
//	// ���֮ǰ�������Ŀ�Ƿ��ڹ������У�����������Ŀ������
//	for (size_t i = 0; i < num_inserted; i++) {
//		assert(filter.Contain(i) == cuckoofilter::Ok);  // ������ʧ�ܣ��������ֹ
//	}
//
//	// ��鲻���ڵ���Ŀ��Ԥ�ƻ���һЩ��
//	size_t total_queries = 0; // �ܲ�ѯ����
//	size_t false_queries = 0; // �󱨴���
//	for (size_t i = total_items; i < 2 * total_items; i++) {
//		if (filter.Contain(i) == cuckoofilter::Ok) {
//			false_queries++;
//		}
//		total_queries++;
//	}
//
//	// �������������
//	std::cout << "false positive rate is "
//		<< 100.0 * false_queries / total_queries << "%\n";
//
//	return 0;
//}
#include <iostream>
#include <chrono>
#include <random>

int main() {
	// ģ�������
	auto start_modulo = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < 1000000; ++i) {
		int result = i % 10;
	}
	auto end_modulo = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_modulo = end_modulo - start_modulo;
	std::cout << "Elapsed time for modulo operation: " << elapsed_modulo.count() << " s\n";

	// ���������
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