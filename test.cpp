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

KSEQ_INIT(gzFile, gzread)

std::vector<std::pair<std::string, uint32_t>> read_target_info(const std::string& filepath) {
	std::ifstream file(filepath);
	std::vector<std::pair<std::string, uint32_t>> targets;
	std::string line;

	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string path;
		uint32_t taxid;
		if (!(iss >> path >> taxid)) { break; }
		targets.emplace_back(path, taxid);
	}

	return targets;
}

std::vector<seqan3::dna4> read_fna_gz(const std::string& filepath) {
	using namespace seqan3::literals;
	gzFile fp = gzopen(filepath.c_str(), "r");
	kseq_t* seq = kseq_init(fp);
	std::vector<seqan3::dna4> sequence;

	while (kseq_read(seq) >= 0) {
		for (size_t i = 0; i < seq->seq.l; ++i) {
			switch (seq->seq.s[i]) {
			case 'A': sequence.push_back('A'_dna4); break;
			case 'C': sequence.push_back('C'_dna4); break;
			case 'G': sequence.push_back('G'_dna4); break;
			case 'T': sequence.push_back('T'_dna4); break;
			default: break;
			}
		}
	}

	kseq_destroy(seq);
	gzclose(fp);

	return sequence;
}

int main() {
	std::string target_info_file = "/mnt/d/code/src/ganon/test_files/build/target_info.tsv";
	auto targets = read_target_info(target_info_file);

	size_t total_items = 100000;
	cuckoofilter::CuckooFilter<size_t, 12> filter(total_items);

	for (const auto& [path, taxid] : targets) {
		auto sequence = read_fna_gz(path);
		auto minimizer_hash = sequence | seqan3::views::minimiser_hash(seqan3::shape{ seqan3::ungapped{4} }, seqan3::window_size{ 8 });
		for (auto hash : minimizer_hash) {
			filter.Add(hash);
		}
	}

	// 序列化 Cuckoo Filter
	{
		std::ofstream os("cuckoofilter.cereal", std::ios::binary);
		cereal::BinaryOutputArchive archive(os);
		archive(filter);
	}

	// 反序列化 Cuckoo Filter
	cuckoofilter::CuckooFilter<size_t, 12> loaded_filter(total_items); // 提供必要的参数
	{
		std::ifstream is("cuckoofilter.cereal", std::ios::binary);
		cereal::BinaryInputArchive archive(is);
		archive(loaded_filter);
	}

	// 检查第一个 hash 是否存在于反序列化后的过滤器中
	if (!targets.empty()) {
		auto first_sequence = read_fna_gz(targets[0].first);
		auto first_minimizer_hash = first_sequence | seqan3::views::minimiser_hash(seqan3::shape{ seqan3::ungapped{4} }, seqan3::window_size{ 8 });
		auto first_hash = *first_minimizer_hash.begin();
		if (loaded_filter.Contain(first_hash) == cuckoofilter::Ok) {
			std::cout << "Hash exists in the loaded filter.\n";
		}
		else {
			std::cout << "Hash does not exist in the loaded filter.\n";
		}
	}

	//输出两个过滤器的各种大小属性
	std::cout << "Original filter: " << filter.SizeInBytes() << " bytes\n";
	std::cout << "Loaded filter: " << loaded_filter.SizeInBytes() << " bytes\n";
	std::cout << "Original filter: " << filter.Size() << " items\n";
	std::cout << "Loaded filter: " << loaded_filter.Size() << " items\n";
	std::cout << "Original filter: " << filter.Info() << "\n";
	std::cout << "Loaded filter: " << loaded_filter.Info() << "\n";

	return 0;
}
