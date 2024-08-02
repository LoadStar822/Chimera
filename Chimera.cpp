/*
 * -----------------------------------------------------------------------------
 * Filename:      Chimera.cpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-07-09
 *
 * Last Modified: 2024-07-09
 *
 * Description:
 *  This is a simple C++ program that outputs "Hello, World!".
 *
 * Version:
 *  1.0
 * -----------------------------------------------------------------------------
 */
#include "Chimera.h"

using cuckoofilter::CuckooFilter;

namespace chimera
{
	//创建一个结构来存储路径和taxid
	typedef struct {
		kstring_t path;
		int taxid;
	} fileInput;

	// 定义别名
	typedef kvec_t(fileInput) fileInputVec;

	// 读取文件并处理最小化器
	// 读取文件并处理最小化器
	void process_files(fileInputVec& test, cuckoofilter::CuckooFilter<uint64_t, 12>& filter, const Config& config) {
		auto minimiser_view = seqan3::views::minimiser_hash(seqan3::shape{ seqan3::ungapped{ config.kmer_size } },
			seqan3::window_size{ config.window_size },
			seqan3::seed{ 0x8F3F73B5CF1C9ADE });

		for (size_t i = 0; i < kv_size(test); ++i) {
			fileInput fi = kv_A(test, i);
			std::cout << "Processing file: " << fi.path.s << std::endl;

			gzFile gzfp = gzopen(fi.path.s, "r");
			if (!gzfp) {
				fprintf(stderr, "Failed to open gzipped file: %s\n", fi.path.s);
				continue;
			}

			kseq_t* seq = kseq_init(gzfp);
			while (kseq_read(seq) >= 0) {
				if (seq->seq.l < config.min_length) {
					std::cout << "Skipped short sequence of length: " << seq->seq.l << std::endl;
					continue;
				}

				std::vector<seqan3::dna4> dna_sequence;
				dna_sequence.reserve(seq->seq.l);
				for (size_t j = 0; j < seq->seq.l; ++j) {
					dna_sequence.push_back(seqan3::assign_char_to(seq->seq.s[j], seqan3::dna4{}));
				}

				std::cout << "Processing sequence of length: " << seq->seq.l << std::endl;
				auto minimisers = dna_sequence | minimiser_view | std::views::common;
				for (auto minimiser : minimisers) {
					uint64_t minimiser_hash = static_cast<uint64_t>(minimiser);
					std::cout << "Minimiser hash: " << minimiser_hash << std::endl;
					if (filter.Add(minimiser_hash) != cuckoofilter::Ok) {
						std::cerr << "Failed to insert item into the filter\n";
					}
					else {
						std::cout << "Inserted minimiser hash: " << minimiser_hash << std::endl;
					}
				}
			}

			kseq_destroy(seq);
			gzclose(gzfp);
		}
	}
}

int main(int argc, char** argv) {
	chimera::fileInputVec test;
	kv_init(test);
	// 打开文件
	FILE* fp = fopen("/mnt/d/code/src/ganon/test_files/build/target_info.tsv", "r");
	if (!fp) {
		perror("Failed to open file");
		return EXIT_FAILURE;
	}

	// 读取文件每一行并解析
	char* line = NULL;
	size_t len = 0;
	ssize_t read;
	while ((read = getline(&line, &len, fp)) != -1) {
		chimera::fileInput fi;
		fi.path.l = fi.path.m = 0;
		fi.path.s = NULL;

		// 解析路径
		char* token = strtok(line, "\t");
		if (token) {
			kputs(token, &fi.path);
		}

		// 解析 taxid
		token = strtok(NULL, "\n");
		if (token) {
			fi.taxid = atoi(token);
		}

		// 将解析结果添加到 kvec_t 中
		kv_push(chimera::fileInput, test, fi);
	}
	free(line);
	fclose(fp);

	//配置参数
	chimera::Config config;

	// 初始化 Cuckoo Filter
	size_t total_items = 100000; // 根据需要设置大小
	cuckoofilter::CuckooFilter<uint64_t, 12> filter(total_items);

	// 处理文件
	chimera::process_files(test, filter, config);

	// 打印Cuckoo Filter中的项目数量
	size_t num_inserted = filter.Size();
	std::cerr << "Inserted " << num_inserted << " items into the filter.\n";

	// 释放 kvec_t 中的字符串内存
	for (size_t i = 0; i < kv_size(test); ++i) {
		free(kv_A(test, i).path.s);
	}
	kv_destroy(test);

	return 0;
}