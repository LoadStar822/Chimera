#include <buildConfig.hpp>
#include <ChimeraBuild.hpp>
#include <iostream>
#include <chrono>
#include <interleaved-cuckoo-filter.h>
#include <iomanip>
#include "khash.h"

namespace ChimeraBuild {
	KHASH_MAP_INIT_STR(str, uint64_t)
		// �������
		void khash_insert(khash_t(str)* h, const std::string& key, uint64_t value) {
		int ret;
		khiter_t k = kh_put(str, h, key.c_str(), &ret);
		if (ret == -1) {
			std::cerr << "Error inserting key: " << key << std::endl;
			return;
		}
		kh_value(h, k) = value;
	}

	// ��������
	bool khash_find(khash_t(str)* h, const std::string& key, uint64_t& value) {
		khiter_t k = kh_get(str, h, key.c_str());
		if (k == kh_end(h)) {
			return false; // δ�ҵ�
		}
		value = kh_value(h, k);
		return true;
	}

	// ��������
	void khash_iterate(khash_t(str)* h) {
		for (khiter_t k = kh_begin(h); k != kh_end(h); ++k) {
			if (kh_exist(h, k)) {
				const char* key = kh_key(h, k);
				uint64_t value = kh_value(h, k);
				std::cout << "Key: " << key << ", Value: " << value << std::endl;
			}
		}
	}

	void print_build_time(long long milliseconds) {
		// �����롢���Ӻ�Сʱ
		long long total_seconds = milliseconds / 1000;
		long long seconds = total_seconds % 60;
		long long total_minutes = total_seconds / 60;
		long long minutes = total_minutes % 60;
		long long hours = total_minutes / 60;

		// ����ʱ�䳤��ѡ��ͬ�ĸ�ʽ���
		if (hours > 0) {
			std::cout << hours << "h " << minutes << "min " << seconds << "s" << std::endl;
		}
		else {
			std::cout << minutes << "min " << seconds << "s" << std::endl;
		}
	}

	void run(BuildConfig config)
	{
		if (config.verbose)
		{
			std::cout << config << std::endl;
		}
		auto build_start = std::chrono::high_resolution_clock::now();

		ICFConfig icfConfig;
		icfConfig.kmer_size = config.kmer_size;
		icfConfig.window_size = config.window_size;

		auto build_end = std::chrono::high_resolution_clock::now();

		// �����ܺ�ʱ���Ժ���Ϊ��λ
		auto build_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(build_end - build_start).count();
	}
}