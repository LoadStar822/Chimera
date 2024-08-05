#include <buildConfig.hpp>
#include <ChimeraBuild.hpp>
#include <iostream>
#include <chrono>
#include <interleaved-cuckoo-filter.h>
#include <iomanip>
#include "khash.h"

namespace ChimeraBuild {
	KHASH_MAP_INIT_STR(str, uint64_t)
		// 插入操作
		void khash_insert(khash_t(str)* h, const std::string& key, uint64_t value) {
		int ret;
		khiter_t k = kh_put(str, h, key.c_str(), &ret);
		if (ret == -1) {
			std::cerr << "Error inserting key: " << key << std::endl;
			return;
		}
		kh_value(h, k) = value;
	}

	// 检索操作
	bool khash_find(khash_t(str)* h, const std::string& key, uint64_t& value) {
		khiter_t k = kh_get(str, h, key.c_str());
		if (k == kh_end(h)) {
			return false; // 未找到
		}
		value = kh_value(h, k);
		return true;
	}

	// 遍历操作
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
		// 计算秒、分钟和小时
		long long total_seconds = milliseconds / 1000;
		long long seconds = total_seconds % 60;
		long long total_minutes = total_seconds / 60;
		long long minutes = total_minutes % 60;
		long long hours = total_minutes / 60;

		// 根据时间长度选择不同的格式输出
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

		// 计算总耗时，以毫秒为单位
		auto build_total_time = std::chrono::duration_cast<std::chrono::milliseconds>(build_end - build_start).count();
	}
}