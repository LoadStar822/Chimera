// Chimera.cpp: 定义应用程序的入口点。
//

#include "Chimera.h"

#include "src/build/filter/cuckoofilter.h"  // 引入Cuckoo Filter的头文件

#include <assert.h>        // 引入断言库，用于调试
#include <math.h>          // 引入数学库
#include <iostream>        // 引入输入输出流库
#include <vector>          // 引入向量库

using cuckoofilter::CuckooFilter;   // 使用cuckoofilter命名空间中的CuckooFilter类

int main(int argc, char** argv) {
	size_t total_items = 58546453;  // 要插入过滤器的总项目数


	// 创建一个Cuckoo Filter，项目类型为size_t，每个项目使用12位：
	// CuckooFilter<size_t, 12> filter(total_items);
	// 要启用半排序，可以将Cuckoo Filter的存储定义为PackedTable，接受size_t类型的键，并为每个键分配13位：
	// CuckooFilter<size_t, 13, cuckoofilter::PackedTable> filter(total_items);
	CuckooFilter<size_t, 12> filter(total_items); // 实例化一个Cuckoo Filter对象

	// 向这个Cuckoo Filter中插入项目
	size_t num_inserted = 0;  // 记录成功插入的项目数量
	for (size_t i = 0; i < total_items; i++, num_inserted++) {
		if (filter.Add(i) != cuckoofilter::Ok) {    // 尝试添加项目，如果添加失败则跳出循环
			break;
		}
	}
	std::cerr << "Inserted " << num_inserted << " items into the filter.\n";

	// 检查之前插入的项目是否在过滤器中，期望所有项目都存在
	for (size_t i = 0; i < num_inserted; i++) {
		assert(filter.Contain(i) == cuckoofilter::Ok);  // 如果检查失败，程序会终止
	}
	std::cerr << "All inserted items are in the filter.\n";

	// 检查不存在的项目，预计会有一些误报
	size_t total_queries = 0; // 总查询次数
	size_t false_queries = 0; // 误报次数
	for (size_t i = total_items; i < 2 * total_items; i++) {
		if (filter.Contain(i) == cuckoofilter::Ok) {
			false_queries++;
		}
		total_queries++;
	}

	// 输出测量的误报率
	std::cerr << "false positive rate is "
		<< 100.0 * false_queries / total_queries << "%\n";

	
	return 0;
}
