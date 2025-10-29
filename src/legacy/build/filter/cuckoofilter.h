#ifndef CUCKOO_FILTER_CUCKOO_FILTER_H_
#define CUCKOO_FILTER_CUCKOO_FILTER_H_

#include <assert.h>
#include <algorithm>

#include "debug.h"
#include "hashutil.h"
#include "packedtable.h"
#include "printutil.h"
#include "singletable.h"

namespace cuckoofilter {
// 定义枚举类型Status表示布谷鸟过滤器操作的状态
enum Status {
	Ok = 0,             // 操作成功
	NotFound = 1,       // 项未找到
	NotEnoughSpace = 2, // 空间不足
	NotSupported = 3,   // 操作不支持
};

// 定义布谷鸟过滤器操作失败前最大的布谷鸟踢出次数
const size_t kMaxCuckooCount = 500;

// 布谷鸟过滤器类，提供Bloomier过滤器接口，支持Add, Delete, Contain操作
// 该类接受三个模板参数：
//   ItemType: 要插入项的类型
//   bits_per_item: 每个项被哈希到的位数
//   TableType: 表的存储类型，默认为SingleTable，可选PackedTable以启用半排序
template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType = SingleTable,
          typename HashFamily = TwoIndependentMultiplyShift>
class CuckooFilter {
  // 存储项的表
  TableType<bits_per_item> *table_;

  // 存储的项数量
  size_t num_items_;

  // 受害者缓存结构体，记录暂时无法插入的项
  typedef struct {
	  size_t index; // 项的索引
	  uint32_t tag; // 项的标签
	  bool used;    // 标记该缓存是否被使用
  } VictimCache;


  // 受害者缓存
  VictimCache victim_;

  // 哈希函数族
  HashFamily hasher_;

  // 内联函数，计算索引哈希
  inline size_t IndexHash(uint32_t hv) const {
	  // table_->num_buckets始终是2的幂，所以模运算可以用按位与代替：
	  return hv & (table_->NumBuckets() - 1);
  }

  // 内联函数，计算标签哈希
  inline uint32_t TagHash(uint32_t hv) const {
	  uint32_t tag;
	  // 将哈希值与bits_per_item的掩码进行按位与运算，得到标签
	  tag = hv & ((1ULL << bits_per_item) - 1);
	  // 如果标签为0，增加1以避免使用0标签
	  tag += (tag == 0);
	  return tag;
  }

  // 内联函数，生成项的索引和标签哈希
  inline void GenerateIndexTagHash(const ItemType& item, size_t* index, uint32_t* tag) const {
	  // 使用哈希函数生成64位哈希值
	  const uint64_t hash = hasher_(item);
	  // 高32位用于计算索引
	  *index = IndexHash(hash >> 32);
	  // 低32位用于计算标签
	  *tag = TagHash(hash);
  }


  // 内联函数，计算替代索引
  inline size_t AltIndex(const size_t index, const uint32_t tag) const {
	  // 使用MurmurHash2的哈希常量0x5bd1e995进行快速哈希
	  return IndexHash((uint32_t)(index ^ (tag * 0x5bd1e995)));
  }

  // 添加项的内部实现，传入索引和标签
  Status AddImpl(const size_t i, const uint32_t tag);

  // 计算装载因子，表示占用率
  double LoadFactor() const {
	  return 1.0 * Size() / table_->SizeInTags();
  }

  // 计算每项的平均位数
  double BitsPerItem() const {
	  return 8.0 * table_->SizeInBytes() / Size();
  }

 public:

  // 构造函数，接收最大键数作为参数
  explicit CuckooFilter(const size_t max_num_keys) : num_items_(0), victim_(), hasher_() {
	  size_t assoc = 4; // 每个桶的关联数
	  size_t num_buckets = upperpower2(std::max<uint64_t>(1, max_num_keys / assoc)); // 计算桶的数量
	  double frac = (double)max_num_keys / num_buckets / assoc; // 计算每个桶的装载因子
	  if (frac > 0.96) { // 如果装载因子超过0.96，桶的数量加倍
		  num_buckets <<= 1;
	  }
	  victim_.used = false; // 初始化受害者缓存为未使用
	  table_ = new TableType<bits_per_item>(num_buckets); // 初始化表
  }

  // 析构函数，释放表的内存
  ~CuckooFilter() { delete table_; }

  // 向过滤器添加一个项
  Status Add(const ItemType& item);

  // 检查一个项是否存在（可能有假阳性率）
  Status Contain(const ItemType& item) const;

  // 从过滤器中删除一个项
  Status Delete(const ItemType& item);

  /* 提供统计信息的方法 */
  // 返回摘要信息
  std::string Info() const;

  // 返回当前插入的项数
  size_t Size() const { return num_items_; }

  // 返回过滤器的字节大小
  size_t SizeInBytes() const { return table_->SizeInBytes(); }

  // 序列化方法，用于持久化过滤器
  template <class Archive>
  void serialize(Archive& ar) {
	  ar(num_items_, victim_.index, victim_.tag, victim_.used);
	  ar(*table_);
      ar(hasher_);
  }
};

// 模板类CuckooFilter的成员函数Add的实现
template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::Add(
    const ItemType &item) {
  size_t i;
  uint32_t tag;

  // 如果受害者缓存已被使用，返回空间不足状态
  if (victim_.used) {
	  return NotEnoughSpace;
  }

  // 生成项的索引和标签哈希值
  GenerateIndexTagHash(item, &i, &tag);
  // 调用内部实现方法AddImpl
  return AddImpl(i, tag);
}

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::AddImpl(
    const size_t i, const uint32_t tag) {
  size_t curindex = i; // 当前索引
  uint32_t curtag = tag; // 当前标签
  uint32_t oldtag;

  // 尝试最多kMaxCuckooCount次
  for (uint32_t count = 0; count < kMaxCuckooCount; count++) {
	  bool kickout = count > 0; // 如果不是第一次尝试，则允许踢出
	  oldtag = 0;
	  // 尝试将标签插入到当前桶中
	  if (table_->InsertTagToBucket(curindex, curtag, kickout, oldtag)) {
		  num_items_++; // 成功插入后，增加项数
		  return Ok; // 返回操作成功状态
	  }
	  // 如果踢出了一个标签，将其设为当前标签
	  if (kickout) {
		  curtag = oldtag;
	  }
	  // 计算替代索引
	  curindex = AltIndex(curindex, curtag);
  }

  // 如果所有尝试均失败，将当前索引和标签存入受害者缓存
  victim_.index = curindex;
  victim_.tag = curtag;
  victim_.used = true;
  return Ok; 
}

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::Contain(
    const ItemType &key) const {
  bool found = false;
  size_t i1, i2;
  uint32_t tag;

  // 生成键的索引和标签哈希值
  GenerateIndexTagHash(key, &i1, &tag);
  // 计算替代索引
  i2 = AltIndex(i1, tag);

  // 确保替代索引计算正确
  assert(i1 == AltIndex(i2, tag));

  // 检查受害者缓存中是否包含该项
  found = victim_.used && (tag == victim_.tag) &&
	  (i1 == victim_.index || i2 == victim_.index);

  // 如果在受害者缓存中找到或者在桶中找到，返回Ok状态
  if (found || table_->FindTagInBuckets(i1, i2, tag)) {
	  return Ok;
  }
  else {
	  return NotFound; // 否则返回未找到状态
  }
}

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::Delete(
    const ItemType &key) {
  size_t i1, i2;
  uint32_t tag;
  // 生成键的索引和标签哈希值
  GenerateIndexTagHash(key, &i1, &tag);
  // 计算替代索引
  i2 = AltIndex(i1, tag);

  // 尝试从第一个桶中删除标签
  if (table_->DeleteTagFromBucket(i1, tag)) {
	  num_items_--; // 删除成功，减少项数
	  goto TryEliminateVictim;
  }
  // 尝试从第二个桶中删除标签
  else if (table_->DeleteTagFromBucket(i2, tag)) {
	  num_items_--; // 删除成功，减少项数
	  goto TryEliminateVictim;
  }
  // 检查受害者缓存中是否包含该项
  else if (victim_.used && tag == victim_.tag &&
	  (i1 == victim_.index || i2 == victim_.index)) {
	  victim_.used = false; // 删除成功，标记受害者缓存未使用
	  return Ok;
  }
  // 如果以上删除操作均未成功，返回未找到状态
  else {
	  return NotFound;
  }

TryEliminateVictim:
  // 如果受害者缓存中有未处理的项，尝试将其插入
  if (victim_.used) {
	  victim_.used = false;
	  size_t i = victim_.index;
	  uint32_t tag = victim_.tag;
	  AddImpl(i, tag); // 尝试重新插入受害者缓存中的项
  }
  return Ok; // 返回操作成功状态
}

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
std::string CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::Info() const {
  std::stringstream ss;
  // 输出布谷鸟过滤器的状态信息
  ss << "CuckooFilter Status:\n"
	  << "\t\t" << table_->Info() << "\n" // 输出表的信息
	  << "\t\tKeys stored: " << Size() << "\n" // 输出存储的键数
	  << "\t\tLoad factor: " << LoadFactor() << "\n" // 输出负载因子
	  << "\t\tHashtable size: " << (table_->SizeInBytes() >> 10) << " KB\n"; // 输出哈希表的大小（以KB为单位）
  if (Size() > 0) {
	  ss << "\t\tbit/key:   " << BitsPerItem() << "\n"; // 如果有存储项，输出每个键的平均位数
  }
  else {
	  ss << "\t\tbit/key:   N/A\n"; // 如果没有存储项，输出N/A
  }
  return ss.str();
}
}  // namespace cuckoofilter
#endif  // CUCKOO_FILTER_CUCKOO_FILTER_H_
