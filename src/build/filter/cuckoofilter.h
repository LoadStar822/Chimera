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
// ����ö������Status��ʾ�����������������״̬
enum Status {
	Ok = 0,             // �����ɹ�
	NotFound = 1,       // ��δ�ҵ�
	NotEnoughSpace = 2, // �ռ䲻��
	NotSupported = 3,   // ������֧��
};

// ���岼�������������ʧ��ǰ���Ĳ������߳�����
const size_t kMaxCuckooCount = 500;

// ������������࣬�ṩBloomier�������ӿڣ�֧��Add, Delete, Contain����
// �����������ģ�������
//   ItemType: Ҫ�����������
//   bits_per_item: ÿ�����ϣ����λ��
//   TableType: ��Ĵ洢���ͣ�Ĭ��ΪSingleTable����ѡPackedTable�����ð�����
template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType = SingleTable,
          typename HashFamily = TwoIndependentMultiplyShift>
class CuckooFilter {
  // �洢��ı�
  TableType<bits_per_item> *table_;

  // �洢��������
  size_t num_items_;

  // �ܺ��߻���ṹ�壬��¼��ʱ�޷��������
  typedef struct {
	  size_t index; // �������
	  uint32_t tag; // ��ı�ǩ
	  bool used;    // ��Ǹû����Ƿ�ʹ��
  } VictimCache;


  // �ܺ��߻���
  VictimCache victim_;

  // ��ϣ������
  HashFamily hasher_;

  // ��������������������ϣ
  inline size_t IndexHash(uint32_t hv) const {
	  // table_->num_bucketsʼ����2���ݣ�����ģ��������ð�λ����棺
	  return hv & (table_->NumBuckets() - 1);
  }

  // ���������������ǩ��ϣ
  inline uint32_t TagHash(uint32_t hv) const {
	  uint32_t tag;
	  // ����ϣֵ��bits_per_item��������а�λ�����㣬�õ���ǩ
	  tag = hv & ((1ULL << bits_per_item) - 1);
	  // �����ǩΪ0������1�Ա���ʹ��0��ǩ
	  tag += (tag == 0);
	  return tag;
  }

  // ����������������������ͱ�ǩ��ϣ
  inline void GenerateIndexTagHash(const ItemType& item, size_t* index, uint32_t* tag) const {
	  // ʹ�ù�ϣ��������64λ��ϣֵ
	  const uint64_t hash = hasher_(item);
	  // ��32λ���ڼ�������
	  *index = IndexHash(hash >> 32);
	  // ��32λ���ڼ����ǩ
	  *tag = TagHash(hash);
  }


  // ���������������������
  inline size_t AltIndex(const size_t index, const uint32_t tag) const {
	  // ʹ��MurmurHash2�Ĺ�ϣ����0x5bd1e995���п��ٹ�ϣ
	  return IndexHash((uint32_t)(index ^ (tag * 0x5bd1e995)));
  }

  // �������ڲ�ʵ�֣����������ͱ�ǩ
  Status AddImpl(const size_t i, const uint32_t tag);

  // ����װ�����ӣ���ʾռ����
  double LoadFactor() const {
	  return 1.0 * Size() / table_->SizeInTags();
  }

  // ����ÿ���ƽ��λ��
  double BitsPerItem() const {
	  return 8.0 * table_->SizeInBytes() / Size();
  }

 public:

  // ���캯����������������Ϊ����
  explicit CuckooFilter(const size_t max_num_keys) : num_items_(0), victim_(), hasher_() {
	  size_t assoc = 4; // ÿ��Ͱ�Ĺ�����
	  size_t num_buckets = upperpower2(std::max<uint64_t>(1, max_num_keys / assoc)); // ����Ͱ������
	  double frac = (double)max_num_keys / num_buckets / assoc; // ����ÿ��Ͱ��װ������
	  if (frac > 0.96) { // ���װ�����ӳ���0.96��Ͱ�������ӱ�
		  num_buckets <<= 1;
	  }
	  victim_.used = false; // ��ʼ���ܺ��߻���Ϊδʹ��
	  table_ = new TableType<bits_per_item>(num_buckets); // ��ʼ����
  }

  // �����������ͷű���ڴ�
  ~CuckooFilter() { delete table_; }

  // ����������һ����
  Status Add(const ItemType& item);

  // ���һ�����Ƿ���ڣ������м������ʣ�
  Status Contain(const ItemType& item) const;

  // �ӹ�������ɾ��һ����
  Status Delete(const ItemType& item);

  /* �ṩͳ����Ϣ�ķ��� */
  // ����ժҪ��Ϣ
  std::string Info() const;

  // ���ص�ǰ���������
  size_t Size() const { return num_items_; }

  // ���ع��������ֽڴ�С
  size_t SizeInBytes() const { return table_->SizeInBytes(); }

  // ���л����������ڳ־û�������
  template <class Archive>
  void serialize(Archive& ar) {
	  ar(num_items_, victim_.index, victim_.tag, victim_.used);
	  ar(*table_);
      ar(hasher_);
  }
};

// ģ����CuckooFilter�ĳ�Ա����Add��ʵ��
template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::Add(
    const ItemType &item) {
  size_t i;
  uint32_t tag;

  // ����ܺ��߻����ѱ�ʹ�ã����ؿռ䲻��״̬
  if (victim_.used) {
	  return NotEnoughSpace;
  }

  // ������������ͱ�ǩ��ϣֵ
  GenerateIndexTagHash(item, &i, &tag);
  // �����ڲ�ʵ�ַ���AddImpl
  return AddImpl(i, tag);
}

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::AddImpl(
    const size_t i, const uint32_t tag) {
  size_t curindex = i; // ��ǰ����
  uint32_t curtag = tag; // ��ǰ��ǩ
  uint32_t oldtag;

  // �������kMaxCuckooCount��
  for (uint32_t count = 0; count < kMaxCuckooCount; count++) {
	  bool kickout = count > 0; // ������ǵ�һ�γ��ԣ��������߳�
	  oldtag = 0;
	  // ���Խ���ǩ���뵽��ǰͰ��
	  if (table_->InsertTagToBucket(curindex, curtag, kickout, oldtag)) {
		  num_items_++; // �ɹ��������������
		  return Ok; // ���ز����ɹ�״̬
	  }
	  // ����߳���һ����ǩ��������Ϊ��ǰ��ǩ
	  if (kickout) {
		  curtag = oldtag;
	  }
	  // �����������
	  curindex = AltIndex(curindex, curtag);
  }

  // ������г��Ծ�ʧ�ܣ�����ǰ�����ͱ�ǩ�����ܺ��߻���
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

  // ���ɼ��������ͱ�ǩ��ϣֵ
  GenerateIndexTagHash(key, &i1, &tag);
  // �����������
  i2 = AltIndex(i1, tag);

  // ȷ���������������ȷ
  assert(i1 == AltIndex(i2, tag));

  // ����ܺ��߻������Ƿ��������
  found = victim_.used && (tag == victim_.tag) &&
	  (i1 == victim_.index || i2 == victim_.index);

  // ������ܺ��߻������ҵ�������Ͱ���ҵ�������Ok״̬
  if (found || table_->FindTagInBuckets(i1, i2, tag)) {
	  return Ok;
  }
  else {
	  return NotFound; // ���򷵻�δ�ҵ�״̬
  }
}

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
Status CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::Delete(
    const ItemType &key) {
  size_t i1, i2;
  uint32_t tag;
  // ���ɼ��������ͱ�ǩ��ϣֵ
  GenerateIndexTagHash(key, &i1, &tag);
  // �����������
  i2 = AltIndex(i1, tag);

  // ���Դӵ�һ��Ͱ��ɾ����ǩ
  if (table_->DeleteTagFromBucket(i1, tag)) {
	  num_items_--; // ɾ���ɹ�����������
	  goto TryEliminateVictim;
  }
  // ���Դӵڶ���Ͱ��ɾ����ǩ
  else if (table_->DeleteTagFromBucket(i2, tag)) {
	  num_items_--; // ɾ���ɹ�����������
	  goto TryEliminateVictim;
  }
  // ����ܺ��߻������Ƿ��������
  else if (victim_.used && tag == victim_.tag &&
	  (i1 == victim_.index || i2 == victim_.index)) {
	  victim_.used = false; // ɾ���ɹ�������ܺ��߻���δʹ��
	  return Ok;
  }
  // �������ɾ��������δ�ɹ�������δ�ҵ�״̬
  else {
	  return NotFound;
  }

TryEliminateVictim:
  // ����ܺ��߻�������δ���������Խ������
  if (victim_.used) {
	  victim_.used = false;
	  size_t i = victim_.index;
	  uint32_t tag = victim_.tag;
	  AddImpl(i, tag); // �������²����ܺ��߻����е���
  }
  return Ok; // ���ز����ɹ�״̬
}

template <typename ItemType, size_t bits_per_item,
          template <size_t> class TableType, typename HashFamily>
std::string CuckooFilter<ItemType, bits_per_item, TableType, HashFamily>::Info() const {
  std::stringstream ss;
  // ����������������״̬��Ϣ
  ss << "CuckooFilter Status:\n"
	  << "\t\t" << table_->Info() << "\n" // ��������Ϣ
	  << "\t\tKeys stored: " << Size() << "\n" // ����洢�ļ���
	  << "\t\tLoad factor: " << LoadFactor() << "\n" // �����������
	  << "\t\tHashtable size: " << (table_->SizeInBytes() >> 10) << " KB\n"; // �����ϣ��Ĵ�С����KBΪ��λ��
  if (Size() > 0) {
	  ss << "\t\tbit/key:   " << BitsPerItem() << "\n"; // ����д洢����ÿ������ƽ��λ��
  }
  else {
	  ss << "\t\tbit/key:   N/A\n"; // ���û�д洢����N/A
  }
  return ss.str();
}
}  // namespace cuckoofilter
#endif  // CUCKOO_FILTER_CUCKOO_FILTER_H_
