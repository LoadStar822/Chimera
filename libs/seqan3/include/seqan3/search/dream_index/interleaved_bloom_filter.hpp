// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan3::interleaved_bloom_filter.
 */

#pragma once

#include <algorithm>
#include <bit>

#include <sdsl/bit_vectors.hpp>

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/detail/strong_type.hpp>

namespace seqan3
{
//!\brief Determines if the Interleaved Bloom Filter is compressed.
//!\ingroup search_dream_index
enum data_layout : bool
{
	uncompressed, //!< 交错布隆过滤器未压缩。
	compressed    //!< 交错布隆过滤器已压缩。
};

//!\brief 表示seqan3::interleaved_bloom_filter的bin数量的强类型。
//!\ingroup search_dream_index
struct bin_count : public detail::strong_type<size_t, bin_count, detail::strong_type_skill::convert>
{
    using detail::strong_type<size_t, bin_count, detail::strong_type_skill::convert>::strong_type;
};

//!\brief 表示seqan3::interleaved_bloom_filter中每个bin的比特数的强类型。
//!\ingroup search_dream_index
struct bin_size : public detail::strong_type<size_t, bin_size, detail::strong_type_skill::convert>
{
    using detail::strong_type<size_t, bin_size, detail::strong_type_skill::convert>::strong_type;
};

//!\brief 表示seqan3::interleaved_bloom_filter的哈希函数数量的强类型。
//!\ingroup search_dream_index
struct hash_function_count : public detail::strong_type<size_t, hash_function_count, detail::strong_type_skill::convert>
{
    using detail::strong_type<size_t, hash_function_count, detail::strong_type_skill::convert>::strong_type;
};

//!\brief 表示seqan3::interleaved_bloom_filter的bin索引的强类型。
//!\ingroup search_dream_index
struct bin_index : public detail::strong_type<size_t, bin_index, detail::strong_type_skill::convert>
{
    using detail::strong_type<size_t, bin_index, detail::strong_type_skill::convert>::strong_type;
};

/*!\brief IBF分箱目录。一种能够高效回答多个bin集合成员查询的数据结构。
 * \ingroup search_dream_index
 * \tparam data_layout_mode_ 指示底层数据类型是否压缩。参见seqan3::data_layout。
 * \implements seqan3::cerealisable
 *
 * \details
 *
 * ### 分箱目录
 *
 * 分箱目录是一种数据结构，可以用来确定元素的集合成员关系。
 * 例如，一个常见的使用案例是通过某种聚类方法（例如分类分箱或基于k-mer相似性聚类基因组序列）将数据库划分为固定数量（例如1024个）的bin。
 * 对于查询，分箱目录可以回答查询可能出现在哪些bin中。
 * 在SeqAn中，我们提供了交错布隆过滤器（IBF），可以高效地回答这些查询。
 *
 * ### 交错布隆过滤器（IBF）
 *
 * 交错布隆过滤器是一种概率数据结构，扩展了[布隆过滤器](https://en.wikipedia.org/wiki/Bloom_filter)。
 * 布隆过滤器可以被认为是一个长度为`n`的位向量和`h`个哈希函数，用于确定集合成员关系。
 * 插入数据时，数据通过`h`个哈希函数进行哈希（返回值在`[0, n)`范围内），并将位向量中对应的`h`个位置设为`1`。
 * 查询数据时，即确定查询是否属于布隆过滤器构建的集合时，通过相同的`h`个哈希函数对查询进行哈希，并检查对应的位置。
 * 如果所有`h`个位置都为`1`，则查询（可能）在数据集中。
 * 由于布隆过滤器具有可变长度，哈希不是双射的，即使查询从未插入到布隆过滤器中，它也可能返回真。
 * 注意，布隆过滤器总是会返回真，如果查询被插入，则可能有假阳性，但没有假阴性。
 *
 * 交错布隆过滤器现在将布隆过滤器的概念应用于多个集合，并提供一个*全局*数据结构，以确定查询在`b`个数据集/bin中的集合成员关系。
 * 从概念上讲，为每个bin创建一个布隆过滤器，使用相同的固定长度和固定哈希函数。
 * 然后将生成的`b`个布隆过滤器交错排列，使得每个布隆过滤器的第`i`位相邻：
 * ```
 * 布隆过滤器0       布隆过滤器1      布隆过滤器2      布隆过滤器3
 * |0.0|0.1|0.2|0.3|    |1.0|1.1|1.2|1.3|   |2.0|2.1|2.2|2.3|   |3.0|3.1|3.2|3.3|
 * ```
 * 其中`x.y`表示第`x`个布隆过滤器的第`y`位。
 * ```
 * 交错布隆过滤器
 * |0.0|1.0|2.0|3.0|0.1|1.1|2.1|3.1|0.2|1.2|2.2|3.2|0.3|1.3|2.3|3.3|
 * ```
 * 现在可以通过计算`h`个哈希函数，在哈希函数指示的位置检索长度为`b`的`h`个子位向量，在所有`b`个bin中搜索查询。
 * 这些子位向量的按位与操作产生了binning向量，一个长度为`b`的位向量，其中第`i`位表示在第`i`个bin中的集合成员关系。
 *
 * ### 查询
 * 要查询交错布隆过滤器中的值，调用seqan3::interleaved_bloom_filter::membership_agent()并使用返回的seqan3::interleaved_bloom_filter::membership_agent_type。
 *
 * 要计算交错布隆过滤器中一系列值的出现次数，调用seqan3::interleaved_bloom_filter::counting_agent()并使用
 * 返回的seqan3::interleaved_bloom_filter::counting_agent_type。
 *
 * ### 压缩
 *
 * 交错布隆过滤器可以通过传递`data_layout::compressed`作为模板参数进行压缩。
 * 压缩的`seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>`只能从`seqan3::interleaved_bloom_filter`构建，在这种情况下，底层位向量被压缩。
 * 压缩的交错布隆过滤器是不可变的，即只支持查询操作。
 *
 * ### 线程安全
 *
 * 交错布隆过滤器保证STL的基本线程安全，即所有对`const`成员函数的调用在多个线程中是安全的（只要没有线程同时调用非`const`成员函数）。
 *
 * 此外，当每个线程处理的bin数是wordsize（=64）的整数倍时，多个线程并发调用`emplace`是安全的。
 * 例如，如果`thread_1`访问bins 0-63，`thread_2`访问bins 64-127，依此类推，则可以安全地从多个线程调用`emplace`。
 */

template <data_layout data_layout_mode_ = data_layout::uncompressed>
class interleaved_bloom_filter
{
private:
    //!\cond
    template <data_layout data_layout_mode>
    friend class interleaved_bloom_filter;
    //!\endcond

    //!\brief 底层数据类型的定义。
    using data_type =
        std::conditional_t<data_layout_mode_ == data_layout::uncompressed, sdsl::bit_vector, sdsl::sd_vector<>>;

	//!\brief 用户指定的bin数量。
	size_t bins{};
	//!\brief IBF中存储的bin数量（`bins`的下一个64的倍数）。
	size_t technical_bins{};
	//!\brief 每个bin的大小，以位为单位。
	size_t bin_size_{};
	//!\brief 进行乘法哈希前需要移位的哈希值位数。
	size_t hash_shift{};
	//!\brief 存储`bins`个bit所需的64位整数的数量（例如，`bins = 50` -> `bin_words = 1`）。
	size_t bin_words{};
	//!\brief 哈希函数的数量。
	size_t hash_funs{};
	//!\brief 位向量。
	data_type data{};
      //!\brief 预计算的用于乘法哈希的种子。我们使用大数以保证均匀分布的哈希。
    static constexpr std::array<size_t, 5> hash_seeds{13572355802537770549ULL, // 2**64 / (e/2)
                                                      13043817825332782213ULL, // 2**64 / sqrt(2)
                                                      10650232656628343401ULL, // 2**64 / sqrt(3)
                                                      16499269484942379435ULL, // 2**64 / (sqrt(5)/2)
                                                      4893150838803335377ULL}; // 2**64 / (3*pi/5)

    /*!\brief 扰乱一个值并将其适配到向量中。
     * \param h 要处理的值。
     * \param seed 使用的种子。
     * \returns 一个表示在`data`范围内的位置的哈希值。
     * \sa https://probablydance.com/2018/06/16/
     * \sa https://lemire.me/blog/2016/06/27
     * 对一个输入值进行哈希处理，并将处理后的值映射到交错布隆过滤器的位向量中的一个位置。
     */
    inline constexpr size_t hash_and_fit(size_t h, size_t const seed) const
    {
        h *= seed;
        assert(hash_shift < 64);
        h ^= h >> hash_shift;         // 将高位通过异或移位到低位
        h *= 11400714819323198485ULL; // 乘以黄金分割比，以扩展h到64位范围
                                      // 如果可能，使用fastrange（不进行除法的整数取模）。
        return h % (technical_bins * bin_size_); // 计算哈希值在data范围内的位置
    }
#ifdef __SIZEOF_INT128__
	    // 如果编译器支持128位整数类型
        h = static_cast<uint64_t>((static_cast<__uint128_t>(h) * static_cast<__uint128_t>(bin_size_)) >> 64);
#else
	    // 如果编译器不支持128位整数类型
        h %= bin_size_;
#endif
        //将哈希值扩大到实际存储的bin数量范围内。这个操作将哈希值映射到整个交错布隆过滤器的位向量范围内。
        h *= technical_bins;
        return h;
    }

public:
    //!\brief 指示交错布隆过滤器是否压缩。
    static constexpr data_layout data_layout_mode = data_layout_mode_;

    class membership_agent_type; // documented upon definition below

    template <std::integral value_t>
    class counting_agent_type; // documented upon definition below

    /*!\name Constructors, destructor and assignment
     * \{
     */
	interleaved_bloom_filter() = default;                                             //!< 默认构造函数
	interleaved_bloom_filter(interleaved_bloom_filter const&) = default;             //!< 默认拷贝构造函数
	interleaved_bloom_filter& operator=(interleaved_bloom_filter const&) = default; //!< 默认拷贝赋值操作符
	interleaved_bloom_filter(interleaved_bloom_filter&&) = default;                  //!< 默认移动构造函数
	interleaved_bloom_filter& operator=(interleaved_bloom_filter&&) = default;      //!< 默认移动赋值操作符
	~interleaved_bloom_filter() = default;                                            //!< 默认析构函数

    /*!\brief 构建一个未压缩的交错布隆过滤器。
     * \param bins_ bin 的数量。
     * \param size 位向量的大小。
     * \param funs 哈希函数的数量。默认值为 2。最少为 1，最多为 5。
     *
     * \attention 该构造函数只能用于构建**未压缩**的交错布隆过滤器。
     *
     * \details
     *
     * ### 示例
     *
     * \include test/snippet/search/dream_index/interleaved_bloom_filter_constructor.cpp
     */

    interleaved_bloom_filter(seqan3::bin_count bins_,
                             seqan3::bin_size size,
                             seqan3::hash_function_count funs = seqan3::hash_function_count{2u})
        requires (data_layout_mode == data_layout::uncompressed)
    {
		// 从传入参数获取实际值
		bins = bins_.get();
		bin_size_ = size.get();
		hash_funs = funs.get();

		// 验证参数的合法性
		if (bins == 0)
			throw std::logic_error{ "The number of bins must be > 0." }; // bin 数量必须大于 0
		if (hash_funs == 0 || hash_funs > 5)
			throw std::logic_error{ "The number of hash functions must be > 0 and <= 5." }; // 哈希函数的数量必须在 1 到 5 之间
		if (bin_size_ == 0)
			throw std::logic_error{ "The size of a bin must be > 0." }; // 每个 bin 的大小必须大于 0

		// 计算哈希移位值
		hash_shift = std::countl_zero(bin_size_); // 计算 bin_size_ 的前导零位数，用于后续哈希计算中的移位操作

		// 计算技术参数
		bin_words = (bins + 63) >> 6;    // 计算存储 bins 个 bin 所需的 64 位整数的数量，向上取整
		technical_bins = bin_words << 6; // 确保技术上的 bin 数量是 64 的倍数

		// 初始化位向量，大小为 technical_bins * bin_size_
		data = sdsl::bit_vector(technical_bins * bin_size_);
    }

    /*!\brief 从压缩的交错布隆过滤器构建一个未压缩的交错布隆过滤器。
     * \param[in] ibf 压缩的 seqan3::interleaved_bloom_filter。
     * \details
     *
     * ### 示例
     *
     * \include test/snippet/search/dream_index/interleaved_bloom_filter_constructor_uncompress.cpp
     */
    interleaved_bloom_filter(interleaved_bloom_filter<data_layout::compressed> const& ibf)
        requires (data_layout_mode == data_layout::uncompressed)
    {
        // 解压缩时，复制所有技术参数
        std::tie(bins, technical_bins, bin_size_, hash_shift, bin_words, hash_funs) =
            std::tie(ibf.bins, ibf.technical_bins, ibf.bin_size_, ibf.hash_shift, ibf.bin_words, ibf.hash_funs);

        // 将压缩的位向量数据解压到未压缩的位向量数据
        data = sdsl::bit_vector{ ibf.data.begin(), ibf.data.end() };
    }

    /*!\brief 构建一个压缩的交错布隆过滤器。
     * \param[in] ibf 未压缩的 seqan3::interleaved_bloom_filter。
     *
     * \attention 该构造函数只能用于构建**压缩**的交错布隆过滤器。
     *
     * \details
     *
     * ### 示例
     *
     * \include test/snippet/search/dream_index/interleaved_bloom_filter_constructor_compressed.cpp
     */
    interleaved_bloom_filter(interleaved_bloom_filter<data_layout::uncompressed> const& ibf)
        requires (data_layout_mode == data_layout::compressed)
    {
        // 压缩时，复制所有技术参数
        std::tie(bins, technical_bins, bin_size_, hash_shift, bin_words, hash_funs) =
            std::tie(ibf.bins, ibf.technical_bins, ibf.bin_size_, ibf.hash_shift, ibf.bin_words, ibf.hash_funs);

        // 将未压缩的位向量数据压缩到压缩的位向量数据
        data = sdsl::sd_vector<>{ ibf.data };
    }
    //!\}

    /*!\brief 将一个值插入到特定的 bin 中。
	 * \param[in] value 要处理的原始数值。
	 * \param[in] bin 要插入的 bin 索引。
	 *
	 * \attention 该函数仅适用于**未压缩**的交错布隆过滤器。
	 *
	 * \details
	 *
	 * ### 示例
	 *
	 * \include test/snippet/search/dream_index/interleaved_bloom_filter_emplace.cpp
	 */
    void emplace(size_t const value, bin_index const bin) noexcept
        requires (data_layout_mode == data_layout::uncompressed)
    {
        // 确保 bin 索引在有效范围内
        assert(bin.get() < bins);
        // 对于每个哈希函数
        for (size_t i = 0; i < hash_funs; ++i)
        {
            // 使用哈希函数计算值的索引
            size_t idx = hash_and_fit(value, hash_seeds[i]);
            // 加上 bin 的偏移量
 
            idx += bin.get();
            // 确保计算的索引在数据范围内
            assert(idx < data.size());
            // 将对应位置的位设为 1
            data[idx] = 1;
        }
    }


    /*!\brief 清除特定的 bin。
	 * \param[in] bin 要清除的 bin 索引。
	 *
	 * \attention 该函数仅适用于**未压缩**的交错布隆过滤器。
	 *
	 * \details
	 *
	 * ### 示例
	 *
	 * \include test/snippet/search/dream_index/interleaved_bloom_filter_clear.cpp
	 */
    void clear(bin_index const bin) noexcept
        requires (data_layout_mode == data_layout::uncompressed)
    {
        // 确保 bin 索引在有效范围内
        assert(bin.get() < bins);
        // 清除 bin 中的所有位
        for (size_t idx = bin.get(), i = 0; i < bin_size_; idx += technical_bins, ++i)
            data[idx] = 0;
    }

    /*!\brief 清除一组 bin 的范围。
	 * \tparam rng_t 范围的类型。必须符合 std::ranges::forward_range，并且引用类型必须是 seqan3::bin_index。
	 * \param[in] bin_range 要清除的 bin 范围。
	 *
	 * \attention 该函数仅适用于**未压缩**的交错布隆过滤器。
	 *
	 * \details
	 *
	 * ### 示例
	 *
	 * \include test/snippet/search/dream_index/interleaved_bloom_filter_clear.cpp
	 */
    template <typename rng_t>
        requires (data_layout_mode == data_layout::uncompressed)
    void clear(rng_t&& bin_range) noexcept
    {
        // 静态断言，确保传入的范围符合 forward_range 概念，并且其引用类型是 seqan3::bin_index
        static_assert(std::ranges::forward_range<rng_t>, "The range of bins to clear must model a forward_range.");
        static_assert(std::same_as<std::remove_cvref_t<std::ranges::range_reference_t<rng_t>>, bin_index>,
            "The reference type of the range to clear must be seqan3::bin_index.");

        // 在非调试模式下，检查每个 bin 的索引是否在有效范围内
#ifndef NDEBUG
        for (auto&& bin : bin_range)
            assert(bin.get() < bins);
#endif // NDEBUG

        // 清除 bin_range 范围内的每个 bin
        for (size_t offset = 0, i = 0; i < bin_size_; offset += technical_bins, ++i)
            for (auto&& bin : bin_range)
                data[bin.get() + offset] = 0;
    }


    /*!\brief 增加存储在交错布隆过滤器中的 bin 数量。
     * \param[in] new_bins_ 新的 bin 数量。
     * \throws std::invalid_argument 如果传入的 bin 数量小于当前 bin 数量。
     *
     * \attention 该函数仅适用于**未压缩**的交错布隆过滤器。
     * \attention 新的 bin 数量必须大于或等于当前 bin 数量。
     * \attention 该函数会使所有为此交错布隆过滤器构建的 seqan3::interleaved_bloom_filter::membership_agent_type 失效。
     *
     * \details
     *
     * 结果 `seqan3::interleaved_bloom_filter` 的大小与 `bin_words` 的增加成比例增加（`bin_words` 是表示 `bins` 的 64 位字的数量），例如，
     * 将具有 40 个 bin 的 `seqan3::interleaved_bloom_filter` 调整为 73 个 bin 也会将 `bin_words` 从 1 增加到 2，因此新的 `seqan3::interleaved_bloom_filter`
     * 的大小将是原来的两倍。这种大小的增加是必要的，以避免使所有计算的哈希函数失效。如果您希望在保持大小不变的情况下添加更多的 bin，则需要重建
     * `seqan3::interleaved_bloom_filter`。
     *
     * ### 示例
     *
     * \include test/snippet/search/dream_index/interleaved_bloom_filter_increase_bin_number_to.cpp
     */
    void increase_bin_number_to(bin_count const new_bins_)
        requires (data_layout_mode == data_layout::uncompressed)
    {
        size_t new_bins = new_bins_.get();

        // 检查新的 bin 数量是否小于当前 bin 数量
        if (new_bins < bins)
            throw std::invalid_argument{ "The number of new bins must be >= the current number of bins." };

        // 计算新的 bin 所需的 64 位字的数量，相当于 ceil(new_bins / 64)
        size_t new_bin_words = (new_bins + 63) >> 6;

        bins = new_bins;

        // 如果 bin_words 没有变化，则不需要内部调整大小
        if (new_bin_words == bin_words)
            return;

        // 计算新的技术 bin 数量和新的位向量大小
        size_t new_technical_bins = new_bin_words << 6;
        size_t new_bits = bin_size_ * new_technical_bins;

        size_t idx_{ new_bits }, idx{ data.size() };
        size_t delta = new_technical_bins - technical_bins + 64;

        // 调整位向量的大小
        data.resize(new_bits);

        // 重新排列位向量中的数据
        for (size_t i = idx_, j = idx; j > 0; i -= new_technical_bins, j -= technical_bins)
        {
            size_t stop = i - new_technical_bins;

            for (size_t ii = i - delta, jj = j - 64; stop && ii >= stop; ii -= 64, jj -= 64)
            {
                uint64_t old = data.get_int(jj);
                data.set_int(jj, 0);
                data.set_int(ii, old);
            }
        }

        // 更新 bin_words 和 technical_bins 的值
        bin_words = new_bin_words;
        technical_bins = new_technical_bins;
    }
	//!\}

    /*!\name Lookup
     * \{
     */
     /*!\brief 返回一个 seqan3::interleaved_bloom_filter::membership_agent_type，用于查找。
      * \attention 调用 seqan3::interleaved_bloom_filter::increase_bin_number_to 会使所有为此交错布隆过滤器构建的
      * `seqan3::interleaved_bloom_filter::membership_agent_type` 失效。
      *
      * \details
      *
      * ### 示例
      *
      * \include test/snippet/search/dream_index/membership_agent_construction.cpp
      * \sa seqan3::interleaved_bloom_filter::membership_agent_type::bulk_contains
      */
    membership_agent_type membership_agent() const
    {
        return membership_agent_type{ *this };
    }


    /*!\brief 返回一个 seqan3::interleaved_bloom_filter::counting_agent_type，用于计数。
     * \attention 调用 seqan3::interleaved_bloom_filter::increase_bin_number_to 会使所有为此交错布隆过滤器构建的
     * `seqan3::interleaved_bloom_filter::counting_agent_type` 失效。
     *
     * \details
     *
     * ### 示例
     *
     * \include test/snippet/search/dream_index/counting_agent_construction.cpp
     * \sa seqan3::interleaved_bloom_filter::counting_agent_type::bulk_count
     */
    template <typename value_t = uint16_t>
    counting_agent_type<value_t> counting_agent() const
    {
        return counting_agent_type<value_t>{*this};
    }
    //!\}

    /*!\brief 返回交错布隆过滤器中使用的哈希函数数量。
     * \returns 哈希函数的数量。
     */
    size_t hash_function_count() const noexcept
    {
        return hash_funs;
    }


    /*!\brief 返回交错布隆过滤器管理的 bin 数量。
     * \returns bin 的数量。
     */
    size_t bin_count() const noexcept
    {
        return bins;
    }


    /*!\brief 返回交错布隆过滤器管理的单个 bin 的大小。
     * \returns 单个 bin 的大小（以位为单位）。
     */
    size_t bin_size() const noexcept
    {
        return bin_size_;
    }

    /*!\brief 返回底层位向量的大小。
     * \returns 底层位向量的大小（以位为单位）。
     */
    size_t bit_size() const noexcept
    {
        return data.size();
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
     /*!\brief 测试两个交错布隆过滤器是否相等。
      * \param[in] lhs 一个 `seqan3::interleaved_bloom_filter` 实例。
      * \param[in] rhs 要比较的另一个 `seqan3::interleaved_bloom_filter` 实例。
      * \returns 如果相等返回 `true`，否则返回 `false`。
      */
    friend bool operator==(interleaved_bloom_filter const& lhs, interleaved_bloom_filter const& rhs) noexcept
    {
        return std::tie(lhs.bins,
            lhs.technical_bins,
            lhs.bin_size_,
            lhs.hash_shift,
            lhs.bin_words,
            lhs.hash_funs,
            lhs.data)
            == std::tie(rhs.bins,
                rhs.technical_bins,
                rhs.bin_size_,
                rhs.hash_shift,
                rhs.bin_words,
                rhs.hash_funs,
                rhs.data);
    }


    /*!\brief 测试两个交错布隆过滤器是否不相等。
     * \param[in] lhs 一个 `seqan3::interleaved_bloom_filter` 实例。
     * \param[in] rhs 要比较的另一个 `seqan3::interleaved_bloom_filter` 实例。
     * \returns 如果不相等返回 `true`，否则返回 `false`。
     */
    friend bool operator!=(interleaved_bloom_filter const& lhs, interleaved_bloom_filter const& rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}

    /*!\name Access
     * \{
     */
     /*!\brief 提供对底层数据结构的直接、不安全的访问。
      * \returns 对 SDSL 位向量的引用。
      *
      * \details
      *
      * \noapi{数据的确切表示是实现定义的。}
      */
    constexpr data_type& raw_data() noexcept
    {
        return data;
    }

    //!\copydoc raw_data()
    constexpr data_type const& raw_data() const noexcept
    {
        return data;
    }
    //!\}

    /*!\cond DEV
     * \brief 序列化支持函数。
     * \tparam archive_t `archive` 的类型；必须满足 seqan3::cereal_archive。
     * \param[in] archive 正在从/向其序列化的存档。
     *
     * \attention 这些函数从不直接调用，详见 \ref serialisation 了解更多细节。
     */
    template <cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t& archive)
    {
        archive(bins);
        archive(technical_bins);
        archive(bin_size_);
        archive(hash_shift);
        archive(bin_words);
        archive(hash_funs);
        archive(data);
    }
    //!\endcond

};

/*!\brief 管理 seqan3::interleaved_bloom_filter 的成员查询。
 * \attention 调用 `seqan3::interleaved_bloom_filter::increase_bin_number_to` 会使 membership_agent 失效。
 *
 * \details
 *
 * ### 示例
 *
 * \include test/snippet/search/dream_index/membership_agent_construction.cpp
 */
template <data_layout data_layout_mode>
class interleaved_bloom_filter<data_layout_mode>::membership_agent_type
{
private:
	//!\brief 增强型 seqan3::interleaved_bloom_filter 的类型。
	using ibf_t = interleaved_bloom_filter<data_layout_mode>;

	//!\brief 增强型 seqan3::interleaved_bloom_filter 的指针。
	ibf_t const* ibf_ptr{ nullptr };

public:
    class binning_bitvector;

    /*!\name Constructors, destructor and assignment
     * \{
     */
	membership_agent_type() = default;                                          //!< 默认构造函数
	membership_agent_type(membership_agent_type const&) = default;             //!< 默认拷贝构造函数
	membership_agent_type& operator=(membership_agent_type const&) = default; //!< 默认拷贝赋值操作符
	membership_agent_type(membership_agent_type&&) = default;                  //!< 默认移动构造函数
	membership_agent_type& operator=(membership_agent_type&&) = default;      //!< 默认移动赋值操作符
	~membership_agent_type() = default;                                         //!< 默认析构函数

    /*!\brief 从一个 seqan3::interleaved_bloom_filter 构造一个 membership_agent_type。
     * \private
     * \param ibf 增强型的 seqan3::interleaved_bloom_filter。
     */ 
    explicit membership_agent_type(ibf_t const & ibf) : ibf_ptr(std::addressof(ibf)), result_buffer(ibf.bin_count())
    {}
    //!\}

	//!\brief 存储 bulk_contains() 的结果。
    binning_bitvector result_buffer;

    /*!\name Lookup
     * \{
     */
     /*!\brief 确定给定值的集合成员资格。
      * \param[in] value 要处理的原始值。
      *
      * \attention 该函数的返回结果必须始终通过引用绑定，例如 `auto &`，以防止复制。
      * \attention 连续调用该函数会使先前返回的引用失效。
      *
      * \details
      *
      * ### 示例
      *
      * \include test/snippet/search/dream_index/membership_agent_bulk_contains.cpp
      *
      * ### 线程安全
      *
      * 并发调用该函数不是线程安全的，请为每个线程创建一个 `seqan3::interleaved_bloom_filter::membership_agent_type`。
      */
    [[nodiscard]] binning_bitvector const & bulk_contains(size_t const value) & noexcept
    {
		// 确保过滤器指针不为空
		assert(ibf_ptr != nullptr);
		// 确保结果缓冲区的大小与 bin 的数量一致
		assert(result_buffer.size() == ibf_ptr->bin_count());

		// 初始化布隆过滤器索引数组，大小为 5（最大哈希函数数量）
		std::array<size_t, 5> bloom_filter_indices;
		// 将哈希种子复制到布隆过滤器索引数组中
		std::memcpy(&bloom_filter_indices, &ibf_ptr->hash_seeds, sizeof(size_t) * ibf_ptr->hash_funs);

		// 计算每个哈希函数的索引
		for (size_t i = 0; i < ibf_ptr->hash_funs; ++i)
			bloom_filter_indices[i] = ibf_ptr->hash_and_fit(value, bloom_filter_indices[i]);

		// 对每个批次进行处理
		for (size_t batch = 0; batch < ibf_ptr->bin_words; ++batch)
		{
			// 初始化一个临时变量 tmp 为全 1
			size_t tmp{ -1ULL };
			// 对每个哈希函数进行按位与操作
			for (size_t i = 0; i < ibf_ptr->hash_funs; ++i)
			{
				// 确保索引在数据范围内
				assert(bloom_filter_indices[i] < ibf_ptr->data.size());
				// 获取当前索引位置的值并与 tmp 进行按位与操作
				tmp &= ibf_ptr->data.get_int(bloom_filter_indices[i]);
				// 将索引加 64，移动到下一个 64 位块
				bloom_filter_indices[i] += 64;
			}
			// 将按位与的结果存储在结果缓冲区中
			result_buffer.data.set_int(batch << 6, tmp);
		}

		// 返回结果缓冲区的引用
		return result_buffer;
    }

	// `bulk_contains` 不能在临时对象上调用，因为返回的引用指向的对象会立即被销毁。
    [[nodiscard]] binning_bitvector const & bulk_contains(size_t const value) && noexcept = delete;
    //!\}
};

//!\brief 表示 `seqan3::interleaved_bloom_filter` 的 `bulk_contains` 调用结果的位向量。
template <data_layout data_layout_mode>
class interleaved_bloom_filter<data_layout_mode>::membership_agent_type::binning_bitvector
{
private:
	//!\brief 底层使用的数据类型。
	using data_type = sdsl::bit_vector;
	//!\brief 位向量。
	data_type data{};

	// 允许 membership_agent_type 访问私有成员
	friend class membership_agent_type;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
	binning_bitvector() = default;                                      //!< 默认构造函数
	binning_bitvector(binning_bitvector const&) = default;             //!< 默认拷贝构造函数
	binning_bitvector& operator=(binning_bitvector const&) = default; //!< 默认拷贝赋值操作符
	binning_bitvector(binning_bitvector&&) = default;                  //!< 默认移动构造函数
	binning_bitvector& operator=(binning_bitvector&&) = default;      //!< 默认移动赋值操作符
	~binning_bitvector() = default;                                     //!< 默认析构函数


	//!\brief 使用给定大小构造位向量。
    explicit binning_bitvector(size_t const size) : data(size)
    {}
    //!\}

	//!\brief 返回元素数量。
    size_t size() const noexcept
    {
        return data.size();
    }

    /*!\name Iterators
     * \{
     */
     //!\brief 返回指向容器第一个元素的迭代器。
    auto begin() noexcept
    {
        return data.begin();
    }

    //!\copydoc begin()
    auto begin() const noexcept
    {
        return data.begin();
    }

    //!\brief 返回指向容器最后一个元素之后的迭代器。
    auto end() noexcept
    {
        return data.end();
    }

    //!\copydoc end()
    auto end() const noexcept
    {
        return data.end();
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
     //!\brief 测试两个 `binning_bitvector` 是否相等。
    friend bool operator==(binning_bitvector const& lhs, binning_bitvector const& rhs) noexcept
    {
        return lhs.data == rhs.data;
    }

    //!\brief 测试两个 `binning_bitvector` 是否不相等。
    friend bool operator!=(binning_bitvector const& lhs, binning_bitvector const& rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}

    /*!\name Access
     * \{
     */
     //!\brief 返回第 i 个元素。
    auto operator[](size_t const i) noexcept
    {
        assert(i < size());
        return data[i];
    }

    //!\copydoc operator[]()
    auto operator[](size_t const i) const noexcept
    {
        assert(i < size());
        return data[i];
    }

    /*!\brief 提供对底层数据结构的直接、不安全的访问。
     * \returns 对 SDSL 位向量的引用。
     *
     * \details
     *
     * \noapi{数据的确切表示是实现定义的。}
     */
    constexpr data_type& raw_data() noexcept
    {
        return data;
    }

    //!\copydoc raw_data()
    constexpr data_type const& raw_data() const noexcept
    {
        return data;
    }
    //!\}
};

/*!\brief 一个数据结构，行为类似于 `std::vector`，可用于整合多次调用
 *        `seqan3::interleaved_bloom_filter::membership_agent_type::bulk_contains` 的结果。
 * \ingroup search_dream_index
 * \tparam value_t 计数的类型。必须满足 `std::integral` 概念。
 *
 * \details
 *
 * 使用 `seqan3::interleaved_bloom_filter::membership_agent_type::bulk_contains` 操作时，一个常见的用例是汇总查询中所有 k-mers 的结果。
 * 这将为每个 bin 提供查询中 k-mers 的数量，这些信息可以用于进一步的过滤或基于 k-mer 计数的丰度估计。
 *
 * `seqan3::counting_vector` 提供了一种简单的方法，通过 `+=` 操作符来汇总各个
 * `seqan3::interleaved_bloom_filter::membership_agent_type::binning_bitvector`。
 *
 * 应根据 `bulk_contains` 的所有调用为特定 bin 返回命中的情况选择 `value_t` 模板参数，以避免溢出。例如，处理短的 Illumina 读数时，
 * `uint8_t` 就足够了，而长读数则需要至少 `uint32_t`。
 *
 * ### 示例
 *
 * \include test/snippet/search/dream_index/counting_vector.cpp
 */
template <std::integral value_t>
class counting_vector : public std::vector<value_t>
{
private:
	//!\brief 基类类型。
    using base_t = std::vector<value_t>;

	//!\brief 判断 binning_bitvector_t 是否为 `seqan3::interleaved_bloom_filter::membership_agent_type::binning_bitvector`。
    template <typename binning_bitvector_t>
    static constexpr bool is_binning_bitvector =
        std::same_as<binning_bitvector_t,
                     interleaved_bloom_filter<data_layout::uncompressed>::membership_agent_type::binning_bitvector>
        || std::same_as<binning_bitvector_t,
                        interleaved_bloom_filter<data_layout::compressed>::membership_agent_type::binning_bitvector>;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
	counting_vector() = default;                                    //!< 默认构造函数。
	counting_vector(counting_vector const&) = default;             //!< 默认拷贝构造函数。
	counting_vector& operator=(counting_vector const&) = default; //!< 默认拷贝赋值操作符。
	counting_vector(counting_vector&&) = default;                  //!< 默认移动构造函数。
	counting_vector& operator=(counting_vector&&) = default;      //!< 默认移动赋值操作符。
	~counting_vector() = default;                                   //!< 默认析构函数。

    using base_t::base_t;
    //!\}

    /*!\brief 对 `seqan3::interleaved_bloom_filter::membership_agent_type::binning_bitvector` 的位进行逐 bin 加法操作。
     * \tparam binning_bitvector_t 右操作数的类型。必须是 `seqan3::interleaved_bloom_filter::membership_agent_type::binning_bitvector`。
     * \param binning_bitvector `seqan3::interleaved_bloom_filter::membership_agent_type::binning_bitvector`。
     * \attention `counting_vector` 必须至少与 `binning_bitvector` 一样大。
     *
     * \details
     *
     * ### 示例
     *
     * \include test/snippet/search/dream_index/counting_vector.cpp
     */
    template <typename binning_bitvector_t>
        requires is_binning_bitvector<binning_bitvector_t>
    counting_vector & operator+=(binning_bitvector_t const & binning_bitvector)
    {
        for_each_set_bin(binning_bitvector,
                         [this](size_t const bin)
                         {
                             ++(*this)[bin];
                         });
        return *this;
    }

    /*!\brief 对 `seqan3::interleaved_bloom_filter::membership_agent_type::binning_bitvector` 的位进行逐 bin 减法操作。
     * \tparam binning_bitvector_t 右操作数的类型。必须是 `seqan3::interleaved_bloom_filter::membership_agent_type::binning_bitvector`。
     * \param binning_bitvector `seqan3::interleaved_bloom_filter::membership_agent_type::binning_bitvector`。
     * \attention `counting_vector` 必须至少与 `binning_bitvector` 一样大。
     */
    template <typename binning_bitvector_t>
        requires is_binning_bitvector<binning_bitvector_t>
    counting_vector & operator-=(binning_bitvector_t const & binning_bitvector)
    {
        for_each_set_bin(binning_bitvector,
                         [this](size_t const bin)
                         {
                             assert((*this)[bin] > 0);
                             --(*this)[bin];
                         });
        return *this;
    }

    /*!\brief 对两个 `seqan3::counting_vector` 进行逐 bin 加法操作。
     * \param rhs 另一个 `seqan3::counting_vector`。
     * \attention `counting_vector` 必须至少与 `rhs` 一样大。
     *
     * \details
     *
     * ### 示例
     *
     * \include test/snippet/search/dream_index/counting_vector.cpp
     */
    counting_vector & operator+=(counting_vector const & rhs)
    {
		assert(this->size() >= rhs.size()); // 确保计数向量的大小不小于 rhs 的大小

		// 使用 std::transform 对两个计数向量进行逐元素加法操作
        std::transform(this->begin(), this->end(), rhs.begin(), this->begin(), std::plus<value_t>());

        return *this;
    }

    /*!\brief 对两个 `seqan3::counting_vector` 进行逐 bin 减法操作。
     * \param rhs 另一个 `seqan3::counting_vector`。
     * \attention `counting_vector` 必须至少与 `rhs` 一样大。
     */
    counting_vector & operator-=(counting_vector const & rhs)
    {
		assert(this->size() >= rhs.size()); // 确保计数向量的大小不小于 rhs 的大小

		// 使用 std::transform 对两个计数向量进行逐元素减法操作，并使用 lambda 函数确保无溢出
        std::transform(this->begin(),
                       this->end(),
                       rhs.begin(),
                       this->begin(),
                       [](auto a, auto b)
                       {
                           assert(a >= b);
                           return a - b;
                       });

        return *this;
    }

private:
    //!\brief Enumerates all bins of a seqan3::interleaved_bloom_filter::membership_agent_type::binning_bitvector.
    template <typename binning_bitvector_t, typename on_bin_fn_t>
    void for_each_set_bin(binning_bitvector_t && binning_bitvector, on_bin_fn_t && on_bin_fn)
    {
        assert(this->size() >= binning_bitvector.size()); // The counting vector may be bigger than what we need.

        // Jump to the next 1 and return the number of places jumped in the bit_sequence
        auto jump_to_next_1bit = [](size_t & x)
        {
            auto const zeros = std::countr_zero(x);
            x >>= zeros; // skip number of zeros
            return zeros;
        };

        // Each iteration can handle 64 bits
        for (size_t bit_pos = 0; bit_pos < binning_bitvector.size(); bit_pos += 64)
        {
            // get 64 bits starting at position `bit_pos`
            size_t bit_sequence = binning_bitvector.raw_data().get_int(bit_pos);

            // process each relative bin inside the bit_sequence
            for (size_t bin = bit_pos; bit_sequence != 0u; ++bin, bit_sequence >>= 1)
            {
                // Jump to the next 1 and
                bin += jump_to_next_1bit(bit_sequence);

                on_bin_fn(bin);
            }
        }
    }
};

/*!\brief Manages counting ranges of values for the seqan3::interleaved_bloom_filter.
 * \attention Calling seqan3::interleaved_bloom_filter::increase_bin_number_to invalidates the counting_agent_type.
 *
 * \details
 *
 * ### Example
 *
 * \include test/snippet/search/dream_index/counting_agent.cpp
 */
template <data_layout data_layout_mode>
template <std::integral value_t>
class interleaved_bloom_filter<data_layout_mode>::counting_agent_type
{
private:
    //!\brief The type of the augmented seqan3::interleaved_bloom_filter.
    using ibf_t = interleaved_bloom_filter<data_layout_mode>;

    //!\brief A pointer to the augmented seqan3::interleaved_bloom_filter.
    ibf_t const * ibf_ptr{nullptr};

    //!\brief Store a seqan3::interleaved_bloom_filter::membership_agent to call `bulk_contains`.
    membership_agent_type membership_agent;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_agent_type() = default;                                        //!< Defaulted.
    counting_agent_type(counting_agent_type const &) = default;             //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type(counting_agent_type &&) = default;                  //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type &&) = default;      //!< Defaulted.
    ~counting_agent_type() = default;                                       //!< Defaulted.

    /*!\brief Construct a counting_agent_type for an existing seqan3::interleaved_bloom_filter.
     * \private
     * \param ibf The seqan3::interleaved_bloom_filter.
     */
    explicit counting_agent_type(ibf_t const & ibf) :
        ibf_ptr(std::addressof(ibf)),
        membership_agent(ibf),
        result_buffer(ibf.bin_count())
    {}
    //!\}

    //!\brief Stores the result of bulk_count().
    counting_vector<value_t> result_buffer;

    /*!\name Counting
     * \{
     */
     /*!\brief 统计一个范围内每个 bin 的出现次数。
      * \tparam value_range_t 值范围的类型。必须符合 `std::ranges::input_range` 概念。引用类型必须符合 `std::unsigned_integral`。
      * \param[in] values 要处理的值范围。
      *
      * \attention 该函数的返回结果必须始终通过引用绑定，例如 `auto &`，以防止复制。
      * \attention 连续调用该函数会使先前返回的引用失效。
      *
      * \details
      *
      * ### 示例
      *
      * \include test/snippet/search/dream_index/counting_agent.cpp
      *
      * ### 线程安全
      *
      * 并发调用该函数不是线程安全的，请为每个线程创建一个 `seqan3::interleaved_bloom_filter::counting_agent_type` 实例。
      */
	template <std::ranges::range value_range_t>
	[[nodiscard]] counting_vector<value_t> const& bulk_count(value_range_t&& values) & noexcept
	{
		assert(ibf_ptr != nullptr);
		assert(result_buffer.size() == ibf_ptr->bin_count());

		// 静态断言，确保 values 符合 input_range 概念，且其中的每个值都是无符号整数类型
		static_assert(std::ranges::input_range<value_range_t>, "The values must model input_range.");
		static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
			"An individual value must be an unsigned integral.");

		// 将 result_buffer 中的所有值初始化为 0
		std::ranges::fill(result_buffer, 0);

		// 遍历所有值，并调用 bulk_contains 进行查询和累加结果
		for (auto&& value : values)
			result_buffer += membership_agent.bulk_contains(value);

		// 返回 result_buffer 的常量引用
		return result_buffer;
	}

    // `bulk_count` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <std::ranges::range value_range_t>
    [[nodiscard]] counting_vector<value_t> const & bulk_count(value_range_t && values) && noexcept = delete;
    //!\}
};

} // namespace seqan3
