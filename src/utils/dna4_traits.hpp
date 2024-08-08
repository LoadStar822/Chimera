// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/io/sequence_file/input.hpp>

namespace raptor
{
/**
 * \brief 定义一个新的 traits 结构体，以便在 SeqAn3 的 sequence_file_input 中使用。
 *
 * 这个结构体继承自 seqan3::sequence_file_input_default_traits_dna，并将序列字母表设置为 seqan3::dna4。
 * 这意味着输入的序列将被视为 DNA 序列，并且只包含 A, C, G, T 这四种碱基。
 */
struct dna4_traits : seqan3::sequence_file_input_default_traits_dna
{
    // 定义序列字母表类型为 seqan3::dna4，只包含 A, C, G, T 四种碱基。
    using sequence_alphabet = seqan3::dna4;
};

} // namespace raptor
