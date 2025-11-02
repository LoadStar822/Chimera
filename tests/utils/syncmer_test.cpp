#include <array>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "utils/Syncmer.hpp"
#include "utils/FeatureHasher.hpp"

namespace {

using seqan3::operator""_dna4;

seqan3::dna4_vector make_sequence(std::string_view pattern, size_t repeats)
{
    seqan3::dna4_vector result;
    result.reserve(pattern.size() * repeats);
    for (size_t i = 0; i < repeats; ++i)
    {
        for (char c : pattern)
        {
            result.push_back(seqan3::assign_char_to(c, seqan3::dna4{}));
        }
    }
    return result;
}

bool expect_vector_equal(std::string const & name,
                         std::vector<uint64_t> const & got,
                         std::vector<uint64_t> const & expected,
                         std::string & message)
{
    if (got == expected)
        return true;

    message = name + " 结果不一致\n  实际: ";
    for (auto v : got)
        message += std::to_string(v) + " ";
    message += "\n  期望: ";
    for (auto v : expected)
        message += std::to_string(v) + " ";
    message += '\n';
    return false;
}

bool expect_throw(std::string const & name, std::function<void()> const & fn, std::string & message)
{
    try
    {
        fn();
    }
    catch (std::invalid_argument const &)
    {
        return true;
    }
    catch (...)
    {
        message = name + " 抛出了非 std::invalid_argument 异常";
        return false;
    }

    message = name + " 未抛出 std::invalid_argument";
    return false;
}

} // namespace

int main()
{
    int failures = 0;
    std::vector<std::string> failure_messages;

    const seqan3::dna4_vector sequence4{"ACGGCGACGTTTAG"_dna4};
    const seqan3::dna5_vector sequence(sequence4.begin(), sequence4.end());

    {
        const std::vector<uint64_t> expected{186, 196, 39, 7, 1870};
        auto result = chimera::syncmer::compute_hashes(sequence, 2, 5, std::array<size_t, 1>{0}, 0, true);
        std::string message;
        if (!expect_vector_equal("open_syncmer", result, expected, message))
        {
            ++failures;
            failure_messages.push_back(std::move(message));
        }
    }

    {
        const std::vector<uint64_t> expected{186, 1426, 196, 39, 7, 1870};
        auto result = chimera::syncmer::compute_hashes(sequence, 2, 5, std::array<size_t, 2>{0, 3}, 0, true);
        std::string message;
        if (!expect_vector_equal("closed_syncmer", result, expected, message))
        {
            ++failures;
            failure_messages.push_back(std::move(message));
        }
    }

    {
        std::string message;
        if (!expect_throw("invalid_positions", [&] {
                chimera::syncmer::compute_hashes(sequence, 2, 5, std::array<size_t, 1>{4});
            }, message))
        {
            ++failures;
            failure_messages.push_back(std::move(message));
        }
    }

    {
        auto params = chimera::feature::auto_params_from_readlen(60);
        if (params.method != chimera::feature::Method::Syncmer)
        {
            ++failures;
            failure_messages.emplace_back("auto_params 在短读下应回退到 syncmer");
        }
    }

    {
        auto params = chimera::feature::auto_params_from_readlen(600);
        chimera::feature::Method expected =
            chimera::feature::strobemer_available() ? chimera::feature::Method::Strobemer
                                                    : chimera::feature::Method::Syncmer;
        if (params.method != expected)
        {
            ++failures;
            failure_messages.emplace_back("auto_params 在长读下应选择正确的方法");
        }
    }

    const auto long_seq = make_sequence("ACGTTGCATGTCAGTCA", 12);   // 192 bp
    const auto short_seq = make_sequence("ACGTACGT", 5);       // 40 bp

    {
        chimera::feature::Params sync_params{};
        sync_params.method = chimera::feature::Method::Syncmer;
        sync_params.sync = {.k = 7, .s = 4, .pos = 1, .seed = 0, .canonical = true};
        auto canonical_hashes = chimera::feature::compute_hashes(long_seq, sync_params);
        sync_params.sync.canonical = false;
        auto non_canonical_hashes = chimera::feature::compute_hashes(long_seq, sync_params);
        if (canonical_hashes.empty() || non_canonical_hashes.empty())
        {
            ++failures;
            failure_messages.emplace_back("syncmer canonical/non-canonical 结果不应为空");
        }
    }

    {
        auto sync_only = chimera::feature::Params{.method = chimera::feature::Method::Syncmer,
                                                  .sync = {.k = 31, .s = 16, .pos = 7, .seed = 0, .canonical = true},
                                                  .strobe = {}};

        chimera::feature::Params auto_params_short{};
        auto_params_short.method = chimera::feature::Method::Auto;
        auto_params_short.sync = sync_only.sync;
        auto_params_short.strobe = {.k = 28, .order = 3, .w_min = 64, .w_max = 80, .seed = 0, .canonical = true};

        auto sync_short_hashes = chimera::feature::compute_hashes(short_seq, sync_only);
        auto auto_short_hashes = chimera::feature::compute_hashes(short_seq, auto_params_short);
        if (sync_short_hashes != auto_short_hashes)
        {
            ++failures;
            failure_messages.emplace_back("AUTO 在短读上应与 syncmer 结果一致");
        }

        chimera::feature::Params auto_params_long{};
        auto_params_long.method = chimera::feature::Method::Auto;
        auto_params_long.sync = sync_only.sync;
        auto_params_long.strobe = {.k = 28, .order = 2, .w_min = 12, .w_max = 32, .seed = 0, .canonical = true};
        auto auto_long_hashes = chimera::feature::compute_hashes(long_seq, auto_params_long);

        if (chimera::feature::strobemer_available())
        {
            chimera::feature::Params strobe_params{};
            strobe_params.method = chimera::feature::Method::Strobemer;
            strobe_params.strobe = {.k = 28, .order = 2, .w_min = 12, .w_max = 32, .seed = 0, .canonical = true};
            auto strobe_hashes = chimera::feature::compute_hashes(long_seq, strobe_params);
            if (strobe_hashes.empty())
            {
                ++failures;
                failure_messages.emplace_back("strobemer 哈希结果不应为空");
            }
            if (auto_long_hashes != strobe_hashes)
            {
                ++failures;
                failure_messages.emplace_back("AUTO 在长读上应与 strobemer 结果一致");
            }
        }
        else
        {
            if (auto_long_hashes.empty())
            {
                ++failures;
                failure_messages.emplace_back("AUTO 在无 strobemer 支持时仍应产出 syncmer 哈希");
            }
        }
    }

    if (failures == 0)
    {
        std::cout << "syncmer tests: 全部通过" << std::endl;
        return 0;
    }

    std::cout << "syncmer tests: 失败 " << failures << " 项" << std::endl;
    for (auto const & msg : failure_messages)
    {
        std::cout << " - " << msg << std::endl;
    }
    return 1;
}
