#include <array>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "utils/Syncmer.hpp"

namespace {

using seqan3::operator""_dna4;

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
