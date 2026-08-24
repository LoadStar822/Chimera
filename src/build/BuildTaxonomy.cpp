#include "BuildTaxonomy.hpp"

#include <cctype>
#include <charconv>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ChimeraBuild
{
	namespace
	{

		uint32_t parse_taxid(std::string_view text) noexcept
		{
			while (!text.empty() && std::isspace(static_cast<unsigned char>(text.front())))
			{
				text.remove_prefix(1);
			}
			while (!text.empty() && std::isspace(static_cast<unsigned char>(text.back())))
			{
				text.remove_suffix(1);
			}
			uint32_t value = 0;
			const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
			return error == std::errc{} && end == text.data() + text.size() ? value : 0;
		}

	} // namespace

	BuildTaxonomy BuildTaxonomy::load_required(const BuildConfig &config, std::string_view purpose)
	{
		std::vector<std::filesystem::path> roots;
		if (!config.taxonomy_dir.empty())
		{
			roots.emplace_back(config.taxonomy_dir);
		}
		if (const char *environment = std::getenv("CHIMERA_NCBI_TAXDUMP_DIR"))
		{
			if (*environment != '\0')
			{
				roots.emplace_back(environment);
			}
		}
		std::filesystem::path nodes;
		for (const auto &root : roots)
		{
			const auto candidate = root / "nodes.dmp";
			if (std::filesystem::exists(candidate))
			{
				nodes = candidate;
				break;
			}
		}
		if (nodes.empty())
		{
			throw std::runtime_error(std::string(purpose) + " requires taxonomy nodes.dmp; provide "
			                                                "--taxonomy-dir or CHIMERA_NCBI_TAXDUMP_DIR");
		}
		std::ifstream input(nodes);
		if (!input)
		{
			throw std::runtime_error("failed to open taxonomy nodes.dmp: " + nodes.string());
		}

		BuildTaxonomy taxonomy;
		std::string line;
		while (std::getline(input, line))
		{
			const size_t first = line.find('|');
			const size_t second = first == std::string::npos ? first : line.find('|', first + 1u);
			const size_t third = second == std::string::npos ? second : line.find('|', second + 1u);
			if (first == std::string::npos || second == std::string::npos || third == std::string::npos)
			{
				continue;
			}
			const uint32_t taxid = parse_taxid(std::string_view(line).substr(0, first));
			const uint32_t parent = parse_taxid(std::string_view(line).substr(first + 1u, second - first - 1u));
			std::string_view rank = std::string_view(line).substr(second + 1u, third - second - 1u);
			while (!rank.empty() && std::isspace(static_cast<unsigned char>(rank.front())))
			{
				rank.remove_prefix(1);
			}
			while (!rank.empty() && std::isspace(static_cast<unsigned char>(rank.back())))
			{
				rank.remove_suffix(1);
			}
			if (taxid == 0)
			{
				continue;
			}
			if (taxid >= taxonomy.parent_.size())
			{
				taxonomy.parent_.resize(static_cast<size_t>(taxid) + 1u, 0);
				taxonomy.is_species_.resize(static_cast<size_t>(taxid) + 1u, 0);
				taxonomy.is_genus_.resize(static_cast<size_t>(taxid) + 1u, 0);
			}
			taxonomy.parent_[taxid] = parent;
			taxonomy.is_species_[taxid] = rank == "species";
			taxonomy.is_genus_[taxid] = rank == "genus";
		}
		const auto merged_path = nodes.parent_path() / "merged.dmp";
		if (std::filesystem::exists(merged_path))
		{
			std::ifstream merged_input(merged_path);
			if (!merged_input)
			{
				throw std::runtime_error("failed to open taxonomy merged.dmp: " + merged_path.string());
			}
			while (std::getline(merged_input, line))
			{
				const size_t first = line.find('|');
				const size_t second = first == std::string::npos ? first : line.find('|', first + 1u);
				if (first == std::string::npos || second == std::string::npos)
				{
					continue;
				}
				const uint32_t old_taxid = parse_taxid(std::string_view(line).substr(0, first));
				const uint32_t new_taxid = parse_taxid(std::string_view(line).substr(first + 1u, second - first - 1u));
				if (old_taxid == 0 || new_taxid == 0)
				{
					continue;
				}
				if (old_taxid >= taxonomy.merged_.size())
				{
					taxonomy.merged_.resize(static_cast<size_t>(old_taxid) + 1u, 0);
				}
				taxonomy.merged_[old_taxid] = new_taxid;
			}
		}
		if (taxonomy.parent_.empty())
		{
			throw std::runtime_error("taxonomy nodes.dmp contains no usable nodes: " + nodes.string());
		}
		return taxonomy;
	}

	uint32_t BuildTaxonomy::to_species(uint32_t taxid) const noexcept
	{
		if (taxid < merged_.size() && merged_[taxid] != 0)
		{
			taxid = merged_[taxid];
		}
		if (taxid == 0 || taxid >= is_species_.size())
		{
			return taxid;
		}
		if (is_species_[taxid])
		{
			return taxid;
		}
		uint32_t current = taxid;
		for (uint32_t steps = 0; steps < 128; ++steps)
		{
			if (current == 0 || current >= parent_.size())
			{
				break;
			}
			const uint32_t parent = parent_[current];
			if (parent == 0 || parent == current)
			{
				break;
			}
			current = parent;
			if (current < is_species_.size() && is_species_[current])
			{
				return current;
			}
		}
		return taxid;
	}

	uint32_t BuildTaxonomy::to_genus(uint32_t taxid) const noexcept
	{
		if (taxid < merged_.size() && merged_[taxid] != 0)
		{
			taxid = merged_[taxid];
		}
		uint32_t current = taxid;
		for (uint32_t steps = 0; steps < 128; ++steps)
		{
			if (current == 0 || current >= parent_.size())
			{
				return 0;
			}
			if (current < is_genus_.size() && is_genus_[current])
			{
				return current;
			}
			const uint32_t parent = parent_[current];
			if (parent == 0 || parent == current)
			{
				return 0;
			}
			current = parent;
		}
		return 0;
	}

} // namespace ChimeraBuild
