#pragma once

#include "buildConfig.hpp"

#include <cstdint>
#include <string_view>
#include <vector>

namespace ChimeraBuild
{

	class BuildTaxonomy
	{
	  public:
		static BuildTaxonomy load_required(const BuildConfig &config, std::string_view purpose);

		uint32_t to_species(uint32_t taxid) const noexcept;
		uint32_t to_genus(uint32_t taxid) const noexcept;

	  private:
		std::vector<uint32_t> parent_;
		std::vector<uint32_t> merged_;
		std::vector<uint8_t> is_species_;
		std::vector<uint8_t> is_genus_;
	};

} // namespace ChimeraBuild
