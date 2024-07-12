#pragma once

namespace chimera
{
	struct Config
	{
		uint8_t kmer_size = 19;
		uint16_t window_size = 31;
		uint64_t min_length = 0;
		kstring_t tmp_output_folder = { 0, 0, NULL };
		kstring_t input_file = { 0, 0, NULL };
		kstring_t output_file = { 0, 0, NULL };
		Config() {
			kputs("/tmp", &tmp_output_folder);
		}
	};
}