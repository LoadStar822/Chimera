/*
 * -----------------------------------------------------------------------------
 * Filename:      khashOperation.hpp
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-07-30
 *
 * Last Modified: 2024-08-06
 *
 * Description:
 *  This is a simple C++ program that demonstrates the usage of a hash table.
 *
 * Version:
 *  1.0
 * -----------------------------------------------------------------------------
 */
#include <iostream>
#include <string>
#include <vector>
#include "khash.h"

 // Define a hash table from string to std::vector<std::string>
KHASH_MAP_INIT_STR(str_vec_str, std::vector<std::string>)

// Insert or update operation: handle std::vector<std::string> and std::string
template <typename T>
void khash_insert(khash_t(str_vec_str)* h, const std::string& key, const T& value) {
	static_assert(std::is_same<T, std::vector<std::string>>::value || std::is_same<T, std::string>::value,
		"Value must be either std::vector<std::string> or std::string");

	int ret;
	// Create a copy of the key
	char* key_copy = strdup(key.c_str()); // Duplicate the string to ensure each key has independent memory space

	// Insert the new key into the hash table
	khiter_t k = kh_put(str_vec_str, h, key_copy, &ret);

	if (ret == -1) { // Insertion failed, free memory
		std::cerr << "Failed to insert key: " << key << std::endl;
		free(key_copy);
		return;
	}

	if (ret == 0) { // Key already exists
		free(key_copy); // Free key_copy if the key already exists, use the existing one
	}

	// Insert or update the value based on the type
	if constexpr (std::is_same<T, std::vector<std::string>>::value) {
		// Insert or update vector
		if (ret != 0) {
			// If it's a new key, insert the entire vector directly
			kh_value(h, k) = value;
			std::cout << "Inserted key: " << key << " with " << value.size() << " values" << std::endl;
		}
		else {
			// If it's an existing key, append the vector content
			std::vector<std::string>& existing_values = kh_value(h, k);
			existing_values.insert(existing_values.end(), value.begin(), value.end());
			std::cout << "Updated key: " << key << " with additional " << value.size() << " values" << std::endl;
		}
	}
	else {
		// Insert or update a single string
		if (ret != 0) {
			// If it's a new key, initialize a vector and insert the single string
			kh_value(h, k) = std::vector<std::string>{ value };
			std::cout << "Inserted key: " << key << " with a single value: " << value << std::endl;
		}
		else {
			// If it's an existing key, append the single string
			std::vector<std::string>& existing_values = kh_value(h, k);
			existing_values.push_back(value);
			std::cout << "Added single value '" << value << "' to key: " << key << std::endl;
		}
	}
}

// Find operation
bool khash_find_vec(khash_t(str_vec_str)* h, const std::string& key, std::vector<std::string>& values) {
	khiter_t k = kh_get(str_vec_str, h, key.c_str()); // Retrieve the key
	if (k != kh_end(h)) {
		values = kh_value(h, k);
		return true;
	}
	return false;
}

// Iterate operation
void khash_iterate_vec(khash_t(str_vec_str)* h) {
	std::cout << "Iterating over hash table:" << std::endl;
	for (khiter_t k = kh_begin(h); k != kh_end(h); ++k) {
		if (kh_exist(h, k)) {
			const char* key = kh_key(h, k);
			const std::vector<std::string>& values = kh_value(h, k);
			std::cout << "Key: " << key << ", Values: [";
			for (const auto& val : values) {
				std::cout << val << " ";
			}
			std::cout << "]" << std::endl;
		}
	}
}

// Delete operation
void khash_delete_vec(khash_t(str_vec_str)* h, const std::string& key) {
	khiter_t k = kh_get(str_vec_str, h, key.c_str()); // Find the key to delete
	if (k != kh_end(h)) {
		free((char*)kh_key(h, k)); // Free the memory of the key
		kh_del(str_vec_str, h, k); // Delete the key
		std::cout << "Deleted key: " << key << std::endl;
	}
	else {
		std::cout << "Key not found, cannot delete: " << key << std::endl;
	}
}