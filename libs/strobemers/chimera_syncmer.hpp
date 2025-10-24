/*
 * Chimera-specific syncmer hashing built on top of the strobemers sliding-window
 * primitives.
 *
 * This implementation mirrors the semantics of SeqAn3's syncmer_hash view, while
 * reusing the lightweight rolling hash logic from strobemers. It supports
 * configurable syncmer positions, optional canonical mode, and deterministic
 * seeding.
 */

#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace strobemers::chimera {

struct SyncmerConfig {
    size_t k;
    size_t s;
    std::vector<size_t> positions;
    uint64_t seed;
    bool canonical;
};

class SyncmerHasher {
public:
    explicit SyncmerHasher(SyncmerConfig config)
        : config_(std::move(config))
    {
        if (config_.k == 0 || config_.s == 0) {
            throw std::invalid_argument("syncmer k-mer/s-mer length must be positive");
        }
        if (config_.s > config_.k) {
            throw std::invalid_argument("syncmer s-mer length must not exceed k-mer length");
        }
        if (config_.k > 32 || config_.s > 32) {
            throw std::invalid_argument("syncmer k-mer/s-mer length must be <= 32 for 64-bit packing");
        }
        window_size_ = config_.k - config_.s + 1;
        if (config_.positions.empty()) {
            throw std::invalid_argument("syncmer positions must not be empty");
        }
        for (size_t pos : config_.positions) {
            if (pos >= window_size_) {
                throw std::invalid_argument("syncmer position out of window range");
            }
        }
        kmask_ = mask_for(config_.k);
        smask_ = mask_for(config_.s);
        kshift_ = static_cast<uint32_t>((config_.k - 1) * 2);
        sshift_ = static_cast<uint32_t>((config_.s - 1) * 2);
        reset();
    }

    void reset() {
        forward_k_ = 0;
        reverse_k_ = 0;
        forward_s_ = 0;
        reverse_s_ = 0;
        symbols_processed_ = 0;
        s_index_ = 0;
        forward_queue_.clear();
        reverse_queue_.clear();
    }

    template <typename InputIt, typename Sentinel, typename RankFn, typename Callback>
    void process(InputIt first, Sentinel last, RankFn rank_of, Callback&& callback) {
        reset();
        for (; first != last; ++first) {
            uint8_t rank = static_cast<uint8_t>(rank_of(*first));
            if (rank > 3) {
                reset();
                continue;
            }
            consume_rank(rank, callback);
        }
    }

private:
    struct MinQueue {
        void push(uint64_t value, size_t index) {
            while (!data_.empty() && value < data_.back().first) {
                data_.pop_back();
            }
            data_.emplace_back(value, index);
        }

        void pop_older(size_t min_index) {
            while (!data_.empty() && data_.front().second < min_index) {
                data_.pop_front();
            }
        }

        [[nodiscard]] bool empty() const {
            return data_.empty();
        }

        [[nodiscard]] uint64_t min_value() const {
            return data_.front().first;
        }

        [[nodiscard]] size_t min_index() const {
            return data_.front().second;
        }

        void clear() {
            data_.clear();
        }

    private:
        std::deque<std::pair<uint64_t, size_t>> data_;
    };

    static uint64_t mask_for(size_t length) {
        if (length == 0) {
            return 0;
        }
        if (length >= 32) {
            return std::numeric_limits<uint64_t>::max();
        }
        return (static_cast<uint64_t>(1) << (length * 2)) - 1;
    }

    [[nodiscard]] bool matches_position(size_t window_start, size_t absolute_index, bool reverse) const {
        if (absolute_index < window_start) {
            return false;
        }
        size_t offset = absolute_index - window_start;
        if (offset >= window_size_) {
            return false;
        }
        if (reverse) {
            offset = window_size_ - 1 - offset;
        }
        for (size_t pos : config_.positions) {
            if (pos == offset) {
                return true;
            }
        }
        return false;
    }

    template <typename Callback>
    void consume_rank(uint8_t rank, Callback& callback) {
        forward_k_ = ((forward_k_ << 2) | rank) & kmask_;
        forward_s_ = ((forward_s_ << 2) | rank) & smask_;

        if (config_.canonical) {
            uint8_t rc_rank = static_cast<uint8_t>(rank ^ 0x3u);
            reverse_k_ = ((reverse_k_ >> 2) | (static_cast<uint64_t>(rc_rank) << kshift_)) & kmask_;
            reverse_s_ = ((reverse_s_ >> 2) | (static_cast<uint64_t>(rc_rank) << sshift_)) & smask_;
        }

        ++symbols_processed_;
        if (symbols_processed_ < config_.s) {
            return;
        }

        uint64_t const hash_seed = config_.canonical ? config_.seed : 0ULL;

        uint64_t hash_s_forward = (forward_s_ & smask_) ^ hash_seed;
        forward_queue_.push(hash_s_forward, s_index_);

        if (config_.canonical) {
            uint64_t hash_s_reverse = (reverse_s_ & smask_) ^ hash_seed;
            reverse_queue_.push(hash_s_reverse, s_index_);
        }

        ++s_index_;

        if (s_index_ < window_size_) {
            return;
        }

        size_t window_start = s_index_ - window_size_;
        forward_queue_.pop_older(window_start);
        if (config_.canonical) {
            reverse_queue_.pop_older(window_start);
        }

        uint64_t hash_k_forward = (forward_k_ & kmask_) ^ hash_seed;

        if (!config_.canonical) {
            if (!forward_queue_.empty() && matches_position(window_start, forward_queue_.min_index(), false)) {
                callback(hash_k_forward);
            }
            return;
        }

        uint64_t hash_k_reverse = (reverse_k_ & kmask_) ^ hash_seed;
        bool choose_forward = hash_k_forward < hash_k_reverse;
        if (hash_k_forward == hash_k_reverse) {
            choose_forward = false;
        }

        if (choose_forward) {
            if (!forward_queue_.empty() && matches_position(window_start, forward_queue_.min_index(), false)) {
                callback(hash_k_forward);
            }
            return;
        }

        if (!reverse_queue_.empty() && matches_position(window_start, reverse_queue_.min_index(), true)) {
            callback(hash_k_reverse);
        }
    }

    SyncmerConfig config_;
    size_t window_size_{0};
    uint64_t kmask_{0};
    uint64_t smask_{0};
    uint32_t kshift_{0};
    uint32_t sshift_{0};

    uint64_t forward_k_{0};
    uint64_t reverse_k_{0};
    uint64_t forward_s_{0};
    uint64_t reverse_s_{0};
    size_t symbols_processed_{0};
    size_t s_index_{0};
    MinQueue forward_queue_;
    MinQueue reverse_queue_;
};

} // namespace strobemers::chimera
