/***********************************************************************
Copyright (c) 2014-2020, Jan Elffers
Copyright (c) 2019-2021, Jo Devriendt
Copyright (c) 2020-2021, Stephan Gocht
Copyright (c) 2014-2021, Jakob Nordström

Parts of the code were copied or adapted from MiniSat.

MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010  Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***********************************************************************/

// NOTE: This file may not be called auxiliary.hpp, because this is a forbidden file name on windows.

#pragma once

#define EXPANDED(x) STR(x)
#define STR(x) #x

#define _unused(x) ((void)(x))  // marks variables unused in release mode, use [[maybe_unused]] where possible

#define INSTANTIATE_CLASS(x)           \
  template class x<int, long long>;    \
  template class x<long long, int128>; \
  template class x<int128, int128>;    \
  template class x<int128, int256>;    \
  template class x<bigint, bigint>;

#define INSTANTIATE_STRUCT(x)           \
  template struct x<int, long long>;    \
  template struct x<long long, int128>; \
  template struct x<int128, int128>;    \
  template struct x<int128, int256>;    \
  template struct x<bigint, bigint>;

#include <sys/resource.h>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <random>
#include <unordered_map>
#include <vector>
#include "typedefs.hpp"

namespace rs {

template <typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::pair<T, U>& p) {
  os << p.first << "," << p.second;
  return os;
}
template <typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<T, U>& m) {
  for (const auto& e : m) os << e << ";";
  return os;
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& m) {
  for (const auto& e : m) os << e << " ";
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const __int128& x) {
  if (x < 0) return os << "-" << -x;
  if (x < 10) return os << (char)(x + '0');
  return os << x / 10 << (char)(x % 10 + '0');
}

namespace aux {

extern std::default_random_engine rng;

template <typename T>
T sto(const std::string& s) {
  return std::stold(s);
}
template <>
inline double sto(const std::string& s) {
  return std::stod(s);
}
template <>
inline std::string sto(const std::string& s) {
  return s;
}

template <typename T>
std::string tos(const T& t) {
  return std::to_string(t);
}
template <>
inline std::string tos(const std::string& s) {
  return s;
}

template <typename T>
void swapErase(T& indexable, size_t index) {
  indexable[index] = std::move(indexable.back());
  indexable.pop_back();
}

template <typename T, typename U>
bool contains(const T& v, const U& x) {
  return std::find(v.cbegin(), v.cend(), x) != v.cend();
}

template <typename T>
T ceildiv(const T& p, const T& q) {
  assert(q > 0);
  assert(p >= 0);
  return (p + q - 1) / q;
}  // NOTE: potential overflow
template <typename T>
T floordiv(const T& p, const T& q) {
  assert(q > 0);
  assert(p >= 0);
  return p / q;
}
template <typename T>
T ceildiv_safe(const T& p, const T& q) {
  assert(q > 0);
  return (p < 0) ? (-floordiv<T>(-p, q)) : ceildiv(p, q);
}
template <typename T>
T floordiv_safe(const T& p, const T& q) {
  assert(q > 0);
  return (p < 0) ? (-ceildiv<T>(-p, q)) : floordiv(p, q);
}
template <typename T>
T mod_safe(const T& p, const T& q) {
  assert(q > 0);
  if (p < 0)
    return q - (-p % q);
  else
    return p % q;
}

// Minisat cpuTime function
static inline double cpuTime() {
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000;
}

template <typename T>
void resizeIntMap(std::vector<T>& _map, typename std::vector<T>::iterator& map, int size, int resize_factor, T init) {
  assert(size < (1 << 28));
  int oldsize = (_map.size() - 1) / 2;
  if (oldsize >= size) return;
  int newsize = oldsize;
  while (newsize < size) newsize = newsize * resize_factor + 1;
  _map.resize(2 * newsize + 1);
  map = _map.begin() + newsize;
  int i = _map.size() - 1;
  for (; i > newsize + oldsize; --i) _map[i] = init;
  for (; i >= newsize - oldsize; --i) _map[i] = _map[i - newsize + oldsize];
  for (; i >= 0; --i) _map[i] = init;
}

template <typename T>
T median(std::vector<T>& v) {
  assert(v.size() > 0);
  size_t n = v.size() / 2;
  std::nth_element(v.begin(), v.begin() + n, v.end());
  return v[n];
}

template <typename T>
double average(const std::vector<T>& v) {
  assert(v.size() > 0);
  return std::accumulate(v.begin(), v.end(), 0.0) / (double)v.size();
}

template <typename T>
T min(const std::vector<T>& v) {
  return *std::min_element(v.begin(), v.end());
}

template <typename T>
T max(const std::vector<T>& v) {
  return *std::max_element(v.begin(), v.end());
}

template <typename T>
T abs(const T& x) {
  return std::abs(x);
}
template <>
inline bigint abs(const bigint& x) {
  return boost::multiprecision::abs(x);
}
template <>
inline boost::multiprecision::int128_t abs(const boost::multiprecision::int128_t& x) {
  return boost::multiprecision::abs(x);
}
template <>
inline int256 abs(const int256& x) {
  return boost::multiprecision::abs(x);
}

template <typename T>
T gcd(const T& x, const T& y) {
  return std::gcd(x, y);
}
template <>
inline bigint gcd(const bigint& x, const bigint& y) {
  return boost::multiprecision::gcd(x, y);
}
template <>
inline boost::multiprecision::int128_t gcd(const boost::multiprecision::int128_t& x,
                                           const boost::multiprecision::int128_t& y) {
  return boost::multiprecision::gcd(x, y);
}
template <>
inline int256 gcd(const int256& x, const int256& y) {
  return boost::multiprecision::gcd(x, y);
}

template <typename T>
T lcm(const T& x, const T& y) {
  return std::lcm(x, y);
}
template <>
inline bigint lcm(const bigint& x, const bigint& y) {
  return boost::multiprecision::lcm(x, y);
}
template <>
inline boost::multiprecision::int128_t lcm(const boost::multiprecision::int128_t& x,
                                           const boost::multiprecision::int128_t& y) {
  return boost::multiprecision::lcm(x, y);
}
template <>
inline int256 lcm(const int256& x, const int256& y) {
  return boost::multiprecision::lcm(x, y);
}

template <typename T>
unsigned msb(const T& x) {
  assert(x > 0);
  // return std::bit_floor(x); // C++20
  return boost::multiprecision::msb(boost::multiprecision::uint128_t(x));
}
template <>
inline unsigned msb(const bigint& x) {
  assert(x > 0);
  return boost::multiprecision::msb(x);
}
template <>
inline unsigned msb(const boost::multiprecision::int128_t& x) {
  assert(x > 0);
  return boost::multiprecision::msb(x);
}
template <>
inline unsigned msb(const int256& x) {
  assert(x > 0);
  return boost::multiprecision::msb(x);
}

template <typename T>
T pow(const T& x, unsigned y) {
  return std::pow(x, y);
}
template <>
inline bigint pow(const bigint& x, unsigned y) {
  return boost::multiprecision::pow(x, y);
}
template <>
inline boost::multiprecision::int128_t pow(const boost::multiprecision::int128_t& x, unsigned y) {
  return boost::multiprecision::pow(x, y);
}
template <>
inline int256 pow(const int256& x, unsigned y) {
  return boost::multiprecision::pow(x, y);
}

template <typename T>
T timeCall(const std::function<T(void)>& f, double& to) {
  double start = cpuTime();
  T result = f();
  to += cpuTime() - start;
  return result;
}
template <>
inline void timeCall(const std::function<void(void)>& f, double& to) {
  double start = cpuTime();
  f();
  to += cpuTime() - start;
}

template <typename T>
bool fits([[maybe_unused]] const bigint& x) {
  return false;
}
template <>
inline bool fits<int>(const bigint& x) {
  return aux::abs(x) <= bigint(limit32);
}
template <>
inline bool fits<long long>(const bigint& x) {
  return aux::abs(x) <= bigint(limit64);
}
template <>
inline bool fits<int128>(const bigint& x) {
  return aux::abs(x) <= bigint(limit128);
}
template <>
inline bool fits<int256>(const bigint& x) {
  return aux::abs(x) <= bigint(limit256);
}
template <>
inline bool fits<bigint>([[maybe_unused]] const bigint& x) {
  return true;
}
template <typename T, typename S>
bool fitsIn([[maybe_unused]] const S& x) {
  return fits<T>(bigint(x));
}

#ifdef __APPLE__
template <>
inline int128 abs(const int128& x) {
  return x < 0 ? -x : x;
}
template <>
inline int128 gcd(const int128& x, const int128& y) {
  return static_cast<int128>(
      boost::multiprecision::gcd(boost::multiprecision::int128_t(x), boost::multiprecision::int128_t(y)));
}
template <>
inline unsigned msb(const int128& x) {
  return boost::multiprecision::msb(boost::multiprecision::uint128_t(x));
}
template <>
inline int128 pow(const int128& x, unsigned y) {
  return static_cast<int128>(boost::multiprecision::pow(boost::multiprecision::int128_t(x), y));
}
#endif

namespace subsum {
// subset sum algorithms used for (incremental) pb reformulations

// compute _all_ subset sums. note that this has exponential space and probably
// shouldn't be used for that reason.
// https://codereview.stackexchange.com/questions/224964/print-sums-of-all-subsets
template <typename SMALL, typename LARGE>
std::vector<LARGE> all(const std::vector<SMALL>& set) {
  std::vector<LARGE> sums{0};
  sums.reserve(std::size_t(1) << set.size());
  for (SMALL val : set)
    std::transform(sums.begin(), sums.end(), std::back_inserter(sums), [=](LARGE sum) { return sum + val; });
  return sums;
}

// incremental subset sum computer, computes the next larger subset sum with
// smaller memory usage
// https://math.stackexchange.com/questions/89419/algorithm-wanted-enumerate-all-subsets-of-a-set-in-order-of-increasing-sums
template <typename SMALL, typename LARGE>
struct IncComp {
  // the input set and the index(es) of the next considered subset
  std::vector<std::pair<SMALL, std::vector<int64_t>>> set;
  // the computed subset sums and a bitmap representing the set
  std::vector<std::pair<LARGE, uint64_t>> sums = {{0, 0}};
  // the threshold of number of sums cached at which cleanup will be triggered
  unsigned int CLEAN_THR;

  IncComp(const std::vector<SMALL>& _set, int cthr = 64) : CLEAN_THR(cthr) {
    assert(_set.size() <= 64);
    std::vector<SMALL> set_copy = _set;
    std::sort(set_copy.begin(), set_copy.end());
    SMALL last = -1;
    for (SMALL val : set_copy) {
      assert(val > 0);
      if (val != last) {
        assert(val > last);
        std::vector<int64_t> vec = {};
        set.emplace_back(std::make_pair(val, vec));
      }
      set.back().second.push_back(0);
      last = val;
    }
  }

  // computes the mask for the selected element
  inline int find_idx(unsigned primary, unsigned secondary) {
    assert(primary < set.size());
    assert(secondary < set[primary].second.size());
    unsigned cnt = 0;
    for (unsigned i = 0; i < primary; i++) cnt += set[i].second.size();
    cnt += secondary;
    return cnt;
  }

  // computes the next subset sum
  bool _step() {
    // find candidate to add to a subset
    int minIdx = -1;
    unsigned secIdx = 0;
    LARGE minVal = -1;
    for (unsigned i = 0; i < set.size(); i++) {
      for (unsigned j = 0; j < set[i].second.size(); j++) {
        int idx = set[i].second[j];
        if (idx >= 0 && (j == 0 || sums[idx].second & (1ul << find_idx(i, j - 1))) &&
            (sums[idx].first + set[i].first < minVal || minVal < 0)) {
          minIdx = i;
          secIdx = j;
          minVal = sums[idx].first + set[i].first;
        }
      }
    }
    if (minVal < 0) return false;
    assert(minIdx >= 0);
    assert(static_cast<size_t>(minIdx) < set.size());
    assert(secIdx < set[minIdx].second.size());
    int trueIdx = find_idx(minIdx, secIdx);
    uint64_t negMask = ~0ul << trueIdx;
    uint64_t posMask = 1ul << trueIdx;
    uint64_t setRep = sums[set[minIdx].second[secIdx]].second;
    // sanity check: neither the select element nor anything bigger is in the selected sum already
    assert(!(setRep & negMask));
    assert(secIdx == 0 || (setRep & (1ul << (trueIdx - 1))));
    // generate new subset, setting the minIdx bit in the
    sums.emplace_back(std::make_pair(minVal, setRep | posMask));
    // update pointer of next-higher cardinality version of current element
    // this iteratively initializes the next-higher cardinality version correctly
    if (secIdx + 1 < set[minIdx].second.size() && !(sums[set[minIdx].second[secIdx + 1]].second & posMask))
      set[minIdx].second[secIdx + 1] = sums.size() - 1;
    // compute next candidate set for the just used element. the next candidate is the next sum that does not include an
    // element with higher index (i.e., higher value or higher precedence)
    int64_t& nextSum = set[minIdx].second[secIdx];
    nextSum++;
    while (static_cast<size_t>(nextSum) < sums.size() &&
           (sums[nextSum].second & negMask || (secIdx > 0 && !(sums[nextSum].second & (1ul << (trueIdx - 1))))))
      nextSum++;
    // finished with this element if no new candidate found
    if (static_cast<size_t>(nextSum) >= sums.size()) nextSum = -1;
    assert(nextSum < 0 || secIdx == 0 || sums[nextSum].second & (1ul << (trueIdx - 1)));
    return true;
  }

  // clean unneeded data from the datastructure
  void _cleanup() {
    // remove all sums from the beginning that are not needed anymore
    int64_t minIdx = sums.size();
    for (auto& elem : set) {
      for (int64_t idx : elem.second)
        if (idx >= 0 && idx < minIdx) minIdx = idx;
    }
    if (minIdx > 0 && minIdx < (int64_t)sums.size()) {
      sums.erase(sums.begin(), sums.begin() + minIdx);
      for (auto& elem : set) {
        for (int64_t& idx : elem.second)
          if (idx > 0) idx -= minIdx;
      }
    }
    // potentially one could try to remove sums in the middle of the references as well, but that seems like quite a lot
    // more work for probably not that much reward
    // removing elements that are done is not feasible since all bitmaps would
    // need to be recomputed
  }

  // gets the next subset sum larger than thr
  LARGE nextLarger(LARGE thr = -1, bool equal = false) {
    if (thr < 0) thr = sums.back().first;
    if (equal) thr--;
    while (sums.back().first <= thr && _step()) {
    };
    if (sums.size() > CLEAN_THR) _cleanup();
    return (sums.back().first > thr) ? sums.back().first : -1;
  }
};

}  // namespace subsum

}  // namespace aux

}  // namespace rs
