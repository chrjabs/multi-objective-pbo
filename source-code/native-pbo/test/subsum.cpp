#include <vector>
#include "auxiliary.hpp"

std::vector<std::pair<std::vector<int>, std::vector<int>>> testSets = {
    {{1, 1, 1, 1}, {0, 1, 2, 3, 4}},
    {{1, 1, 2, 2}, {0, 1, 2, 3, 4, 5, 6}},
    {{1, 4, 5, 9}, {0, 1, 4, 5, 6, 9, 10, 13, 14, 15, 18, 19}},
    {{1, 1, 1, 1, 1}, {0, 1, 2, 3, 4, 5}},
    {{1, 9, 10, 16, 16, 18, 20, 36, 72},
     {0,   1,   9,   10,  11,  16,  17,  18,  19,  20,  21,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,
      36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,
      58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
      80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,  99,  100, 101,
      102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123,
      124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145,
      146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
      168, 169, 170, 171, 172, 173, 177, 178, 179, 180, 181, 182, 187, 188, 189, 197, 198}}};

int test_all() {
  for (auto& test : testSets) {
    auto res = rs::aux::subsum::all<int, int>(test.first);
    std::sort(res.begin(), res.end());
    res.erase(std::unique(res.begin(), res.end()), res.end());
    if (res != test.second) return 1;
  }
  return 0;
}

int test_inc() {
  for (auto& test : testSets) {
    auto incComp = rs::aux::subsum::IncComp<int, int>(test.first, 16);
    for (size_t i = 1; i < test.second.size(); i++) {
      int nl = incComp.nextLarger();
      if (nl != test.second[i]) {
        return 1;
      }
    }
  }
  return 0;
}

int main() {
  if (test_all()) return 1;
  if (test_inc()) return 1;
  return 0;
}
