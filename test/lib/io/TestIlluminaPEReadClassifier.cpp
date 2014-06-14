#include "io/IlluminaPEReadClassifier.hpp"

#include <gtest/gtest.h>

#include <vector>
#include <string>

TEST(TestIlluminaPE, classify) {
    std::vector<std::string> param_names{
        "rev",
        "mrev",
        "left",
        "ppair",
        "big",
        "sml"};

    int const n = 1 << param_names.size();

    std::cout << param_names[0];
    for (size_t j = 1; j < param_names.size(); ++j) {
        std::cout << "\t" << param_names[j];
    }
    std::cout << "\tflag\n";


    std::vector<ReadFlag> values(n);
    for (int i = 0; i < n; ++i) {
        values[i] = pe_classify(
            bool(i & 1),
            bool(i & 2),
            bool(i & 4),
            bool(i & 8),
            bool(i & 16),
            bool(i & 32)
            );

        for (size_t j = 0; j < param_names.size(); ++j) {
            std::cout << bool(i & (1 << j)) << "\t";
        }
        std::cout << FLAG_VALUES.string_name(values[i]) << "\n";
    }
}
