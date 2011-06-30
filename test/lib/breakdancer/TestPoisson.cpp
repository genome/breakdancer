#include "breakdancer/Poisson.h"

#include <cmath>
#include <gtest/gtest.h>

TEST(Poisson, LAdd) {
    double a = log(30);
    double b = log(20);
    double expected = log(50);
    ASSERT_NEAR(expected, LAdd(a, b), 1e-6);
}
