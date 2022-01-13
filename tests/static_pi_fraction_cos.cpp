#include <catch2/catch.hpp>
#include <cdm_maths.hpp>

#include <array>
#include <iostream>

using namespace cdm;

#define FRACTION_TEST_CASE_TYPE(T, NUMERATOR, DENOMINATOR)            \
	TEST_CASE("cos(static_pi_fraction) " #T " (PI " #NUMERATOR        \
	          "/" #DENOMINATOR ")",                                   \
	          "[working][unittest]")                                  \
	{                                                                 \
		static_pi_fraction<NUMERATOR, DENOMINATOR> f;                 \
		radian_t<T> r0 = f;                                           \
		radian_t<T> r1{(pi * T(NUMERATOR)) / T(DENOMINATOR)};         \
		degree_t<T> d0 = f;                                           \
		degree_t<T> d1{r0};                                           \
		degree_t<T> d2{r1};                                           \
                                                                      \
		CHECK(cdm::cos<T>(f) == Approx(cdm::cos(r0)).margin(1.0e-6)); \
		CHECK(cdm::cos<T>(f) == Approx(cdm::cos(r1)).margin(1.0e-6)); \
		CHECK(cdm::cos<T>(f) == Approx(cdm::cos(d0)).margin(1.0e-6)); \
		CHECK(cdm::cos<T>(f) == Approx(cdm::cos(d1)).margin(1.0e-6)); \
		CHECK(cdm::cos<T>(f) == Approx(cdm::cos(d2)).margin(1.0e-6)); \
		CHECK(cdm::cos<T>(f) ==                                       \
		      Approx(std::cos(T(pi * NUMERATOR) / T(DENOMINATOR)))    \
		          .margin(1.0e-6));                                   \
	}
#define FRACTION_TEST_CASE(NUMERATOR, DENOMINATOR)         \
	FRACTION_TEST_CASE_TYPE(float, NUMERATOR, DENOMINATOR) \
	FRACTION_TEST_CASE_TYPE(double, NUMERATOR, DENOMINATOR)

#define INSTANCIATE_TEST_CASES \
	FRACTION_TEST_CASE(0, 1)   \
	FRACTION_TEST_CASE(0, 2)   \
	FRACTION_TEST_CASE(0, 3)   \
	FRACTION_TEST_CASE(0, 4)   \
	FRACTION_TEST_CASE(0, 5)   \
	FRACTION_TEST_CASE(0, 6)   \
	FRACTION_TEST_CASE(0, 7)   \
	FRACTION_TEST_CASE(0, 8)   \
	FRACTION_TEST_CASE(0, 9)   \
	FRACTION_TEST_CASE(1, 1)   \
	FRACTION_TEST_CASE(1, 2)   \
	FRACTION_TEST_CASE(1, 3)   \
	FRACTION_TEST_CASE(1, 4)   \
	FRACTION_TEST_CASE(1, 5)   \
	FRACTION_TEST_CASE(1, 6)   \
	FRACTION_TEST_CASE(1, 7)   \
	FRACTION_TEST_CASE(1, 8)   \
	FRACTION_TEST_CASE(1, 9)   \
	FRACTION_TEST_CASE(2, 1)   \
	FRACTION_TEST_CASE(2, 2)   \
	FRACTION_TEST_CASE(2, 3)   \
	FRACTION_TEST_CASE(2, 4)   \
	FRACTION_TEST_CASE(2, 5)   \
	FRACTION_TEST_CASE(2, 6)   \
	FRACTION_TEST_CASE(2, 7)   \
	FRACTION_TEST_CASE(2, 8)   \
	FRACTION_TEST_CASE(2, 9)   \
	FRACTION_TEST_CASE(3, 1)   \
	FRACTION_TEST_CASE(3, 2)   \
	FRACTION_TEST_CASE(3, 3)   \
	FRACTION_TEST_CASE(3, 4)   \
	FRACTION_TEST_CASE(3, 5)   \
	FRACTION_TEST_CASE(3, 6)   \
	FRACTION_TEST_CASE(3, 7)   \
	FRACTION_TEST_CASE(3, 8)   \
	FRACTION_TEST_CASE(3, 9)   \
	FRACTION_TEST_CASE(4, 1)   \
	FRACTION_TEST_CASE(4, 2)   \
	FRACTION_TEST_CASE(4, 3)   \
	FRACTION_TEST_CASE(4, 4)   \
	FRACTION_TEST_CASE(4, 5)   \
	FRACTION_TEST_CASE(4, 6)   \
	FRACTION_TEST_CASE(4, 7)   \
	FRACTION_TEST_CASE(4, 8)   \
	FRACTION_TEST_CASE(4, 9)   \
	FRACTION_TEST_CASE(5, 1)   \
	FRACTION_TEST_CASE(5, 2)   \
	FRACTION_TEST_CASE(5, 3)   \
	FRACTION_TEST_CASE(5, 4)   \
	FRACTION_TEST_CASE(5, 5)   \
	FRACTION_TEST_CASE(5, 6)   \
	FRACTION_TEST_CASE(5, 7)   \
	FRACTION_TEST_CASE(5, 8)   \
	FRACTION_TEST_CASE(5, 9)   \
	FRACTION_TEST_CASE(6, 1)   \
	FRACTION_TEST_CASE(6, 2)   \
	FRACTION_TEST_CASE(6, 3)   \
	FRACTION_TEST_CASE(6, 4)   \
	FRACTION_TEST_CASE(6, 5)   \
	FRACTION_TEST_CASE(6, 6)   \
	FRACTION_TEST_CASE(6, 7)   \
	FRACTION_TEST_CASE(6, 8)   \
	FRACTION_TEST_CASE(6, 9)   \
	FRACTION_TEST_CASE(7, 1)   \
	FRACTION_TEST_CASE(7, 2)   \
	FRACTION_TEST_CASE(7, 3)   \
	FRACTION_TEST_CASE(7, 4)   \
	FRACTION_TEST_CASE(7, 5)   \
	FRACTION_TEST_CASE(7, 6)   \
	FRACTION_TEST_CASE(7, 7)   \
	FRACTION_TEST_CASE(7, 8)   \
	FRACTION_TEST_CASE(7, 9)   \
	FRACTION_TEST_CASE(8, 1)   \
	FRACTION_TEST_CASE(8, 2)   \
	FRACTION_TEST_CASE(8, 3)   \
	FRACTION_TEST_CASE(8, 4)   \
	FRACTION_TEST_CASE(8, 5)   \
	FRACTION_TEST_CASE(8, 6)   \
	FRACTION_TEST_CASE(8, 7)   \
	FRACTION_TEST_CASE(8, 8)   \
	FRACTION_TEST_CASE(8, 9)   \
	FRACTION_TEST_CASE(9, 1)   \
	FRACTION_TEST_CASE(9, 2)   \
	FRACTION_TEST_CASE(9, 3)   \
	FRACTION_TEST_CASE(9, 4)   \
	FRACTION_TEST_CASE(9, 5)   \
	FRACTION_TEST_CASE(9, 6)   \
	FRACTION_TEST_CASE(9, 7)   \
	FRACTION_TEST_CASE(9, 8)   \
	FRACTION_TEST_CASE(9, 9)   \
	FRACTION_TEST_CASE(-1, 1)  \
	FRACTION_TEST_CASE(-1, 2)  \
	FRACTION_TEST_CASE(-1, 3)  \
	FRACTION_TEST_CASE(-1, 4)  \
	FRACTION_TEST_CASE(-1, 5)  \
	FRACTION_TEST_CASE(-1, 6)  \
	FRACTION_TEST_CASE(-1, 7)  \
	FRACTION_TEST_CASE(-1, 8)  \
	FRACTION_TEST_CASE(-1, 9)  \
	FRACTION_TEST_CASE(-2, 1)  \
	FRACTION_TEST_CASE(-2, 2)  \
	FRACTION_TEST_CASE(-2, 3)  \
	FRACTION_TEST_CASE(-2, 4)  \
	FRACTION_TEST_CASE(-2, 5)  \
	FRACTION_TEST_CASE(-2, 6)  \
	FRACTION_TEST_CASE(-2, 7)  \
	FRACTION_TEST_CASE(-2, 8)  \
	FRACTION_TEST_CASE(-2, 9)  \
	FRACTION_TEST_CASE(-3, 1)  \
	FRACTION_TEST_CASE(-3, 2)  \
	FRACTION_TEST_CASE(-3, 3)  \
	FRACTION_TEST_CASE(-3, 4)  \
	FRACTION_TEST_CASE(-3, 5)  \
	FRACTION_TEST_CASE(-3, 6)  \
	FRACTION_TEST_CASE(-3, 7)  \
	FRACTION_TEST_CASE(-3, 8)  \
	FRACTION_TEST_CASE(-3, 9)  \
	FRACTION_TEST_CASE(-4, 1)  \
	FRACTION_TEST_CASE(-4, 2)  \
	FRACTION_TEST_CASE(-4, 3)  \
	FRACTION_TEST_CASE(-4, 4)  \
	FRACTION_TEST_CASE(-4, 5)  \
	FRACTION_TEST_CASE(-4, 6)  \
	FRACTION_TEST_CASE(-4, 7)  \
	FRACTION_TEST_CASE(-4, 8)  \
	FRACTION_TEST_CASE(-4, 9)  \
	FRACTION_TEST_CASE(-5, 1)  \
	FRACTION_TEST_CASE(-5, 2)  \
	FRACTION_TEST_CASE(-5, 3)  \
	FRACTION_TEST_CASE(-5, 4)  \
	FRACTION_TEST_CASE(-5, 5)  \
	FRACTION_TEST_CASE(-5, 6)  \
	FRACTION_TEST_CASE(-5, 7)  \
	FRACTION_TEST_CASE(-5, 8)  \
	FRACTION_TEST_CASE(-5, 9)  \
	FRACTION_TEST_CASE(-6, 1)  \
	FRACTION_TEST_CASE(-6, 2)  \
	FRACTION_TEST_CASE(-6, 3)  \
	FRACTION_TEST_CASE(-6, 4)  \
	FRACTION_TEST_CASE(-6, 5)  \
	FRACTION_TEST_CASE(-6, 6)  \
	FRACTION_TEST_CASE(-6, 7)  \
	FRACTION_TEST_CASE(-6, 8)  \
	FRACTION_TEST_CASE(-6, 9)  \
	FRACTION_TEST_CASE(-7, 1)  \
	FRACTION_TEST_CASE(-7, 2)  \
	FRACTION_TEST_CASE(-7, 3)  \
	FRACTION_TEST_CASE(-7, 4)  \
	FRACTION_TEST_CASE(-7, 5)  \
	FRACTION_TEST_CASE(-7, 6)  \
	FRACTION_TEST_CASE(-7, 7)  \
	FRACTION_TEST_CASE(-7, 8)  \
	FRACTION_TEST_CASE(-7, 9)  \
	FRACTION_TEST_CASE(-8, 1)  \
	FRACTION_TEST_CASE(-8, 2)  \
	FRACTION_TEST_CASE(-8, 3)  \
	FRACTION_TEST_CASE(-8, 4)  \
	FRACTION_TEST_CASE(-8, 5)  \
	FRACTION_TEST_CASE(-8, 6)  \
	FRACTION_TEST_CASE(-8, 7)  \
	FRACTION_TEST_CASE(-8, 8)  \
	FRACTION_TEST_CASE(-8, 9)  \
	FRACTION_TEST_CASE(-9, 1)  \
	FRACTION_TEST_CASE(-9, 2)  \
	FRACTION_TEST_CASE(-9, 3)  \
	FRACTION_TEST_CASE(-9, 4)  \
	FRACTION_TEST_CASE(-9, 5)  \
	FRACTION_TEST_CASE(-9, 6)  \
	FRACTION_TEST_CASE(-9, 7)  \
	FRACTION_TEST_CASE(-9, 8)  \
	FRACTION_TEST_CASE(-9, 9)  \
	FRACTION_TEST_CASE(-1, -1) \
	FRACTION_TEST_CASE(-1, -2) \
	FRACTION_TEST_CASE(-1, -3) \
	FRACTION_TEST_CASE(-1, -4) \
	FRACTION_TEST_CASE(-1, -5) \
	FRACTION_TEST_CASE(-1, -6) \
	FRACTION_TEST_CASE(-1, -7) \
	FRACTION_TEST_CASE(-1, -8) \
	FRACTION_TEST_CASE(-1, -9) \
	FRACTION_TEST_CASE(-2, -1) \
	FRACTION_TEST_CASE(-2, -2) \
	FRACTION_TEST_CASE(-2, -3) \
	FRACTION_TEST_CASE(-2, -4) \
	FRACTION_TEST_CASE(-2, -5) \
	FRACTION_TEST_CASE(-2, -6) \
	FRACTION_TEST_CASE(-2, -7) \
	FRACTION_TEST_CASE(-2, -8) \
	FRACTION_TEST_CASE(-2, -9) \
	FRACTION_TEST_CASE(-3, -1) \
	FRACTION_TEST_CASE(-3, -2) \
	FRACTION_TEST_CASE(-3, -3) \
	FRACTION_TEST_CASE(-3, -4) \
	FRACTION_TEST_CASE(-3, -5) \
	FRACTION_TEST_CASE(-3, -6) \
	FRACTION_TEST_CASE(-3, -7) \
	FRACTION_TEST_CASE(-3, -8) \
	FRACTION_TEST_CASE(-3, -9) \
	FRACTION_TEST_CASE(-4, -1) \
	FRACTION_TEST_CASE(-4, -2) \
	FRACTION_TEST_CASE(-4, -3) \
	FRACTION_TEST_CASE(-4, -4) \
	FRACTION_TEST_CASE(-4, -5) \
	FRACTION_TEST_CASE(-4, -6) \
	FRACTION_TEST_CASE(-4, -7) \
	FRACTION_TEST_CASE(-4, -8) \
	FRACTION_TEST_CASE(-4, -9) \
	FRACTION_TEST_CASE(-5, -1) \
	FRACTION_TEST_CASE(-5, -2) \
	FRACTION_TEST_CASE(-5, -3) \
	FRACTION_TEST_CASE(-5, -4) \
	FRACTION_TEST_CASE(-5, -5) \
	FRACTION_TEST_CASE(-5, -6) \
	FRACTION_TEST_CASE(-5, -7) \
	FRACTION_TEST_CASE(-5, -8) \
	FRACTION_TEST_CASE(-5, -9) \
	FRACTION_TEST_CASE(-6, -1) \
	FRACTION_TEST_CASE(-6, -2) \
	FRACTION_TEST_CASE(-6, -3) \
	FRACTION_TEST_CASE(-6, -4) \
	FRACTION_TEST_CASE(-6, -5) \
	FRACTION_TEST_CASE(-6, -6) \
	FRACTION_TEST_CASE(-6, -7) \
	FRACTION_TEST_CASE(-6, -8) \
	FRACTION_TEST_CASE(-6, -9) \
	FRACTION_TEST_CASE(-7, -1) \
	FRACTION_TEST_CASE(-7, -2) \
	FRACTION_TEST_CASE(-7, -3) \
	FRACTION_TEST_CASE(-7, -4) \
	FRACTION_TEST_CASE(-7, -5) \
	FRACTION_TEST_CASE(-7, -6) \
	FRACTION_TEST_CASE(-7, -7) \
	FRACTION_TEST_CASE(-7, -8) \
	FRACTION_TEST_CASE(-7, -9) \
	FRACTION_TEST_CASE(-8, -1) \
	FRACTION_TEST_CASE(-8, -2) \
	FRACTION_TEST_CASE(-8, -3) \
	FRACTION_TEST_CASE(-8, -4) \
	FRACTION_TEST_CASE(-8, -5) \
	FRACTION_TEST_CASE(-8, -6) \
	FRACTION_TEST_CASE(-8, -7) \
	FRACTION_TEST_CASE(-8, -8) \
	FRACTION_TEST_CASE(-8, -9) \
	FRACTION_TEST_CASE(-9, -1) \
	FRACTION_TEST_CASE(-9, -2) \
	FRACTION_TEST_CASE(-9, -3) \
	FRACTION_TEST_CASE(-9, -4) \
	FRACTION_TEST_CASE(-9, -5) \
	FRACTION_TEST_CASE(-9, -6) \
	FRACTION_TEST_CASE(-9, -7) \
	FRACTION_TEST_CASE(-9, -8) \
	FRACTION_TEST_CASE(-9, -9) \
	FRACTION_TEST_CASE(1, -1)  \
	FRACTION_TEST_CASE(1, -2)  \
	FRACTION_TEST_CASE(1, -3)  \
	FRACTION_TEST_CASE(1, -4)  \
	FRACTION_TEST_CASE(1, -5)  \
	FRACTION_TEST_CASE(1, -6)  \
	FRACTION_TEST_CASE(1, -7)  \
	FRACTION_TEST_CASE(1, -8)  \
	FRACTION_TEST_CASE(1, -9)  \
	FRACTION_TEST_CASE(2, -1)  \
	FRACTION_TEST_CASE(2, -2)  \
	FRACTION_TEST_CASE(2, -3)  \
	FRACTION_TEST_CASE(2, -4)  \
	FRACTION_TEST_CASE(2, -5)  \
	FRACTION_TEST_CASE(2, -6)  \
	FRACTION_TEST_CASE(2, -7)  \
	FRACTION_TEST_CASE(2, -8)  \
	FRACTION_TEST_CASE(2, -9)  \
	FRACTION_TEST_CASE(3, -1)  \
	FRACTION_TEST_CASE(3, -2)  \
	FRACTION_TEST_CASE(3, -3)  \
	FRACTION_TEST_CASE(3, -4)  \
	FRACTION_TEST_CASE(3, -5)  \
	FRACTION_TEST_CASE(3, -6)  \
	FRACTION_TEST_CASE(3, -7)  \
	FRACTION_TEST_CASE(3, -8)  \
	FRACTION_TEST_CASE(3, -9)  \
	FRACTION_TEST_CASE(4, -1)  \
	FRACTION_TEST_CASE(4, -2)  \
	FRACTION_TEST_CASE(4, -3)  \
	FRACTION_TEST_CASE(4, -4)  \
	FRACTION_TEST_CASE(4, -5)  \
	FRACTION_TEST_CASE(4, -6)  \
	FRACTION_TEST_CASE(4, -7)  \
	FRACTION_TEST_CASE(4, -8)  \
	FRACTION_TEST_CASE(4, -9)  \
	FRACTION_TEST_CASE(5, -1)  \
	FRACTION_TEST_CASE(5, -2)  \
	FRACTION_TEST_CASE(5, -3)  \
	FRACTION_TEST_CASE(5, -4)  \
	FRACTION_TEST_CASE(5, -5)  \
	FRACTION_TEST_CASE(5, -6)  \
	FRACTION_TEST_CASE(5, -7)  \
	FRACTION_TEST_CASE(5, -8)  \
	FRACTION_TEST_CASE(5, -9)  \
	FRACTION_TEST_CASE(6, -1)  \
	FRACTION_TEST_CASE(6, -2)  \
	FRACTION_TEST_CASE(6, -3)  \
	FRACTION_TEST_CASE(6, -4)  \
	FRACTION_TEST_CASE(6, -5)  \
	FRACTION_TEST_CASE(6, -6)  \
	FRACTION_TEST_CASE(6, -7)  \
	FRACTION_TEST_CASE(6, -8)  \
	FRACTION_TEST_CASE(6, -9)  \
	FRACTION_TEST_CASE(7, -1)  \
	FRACTION_TEST_CASE(7, -2)  \
	FRACTION_TEST_CASE(7, -3)  \
	FRACTION_TEST_CASE(7, -4)  \
	FRACTION_TEST_CASE(7, -5)  \
	FRACTION_TEST_CASE(7, -6)  \
	FRACTION_TEST_CASE(7, -7)  \
	FRACTION_TEST_CASE(7, -8)  \
	FRACTION_TEST_CASE(7, -9)  \
	FRACTION_TEST_CASE(8, -1)  \
	FRACTION_TEST_CASE(8, -2)  \
	FRACTION_TEST_CASE(8, -3)  \
	FRACTION_TEST_CASE(8, -4)  \
	FRACTION_TEST_CASE(8, -5)  \
	FRACTION_TEST_CASE(8, -6)  \
	FRACTION_TEST_CASE(8, -7)  \
	FRACTION_TEST_CASE(8, -8)  \
	FRACTION_TEST_CASE(8, -9)  \
	FRACTION_TEST_CASE(9, -1)  \
	FRACTION_TEST_CASE(9, -2)  \
	FRACTION_TEST_CASE(9, -3)  \
	FRACTION_TEST_CASE(9, -4)  \
	FRACTION_TEST_CASE(9, -5)  \
	FRACTION_TEST_CASE(9, -6)  \
	FRACTION_TEST_CASE(9, -7)  \
	FRACTION_TEST_CASE(9, -8)  \
	FRACTION_TEST_CASE(9, -9)

INSTANCIATE_TEST_CASES
