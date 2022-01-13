#include <catch2/catch.hpp>
#include <cdm_maths.hpp>

#include <array>
#include <iostream>

using namespace cdm;

TEST_CASE("matrix3::rotation_around_x()", "[inProgress]")
{
	radian angle = 90_deg;
	matrix3 m1 = matrix3::rotation_around_x(angle);
	matrix3 m2 = matrix3(std::array<float, 9>{1.0f, 0.0f, 0.0f, 0.0f,
	                                          cos(angle), sin(angle), 0.0f,
	                                          -sin(angle), cos(angle)});

	CHECK(m1.column(0).row(0) == Approx(m2.column(0).row(0)).margin(5.0e-8));
	CHECK(m1.column(0).row(1) == Approx(m2.column(0).row(1)).margin(5.0e-8));
	CHECK(m1.column(0).row(2) == Approx(m2.column(0).row(2)).margin(5.0e-8));
	CHECK(m1.column(1).row(0) == Approx(m2.column(1).row(0)).margin(5.0e-8));
	CHECK(m1.column(1).row(1) == Approx(m2.column(1).row(1)).margin(5.0e-8));
	CHECK(m1.column(1).row(2) == Approx(m2.column(1).row(2)).margin(5.0e-8));
	CHECK(m1.column(2).row(0) == Approx(m2.column(2).row(0)).margin(5.0e-8));
	CHECK(m1.column(2).row(1) == Approx(m2.column(2).row(1)).margin(5.0e-8));
	CHECK(m1.column(2).row(2) == Approx(m2.column(2).row(2)).margin(5.0e-8));
}
