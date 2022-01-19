#include <common.hpp>

TEST_CASE("quaternion::operator*(vector3) and quaternion::identity()",
          "[working][unittest][quaternion]")
{
	const quaternion q1 = quaternion::identity();

	float x = GENERATE(-1.0f, 0.0f, 1.0f);
	float y = GENERATE(-1.0f, 0.0f, 1.0f);
	float z = GENERATE(-1.0f, 0.0f, 1.0f);

	CHECK((q1 * vector3{x, y, z}) == vector3{x, y, z});
}

TEST_CASE(
    "quaternion::operator*(vector3) and quaternion::quaternion(direction, "
    "radian)",
    "[working][unittest][quaternion]")
{
	const quaternion q1 = quaternion(direction(0.0f, 0.0f, 1.0f), 90_deg);

	const vector3 v1{1, 0, 0};
	const vector3 v2 = q1 * v1;

	CHECK_THAT(v2, Vector3Matcher({0, 1, 0}, 1.0e-6));
}

TEST_CASE(
    "quaternion::operator*(vector3) and quaternion::quaternion(direction, "
    "static_pi_fraction)",
    "[working][unittest][quaternion]")
{
	using HalfPi = static_pi_fraction<1, 2>;
	const quaternion q1 = quaternion(direction(0.0f, 0.0f, 1.0f), HalfPi{});

	const vector3 v1{1, 0, 0};
	const vector3 v2 = q1 * v1;
	
	CHECK_THAT(v2, Vector3Matcher({0, 1, 0}, 1.0e-6));
}
