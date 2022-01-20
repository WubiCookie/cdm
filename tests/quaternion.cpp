#include <common.hpp>

INFO_BEGIN(quaternion)

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
	{
		const quaternion q = quaternion(direction::posX(), 90_deg);
		CHECK_THAT(q * vector3(0, 1, 0), Vector3Matcher({0, 0, 1}));
		CHECK_THAT(q * vector3(0, 0, 1), Vector3Matcher({0, -1, 0}));
		CHECK_THAT(q * vector3(0, -1, 0), Vector3Matcher({0, 0, -1}));
		CHECK_THAT(q * vector3(0, 0, -1), Vector3Matcher({0, 1, 0}));
	}
	{
		const quaternion q = quaternion(direction::posY(), 90_deg);
		CHECK_THAT(q * vector3(1, 0, 0), Vector3Matcher({0, 0, -1}));
		CHECK_THAT(q * vector3(0, 0, -1), Vector3Matcher({-1, 0, 0}));
		CHECK_THAT(q * vector3(-1, 0, 0), Vector3Matcher({0, 0, 1}));
		CHECK_THAT(q * vector3(0, 0, 1), Vector3Matcher({1, 0, 0}));
	}
	{
		const quaternion q = quaternion(direction::posZ(), 90_deg);
		CHECK_THAT(q * vector3(1, 0, 0), Vector3Matcher({0, 1, 0}));
		CHECK_THAT(q * vector3(0, 1, 0), Vector3Matcher({-1, 0, 0}));
		CHECK_THAT(q * vector3(-1, 0, 0), Vector3Matcher({0, -1, 0}));
		CHECK_THAT(q * vector3(0, -1, 0), Vector3Matcher({1, 0, 0}));
	}
}

TEST_CASE(
    "quaternion::operator*(vector3) and quaternion::quaternion(direction, "
    "static_pi_fraction)",
    "[working][unittest][quaternion]")
{
	using HalfPi = static_pi_fraction<1, 2>;

	{
		const quaternion q = quaternion(direction(1.0f, 0.0f, 0.0f), HalfPi{});
		CHECK_THAT(q * vector3(0, 1, 0), Vector3Matcher({0, 0, 1}));
		CHECK_THAT(q * vector3(0, 0, 1), Vector3Matcher({0, -1, 0}));
		CHECK_THAT(q * vector3(0, -1, 0), Vector3Matcher({0, 0, -1}));
		CHECK_THAT(q * vector3(0, 0, -1), Vector3Matcher({0, 1, 0}));
	}
	{
		const quaternion q = quaternion(direction(0.0f, 1.0f, 0.0f), HalfPi{});
		CHECK_THAT(q * vector3(1, 0, 0), Vector3Matcher({0, 0, -1}));
		CHECK_THAT(q * vector3(0, 0, -1), Vector3Matcher({-1, 0, 0}));
		CHECK_THAT(q * vector3(-1, 0, 0), Vector3Matcher({0, 0, 1}));
		CHECK_THAT(q * vector3(0, 0, 1), Vector3Matcher({1, 0, 0}));
	}
	{
		const quaternion q = quaternion(direction(0.0f, 0.0f, 1.0f), HalfPi{});
		CHECK_THAT(q * vector3(1, 0, 0), Vector3Matcher({0, 1, 0}));
		CHECK_THAT(q * vector3(0, 1, 0), Vector3Matcher({-1, 0, 0}));
		CHECK_THAT(q * vector3(-1, 0, 0), Vector3Matcher({0, -1, 0}));
		CHECK_THAT(q * vector3(0, -1, 0), Vector3Matcher({1, 0, 0}));
	}
}

INFO_END(quaternion)
