#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <common.hpp>

TEST_CASE("quaternion::identity()", "[working][unittest][quaternion]")
{
	const quaternion q1 = quaternion::identity();

	float x = GENERATE(-1.0f, 0.0f, 1.0f);
	float y = GENERATE(-1.0f, 0.0f, 1.0f);
	float z = GENERATE(-1.0f, 0.0f, 1.0f);

	CHECK((q1 * vector3{x, y, z}) == vector3{x, y, z});
}

TEST_CASE("quaternion::constructors", "[working][unittest][quaternion]")
{
	const quaternion q1{1, 2, 3, 4};
	auto a1 = q1.to_array();
	const quaternion q2 = a1;

	CHECK(q1 == q2);
}

TEST_CASE("quaternion::quaternion(direction, radian)",
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

TEST_CASE("quaternion::quaternion(direction, static_pi_fraction)",
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

TEST_CASE("quaternion::inverse()", "[working][unittest][quaternion]")
{
	using HalfPi = static_pi_fraction<1, 2>;

	{
		quaternion q = quaternion(direction(1.0f, 0.0f, 0.0f), HalfPi{});
		q.inverse();
		CHECK_THAT(q * vector3(0, 0, 1), Vector3Matcher({0, 1, 0}));
		CHECK_THAT(q * vector3(0, -1, 0), Vector3Matcher({0, 0, 1}));
		CHECK_THAT(q * vector3(0, 0, -1), Vector3Matcher({0, -1, 0}));
		CHECK_THAT(q * vector3(0, 1, 0), Vector3Matcher({0, 0, -1}));
	}
	{
		quaternion q = quaternion(direction(0.0f, 1.0f, 0.0f), HalfPi{});
		q.inverse();
		CHECK_THAT(q * vector3(0, 0, -1), Vector3Matcher({1, 0, 0}));
		CHECK_THAT(q * vector3(-1, 0, 0), Vector3Matcher({0, 0, -1}));
		CHECK_THAT(q * vector3(0, 0, 1), Vector3Matcher({-1, 0, 0}));
		CHECK_THAT(q * vector3(1, 0, 0), Vector3Matcher({0, 0, 1}));
	}
	{
		quaternion q = quaternion(direction(0.0f, 0.0f, 1.0f), HalfPi{});
		q.inverse();
		CHECK_THAT(q * vector3(0, 1, 0), Vector3Matcher({1, 0, 0}));
		CHECK_THAT(q * vector3(-1, 0, 0), Vector3Matcher({0, 1, 0}));
		CHECK_THAT(q * vector3(0, -1, 0), Vector3Matcher({-1, 0, 0}));
		CHECK_THAT(q * vector3(1, 0, 0), Vector3Matcher({0, -1, 0}));
	}
}

TEST_CASE("quaternion::inverse()_2", "[working][unittest][quaternion]")
{
	using Catch::Generators::range;

	const float dx = GENERATE(range(-1.0f, 1.0f, 0.2f));
	const float dy = GENERATE(range(-1.0f, 1.0f, 0.2f));
	const float dz = GENERATE(range(-1.0f, 1.0f, 0.2f));
	const float dw = GENERATE(range(-1.0f, 1.0f, 0.2f));

	if (dx != 0 && dy != 0 && dz != 0 && dw != 0)
	{
		auto q = quaternion{
		    dx,
		    dy,
		    dz,
		    dw,
		};
		q.normalize();

		const vector3 v{
		    GENERATE(-2.0f, 0.0f, 1.0f),
		    GENERATE(-2.0f, 0.0f, 1.0f),
		    GENERATE(-2.0f, 0.0f, 1.0f),
		};

		const auto v0 = q * v;

		q.inverse();
		q.inverse();

		const auto v1 = q * v;

		CHECK(v0 == v1);
	}
}
