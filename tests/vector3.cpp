#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <common.hpp>

TEST_CASE("vector3::vector3(std::array<float, 3>)",
          "[working][unittest][vector3]")
{
	const std::array<float, 3> a{0, 0, 0};

	const vector3 v1 = a;

	CHECK(v1.x == 0);
	CHECK(v1.y == 0);
	CHECK(v1.z == 0);

	const std::array<float, 3> b{-89453.6654f, 3.14159f, 5566656.66656f};

	const vector3 v2 = b;

	CHECK(v2.x == -89453.6654f);
	CHECK(v2.y == 3.14159f);
	CHECK(v2.z == 5566656.66656f);
}

TEST_CASE("vector3::vector3(float, float)", "[working][unittest][vector3]")
{
	const vector3 v1{0, 0, 0};

	CHECK(v1.x == 0);
	CHECK(v1.y == 0);
	CHECK(v1.z == 0);

	const vector3 v2{-89453.6654f, 3.14159f, 5566656.66656f};

	CHECK(v2.x == -89453.6654f);
	CHECK(v2.y == 3.14159f);
	CHECK(v2.z == 5566656.66656f);
}

TEST_CASE("vector3::norm_squared()", "[working][unittest][vector3]")
{
	CHECK(vector3{1, 1, 1}.norm_squared() == 3);
	CHECK(vector3{2, 2, 2}.norm_squared() == 12);
	CHECK(vector3{-1, -1, 0}.norm_squared() == 2);
	CHECK(vector3{-1, -1, -1}.norm_squared() == 3);
	CHECK(vector3{-1, 1, 0}.norm_squared() == 2);
}

TEST_CASE("vector3::norm()", "[working][unittest][vector3]")
{
	CHECK(vector3{1, 1, 0}.norm() == Approx(sqrt2));
	CHECK(vector3{1, 1, 1}.norm() == Approx(sqrt3));
	CHECK(vector3{2, 2, 0}.norm() == Approx(std::sqrt(8)));
	CHECK(vector3{-1, -1, 0}.norm() == Approx(sqrt2));
	CHECK(vector3{-1, -1, -1}.norm() == Approx(sqrt3));
	CHECK(vector3{-1, 1, 0}.norm() == Approx(sqrt2));
}

TEST_CASE("vector3::normalize()", "[working][unittest][vector3]")
{
	CHECK(vector3{1, 1, 1}.normalize().norm() == Approx(1));
	CHECK(vector3{2, 2, 2}.normalize().norm() == Approx(1));
	CHECK(vector3{-1, -1, -1}.normalize().norm() == Approx(1));
	CHECK(vector3{-1, 1, 1}.normalize().norm() == Approx(1));
	CHECK(vector3{-89453.6654f, 3.14159f, 5566656.66656f}.normalize().norm() ==
	      Approx(1));
}

TEST_CASE("vector3::get_normalized()", "[working][unittest][vector3]")
{
	CHECK(vector3{1, 1, 1}.get_normalized().norm() == Approx(1));
	CHECK(vector3{2, 2, 2}.get_normalized().norm() == Approx(1));
	CHECK(vector3{-1, -1, -1}.get_normalized().norm() == Approx(1));
	CHECK(vector3{-1, 1, 1}.get_normalized().norm() == Approx(1));
	CHECK(vector3{-89453.6654f, 3.14159f, 5566656.66656f}
	          .get_normalized()
	          .norm() == Approx(1));

	CHECK(vector3{1, 1, 1}.get_normalized() == vector3{1, 1, 1}.normalize());
	CHECK(vector3{2, 2, 2}.get_normalized() == vector3{2, 2, 2}.normalize());
	CHECK(vector3{-1, -1, -1}.get_normalized() ==
	      vector3{-1, -1, -1}.normalize());
	CHECK(vector3{-1, 1, 1}.get_normalized() == vector3{-1, 1, 1}.normalize());
	CHECK(vector3{-89453.6654f, 3.14159f, 5566656.66656f}.get_normalized() ==
	      vector3{-89453.6654f, 3.14159f, 5566656.66656f}.normalize());
}

TEST_CASE("vector3::clamp(vector3, vector3)", "[working][unittest][vector3]")
{
	CHECK(vector3{1, 1, 1}.clamp({0, 0, 0}, {1, 1, 1}) == vector3{1, 1, 1});
	CHECK(vector3{2, 2, 2}.clamp({0, 0, 0}, {1, 1, 1}) == vector3{1, 1, 1});
	CHECK(vector3{-1, -1, -1}.clamp({0, 0, 0}, {1, 1, 1}) == vector3{0, 0, 0});
	CHECK(vector3{-1, 1, 0}.clamp({0, 0, 0}, {1, 1, 1}) == vector3{0, 1, 0});
	CHECK(vector3{-89453.6654f, 3.14159f, 5566656.66656f}.clamp(
	          {0, 0, 0}, {1, 1, 1}) == vector3{0, 1, 1});
}

TEST_CASE("vector3::get_clamped(vector3, vector3)",
          "[working][unittest][vector3]")
{
	CHECK(vector3{1, 1, 1}.clamp({0, 0, 0}, {1, 1, 1}) ==
	      vector3{1, 1, 1}.get_clamped({0, 0, 0}, {1, 1, 1}));
	CHECK(vector3{2, 2, 2}.clamp({0, 0, 0}, {1, 1, 1}) ==
	      vector3{2, 2, 2}.get_clamped({0, 0, 0}, {1, 1, 1}));
	CHECK(vector3{-1, -1, -1}.clamp({0, 0, 0}, {1, 1, 1}) ==
	      vector3{-1, -1, -1}.get_clamped({0, 0, 0}, {1, 1, 1}));
	CHECK(vector3{-1, 1, 0}.clamp({0, 0, 0}, {1, 1, 1}) ==
	      vector3{-1, 1, 0}.get_clamped({0, 0, 0}, {1, 1, 1}));
	CHECK(vector3{-89453.6654f, 3.14159f, 5566656.66656f}.clamp({0, 0, 0},
	                                                            {1, 1, 1}) ==
	      vector3{-89453.6654f, 3.14159f, 5566656.66656f}.get_clamped(
	          {0, 0, 0}, {1, 1, 1}));
}

TEST_CASE("vector3::operator+(vector3)", "[working][unittest][vector3]")
{
	{
		const vector3 v1{0, 0, 0};
		const vector3 v2{-89453.6654f, 3.14159f, 5566656.66656f};
		const vector3 v3 = v1 + v2;

		CHECK(v3.x == -89453.6654f);
		CHECK(v3.y == 3.14159f);
		CHECK(v3.z == 5566656.66656f);
	}
	{
		const vector3 v1{1, 1, 1};
		const vector3 v2{-89453.6654f, 3.14159f, 5566656.66656f};
		const vector3 v3 = v1 + v2;

		CHECK(v3.x == -89453.6654f + 1);
		CHECK(v3.y == 3.14159f + 1);
		CHECK(v3.z == 5566656.66656f + 1);
	}
}

TEST_CASE("vector3::operator-(vector3)", "[working][unittest][vector3]")
{
	{
		const vector3 v1{0, 0, 0};
		const vector3 v2{-89453.6654f, 3.14159f, 5566656.66656f};
		const vector3 v3 = v1 - v2;

		CHECK(v3.x == 89453.6654f);
		CHECK(v3.y == -3.14159f);
		CHECK(v3.z == -5566656.66656f);
	}
	{
		const vector3 v1{1, 1, 1};
		const vector3 v2{-89453.6654f, 3.14159f, 5566656.66656f};
		const vector3 v3 = v1 - v2;

		CHECK(v3.x == 1 - -89453.6654f);
		CHECK(v3.y == 1 - 3.14159f);
		CHECK(v3.z == 1 - 5566656.66656f);
	}
}

TEST_CASE("vector3::operator*(float)", "[working][unittest][vector3]")
{
	{
		const vector3 v1{0, 0, 0};
		const vector3 v2 = v1 * 2.0f;

		CHECK(v2.x == 0);
		CHECK(v2.y == 0);
		CHECK(v2.z == 0);
	}
	{
		const vector3 v1{0, 0, 0};
		const vector3 v2 = 2.0f * v1;

		CHECK(v2.x == 0);
		CHECK(v2.y == 0);
		CHECK(v2.z == 0);
	}
	{
		const vector3 v1{1, 1, 1};
		const vector3 v2 = v1 * 2.0f;

		CHECK(v2.x == 2);
		CHECK(v2.y == 2);
		CHECK(v2.z == 2);
	}
	{
		const vector3 v1{1, 1, 1};
		const vector3 v2 = 2.0f * v1;

		CHECK(v2.x == 2);
		CHECK(v2.y == 2);
		CHECK(v2.z == 2);
	}
	{
		const vector3 v1{3, -11, 555.545f};
		const vector3 v2 = v1 * 5.0f;

		CHECK(v2.x == v1.x * 5);
		CHECK(v2.y == v1.y * 5);
		CHECK(v2.z == v1.z * 5);
	}
	{
		const vector3 v1{3, -11, 555.545f};
		const vector3 v2 = 5.0f * v1;

		CHECK(v2.x == v1.x * 5);
		CHECK(v2.y == v1.y * 5);
		CHECK(v2.z == v1.z * 5);
	}
	{
		const vector3 v1{3, -11, 555.545f};
		const vector3 v2 = v1 * 0.0f;

		CHECK(v2.x == 0);
		CHECK(v2.y == 0);
		CHECK(v2.z == 0);
	}
	{
		const vector3 v1{3, -11, 555.545f};
		const vector3 v2 = 0.0f * v1;

		CHECK(v2.x == 0);
		CHECK(v2.y == 0);
		CHECK(v2.z == 0);
	}
}

TEST_CASE("vector3::operator/(float)", "[working][unittest][vector3]")
{
	{
		const vector3 v1{0, 0, 0};
		const vector3 v2 = v1 / 2.0f;

		CHECK(v2.x == 0);
		CHECK(v2.y == 0);
		CHECK(v2.z == 0);
	}
	{
		const vector3 v1{1, 1, 1};
		const vector3 v2 = v1 / 2.0f;

		CHECK(v2.x == 0.5f);
		CHECK(v2.y == 0.5f);
		CHECK(v2.z == 0.5f);
	}
	{
		const vector3 v1{3.0f, -11.0f, 555.545f};
		const vector3 v2 = v1 / 5.0f;

		CHECK(v2.x == v1.x / 5.0f);
		CHECK(v2.y == v1.y / 5.0f);
		CHECK(v2.z == v1.z / 5.0f);
	}
}

TEST_CASE("vector3::operator+=(vector3)", "[working][unittest][vector3]")
{
	{
		vector3 v1{0, 0, 0};
		const vector3 v2{-89453.6654f, 3.14159f, 5566656.66656f};
		v1 += v2;

		CHECK(v1.x == -89453.6654f);
		CHECK(v1.y == 3.14159f);
		CHECK(v1.z == 5566656.66656f);
	}
	{
		vector3 v1{1, 1, 1};
		const vector3 v2{-89453.6654f, 3.14159f, 5566656.66656f};
		v1 += v2;

		CHECK(v1.x == -89453.6654f + 1);
		CHECK(v1.y == 3.14159f + 1);
		CHECK(v1.z == 5566656.66656f + 1);
	}
}

TEST_CASE("vector3::operator-=(vector3)", "[working][unittest][vector3]")
{
	{
		vector3 v1{0, 0, 0};
		const vector3 v2{-89453.6654f, 3.14159f, 5566656.66656f};
		v1 -= v2;

		CHECK(v1.x == 89453.6654f);
		CHECK(v1.y == -3.14159f);
		CHECK(v1.z == -5566656.66656f);
	}
	{
		vector3 v1{1, 1, 1};
		const vector3 v2{-89453.6654f, 3.14159f, 5566656.66656f};
		v1 -= v2;

		CHECK(v1.x == 1 - -89453.6654f);
		CHECK(v1.y == 1 - 3.14159f);
		CHECK(v1.z == 1 - 5566656.66656f);
	}
}

TEST_CASE("vector3::operator*=(float)", "[working][unittest][vector3]")
{
	{
		vector3 v1{0, 0, 0};
		v1 *= 2.0f;

		CHECK(v1.x == 0);
		CHECK(v1.y == 0);
		CHECK(v1.z == 0);
	}
	{
		vector3 v1{1, 1, 1};
		v1 *= 2.0f;

		CHECK(v1.x == 2);
		CHECK(v1.y == 2);
		CHECK(v1.z == 2);
	}
	{
		vector3 v1{3, -11, 555.545f};
		v1 *= 5.0f;

		CHECK(v1.x == 3.0f * 5.0f);
		CHECK(v1.y == -11.0f * 5.0f);
		CHECK(v1.z == 555.545f * 5.0f);
	}
	{
		vector3 v1{3, -11, 555.545f};
		v1 *= 0.0f;

		CHECK(v1.x == 0);
		CHECK(v1.y == 0);
		CHECK(v1.z == 0);
	}
}

TEST_CASE("vector3::operator/=(float)", "[working][unittest][vector3]")
{
	{
		vector3 v1{0, 0, 0};
		v1 /= 2.0f;

		CHECK(v1.x == 0);
		CHECK(v1.y == 0);
		CHECK(v1.z == 0);
	}
	{
		vector3 v1{1, 1, 1};
		v1 /= 2.0f;

		CHECK(v1.x == 0.5f);
		CHECK(v1.y == 0.5f);
		CHECK(v1.z == 0.5f);
	}
	{
		vector3 v1{3.0f, -11.0f, 555.545f};
		v1 /= 5.0f;

		CHECK(v1.x == 3.0f / 5.0f);
		CHECK(v1.y == -11.0f / 5.0f);
		CHECK(v1.z == Approx(555.545f / 5.0f));
	}
}

TEST_CASE("vector3::operator-()", "[working][unittest][vector3]")
{
	const vector3 v1{0, 0, 0};
	const vector3 v12 = -v1;

	CHECK(v12.x == -0);
	CHECK(v12.y == -0);
	CHECK(v12.z == -0);

	const vector3 v2{-89453.6654f, 3.14159f, 5566656.66656f};
	const vector3 v22 = -v2;

	CHECK(v22.x == 89453.6654f);
	CHECK(v22.y == -3.14159f);
	CHECK(v22.z == -5566656.66656f);
}

TEST_CASE("vector3::vector3::operator==()", "[working][unittest][vector3]")
{
	const vector3 v1{0, 0, 0};
	const vector3 v2{0, 0, 0};

	CHECK(v1 == v2);
	CHECK(v1.x == v2.x);
	CHECK(v1.y == v2.y);
	CHECK(v1.z == v2.z);

	const vector3 v3{-89453.6654f, 3.14159f, 5566656.66656f};
	const vector3 v4{-89453.6654f, 3.14159f, 5566656.66656f};

	CHECK(v3 == v4);
	CHECK(v3.x == v4.x);
	CHECK(v3.y == v4.y);
	CHECK(v3.z == v4.z);

	CHECK_FALSE(v1 == v3);
}

TEST_CASE("vector3::vector3::operator!=()", "[working][unittest][vector3]")
{
	const vector3 v1{0, 0, 0};
	const vector3 v2{0, 0, 0};

	CHECK(!(v1 != v2));
	CHECK_FALSE(v1 != v2);
	CHECK(v1.x == v2.x);
	CHECK(v1.y == v2.y);
	CHECK(v1.z == v2.z);

	const vector3 v3{-89453.6654f, 3.14159f, 5566656.66656f};
	const vector3 v4{-89453.6654f, 3.14159f, 5566656.66656f};

	CHECK(!(v3 != v4));
	CHECK_FALSE(v3 != v4);
	CHECK(v3.x == v4.x);
	CHECK(v3.y == v4.y);
	CHECK(v3.z == v4.z);

	CHECK(v1 != v3);
}

TEST_CASE("vector3::dot(vector3, vector3)", "[working][unittest][vector3]")
{
	CHECK(dot(vector3{0, 0, 0}, vector3{0, 0, 0}) == 0);
	CHECK(dot(vector3{1, 0, 0}, vector3{0, 1, 0}) == 0);
	CHECK(dot(vector3{1, 0, 0}, vector3{1, 0, 0}) == 1);
	CHECK(dot(vector3{0, 1, 0}, vector3{0, 1, 0}) == 1);
	CHECK(dot(vector3{1, 1, 0}, vector3{1, 0, 0}) == 1);
	CHECK(dot(vector3{0.5f, 1, 0}, vector3{1, 0, 0}) == 0.5f);
	CHECK(dot(vector3{1, 0, 0}, vector3{0.5f, 1, 0}) == 0.5f);
}

TEST_CASE("vector3::cross(vector3, vector3)", "[working][unittest][vector3]")
{
	CHECK(cross(vector3{0, 0, 0}, vector3{0, 0, 0}) == vector3{0, 0, 0});

	CHECK(cross(vector3{1, 0, 0}, vector3{0, 1, 0}) == vector3{0, 0, 1});
	CHECK(cross(vector3{0, 1, 0}, vector3{0, 0, 1}) == vector3{1, 0, 0});
	CHECK(cross(vector3{0, 0, 1}, vector3{1, 0, 0}) == vector3{0, 1, 0});

	CHECK(cross(vector3{0, 1, 0}, vector3{1, 0, 0}) == vector3{0, 0, -1});
	CHECK(cross(vector3{0, 0, 1}, vector3{0, 1, 0}) == vector3{-1, 0, 0});
	CHECK(cross(vector3{1, 0, 0}, vector3{0, 0, 1}) == vector3{0, -1, 0});
}

TEST_CASE("vector3::distance_between(vector3, vector3)",
          "[working][unittest][vector3]")
{
	CHECK(distance_between({0, 0, 0}, vector3{0, 0, 0}) == 0);
	CHECK(distance_between({1, 0, 0}, vector3{0, 1, 0}) == Approx(sqrt2));
	CHECK(distance_between({0, 1, 0}, vector3{1, 0, 0}) == Approx(sqrt2));
	CHECK(distance_between({1, 0, 0}, vector3{1, 0, 0}) == 0);
	CHECK(distance_between({0, 1, 0}, vector3{0, 1, 0}) == 0);
	CHECK(distance_between({1, 0, 0}, vector3{0, 0.5f, 0}.get_normalized()) ==
	      Approx(sqrt2));
	CHECK(distance_between({1, 0, 0}, vector3{-1, 0, 0}) == 2);
	CHECK(distance_between({2, 0, 0}, vector3{-1, 0, 0}) == 3);
}

TEST_CASE("vector3::distance_squared_between(vector3, vector3)",
          "[working][unittest][vector3]")
{
	CHECK(distance_squared_between({0, 0, 0}, vector3{0, 0, 0}) == 0);
	CHECK(distance_squared_between({1, 0, 0}, vector3{0, 1, 0}) == 2);
	CHECK(distance_squared_between({0, 1, 0}, vector3{1, 0, 0}) == 2);
	CHECK(distance_squared_between({1, 0, 0}, vector3{1, 0, 0}) == 0);
	CHECK(distance_squared_between({0, 1, 0}, vector3{0, 1, 0}) == 0);
	CHECK(distance_squared_between({1, 0, 0},
	                               vector3{0, 0.5f, 0}.get_normalized()) == 2);
	CHECK(distance_squared_between({1, 0, 0}, vector3{-1, 0, 0}) == 4);
	CHECK(distance_squared_between({2, 0, 0}, vector3{-1, 0, 0}) == 9);
}

TEST_CASE("vector3::angle_between(vector3, vector3)",
          "[working][unittest][vector3]")
{
	CHECK(angle_between({0, 0, 0}, vector3{0, 0, 0}) == 0_rad);
	CHECK(angle_between({1, 0, 0}, vector3{0, 0, 0}) == 0_rad);
	CHECK(angle_between({0, 0, 0}, vector3{1, 0, 0}) == 0_rad);

	CHECK(angle_between({1, 0, 0}, vector3{0, 1, 0}) == radian(90_deg));
	CHECK(angle_between({1, 0, 0}, vector3{1, 1, 0}) == radian(45_deg));
	CHECK(angle_between({0, 1, 0}, vector3{1, 0, 0}) == radian(90_deg));
	CHECK(angle_between({1, 0, 0}, vector3{1, 0, 0}) == radian(0_deg));
	CHECK(angle_between({0, 1, 0}, vector3{0, 1, 0}) == radian(0_deg));
	CHECK(angle_between({1, 0, 0}, vector3{0, 0.5f, 0}.get_normalized()) ==
	      1_pi / 2.0f);
	CHECK(angle_between({1, 0, 0}, vector3{-1, 0, 0}) == 1_pi);
	CHECK(angle_between({2, 0, 0}, vector3{-1, 0, 0}) == 1_pi);

	CHECK(angle_between({1, 0, 0}, vector3{0, 0, 1}) == radian(90_deg));
	CHECK(angle_between({1, 0, 0}, vector3{1, 0, 1}) == radian(45_deg));
	CHECK(angle_between({0, 0, 1}, vector3{1, 0, 0}) == radian(90_deg));
	CHECK(angle_between({1, 0, 0}, vector3{1, 0, 0}) == radian(0_deg));
	CHECK(angle_between({0, 0, 1}, vector3{0, 0, 1}) == radian(0_deg));
	CHECK(angle_between({1, 0, 0}, vector3{0, 0, 0.5f}.get_normalized()) ==
	      1_pi / 2.0f);
	CHECK(angle_between({1, 0, 0}, vector3{-1, 0, 0}) == 1_pi);
	CHECK(angle_between({2, 0, 0}, vector3{-1, 0, 0}) == 1_pi);

	CHECK(angle_between({0, 1, 0}, vector3{0, 0, 1}) == radian(90_deg));
	CHECK(angle_between({0, 1, 0}, vector3{0, 1, 1}) == radian(45_deg));
	CHECK(angle_between({0, 0, 1}, vector3{0, 1, 0}) == radian(90_deg));
	CHECK(angle_between({0, 1, 0}, vector3{0, 1, 0}) == radian(0_deg));
	CHECK(angle_between({0, 0, 1}, vector3{0, 0, 1}) == radian(0_deg));
	CHECK(angle_between({0, 1, 0}, vector3{0, 0, 0.5f}.get_normalized()) ==
	      1_pi / 2.0f);
	CHECK(angle_between({0, 1, 0}, vector3{0, -1, 0}) == 1_pi);
	CHECK(angle_between({0, 2, 0}, vector3{0, -1, 0}) == 1_pi);
}

TEST_CASE("vector3::element_wise_min(vector3, vector3)",
          "[working][unittest][vector3]")
{
	CHECK(element_wise_min({-89453.6654f, 3.14159f, 5566656.66656f},
	                       vector3{0, 0, 0}) == vector3{-89453.6654f, 0, 0});
}

TEST_CASE("vector3::element_wise_max(vector3, vector3)",
          "[working][unittest][vector3]")
{
	CHECK(element_wise_max({-89453.6654f, 3.14159f, 5566656.66656f},
	                       vector3{0, 0, 0}) ==
	      vector3{0, 3.14159f, 5566656.66656f});
}

TEST_CASE("vector3::angle_around_axis(vector3, vector3, direction)",
          "[working][unittest][vector3]")
{
	CHECK(angle_around_axis(vector3{1, 0, 0}, vector3{0, 1, 0},
	                        direction3::posZ()) == radian(90_deg));
	CHECK(angle_around_axis(vector3{1, 0, 0}, vector3{0, 1, 0},
	                        direction3::negZ()) == radian(-90_deg));
	CHECK(angle_around_axis(vector3{0, 1, 0}, vector3{1, 0, 0},
	                        direction3::posZ()) == radian(-90_deg));
	CHECK(angle_around_axis(vector3{0, 1, 0}, vector3{1, 0, 0},
	                        direction3::negZ()) == radian(90_deg));
	CHECK(angle_around_axis(vector3{1, 0, 0}, vector3{1, 0, 0},
	                        direction3::posZ()) == radian(0_deg));
	CHECK(angle_around_axis(vector3{1, 0, 0}, vector3{-1, 0, 0},
	                        direction3::posZ()) == radian(180_deg));
	CHECK(angle_around_axis(vector3{1, 0, 0}, vector3{0, -1, 0},
	                        direction3::posZ()) == radian(-90_deg));
	CHECK(angle_around_axis(vector3{1, 0, 0}, vector3{1, 1, 0},
	                        direction3::posZ()) == radian(45_deg));
}
