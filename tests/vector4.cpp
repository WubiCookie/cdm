#include <common.hpp>

INFO_BEGIN(vector4)

TEST_CASE("vector4::vector4(std::array<float, 4>)",
          "[working][unittest][vector4]")
{
	const std::array<float, 4> a{0, 0, 0, 0};

	const vector4 v1 = a;

	CHECK(v1.x == 0);
	CHECK(v1.y == 0);
	CHECK(v1.z == 0);
	CHECK(v1.w == 0);

	const std::array<float, 4> b{-89453.6654f, 3.14159f, 5566656.66656f,
	                             -0.001f};

	const vector4 v2 = b;

	CHECK(v2.x == -89453.6654f);
	CHECK(v2.y == 3.14159f);
	CHECK(v2.z == 5566656.66656f);
	CHECK(v2.w == -0.001f);
}

TEST_CASE("vector4::vector4(float, float)", "[working][unittest][vector4]")
{
	const vector4 v1{0, 0, 0, 0};

	CHECK(v1.x == 0);
	CHECK(v1.y == 0);
	CHECK(v1.z == 0);

	const vector4 v2{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f};

	CHECK(v2.x == -89453.6654f);
	CHECK(v2.y == 3.14159f);
	CHECK(v2.z == 5566656.66656f);
	CHECK(v2.w == -0.001f);
}

TEST_CASE("vector4::norm_squared()", "[working][unittest][vector4]")
{
	CHECK(vector4{0, 0, 0, 0}.norm_squared() == 0);
	CHECK(vector4{1, 1, 1, 1}.norm_squared() == 4);
	CHECK(vector4{2, 2, 2, 2}.norm_squared() == 16);
	CHECK(vector4{-1, -1, 0, 0}.norm_squared() == 2);
	CHECK(vector4{-1, -1, -1, -1}.norm_squared() == 4);
	CHECK(vector4{-1, 1, 0, 0}.norm_squared() == 2);
}

TEST_CASE("vector4::norm()", "[working][unittest][vector4]")
{
	CHECK(vector4{0, 0, 0, 0}.norm() == 0);
	CHECK(vector4{1, 1, 0, 0}.norm() == Approx(sqrt2));
	CHECK(vector4{1, 1, 1, 1}.norm() == Approx(2.0f));
	CHECK(vector4{2, 2, 0, 0}.norm() == Approx(std::sqrt(8)));
	CHECK(vector4{-1, -1, 0, 0}.norm() == Approx(sqrt2));
	CHECK(vector4{-1, -1, -1, -1}.norm() == Approx(2.0f));
	CHECK(vector4{-1, 1, 0, 0}.norm() == Approx(sqrt2));
}

TEST_CASE("vector4::normalize()", "[working][unittest][vector4]")
{
	CHECK(vector4{1, 1, 1, 1}.normalize().norm() == Approx(1));
	CHECK(vector4{2, 2, 2, 2}.normalize().norm() == Approx(1));
	CHECK(vector4{-1, -1, -1, -1}.normalize().norm() == Approx(1));
	CHECK(vector4{-1, 1, 1, -1}.normalize().norm() == Approx(1));
	CHECK(vector4{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f}
	          .normalize()
	          .norm() == Approx(1));
}

TEST_CASE("vector4::get_normalized()", "[working][unittest][vector4]")
{
	CHECK(vector4{1, 1, 1, 1}.get_normalized().norm() == Approx(1));
	CHECK(vector4{2, 2, 2, 2}.get_normalized().norm() == Approx(1));
	CHECK(vector4{-1, -1, -1, -1}.get_normalized().norm() == Approx(1));
	CHECK(vector4{-1, 1, 1, -1}.get_normalized().norm() == Approx(1));
	CHECK(vector4{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f}
	          .get_normalized()
	          .norm() == Approx(1));

	CHECK(vector4{1, 1, 1, 1}.get_normalized() ==
	      vector4{1, 1, 1, 1}.normalize());
	CHECK(vector4{2, 2, 2, 2}.get_normalized() ==
	      vector4{2, 2, 2, 2}.normalize());
	CHECK(vector4{-1, -1, -1, -1}.get_normalized() ==
	      vector4{-1, -1, -1, -1}.normalize());
	CHECK(vector4{-1, 1, 1, -1}.get_normalized() ==
	      vector4{-1, 1, 1, -1}.normalize());
	CHECK(
	    vector4{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f}
	        .get_normalized() ==
	    vector4{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f}.normalize());
}

TEST_CASE("vector4::clamp(vector4, vector4)", "[working][unittest][vector4]")
{
	const vector4 zero{0, 0, 0, 0};
	const vector4 one{1, 1, 1, 1};

	CHECK(vector4{1, 1, 1, 1}.clamp(zero, one) == vector4{1, 1, 1, 1});
	CHECK(vector4{2, 2, 2, 2}.clamp(zero, one) == vector4{1, 1, 1, 1});
	CHECK(vector4{-1, -1, -1, -1}.clamp(zero, one) == vector4{0, 0, 0, 0});
	CHECK(vector4{-1, 1, 0, -1}.clamp(zero, one) == vector4{0, 1, 0, 0});
	CHECK(vector4{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f}.clamp(
	          {0, 0, 0, 0}, {1, 1, 1, 1}) == vector4{0, 1, 1, 0});
}

TEST_CASE("vector4::get_clamped(vector4, vector4)",
          "[working][unittest][vector4]")
{
	const vector4 zero{0, 0, 0, 0};
	const vector4 one{1, 1, 1, 1};

	CHECK(zero.get_clamped(zero, one) == zero);
	CHECK(one.get_clamped(zero, one) == one);
	CHECK(vector4{1, 1, 1, 1}.clamp(zero, one) ==
	      vector4{1, 1, 1, 1}.get_clamped(zero, one));
	CHECK(vector4{2, 2, 2, 2}.clamp(zero, one) ==
	      vector4{2, 2, 2, 2}.get_clamped(zero, one));
	CHECK(vector4{-1, -1, -1, -1}.clamp(zero, one) ==
	      vector4{-1, -1, -1, -1}.get_clamped(zero, one));
	CHECK(vector4{-1, 1, 0, -1}.clamp(zero, one) ==
	      vector4{-1, 1, 0, -1}.get_clamped(zero, one));
	CHECK(vector4{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f}.clamp(
	          zero, one) ==
	      vector4{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f}.get_clamped(
	          zero, one));
}

TEST_CASE("vector4::operator+(vector4)", "[working][unittest][vector4]")
{
	{
		const vector4 v1{0, 0, 0, 0};
		const vector4 v2{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f};
		const vector4 v3 = v1 + v2;

		CHECK(v3.x == -89453.6654f);
		CHECK(v3.y == 3.14159f);
		CHECK(v3.z == 5566656.66656f);
		CHECK(v3.w == -0.001f);
	}
	{
		const vector4 v1{1, 1, 1, 1};
		const vector4 v2{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f};
		const vector4 v3 = v1 + v2;

		CHECK(v3.x == -89453.6654f + 1);
		CHECK(v3.y == 3.14159f + 1);
		CHECK(v3.z == 5566656.66656f + 1);
		CHECK(v3.w == -0.001f + 1);
	}
}

TEST_CASE("vector4::operator-(vector4)", "[working][unittest][vector4]")
{
	{
		const vector4 v1{0, 0, 0, 0};
		const vector4 v2{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f};
		const vector4 v3 = v1 - v2;

		CHECK(v3.x == 89453.6654f);
		CHECK(v3.y == -3.14159f);
		CHECK(v3.z == -5566656.66656f);
		CHECK(v3.w == 0.001f);
	}
	{
		const vector4 v1{1, 1, 1, 1};
		const vector4 v2{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f};
		const vector4 v3 = v1 - v2;

		CHECK(v3.x == 1 - -89453.6654f);
		CHECK(v3.y == 1 - 3.14159f);
		CHECK(v3.z == 1 - 5566656.66656f);
		CHECK(v3.w == 1 - -0.001f);
	}
}

TEST_CASE("vector4::operator*(float)", "[working][unittest][vector4]")
{
	{
		const vector4 v1{0, 0, 0, 0};
		const vector4 v2 = v1 * 2.0f;

		CHECK(v2.x == 0);
		CHECK(v2.y == 0);
		CHECK(v2.z == 0);
		CHECK(v2.w == 0);
	}
	{
		const vector4 v1{0, 0, 0, 0};
		const vector4 v2 = 2.0f * v1;

		CHECK(v2.x == 0);
		CHECK(v2.y == 0);
		CHECK(v2.z == 0);
		CHECK(v2.w == 0);
	}
	{
		const vector4 v1{1, 1, 1, 1};
		const vector4 v2 = v1 * 2.0f;

		CHECK(v2.x == 2);
		CHECK(v2.y == 2);
		CHECK(v2.z == 2);
		CHECK(v2.w == 2);
	}
	{
		const vector4 v1{1, 1, 1, 1};
		const vector4 v2 = 2.0f * v1;

		CHECK(v2.x == 2);
		CHECK(v2.y == 2);
		CHECK(v2.z == 2);
		CHECK(v2.w == 2);
	}
	{
		const vector4 v1{3, -11, 555.545f, -0.001f};
		const vector4 v2 = v1 * 5.0f;

		CHECK(v2.x == v1.x * 5);
		CHECK(v2.y == v1.y * 5);
		CHECK(v2.z == v1.z * 5);
		CHECK(v2.w == v1.w * 5);
	}
	{
		const vector4 v1{3, -11, 555.545f, -0.001f};
		const vector4 v2 = 5.0f * v1;

		CHECK(v2.x == v1.x * 5);
		CHECK(v2.y == v1.y * 5);
		CHECK(v2.z == v1.z * 5);
		CHECK(v2.w == v1.w * 5);
	}
	{
		const vector4 v1{3, -11, 555.545f, -0.001f};
		const vector4 v2 = v1 * 0.0f;

		CHECK(v2.x == 0);
		CHECK(v2.y == 0);
		CHECK(v2.z == 0);
		CHECK(v2.w == 0);
	}
	{
		const vector4 v1{3, -11, 555.545f, -0.001f};
		const vector4 v2 = 0.0f * v1;

		CHECK(v2.x == 0);
		CHECK(v2.y == 0);
		CHECK(v2.z == 0);
		CHECK(v2.w == 0);
	}
}

TEST_CASE("vector4::operator/(float)", "[working][unittest][vector4]")
{
	{
		const vector4 v1{0, 0, 0, 0};
		const vector4 v2 = v1 / 2.0f;

		CHECK(v2.x == 0);
		CHECK(v2.y == 0);
		CHECK(v2.z == 0);
		CHECK(v2.w == 0);
	}
	{
		const vector4 v1{1, 1, 1, 1};
		const vector4 v2 = v1 / 2.0f;

		CHECK(v2.x == 0.5f);
		CHECK(v2.y == 0.5f);
		CHECK(v2.z == 0.5f);
		CHECK(v2.w == 0.5f);
	}
	{
		const vector4 v1{3.0f, -11.0f, 555.545f, -0.001f};
		const vector4 v2 = v1 / 5.0f;

		CHECK(v2.x == v1.x / 5.0f);
		CHECK(v2.y == v1.y / 5.0f);
		CHECK(v2.z == v1.z / 5.0f);
		CHECK(v2.w == v1.w / 5.0f);
	}
}

TEST_CASE("vector4::operator+=(vector4)", "[working][unittest][vector4]")
{
	{
		vector4 v1{0, 0, 0, 0};
		const vector4 v2{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f};
		v1 += v2;

		CHECK(v1.x == -89453.6654f);
		CHECK(v1.y == 3.14159f);
		CHECK(v1.z == 5566656.66656f);
		CHECK(v1.w == -0.001f);
	}
	{
		vector4 v1{1, 1, 1, 1};
		const vector4 v2{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f};
		v1 += v2;

		CHECK(v1.x == -89453.6654f + 1);
		CHECK(v1.y == 3.14159f + 1);
		CHECK(v1.z == 5566656.66656f + 1);
		CHECK(v1.w == -0.001f + 1);
	}
}

TEST_CASE("vector4::operator-=(vector4)", "[working][unittest][vector4]")
{
	{
		vector4 v1{0, 0, 0, 0};
		const vector4 v2{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f};
		v1 -= v2;

		CHECK(v1.x == 89453.6654f);
		CHECK(v1.y == -3.14159f);
		CHECK(v1.z == -5566656.66656f);
		CHECK(v1.w == 0.001f);
	}
	{
		vector4 v1{1, 1, 1, 1};
		const vector4 v2{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f};
		v1 -= v2;

		CHECK(v1.x == 1 - -89453.6654f);
		CHECK(v1.y == 1 - 3.14159f);
		CHECK(v1.z == 1 - 5566656.66656f);
		CHECK(v1.w == 1 - -0.001f);
	}
}

TEST_CASE("vector4::operator*=(float)", "[working][unittest][vector4]")
{
	{
		vector4 v1{0, 0, 0, 0};
		v1 *= 2.0f;

		CHECK(v1.x == 0);
		CHECK(v1.y == 0);
		CHECK(v1.z == 0);
		CHECK(v1.w == 0);
	}
	{
		vector4 v1{1, 1, 1, 1};
		v1 *= 2.0f;

		CHECK(v1.x == 2);
		CHECK(v1.y == 2);
		CHECK(v1.z == 2);
		CHECK(v1.w == 2);
	}
	{
		vector4 v1{3, -11, 555.545f, -0.001f};
		v1 *= 5.0f;

		CHECK(v1.x == 3.0f * 5.0f);
		CHECK(v1.y == -11.0f * 5.0f);
		CHECK(v1.z == 555.545f * 5.0f);
		CHECK(v1.w == -0.001f * 5.0f);
	}
	{
		vector4 v1{3, -11, 555.545f, -0.001f};
		v1 *= 0.0f;

		CHECK(v1.x == 0);
		CHECK(v1.y == 0);
		CHECK(v1.z == 0);
		CHECK(v1.w == 0);
	}
}

TEST_CASE("vector4::operator/=(float)", "[working][unittest][vector4]")
{
	{
		vector4 v1{0, 0, 0, 0};
		v1 /= 2.0f;

		CHECK(v1.x == 0);
		CHECK(v1.y == 0);
		CHECK(v1.z == 0);
		CHECK(v1.w == 0);
	}
	{
		vector4 v1{1, 1, 1, 1};
		v1 /= 2.0f;

		CHECK(v1.x == 0.5f);
		CHECK(v1.y == 0.5f);
		CHECK(v1.z == 0.5f);
		CHECK(v1.w == 0.5f);
	}
	{
		vector4 v1{3.0f, -11.0f, 555.545f, -0.001f};
		v1 /= 5.0f;

		CHECK(v1.x == 3.0f / 5.0f);
		CHECK(v1.y == -11.0f / 5.0f);
		CHECK(v1.z == Approx(555.545f / 5.0f));
		CHECK(v1.w == -0.001f / 5.0f);
	}
}

TEST_CASE("vector4::operator-()", "[working][unittest][vector4]")
{
	const vector4 v1{0, 0, 0, 0};
	const vector4 v12 = -v1;

	CHECK(v12.x == -0);
	CHECK(v12.y == -0);
	CHECK(v12.z == -0);
	CHECK(v12.w == -0);

	const vector4 v2{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f};
	const vector4 v22 = -v2;

	CHECK(v22.x == 89453.6654f);
	CHECK(v22.y == -3.14159f);
	CHECK(v22.z == -5566656.66656f);
	CHECK(v22.w == 0.001f);
}

TEST_CASE("vector4::operator==()", "[working][unittest][vector4]")
{
	const vector4 v1{0, 0, 0, 0};
	const vector4 v2{0, 0, 0, 0};

	CHECK(v1 == v2);
	CHECK(v1.x == v2.x);
	CHECK(v1.y == v2.y);
	CHECK(v1.z == v2.z);
	CHECK(v1.w == v2.w);

	const vector4 v3{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f};
	const vector4 v4{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f};

	CHECK(v3 == v4);
	CHECK(v3.x == v4.x);
	CHECK(v3.y == v4.y);
	CHECK(v3.z == v4.z);
	CHECK(v3.w == v4.w);

	CHECK_FALSE(v1 == v3);
}

TEST_CASE("vector4::operator!=()", "[working][unittest][vector4]")
{
	const vector4 v1{0, 0, 0, 0};
	const vector4 v2{0, 0, 0, 0};

	CHECK(!(v1 != v2));
	CHECK_FALSE(v1 != v2);
	CHECK(v1.x == v2.x);
	CHECK(v1.y == v2.y);
	CHECK(v1.z == v2.z);
	CHECK(v1.w == v2.w);

	const vector4 v3{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f};
	const vector4 v4{-89453.6654f, 3.14159f, 5566656.66656f, -0.001f};

	CHECK(!(v3 != v4));
	CHECK_FALSE(v3 != v4);
	CHECK(v3.x == v4.x);
	CHECK(v3.y == v4.y);
	CHECK(v3.z == v4.z);
	CHECK(v3.w == v4.w);

	CHECK(v1 != v3);
}

TEST_CASE("dot(vector4, vector4)", "[working][unittest][vector4]")
{
	CHECK(dot(vector4{0, 0, 0, 0}, {0, 0, 0, 0}) == 0);
	CHECK(dot(vector4{1, 0, 0, 0}, {0, 1, 0, 0}) == 0);
	CHECK(dot(vector4{1, 0, 0, 0}, {1, 0, 0, 0}) == 1);
	CHECK(dot(vector4{0, 1, 0, 0}, {0, 1, 0, 0}) == 1);
	CHECK(dot(vector4{1, 1, 0, 0}, {1, 0, 0, 0}) == 1);
	CHECK(dot(vector4{0.5f, 1, 0, 0}, {1, 0, 0, 0}) == 0.5f);
	CHECK(dot(vector4{1, 0, 0, 0}, {0.5f, 1, 0, 0}) == 0.5f);
}

TEST_CASE("distance_between(vector4, vector4)", "[working][unittest][vector4]")
{
	CHECK(distance_between({0, 0, 0, 0}, vector4{0, 0, 0, 0}) == 0);
	CHECK(distance_between({1, 0, 0, 0}, vector4{0, 1, 0, 0}) ==
	      Approx(sqrt2));
	CHECK(distance_between({0, 1, 0, 0}, vector4{1, 0, 0, 0}) ==
	      Approx(sqrt2));
	CHECK(distance_between({1, 0, 0, 0}, vector4{1, 0, 0, 0}) == 0);
	CHECK(distance_between({0, 1, 0, 0}, vector4{0, 1, 0, 0}) == 0);
	CHECK(distance_between({1, 0, 0, 0},
	                       vector4{0, 0.5f, 0, 0}.get_normalized()) ==
	      Approx(sqrt2));
	CHECK(distance_between({1, 0, 0, 0}, vector4{-1, 0, 0, 0}) == 2);
	CHECK(distance_between({2, 0, 0, 0}, vector4{-1, 0, 0, 0}) == 3);
}

TEST_CASE("distance_squared_between(vector4, vector4)",
          "[working][unittest][vector4]")
{
	CHECK(distance_squared_between({0, 0, 0, 0}, vector4{0, 0, 0, 0}) == 0);
	CHECK(distance_squared_between({1, 0, 0, 0}, vector4{0, 1, 0, 0}) == 2);
	CHECK(distance_squared_between({0, 1, 0, 0}, vector4{1, 0, 0, 0}) == 2);
	CHECK(distance_squared_between({1, 0, 0, 0}, vector4{1, 0, 0, 0}) == 0);
	CHECK(distance_squared_between({0, 1, 0, 0}, vector4{0, 1, 0, 0}) == 0);
	CHECK(distance_squared_between(
	          {1, 0, 0, 0}, vector4{0, 0.5f, 0, 0}.get_normalized()) == 2);
	CHECK(distance_squared_between({1, 0, 0, 0}, vector4{-1, 0, 0, 0}) == 4);
	CHECK(distance_squared_between({2, 0, 0, 0}, vector4{-1, 0, 0, 0}) == 9);
}

TEST_CASE("element_wise_min(vector4, vector4)", "[working][unittest][vector4]")
{
	CHECK(element_wise_min({-89453.6654f, 3.14159f, 5566656.66656f, -0.001f},
	                       vector4{0, 0, 0, 0}) ==
	      vector4{-89453.6654f, 0, 0, -0.001f});
}

TEST_CASE("element_wise_max(vector4, vector4)", "[working][unittest][vector4]")
{
	CHECK(element_wise_max({-89453.6654f, 3.14159f, 5566656.66656f, -0.001f},
	                       vector4{0, 0, 0, 0}) ==
	      vector4{0, 3.14159f, 5566656.66656f, 0});
}

INFO_END(vector4)
