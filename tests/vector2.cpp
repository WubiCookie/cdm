#include <common.hpp>

TEST_CASE("vector2::vector2(std::array<float, 2>)", "[working][unittest]")
{
	const std::array<float, 2> a{0.0f, 0.0f};

	const vector2 v1 = a;

	CHECK(v1.x == 0.0f);
	CHECK(v1.y == 0.0f);

	const std::array<float, 2> b{-89453.6654f, 3.14159f};

	const vector2 v2 = b;

	CHECK(v2.x == -89453.6654f);
	CHECK(v2.y == 3.14159f);
}

TEST_CASE("vector2::vector2(float, float)", "[working][unittest]")
{
	const vector2 v1{0.0f, 0.0f};

	CHECK(v1.x == 0.0f);
	CHECK(v1.y == 0.0f);

	const vector2 v2{-89453.6654f, 3.14159f};

	CHECK(v2.x == -89453.6654f);
	CHECK(v2.y == 3.14159f);
}

TEST_CASE("vector2::norm_squared()", "[working][unittest]")
{
	CHECK(vector2{1.0f, 1.0f}.norm_squared() == 2.0f);
	CHECK(vector2{2.0f, 2.0f}.norm_squared() == 8.0f);
	CHECK(vector2{-1.0f, -1.0f}.norm_squared() == 2.0f);
	CHECK(vector2{-1.0f, 1.0f}.norm_squared() == 2.0f);
}

TEST_CASE("vector2::norm()", "[working][unittest]")
{
	CHECK(vector2{1.0f, 1.0f}.norm() == Approx(sqrt2));
	CHECK(vector2{2.0f, 2.0f}.norm() == Approx(std::sqrt(8.0f)));
	CHECK(vector2{-1.0f, -1.0f}.norm() == Approx(sqrt2));
	CHECK(vector2{-1.0f, 1.0f}.norm() == Approx(sqrt2));
}

TEST_CASE("vector2::normalize()", "[working][unittest]")
{
	CHECK(vector2{1.0f, 1.0f}.normalize().norm() == Approx(1.0f));
	CHECK(vector2{2.0f, 2.0f}.normalize().norm() == Approx(1.0f));
	CHECK(vector2{-1.0f, -1.0f}.normalize().norm() == Approx(1.0f));
	CHECK(vector2{-1.0f, 1.0f}.normalize().norm() == Approx(1.0f));
	CHECK(vector2{-89453.6654f, 3.14159f}.normalize().norm() == Approx(1.0f));
}

TEST_CASE("vector2::get_normalized()", "[working][unittest]")
{
	CHECK(vector2{1.0f, 1.0f}.get_normalized().norm() == Approx(1.0f));
	CHECK(vector2{2.0f, 2.0f}.get_normalized().norm() == Approx(1.0f));
	CHECK(vector2{-1.0f, -1.0f}.get_normalized().norm() == Approx(1.0f));
	CHECK(vector2{-1.0f, 1.0f}.get_normalized().norm() == Approx(1.0f));
	CHECK(vector2{-89453.6654f, 3.14159f}.get_normalized().norm() ==
	      Approx(1.0f));

	CHECK(vector2{1.0f, 1.0f}.get_normalized() ==
	      vector2{1.0f, 1.0f}.normalize());
	CHECK(vector2{2.0f, 2.0f}.get_normalized() ==
	      vector2{2.0f, 2.0f}.normalize());
	CHECK(vector2{-1.0f, -1.0f}.get_normalized() ==
	      vector2{-1.0f, -1.0f}.normalize());
	CHECK(vector2{-1.0f, 1.0f}.get_normalized() ==
	      vector2{-1.0f, 1.0f}.normalize());
	CHECK(vector2{-89453.6654f, 3.14159f}.get_normalized() ==
	      vector2{-89453.6654f, 3.14159f}.normalize());
}

TEST_CASE("vector2::clamp(vector2, vector2)", "[working][unittest]")
{
	CHECK(vector2{1.0f, 1.0f}.clamp({0.0f, 0.0f}, {1.0f, 1.0f}) ==
	      vector2{1.0f, 1.0f});
	CHECK(vector2{2.0f, 2.0f}.clamp({0.0f, 0.0f}, {1.0f, 1.0f}) ==
	      vector2{1.0f, 1.0f});
	CHECK(vector2{-1.0f, -1.0f}.clamp({0.0f, 0.0f}, {1.0f, 1.0f}) ==
	      vector2{0.0f, 0.0f});
	CHECK(vector2{-1.0f, 1.0f}.clamp({0.0f, 0.0f}, {1.0f, 1.0f}) ==
	      vector2{0.0f, 1.0f});
	CHECK(vector2{-89453.6654f, 3.14159f}.clamp({0.0f, 0.0f}, {1.0f, 1.0f}) ==
	      vector2{0.0f, 1.0f});
}

TEST_CASE("vector2::get_clamped(vector2, vector2)", "[working][unittest]")
{
	CHECK(vector2{1.0f, 1.0f}.clamp({0.0f, 0.0f}, {1.0f, 1.0f}) ==
	      vector2{1.0f, 1.0f}.get_clamped({0.0f, 0.0f}, {1.0f, 1.0f}));
	CHECK(vector2{2.0f, 2.0f}.clamp({0.0f, 0.0f}, {1.0f, 1.0f}) ==
	      vector2{2.0f, 2.0f}.get_clamped({0.0f, 0.0f}, {1.0f, 1.0f}));
	CHECK(vector2{-1.0f, -1.0f}.clamp({0.0f, 0.0f}, {1.0f, 1.0f}) ==
	      vector2{-1.0f, -1.0f}.get_clamped({0.0f, 0.0f}, {1.0f, 1.0f}));
	CHECK(vector2{-1.0f, 1.0f}.clamp({0.0f, 0.0f}, {1.0f, 1.0f}) ==
	      vector2{-1.0f, 1.0f}.get_clamped({0.0f, 0.0f}, {1.0f, 1.0f}));
	CHECK(vector2{-89453.6654f, 3.14159f}.clamp({0.0f, 0.0f}, {1.0f, 1.0f}) ==
	      vector2{-89453.6654f, 3.14159f}.get_clamped({0.0f, 0.0f},
	                                                  {1.0f, 1.0f}));
}

TEST_CASE("vector2::operator[](size_t)", "[working][unittest]")
{
	{
		const vector2 v1{0.0f, 0.0f};

		CHECK(v1.x == 0.0f);
		CHECK(v1.y == 0.0f);
		CHECK(v1.x == v1[0]);
		CHECK(v1.y == v1[1]);

		const vector2 v2{-89453.6654f, 3.14159f};

		CHECK(v2.x == -89453.6654f);
		CHECK(v2.y == 3.14159f);
		CHECK(v2.x == v2[0]);
		CHECK(v2.y == v2[1]);
	}
	{
		vector2 v1{0.0f, 0.0f};

		CHECK(v1.x == 0.0f);
		CHECK(v1.y == 0.0f);
		CHECK(v1.x == v1[0]);
		CHECK(v1.y == v1[1]);

		vector2 v2{-89453.6654f, 3.14159f};

		CHECK(v2.x == -89453.6654f);
		CHECK(v2.y == 3.14159f);
		CHECK(v2.x == v2[0]);
		CHECK(v2.y == v2[1]);
	}
}

TEST_CASE("vector2::operator+(vector2)", "[working][unittest]")
{
	{
		const vector2 v1{0.0f, 0.0f};
		const vector2 v2{-89453.6654f, 3.14159f};
		const vector2 v3 = v1 + v2;

		CHECK(v3.x == -89453.6654f);
		CHECK(v3.y == 3.14159f);
	}
	{
		const vector2 v1{1.0f, 1.0f};
		const vector2 v2{-89453.6654f, 3.14159f};
		const vector2 v3 = v1 + v2;

		CHECK(v3.x == -89453.6654f + 1.0f);
		CHECK(v3.y == 3.14159f + 1.0f);
	}
}

TEST_CASE("vector2::operator-(vector2)", "[working][unittest]")
{
	{
		const vector2 v1{0.0f, 0.0f};
		const vector2 v2{-89453.6654f, 3.14159f};
		const vector2 v3 = v1 - v2;

		CHECK(v3.x == 89453.6654f);
		CHECK(v3.y == -3.14159f);
	}
	{
		const vector2 v1{1.0f, 1.0f};
		const vector2 v2{-89453.6654f, 3.14159f};
		const vector2 v3 = v1 - v2;

		CHECK(v3.x == 1.0f - -89453.6654f);
		CHECK(v3.y == 1.0f - 3.14159f);
	}
}

TEST_CASE("vector2::operator*(float)", "[working][unittest]")
{
	{
		const vector2 v1{0.0f, 0.0f};
		const vector2 v2 = v1 * 2.0f;

		CHECK(v2.x == 0.0f);
		CHECK(v2.y == 0.0f);
	}
	{
		const vector2 v1{0.0f, 0.0f};
		const vector2 v2 = 2.0f * v1;

		CHECK(v2.x == 0.0f);
		CHECK(v2.y == 0.0f);
	}
	{
		const vector2 v1{1.0f, 1.0f};
		const vector2 v2 = v1 * 2.0f;

		CHECK(v2.x == 2.0f);
		CHECK(v2.y == 2.0f);
	}
	{
		const vector2 v1{1.0f, 1.0f};
		const vector2 v2 = 2.0f * v1;

		CHECK(v2.x == 2.0f);
		CHECK(v2.y == 2.0f);
	}
	{
		const vector2 v1{3.0f, -11.0f};
		const vector2 v2 = v1 * 5.0f;

		CHECK(v2.x == v1.x * 5.0f);
		CHECK(v2.y == v1.y * 5.0f);
	}
	{
		const vector2 v1{3.0f, -11.0f};
		const vector2 v2 = 5.0f * v1;

		CHECK(v2.x == v1.x * 5.0f);
		CHECK(v2.y == v1.y * 5.0f);
	}
	{
		const vector2 v1{3.0f, -11.0f};
		const vector2 v2 = v1 * 0.0f;

		CHECK(v2.x == 0.0f);
		CHECK(v2.y == 0.0f);
	}
	{
		const vector2 v1{3.0f, -11.0f};
		const vector2 v2 = 0.0f * v1;

		CHECK(v2.x == 0.0f);
		CHECK(v2.y == 0.0f);
	}
}

TEST_CASE("vector2::operator/(float)", "[working][unittest]")
{
	{
		const vector2 v1{0.0f, 0.0f};
		const vector2 v2 = v1 / 2.0f;

		CHECK(v2.x == 0.0f);
		CHECK(v2.y == 0.0f);
	}
	{
		const vector2 v1{1.0f, 1.0f};
		const vector2 v2 = v1 / 2.0f;

		CHECK(v2.x == 0.5f);
		CHECK(v2.y == 0.5f);
	}
	{
		const vector2 v1{3.0f, -11.0f};
		const vector2 v2 = v1 / 5.0f;

		CHECK(v2.x == v1.x / 5.0f);
		CHECK(v2.y == v1.y / 5.0f);
	}
}

TEST_CASE("vector2::operator+=(vector2)", "[working][unittest]")
{
	{
		vector2 v1{0.0f, 0.0f};
		const vector2 v2{-89453.6654f, 3.14159f};
		v1 += v2;

		CHECK(v1.x == -89453.6654f);
		CHECK(v1.y == 3.14159f);
	}
	{
		vector2 v1{1.0f, 1.0f};
		const vector2 v2{-89453.6654f, 3.14159f};
		v1 += v2;

		CHECK(v1.x == -89453.6654f + 1.0f);
		CHECK(v1.y == 3.14159f + 1.0f);
	}
}

TEST_CASE("vector2::operator-=(vector2)", "[working][unittest]")
{
	{
		vector2 v1{0.0f, 0.0f};
		const vector2 v2{-89453.6654f, 3.14159f};
		v1 -= v2;

		CHECK(v1.x == 89453.6654f);
		CHECK(v1.y == -3.14159f);
	}
	{
		vector2 v1{1.0f, 1.0f};
		const vector2 v2{-89453.6654f, 3.14159f};
		v1 -= v2;

		CHECK(v1.x == 1.0f - -89453.6654f);
		CHECK(v1.y == 1.0f - 3.14159f);
	}
}

TEST_CASE("vector2::operator*=(float)", "[working][unittest]")
{
	{
		vector2 v1{0.0f, 0.0f};
		v1 *= 2.0f;

		CHECK(v1.x == 0.0f);
		CHECK(v1.y == 0.0f);
	}
	{
		vector2 v1{1.0f, 1.0f};
		v1 *= 2.0f;

		CHECK(v1.x == 2.0f);
		CHECK(v1.y == 2.0f);
	}
	{
		vector2 v1{3.0f, -11.0f};
		v1 *= 5.0f;

		CHECK(v1.x == 3.0f * 5.0f);
		CHECK(v1.y == -11.0f * 5.0f);
	}
	{
		vector2 v1{3.0f, -11.0f};
		v1 *= 0.0f;

		CHECK(v1.x ==  0.0f);
		CHECK(v1.y ==  0.0f);
	}
}

TEST_CASE("vector2::operator/=(float)", "[working][unittest]")
{
	{
		vector2 v1{0.0f, 0.0f};
		v1 /= 2.0f;

		CHECK(v1.x == 0.0f);
		CHECK(v1.y == 0.0f);
	}
	{
		vector2 v1{1.0f, 1.0f};
		v1 /= 2.0f;

		CHECK(v1.x == 0.5f);
		CHECK(v1.y == 0.5f);
	}
	{
		vector2 v1{3.0f, -11.0f};
		v1 /= 5.0f;

		CHECK(v1.x == 3.0f / 5.0f);
		CHECK(v1.y == -11.0f / 5.0f);
	}
}

TEST_CASE("vector2::operator-()", "[working][unittest]")
{
	const vector2 v1{0.0f, 0.0f};
	const vector2 v12 = -v1;

	CHECK(v12.x == -0.0f);
	CHECK(v12.y == -0.0f);

	const vector2 v2{-89453.6654f, 3.14159f};
	const vector2 v22 = -v2;

	CHECK(v22.x == 89453.6654f);
	CHECK(v22.y == -3.14159f);
}

TEST_CASE("vector2::operator==()", "[working][unittest]")
{
	const vector2 v1{0.0f, 0.0f};
	const vector2 v2{0.0f, 0.0f};

	CHECK(v1 == v2);
	CHECK(v1.x == v2.x);
	CHECK(v1.y == v2.y);

	const vector2 v3{-89453.6654f, 3.14159f};
	const vector2 v4{-89453.6654f, 3.14159f};
	
	CHECK(v3 == v4);
	CHECK(v3.x == v4.x);
	CHECK(v3.y == v4.y);
	
	CHECK_FALSE(v1 == v3);
}

TEST_CASE("vector2::operator!=()", "[working][unittest]")
{
	const vector2 v1{0.0f, 0.0f};
	const vector2 v2{0.0f, 0.0f};

	CHECK(!(v1 != v2));
	CHECK_FALSE(v1 != v2);
	CHECK(v1.x == v2.x);
	CHECK(v1.y == v2.y);

	const vector2 v3{-89453.6654f, 3.14159f};
	const vector2 v4{-89453.6654f, 3.14159f};
	
	CHECK(!(v3 != v4));
	CHECK_FALSE(v3 != v4);
	CHECK(v3.x == v4.x);
	CHECK(v3.y == v4.y);
	
	CHECK(v1 != v3);
}

TEST_CASE("dot(vector2, vector2)", "[working][unittest]")
{
	CHECK(dot(vector2{ 0.0f, 0.0f }, { 0.0f, 0.0f } ) == 0.0f);
	CHECK(dot(vector2{ 1.0f, 0.0f }, { 0.0f, 1.0f } ) == 0.0f);
	CHECK(dot(vector2{ 1.0f, 0.0f }, { 1.0f, 0.0f } ) == 1.0f);
	CHECK(dot(vector2{ 0.0f, 1.0f }, { 0.0f, 1.0f } ) == 1.0f);
	CHECK(dot(vector2{ 1.0f, 1.0f }, { 1.0f, 0.0f } ) == 1.0f);
	CHECK(dot(vector2{ 0.5f, 1.0f }, { 1.0f, 0.0f } ) == 0.5f);
	CHECK(dot(vector2{ 1.0f, 0.0f }, { 0.5f, 1.0f } ) == 0.5f);
}

TEST_CASE("cross(vector2, vector2)", "[working][unittest]")
{
	CHECK(cross({ 0.0f, 0.0f }, vector2{ 0.0f, 0.0f } ) == 0.0f);
	CHECK(cross({ 1.0f, 0.0f }, vector2{ 0.0f, 1.0f } ) == 1.0f);
	CHECK(cross({ 0.0f, 1.0f }, vector2{ 1.0f, 0.0f } ) == -1.0f);
	CHECK(cross({ 1.0f, 0.0f }, vector2{ 1.0f, 0.0f } ) == 0.0f);
	CHECK(cross({ 0.0f, 1.0f }, vector2{ 0.0f, 1.0f } ) == 0.0f);
	CHECK(cross({ 1.0f, 0.0f }, vector2{ 0.0f, 0.5f }.get_normalized() ) == 1.0f);
	CHECK(cross({ 1.0f, 0.0f }, vector2{ -1.0f, 0.0f } ) == 0.0f);
	CHECK(cross({ 1.0f, 0.0f }, vector2{ 1.0f, 1.0f }.get_normalized() ) == Approx(sqrt2 / 2.0f));
	CHECK(cross({ 2.0f, 0.0f }, vector2{ -1.0f, 0.0f } ) == 0.0f);
	CHECK(cross({ 2.0f, 0.0f }, vector2{ 0.0f, 1.0f } ) == 2.0f);
	CHECK(cross({ 0.0f, 1.0f }, vector2{ 2.0f, 0.0f } ) == -2.0f);
}

TEST_CASE("distance_between(vector2, vector2)", "[working][unittest]")
{
	CHECK(distance_between({ 0.0f, 0.0f }, vector2{ 0.0f, 0.0f }) == 0.0f);
	CHECK(distance_between({ 1.0f, 0.0f }, vector2{ 0.0f, 1.0f }) == Approx(sqrt2));
	CHECK(distance_between({ 0.0f, 1.0f }, vector2{ 1.0f, 0.0f }) == Approx(sqrt2));
	CHECK(distance_between({ 1.0f, 0.0f }, vector2{ 1.0f, 0.0f }) == 0.0f);
	CHECK(distance_between({ 0.0f, 1.0f }, vector2{ 0.0f, 1.0f }) == 0.0f);
	CHECK(distance_between({ 1.0f, 0.0f }, vector2{ 0.0f, 0.5f }.get_normalized()) == Approx(sqrt2));
	CHECK(distance_between({ 1.0f, 0.0f }, vector2{ -1.0f, 0.0f }) == 2.0f);
	CHECK(distance_between({ 2.0f, 0.0f }, vector2{ -1.0f, 0.0f }) == 3.0f);
}

TEST_CASE("distance_squared_between(vector2, vector2)", "[working][unittest]")
{
	CHECK(distance_squared_between({ 0.0f, 0.0f }, vector2{ 0.0f, 0.0f }) == 0.0f);
	CHECK(distance_squared_between({ 1.0f, 0.0f }, vector2{ 0.0f, 1.0f }) == 2.0f);
	CHECK(distance_squared_between({ 0.0f, 1.0f }, vector2{ 1.0f, 0.0f }) == 2.0f);
	CHECK(distance_squared_between({ 1.0f, 0.0f }, vector2{ 1.0f, 0.0f }) == 0.0f);
	CHECK(distance_squared_between({ 0.0f, 1.0f }, vector2{ 0.0f, 1.0f }) == 0.0f);
	CHECK(distance_squared_between({ 1.0f, 0.0f }, vector2{ 0.0f, 0.5f }.get_normalized()) == 2.0f);
	CHECK(distance_squared_between({ 1.0f, 0.0f }, vector2{ -1.0f, 0.0f }) == 4.0f);
	CHECK(distance_squared_between({ 2.0f, 0.0f }, vector2{ -1.0f, 0.0f }) == 9.0f);
}

TEST_CASE("angle_between(vector2, vector2)", "[working][unittest]")
{
	CHECK(angle_between({ 0.0f, 0.0f }, vector2{ 0.0f, 0.0f }) == 0_rad);
	CHECK(angle_between({ 1.0f, 0.0f }, vector2{ 0.0f, 1.0f }) == radian(90_deg));
	CHECK(angle_between({ 1.0f, 0.0f }, vector2{ 1.0f, 1.0f }) == radian(45_deg));
	CHECK(angle_between({ 0.0f, 1.0f }, vector2{ 1.0f, 0.0f }) == radian(-90_deg));
	CHECK(angle_between({ 1.0f, 0.0f }, vector2{ 1.0f, 0.0f }) == radian(0_deg));
	CHECK(angle_between({ 0.0f, 1.0f }, vector2{ 0.0f, 1.0f }) == radian(0_deg));
	CHECK(angle_between({ 1.0f, 0.0f }, vector2{ 0.0f, 0.5f }.get_normalized()) == 1_pi / 2.0f);
	CHECK(angle_between({ 1.0f, 0.0f }, vector2{ -1.0f, 0.0f }) == 1_pi);
	CHECK(angle_between({ 2.0f, 0.0f }, vector2{ -1.0f, 0.0f }) == 1_pi);
}
