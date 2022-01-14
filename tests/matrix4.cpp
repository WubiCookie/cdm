#include <catch2/catch.hpp>
#include <cdm_maths.hpp>

#include <array>
#include <iostream>

using namespace cdm;

#define REQUIRE_MATRICES_EQUAL(m1, m2)                   \
	REQUIRE(m1.column(0).row(0) == m2.column(0).row(0)); \
	REQUIRE(m1.column(0).row(1) == m2.column(0).row(1)); \
	REQUIRE(m1.column(0).row(2) == m2.column(0).row(2)); \
	REQUIRE(m1.column(0).row(3) == m2.column(0).row(3)); \
	REQUIRE(m1.column(1).row(0) == m2.column(1).row(0)); \
	REQUIRE(m1.column(1).row(1) == m2.column(1).row(1)); \
	REQUIRE(m1.column(1).row(2) == m2.column(1).row(2)); \
	REQUIRE(m1.column(1).row(3) == m2.column(1).row(3)); \
	REQUIRE(m1.column(2).row(0) == m2.column(2).row(0)); \
	REQUIRE(m1.column(2).row(1) == m2.column(2).row(1)); \
	REQUIRE(m1.column(2).row(2) == m2.column(2).row(2)); \
	REQUIRE(m1.column(2).row(3) == m2.column(2).row(3)); \
	REQUIRE(m1.column(3).row(0) == m2.column(3).row(0)); \
	REQUIRE(m1.column(3).row(1) == m2.column(3).row(1)); \
	REQUIRE(m1.column(3).row(2) == m2.column(3).row(2)); \
	REQUIRE(m1.column(3).row(3) == m2.column(3).row(3));

#define REQUIRE_MATRICES_EQUAL_APPROX(m1, m2)                    \
	REQUIRE(m1.column(0).row(0) == Approx(m2.column(0).row(0))); \
	REQUIRE(m1.column(0).row(1) == Approx(m2.column(0).row(1))); \
	REQUIRE(m1.column(0).row(2) == Approx(m2.column(0).row(2))); \
	REQUIRE(m1.column(0).row(3) == Approx(m2.column(0).row(3))); \
	REQUIRE(m1.column(1).row(0) == Approx(m2.column(1).row(0))); \
	REQUIRE(m1.column(1).row(1) == Approx(m2.column(1).row(1))); \
	REQUIRE(m1.column(1).row(2) == Approx(m2.column(1).row(2))); \
	REQUIRE(m1.column(1).row(3) == Approx(m2.column(1).row(3))); \
	REQUIRE(m1.column(2).row(0) == Approx(m2.column(2).row(0))); \
	REQUIRE(m1.column(2).row(1) == Approx(m2.column(2).row(1))); \
	REQUIRE(m1.column(2).row(2) == Approx(m2.column(2).row(2))); \
	REQUIRE(m1.column(2).row(3) == Approx(m2.column(2).row(3))); \
	REQUIRE(m1.column(3).row(0) == Approx(m2.column(3).row(0))); \
	REQUIRE(m1.column(3).row(1) == Approx(m2.column(3).row(1))); \
	REQUIRE(m1.column(3).row(2) == Approx(m2.column(3).row(2))); \
	REQUIRE(m1.column(3).row(3) == Approx(m2.column(3).row(3)));

#define REQUIRE_MATRICES_EQUAL_APPROX_MARGIN(m1, m2, MARGIN) \
	REQUIRE(m1.column(0).row(0) ==                           \
	        Approx(m2.column(0).row(0)).margin(MARGIN));     \
	REQUIRE(m1.column(0).row(1) ==                           \
	        Approx(m2.column(0).row(1)).margin(MARGIN));     \
	REQUIRE(m1.column(0).row(2) ==                           \
	        Approx(m2.column(0).row(2)).margin(MARGIN));     \
	REQUIRE(m1.column(0).row(3) ==                           \
	        Approx(m2.column(0).row(3)).margin(MARGIN));     \
	REQUIRE(m1.column(1).row(0) ==                           \
	        Approx(m2.column(1).row(0)).margin(MARGIN));     \
	REQUIRE(m1.column(1).row(1) ==                           \
	        Approx(m2.column(1).row(1)).margin(MARGIN));     \
	REQUIRE(m1.column(1).row(2) ==                           \
	        Approx(m2.column(1).row(2)).margin(MARGIN));     \
	REQUIRE(m1.column(1).row(3) ==                           \
	        Approx(m2.column(1).row(3)).margin(MARGIN));     \
	REQUIRE(m1.column(2).row(0) ==                           \
	        Approx(m2.column(2).row(0)).margin(MARGIN));     \
	REQUIRE(m1.column(2).row(1) ==                           \
	        Approx(m2.column(2).row(1)).margin(MARGIN));     \
	REQUIRE(m1.column(2).row(2) ==                           \
	        Approx(m2.column(2).row(2)).margin(MARGIN));     \
	REQUIRE(m1.column(2).row(3) ==                           \
	        Approx(m2.column(2).row(3)).margin(MARGIN));     \
	REQUIRE(m1.column(3).row(0) ==                           \
	        Approx(m2.column(3).row(0)).margin(MARGIN));     \
	REQUIRE(m1.column(3).row(1) ==                           \
	        Approx(m2.column(3).row(1)).margin(MARGIN));     \
	REQUIRE(m1.column(3).row(2) ==                           \
	        Approx(m2.column(3).row(2)).margin(MARGIN));     \
	REQUIRE(m1.column(3).row(3) == Approx(m2.column(3).row(3)).margin(MARGIN));

TEST_CASE("matrix4::constructor", "[working][unittest]")
{
	std::array<float, 16> a{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
	                        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

	matrix4 m1 = matrix4(a);

	REQUIRE(m1.column(0).row(0) == 0.0f);
	REQUIRE(m1.column(0).row(1) == 0.0f);
	REQUIRE(m1.column(0).row(2) == 0.0f);
	REQUIRE(m1.column(0).row(3) == 0.0f);
	REQUIRE(m1.column(1).row(0) == 0.0f);
	REQUIRE(m1.column(1).row(1) == 0.0f);
	REQUIRE(m1.column(1).row(2) == 0.0f);
	REQUIRE(m1.column(1).row(3) == 0.0f);
	REQUIRE(m1.column(2).row(0) == 0.0f);
	REQUIRE(m1.column(2).row(1) == 0.0f);
	REQUIRE(m1.column(2).row(2) == 0.0f);
	REQUIRE(m1.column(2).row(3) == 0.0f);
	REQUIRE(m1.column(3).row(0) == 0.0f);
	REQUIRE(m1.column(3).row(1) == 0.0f);
	REQUIRE(m1.column(3).row(2) == 0.0f);
	REQUIRE(m1.column(3).row(3) == 0.0f);

	std::array<float, 16> b{
	    //
	    1.0f,         2.5f,     0.1f,        35.9f,           //
	    112.0f,       51.0f,    -660.5f,     232323.23f,      //
	    -89453.6654f, 3.14159f, -3.23123f,   999999999.999f,  //
	    -1.0f,        0.001f,   0.01203212f, 111100.0f        //
	};

	matrix4 m2 = matrix4(b);

	REQUIRE(m2.column(0).row(0) == 1.0f);
	REQUIRE(m2.column(0).row(1) == 112.0f);
	REQUIRE(m2.column(0).row(2) == -89453.6654f);
	REQUIRE(m2.column(0).row(3) == -1.0f);
	REQUIRE(m2.column(1).row(0) == 2.5f);
	REQUIRE(m2.column(1).row(1) == 51.0f);
	REQUIRE(m2.column(1).row(2) == 3.14159f);
	REQUIRE(m2.column(1).row(3) == 0.001f);
	REQUIRE(m2.column(2).row(0) == 0.1f);
	REQUIRE(m2.column(2).row(1) == -660.5f);
	REQUIRE(m2.column(2).row(2) == -3.23123f);
	REQUIRE(m2.column(2).row(3) == 0.01203212f);
	REQUIRE(m2.column(3).row(0) == 35.9f);
	REQUIRE(m2.column(3).row(1) == 232323.23f);
	REQUIRE(m2.column(3).row(2) == 999999999.999f);
	REQUIRE(m2.column(3).row(3) == 111100.0f);
}

TEST_CASE("matrix4::zero()", "[working][unittest]")
{
	matrix4 m = matrix4::zero();
	REQUIRE(m.column(0).row(0) == 0.0f);
	REQUIRE(m.column(0).row(1) == 0.0f);
	REQUIRE(m.column(0).row(2) == 0.0f);
	REQUIRE(m.column(0).row(3) == 0.0f);
	REQUIRE(m.column(1).row(0) == 0.0f);
	REQUIRE(m.column(1).row(1) == 0.0f);
	REQUIRE(m.column(1).row(2) == 0.0f);
	REQUIRE(m.column(1).row(3) == 0.0f);
	REQUIRE(m.column(2).row(0) == 0.0f);
	REQUIRE(m.column(2).row(1) == 0.0f);
	REQUIRE(m.column(2).row(2) == 0.0f);
	REQUIRE(m.column(2).row(3) == 0.0f);
	REQUIRE(m.column(3).row(0) == 0.0f);
	REQUIRE(m.column(3).row(1) == 0.0f);
	REQUIRE(m.column(3).row(2) == 0.0f);
	REQUIRE(m.column(3).row(3) == 0.0f);
}

TEST_CASE("matrix4::identity()", "[working][unittest]")
{
	matrix4 m = matrix4::identity();

	REQUIRE(m.column(0).row(0) == 1.0f);
	REQUIRE(m.column(1).row(1) == 1.0f);
	REQUIRE(m.column(2).row(2) == 1.0f);
	REQUIRE(m.column(3).row(3) == 1.0f);

	REQUIRE(m.column(0).row(1) == 0.0f);
	REQUIRE(m.column(0).row(2) == 0.0f);
	REQUIRE(m.column(0).row(3) == 0.0f);

	REQUIRE(m.column(1).row(0) == 0.0f);
	REQUIRE(m.column(1).row(2) == 0.0f);
	REQUIRE(m.column(1).row(3) == 0.0f);

	REQUIRE(m.column(2).row(0) == 0.0f);
	REQUIRE(m.column(2).row(1) == 0.0f);
	REQUIRE(m.column(2).row(3) == 0.0f);

	REQUIRE(m.column(3).row(0) == 0.0f);
	REQUIRE(m.column(3).row(1) == 0.0f);
	REQUIRE(m.column(3).row(2) == 0.0f);
}

TEST_CASE("matrix4::operator==(const matrix4&)", "[working][unittest]")
{
	matrix4 m1 = matrix4::zero();
	matrix4 m2 = matrix4::zero();

	REQUIRE(m1 == m2);

	m1 = matrix4::identity();
	m2 = matrix4::identity();

	REQUIRE(m1 == m2);

	std::array<float, 16> a{
	    //
	    1.0f,         2.5f,     0.1f,        35.9f,           //
	    112.0f,       51.0f,    -660.5f,     232323.23f,      //
	    -89453.6654f, 3.14159f, -3.23123f,   999999999.999f,  //
	    -1.0f,        0.001f,   0.01203212f, 111100.0f        //
	};

	m1 = matrix4(a);
	m2 = matrix4(a);

	REQUIRE(m1 == m2);

	m2.column(0).row(0) = 2.0f;

	REQUIRE_FALSE(m1 == m2);

	m1.column(0).row(0) = 1.0f;
	m1.column(0).row(1) = 2.0f;
	m1.column(0).row(2) = 3.0f;
	m1.column(0).row(3) = 4.0f;
	m1.column(1).row(0) = 5.0f;
	m1.column(1).row(1) = 6.0f;
	m1.column(1).row(2) = 7.0f;
	m1.column(1).row(3) = 8.0f;
	m1.column(2).row(0) = 9.0f;
	m1.column(2).row(1) = 10.0f;
	m1.column(2).row(2) = 11.0f;
	m1.column(2).row(3) = 12.0f;
	m1.column(3).row(0) = 13.0f;
	m1.column(3).row(1) = 14.0f;
	m1.column(3).row(2) = 15.0f;
	m1.column(3).row(3) = 16.0f;

	m2.column(0).row(0) = 1.0f;
	m2.column(0).row(1) = 2.0f;
	m2.column(0).row(2) = 3.0f;
	m2.column(0).row(3) = 4.0f;
	m2.column(1).row(0) = 5.0f;
	m2.column(1).row(1) = 6.0f;
	m2.column(1).row(2) = 7.0f;
	m2.column(1).row(3) = 8.0f;
	m2.column(2).row(0) = 9.0f;
	m2.column(2).row(1) = 10.0f;
	m2.column(2).row(2) = 11.0f;
	m2.column(2).row(3) = 12.0f;
	m2.column(3).row(0) = 13.0f;
	m2.column(3).row(1) = 14.0f;
	m2.column(3).row(2) = 15.0f;
	m2.column(3).row(3) = 16.0f;

	REQUIRE(m1 == m2);
	REQUIRE_MATRICES_EQUAL(m1, m2);
}

TEST_CASE("matrix4::to_array() and matrix4::matrix4(array)",
          "[working][unittest]")
{
	matrix4 m1 = matrix4::zero();
	matrix4 m2 = matrix4::zero();

	REQUIRE(m1 == m2);
	REQUIRE(m1.to_array() == m2.to_array());

	m1 = matrix4::identity();
	m2 = matrix4::identity();

	REQUIRE(m1 == m2);
	REQUIRE(m1.to_array() == m2.to_array());

	std::array<float, 16> a{
	    1.0f,         2.5f,     0.1f,        35.9f,           //
	    112.0f,       51.0f,    -660.5f,     232323.23f,      //
	    -89453.6654f, 3.14159f, -3.23123f,   999999999.999f,  //
	    -1.0f,        0.001f,   0.01203212f, 111100.0f        //
	};

	m1 = matrix4(a);
	m2 = matrix4(a);

	REQUIRE(m1.to_array() == m2.to_array());
	REQUIRE(matrix4(m1.to_array()) == matrix4(m2.to_array()));
	REQUIRE(matrix4(m1.to_array()).to_array() ==
	        matrix4(m2.to_array()).to_array());

	m2.column(0).row(0) = 2.0f;

	REQUIRE_FALSE(m1 == m2);
	REQUIRE_FALSE(m1.to_array() == m2.to_array());
}

TEST_CASE("matrix4::rotation(euler_angles)", "[working][unittest]")
{
	euler_angles r{radian(float(pi / 2.0)), 0_rad, 0_rad};

	matrix4 m2 = matrix4::rotation(r);
	matrix4 m4 = matrix4(matrix3::rotation(r));
	std::array<float, 16> a = m2.to_array();
	std::array<float, 16> g = m4.to_array();

	CHECK(a == g);

	matrix3 RX =
	    matrix3(std::array<float, 9>{1.0f, 0.0f, 0.0f, 0.0f, cos(r.x),
	                                 sin(r.x), 0.0f, -sin(r.x), cos(r.x)});
	matrix3 RY =
	    matrix3(std::array<float, 9>{cos(r.y), 0.0f, -sin(r.y), 0.0f, 1.0f,
	                                 0.0f, sin(r.y), 0.0f, cos(r.y)});
	matrix3 RZ =
	    matrix3(std::array<float, 9>{cos(r.z), sin(r.z), 0.0f, -sin(r.z),
	                                 cos(r.z), 0.0f, 0.0f, 0.0f, 1.0f});
	matrix3 m5 = RY * RX * RZ;
	matrix4 m6 = matrix4(m5);
	std::array<float, 16> m = m6.to_array();
	for (size_t i = 0; i < a.size(); i++)
		REQUIRE(a[i] == Approx(m[i]).epsilon(0.1));
}

TEST_CASE("matrix4::rotation(quaternion)", "[working][unittest]")
{
	quaternion q{1.0f, 2.0f, 3.0f, 4.0f};
	q.normalize();

	matrix4 m2 = matrix4::rotation(q);
	matrix3 m3;
	matrix4 m4 = matrix4(m3.rotation(q));
	std::array<float, 16> a = m2.to_array();
	std::array<float, 16> g = m4.to_array();
	CHECK(a == g);
}

TEST_CASE("matrix4::translation(vector3)", "[working][unittest]")
{
	vector3 t{1.0f, 2.0f, 3.0f};
	matrix4 m = matrix4::translation(t);

	REQUIRE(m.column(0).row(0) == 1.0f);
	REQUIRE(m.column(0).row(1) == 0.0f);
	REQUIRE(m.column(0).row(2) == 0.0f);
	REQUIRE(m.column(0).row(3) == 0.0f);

	REQUIRE(m.column(1).row(0) == 0.0f);
	REQUIRE(m.column(1).row(1) == 1.0f);
	REQUIRE(m.column(1).row(2) == 0.0f);
	REQUIRE(m.column(1).row(3) == 0.0f);

	REQUIRE(m.column(2).row(0) == 0.0f);
	REQUIRE(m.column(2).row(1) == 0.0f);
	REQUIRE(m.column(2).row(2) == 1.0f);
	REQUIRE(m.column(2).row(3) == 0.0f);

	REQUIRE(m.column(3).row(0) == t.x);
	REQUIRE(m.column(3).row(1) == t.y);
	REQUIRE(m.column(3).row(2) == t.z);
	REQUIRE(m.column(3).row(3) == 1.0f);
}

TEST_CASE("matrix4::translation(float, float, float)", "[working][unittest]")
{
	float x{0.02f};
	float y{0.01f};
	float z{0.03f};

	matrix4 m = matrix4::translation(x, y, z);

	REQUIRE(m.column(0).row(0) == 1.0f);
	REQUIRE(m.column(0).row(1) == 0.0f);
	REQUIRE(m.column(0).row(2) == 0.0f);
	REQUIRE(m.column(0).row(3) == 0.0f);

	REQUIRE(m.column(1).row(0) == 0.0f);
	REQUIRE(m.column(1).row(1) == 1.0f);
	REQUIRE(m.column(1).row(2) == 0.0f);
	REQUIRE(m.column(1).row(3) == 0.0f);

	REQUIRE(m.column(2).row(0) == 0.0f);
	REQUIRE(m.column(2).row(1) == 0.0f);
	REQUIRE(m.column(2).row(2) == 1.0f);
	REQUIRE(m.column(2).row(3) == 0.0f);

	REQUIRE(m.column(3).row(0) == x);
	REQUIRE(m.column(3).row(1) == y);
	REQUIRE(m.column(3).row(2) == z);
	REQUIRE(m.column(3).row(3) == 1.0f);
}

TEST_CASE("matrix4::scale(vector3)", "[working][unittest]")
{
	vector3 s{0.5f, 0.4f, 0.3f};
	matrix4 m = matrix4::scale(s);

	REQUIRE(m.column(0).row(0) == s.x);
	REQUIRE(m.column(1).row(1) == s.y);
	REQUIRE(m.column(2).row(2) == s.z);
	REQUIRE(m.column(3).row(3) == 1.0f);

	REQUIRE(m.column(0).row(1) == 0.0f);
	REQUIRE(m.column(0).row(2) == 0.0f);
	REQUIRE(m.column(0).row(3) == 0.0f);

	REQUIRE(m.column(1).row(0) == 0.0f);
	REQUIRE(m.column(1).row(2) == 0.0f);
	REQUIRE(m.column(1).row(3) == 0.0f);

	REQUIRE(m.column(2).row(0) == 0.0f);
	REQUIRE(m.column(2).row(1) == 0.0f);
	REQUIRE(m.column(2).row(3) == 0.0f);

	REQUIRE(m.column(3).row(0) == 0.0f);
	REQUIRE(m.column(3).row(1) == 0.0f);
	REQUIRE(m.column(3).row(2) == 0.0f);
}

TEST_CASE("matrix4::scale(float, float, float)", "[working][unittest]")
{
	float x{0.5f};
	float y{0.4f};
	float z{0.3f};
	matrix4 m = matrix4::scale(x, y, z);

	REQUIRE(m.column(0).row(0) == x);
	REQUIRE(m.column(1).row(1) == y);
	REQUIRE(m.column(2).row(2) == z);
	REQUIRE(m.column(3).row(3) == 1.0f);

	REQUIRE(m.column(0).row(1) == 0.0f);
	REQUIRE(m.column(0).row(2) == 0.0f);
	REQUIRE(m.column(0).row(3) == 0.0f);

	REQUIRE(m.column(1).row(0) == 0.0f);
	REQUIRE(m.column(1).row(2) == 0.0f);
	REQUIRE(m.column(1).row(3) == 0.0f);

	REQUIRE(m.column(2).row(0) == 0.0f);
	REQUIRE(m.column(2).row(1) == 0.0f);
	REQUIRE(m.column(2).row(3) == 0.0f);

	REQUIRE(m.column(3).row(0) == 0.0f);
	REQUIRE(m.column(3).row(1) == 0.0f);
	REQUIRE(m.column(3).row(2) == 0.0f);
}

TEST_CASE("matrix4::scale(float)", "[working][unittest]")
{
	float s{0.5f};
	matrix4 m = matrix4::scale(s);

	REQUIRE(m.column(0).row(0) == s);
	REQUIRE(m.column(1).row(1) == s);
	REQUIRE(m.column(2).row(2) == s);
	REQUIRE(m.column(3).row(3) == 1.0f);

	REQUIRE(m.column(0).row(1) == 0.0f);
	REQUIRE(m.column(0).row(2) == 0.0f);
	REQUIRE(m.column(0).row(3) == 0.0f);

	REQUIRE(m.column(1).row(0) == 0.0f);
	REQUIRE(m.column(1).row(2) == 0.0f);
	REQUIRE(m.column(1).row(3) == 0.0f);

	REQUIRE(m.column(2).row(0) == 0.0f);
	REQUIRE(m.column(2).row(1) == 0.0f);
	REQUIRE(m.column(2).row(3) == 0.0f);

	REQUIRE(m.column(3).row(0) == 0.0f);
	REQUIRE(m.column(3).row(1) == 0.0f);
	REQUIRE(m.column(3).row(2) == 0.0f);
}

TEST_CASE("matrix4::rotation_around_x(radian)", "[working][unittest]")
{
	const vector4 v0{0.0f, 0.0f, 1.0f, 0.0f};

	{
		radian rotation{0_rad};
		matrix4 m = matrix4::rotation_around_x(rotation);
		vector4 v1 = m * v0;
		CHECK(v1 == v1);
	}

	const auto almost_0 = Approx(0.0f).margin(1.0e-7);

	{
		radian rotation{90_deg};
		//matrix4 m = matrix4::rotation_around_x(rotation);
		matrix4 m = matrix4::rotation(quaternion{ vector3{1.0f, 0.0f, 0.0f}, static_pi_fraction<1,2>{} });
		vector4 v1 = m * v0;
		CHECK(std::abs(v1.x) == almost_0);
		CHECK(v1.y == Approx(1.0f).margin(1.0e-7));
		CHECK(std::abs(v1.z) == almost_0);
		CHECK(std::abs(v1.w) == almost_0);
	}

	{
		radian rotation{pi};
		matrix4 m = matrix4::rotation_around_x(rotation);
		vector4 v1 = m * v0;
		CHECK(v1.x == Approx(-v0.x).margin(1.0e-7));
		CHECK(v1.y == Approx(-v0.y).margin(1.0e-7));
		CHECK(v1.z == Approx(-v0.z).margin(1.0e-7));
		CHECK(v1.w == Approx(-v0.w).margin(1.0e-7));
	}
}

TEST_CASE("matrix4::rotation_around_y(radian)", "[working][unittest]")
{
	const vector4 v0{1.0f, 0.0f, 0.0f, 0.0f};

	{
		radian rotation{0_rad};
		matrix4 m = matrix4::rotation_around_y(rotation);
		vector4 v1 = m * v0;
		CHECK(v1 == v1);
	}

	const auto almost_0 = Approx(0.0f).margin(1.0e-7);

	{
		radian rotation{90_deg};
		matrix4 m = matrix4::rotation_around_y(rotation);
		vector4 v1 = m * v0;
		CHECK(std::abs(v1.x) == almost_0);
		CHECK(std::abs(v1.y) == almost_0);
		CHECK(v1.z == Approx(1.0f).margin(1.0e-7));
		CHECK(std::abs(v1.w) == almost_0);
	}

	{
		radian rotation{pi};
		matrix4 m = matrix4::rotation_around_y(rotation);
		vector4 v1 = m * v0;
		CHECK(v1.x == Approx(-v0.x).margin(1.0e-7));
		CHECK(v1.y == Approx(-v0.y).margin(1.0e-7));
		CHECK(v1.z == Approx(-v0.z).margin(1.0e-7));
		CHECK(v1.w == Approx(-v0.w).margin(1.0e-7));
	}
}

TEST_CASE("matrix4::rotation_around_z(radian)", "[working][unittest]")
{
	const vector4 v0{0.0f, 1.0f, 0.0f, 0.0f};

	{
		radian rotation{0_rad};
		matrix4 m = matrix4::rotation_around_z(rotation);
		vector4 v1 = m * v0;
		CHECK(v1 == v1);
	}

	const auto almost_0 = Approx(0.0f).margin(1.0e-7);

	{
		radian rotation{90_deg};
		matrix4 m = matrix4::rotation_around_z(rotation);
		vector4 v1 = m * v0;
		CHECK(v1.x == Approx(1.0f).margin(1.0e-7));
		CHECK(std::abs(v1.y) == almost_0);
		CHECK(std::abs(v1.z) == almost_0);
		CHECK(std::abs(v1.w) == almost_0);
	}
	
	{
		radian rotation{pi};
		matrix4 m = matrix4::rotation_around_z(rotation);
		vector4 v1 = m * v0;
		CHECK(v1.x == Approx(-v0.x).margin(1.0e-7));
		CHECK(v1.y == Approx(-v0.y).margin(1.0e-7));
		CHECK(v1.z == Approx(-v0.z).margin(1.0e-7));
		CHECK(v1.w == Approx(-v0.w).margin(1.0e-7));
	}
}

TEST_CASE("matrix4::rotation_around_z(complex)", "[working][unittest]")
{
	const vector4 v0{0.0f, 1.0f, 0.0f, 0.0f};

	{
		complex rotation{0_rad};
		matrix4 m = matrix4::rotation_around_z(rotation);
		vector4 v1 = m * v0;
		CHECK(v1 == v1);
	}

	const auto almost_0 = Approx(0.0f).margin(1.0e-7);

	{
		complex rotation{90_deg};
		//matrix4 m = matrix4::rotation_around_z(rotation);
		matrix4 m = matrix4::rotation(quaternion{ vector3{0.0f, 0.0f, 1.0f}, static_pi_fraction<1,2>{} });
		vector4 v1 = m * v0;
		CHECK(v1.x == Approx(1.0f).margin(1.0e-7));
		CHECK(std::abs(v1.y) == almost_0);
		CHECK(std::abs(v1.z) == almost_0);
		CHECK(std::abs(v1.w) == almost_0);
	}
	
	{
		complex rotation{1_pi};
		matrix4 m = matrix4::rotation_around_z(rotation);
		vector4 v1 = m * v0;
		CHECK(v1.x == Approx(-v0.x).margin(1.0e-7));
		CHECK(v1.y == Approx(-v0.y).margin(1.0e-7));
		CHECK(v1.z == Approx(-v0.z).margin(1.0e-7));
		CHECK(v1.w == Approx(-v0.w).margin(1.0e-7));
	}
}

TEST_CASE("matrix4::rows(vector4, vector4, vector4, vector4)",
          "[working][unittest]")
{
	matrix4 m = matrix4::rows(
	    {1.0f, 2.0f, 3.0f, 4.0f}, {5.0f, 6.0f, 7.0f, 8.0f},
	    {9.0f, 10.0f, 11.0f, 12.0f}, {13.0f, 14.0f, 15.0f, 16.0f});
	REQUIRE(m.row(0).column(0) == 1.0f);
	REQUIRE(m.row(0).column(1) == 2.0f);
	REQUIRE(m.row(0).column(2) == 3.0f);
	REQUIRE(m.row(0).column(3) == 4.0f);
	REQUIRE(m.row(1).column(0) == 5.0f);
	REQUIRE(m.row(1).column(1) == 6.0f);
	REQUIRE(m.row(1).column(2) == 7.0f);
	REQUIRE(m.row(1).column(3) == 8.0f);
	REQUIRE(m.row(2).column(0) == 9.0f);
	REQUIRE(m.row(2).column(1) == 10.0f);
	REQUIRE(m.row(2).column(2) == 11.0f);
	REQUIRE(m.row(2).column(3) == 12.0f);
	REQUIRE(m.row(3).column(0) == 13.0f);
	REQUIRE(m.row(3).column(1) == 14.0f);
	REQUIRE(m.row(3).column(2) == 15.0f);
	REQUIRE(m.row(3).column(3) == 16.0f);
}

TEST_CASE("matrix4::row(int)", "[working][unittest]")
{
	vector4 v0{1.0f, 2.0f, 3.0f, 4.0f};
	vector4 v1{5.0f, 6.0f, 7.0f, 8.0f};
	vector4 v2{9.0f, 10.0f, 11.0f, 12.0f};
	vector4 v3{13.0f, 14.0f, 15.0f, 16.0f};

	matrix4 m = matrix4::rows(v0, v1, v2, v3);

	REQUIRE(v0 == m.row(0));
	REQUIRE(v1 == m.row(1));
	REQUIRE(v2 == m.row(2));
	REQUIRE(v3 == m.row(3));
}

TEST_CASE("matrix4::column(int)", "[working][unittest]")
{
	vector4 v0{1.0f, 2.0f, 3.0f, 4.0f};
	vector4 v1{5.0f, 6.0f, 7.0f, 8.0f};
	vector4 v2{9.0f, 10.0f, 11.0f, 12.0f};
	vector4 v3{13.0f, 14.0f, 15.0f, 16.0f};

	matrix4 m = matrix4::columns(v0, v1, v2, v3);

	REQUIRE(v0 == m.column(0));
	REQUIRE(v1 == m.column(1));
	REQUIRE(v2 == m.column(2));
	REQUIRE(v3 == m.column(3));
}

TEST_CASE("matrix4::columns(vector4, vector4, vector4, vector4)",
          "[working][unittest]")
{
	matrix4 m = matrix4::columns(
	    {1.0f, 2.0f, 3.0f, 4.0f}, {5.0f, 6.0f, 7.0f, 8.0f},
	    {9.0f, 10.0f, 11.0f, 12.0f}, {13.0f, 14.0f, 15.0f, 16.0f});
	REQUIRE(m.column(0).row(0) == 1.0f);
	REQUIRE(m.column(0).row(1) == 2.0f);
	REQUIRE(m.column(0).row(2) == 3.0f);
	REQUIRE(m.column(0).row(3) == 4.0f);
	REQUIRE(m.column(1).row(0) == 5.0f);
	REQUIRE(m.column(1).row(1) == 6.0f);
	REQUIRE(m.column(1).row(2) == 7.0f);
	REQUIRE(m.column(1).row(3) == 8.0f);
	REQUIRE(m.column(2).row(0) == 9.0f);
	REQUIRE(m.column(2).row(1) == 10.0f);
	REQUIRE(m.column(2).row(2) == 11.0f);
	REQUIRE(m.column(2).row(3) == 12.0f);
	REQUIRE(m.column(3).row(0) == 13.0f);
	REQUIRE(m.column(3).row(1) == 14.0f);
	REQUIRE(m.column(3).row(2) == 15.0f);
	REQUIRE(m.column(3).row(3) == 16.0f);
}

TEST_CASE("matrix4::is_orthogonal()", "[working][unittest]")
{
	CHECK(false == matrix4::zero().is_orthogonal());
	matrix4 m2 = matrix4(
	    std::array<float, 16>{1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
	                          0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f});
	CHECK(true == m2.is_orthogonal());
}

TEST_CASE("matrix4::is_homogenous()", "[working][unittest]")
{
	CHECK(false == matrix4::zero().is_homogenous());
	matrix4 m2 = matrix4(
	    std::array<float, 16>{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
	                          0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f});
	CHECK(true == m2.is_homogenous());
}

TEST_CASE("matrix4::inverse()", "[working][unittest]")
{
	matrix4 m1;
	matrix4 m2;

	// translation
	{
		vector3 translation{1.0f, 513.0f, -984.5f};

		m1 = matrix4::translation(translation);
		m1.inverse();
		m2 = matrix4::translation(-translation);

		CHECK(m1 == m2);
	}

	// scale
	{
		vector3 scale{4589.0f, 132.015f, 0.00125f};
		vector3 invScale{1.0f / scale.x, 1.0f / scale.y, 1.0f / scale.z};

		m1 = matrix4::scale(scale);
		m1.inverse();
		m2 = matrix4::scale(invScale);

		REQUIRE_MATRICES_EQUAL_APPROX(m1, m2);
	}

	const radian rotation = 90_deg;

	// rotation around x
	{
		m1 = matrix4::rotation_around_x(rotation);
		m1.inverse();
		m2 = m1;

		const vector4 p0{1.f, 2.0f, 3.0f, 1.0f};
		vector4 p1 = m1 * p0;
		vector4 p2 = m1 * p1;

		CHECK(p0.x == Approx(p2.x));
		CHECK(p0.y == Approx(p2.y));
		CHECK(p0.z == Approx(p2.z));
	}

	// rotation around y
	{
		m1 = matrix4::rotation_around_y(rotation);
		m1.inverse();
		m2 = m1;

		const vector4 p0{1.f, 2.0f, 3.0f, 1.0f};
		vector4 p1 = m1 * p0;
		vector4 p2 = m1 * p1;

		CHECK(p0.x == Approx(p2.x));
		CHECK(p0.y == Approx(p2.y));
		CHECK(p0.z == Approx(p2.z));
	}

	// rotation around z
	{
		m1 = matrix4::rotation_around_z(rotation);
		m1.inverse();
		m2 = m1;

		const vector4 p0{1.f, 2.0f, 3.0f, 1.0f};
		vector4 p1 = m1 * p0;
		vector4 p2 = m1 * p1;

		CHECK(p0.x == Approx(p2.x));
		CHECK(p0.y == Approx(p2.y));
		CHECK(p0.z == Approx(p2.z));
	}
}

TEST_CASE("matrix4::get_inversed()", "[working][unittest]")
{
	std::array<float, 16> a{1.0f, 4.0f, 2.0f, 3.0f, 8.0f, 10.0f, 2.0f, 6.0f,
	                        1.0f, 7.0f, 2.0f, 2.0f, 2.0f, 2.0f,  2.0f, 2.0f};
	std::array<float, 16> b{-0.4f, 0.1f,   -0.0f,  0.3f,   -0.08f, 0.02f,
	                        0.2f,  -0.14f, -0.28f, -0.18f, 0.2f,   0.76f,
	                        0.76f, 0.06f,  -0.4f,  -0.42f};

	matrix4 m1{a};
	matrix4 m2{m1.get_inversed()};
	matrix4 m3{b};

	REQUIRE_MATRICES_EQUAL_APPROX(m2, m3);

	m1.inverse();
	REQUIRE(m1 == m2);
}

TEST_CASE("matrix4::transpose()", "[working][unittest]")
{
	matrix4 m1;

	m1.column(0).row(0) = 1.0f;
	m1.column(0).row(1) = 2.0f;
	m1.column(0).row(2) = 3.0f;
	m1.column(0).row(3) = 4.0f;
	m1.column(1).row(0) = 5.0f;
	m1.column(1).row(1) = 6.0f;
	m1.column(1).row(2) = 7.0f;
	m1.column(1).row(3) = 8.0f;
	m1.column(2).row(0) = 9.0f;
	m1.column(2).row(1) = 10.0f;
	m1.column(2).row(2) = 11.0f;
	m1.column(2).row(3) = 12.0f;
	m1.column(3).row(0) = 13.0f;
	m1.column(3).row(1) = 14.0f;
	m1.column(3).row(2) = 15.0f;
	m1.column(3).row(3) = 16.0f;

	matrix4 m2 = m1;
	m2.transpose();

	REQUIRE(m2.row(0).column(0) == 1.0f);
	REQUIRE(m2.row(0).column(1) == 2.0f);
	REQUIRE(m2.row(0).column(2) == 3.0f);
	REQUIRE(m2.row(0).column(3) == 4.0f);
	REQUIRE(m2.row(1).column(0) == 5.0f);
	REQUIRE(m2.row(1).column(1) == 6.0f);
	REQUIRE(m2.row(1).column(2) == 7.0f);
	REQUIRE(m2.row(1).column(3) == 8.0f);
	REQUIRE(m2.row(2).column(0) == 9.0f);
	REQUIRE(m2.row(2).column(1) == 10.0f);
	REQUIRE(m2.row(2).column(2) == 11.0f);
	REQUIRE(m2.row(2).column(3) == 12.0f);
	REQUIRE(m2.row(3).column(0) == 13.0f);
	REQUIRE(m2.row(3).column(1) == 14.0f);
	REQUIRE(m2.row(3).column(2) == 15.0f);
	REQUIRE(m2.row(3).column(3) == 16.0f);
}

TEST_CASE("matrix4::get_transposed()", "[working][unittest]")
{
	matrix4 m1;

	m1.column(0).row(0) = 1.0f;
	m1.column(0).row(1) = 2.0f;
	m1.column(0).row(2) = 3.0f;
	m1.column(0).row(3) = 4.0f;
	m1.column(1).row(0) = 5.0f;
	m1.column(1).row(1) = 6.0f;
	m1.column(1).row(2) = 7.0f;
	m1.column(1).row(3) = 8.0f;
	m1.column(2).row(0) = 9.0f;
	m1.column(2).row(1) = 10.0f;
	m1.column(2).row(2) = 11.0f;
	m1.column(2).row(3) = 12.0f;
	m1.column(3).row(0) = 13.0f;
	m1.column(3).row(1) = 14.0f;
	m1.column(3).row(2) = 15.0f;
	m1.column(3).row(3) = 16.0f;

	matrix4 m2 = m1;
	m2.transpose();

	CHECK(m1.get_transposed() == m2);
}

TEST_CASE("matrix4::determinant()", "[working][unittest]")
{
	matrix4 m2 = matrix4(
	    std::array<float, 16>{1.0f, 4.0f, 2.0f, 3.0f, 8.0f, 10.0f, 2.0f, 6.0f,
	                          1.0f, 7.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f});
	float d = m2.determinant();
	CHECK(d == -100);
}

TEST_CASE("matrix4::operator*(const matrix4&)", "[working][unittest]")
{
	matrix4 m1;

	m1.column(0).row(0) = 1.0f;
	m1.column(0).row(1) = 2.0f;
	m1.column(0).row(2) = 3.0f;
	m1.column(0).row(3) = 4.0f;
	m1.column(1).row(0) = 5.0f;
	m1.column(1).row(1) = 6.0f;
	m1.column(1).row(2) = 7.0f;
	m1.column(1).row(3) = 8.0f;
	m1.column(2).row(0) = 9.0f;
	m1.column(2).row(1) = 10.0f;
	m1.column(2).row(2) = 11.0f;
	m1.column(2).row(3) = 12.0f;
	m1.column(3).row(0) = 13.0f;
	m1.column(3).row(1) = 14.0f;
	m1.column(3).row(2) = 15.0f;
	m1.column(3).row(3) = 16.0f;

	matrix4 m2 = m1 * m1;

	REQUIRE(m2.column(0).row(0) == 90.0f);
	REQUIRE(m2.column(0).row(1) == 100.0f);
	REQUIRE(m2.column(0).row(2) == 110.0f);
	REQUIRE(m2.column(0).row(3) == 120.0f);
	REQUIRE(m2.column(1).row(0) == 202.0f);
	REQUIRE(m2.column(1).row(1) == 228.0f);
	REQUIRE(m2.column(1).row(2) == 254.0f);
	REQUIRE(m2.column(1).row(3) == 280.0f);
	REQUIRE(m2.column(2).row(0) == 314.0f);
	REQUIRE(m2.column(2).row(1) == 356.0f);
	REQUIRE(m2.column(2).row(2) == 398.0f);
	REQUIRE(m2.column(2).row(3) == 440.0f);
	REQUIRE(m2.column(3).row(0) == 426.0f);
	REQUIRE(m2.column(3).row(1) == 484.0f);
	REQUIRE(m2.column(3).row(2) == 542.0f);
	REQUIRE(m2.column(3).row(3) == 600.0f);

	vector3 t{123.0f, 456.0f, 789.0f};
	m1 = matrix4::translation(t);
	m2 = m1 * m1;

	REQUIRE(m2.column(0).row(0) == 1.0f);
	REQUIRE(m2.column(0).row(1) == 0.0f);
	REQUIRE(m2.column(0).row(2) == 0.0f);
	REQUIRE(m2.column(0).row(3) == 0.0f);
	REQUIRE(m2.column(1).row(0) == 0.0f);
	REQUIRE(m2.column(1).row(1) == 1.0f);
	REQUIRE(m2.column(1).row(2) == 0.0f);
	REQUIRE(m2.column(1).row(3) == 0.0f);
	REQUIRE(m2.column(2).row(0) == 0.0f);
	REQUIRE(m2.column(2).row(1) == 0.0f);
	REQUIRE(m2.column(2).row(2) == 1.0f);
	REQUIRE(m2.column(2).row(3) == 0.0f);
	REQUIRE(m2.column(3).row(0) == (t.x * 2.0f));
	REQUIRE(m2.column(3).row(1) == (t.y * 2.0f));
	REQUIRE(m2.column(3).row(2) == (t.z * 2.0f));
	REQUIRE(m2.column(3).row(3) == 1.0f);
}

TEST_CASE("matrix4::operator*(vector4)", "[working][unittest]")
{
	const vector4 v0{1.0f, 2.0f, 3.0f, 1.0f};

	const vector3 t{123.0f, 456.0f, 789.0f};
	matrix4 m = matrix4::translation(t);

	vector4 v1 = m * v0;

	REQUIRE(v1.x == (v0.x + t.x));
	REQUIRE(v1.y == (v0.y + t.y));
	REQUIRE(v1.z == (v0.z + t.z));

	float scale{20.01f};
	m = matrix4::scale(scale);

	v1 = m * v0;

	REQUIRE(v1.x == (v0.x * scale));
	REQUIRE(v1.y == (v0.y * scale));
	REQUIRE(v1.z == (v0.z * scale));
}
