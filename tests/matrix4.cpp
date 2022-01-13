#include <catch2/catch.hpp>
#include <cdm_maths.hpp>

#include <array>
#include <iostream>

using namespace cdm;

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
}

TEST_CASE("matrix4::to_array() and matrix4::matrix4(array)", "[working][unittest]")
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
	    //
	    1.0f,         2.5f,     0.1f,        35.9f,           //
	    112.0f,       51.0f,    -660.5f,     232323.23f,      //
	    -89453.6654f, 3.14159f, -3.23123f,   999999999.999f,  //
	    -1.0f,        0.001f,   0.01203212f, 111100.0f        //
	};

	m1 = matrix4(a);
	m2 = matrix4(a);

	REQUIRE(m1.to_array() == m2.to_array());
	REQUIRE(matrix4(m1.to_array()) == matrix4(m2.to_array()));
	REQUIRE(matrix4(m1.to_array()).to_array() == matrix4(m2.to_array()).to_array());

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
	matrix4 m2 = matrix4::translation(t);
	std::array<float, 16> a = m2.to_array();
	std::array<float, 16> g = {
		1.0f, 0.0f, 0.0f, t.x,
		0.0f, 1.0f, 0.0f, t.y,
		0.0f, 0.0f, 1.0f, t.z,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	
	CHECK(a == g);
}

TEST_CASE("matrix4::translation(float, float, float)", "[working][unittest]")
{
	matrix4 m2 = matrix4::translation(0.02f, 0.01f, 0.01f);
	matrix4 m4 = matrix4::translation({0.02f, 0.01f, 0.01f});
	std::array<float, 16> a = m2.to_array();
	std::array<float, 16> g = m4.to_array();
	CHECK(a == g);
}

TEST_CASE("matrix4::scale(vector3)", "[working][unittest]")
{
	vector3 s;
	s = {0.5f, 0.0f, 0.0f};
	matrix4 m2 = matrix4::scale(s);
	std::array<float, 16> a = m2.to_array();
	std::array<float, 16> g = {s.x,  0.0f, 0.0f, 0.0f, 0.0f, s.y,  0.0f, 0.0f,
	                           0.0f, 0.0f, s.z,  0.0f, 0.0f, 0.0f, 0.0f, 1.0f};
	CHECK(s.x == 0.5f);
	CHECK(a == g);
}

TEST_CASE("matrix4::get_determinant()", "[working][unittest]")
{
	matrix4 m2 = matrix4(
	    std::array<float, 16>{1.0f, 4.0f, 2.0f, 3.0f, 8.0f, 10.0f, 2.0f, 6.0f,
	                          1.0f, 7.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f});
	float d = m2.determinant();
	CHECK(d == -100);
}

TEST_CASE("matrix4::get_transposed()", "[working][unittest]")
{
	matrix4 m = matrix4(
	    std::array<float, 16>{1.0f, 4.0f, 2.0f, 3.0f, 8.0f, 10.0f, 2.0f, 6.0f,
	                          1.0f, 7.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f});
	matrix4 t = matrix4(
	    std::array<float, 16>{1.0f, 8.0f, 1.0f, 2.0f, 4.0f, 10.0f, 7.0f, 2.0f,
	                          2.0f, 2.0f, 2.0f, 2.0f, 3.0f, 6.0f, 2.0f, 2.0f});
	matrix4 m2 = m.get_transposed();
	std::array<float, 16> a = t.to_array();
	std::array<float, 16> g = m2.to_array();
	CHECK(a == g);
}

TEST_CASE("matrix4::get_inversed()", "[working][unittest]")
{
	matrix4 m = matrix4(
	    std::array<float, 16>{1.0f, 4.0f, 2.0f, 3.0f, 8.0f, 10.0f, 2.0f, 6.0f,
	                          1.0f, 7.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f});
	matrix4 i = matrix4(std::array<float, 16>{
	    -0.4f, 0.1f, -0.0f, 0.3f, -0.08f, 0.02f, 0.2f, -0.14f, -0.28f, -0.18f,
	    0.2f, 0.76f, 0.76f, 0.06f, -0.4f, -0.42f});
	matrix4 m2 = m.get_inversed();
	std::array<float, 16> a = i.to_array();
	std::array<float, 16> g = m2.to_array();
	REQUIRE(a.size() == g.size());
	for (size_t i = 0; i < a.size(); i++)
		REQUIRE(a[i] == Approx(g[i]));
}

TEST_CASE("matrix4::is_homogenous()", "[working][unittest]")
{
	CHECK(false == matrix4::zero().is_homogenous());
	matrix4 m2 = matrix4(
	    std::array<float, 16>{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
	                          0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f});
	CHECK(true == m2.is_homogenous());
}

TEST_CASE("matrix4::is_orthogonal()", "[working][unittest]")
{
	CHECK(false == matrix4::zero().is_orthogonal());
	matrix4 m2 = matrix4(
	    std::array<float, 16>{1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
	                          0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f});
	CHECK(true == m2.is_orthogonal());
}
