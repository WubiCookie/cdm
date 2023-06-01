#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <common.hpp>

TEST_CASE("matrix3::matrix3(std::array<float, 9>)",
          "[working][unittest][matrix3]")
{
	const std::array<float, 9> a{0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
	                             0.0f, 0.0f, 0.0f, 0.0f};

	matrix3 m1 = matrix3(a);

	REQUIRE(m1.column(0).row(0) == 0.0f);
	REQUIRE(m1.column(0).row(1) == 0.0f);
	REQUIRE(m1.column(0).row(2) == 0.0f);
	REQUIRE(m1.column(1).row(0) == 0.0f);
	REQUIRE(m1.column(1).row(1) == 0.0f);
	REQUIRE(m1.column(1).row(2) == 0.0f);
	REQUIRE(m1.column(2).row(0) == 0.0f);
	REQUIRE(m1.column(2).row(1) == 0.0f);
	REQUIRE(m1.column(2).row(2) == 0.0f);

	const std::array<float, 9> b{
	    //
	    1.0f,         2.5f,     0.1f,      //
	    112.0f,       51.0f,    -660.5f,   //
	    -89453.6654f, 3.14159f, -3.23123f  //
	};

	matrix3 m2 = matrix3(b);

	REQUIRE(m2.column(0).row(0) == 1.0f);
	REQUIRE(m2.column(0).row(1) == 112.0f);
	REQUIRE(m2.column(0).row(2) == -89453.6654f);
	REQUIRE(m2.column(1).row(0) == 2.5f);
	REQUIRE(m2.column(1).row(1) == 51.0f);
	REQUIRE(m2.column(1).row(2) == 3.14159f);
	REQUIRE(m2.column(2).row(0) == 0.1f);
	REQUIRE(m2.column(2).row(1) == -660.5f);
	REQUIRE(m2.column(2).row(2) == -3.23123f);
}

TEST_CASE("matrix3::matrix3(matrix3)", "[working][unittest][matrix3]")
{
	const std::array<float, 9> a{1.0f, 2.0f, 3.0f, 4.0f, 5.0f,
	                             6.0f, 7.0f, 8.0f, 9.0f};

	const matrix3 m1(a);

	REQUIRE(m1.column(0).row(0) == 1.0f);
	REQUIRE(m1.column(1).row(0) == 2.0f);
	REQUIRE(m1.column(2).row(0) == 3.0f);
	REQUIRE(m1.column(0).row(1) == 4.0f);
	REQUIRE(m1.column(1).row(1) == 5.0f);
	REQUIRE(m1.column(2).row(1) == 6.0f);
	REQUIRE(m1.column(0).row(2) == 7.0f);
	REQUIRE(m1.column(1).row(2) == 8.0f);
	REQUIRE(m1.column(2).row(2) == 9.0f);

	const matrix3 m2(m1);

	CHECK(m2.column(0).row(0) == m1.column(0).row(0));
	CHECK(m2.column(1).row(0) == m1.column(1).row(0));
	CHECK(m2.column(2).row(0) == m1.column(2).row(0));
	CHECK(m2.column(0).row(1) == m1.column(0).row(1));
	CHECK(m2.column(1).row(1) == m1.column(1).row(1));
	CHECK(m2.column(2).row(1) == m1.column(2).row(1));
	CHECK(m2.column(0).row(2) == m1.column(0).row(2));
	CHECK(m2.column(1).row(2) == m1.column(1).row(2));
	CHECK(m2.column(2).row(2) == m1.column(2).row(2));
}

TEST_CASE("matrix3::zero()", "[working][unittest][matrix3]")
{
	matrix3 m = matrix3::zero();
	REQUIRE(m.column(0).row(0) == 0.0f);
	REQUIRE(m.column(0).row(1) == 0.0f);
	REQUIRE(m.column(0).row(2) == 0.0f);
	REQUIRE(m.column(1).row(0) == 0.0f);
	REQUIRE(m.column(1).row(1) == 0.0f);
	REQUIRE(m.column(1).row(2) == 0.0f);
	REQUIRE(m.column(2).row(0) == 0.0f);
	REQUIRE(m.column(2).row(1) == 0.0f);
	REQUIRE(m.column(2).row(2) == 0.0f);
}

TEST_CASE("matrix3::identity()", "[working][unittest][matrix3]")
{
	matrix3 m = matrix3::identity();

	REQUIRE(m.column(0).row(0) == 1.0f);
	REQUIRE(m.column(1).row(1) == 1.0f);
	REQUIRE(m.column(2).row(2) == 1.0f);

	REQUIRE(m.column(0).row(1) == 0.0f);
	REQUIRE(m.column(0).row(2) == 0.0f);

	REQUIRE(m.column(1).row(0) == 0.0f);
	REQUIRE(m.column(1).row(2) == 0.0f);

	REQUIRE(m.column(2).row(0) == 0.0f);
	REQUIRE(m.column(2).row(1) == 0.0f);
}

TEST_CASE("matrix3::operator==(const matrix3&)",
          "[working][unittest][matrix3]")
{
	matrix3 m1 = matrix3::zero();
	matrix3 m2 = matrix3::zero();

	REQUIRE(m1 == m2);

	m1 = matrix3::identity();
	m2 = matrix3::identity();

	REQUIRE(m1 == m2);

	const std::array<float, 9> a{
	    1.0f,         2.5f,     0.1f,      //
	    112.0f,       51.0f,    -660.5f,   //
	    -89453.6654f, 3.14159f, -3.23123f  //
	};

	m1 = matrix3(a);
	m2 = matrix3(a);

	REQUIRE(m1 == m2);

	m2.column(0).row(0) = 2.0f;

	REQUIRE_FALSE(m1 == m2);

	m1.column(0).row(0) = 1.0f;
	m1.column(0).row(1) = 2.0f;
	m1.column(0).row(2) = 3.0f;
	m1.column(1).row(0) = 5.0f;
	m1.column(1).row(1) = 6.0f;
	m1.column(1).row(2) = 7.0f;
	m1.column(2).row(0) = 9.0f;
	m1.column(2).row(1) = 10.0f;
	m1.column(2).row(2) = 11.0f;

	m2.column(0).row(0) = 1.0f;
	m2.column(0).row(1) = 2.0f;
	m2.column(0).row(2) = 3.0f;
	m2.column(1).row(0) = 5.0f;
	m2.column(1).row(1) = 6.0f;
	m2.column(1).row(2) = 7.0f;
	m2.column(2).row(0) = 9.0f;
	m2.column(2).row(1) = 10.0f;
	m2.column(2).row(2) = 11.0f;

	REQUIRE(m1 == m2);
	REQUIRE_THAT(m1, Matrix3Matcher(m2, 0.0));
}

TEST_CASE("matrix3::to_array() and matrix3::matrix3(array)",
          "[working][unittest][matrix3]")
{
	matrix3 m1 = matrix3::zero();
	matrix3 m2 = matrix3::zero();

	REQUIRE(m1 == m2);
	REQUIRE(m1.to_array() == m2.to_array());

	m1 = matrix3::identity();
	m2 = matrix3::identity();

	REQUIRE(m1 == m2);
	REQUIRE(m1.to_array() == m2.to_array());

	const std::array<float, 9> a{
	    1.0f,         2.5f,     0.1f,      //
	    112.0f,       51.0f,    -660.5f,   //
	    -89453.6654f, 3.14159f, -3.23123f  //
	};

	m1 = matrix3(a);
	m2 = matrix3(a);

	REQUIRE(m1.to_array() == m2.to_array());
	REQUIRE(matrix3(m1.to_array()) == matrix3(m2.to_array()));
	REQUIRE(matrix3(m1.to_array()).to_array() ==
	        matrix3(m2.to_array()).to_array());

	m2.column(0).row(0) = 2.0f;

	REQUIRE_FALSE(m1 == m2);
	REQUIRE_FALSE(m1.to_array() == m2.to_array());
}

TEST_CASE("matrix3::rotation(euler_angles)_X", "[working][unittest][matrix3]")
{
	const auto angle = radian(float(pi / 2.0));

	const euler_angles r{angle, 0_rad, 0_rad};

	const matrix3 m = matrix3::rotation(r);
	const std::array<float, 9> a = m.to_array();

	const auto rotx = m * vector3(1.0f, 0.0f, 0.0f);
	CHECK_THAT(rotx, Vector3Matcher({1.0f, 0.0f, 0.0f}));
	const auto roty = m * vector3(0.0f, 1.0f, 0.0f);
	CHECK_THAT(roty, Vector3Matcher({0.0f, 0.0f, 1.0f}));
	const auto rotz = m * vector3(0.0f, 0.0f, 1.0f);
	CHECK_THAT(rotz, Vector3Matcher({0.0f, -1.0f, 0.0f}));

	{
		const matrix3 m2 = matrix3::rotation_around_x(angle);

		CHECK_THAT(m, Matrix3Matcher(m2));
	}
}

TEST_CASE("matrix3::rotation(euler_angles)_Y", "[working][unittest][matrix3]")
{
	const auto angle = radian(float(pi / 2.0));

	const euler_angles r{0_rad, angle, 0_rad};

	const matrix3 m = matrix3::rotation(r);
	const std::array<float, 9> a = m.to_array();

	const auto rotx = m * vector3(1.0f, 0.0f, 0.0f);
	CHECK_THAT(rotx, Vector3Matcher({0.0f, 0.0f, -1.0f}));
	const auto roty = m * vector3(0.0f, 1.0f, 0.0f);
	CHECK_THAT(roty, Vector3Matcher({0.0f, 1.0f, 0.0f}));
	const auto rotz = m * vector3(0.0f, 0.0f, 1.0f);
	CHECK_THAT(rotz, Vector3Matcher({1.0f, 0.0f, 0.0f}));

	{
		const matrix3 m2 = matrix3::rotation_around_y(angle);

		CHECK_THAT(m, Matrix3Matcher(m2));
	}
}

TEST_CASE("matrix3::rotation(euler_angles)_Z", "[working][unittest][matrix3]")
{
	const auto angle = radian(float(pi / 2.0));

	const euler_angles r{0_rad, 0_rad, angle};

	const matrix3 m = matrix3::rotation(r);
	const std::array<float, 9> a = m.to_array();

	const auto rotx = m * vector3(1.0f, 0.0f, 0.0f);
	CHECK_THAT(rotx, Vector3Matcher({0.0f, 1.0f, 0.0f}));
	const auto roty = m * vector3(0.0f, 1.0f, 0.0f);
	CHECK_THAT(roty, Vector3Matcher({-1.0f, 0.0f, 0.0f}));
	const auto rotz = m * vector3(0.0f, 0.0f, 1.0f);
	CHECK_THAT(rotz, Vector3Matcher({0.0f, 0.0f, 1.0f}));

	{
		const matrix3 m2 = matrix3::rotation_around_z(angle);

		CHECK_THAT(m, Matrix3Matcher(m2));
	}
}

TEST_CASE("matrix3::rotation(quaternion)", "[working][unittest][matrix3]")
{
	quaternion q{1.0f, 2.0f, 3.0f, 4.0f};
	q.normalize();

	const matrix3 m2 = matrix3::rotation(q);
	const matrix4 m4 = matrix4::rotation(q);

	matrix3 m3;
	m3.row(0) = m4.row(0).xyz();
	m3.row(1) = m4.row(1).xyz();
	m3.row(2) = m4.row(2).xyz();

	const std::array<float, 9> a = m2.to_array();
	const std::array<float, 9> b = m3.to_array();

	REQUIRE(a == b);
}

TEST_CASE("matrix3::scale(vector3)", "[working][unittest][matrix3]")
{
	vector3 s{0.5f, 0.4f, 0.3f};
	matrix3 m = matrix3::scale(s);

	REQUIRE(m.column(0).row(0) == s.x);
	REQUIRE(m.column(1).row(1) == s.y);
	REQUIRE(m.column(2).row(2) == s.z);

	REQUIRE(m.column(0).row(1) == 0.0f);
	REQUIRE(m.column(0).row(2) == 0.0f);

	REQUIRE(m.column(1).row(0) == 0.0f);
	REQUIRE(m.column(1).row(2) == 0.0f);

	REQUIRE(m.column(2).row(0) == 0.0f);
	REQUIRE(m.column(2).row(1) == 0.0f);
}

TEST_CASE("matrix3::scale(float, float, float)",
          "[working][unittest][matrix3]")
{
	float x{0.5f};
	float y{0.4f};
	float z{0.3f};
	matrix3 m = matrix3::scale(x, y, z);

	REQUIRE(m.column(0).row(0) == x);
	REQUIRE(m.column(1).row(1) == y);
	REQUIRE(m.column(2).row(2) == z);

	REQUIRE(m.column(0).row(1) == 0.0f);
	REQUIRE(m.column(0).row(2) == 0.0f);

	REQUIRE(m.column(1).row(0) == 0.0f);
	REQUIRE(m.column(1).row(2) == 0.0f);

	REQUIRE(m.column(2).row(0) == 0.0f);
	REQUIRE(m.column(2).row(1) == 0.0f);
}

TEST_CASE("matrix3::scale(float)", "[working][unittest][matrix3]")
{
	float s{0.5f};
	matrix3 m = matrix3::scale(s);

	REQUIRE(m.column(0).row(0) == s);
	REQUIRE(m.column(1).row(1) == s);
	REQUIRE(m.column(2).row(2) == s);

	REQUIRE(m.column(0).row(1) == 0.0f);
	REQUIRE(m.column(0).row(2) == 0.0f);

	REQUIRE(m.column(1).row(0) == 0.0f);
	REQUIRE(m.column(1).row(2) == 0.0f);

	REQUIRE(m.column(2).row(0) == 0.0f);
	REQUIRE(m.column(2).row(1) == 0.0f);
}

TEST_CASE("matrix3::rotation_around_x(radian)", "[working][unittest][matrix3]")
{
	using ::Catch::Matchers::WithinAbs;

	const direction3 axis = direction3::posX();

	const vector3 v0{0.0f, 1.0f, 0.0f};
	const vector3 vE{0.0f, 0.0f, 1.0f};

	{
		const radian rotation{0_rad};
		const matrix3 m = matrix3::rotation_around_x(rotation);
		const vector3 v1 = m * v0;
		CHECK(v1 == v0);
	}

	{
		const radian rotation{90_deg};
		const matrix3 m = matrix3::rotation_around_x(rotation);
		const vector3 v1 = m * v0;
		CHECK_THAT(v1, Vector3Matcher(vE));

		const matrix3 m2 =
		    matrix3::rotation(quaternion{axis, static_pi_fraction<1, 2>{}});
		const vector3 v2 = m2 * v0;
		CHECK_THAT(v2, Vector3Matcher(vE));

		CHECK_THAT(m, Matrix3Matcher(m2));

		const matrix3 m3 = matrix3::rotation(axis, rotation);
		const vector3 v3 = m3 * v0;
		CHECK_THAT(v3, Vector3Matcher(vE));

		CHECK_THAT(m2, Matrix3Matcher(m3));

		const matrix3 m4 = matrix3::rotation(axis, static_pi_fraction<1, 2>{});
		const vector3 v4 = m4 * v0;
		CHECK_THAT(v4, Vector3Matcher(vE));

		CHECK_THAT(m3, Matrix3Matcher(m4));
	}

	{
		const radian rotation{15_deg};
		const matrix3 m = matrix3::rotation_around_x(rotation);
		const matrix3 m2 = matrix3::rotation(quaternion{axis, rotation});
		CHECK_THAT(m, Matrix3Matcher(m2));
	}

	{
		const radian rotation{45_deg};
		const matrix3 m = matrix3::rotation_around_x(rotation);
		const matrix3 m2 = matrix3::rotation(quaternion{axis, rotation});
		CHECK_THAT(m, Matrix3Matcher(m2));
	}

	{
		const radian rotation{90_deg};
		const matrix3 m = matrix3::rotation_around_x(rotation);
		const matrix3 m2 = matrix3::rotation(quaternion{axis, rotation});
		CHECK_THAT(m, Matrix3Matcher(m2));
	}

	{
		const radian rotation{180_deg};
		const matrix3 m = matrix3::rotation_around_x(rotation);
		const matrix3 m2 = matrix3::rotation(quaternion{axis, rotation});
		CHECK_THAT(m, Matrix3Matcher(m2));
	}

	{
		const radian rotation{270_deg};
		const matrix3 m = matrix3::rotation_around_x(rotation);
		const matrix3 m2 = matrix3::rotation(quaternion{axis, rotation});
		CHECK_THAT(m, Matrix3Matcher(m2));
	}

	{
		const radian rotation{360_deg};
		const matrix3 m = matrix3::rotation_around_x(rotation);
		const matrix3 m2 = matrix3::rotation(quaternion{axis, rotation});
		CHECK_THAT(m, Matrix3Matcher(m2));
	}

	{
		const radian rotation{pi};
		const matrix3 m = matrix3::rotation_around_x(rotation);

		const vector3 v1 = m * v0;
		CHECK_THAT(v1, Vector3Matcher(-v0));
	}
}

TEST_CASE("matrix3::rotation_around_y(radian)", "[working][unittest][matrix3]")
{
	using ::Catch::Matchers::WithinAbs;

	const vector3 v0{0.0f, 0.0f, 1.0f};
	const vector3 vE{1.0f, 0.0f, 0.0f};

	{
		radian rotation{0_rad};
		matrix3 m = matrix3::rotation_around_y(rotation);
		vector3 v1 = m * v0;
		CHECK(v1 == v0);
	}

	const auto almost_0 = Approx(0.0f).margin(1.0e-6);

	{
		radian rotation{90_deg};
		matrix3 m = matrix3::rotation_around_y(rotation);
		vector3 v1 = m * v0;
		CHECK_THAT(v1, Vector3Matcher(vE));

		matrix3 m2 = matrix3::rotation(
		    quaternion{direction3::posY(), static_pi_fraction<1, 2>{}});
		CHECK_THAT(m, Matrix3Matcher(m2));
	}

	{
		radian rotation{pi};
		matrix3 m = matrix3::rotation_around_y(rotation);
		vector3 v1 = m * v0;
		CHECK_THAT(v1, Vector3Matcher(-v0));
	}
}

TEST_CASE("matrix3::rotation_around_z(radian)", "[working][unittest][matrix3]")
{
	using ::Catch::Matchers::WithinAbs;

	const vector3 v0{1.0f, 0.0f, 0.0f};
	const vector3 vE{0.0f, 1.0f, 0.0f};

	{
		const radian rotation{0_rad};
		const matrix3 m = matrix3::rotation_around_z(rotation);
		const vector3 v1 = m * v0;
		CHECK(v1 == v0);
	}

	{
		const radian rotation{90_deg};
		const matrix3 m = matrix3::rotation_around_z(rotation);
		const vector3 v1 = m * v0;
		CHECK_THAT(v1, Vector3Matcher(vE));

		const matrix3 m2 = matrix3::rotation(
		    quaternion{direction3::posZ(), static_pi_fraction<1, 2>{}});
		const vector3 v2 = m2 * v0;
		CHECK_THAT(v2, Vector3Matcher(vE));

		CHECK_THAT(m, Matrix3Matcher(m2));
	}

	{
		const radian rotation{pi};
		const matrix3 m = matrix3::rotation_around_z(rotation);
		const vector3 v1 = m * v0;
		CHECK_THAT(v1, Vector3Matcher(-v0));
	}
}

TEST_CASE("matrix3::rotation_around_z(complex)",
          "[working][unittest][matrix3]")
{
	const vector3 v0{1.0f, 0.0f, 0.0f};

	{
		normalized<complex> rotation{0_rad};
		matrix3 m = matrix3::rotation_around_z(rotation);
		vector3 v1 = m * v0;
		CHECK(v1 == v1);
	}

	const auto almost_0 = Approx(0.0f).margin(1.0e-6);

	{
		normalized<complex> rotation{90_deg};
		matrix3 m = matrix3::rotation_around_z(rotation);
		vector3 v1 = m * v0;
		CHECK(std::abs(v1.x) == almost_0);
		CHECK(v1.y == Approx(1.0f).margin(1.0e-6));
		CHECK(std::abs(v1.z) == almost_0);

		matrix3 m2 = matrix3::rotation(
		    quaternion{direction3::posZ(), static_pi_fraction<1, 2>{}});
		REQUIRE_THAT(m, Matrix3Matcher(m2));
	}

	{
		normalized<complex> rotation{1_pi};
		matrix3 m = matrix3::rotation_around_z(rotation);
		vector3 v1 = m * v0;
		CHECK(v1.x == Approx(-v0.x).margin(1.0e-6));
		CHECK(v1.y == Approx(-v0.y).margin(1.0e-6));
		CHECK(v1.z == Approx(-v0.z).margin(1.0e-6));
	}
}

TEST_CASE("matrix3::rows(vector3, vector3, vector3)",
          "[working][unittest][matrix3]")
{
	matrix3 m = matrix3::rows({1.0f, 2.0f, 3.0f}, {5.0f, 6.0f, 7.0f},
	                          {9.0f, 10.0f, 11.0f});
	REQUIRE(m.row(0).column(0) == 1.0f);
	REQUIRE(m.row(0).column(1) == 2.0f);
	REQUIRE(m.row(0).column(2) == 3.0f);
	REQUIRE(m.row(1).column(0) == 5.0f);
	REQUIRE(m.row(1).column(1) == 6.0f);
	REQUIRE(m.row(1).column(2) == 7.0f);
	REQUIRE(m.row(2).column(0) == 9.0f);
	REQUIRE(m.row(2).column(1) == 10.0f);
	REQUIRE(m.row(2).column(2) == 11.0f);
}

TEST_CASE("matrix3::columns(vector3, vector3, vector3)",
          "[working][unittest][matrix3]")
{
	matrix3 m = matrix3::columns({1.0f, 2.0f, 3.0f}, {5.0f, 6.0f, 7.0f},
	                             {9.0f, 10.0f, 11.0f});
	REQUIRE(m.column(0).row(0) == 1.0f);
	REQUIRE(m.column(0).row(1) == 2.0f);
	REQUIRE(m.column(0).row(2) == 3.0f);
	REQUIRE(m.column(1).row(0) == 5.0f);
	REQUIRE(m.column(1).row(1) == 6.0f);
	REQUIRE(m.column(1).row(2) == 7.0f);
	REQUIRE(m.column(2).row(0) == 9.0f);
	REQUIRE(m.column(2).row(1) == 10.0f);
	REQUIRE(m.column(2).row(2) == 11.0f);
}

TEST_CASE("matrix3::row(int)", "[working][unittest][matrix3]")
{
	vector3 v0{1.0f, 2.0f, 3.0f};
	vector3 v1{5.0f, 6.0f, 7.0f};
	vector3 v2{9.0f, 10.0f, 11.0f};

	matrix3 m = matrix3::rows(v0, v1, v2);

	REQUIRE(v0 == m.row(0));
	REQUIRE(v1 == m.row(1));
	REQUIRE(v2 == m.row(2));

	for (int i = 0; i < 3; ++i)
	{
		m.row(i) = vector3(-1.0f, -2.0f, -3.0f);
		REQUIRE(m.row(i).column(0) == -1.0f);
		REQUIRE(m.row(i).column(1) == -2.0f);
		REQUIRE(m.row(i).column(2) == -3.0f);
		REQUIRE(m.row(i) == vector3(-1.0f, -2.0f, -3.0f));
	}
}

TEST_CASE("matrix3::column(int)", "[working][unittest][matrix3]")
{
	vector3 v0{1.0f, 2.0f, 3.0f};
	vector3 v1{5.0f, 6.0f, 7.0f};
	vector3 v2{9.0f, 10.0f, 11.0f};

	matrix3 m = matrix3::columns(v0, v1, v2);

	REQUIRE(v0 == m.column(0));
	REQUIRE(v1 == m.column(1));
	REQUIRE(v2 == m.column(2));

	for (int i = 0; i < 3; ++i)
	{
		m.column(i) = vector3(-1.0f, -2.0f, -3.0f);
		REQUIRE(m.column(i).row(0) == -1.0f);
		REQUIRE(m.column(i).row(1) == -2.0f);
		REQUIRE(m.column(i).row(2) == -3.0f);
		REQUIRE(m.column(i) == vector3(-1.0f, -2.0f, -3.0f));
	}
}

TEST_CASE("matrix3::diag(int)", "[working][unittest][matrix3]")
{
	vector3 v0{1.0f, 2.0f, 3.0f};
	vector3 v1{5.0f, 6.0f, 7.0f};
	vector3 v2{9.0f, 10.0f, 11.0f};

	matrix3 m = matrix3::columns(v0, v1, v2);

	REQUIRE(v0.x == m.diag(0));
	REQUIRE(v1.y == m.diag(1));
	REQUIRE(v2.z == m.diag(2));

	m.diag(0) = -1.0f;
	m.diag(1) = -2.0f;
	m.diag(2) = -3.0f;

	REQUIRE(m.diag(0) == -1.0f);
	REQUIRE(m.diag(1) == -2.0f);
	REQUIRE(m.diag(2) == -3.0f);
	REQUIRE(m.diag() == vector3(-1.0f, -2.0f, -3.0f));
}

TEST_CASE("matrix3::diag()", "[working][unittest][matrix3]")
{
	vector3 v0{1.0f, 2.0f, 3.0f};
	vector3 v1{5.0f, 6.0f, 7.0f};
	vector3 v2{9.0f, 10.0f, 11.0f};

	matrix3 m = matrix3::columns(v0, v1, v2);

	REQUIRE(v0.x == m.diag().x);
	REQUIRE(v1.y == m.diag().y);
	REQUIRE(v2.z == m.diag().z);

	m.diag() = vector3(-1.0f, -2.0f, -3.0f);

	REQUIRE(m.diag().x == -1.0f);
	REQUIRE(m.diag().y == -2.0f);
	REQUIRE(m.diag().z == -3.0f);
	REQUIRE(m.diag() == vector3(-1.0f, -2.0f, -3.0f));
}

TEST_CASE("matrix3::is_orthogonal()", "[working][unittest][matrix3]")
{
	CHECK(false == matrix3::zero().is_orthogonal());
	const matrix3 m2 = matrix3::identity();
	CHECK(true == m2.is_orthogonal());
}

TEST_CASE("matrix3::inverse()", "[working][unittest][matrix3]")
{
	matrix3 m1;

	// scale
	{
		const vector3 scale{4589.0f, 132.015f, 0.00125f};
		const vector3 invScale{1.0f / scale.x, 1.0f / scale.y, 1.0f / scale.z};

		m1 = matrix3::scale(scale);
		m1.inverse();
		matrix3 m2 = matrix3::scale(invScale);

		CHECK_THAT(m1, Matrix3Matcher(m2));
	}

	const radian rotation = 90_deg;

	// rotation around x
	{
		const vector3 p0{1.f, 2.0f, 3.0f};

		m1 = matrix3::rotation_around_x(rotation);
		const vector3 p1 = m1 * p0;
		m1.inverse();
		const vector3 p2 = m1 * p1;

		CHECK_THAT(p0, Vector3Matcher(p2));
	}

	// rotation around y
	{
		const vector3 p0{1.f, 2.0f, 3.0f};

		m1 = matrix3::rotation_around_y(rotation);
		const vector3 p1 = m1 * p0;
		m1.inverse();
		const vector3 p2 = m1 * p1;

		CHECK_THAT(p0, Vector3Matcher(p2));
	}

	// rotation around z
	{
		const vector3 p0{1.f, 2.0f, 3.0f};

		m1 = matrix3::rotation_around_z(rotation);
		const vector3 p1 = m1 * p0;
		m1.inverse();
		const vector3 p2 = m1 * p1;

		CHECK_THAT(p0, Vector3Matcher(p2));
	}
}

TEST_CASE("matrix3::get_inversed()", "[working][unittest][matrix3]")
{
	const std::array<float, 9> a{1.0f,  4.0f, 2.0f, 3.0f, 8.0f,
	                             10.0f, 2.0f, 6.0f, 1.0f};

	matrix3 m1{a};
	const matrix3 m2{m1.get_inversed()};
	m1.inverse();

	REQUIRE(m1 == m2);
}

TEST_CASE("matrix3::transpose()", "[working][unittest][matrix3]")
{
	matrix3 m1;

	m1.column(0).row(0) = 1.0f;
	m1.column(0).row(1) = 2.0f;
	m1.column(0).row(2) = 3.0f;
	m1.column(1).row(0) = 5.0f;
	m1.column(1).row(1) = 6.0f;
	m1.column(1).row(2) = 7.0f;
	m1.column(2).row(0) = 9.0f;
	m1.column(2).row(1) = 10.0f;
	m1.column(2).row(2) = 11.0f;

	matrix3 m2 = m1;
	m2.transpose();

	REQUIRE(m2.row(0).column(0) == 1.0f);
	REQUIRE(m2.row(0).column(1) == 2.0f);
	REQUIRE(m2.row(0).column(2) == 3.0f);
	REQUIRE(m2.row(1).column(0) == 5.0f);
	REQUIRE(m2.row(1).column(1) == 6.0f);
	REQUIRE(m2.row(1).column(2) == 7.0f);
	REQUIRE(m2.row(2).column(0) == 9.0f);
	REQUIRE(m2.row(2).column(1) == 10.0f);
	REQUIRE(m2.row(2).column(2) == 11.0f);
}

TEST_CASE("matrix3::get_transposed()", "[working][unittest][matrix3]")
{
	matrix3 m1;

	m1.column(0).row(0) = 1.0f;
	m1.column(0).row(1) = 2.0f;
	m1.column(0).row(2) = 3.0f;
	m1.column(1).row(0) = 5.0f;
	m1.column(1).row(1) = 6.0f;
	m1.column(1).row(2) = 7.0f;
	m1.column(2).row(0) = 9.0f;
	m1.column(2).row(1) = 10.0f;
	m1.column(2).row(2) = 11.0f;

	matrix3 m2 = m1;
	m2.transpose();

	CHECK(m1.get_transposed() == m2);
}

TEST_CASE("matrix3::determinant()", "[working][unittest][matrix3]")
{
	const matrix3 m2 = matrix3(std::array<float, 9>{
	    1.0f, 4.0f, 2.0f, 3.0f, 8.0f, 10.0f, 2.0f, 6.0f, 1.0f});
	const float d = m2.determinant();
	CHECK(d == 20.0f);
}

TEST_CASE("matrix3::operator*(const matrix3&)", "[working][unittest][matrix3]")
{
	matrix3 m1;

	m1.column(0).row(0) = 1.0f;
	m1.column(0).row(1) = 2.0f;
	m1.column(0).row(2) = 3.0f;
	m1.column(1).row(0) = 4.0f;
	m1.column(1).row(1) = 5.0f;
	m1.column(1).row(2) = 6.0f;
	m1.column(2).row(0) = 7.0f;
	m1.column(2).row(1) = 8.0f;
	m1.column(2).row(2) = 9.0f;

	matrix3 m2 = m1 * m1;

	REQUIRE(m2.column(0).row(0) == 30.0f);
	REQUIRE(m2.column(0).row(1) == 36.0f);
	REQUIRE(m2.column(0).row(2) == 42.0f);
	REQUIRE(m2.column(1).row(0) == 66.0f);
	REQUIRE(m2.column(1).row(1) == 81.0f);
	REQUIRE(m2.column(1).row(2) == 96.0f);
	REQUIRE(m2.column(2).row(0) == 102.0f);
	REQUIRE(m2.column(2).row(1) == 126.0f);
	REQUIRE(m2.column(2).row(2) == 150.0f);

	vector3 s{123.0f, 456.0f, 789.0f};
	m1 = matrix3::scale(s);
	m2 = m1 * m1;

	REQUIRE(m2.column(0).row(0) == s.x * s.x);
	REQUIRE(m2.column(0).row(1) == 0.0f);
	REQUIRE(m2.column(0).row(2) == 0.0f);
	REQUIRE(m2.column(1).row(0) == 0.0f);
	REQUIRE(m2.column(1).row(1) == s.y * s.y);
	REQUIRE(m2.column(1).row(2) == 0.0f);
	REQUIRE(m2.column(2).row(0) == 0.0f);
	REQUIRE(m2.column(2).row(1) == 0.0f);
	REQUIRE(m2.column(2).row(2) == s.z * s.z);
}

TEST_CASE("matrix3::operator*(vector3)", "[working][unittest][matrix3]")
{
	const vector3 v0{1.0f, 2.0f, 3.0f};

	const vector3 s{123.0f, 456.0f, 789.0f};
	matrix3 m = matrix3::scale(s);

	vector3 v1 = m * v0;

	REQUIRE(v1.x == (v0.x * s.x));
	REQUIRE(v1.y == (v0.y * s.y));
	REQUIRE(v1.z == (v0.z * s.z));

	float scale{20.01f};
	m = matrix3::scale(scale);

	v1 = m * v0;

	REQUIRE(v1.x == (v0.x * scale));
	REQUIRE(v1.y == (v0.y * scale));
	REQUIRE(v1.z == (v0.z * scale));
}
