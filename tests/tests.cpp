#include <catch2/catch.hpp>
#include <cdm_maths.hpp>
#include <iostream>

using namespace cdm;
using namespace std;


TEST_CASE("matrix4", "[working]" "[unittest]") {

	float
		m00 = 0.0f, m10 = 0.0f, m20 = 0.0f, m30 = 0.0f,
		m01 = 0.0f, m11 = 0.0f, m21 = 0.0f, m31 = 0.0f,
		m02 = 0.0f, m12 = 0.0f, m22 = 0.0f, m32 = 0.0f,
		m03 = 0.0f, m13 = 0.0f, m23 = 0.0f, m33 = 0.0f;

	matrix4* m = new matrix4();
	matrix4* m1 = new matrix4(m00, m10, m20, m30,
		m01, m11, m21, m31,
		m02, m12, m22, m32,
		m03, m13, m23, m33);

	SECTION("constructor") {
		array<float, 16> a = m->operator array<float, 16>();
		array<float, 16> b = m1->operator array<float, 16>();
		CHECK(a == b);
	}

	SECTION("zero()") {
		matrix4 m2 = m1->zero();
		array<float, 16> a = m2.operator array<float, 16>();
		CHECK(m2.m00 == 0.0f); 
		array<float, 16> g = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
		CHECK(a == g);
	}

	SECTION("identity()") {
		matrix4 m2 = m1->identity();
		array<float, 16> a = m2.operator array<float, 16>();
		CHECK(m2.m00 == 1.00f);
		array<float, 16> g = { 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f };
		CHECK(a == g);
	}
	
	SECTION("rotation(euler_angles r)") {

		euler_angles r;
		r = { radian(1.5708f), radian(0), radian(0)};
		matrix4 m2 = m1->rotation(r);
		matrix3 m3;
		matrix4 m4 = matrix4(m3.rotation(r));
		array<float, 16> a = m2;
		array<float, 16> g = m4;
		CHECK(a == g);
		matrix3 RX = {1.0f, 0.0f, 0.0f, 0.0f, cos(r.x), sin(r.x), 0.0f, -sin(r.x), cos(r.x)};
		matrix3 RY = { cos(r.y), 0.0f, -sin(r.y), 0.0f, 1.0f, 0.0f, sin(r.y), 0.0f, cos(r.y)};
		matrix3 RZ = { cos(r.z), sin(r.z), 0.0f, -sin(r.z), cos(r.z), 0.0f, 0.0f, 0.0f, 1.0f };
		matrix3 m5 = RY * RX * RZ;
		matrix4 m6 = matrix4(m5);
		array<float, 16> m = m6;
		for (size_t i = 0; i < a.size(); i++)
			REQUIRE(a[i] == Approx(m[i]).epsilon(0.1));
	}

	SECTION("rotation(quaternion q)") {
		quaternion q;
		matrix4 m2 = m1->rotation(q);
		matrix3 m3;
		matrix4 m4 = matrix4(m3.rotation(q));
		array<float, 16> a = m2.operator array<float, 16>();
		array<float, 16> g = m4.operator array<float, 16>();
		CHECK(a == g);
	}

	SECTION("translation(vector3 t)") {
		vector3 t;
		matrix4 m2 = m1->translation(t);
		array<float, 16> a = m2.operator array<float, 16>();
		array<float, 16> g = { 1.0f, 0.0f, 0.0f, t.x, 0.0f, 1.0f, 0.0f, t.y, 0.0f, 0.0f, 1.0f, t.z, 0.0f, 0.0f, 0.0f, 1.0f };
		CHECK(a == g);

	}

	SECTION("translation(float x, float y, float z") {
		matrix4 m2 = m1->translation(0.02f, 0.01f, 0.01f);
		matrix4 m4 = m4.translation({ 0.02f, 0.01f, 0.01f });
		array<float, 16> a = m2.operator array<float, 16>();
		array<float, 16> g = m4.operator array<float, 16>();
		CHECK(a == g);
	}

	SECTION("scale(vector3 s") {
		vector3 s;
		s = { 0.5f, 0.0f, 0.0f };
		matrix4 m2 = m1->scale(s);
		array<float, 16> a = m2.operator array<float, 16>();
		array<float, 16> g = { s.x, 0.0f, 0.0f, 0.0f, 0.0f, s.y, 0.0f, 0.0f, 0.0f, 0.0f, s.z, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f };
		CHECK(s.x == 0.5f);
		CHECK(a == g);
	}

	SECTION("get_determinant()") {
		matrix4 m2 = { 1.0f, 4.0f, 2.0f, 3.0f, 8.0f, 10.0f, 2.0f, 6.0f, 1.0f, 7.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f };
		float d = m2.get_determinant();
		CHECK(d == -100);
	}

	SECTION("get_transposed()") {
		matrix4 m = { 1.0f, 4.0f, 2.0f, 3.0f, 8.0f, 10.0f, 2.0f, 6.0f, 1.0f, 7.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f };
		matrix4 t = { 1.0f, 8.0f, 1.0f, 2.0f, 4.0f, 10.0f, 7.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 3.0f, 6.0f, 2.0f, 2.0f };
		matrix4 m2 = m.get_transposed();
		array<float, 16> a = t.operator array<float, 16>();
		array<float, 16> g = m2.operator array<float, 16>();
		CHECK(a == g);
	}

	SECTION("get_inversed()") {
		matrix4 m = { 1.0f, 4.0f, 2.0f, 3.0f, 8.0f, 10.0f, 2.0f, 6.0f, 1.0f, 7.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f };
		matrix4 i = { -0.4f, 0.1f, -0.0f, 0.3f, -0.08f, 0.02f, 0.2f, -0.14f, -0.28f, -0.18f, 0.2f, 0.76f, 0.76f, 0.06f, -0.4f, -0.42f };
		matrix4 m2 = m.get_inversed();
		array<float, 16> a = i.operator array<float, 16>();
		array<float, 16> g = m2.operator array<float, 16>();
		REQUIRE(a.size() == g.size());
		for (size_t i = 0; i < a.size(); i++)
			REQUIRE(a[i] == Approx(g[i])); // Approx: comparing floating point values
	}

	SECTION("is_homogenous()") {
		CHECK(false == m1->is_homogenous());
		matrix4 m2 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f };
		CHECK(true == m2.is_homogenous());
	}

	SECTION("is_orthogonal()") {
		CHECK(false == m1->is_orthogonal());
		matrix4 m2 = {1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f};
		CHECK(true == m2.is_orthogonal());
	}

	SECTION("rotation around x()") {
		radian angle = radian(1.5708f);
		matrix3 m3;
		m3 = m3.rotation_around_x(angle);
		matrix3 m4 = { 1.0f, 0.0f, 0.0f, 0.0f, cos(angle), sin(angle), 0.0f, -sin(angle), cos(angle) };
		matrix4 m5 = matrix4(m4);
		matrix4 m6 = matrix4(m3);
		CHECK(m3.m11 == cos(angle));
		CHECK(m3.m22 == cos(angle));
		CHECK(m3.m12 == -sin(angle));
		CHECK(m3.m21 == sin(angle));
		array<float, 16> a = m6.operator array<float, 16>();
		array<float, 16> g = m5.operator array<float, 16>();
		CHECK(a == g);
	}

	delete m;
	delete m1;
}

TEST_CASE("matrix3", "[inProgress]") {

}