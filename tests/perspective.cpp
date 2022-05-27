#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <common.hpp>

TEST_CASE("operator*(perspective, vector4), perspective::to_matrix",
          "[working][unittest][perspective]")
{
	const perspective p{1_pi / 2.0f, 1.0f, 0.01f, 100.0f};
	const matrix4 pm = p.to_matrix4();

	const float x = GENERATE(-1.0f, 0.0f, 1.0f);
	const float y = GENERATE(-1.0f, 0.0f, 1.0f);
	const float z = GENERATE(-1.0f, 0.0f, 1.0f);
	const float w = GENERATE(-1.0f, 0.0f, 1.0f);
	const vector4 v{x, y, z, w};

	CHECK((p * v) == (pm * v));
}

TEST_CASE("operator*(perspective, vector4), perspective::to_inverse_matrix",
          "[working][unittest][perspective]")
{
	const perspective p{1_pi / 2.0f, 1.0f, 0.01f, 100.0f};
	const matrix4 pim = p.to_matrix4().get_inversed();
	const matrix4 ipm = p.to_inverse_matrix4();

	CHECK(pim == ipm);

	const float x = GENERATE(-1.0f, 0.0f, 1.0f);
	const float y = GENERATE(-1.0f, 0.0f, 1.0f);
	const float z = GENERATE(-1.0f, 0.0f, 1.0f);
	const float w = GENERATE(-1.0f, 0.0f, 1.0f);
	const vector4 v{x, y, z, w};

	CHECK((pim * v) == (ipm * v));
}

TEST_CASE("perspective::set_far(float)", "[working][unittest][perspective]")
{
	perspective p{1_pi / 2.0f, 1.0f, 0.01f, 100.0f};

	auto v0 = p * vector4{0, 0, 100, 1};
	v0 /= v0.w;

	CHECK_THAT(v0, Vector4Matcher({0, 0, 1, 1}, 0.001));

	p.set_far(42.666f);
	
	v0 = p * vector4{0, 0, 42.666f, 1};
	v0 /= v0.w;

	CHECK_THAT(v0, Vector4Matcher({0, 0, 1, 1}, 0.001));
}

TEST_CASE("operator*(perspective, matrix4), perspective::to_inverse_matrix",
          "[working][unittest][perspective]")
{
	const perspective p{1_pi / 2.0f, 1.0f, 0.01f, 100.0f};
	matrix4 m(std::array<float, 16>{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16});

	CHECK(p * m == p.to_matrix4() * m);
}

TEST_CASE("operator*(matrix4, perspective), perspective::to_inverse_matrix",
          "[working][unittest][perspective]")
{
	const perspective p{1_pi / 2.0f, 1.0f, 0.01f, 100.0f};
	matrix4 m(std::array<float, 16>{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16});

	CHECK(m * p == m * p.to_matrix4());
}
