#include <common.hpp>

TEST_CASE("transform3d::to_matrix()", "[working][unittest][transform3d]")
{
	const transform3d tr{
	    .position = vector3{1, 2, 3},
	    .rotation = quaternion{direction{0.0f, 1.0f, 1.0f}, 74_deg},
	    .scale = vector3{4, 5, 6},
	};

	matrix4 Tr{matrix4::translation(tr.position)};
	matrix4 R{matrix4::rotation(tr.rotation)};
	matrix4 S{matrix4::scale(tr.scale)};
	const matrix4 m = Tr * S * R;
	const matrix4 trm = tr.to_matrix4();

	CHECK(trm == m);
}

TEST_CASE("transform3d::operator*(vector3)",
          "[working][unittest][transform3d]")
{
	const transform3d tr{
	    .position = vector3{1, 2, 3},
	    .rotation = quaternion{direction{0.0f, 1.0f, 1.0f}, 74_deg},
	    .scale = vector3{4, 5, 6},
	};
	const matrix4 trm = tr.to_matrix4();

	const float x = GENERATE(-10, -1, 0, 1, 10);
	const float y = GENERATE(-10, -1, 0, 1, 10);
	const float z = GENERATE(-10, -1, 0, 1, 10);

	const vector3 v0{x, y, z};
	const vector4 v1 = vector4(tr * v0, 1.0f);
	const vector4 v2 = trm * vector4(v0, 1.0f);

	CHECK_THAT(v1, Vector4Matcher(v2, 1.0e-5));
}
