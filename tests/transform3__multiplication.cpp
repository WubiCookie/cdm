#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <common.hpp>

TEST_CASE("transform3::operator*(vector3)", "[working][unittest][transform3]")
{
	const float dx = GENERATE(-1, 0, 1);
	const float dy = GENERATE(-1, 0, 1);
	const float dz = GENERATE(-1, 0, 1);

	if (dx != 0 && dy != 0 && dz != 0)
	{
		const transform3 tr{
		    .position = vector3(GENERATE(0, 1),  //
		                        GENERATE(0, 1),  //
		                        GENERATE(0, 1)   //
		                        ),
		    .rotation =
		        quaternion{
		            direction3(dx, dy, dz),
		            GENERATE(-180_deg, -145_deg, -90_deg, -0_deg, 0_deg,
		                     15_deg, 45_deg, 90_deg, 180_deg, 360_deg),
		        },
		    .scale = vector3(GENERATE(-1, 0, 0.1, 1),  //
		                     GENERATE(-1, 0, 0.1, 1),  //
		                     GENERATE(-1, 0, 0.1, 1)   //
		                     ),
		};
		const matrix4 trm = tr.to_matrix4();

		const float x = GENERATE(-1, 0, 1, 10);
		const float y = GENERATE(-1, 0, 1, 10);
		const float z = GENERATE(-1, 0, 1, 10);

		const vector3 v0{x, y, z};
		const vector4 v1 = vector4(tr * v0, 1.0f);
		const vector4 v2 = trm * vector4(v0, 1.0f);

		CHECK_THAT(v1, Vector4Matcher(v2, 1.0e-5));
	}
}
