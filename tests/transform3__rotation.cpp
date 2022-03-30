#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <common.hpp>

TEST_CASE("transform3::rotate(quaternion)",
          "[not-working][unittest][transform3]")
{
	const float dx = GENERATE(-1, 0, 1);
	const float dy = GENERATE(-1, 0, 1);
	const float dz = GENERATE(-1, 0, 1);

	const float dx2 = GENERATE(-1, 0, 1);
	const float dy2 = GENERATE(-1, 0, 1);
	const float dz2 = GENERATE(-1, 0, 1);

	if ((dx != 0 && dy != 0 && dz != 0) && (dx2 != 0 && dy2 != 0 && dz2 != 0))
	{
		transform3 tr1{
		    .position = vector3(GENERATE(0, 1),  //
		                        GENERATE(0, 1),  //
		                        GENERATE(0, 1)   //
		                        ),
		    .rotation =
		        quaternion{
		            direction(dx, dy, dz),
		            // GENERATE(-180_deg, -145_deg, -90_deg, -0_deg, 0_deg,
		            //         15_deg, 45_deg, 90_deg, 180_deg, 360_deg),
		            GENERATE(-180_deg, -90_deg, -0_deg, 0_deg, 90_deg, 180_deg,
		                     360_deg),
		        },
		    .scale = vector3(GENERATE(-1, 0, 0.1, 1),  //
		                     GENERATE(-1, 0, 0.1, 1),  //
		                     GENERATE(-1, 0, 0.1, 1))  //
		};

		const quaternion rotation =
		    quaternion(direction(dx2, dy2, dz2),
		               GENERATE(-180_deg, -90_deg, -0_deg, 0_deg, 90_deg,
		                        180_deg, 360_deg));

		// matrix4 rotation = matrix4::rotation(tr1.rotation);

		matrix4 m1 = matrix4::rotation(rotation) * tr1.to_matrix4();

		tr1.rotate(rotation);
		matrix4 m2 = tr1.to_matrix4();

		// CHECK_THAT(m1, Matrix4Matcher(m2, 1.0e-5));

		// CHECK_THAT(m1 * vector4(0, 0, 0, 1),
		//           Vector4Matcher(m2 * vector4{0, 0, 0, 1}, 1.0e-5));
		// CHECK_THAT(m1 * vector4(0, 0, 0, 1),
		//           Vector4Matcher(vector4(tr1 * vector3{0, 0, 0},
		//           1), 1.0e-5));
		// CHECK_THAT(m2 * vector4(0, 0, 0, 1),
		//           Vector4Matcher(vector4(tr1 * vector3{0, 0, 0},
		//           1), 1.0e-5));
		// CHECK_THAT(m1 * vector4(1, 1, 1, 1),
		//           Vector4Matcher(m2 * vector4{1, 1, 1, 1}, 1.0e-5));
		// CHECK_THAT(m1 * vector4(1, 1, 1, 1),
		//           Vector4Matcher(vector4(tr1 * vector3{1, 1, 1},
		//           1), 1.0e-5));
		// CHECK_THAT(m2 * vector4(1, 1, 1, 1),
		//           Vector4Matcher(vector4(tr1 * vector3{1, 1, 1},
		//           1), 1.0e-5));
	}
}
