#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <common.hpp>

TEST_CASE("translate_absolute(vector3), translate_relative(vector3)",
          "[working][unittest][transform3]")
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
		            direction(dx, dy, dz),
		            GENERATE(-180_deg, -145_deg, -90_deg, -0_deg, 0_deg,
		                     15_deg, 45_deg, 90_deg, 180_deg, 360_deg),
		        },
		    .scale = vector3(GENERATE(-1, 0, 0.1, 1),  //
		                     GENERATE(-1, 0, 0.1, 1),  //
		                     GENERATE(-1, 0, 0.1, 1))  //
		};

		const vector3 translation(GENERATE(-1, 0, 0.1f, 1),   //
		                          GENERATE(-1, 0, 0.1f, 1),   //
		                          GENERATE(-1, 0, 0.1f, 1));  //

		{
			matrix4 m1 = matrix4::translation(translation) * tr.to_matrix4();

			transform3 tr1{tr};
			tr1.translate_absolute(translation);
			matrix4 m2 = tr1.to_matrix4();

			// transform3 tr2 = transform3::identity();
			// tr2.position = vector3(GENERATE(0, 1),  //
			//                       GENERATE(0, 1),  //
			//                       GENERATE(0, 1)); //
			// matrix4 m3 = (tr2 * tr).to_matrix4();

			CHECK(m1 == m2);

			CHECK_THAT(m1 * vector4(0, 0, 0, 1),
			           Vector4Matcher(m2 * vector4{0, 0, 0, 1}, 1.0e-5));
			CHECK_THAT(
			    m1 * vector4(0, 0, 0, 1),
			    Vector4Matcher(vector4(tr1 * vector3{0, 0, 0}, 1), 1.0e-5));
			CHECK_THAT(
			    m2 * vector4(0, 0, 0, 1),
			    Vector4Matcher(vector4(tr1 * vector3{0, 0, 0}, 1), 1.0e-5));
			CHECK_THAT(m1 * vector4(1, 1, 1, 1),
			           Vector4Matcher(m2 * vector4{1, 1, 1, 1}, 1.0e-5));
			CHECK_THAT(
			    m1 * vector4(1, 1, 1, 1),
			    Vector4Matcher(vector4(tr1 * vector3{1, 1, 1}, 1), 1.0e-5));
			CHECK_THAT(
			    m2 * vector4(1, 1, 1, 1),
			    Vector4Matcher(vector4(tr1 * vector3{1, 1, 1}, 1), 1.0e-5));
		}

		{
			transform3 tr1{tr};
			matrix4 rotation = matrix4::rotation(tr1.rotation);

			matrix4 m1 = rotation * matrix4::translation(translation) *
			             rotation.get_inversed() * tr1.to_matrix4();

			tr1.translate_relative(translation);
			matrix4 m2 = tr1.to_matrix4();

			// transform3 tr2 = transform3::identity();
			// tr2.position = vector3(GENERATE(0, 1),  //
			//                       GENERATE(0, 1),  //
			//                       GENERATE(0, 1)); //
			// matrix4 m3 = (tr2 * tr).to_matrix4();

			CHECK_THAT(m1, Matrix4Matcher(m2, 1.0e-5));

			CHECK_THAT(m1 * vector4(0, 0, 0, 1),
			           Vector4Matcher(m2 * vector4{0, 0, 0, 1}, 1.0e-5));
			CHECK_THAT(
			    m1 * vector4(0, 0, 0, 1),
			    Vector4Matcher(vector4(tr1 * vector3{0, 0, 0}, 1), 1.0e-5));
			CHECK_THAT(
			    m2 * vector4(0, 0, 0, 1),
			    Vector4Matcher(vector4(tr1 * vector3{0, 0, 0}, 1), 1.0e-5));
			CHECK_THAT(m1 * vector4(1, 1, 1, 1),
			           Vector4Matcher(m2 * vector4{1, 1, 1, 1}, 1.0e-5));
			CHECK_THAT(
			    m1 * vector4(1, 1, 1, 1),
			    Vector4Matcher(vector4(tr1 * vector3{1, 1, 1}, 1), 1.0e-5));
			CHECK_THAT(
			    m2 * vector4(1, 1, 1, 1),
			    Vector4Matcher(vector4(tr1 * vector3{1, 1, 1}, 1), 1.0e-5));
		}
	}
}
