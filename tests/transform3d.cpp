#include <common.hpp>

TEST_CASE("transform3d::rotate(quaternion)",
          "[working][unittest][transform3d]")
{
	const float dx = GENERATE(-1, 0, 1);
	const float dy = GENERATE(-1, 0, 1);
	const float dz = GENERATE(-1, 0, 1);

	const float dx2 = GENERATE(-1, 0, 1);
	const float dy2 = GENERATE(-1, 0, 1);
	const float dz2 = GENERATE(-1, 0, 1);

	if ((dx != 0 && dy != 0 && dz != 0) && (dx2 != 0 && dy2 != 0 && dz2 != 0))
	{
		transform3d tr1{
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

TEST_CASE("translate_absolute(vector3), translate_relative(vector3)",
          "[working][unittest][transform3d]")
{
	const float dx = GENERATE(-1, 0, 1);
	const float dy = GENERATE(-1, 0, 1);
	const float dz = GENERATE(-1, 0, 1);

	if (dx != 0 && dy != 0 && dz != 0)
	{
		const transform3d tr{
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

			transform3d tr1{tr};
			tr1.translate_absolute(translation);
			matrix4 m2 = tr1.to_matrix4();

			// transform3d tr2 = transform3d::identity();
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
			transform3d tr1{tr};
			matrix4 rotation = matrix4::rotation(tr1.rotation);

			matrix4 m1 = rotation * matrix4::translation(translation) *
			             rotation.get_inversed() * tr1.to_matrix4();

			tr1.translate_relative(translation);
			matrix4 m2 = tr1.to_matrix4();

			// transform3d tr2 = transform3d::identity();
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

/*
TEST_CASE("transform3d::operator*(transform3d)",
          "[working][unittest][transform3d]")
{
    // const transform3d tr1{
    //    .position = vector3(GENERATE(-1, 0, 1),  //
    //                        GENERATE(-1, 0, 1),  //
    //                        GENERATE(-1, 0, 1)   //
    //                        ),
    //    .rotation =
    //        quaternion{
    //            direction(GENERATE(-1, 0, 1),  //
    //                      GENERATE(-1, 0, 1),  //
    //                      GENERATE(-1, 0, 1)   //
    //                      ),
    //            GENERATE(-180_deg, -145_deg, -90_deg, -45_deg, -15_deg,
    //            -0_deg,
    //                     0_deg, 15_deg, 45_deg, 90_deg, 145_deg, 180_deg,
    //                     200_deg, 280_deg, 360_deg),
    //        },
    //    .scale = vector3(GENERATE(-0.1, 0, 0.1, 0.2, 0.3, 1),  //
    //                     GENERATE(-0.1, 0, 0.1, 0.2, 0.3, 1),  //
    //                     GENERATE(-0.1, 0, 0.1, 0.2, 0.3, 1)   //
    //                     ),
    //};

    float dx = GENERATE(0);
    float dy = GENERATE(0, 1);
    float dz = GENERATE(1);
    degree a = GENERATE(-180_deg, -145_deg, -90_deg, -0_deg, 0_deg, 15_deg,
                        45_deg, 90_deg, 180_deg, 360_deg);

    const transform3d tr1{
        .position = vector3(1,  //
                            1,  //
                            1   //
                            ),
        .rotation =
            quaternion{
                direction(dx,  //
                          dy,  //
                          dz   //
                          ),
                a,
            },
        .scale = vector3(0,  //
                         1,  //
                         1   //
                         ),
    };

    // const transform3d tr2{
    //	.position = vector3(
    //		GENERATE(-10, -1, 0, 1, 10),//
    //		GENERATE(-10, -1, 0, 1, 10),//
    //		GENERATE(-10, -1, 0, 1, 10)	//
    //	),
    //	.rotation = quaternion{
    //		direction(
    //			GENERATE(-10, -1, 0, 1, 10),//
    //			GENERATE(-10, -1, 0, 1, 10),//
    //			GENERATE(-10, -1, 0, 1, 10)//
    //		),
    //		GENERATE(-180_deg, -145_deg, -90_deg, -45_deg, -15_deg, -0_deg,
    // 0_deg, 15_deg, 45_deg, 90_deg, 145_deg, 180_deg, 200_deg, 280_deg,
    // 360_deg),
    //	},
    //	.scale = vector3(
    //		GENERATE(-100, -1, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 1, 100),//
    //		GENERATE(-100, -1, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 1, 100),//
    //		GENERATE(-100, -1, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 1, 100) //
    //	),
    //};
    const transform3d tr2{tr1};

    // const transform3d tr1{
    //	.position = vector3(
    //		GENERATE(0, 1, 2),//
    //		GENERATE(0),//
    //		GENERATE(0)	//
    //	),
    //	.rotation = quaternion{
    //		direction(
    //			GENERATE(1),//
    //			GENERATE(0),//
    //			GENERATE(0)//
    //		),
    //		GENERATE(0_deg),
    //	},
    //	.scale = vector3(
    //		GENERATE(-100, -1, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 1, 100),//
    //		GENERATE(-100, -1, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 1, 100),//
    //		GENERATE(-100, -1, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 1, 100) //
    //	),
    //};
    //
    // const transform3d tr2{
    //	.position = vector3(
    //		GENERATE(0, 1, 2),//
    //		GENERATE(0),//
    //		GENERATE(0)	//
    //	),
    //	.rotation = quaternion{
    //		direction(
    //			GENERATE(1),//
    //			GENERATE(0),//
    //			GENERATE(0)//
    //		),
    //		GENERATE(0_deg),
    //	},
    //	.scale = vector3(
    //		GENERATE(-100, -1, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 1, 100),//
    //		GENERATE(-100, -1, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 1, 100),//
    //		GENERATE(-100, -1, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 1, 100) //
    //	),
    //};

    CHECK_THAT((tr1 * tr2).to_matrix4(),
               Matrix4Matcher((tr1.to_matrix4() * tr2.to_matrix4()), 1.0e-2));

    if (not Matrix4Matcher((tr1.to_matrix4() * tr2.to_matrix4()), 1.0e-2)
                .match((tr1 * tr2).to_matrix4()))
    {
        std::cout << dx << ", " << dy << ", " << dz << ", " << a << "\n";
        std::cout << tr1 << "\n";
        std::cout << tr1.to_matrix4() << "\n\n";
        std::cout << tr2 << "\n";
        std::cout << tr2.to_matrix4() << "\n\n";
        std::cout << tr1 * tr2 << "\n";
        std::cout << "fail\n";
    }
}
//*/

TEST_CASE("transform3d::to_matrix()", "[working][unittest][transform3d]")
{
	const float dx = GENERATE(-1, 0, 1);
	const float dy = GENERATE(-1, 0, 1);
	const float dz = GENERATE(-1, 0, 1);

	if (dx != 0 && dy != 0 && dz != 0)
	{
		const transform3d tr{
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
		                     GENERATE(-1, 0, 0.1, 1)   //
		                     ),
		};

		matrix4 Tr{matrix4::translation(tr.position)};
		matrix4 R{matrix4::rotation(tr.rotation)};
		matrix4 S{matrix4::scale(tr.scale)};
		const matrix4 m = Tr * S * R;
		const matrix4 trm = tr.to_matrix4();

		CHECK(trm == m);
	}
}

TEST_CASE("transform3d::operator*(vector3)",
          "[working][unittest][transform3d]")
{
	const float dx = GENERATE(-1, 0, 1);
	const float dy = GENERATE(-1, 0, 1);
	const float dz = GENERATE(-1, 0, 1);

	if (dx != 0 && dy != 0 && dz != 0)
	{
		const transform3d tr{
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
