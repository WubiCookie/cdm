#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <common.hpp>

/*
TEST_CASE("transform3::operator*(transform3)",
          "[working][unittest][transform3]")
{
    // const transform3 tr1{
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

    const transform3 tr1{
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

    // const transform3 tr2{
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
    const transform3 tr2{tr1};

    // const transform3 tr1{
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
    // const transform3 tr2{
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

TEST_CASE("transform3::to_matrix()", "[working][unittest][transform3]")
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
