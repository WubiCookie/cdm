#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <common.hpp>

TEST_CASE("segment2::length()", "[working][unittest][segment2]")
{
	{
		segment2 s{
		    .origin{0.0f, 0.0f},
		    .end{1.0f, 0.0f},
		};

		CHECK(s.length() == 1.0f);
	}
}

TEST_CASE("collides(segment2, segment2)", "[working][unittest][segment2]")
{
	{
		segment2 s0{
		    .origin{-1.0f, 0.0f},
		    .end{1.0f, 0.0f},
		};
		segment2 s1{
		    .origin{0.0f, -1.0f},
		    .end{0.0f, 1.0f},
		};

		auto [p0, p1] = collides(s0, s1);

		REQUIRE(p0.has_value());
		CHECK(p0.value() == vector2{0.0f, 0.0f});
	}
}
