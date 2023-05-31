#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <common.hpp>

TEST_CASE("collides(segment3, plane)", "[working][unittest][segment3]")
{
	cdm::plane p{
	    {0.0f, 0.0f, 0.0f},
	    cdm::direction3::already_normalized(0.0f, 1.0f, 0.0f),
	};

	{
		cdm::segment3 s{
		    {0.0f, 0.0f, 0.0f},
		    {0.0f, 1.0f, 0.0f},
		};
		CHECK(cdm::collides(s, p).has_value() == true);
	}
	{
		cdm::segment3 s{
		    {0.0f, 0.1f, 0.0f},
		    {0.0f, 1.0f, 0.0f},
		};
		CHECK(cdm::collides(s, p).has_value() == false);
	}
	{
		cdm::segment3 s{
		    {0.0f, -1.0f, 0.0f},
		    {0.0f, 1.0f, 0.0f},
		};
		CHECK(cdm::collides(s, p).has_value() == true);
	}
	{
		cdm::segment3 s{
		    {-1.0f, -1.0f, -1.0f},
		    {1.0f, 1.0f, 1.0f},
		};
		CHECK(cdm::collides(s, p).has_value() == true);
	}
}

TEST_CASE("collides(segment3, plane, vector3)",
          "[working][unittest][segment3]")
{
	cdm::plane p{
	    {0.0f, 0.0f, 0.0f},
	    cdm::direction3::already_normalized(0.0f, 1.0f, 0.0f),
	};

	{
		cdm::segment3 s{
		    {0.0f, 0.0f, 0.0f},
		    {0.0f, 1.0f, 0.0f},
		};
		std::optional<cdm::vector3> v = cdm::collides(s, p);
		REQUIRE(v.has_value() == true);
		CHECK(v.value() == s.origin);
	}
	{
		cdm::segment3 s{
		    {0.0f, 0.1f, 0.0f},
		    {0.0f, 1.0f, 0.0f},
		};
		CHECK(cdm::collides(s, p).has_value() == false);
	}
	{
		cdm::segment3 s{
		    {0.0f, -1.0f, 0.0f},
		    {0.0f, 1.0f, 0.0f},
		};
		std::optional<cdm::vector3> v = cdm::collides(s, p);
		REQUIRE(v.has_value() == true);
		CHECK(v.value() == cdm::vector3{0.0f, 0.0f, 0.0f});
	}
	{
		cdm::segment3 s{
		    {-1.0f, -1.0f, -1.0f},
		    {1.0f, 1.0f, 1.0f},
		};
		std::optional<cdm::vector3> v = cdm::collides(s, p);
		REQUIRE(v.has_value() == true);
		CHECK(v.value() == cdm::vector3{0.0f, 0.0f, 0.0f});
	}
}
