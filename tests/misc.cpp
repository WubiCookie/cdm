#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <common.hpp>

TEST_CASE("sign", "[working][unittest][misc]")
{
	CHECK(cdm::sign(0.0f) == 1);
	CHECK(cdm::sign(-0.0f) == 1);
	CHECK(cdm::sign(1.0f) == 1);
	CHECK(cdm::sign(-1.0f) == -1);
}

TEST_CASE("signnum", "[working][unittest][misc]")
{
	CHECK(cdm::signnum(0.0f) == 0);
	CHECK(cdm::signnum(-0.0f) == 0);
	CHECK(cdm::signnum(1.0f) == 1);
	CHECK(cdm::signnum(-1.0f) == -1);
}
