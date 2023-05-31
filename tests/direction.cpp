#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <common.hpp>

TEST_CASE("direction3::constructors", "[working][unittest][direction3]")
{
	direction3 d1 = direction3(vector3{ 2, 0, 0 });
	CHECK(d1.x() == 1);
	CHECK(d1.y() == 0);
	CHECK(d1.z() == 0);

	direction3 d2{ vector3{ 0, 2, 0 } };
	CHECK(d2.x() == 0);
	CHECK(d2.y() == 1);
	CHECK(d2.z() == 0);

	direction3 d3{ 0, 0, 1 };
	CHECK(d3.x() == 0);
	CHECK(d3.y() == 0);
	CHECK(d3.z() == 1);
}

TEST_CASE("direction3::assignment", "[working][unittest][direction3]")
{
	direction3 d;

	d = direction3(vector3{ 0, 2, 0 });
	CHECK(d.x() == 0);
	CHECK(d.y() == 1);
	CHECK(d.z() == 0);
}

TEST_CASE("direction3::operator*()", "[working][unittest][direction3]")
{
	direction3 d{ 0, 2, 0 };

	const vector3& v = d;
	CHECK(v.x == 0);
	CHECK(v.y == 1);
	CHECK(v.z == 0);
}

TEST_CASE("direction3::operator vector3()", "[working][unittest][direction3]")
{
	direction3 d{ 0, 2, 0 };

	vector3 v = d;
	CHECK(v.x == 0);
	CHECK(v.y == 1);
	CHECK(v.z == 0);
}

TEST_CASE("direction3::operator-()", "[working][unittest][direction3]")
{
	CHECK(-direction3::posX() == direction3::negX());
	CHECK(-direction3::posY() == direction3::negY());
	CHECK(-direction3::posZ() == direction3::negZ());
}
