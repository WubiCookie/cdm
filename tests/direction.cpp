#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <common.hpp>

TEST_CASE("direction::constructors", "[working][unittest][direction]")
{
	direction d1 = direction(vector3{ 2, 0, 0 });
	CHECK(d1->x == 1);
	CHECK(d1->y == 0);
	CHECK(d1->z == 0);

	direction d2{ vector3{ 0, 2, 0 } };
	CHECK(d2->x == 0);
	CHECK(d2->y == 1);
	CHECK(d2->z == 0);

	direction d3{ 0, 0, 1 };
	CHECK(d3->x == 0);
	CHECK(d3->y == 0);
	CHECK(d3->z == 1);
}

TEST_CASE("direction::assignment", "[working][unittest][direction]")
{
	direction d;

	d = direction(vector3{ 0, 2, 0 });
	CHECK(d->x == 0);
	CHECK(d->y == 1);
	CHECK(d->z == 0);
}

TEST_CASE("direction::operator*()", "[working][unittest][direction]")
{
	direction d{ 0, 2, 0 };

	const vector3& v = *d;
	CHECK(v.x == 0);
	CHECK(v.y == 1);
	CHECK(v.z == 0);
}

TEST_CASE("direction::operator vector3()", "[working][unittest][direction]")
{
	direction d{ 0, 2, 0 };

	vector3 v = d;
	CHECK(v.x == 0);
	CHECK(v.y == 1);
	CHECK(v.z == 0);
}

TEST_CASE("direction::operator-()", "[working][unittest][direction]")
{
	CHECK(-direction::posX() == direction::negX());
	CHECK(-direction::posY() == direction::negY());
	CHECK(-direction::posZ() == direction::negZ());
}
