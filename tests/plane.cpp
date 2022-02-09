#include <common.hpp>

TEST_CASE("plane::evaluate(vector3)", "[working][unittest][plane]")
{
	{
		const plane p{
		    {0, 0, 0},
		    direction::posY(),
		};

		CHECK(p.evaluate({0, 0, 0}) == 0.0f);
		CHECK(p.evaluate({1, 0, 0}) == 0.0f);
		CHECK(p.evaluate({0, 1, 0}) == 1.0f);
		CHECK(p.evaluate({0, 0, 1}) == 0.0f);
		CHECK(p.evaluate({1, 1, 1}) == 1.0f);
		CHECK(p.evaluate({0, 0.5f, 0}) == 0.5f);
	}
	{
		const plane p{
		    {0, 0, 0},
		    direction::negY(),
		};

		CHECK(p.evaluate({0, 0, 0}) == 0.0f);
		CHECK(p.evaluate({1, 0, 0}) == 0.0f);
		CHECK(p.evaluate({0, 1, 0}) == -1.0f);
		CHECK(p.evaluate({0, 0, 1}) == 0.0f);
		CHECK(p.evaluate({1, 1, 1}) == -1.0f);
		CHECK(p.evaluate({0, 0.5f, 0}) == -0.5f);
	}
	{
		const plane p{
		    {0, 1, 0},
		    direction::posY(),
		};

		CHECK(p.evaluate({0, 0, 0}) == -1.0f);
		CHECK(p.evaluate({1, 0, 0}) == -1.0f);
		CHECK(p.evaluate({0, 1, 0}) == 0.0f);
		CHECK(p.evaluate({0, 0, 1}) == -1.0f);
		CHECK(p.evaluate({1, 1, 1}) == 0.0f);
		CHECK(p.evaluate({0, 0.5f, 0}) == -0.5f);
	}
}

TEST_CASE("plane::project3d(vector3)", "[working][unittest][plane]")
{
	using namespace Catch::Generators;

	const plane p{
	    {
	        GENERATE(range(-1.0f, 1.0f, 0.2f)),
	        GENERATE(range(-1.0f, 1.0f, 0.2f)),
	        GENERATE(range(-1.0f, 1.0f, 0.2f)),
	    },
	    direction::posY(),
	};

	const vector3 v{
	    GENERATE(range(-2.0f, 2.0f, 0.5f)),
	    GENERATE(range(-2.0f, 2.0f, 0.5f)),
	    GENERATE(range(-2.0f, 2.0f, 0.5f)),
	};

	CHECK(distance_between(v, p.project3d(v)) == std::abs(p.evaluate(v)));
	CHECK(p.evaluate(p.project3d(v)) == Approx(0.0f).margin(1.0e-6));
}

TEST_CASE("plane::project2d(vector3, direction)", "[working][unittest][plane]")
{
	using namespace Catch::Generators;

	const plane p{
	    {0, 0, 0},
	    direction::posZ(),
	};
	const direction tangent = direction::posX();

	const float x = GENERATE(range(-2.0f, 2.0f, 0.1f));
	const float y = GENERATE(range(-2.0f, 2.0f, 0.1f));
	const float z = GENERATE(range(-2.0f, 2.0f, 0.1f));

	CHECK(p.project2d({x, y, z}, tangent) == vector2{x, y});
}

TEST_CASE("plane::unproject(vector2, direction)", "[working][unittest][plane]")
{
	using namespace Catch::Generators;

	const plane p{
	    {
	        0.0f,
	        0.0f,
	        GENERATE(range(-2.0f, 2.0f, 0.1f)),
	    },
	    direction::posZ(),
	};
	const direction tangent = direction::posX();

	const float x = GENERATE(range(-2.0f, 2.0f, 0.1f));
	const float y = GENERATE(range(-2.0f, 2.0f, 0.1f));

	CHECK(p.unproject({x, y}, tangent) == vector3{x, y, p.origin.z});
}
