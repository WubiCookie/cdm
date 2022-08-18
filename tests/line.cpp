#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <common.hpp>

TEST_CASE("line<SlopeIntercept>::line(T, T)", "[working][unittest][line]")
{
	using line_type = line<line_representation::SlopeIntercept>;
	{
		line_type l(0.0f, 0.0f);

		CHECK(l.slope == 0.0f);
		CHECK(l.y_intercept == 0.0f);
	}
	{
		line_type l(1.0f, 0.0f);

		CHECK(l.slope == 1.0f);
		CHECK(l.y_intercept == 0.0f);
	}
	{
		line_type l(0.0f, 1.0f);

		CHECK(l.slope == 0.0f);
		CHECK(l.y_intercept == 1.0f);
	}
	{
		line_type l(-1.0f, -1.0f);

		CHECK(l.slope == -1.0f);
		CHECK(l.y_intercept == -1.0f);
	}
}

TEST_CASE("line<SlopeIntercept>::line(vector2_t<T>)", "[working][unittest][line]")
{
	using line_type = line<line_representation::SlopeIntercept>;
	{
		line_type l(vector2{1.0f, 1.0f});

		CHECK(l.slope == 1.0f);
		CHECK(l.y_intercept == 0.0f);
	}
	{
		line_type l(vector2{2.0f, 2.0f});

		CHECK(l.slope == 1.0f);
		CHECK(l.y_intercept == 0.0f);
	}
	{
		line_type l(vector2{1.0f, 0.0f});

		CHECK(l.slope == 0.0f);
		CHECK(l.y_intercept == 0.0f);
	}
	{
		line_type l(vector2{2.0f, 0.0f});

		CHECK(l.slope == 0.0f);
		CHECK(l.y_intercept == 0.0f);
	}
	{
		line_type l(vector2{2.0f, 1.0f});

		CHECK(l.slope == 0.5f);
		CHECK(l.y_intercept == 0.0f);
	}
	{
		line_type l(vector2{-2.0f, 1.0f});

		CHECK(l.slope == -0.5f);
		CHECK(l.y_intercept == 0.0f);
	}
	{
		line_type l(vector2{1.0f, 2.0f});

		CHECK(l.slope == 2.0f);
		CHECK(l.y_intercept == 0.0f);
	}
	{
		line_type l(vector2{1.0f, -2.0f});

		CHECK(l.slope == -2.0f);
		CHECK(l.y_intercept == 0.0f);
	}
	{
		line_type l(vector2{-1.0f, 2.0f});

		CHECK(l.slope == -2.0f);
		CHECK(l.y_intercept == 0.0f);
	}
}

TEST_CASE("line<SlopeIntercept>::line(vector2_t<T>, vector2_t<T>)",
          "[working][unittest][line]")
{
	using line_type = line<line_representation::SlopeIntercept>;
	{
		line_type l({0.0f, 0.0f}, {1.0f, 1.0f});

		CHECK(l.slope == 1.0f);
		CHECK(l.y_intercept == 0.0f);
	}
	{
		line_type l({0.0f, 1.0f}, {1.0f, 1.0f});

		CHECK(l.slope == 0.0f);
		CHECK(l.y_intercept == 1.0f);
	}
	{
		line_type l({0.0f, 0.0f}, {1.0f, -1.0f});

		CHECK(l.slope == -1.0f);
		CHECK(l.y_intercept == 0.0f);
	}
}

TEST_CASE("line<SlopeIntercept>::resolve_for_x", "[working][unittest][line]")
{
	using line_type = line<line_representation::SlopeIntercept>;
	{
		line_type l(0.0f, 0.0f);

		CHECK(l.resolve_for_x(0.0f) == vector2(0.0f, 0.0f));
		CHECK(l.resolve_for_x(1.0f) == vector2(1.0f, 0.0f));
		CHECK(l.resolve_for_x(-1.0f) == vector2(-1.0f, 0.0f));
	}
	{
		line_type l(1.0f, 0.0f);

		CHECK(l.resolve_for_x(0.0f) == vector2(0.0f, 0.0f));
		CHECK(l.resolve_for_x(1.0f) == vector2(1.0f, 1.0f));
		CHECK(l.resolve_for_x(-1.0f) == vector2(-1.0f, -1.0f));
	}
	{
		line_type l(0.0f, 1.0f);

		CHECK(l.resolve_for_x(0.0f) == vector2(0.0f, 1.0f));
		CHECK(l.resolve_for_x(1.0f) == vector2(1.0f, 1.0f));
		CHECK(l.resolve_for_x(-1.0f) == vector2(-1.0f, 1.0f));
	}
	{
		line_type l(-1.0f, -1.0f);

		CHECK(l.resolve_for_x(0.0f) == vector2(0.0f, -1.0f));
		CHECK(l.resolve_for_x(1.0f) == vector2(1.0f, -2.0f));
		CHECK(l.resolve_for_x(-1.0f) == vector2(-1.0f, 0.0f));
	}
}

TEST_CASE("line<SlopeIntercept>::resolve_for_y", "[working][unittest][line]")
{
	using line_type = line<line_representation::SlopeIntercept>;
	{
		line_type l(1.0f, 0.0f);

		CHECK(l.resolve_for_y(0.0f) == vector2(0.0f, 0.0f));
		CHECK(l.resolve_for_y(1.0f) == vector2(1.0f, 1.0f));
		CHECK(l.resolve_for_y(-1.0f) == vector2(-1.0f, -1.0f));
	}
	{
		line_type l(-1.0f, -1.0f);

		CHECK(l.resolve_for_y(0.0f) == vector2(-1.0f, 0.0f));
		CHECK(l.resolve_for_y(1.0f) == vector2(-2.0f, 1.0f));
		CHECK(l.resolve_for_y(-1.0f) == vector2(0.0f, -1.0f));
	}
}
