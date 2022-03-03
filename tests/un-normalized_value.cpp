#include <common.hpp>

TEST_CASE("normalized_value::normalized_value()", "[working][unittest][normalized_value]")
{
	normalized_value<float> v;

	CHECK(v.value() == 0.0f);
}

TEST_CASE("normalized_value::normalized_value(value)", "[working][unittest][normalized_value]")
{
	CHECK(normalized_value{0.0f}.value() == 0.0f);
	CHECK(normalized_value{0.5f}.value() == 0.5f);
	CHECK(normalized_value{1.0f}.value() == 1.0f);
	CHECK(normalized_value{-1.0f}.value() == 0.0f);
	CHECK(normalized_value{2.0f}.value() == 1.0f);
}

TEST_CASE("normalized_value assignation", "[working][unittest][normalized_value]")
{
	normalized_value<float> v;
	CHECK(v.value() == 0.0f);

	v = 1.0f;
	CHECK(v.value() == 1.0f);

	v = 12.0f;
	CHECK(v.value() == 1.0f);
	
	v = 0.5f;
	CHECK(v.value() == 0.5f);
	
	v = 0.0f;
	CHECK(v.value() == 0.0f);
	
	v = -1000.0f;
	CHECK(v.value() == 0.0f);
}

TEST_CASE("unnormalized_value::unnormalized_value(value_domain[, value])", "[working][unittest][unnormalized_value]")
{
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value v{ d };
		CHECK(v.value() == 1.0f);
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value v{ d, 0.0f };
		CHECK(v.value() == 1.0f);
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value v{ d, 1.0f };
		CHECK(v.value() == 1.0f);
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value v{ d, 2.0f };
		CHECK(v.value() == 2.0f);
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value v{ d, 1.5f };
		CHECK(v.value() == 1.5f);
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value v{ d, 1.984151f };
		CHECK(v.value() == 1.984151f);
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value v{ d, 3.0f };
		CHECK(v.value() == 2.0f);
	}
}

TEST_CASE("unnormalized_value::setValue(value)", "[working][unittest][unnormalized_value]")
{
	value_domain d{ 1.0f, 2.0f };
	unnormalized_value v{ d };
	
	v.setValue(1.0f);
	CHECK(v.value() == 1.0f);
	v.setValue(1.5f);
	CHECK(v.value() == 1.5f);
	v.setValue(2.0f);
	CHECK(v.value() == 2.0f);
	v.setValue(0.0f);
	CHECK(v.value() == 1.0f);
	v.setValue(-50.0f);
	CHECK(v.value() == 1.0f);
	v.setValue(1.984151f);
	CHECK(v.value() == 1.984151f);
}

TEST_CASE("unnormalized_value::unnormalized_value(normalized_value)", "[working][unittest][unnormalized_value]")
{
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value v{ d, normalized_value{0.0f} };
		CHECK(v.value() == 1.0f);
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value v{ d, normalized_value{1.0f} };
		CHECK(v.value() == 2.0f);
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value v{ d, normalized_value{0.5f} };
		CHECK(v.value() == 1.5f);
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value v{ d, normalized_value{1.0f/3.0f} };
		CHECK(v.value() == Approx(1.33333333f));
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value v{ d, normalized_value{10.0f} };
		CHECK(v.value() == 2.0f);
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value v{ d, normalized_value{-10.0f} };
		CHECK(v.value() == 1.0f);
	}
}

TEST_CASE("normalized_value::normalized_value(unnormalized_value)", "[working][unittest][normalized_value]")
{
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value u{ d, 1.0f };
		normalized_value n{ u };
		CHECK(n.value() == 0.0f);
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value u{ d, 2.0f };
		normalized_value n{ u };
		CHECK(n.value() == 1.0f);
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value u{ d, 0.0f };
		normalized_value n{ u };
		CHECK(n.value() == 0.0f);
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value u{ d, 100.0f };
		normalized_value n{ u };
		CHECK(n.value() == 1.0f);
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value u{ d, 1.5f };
		normalized_value n{ u };
		CHECK(n.value() == 0.5f);
	}
	{
		value_domain d{ 1.0f, 2.0f };
		unnormalized_value u{ d, 1.333333f };
		normalized_value n{ u };
		CHECK(n.value() == Approx(0.3333333f));
	}
}

TEST_CASE("un-normalized_value_misc",
          "[working][unittest][normalized_value]")
{
	{
		const value_domain d{1.0f, 2.0f};
		const unnormalized_value u0{d, 0.25f};
		const normalized_value n{u0};
		const unnormalized_value u1{d, n};
		CHECK(u0.value() == Approx(u1.value()).margin(1.0e-5));
	}
	{
		const value_domain d{42.0f, -47.3f};
		const unnormalized_value u0{d, 0.25f};
		const normalized_value n{u0};
		const unnormalized_value u1{d, n};
		CHECK(u0.value() == Approx(u1.value()).margin(1.0e-5));
	}
}
