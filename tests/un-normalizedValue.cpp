#include <common.hpp>

TEST_CASE("normalizedValue::normalizedValue()", "[working][unittest][normalizedValue]")
{
	normalizedValue<float> v;

	CHECK(v.value() == 0.0f);
}

TEST_CASE("normalizedValue::normalizedValue(value)", "[working][unittest][normalizedValue]")
{
	CHECK(normalizedValue<float>{0.0f}.value() == 0.0f);
	CHECK(normalizedValue<float>{0.5f}.value() == 0.5f);
	CHECK(normalizedValue<float>{1.0f}.value() == 1.0f);
	CHECK(normalizedValue<float>{-1.0f}.value() == 0.0f);
	CHECK(normalizedValue<float>{2.0f}.value() == 1.0f);
}

TEST_CASE("normalizedValue assignation", "[working][unittest][normalizedValue]")
{
	normalizedValue<float> v;
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

TEST_CASE("unnormalizedValue::unnormalizedValue(valueDomain[, value])", "[working][unittest][unnormalizedValue]")
{
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> v{ d };
		CHECK(v.value() == 1.0f);
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> v{ d, 0.0f };
		CHECK(v.value() == 1.0f);
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> v{ d, 1.0f };
		CHECK(v.value() == 1.0f);
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> v{ d, 2.0f };
		CHECK(v.value() == 2.0f);
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> v{ d, 1.5f };
		CHECK(v.value() == 1.5f);
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> v{ d, 1.984151f };
		CHECK(v.value() == 1.984151f);
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> v{ d, 3.0f };
		CHECK(v.value() == 2.0f);
	}
}

TEST_CASE("unnormalizedValue::setValue(value)", "[working][unittest][unnormalizedValue]")
{
	valueDomain<float> d{ 1.0f, 2.0f };
	unnormalizedValue<float> v{ d };
	
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

TEST_CASE("unnormalizedValue::unnormalizedValue(normalizedValue)", "[working][unittest][unnormalizedValue]")
{
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> v{ d, normalizedValue<float>{0.0f} };
		CHECK(v.value() == 1.0f);
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> v{ d, normalizedValue<float>{1.0f} };
		CHECK(v.value() == 2.0f);
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> v{ d, normalizedValue<float>{0.5f} };
		CHECK(v.value() == 1.5f);
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> v{ d, normalizedValue<float>{1.0f/3.0f} };
		CHECK(v.value() == Approx(1.33333333f));
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> v{ d, normalizedValue<float>{10.0f} };
		CHECK(v.value() == 2.0f);
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> v{ d, normalizedValue<float>{-10.0f} };
		CHECK(v.value() == 1.0f);
	}
}

TEST_CASE("normalizedValue::normalizedValue(unnormalizedValue)", "[working][unittest][normalizedValue]")
{
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> u{ d, 1.0f };
		normalizedValue<float> n{ u };
		CHECK(n.value() == 0.0f);
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> u{ d, 2.0f };
		normalizedValue<float> n{ u };
		CHECK(n.value() == 1.0f);
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> u{ d, 0.0f };
		normalizedValue<float> n{ u };
		CHECK(n.value() == 0.0f);
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> u{ d, 100.0f };
		normalizedValue<float> n{ u };
		CHECK(n.value() == 1.0f);
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> u{ d, 1.5f };
		normalizedValue<float> n{ u };
		CHECK(n.value() == 0.5f);
	}
	{
		valueDomain<float> d{ 1.0f, 2.0f };
		unnormalizedValue<float> u{ d, 1.333333f };
		normalizedValue<float> n{ u };
		CHECK(n.value() == Approx(0.3333333f));
	}
}
