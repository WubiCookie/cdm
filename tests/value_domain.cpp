#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <common.hpp>

TEST_CASE("value_domain::value_domain()", "[working][unittest][value_domain]")
{
	value_domain<float> domain;

	CHECK(domain.min() == 0.0f);
	CHECK(domain.max() == 1.0f);
}

TEST_CASE("value_domain constructions/assignation",
          "[working][unittest][value_domain]")
{
	{
		value_domain domain{1.2f, 5.25f};
		domain = value_domain{10.2f, 50.25f};

		CHECK(domain.min() == 10.2f);
		CHECK(domain.max() == 50.25f);
	}
	{
		value_domain domain0{1.2f, 5.25f};
		value_domain domain1 = domain0;

		CHECK(domain0.min() == 1.2f);
		CHECK(domain0.max() == 5.25f);
		CHECK(domain1.min() == 1.2f);
		CHECK(domain1.max() == 5.25f);
	}
	{
		value_domain domain0{1.2f, 5.25f};
		value_domain domain1 = std::move(domain0);

		CHECK(domain1.min() == 1.2f);
		CHECK(domain1.max() == 5.25f);
	}
	{
		value_domain domain0{1.2f, 5.25f};
		value_domain domain1{10.2f, 50.25f};
		domain1 = std::move(domain0);

		CHECK(domain1.min() == 1.2f);
		CHECK(domain1.max() == 5.25f);
	}
}

TEST_CASE("value_domain::value_domain(T lim0, T lim1)",
          "[working][unittest][value_domain]")
{
	value_domain domain1{1.2f, 5.25f};

	CHECK(domain1.lim0() == 1.2f);
	CHECK(domain1.lim1() == 5.25f);
	CHECK(domain1.min() == 1.2f);
	CHECK(domain1.max() == 5.25f);

	value_domain domain2{5.25f, 1.2f};

	CHECK(domain2.lim0() == 5.25f);
	CHECK(domain2.lim1() == 1.2f);
	CHECK(domain2.min() == 1.2f);
	CHECK(domain2.max() == 5.25f);
}

TEST_CASE("value_domain::clamp(T value)", "[working][unittest][value_domain]")
{
	value_domain domain{1.2f, 5.25f};

	CHECK(domain.clamp(0.0f) == 1.2f);
	CHECK(domain.clamp(1.0f) == 1.2f);
	CHECK(domain.clamp(1.2f) == 1.2f);
	CHECK(domain.clamp(1.3f) == 1.3f);
	CHECK(domain.clamp(3.0f) == 3.0f);
	CHECK(domain.clamp(5.0f) == 5.0f);
	CHECK(domain.clamp(5.5f) == 5.25f);
	CHECK(domain.clamp(1000.0f) == 5.25f);
}

TEST_CASE("value_domain::range()", "[working][unittest][value_domain]")
{
	CHECK(value_domain{0.0f, 1.0f}.range() == 1.0f);
	CHECK(value_domain{0.0f, 2.0f}.range() == 2.0f);
	CHECK(value_domain{1.0f, 2.0f}.range() == 1.0f);
	CHECK(value_domain{-1.0f, 2.0f}.range() == 3.0f);
	CHECK(value_domain{1.0f, 0.0f}.range() == 1.0f);
	CHECK(value_domain{2.0f, 0.0f}.range() == 2.0f);
	CHECK(value_domain{2.0f, 1.0f}.range() == 1.0f);
	CHECK(value_domain{2.0f, -1.0f}.range() == 3.0f);
}

TEST_CASE("value_domain::lerp(normalized_value value)",
          "[working][unittest][value_domain]")
{
	value_domain domain1{0.0f, 1.0f};

	CHECK(domain1.lerp(normalized_value(0.0f)) == domain1.lim0());
	CHECK(domain1.lerp(normalized_value(1.0f)) == domain1.lim1());
	CHECK(domain1.lerp(normalized_value(0.25f)) == 0.25f);
	CHECK(domain1.lerp(normalized_value(0.5f)) == 0.5f);
	CHECK(domain1.lerp(normalized_value(0.75f)) == 0.75f);
	CHECK(domain1.lerp(normalized_value(-145600.75f)) == domain1.lim0());
	CHECK(domain1.lerp(normalized_value(1000.75f)) == domain1.lim1());

	value_domain domain2{1.0f, 0.0f};

	CHECK(domain2.lerp(normalized_value(0.0f)) == domain2.lim0());
	CHECK(domain2.lerp(normalized_value(1.0f)) == domain2.lim1());
	CHECK(domain2.lerp(normalized_value(0.25f)) == 0.75f);
	CHECK(domain2.lerp(normalized_value(0.5f)) == 0.5f);
	CHECK(domain2.lerp(normalized_value(0.75f)) == 0.25f);
	CHECK(domain2.lerp(normalized_value(-145600.75f)) == domain2.lim0());
	CHECK(domain2.lerp(normalized_value(1000.75f)) == domain2.lim1());
}
