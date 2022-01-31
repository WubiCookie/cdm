#include <common.hpp>

TEST_CASE("valueDomain::valueDomain()", "[working][unittest][valueDomain]")
{
	valueDomain<float> domain;

	CHECK(domain.min() == 0.0f);
	CHECK(domain.max() == 1.0f);
}

TEST_CASE("valueDomain constructions/assignation",
          "[working][unittest][valueDomain]")
{
	{
		valueDomain<float> domain{1.2f, 5.25f};
		domain = valueDomain<float>{10.2f, 50.25f};

		CHECK(domain.min() == 10.2f);
		CHECK(domain.max() == 50.25f);
	}
	{
		valueDomain<float> domain0{1.2f, 5.25f};
		valueDomain<float> domain1 = domain0;

		CHECK(domain0.min() == 1.2f);
		CHECK(domain0.max() == 5.25f);
		CHECK(domain1.min() == 1.2f);
		CHECK(domain1.max() == 5.25f);
	}
	{
		valueDomain<float> domain0{1.2f, 5.25f};
		valueDomain<float> domain1 = std::move(domain0);

		CHECK(domain1.min() == 1.2f);
		CHECK(domain1.max() == 5.25f);
	}
	{
		valueDomain<float> domain0{1.2f, 5.25f};
		valueDomain<float> domain1{10.2f, 50.25f};
		domain1 = std::move(domain0);

		CHECK(domain1.min() == 1.2f);
		CHECK(domain1.max() == 5.25f);
	}
}

TEST_CASE("valueDomain::valueDomain(T min, T max)",
          "[working][unittest][valueDomain]")
{
	valueDomain<float> domain{1.2f, 5.25f};

	CHECK(domain.min() == 1.2f);
	CHECK(domain.max() == 5.25f);
}

TEST_CASE("valueDomain::clamp(T value)", "[working][unittest][valueDomain]")
{
	valueDomain<float> domain{1.2f, 5.25f};

	CHECK(domain.clamp(0.0f) == 1.2f);
	CHECK(domain.clamp(1.0f) == 1.2f);
	CHECK(domain.clamp(1.2f) == 1.2f);
	CHECK(domain.clamp(1.3f) == 1.3f);
	CHECK(domain.clamp(3.0f) == 3.0f);
	CHECK(domain.clamp(5.0f) == 5.0f);
	CHECK(domain.clamp(5.5f) == 5.25f);
	CHECK(domain.clamp(1000.0f) == 5.25f);
}

TEST_CASE("valueDomain::range()", "[working][unittest][valueDomain]")
{
	CHECK(valueDomain<float>{0.0f, 1.0f}.range() == 1.0f);
	CHECK(valueDomain<float>{0.0f, 2.0f}.range() == 2.0f);
	CHECK(valueDomain<float>{1.0f, 2.0f}.range() == 1.0f);
	CHECK(valueDomain<float>{-1.0f, 2.0f}.range() == 3.0f);
}

TEST_CASE("valueDomain::lerp(normalizedValue value)",
          "[working][unittest][valueDomain]")
{
	valueDomain<float> domain{0.0f, 1.0f};

	CHECK(domain.lerp(normalizedValue<float>(0.0f)) == domain.min());
	CHECK(domain.lerp(normalizedValue<float>(1.0f)) == domain.max());
	CHECK(domain.lerp(normalizedValue<float>(0.25f)) == 0.25f);
	CHECK(domain.lerp(normalizedValue<float>(0.5f)) == 0.5f);
	CHECK(domain.lerp(normalizedValue<float>(0.75f)) == 0.75f);
	CHECK(domain.lerp(normalizedValue<float>(-145600.75f)) == domain.min());
	CHECK(domain.lerp(normalizedValue<float>(1000.75f)) == domain.max());
}
