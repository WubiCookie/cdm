#ifndef CDM_TESTS_COMMON_HPP
#define CDM_TESTS_COMMON_HPP 1

#include <catch2/catch.hpp>
#include <cdm_maths.hpp>

#include <array>
#include <iostream>

#define INFO_BEGIN(NAME)                                              \
	TEST_CASE("begin " #NAME, "[working][unittest][" #NAME "][info]") \
	{                                                                 \
		std::cout << "begin " #NAME "\n";                             \
	}

#define INFO_END(NAME)                                              \
	TEST_CASE("end " #NAME, "[working][unittest][" #NAME "][info]") \
	{                                                               \
		std::cout << #NAME " done!\n\n";                            \
	}

using namespace cdm;

struct Matrix3Matcher : Catch::Matchers::Impl::MatcherBase<matrix3>
{
	Matrix3Matcher(matrix3 target, double margin = 1.0e-6)
	    : m_target{target}, m_margin{margin}
	{
		CATCH_ENFORCE(margin >= 0,
		              "Invalid margin: " << margin << '.'
		                                 << " Margin has to be non-negative.");
	}
	bool match(const matrix3& matchee) const override
	{
		for (int c{0}; c < 3; ++c)
			for (int r{0}; r < 3; ++r)
			{
				float m = matchee.column(c).row(r);
				float t = m_target.column(c).row(r);
				bool valid = m + m_margin >= t && t + m_margin >= m;

				if (!valid)
					return false;
			}

		return true;
	}
	std::string describe() const override
	{
		return "\nis within " + ::Catch::Detail::stringify(m_margin) +
		       " of\n" + ::Catch::Detail::stringify(m_target);
	}

private:
	matrix3 m_target;
	double m_margin;
};

struct Matrix4Matcher : Catch::Matchers::Impl::MatcherBase<matrix4>
{
	Matrix4Matcher(matrix4 target, double margin = 1.0e-6)
	    : m_target{target}, m_margin{margin}
	{
		CATCH_ENFORCE(margin >= 0,
		              "Invalid margin: " << margin << '.'
		                                 << " Margin has to be non-negative.");
	}
	bool match(const matrix4& matchee) const override
	{
		for (int c{0}; c < 4; ++c)
			for (int r{0}; r < 4; ++r)
			{
				float m = matchee.column(c).row(r);
				float t = m_target.column(c).row(r);
				bool valid = m + m_margin >= t && t + m_margin >= m;

				if (!valid)
					return false;
			}

		return true;
	}
	std::string describe() const override
	{
		return "\nis within " + ::Catch::Detail::stringify(m_margin) +
		       " of\n" + ::Catch::Detail::stringify(m_target);
	}

private:
	matrix4 m_target;
	double m_margin;
};

struct Vector2Matcher : Catch::Matchers::Impl::MatcherBase<vector2>
{
	Vector2Matcher(vector2 target, double margin = 1.0e-6)
	    : m_target{target}, m_margin{margin}
	{
		CATCH_ENFORCE(margin >= 0,
		              "Invalid margin: " << margin << '.'
		                                 << " Margin has to be non-negative.");
	}
	bool match(const vector2& matchee) const override
	{
		bool valid = true;
		// clang-format off
		valid &= matchee.x + m_margin >= m_target.x && m_target.x + m_margin >= matchee.x;
		valid &= matchee.y + m_margin >= m_target.y && m_target.y + m_margin >= matchee.y;
		// clang-format on

		return valid;
	}
	std::string describe() const override
	{
		return "\nis within " + ::Catch::Detail::stringify(m_margin) +
		       " of\n" + ::Catch::Detail::stringify(m_target);
	}

private:
	vector2 m_target;
	double m_margin;
};

struct Vector3Matcher : Catch::Matchers::Impl::MatcherBase<vector3>
{
	Vector3Matcher(vector3 target, double margin = 1.0e-6)
	    : m_target{target}, m_margin{margin}
	{
		CATCH_ENFORCE(margin >= 0,
		              "Invalid margin: " << margin << '.'
		                                 << " Margin has to be non-negative.");
	}
	bool match(const vector3& matchee) const override
	{
		bool valid = true;
		// clang-format off
		valid &= matchee.x + m_margin >= m_target.x && m_target.x + m_margin >= matchee.x;
		valid &= matchee.y + m_margin >= m_target.y && m_target.y + m_margin >= matchee.y;
		valid &= matchee.z + m_margin >= m_target.z && m_target.z + m_margin >= matchee.z;
		// clang-format on

		return valid;
	}
	std::string describe() const override
	{
		return "\nis within " + ::Catch::Detail::stringify(m_margin) +
		       " of\n" + ::Catch::Detail::stringify(m_target);
	}

private:
	vector3 m_target;
	double m_margin;
};

struct Vector4Matcher : Catch::Matchers::Impl::MatcherBase<vector4>
{
	Vector4Matcher(vector4 target, double margin = 1.0e-6)
	    : m_target{target}, m_margin{margin}
	{
		CATCH_ENFORCE(margin >= 0,
		              "Invalid margin: " << margin << '.'
		                                 << " Margin has to be non-negative.");
	}
	bool match(const vector4& matchee) const override
	{
		bool valid = true;
		// clang-format off
		valid &= matchee.x + m_margin >= m_target.x && m_target.x + m_margin >= matchee.x;
		valid &= matchee.y + m_margin >= m_target.y && m_target.y + m_margin >= matchee.y;
		valid &= matchee.z + m_margin >= m_target.z && m_target.z + m_margin >= matchee.z;
		valid &= matchee.w + m_margin >= m_target.w && m_target.w + m_margin >= matchee.w;
		// clang-format on

		return valid;
	}
	std::string describe() const override
	{
		return "\nis within " + ::Catch::Detail::stringify(m_margin) +
		       " of\n" + ::Catch::Detail::stringify(m_target);
	}

private:
	vector4 m_target;
	double m_margin;
};

struct PerspectiveMatcher : Catch::Matchers::Impl::MatcherBase<perspective>
{
	PerspectiveMatcher(perspective target, double margin = 1.0e-6)
	    : m_target{target}, m_margin{margin}
	{
		CATCH_ENFORCE(margin >= 0,
		              "Invalid margin: " << margin << '.'
		                                 << " Margin has to be non-negative.");
	}
	bool match(const perspective& matchee) const override
	{
		bool valid = true;
		auto validate = [&](float a, float b)
		{ valid &= a + m_margin >= b && b + m_margin >= a; };

		validate(float(matchee.get_angle()), float(m_target.get_angle()));
		validate(matchee.get_far(), m_target.get_far());
		validate(matchee.get_near(), m_target.get_near());
		validate(matchee.get_ratio(), m_target.get_ratio());

		return valid;
	}
	std::string describe() const override
	{
		return "\nis within " + ::Catch::Detail::stringify(m_margin) +
		       " of\n" + ::Catch::Detail::stringify(m_target.to_matrix4());
	}

private:
	perspective m_target;
	double m_margin;
};

#endif  // CDM_TESTS_COMMON_HPP
