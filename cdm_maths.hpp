/* cdm_maths v3.0.0
   C++20 geometric library
   https://github.com/WubiCookie/cdm
   no warranty implied; use at your own risk

LICENSE

       DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
                   Version 2, December 2004

Copyright (C) 2022 Charles Seizilles de Mazancourt <charles DOT de DOT
mazancourt AT hotmail DOT fr>

Everyone is permitted to copy and distribute verbatim or modified
copies of this license document, and changing it is allowed as long
as the name is changed.

           DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
  TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

 0. You just DO WHAT THE FUCK YOU WANT TO.

CREDITS

Written by Charles Seizilles de Mazancourt
*/

#ifndef CDM_MATHS_HPP
#define CDM_MATHS_HPP 1

#include <cdm_concepts.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif
#ifdef near
#undef near
#endif
#ifdef far
#undef far
#endif

namespace cdm
{
#pragma region constants declarations
// clang-format off
constexpr uint64_t version_major = 3;
constexpr uint64_t version_minor = 0;
constexpr uint64_t version_patch = 0;

using std::numbers::pi;
using std::numbers::inv_pi;
using std::numbers::inv_sqrtpi;
constexpr double deg_to_rad = pi / 180.0;
constexpr double rad_to_deg = 180.0 / pi;
using std::numbers::sqrt2;
using std::numbers::sqrt3;
constexpr double inv_sqrt2 = 0.7071067811865475244008443621048490392848359376884740365883398689;
using std::numbers::inv_sqrt3;
constexpr double sqrt3_over_2 = 0.8660254037844386467637231707529361834714026269051903140279034897;
constexpr double epsilon = 1.0e-05;
// clang-format on
#pragma endregion

#pragma region type_traits_and_concepts
// clang-format off
namespace detail
{
template <typename T> struct is_vector : std::false_type {};
template <typename T> inline constexpr bool is_vector_v = is_vector<T>::value;

template <typename T> struct is_matrix : std::false_type {};
template <typename T> inline constexpr bool is_matrix_v = is_matrix<T>::value;
}  // namespace detail

template <typename T> concept vector = detail::is_vector_v<T>;
template <typename T> concept matrix = detail::is_matrix_v<T>;
template <typename T> concept normalizable = requires(T& t, const T& ct)
{
	{ t.normalize() } -> std::convertible_to<T>;
	{ t.get_normalized() } -> std::convertible_to<T>;
	{ ct.get_normalized() } -> std::convertible_to<T>;
};
// clang-format on
#pragma endregion

#pragma region forward_declarations
template <normalizable T>
class normalized;
template <arithmetic T>
struct complex_T;
template <arithmetic T>
struct radian_T;
template <arithmetic T>
struct degree_T;
template <std::signed_integral T>
struct pi_fraction_T;
template <std::signed_integral T, T NumeratorT, T DenominatorT>
struct static_pi_fraction_T;
template <arithmetic T>
struct vector2_T;
template <arithmetic T>
struct vector3_T;
template <arithmetic T>
struct vector4_T;
template <arithmetic T>
struct matrix2_T;
template <arithmetic T>
struct matrix3_T;
template <arithmetic T>
struct matrix4_T;
template <arithmetic T>
struct perspective_T;
template <arithmetic T>
struct direction2_T;
template <arithmetic T>
struct direction3_T;
template <arithmetic T>
struct euler_angles_T;
template <arithmetic T>
struct quaternion_T;
enum class line_representation
{
	SlopeIntercept,
	Points,
	PointAngleDegree,
	PointAngleRadian,
	PointAnglePiFraction,
};
template <arithmetic T, line_representation representation>
struct line_T;
template <arithmetic T>
struct segment2_T;
template <arithmetic T>
struct segment3_T;
template <arithmetic T>
struct plane_T;
template <arithmetic T>
struct oriented_plane_T;
template <arithmetic T>
struct circle_T;
template <arithmetic T>
struct ray2_T;
template <arithmetic T>
struct ray3_T;
template <arithmetic T>
struct aabb2_T;
template <arithmetic T>
struct aabb3_T;
template <arithmetic T>
struct transform2_T;
template <arithmetic T>
struct transform3_T;
template <arithmetic T>
struct uniform_transform2_T;
template <arithmetic T>
struct uniform_transform3_T;
template <arithmetic T>
struct unscaled_transform2_T;
template <arithmetic T>
struct unscaled_transform3_T;
template <typename T>
    requires arithmetic<T> || vector<T>
class value_domain_T;
template <typename T>
class unnormalized_value_T;
template <typename T>
class normalized_value_T;
#pragma endregion

#pragma region misc
namespace detail
{
template <int x, int y, arithmetic T>
T get_quaternion_t_matrix_element(quaternion_T<T> q)
{
	static_assert(x < 3, "x must be [0;2]");
	static_assert(y < 3, "y must be [0;2]");
	if constexpr (y == 0)
	{
		if constexpr (x == 0)
			return T(1) - T(2) * (q.y * q.y + q.z * q.z);
		if constexpr (x == 1)
			return T(2) * (q.x * q.y - q.z * q.w);
		if constexpr (x == 2)
			return T(2) * (q.x * q.z + q.y * q.w);
	}
	if constexpr (y == 1)
	{
		if constexpr (x == 0)
			return T(2) * (q.x * q.y + q.z * q.w);
		if constexpr (x == 1)
			return T(1) - T(2) * (q.x * q.x + q.z * q.z);
		if constexpr (x == 2)
			return T(2) * (q.y * q.z - q.x * q.w);
	}
	if constexpr (y == 2)
	{
		if constexpr (x == 0)
			return T(2) * (q.x * q.z - q.y * q.w);
		if constexpr (x == 1)
			return T(2) * (q.y * q.z + q.x * q.w);
		if constexpr (x == 2)
			return T(1) - T(2) * (q.x * q.x + q.y * q.y);
	}
}
}  // namespace detail

template <typename Functor, arithmetic T>
std::vector<vector3_T<T>> function2D_sampler(const Functor& functor,
                                             T min,
                                             T max,
                                             T step)
{
	std::vector<vector3_T<T>> res;

	if (min > max)
		std::swap(min, max);

	if (step < epsilon)
		step = epsilon;

	for (T f = min; f < max; f += step)
	{
		res.push_back(functor(f));
	}

	return res;
}

template <typename Functor, arithmetic T>
std::vector<std::vector<vector3_T<T>>> function3D_sampler(
    const Functor& functor,
    vector2_T<T> min,
    vector2_T<T> max,
    vector2_T<T> step)
{
	std::vector<std::vector<vector3_T<T>>> res;

	if (min.x > max.x)
		std::swap(min.x, max.x);
	if (min.y > max.y)
		std::swap(min.y, max.y);

	if (step.x < epsilon)
		step.x = epsilon;
	if (step.y < epsilon)
		step.y = epsilon;

	size_t xCount = 0;
	size_t yCount = 0;

	for (T x = min.x; x < max.x; x += step.x)
		xCount++;
	for (T y = min.y; y < max.y; y += step.y)
		yCount++;
	res.reserve(yCount);

	for (T y = min.y; y < max.y; y += step.y)
	{
		std::vector<vector3_T<T>> row;
		row.reserve(xCount);
		for (T x = min.x; x < max.x; x += step.x)
		{
			row.push_back(functor(x, y));
		}
		res.push_back(std::move(row));
	}

	return res;
}

// template <arithmetic T>
// constexpr T lerp(T begin, T end, T percent)
//{
//	return std::lerp(begin, end, percent);
// }

template <typename T>
    requires arithmetic<T> || vector<T>
constexpr T clamp(T value, T min, T max)
{
	if constexpr (arithmetic<T>)
		return std::clamp(value, min, max);
	else
		return value.clamp(min, max);
}

template <arithmetic T>
bool nearly_equal(T f1, T f2, T e = T(epsilon))
{
	return std::abs(f1 - f2) < e;
}

template <arithmetic T>
constexpr int sign(T val)
{
	return (T(0) <= val) - (val < T(0));
}

template <arithmetic T>
constexpr int signnum(T val)
{
	return (T(0) < val) - (val < T(0));
}

template <arithmetic T>
constexpr T element_wise_min(T v0, T v1)
{
	return std::min(v0, v1);
}

template <arithmetic T>
constexpr T element_wise_max(T v0, T v1)
{
	return std::max(v0, v1);
}

template <arithmetic T>
constexpr T element_wise_abs(T v)
{
	return std::abs(v);
}

template <arithmetic T>
constexpr T element_wise_lerp(T begin, T end, T percent)
{
	return std::lerp(begin, end, percent);
}

template <typename T>
    requires arithmetic<T> || vector<T> || matrix<T>
constexpr T zero()
{
	if constexpr (arithmetic<T>)
		return T(0);
	else
		return T::zero();
}

template <typename T>
    requires arithmetic<T> || vector<T> || matrix<T>
constexpr T one()
{
	if constexpr (arithmetic<T>)
		return T(1);
	else
		return T::one();
}
#pragma endregion

#pragma region type_traits_and_concepts_details
// clang-format off
namespace detail
{
template <typename T> struct is_vector<vector2_T<T>> : std::true_type {};
template <typename T> struct is_vector<vector4_T<T>> : std::true_type {};
template <typename T> struct is_vector<vector3_T<T>> : std::true_type {};

template <typename T> struct is_vector<direction2_T<T>> : std::true_type {};
template <typename T> struct is_vector<direction3_T<T>> : std::true_type {};

template <typename T> struct is_matrix<matrix2_T<T>> : std::true_type {};
template <typename T> struct is_matrix<matrix3_T<T>> : std::true_type {};
template <typename T> struct is_matrix<matrix4_T<T>> : std::true_type {};
}  // namespace detail
//  clang-format on
#pragma endregion

#pragma region declarations
#pragma region declaration_complex_T
template <arithmetic T>
struct complex_T
{
	T r;
	T i;

	complex_T() = default;
	complex_T(T r_, T i_) : r{r_}, i{i_} {}
	explicit complex_T(radian_T<T> angle);
	complex_T(const complex_T&) = default;
	complex_T(complex_T&&) = default;
	complex_T& operator=(const complex_T&) = default;
	complex_T& operator=(complex_T&&) = default;

	T norm() const;
	T norm_squared() const;

	complex_T& normalize();
	complex_T get_normalized() const;

	complex_T& conjugate();
	complex_T get_conjugated() const;

	complex_T operator+(complex_T c) const;
	complex_T operator-(complex_T c) const;

	complex_T& operator+=(complex_T c);
	complex_T& operator-=(complex_T c);
	complex_T& operator*=(complex_T c);

	static complex_T identity();

	using underlying_type = T;
};

template <arithmetic T>
complex_T<T> operator*(complex_T<T> c1, complex_T<T> c2);
template <arithmetic T>
complex_T<T> operator*(complex_T<T> c, T f);
template <arithmetic T>
complex_T<T> operator*(normalized<complex_T<T>> c1, complex_T<T> c2);
template <arithmetic T>
complex_T<T> operator*(complex_T<T> c1, normalized<complex_T<T>> c2);
template <arithmetic T>
normalized<complex_T<T>> operator*(normalized<complex_T<T>> c1,
                                   normalized<complex_T<T>> c2);
template <arithmetic T>
vector2_T<T> operator*(normalized<complex_T<T>> c, vector2_T<T> v);

template <arithmetic T>
T norm(complex_T<T> c);
template <arithmetic T>
T norm_squared(complex_T<T> c);
template <arithmetic T>
complex_T<T> normalize(complex_T<T> c);
template <arithmetic T>
complex_T<T> conjugate(complex_T<T> c);
#pragma endregion

#pragma region declaration_radian_T
template <arithmetic T>
struct radian_T
{
private:
	T angle;

public:
	radian_T() = default;
	constexpr explicit radian_T(T f);
	radian_T(const radian_T& r) = default;
	radian_T(radian_T&& r) = default;
	constexpr radian_T(degree_T<T> d);

	explicit operator T() const;

	template <arithmetic U>
	explicit operator radian_T<U>() const
	{
		return radian_T<U>{U(angle)};
	}

	radian_T& operator+=(radian_T r);
	radian_T& operator-=(radian_T r);
	radian_T& operator*=(radian_T r);
	radian_T& operator/=(radian_T r);

	radian_T operator-() const;

	radian_T& operator=(const radian_T& r) = default;
	radian_T& operator=(degree_T<T> d);

	using underlying_type = T;
};

template <arithmetic T>
radian_T<T> operator+(radian_T<T>, radian_T<T>);
template <arithmetic T>
radian_T<T> operator-(radian_T<T>, radian_T<T>);
template <arithmetic T>
radian_T<T> operator*(T, radian_T<T>);
template <arithmetic T>
radian_T<T> operator*(radian_T<T>, T);
template <arithmetic T>
radian_T<T> operator/(radian_T<T>, T);
template <arithmetic T>
radian_T<T>& operator+=(radian_T<T>&, T);
template <arithmetic T>
radian_T<T>& operator-=(radian_T<T>&, T);
template <arithmetic T>
radian_T<T>& operator*=(radian_T<T>&, T);
template <arithmetic T>
radian_T<T>& operator/=(radian_T<T>&, T);

template <arithmetic T>
bool operator<(radian_T<T>, radian_T<T>);
template <arithmetic T>
bool operator>(radian_T<T>, radian_T<T>);
template <arithmetic T>
bool operator==(radian_T<T>, radian_T<T>);
template <arithmetic T>
bool operator!=(radian_T<T>, radian_T<T>);
template <arithmetic T>
bool operator>=(radian_T<T>, radian_T<T>);
template <arithmetic T>
bool operator<=(radian_T<T>, radian_T<T>);

template <arithmetic T>
T sin(radian_T<T>);
template <arithmetic T>
T cos(radian_T<T>);
template <arithmetic T>
T tan(radian_T<T>);
template <arithmetic T>
T asin(radian_T<T>);
template <arithmetic T>
T acos(radian_T<T>);
template <arithmetic T>
T atan(radian_T<T>);
template <arithmetic T>
T sinh(radian_T<T>);
template <arithmetic T>
T cosh(radian_T<T>);
template <arithmetic T>
T tanh(radian_T<T>);
template <arithmetic T>
T asinh(radian_T<T>);
template <arithmetic T>
T acosh(radian_T<T>);
template <arithmetic T>
T atanh(radian_T<T>);
#pragma endregion

#pragma region declaration_degree_T
template <arithmetic T>
struct degree_T
{
private:
	T angle;

public:
	degree_T() = default;
	constexpr explicit degree_T(T f);
	degree_T(const degree_T& d) = default;
	degree_T(degree_T&& d) = default;
	constexpr degree_T(radian_T<T> r);

	explicit operator T() const;

	template <arithmetic U>
	explicit operator degree_T<U>() const
	{
		return degree_T<U>{U(angle)};
	}

	degree_T& operator+=(degree_T d);
	degree_T& operator-=(degree_T d);
	degree_T& operator*=(degree_T d);
	degree_T& operator/=(degree_T d);

	degree_T operator-() const;

	degree_T& operator=(const degree_T& d) = default;
	degree_T& operator=(radian_T<T> r);

	using underlying_type = T;
};

template <arithmetic T>
degree_T<T> operator+(degree_T<T>, degree_T<T>);
template <arithmetic T>
degree_T<T> operator-(degree_T<T>, degree_T<T>);
template <arithmetic T>
degree_T<T> operator*(T, degree_T<T>);
template <arithmetic T>
degree_T<T> operator*(degree_T<T>, T);
template <arithmetic T>
degree_T<T> operator/(degree_T<T>, T);
template <arithmetic T>
degree_T<T>& operator+=(degree_T<T>&, T);
template <arithmetic T>
degree_T<T>& operator-=(degree_T<T>&, T);
template <arithmetic T>
degree_T<T>& operator*=(degree_T<T>&, T);
template <arithmetic T>
degree_T<T>& operator/=(degree_T<T>&, T);

template <arithmetic T>
bool operator<(degree_T<T>, degree_T<T>);
template <arithmetic T>
bool operator>(degree_T<T>, degree_T<T>);
template <arithmetic T>
bool operator==(degree_T<T>, degree_T<T>);
template <arithmetic T>
bool operator!=(degree_T<T>, degree_T<T>);
template <arithmetic T>
bool operator>=(degree_T<T>, degree_T<T>);
template <arithmetic T>
bool operator<=(degree_T<T>, degree_T<T>);

template <arithmetic T>
T sin(degree_T<T>);
template <arithmetic T>
T cos(degree_T<T>);
template <arithmetic T>
T tan(degree_T<T>);
template <arithmetic T>
T asin(degree_T<T>);
template <arithmetic T>
T acos(degree_T<T>);
template <arithmetic T>
T atan(degree_T<T>);
template <arithmetic T>
T sinh(degree_T<T>);
template <arithmetic T>
T cosh(degree_T<T>);
template <arithmetic T>
T tanh(degree_T<T>);
template <arithmetic T>
T asinh(degree_T<T>);
template <arithmetic T>
T acosh(degree_T<T>);
template <arithmetic T>
T atanh(degree_T<T>);
#pragma endregion

#pragma region declaration_pi_fraction_T
template <std::signed_integral T>
struct pi_fraction_T
{
	T numerator;
	T denominator;

	pi_fraction_T() = default;
	pi_fraction_T(T num, T den) : numerator{num}, denominator{den} {}
	pi_fraction_T(const pi_fraction_T& d) = default;
	pi_fraction_T(pi_fraction_T&& d) = default;

	pi_fraction_T& operator=(const pi_fraction_T&) = default;
	pi_fraction_T& operator=(pi_fraction_T&&) = default;

	template <arithmetic U>
	operator radian_T<U>() const
	{
		return radian_T<U>((U(numerator) * U(pi)) / U(denominator));
	}
	template <arithmetic U>
	operator degree_T<U>() const
	{
		return operator radian_T<U>();
	}

	pi_fraction_T operator-() const
	{
		return pi_fraction_T{-numerator, denominator};
	}

	using underlying_type = T;
};

template <arithmetic U, arithmetic T>
U sin(pi_fraction_T<T> d)
{
	if constexpr (d.numerator == T(0))
		return U(0);
	else if constexpr (d == static_pi_fraction_T<T, T(1), T(1)>{})
		return U(0);
	else if constexpr (d == static_pi_fraction_T<T, T(1), T(2)>{})
		return U(1);
	else if constexpr (d == static_pi_fraction_T<T, T(1), T(3)>{})
		return U(sqrt3_over_2);
	else if constexpr (d == static_pi_fraction_T<T, T(1), T(4)>{})
		return U(inv_sqrt2);
	else if constexpr (d == static_pi_fraction_T<T, T(1), T(6)>{})
		return U(0.5);

	else if constexpr (d == static_pi_fraction_T<T, T(-1), T(1)>{})
		return U(-0);
	else if constexpr (d == static_pi_fraction_T<T, T(-1), T(2)>{})
		return U(-1);
	else if constexpr (d == static_pi_fraction_T<T, T(-1), T(3)>{})
		return U(-sqrt3_over_2);
	else if constexpr (d == static_pi_fraction_T<T, T(-1), T(4)>{})
		return U(-inv_sqrt2);
	else if constexpr (d == static_pi_fraction_T<T, T(-1), T(6)>{})
		return U(-0.5);

	else if constexpr (d == static_pi_fraction_T<T, T(2), T(1)>{})
		return U(0);
	else if constexpr (d == static_pi_fraction_T<T, T(2), T(3)>{})
		return U(sqrt3_over_2);

	else if constexpr (d == static_pi_fraction_T<T, T(-2), T(1)>{})
		return U(-0);
	else if constexpr (d == static_pi_fraction_T<T, T(-2), T(3)>{})
		return U(-sqrt3_over_2);

	else if constexpr (d == static_pi_fraction_T<T, T(3), T(2)>{})
		return U(-1);
	else if constexpr (d == static_pi_fraction_T<T, T(3), T(4)>{})
		return U(inv_sqrt2);

	else if constexpr (d == static_pi_fraction_T<T, T(-3), T(2)>{})
		return U(1);
	else if constexpr (d == static_pi_fraction_T<T, T(-3), T(4)>{})
		return U(-inv_sqrt2);

	constexpr radian_T<U> r = d;
	return sin(r);
}
template <arithmetic U, arithmetic T>
U cos(pi_fraction_T<T> d)
{
	if constexpr (d.numerator == T(0))
		return U(1);
	else if constexpr (d.numerator == T(1) || d.numerator == T(-1))
	{
		if constexpr (d.denominator == T(1))
			return U(-1);
		else if constexpr (d.denominator == T(2))
			return U(0);
		else if constexpr (d.denominator == T(3))
			return U(0.5);
		else if constexpr (d.denominator == T(4))
			return U(inv_sqrt2);
		else if constexpr (d.denominator == T(6))
			return U(sqrt3_over_2);
	}
	else if constexpr (d.numerator == T(2) || d.numerator == T(-2))
	{
		if constexpr (d.denominator == T(1))
			return U(1);
		else if constexpr (d.denominator == T(3))
			return U(-0.5);
	}
	else if constexpr (d.numerator == T(3) || d.numerator == T(-3))
	{
		if constexpr (d.denominator == T(2))
			return U(0);
		else if constexpr (d.denominator == T(4))
			return U(-inv_sqrt2);
	}

	constexpr radian_T<U> r = d;
	return cos(r);
}
template <arithmetic U, arithmetic T>
U tan(pi_fraction_T<T> d)
{
	if constexpr (d.numerator == T(0))
		return U(0);
	else if constexpr (d.numerator == T(1))
	{
		if constexpr (d.denominator == T(1))
			return U(0);
		else if constexpr (d.denominator == T(2))
			return -std::numeric_limits<U>::infinity();
		else if constexpr (d.denominator == T(3))
			return U(sqrt3);
		else if constexpr (d.denominator == T(4))
			return U(1);
		else if constexpr (d.denominator == T(6))
			return U(inv_sqrt3);
	}
	else if constexpr (d.numerator == T(-1))
	{
		if constexpr (d.denominator == T(1))
			return U(-0);
		else if constexpr (d.denominator == T(2))
			return std::numeric_limits<U>::infinity();
		else if constexpr (d.denominator == T(3))
			return U(-sqrt3);
		else if constexpr (d.denominator == T(4))
			return U(-1);
		else if constexpr (d.denominator == T(6))
			return U(-inv_sqrt3);
	}
	else if constexpr (d.numerator == T(2))
	{
		if constexpr (d.denominator == T(1))
			return U(0);
		else if constexpr (d.denominator == T(3))
			return U(-sqrt3);
	}
	else if constexpr (d.numerator == T(-2))
	{
		if constexpr (d.denominator == T(1))
			return U(0);
		else if constexpr (d.denominator == T(3))
			return U(sqrt3);
	}
	else if constexpr (d.numerator == T(3))
	{
		if constexpr (d.denominator == T(2))
			return -std::numeric_limits<U>::infinity();
		else if constexpr (d.denominator == T(4))
			return U(-1);
	}
	else if constexpr (d.numerator == T(-3))
	{
		if constexpr (d.denominator == T(2))
			return std::numeric_limits<U>::infinity();
		else if constexpr (d.denominator == T(4))
			return U(1);
	}

	constexpr radian_T<U> r = d;
	return tan(r);
}
template <arithmetic U, arithmetic T>
U asin(pi_fraction_T<T> d)
{
	/// TODO: implement
	constexpr radian_T<U> r = d;
	return asin(r);
}
template <arithmetic U, arithmetic T>
U acos(pi_fraction_T<T> d)
{
	/// TODO: implement
	return acos(radian_T<U>(d));
}
template <arithmetic U, arithmetic T>
U atan(pi_fraction_T<T> d)
{
	/// TODO: implement
	return atan(radian_T<U>(d));
}
template <arithmetic U, arithmetic T>
U sinh(pi_fraction_T<T> d)
{
	/// TODO: implement
	constexpr radian_T<U> r = d;
	return sinh(r);
}
template <arithmetic U, arithmetic T>
U cosh(pi_fraction_T<T> d)
{
	/// TODO: implement
	constexpr radian_T<U> r = d;
	return cosh(r);
}
template <arithmetic U, arithmetic T>
U tanh(pi_fraction_T<T> d)
{
	/// TODO: implement
	constexpr radian_T<U> r = d;
	return tanh(r);
}
template <arithmetic U, arithmetic T>
U asinh(pi_fraction_T<T> d)
{
	/// TODO: implement
	constexpr radian_T<U> r = d;
	return asinh(r);
}
template <arithmetic U, arithmetic T>
U acosh(pi_fraction_T<T> d)
{
	/// TODO: implement
	constexpr radian_T<U> r = d;
	return acosh(r);
}
template <arithmetic U, arithmetic T>
U atanh(pi_fraction_T<T> d)
{
	/// TODO: implement
	constexpr radian_T<U> r = d;
	return atanh(r);
}
#pragma endregion

#pragma region declaration_static_pi_fraction_T
template <std::signed_integral T, T NumeratorT, T DenominatorT>
struct static_pi_fraction_T
{
	static_assert(DenominatorT != T(0), "the denominator can not be 0");

	template <T NumT, T DenT>
	static constexpr std::pair<T, T> resign()
	{
		if constexpr (DenT < T(0))
			return {-NumT, -DenT};
		else
			return {NumT, DenT};
	}

	template <T NumT, T DenT>
	static constexpr std::pair<T, T> simplify()
	{
		constexpr auto p = resign<NumT, DenT>();
		constexpr auto numt = p.first;
		constexpr auto dent = p.second;

		if constexpr (numt % dent == T(0))
			return {numt / dent, T(1)};
		else if constexpr (dent % numt == T(0))
		{
			if constexpr (numt > T(0))
				return {T(1), dent / numt};
			else
				return {T(-1), -dent / numt};
		}
		else
			return {numt, dent};
	}

	template <T NumT, T DenT>
	static constexpr std::pair<T, T> wrap()
	{
		constexpr auto p = simplify<NumT, DenT>();
		constexpr auto numt = p.first;
		constexpr auto dent = p.second;

		if constexpr (dent == T(1))
			return {numt % T(2), dent};
		else if constexpr (dent == T(2))
			return {numt % T(4), dent};
		else
			return {numt, dent};
	}

	static constexpr T numerator{wrap<NumeratorT, DenominatorT>().first};
	static constexpr T denominator{wrap<NumeratorT, DenominatorT>().second};

	template <arithmetic U>
	constexpr operator radian_T<U>() const
	{
		if constexpr (numerator == T(0))
			return radian_T<U>(U(0));
		else
			return radian_T<U>((U(numerator) * U(pi)) / U(denominator));
	}
	template <arithmetic U>
	constexpr operator degree_T<U>() const
	{
		if constexpr (numerator == T(0))
			return degree_T<U>(U(0));
		else
			return degree_T<U>(
			    radian_T<U>((U(numerator) * U(pi)) / U(denominator)));
	}
	constexpr operator pi_fraction_T<T>() const
	{
		return {numerator, denominator};
	}

	static_pi_fraction_T<T, -numerator, denominator> operator-() const
	{
		return {};
	}

	using underlying_type = T;
};

template <std::signed_integral T,
          T NumeratorTL,
          T DenominatorTL,
          T NumeratorTR,
          T DenominatorTR>
constexpr bool operator==(
    const static_pi_fraction_T<T, NumeratorTL, DenominatorTL>& lhs,
    const static_pi_fraction_T<T, NumeratorTR, DenominatorTR>& rhs)
{
	return lhs.numerator == rhs.numerator &&
	       lhs.denominator == rhs.denominator;
}
template <std::signed_integral T,
          T NumeratorTL,
          T DenominatorTL,
          T NumeratorTR,
          T DenominatorTR>
constexpr bool operator!=(
    const static_pi_fraction_T<T, NumeratorTL, DenominatorTL>& lhs,
    const static_pi_fraction_T<T, NumeratorTR, DenominatorTR>& rhs)
{
	return !(lhs == rhs);
}

template <arithmetic U, std::signed_integral T, T NumeratorT, T DenominatorT>
U sin(static_pi_fraction_T<T, NumeratorT, DenominatorT> d)
{
	if constexpr (d.numerator == T(0))
		return U(0);
	else if constexpr (d == static_pi_fraction_T<T, T(1), T(1)>{})
		return U(0);
	else if constexpr (d == static_pi_fraction_T<T, T(1), T(2)>{})
		return U(1);
	else if constexpr (d == static_pi_fraction_T<T, T(1), T(3)>{})
		return U(sqrt3_over_2);
	else if constexpr (d == static_pi_fraction_T<T, T(1), T(4)>{})
		return U(inv_sqrt2);
	else if constexpr (d == static_pi_fraction_T<T, T(1), T(6)>{})
		return U(0.5);

	else if constexpr (d == static_pi_fraction_T<T, T(-1), T(1)>{})
		return U(-0);
	else if constexpr (d == static_pi_fraction_T<T, T(-1), T(2)>{})
		return U(-1);
	else if constexpr (d == static_pi_fraction_T<T, T(-1), T(3)>{})
		return U(-sqrt3_over_2);
	else if constexpr (d == static_pi_fraction_T<T, T(-1), T(4)>{})
		return U(-inv_sqrt2);
	else if constexpr (d == static_pi_fraction_T<T, T(-1), T(6)>{})
		return U(-0.5);

	else if constexpr (d == static_pi_fraction_T<T, T(2), T(1)>{})
		return U(0);
	else if constexpr (d == static_pi_fraction_T<T, T(2), T(3)>{})
		return U(sqrt3_over_2);

	else if constexpr (d == static_pi_fraction_T<T, T(-2), T(1)>{})
		return U(-0);
	else if constexpr (d == static_pi_fraction_T<T, T(-2), T(3)>{})
		return U(-sqrt3_over_2);

	else if constexpr (d == static_pi_fraction_T<T, T(3), T(2)>{})
		return U(-1);
	else if constexpr (d == static_pi_fraction_T<T, T(3), T(4)>{})
		return U(inv_sqrt2);

	else if constexpr (d == static_pi_fraction_T<T, T(-3), T(2)>{})
		return U(1);
	else if constexpr (d == static_pi_fraction_T<T, T(-3), T(4)>{})
		return U(-inv_sqrt2);

	constexpr radian_T<U> r = d;
	return sin(r);
}
template <arithmetic U, std::signed_integral T, T NumeratorT, T DenominatorT>
U cos(static_pi_fraction_T<T, NumeratorT, DenominatorT> d)
{
	if constexpr (d.numerator == T(0))
		return U(1);
	else if constexpr (d.numerator == T(1) || d.numerator == T(-1))
	{
		if constexpr (d.denominator == T(1))
			return U(-1);
		else if constexpr (d.denominator == T(2))
			return U(0);
		else if constexpr (d.denominator == T(3))
			return U(0.5);
		else if constexpr (d.denominator == T(4))
			return U(inv_sqrt2);
		else if constexpr (d.denominator == T(6))
			return U(sqrt3_over_2);
	}
	else if constexpr (d.numerator == T(2) || d.numerator == T(-2))
	{
		if constexpr (d.denominator == T(1))
			return U(1);
		else if constexpr (d.denominator == T(3))
			return U(-0.5);
	}
	else if constexpr (d.numerator == T(3) || d.numerator == T(-3))
	{
		if constexpr (d.denominator == T(2))
			return U(0);
		else if constexpr (d.denominator == T(4))
			return U(-inv_sqrt2);
	}

	constexpr radian_T<U> r = d;
	return cos(r);
}
template <arithmetic U, arithmetic T, T NumeratorT, T DenominatorT>
U tan(static_pi_fraction_T<T, NumeratorT, DenominatorT> d)
{
	if constexpr (d.numerator == T(0))
		return U(0);
	else if constexpr (d.numerator == T(1))
	{
		if constexpr (d.denominator == T(1))
			return U(0);
		else if constexpr (d.denominator == T(3))
			return U(sqrt3);
		else if constexpr (d.denominator == T(4))
			return U(1);
		else if constexpr (d.denominator == T(6))
			return U(inv_sqrt3);
	}
	else if constexpr (d.numerator == T(-1))
	{
		if constexpr (d.denominator == T(1))
			return U(-0);
		else if constexpr (d.denominator == T(3))
			return U(-sqrt3);
		else if constexpr (d.denominator == T(4))
			return U(-1);
		else if constexpr (d.denominator == T(6))
			return U(-inv_sqrt3);
	}
	else if constexpr (d.numerator == T(2))
	{
		if constexpr (d.denominator == T(1))
			return U(0);
		else if constexpr (d.denominator == T(3))
			return U(-sqrt3);
	}
	else if constexpr (d.numerator == T(-2))
	{
		if constexpr (d.denominator == T(1))
			return U(0);
		else if constexpr (d.denominator == T(3))
			return U(sqrt3);
	}
	else if constexpr (d.numerator == T(3))
	{
		if constexpr (d.denominator == T(4))
			return U(-1);
	}
	else if constexpr (d.numerator == T(-3))
	{
		if constexpr (d.denominator == T(4))
			return U(1);
	}

	constexpr radian_T<U> r = d;
	return tan(r);
}
template <arithmetic U, arithmetic T, T NumeratorT, T DenominatorT>
U asin(static_pi_fraction_T<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	constexpr radian_T<U> r = d;
	return asin(r);
}
template <arithmetic U, arithmetic T, T NumeratorT, T DenominatorT>
U acos(static_pi_fraction_T<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	return acos(radian_T<U>(d));
}
template <arithmetic U, arithmetic T, T NumeratorT, T DenominatorT>
U atan(static_pi_fraction_T<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	return atan(radian_T<U>(d));
}
template <arithmetic U, arithmetic T, T NumeratorT, T DenominatorT>
U sinh(static_pi_fraction_T<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	constexpr radian_T<U> r = d;
	return sinh(r);
}
template <arithmetic U, arithmetic T, T NumeratorT, T DenominatorT>
U cosh(static_pi_fraction_T<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	constexpr radian_T<U> r = d;
	return cosh(r);
}
template <arithmetic U, arithmetic T, T NumeratorT, T DenominatorT>
U tanh(static_pi_fraction_T<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	constexpr radian_T<U> r = d;
	return tanh(r);
}
template <arithmetic U, arithmetic T, T NumeratorT, T DenominatorT>
U asinh(static_pi_fraction_T<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	constexpr radian_T<U> r = d;
	return asinh(r);
}
template <arithmetic U, arithmetic T, T NumeratorT, T DenominatorT>
U acosh(static_pi_fraction_T<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	constexpr radian_T<U> r = d;
	return acosh(r);
}
template <arithmetic U, arithmetic T, T NumeratorT, T DenominatorT>
U atanh(static_pi_fraction_T<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	constexpr radian_T<U> r = d;
	return atanh(r);
}
#pragma endregion

#pragma region declaration_vector2_t
template <arithmetic T>
struct vector2_T
{
	T x, y;

	vector2_T() = default;
	constexpr vector2_T(T x_, T y_) : x{x_}, y{y_} {}
	constexpr vector2_T(const direction2_T<T>& d);
	constexpr vector2_T(std::array<T, 2> a) : vector2_T(a[0], a[1]) {}

	template <arithmetic U>
	explicit operator vector2_T<U>() const
	{
		return vector2_T<U>{static_cast<U>(x), static_cast<U>(y)};
	}

	operator std::array<T, 2>() const { return {x, y}; }

	template <arithmetic U = T>
	std::array<U, 2> to_array() const
	{
		return {static_cast<U>(x), static_cast<U>(y)};
	}

	T norm() const { return std::sqrt(norm_squared()); }
	T norm_squared() const { return x * x + y * y; }
	vector2_T& normalize()
	{
		T n{norm()};
		x /= n;
		y /= n;
		return *this;
	}
	vector2_T get_normalized() const
	{
		vector2_T res{*this};
		res.normalize();
		return res;
	}
	constexpr vector2_T& clamp(vector2_T min, vector2_T max)
	{
		x = std::clamp(x, min.x, max.x);
		y = std::clamp(y, min.y, max.y);
		return *this;
	}
	constexpr vector2_T get_clamped(vector2_T min, vector2_T max) const
	{
		vector2_T res{*this};
		res.clamp(min, max);
		return res;
	}

	T& operator[](size_t i)
	{
		if (i == 0)
			return x;
		else
			return y;
	}
	const T& operator[](size_t i) const
	{
		if (i == 0)
			return x;
		else
			return y;
	}

	vector2_T operator+(vector2_T v) const;
	vector2_T operator-(vector2_T v) const;

	vector2_T& operator+=(vector2_T v);
	vector2_T& operator-=(vector2_T v);

	vector2_T operator-() const;

	bool operator==(vector2_T v) const;
	bool operator!=(vector2_T v) const;

	static constexpr vector2_T zero() { return { T(0), T(0) }; }
	static constexpr vector2_T one() { return { T(1), T(1) }; }

	using underlying_type = T;
};

template <arithmetic T>
vector2_T<T> operator*(T f, vector2_T<T> v);
template <arithmetic T>
vector2_T<T> operator*(vector2_T<T> v, T f);
template <arithmetic T>
vector2_T<T> operator/(vector2_T<T> v, T f);
template <arithmetic T>
vector2_T<T>& operator*=(vector2_T<T>& v, T f);
template <arithmetic T>
vector2_T<T>& operator/=(vector2_T<T>& v, T f);

template <arithmetic T>
constexpr T dot(vector2_T<T> lhs, vector2_T<T> rhs);
template <arithmetic T>
constexpr T cross(vector2_T<T> lhs, vector2_T<T> rhs);
//template <arithmetic T>
//constexpr vector2_T<T> clamp(vector2_T<T> v,
//                             vector2_T<T> min,
//                             vector2_T<T> max);
template <arithmetic T>
constexpr vector2_T<T> lerp(vector2_T<T> begin, vector2_T<T> end, T percent);
template <arithmetic T>
vector2_T<T> nlerp(vector2_T<T> begin, vector2_T<T> end, T percent);
template <arithmetic T>
vector2_T<T> slerp(vector2_T<T> begin, vector2_T<T> end, T percent);
template <arithmetic T>
T distance_between(vector2_T<T> v1, vector2_T<T> v2);
template <arithmetic T>
T distance_squared_between(vector2_T<T> v1, vector2_T<T> v2);
template <arithmetic T>
vector2_T<T> from_to(vector2_T<T> from, vector2_T<T> to);
template <arithmetic T>
radian_T<T> angle_between(vector2_T<T> v1, vector2_T<T> v2);
template <arithmetic T>
bool nearly_equal(vector2_T<T> v1, vector2_T<T> v2, T e = T(epsilon));
template <arithmetic T>
constexpr vector2_T<T> element_wise_min(vector2_T<T> v0, vector2_T<T> v1);
template <arithmetic T>
constexpr vector2_T<T> element_wise_max(vector2_T<T> v0, vector2_T<T> v1);
template <arithmetic T>
constexpr vector2_T<T> element_wise_abs(vector2_T<T> v);
template <arithmetic T>
constexpr vector2_T<T> element_wise_lerp(vector2_T<T> begin, vector2_T<T> end, vector2_T<T> percent);
//template <typename vector2_T>
//constexpr vector2_T zero();
//template <typename vector2_T>
//constexpr vector2_T one();
#pragma endregion

#pragma region declaration_vector3_t
template <arithmetic T>
struct vector3_T
{
	T x, y, z;

	vector3_T() = default;
	constexpr vector3_T(T x_, T y_, T z_) : x{x_}, y{y_}, z{z_} {}
	constexpr vector3_T(vector2_T<T> v, T z_) : vector3_T(v.x, v.y, z_) {}
	constexpr vector3_T(const direction3_T<T>& d);
	constexpr vector3_T(std::array<T, 3> a) : vector3_T(a[0], a[1], a[2]) {}

	template <arithmetic U>
	explicit operator vector3_T<U>() const
	{
		return vector3_T<U>{static_cast<U>(x), static_cast<U>(y),
		                    static_cast<U>(z)};
	}

	operator std::array<T, 3>() const { return {x, y, z}; }

	template <arithmetic U = T>
	std::array<U, 3> to_array() const
	{
		return {static_cast<U>(x), static_cast<U>(y), static_cast<U>(z)};
	}

	vector2_T<T> xy() const;

	T norm() const { return std::sqrt(norm_squared()); }
	T norm_squared() const { return x * x + y * y + z * z; }
	vector3_T& normalize()
	{
		T n = norm();
		x /= n;
		y /= n;
		z /= n;
		return *this;
	}
	vector3_T get_normalized() const
	{
		vector3_T res{*this};
		res.normalize();
		return res;
	}
	constexpr vector3_T& clamp(vector3_T min, vector3_T max)
	{
		x = std::clamp(x, min.x, max.x);
		y = std::clamp(y, min.y, max.y);
		z = std::clamp(z, min.z, max.z);
		return *this;
	}
	constexpr vector3_T get_clamped(vector3_T min, vector3_T max) const
	{
		vector3_T res{*this};
		res.clamp(min, max);
		return res;
	}

	vector3_T operator+(vector3_T v) const;
	vector3_T operator-(vector3_T v) const;

	vector3_T& operator+=(vector3_T v);
	vector3_T& operator-=(vector3_T v);

	vector3_T operator-() const;

	bool operator==(vector3_T v) const;
	bool operator!=(vector3_T v) const;

	static constexpr vector3_T zero() { return { T(0), T(0), T(0) }; }
	static constexpr vector3_T one() { return { T(1), T(1), T(1) }; }

	using underlying_type = T;
};

template <arithmetic T>
vector3_T<T> operator*(T f, vector3_T<T> v);
template <arithmetic T>
vector3_T<T> operator*(vector3_T<T> v, T f);
template <arithmetic T>
vector3_T<T> operator/(vector3_T<T> v, T f);
template <arithmetic T>
vector3_T<T>& operator*=(vector3_T<T>& v, T f);
template <arithmetic T>
vector3_T<T>& operator/=(vector3_T<T>& v, T f);

template <arithmetic T>
T norm(vector3_T<T> v);
template <arithmetic T>
T norm_squared(vector3_T<T> v);
template <arithmetic T>
vector3_T<T> normalize(vector3_T<T> v);
//template <arithmetic T>
//constexpr vector3_T<T> clamp(vector3_T<T> v,
//                             vector3_T<T> min,
//                             vector3_T<T> max);
template <arithmetic T>
constexpr T dot(vector3_T<T> lhs, vector3_T<T> rhs);
template <arithmetic T>
constexpr vector3_T<T> cross(vector3_T<T> lhs, vector3_T<T> rhs);
template <arithmetic T>
vector3_T<T> lerp(vector3_T<T> begin, vector3_T<T> end, T percent);
template <arithmetic T>
vector3_T<T> nlerp(vector3_T<T> begin, vector3_T<T> end, T percent);
template <arithmetic T>
T distance_between(vector3_T<T> v1, vector3_T<T> v2);
template <arithmetic T>
T distance_squared_between(vector3_T<T> v1, vector3_T<T> v2);
template <arithmetic T>
vector3_T<T> from_to(vector3_T<T> from, vector3_T<T> to);
template <arithmetic T>
radian_T<T> angle_between(vector3_T<T> v1, vector3_T<T> v2);
template <arithmetic T>
bool nearly_equal(vector3_T<T> v1, vector3_T<T> v2, T e = epsilon);
template <arithmetic T>
constexpr vector3_T<T> element_wise_min(vector3_T<T> v0, vector3_T<T> v1);
template <arithmetic T>
constexpr vector3_T<T> element_wise_max(vector3_T<T> v0, vector3_T<T> v1);
template <arithmetic T>
constexpr vector3_T<T> element_wise_abs(vector3_T<T> v);
template <arithmetic T>
constexpr vector3_T<T> element_wise_lerp(vector3_T<T> begin, vector3_T<T> end, vector3_T<T> percent);
template <arithmetic T>
radian_T<T> angle_around_axis(vector3_T<T> v0,
                              vector3_T<T> v1,
                              direction3_T<T> axis);
#pragma endregion

#pragma region declaration_vector4_t
template <arithmetic T>
struct vector4_T
{
	T x, y, z, w;

	vector4_T() = default;
	constexpr vector4_T(T x_, T y_, T z_, T w_) : x{x_}, y{y_}, z{z_}, w{w_} {}
	constexpr vector4_T(vector2_T<T> v, T z_, T w_) : vector4_T(v.x, v.y, z_, w_) {}
	constexpr vector4_T(vector3_T<T> v, T w_) : vector4_T(v.x, v.y, v.z, w_) {}
	constexpr vector4_T(std::array<T, 4> a) : vector4_T(a[0], a[1], a[2], a[3]) {}

	template <arithmetic U>
	explicit operator vector4_T<U>() const
	{
		return vector4_T<U>{static_cast<U>(x), static_cast<U>(y),
		                    static_cast<U>(z), static_cast<U>(w)};
	}

	operator std::array<T, 4>() const { return {x, y, z, w}; }

	template <arithmetic U = T>
	std::array<U, 4> to_array() const
	{
		return {static_cast<U>(x), static_cast<U>(y), static_cast<U>(z),
		        static_cast<U>(w)};
	}

	vector2_T<T> xy() const;
	vector3_T<T> xyz() const;

	T norm() const { return std::sqrt(norm_squared()); }
	T norm_squared() const { return x * x + y * y + z * z + w * w; }
	vector4_T& normalize()
	{
		T n = norm();
		x /= n;
		y /= n;
		z /= n;
		w /= n;
		return *this;
	}
	vector4_T get_normalized() const
	{
		vector4_T res{*this};
		res.normalize();
		return res;
	}
	constexpr vector4_T& clamp(vector4_T min, vector4_T max)
	{
		x = std::clamp(x, min.x, max.x);
		y = std::clamp(y, min.y, max.y);
		z = std::clamp(z, min.z, max.z);
		w = std::clamp(w, min.w, max.w);
		return *this;
	}
	constexpr vector4_T get_clamped(vector4_T min, vector4_T max) const
	{
		vector4_T res{*this};
		res.clamp(min, max);
		return res;
	}

	vector4_T operator+(vector4_T v) const;
	vector4_T operator-(vector4_T v) const;

	vector4_T& operator+=(vector4_T v);
	vector4_T& operator-=(vector4_T v);

	vector4_T operator-() const;

	bool operator==(vector4_T v) const;
	bool operator!=(vector4_T v) const;

	static constexpr vector4_T zero() { return { T(0), T(0), T(0), T(0) }; }
	static constexpr vector4_T one() { return { T(1), T(1), T(1), T(1) }; }

	using underlying_type = T;
};

template <arithmetic T>
vector4_T<T> operator*(T f, vector4_T<T> v);
template <arithmetic T>
vector4_T<T> operator*(vector4_T<T> v, T f);
template <arithmetic T>
vector4_T<T> operator/(vector4_T<T> v, T f);
template <arithmetic T>
vector4_T<T>& operator*=(vector4_T<T>& v, T f);
template <arithmetic T>
vector4_T<T>& operator/=(vector4_T<T>& v, T f);

template <arithmetic T>
T norm(vector4_T<T> v);
template <arithmetic T>
T norm_squared(vector4_T<T> v);
template <arithmetic T>
vector4_T<T> normalize(vector4_T<T> v);
//template <arithmetic T>
//constexpr vector4_T<T> clamp(vector4_T<T> v,
//                             vector4_T<T> min,
//                             vector4_T<T> max);
template <arithmetic T>
T dot(vector4_T<T> lhs, vector4_T<T> rhs);
template <arithmetic T>
vector4_T<T> lerp(vector4_T<T> begin, vector4_T<T> end, T percent);
template <arithmetic T>
vector4_T<T> nlerp(vector4_T<T> begin, vector4_T<T> end, T percent);
template <arithmetic T>
T distance_between(vector4_T<T> v1, vector4_T<T> v2);
template <arithmetic T>
T distance_squared_between(vector4_T<T> v1, vector4_T<T> v2);
template <arithmetic T>
vector4_T<T> from_to(vector4_T<T> from, vector4_T<T> to);
template <arithmetic T>
bool nearly_equal(vector4_T<T> v1, vector4_T<T> v2, T e = epsilon);
template <arithmetic T>
constexpr vector4_T<T> element_wise_min(vector4_T<T> v0, vector4_T<T> v1);
template <arithmetic T>
constexpr vector4_T<T> element_wise_max(vector4_T<T> v0, vector4_T<T> v1);
template <arithmetic T>
constexpr vector4_T<T> element_wise_abs(vector4_T<T> v);
template <arithmetic T>
constexpr vector4_T<T> element_wise_lerp(vector4_T<T> begin, vector4_T<T> end, vector4_T<T> percent);
#pragma endregion

#pragma region declaration_normalized
template <normalizable T>
class normalized
{
	T vector;

public:
	normalized() = default;
	normalized(const normalized&) = default;
	normalized(normalized&&) = default;
	explicit normalized(const T& t);
	template <typename... Args>
	explicit normalized(Args... args) : normalized{T{args...}}
	{
	}

	static normalized already_normalized(const T& t);
	static normalized already_normalized(T&& t) noexcept;

	const T& operator*() const;
	const T* operator->() const;
	operator const T&() const;

	normalized operator-() const;

	bool operator==(const normalized& n) const;
	bool operator!=(const normalized& n) const;
	bool operator==(const T& v) const;
	bool operator!=(const T& v) const;

	normalized& operator=(const normalized&) = default;
	normalized& operator=(normalized&&) = default;

	using underlying_type = T;
};
#pragma endregion

#pragma region declaration_matrix2_t
template <arithmetic T>
struct matrix2_T
{
	template <arithmetic U>
	friend class matrix2_T;
	template <arithmetic U>
	friend class matrix3_T;
	template <arithmetic U>
	friend class matrix4_T;

private:
	T m00, m01, m10, m11;

	matrix2_T(T m00, T m10, T m01, T m11);

public:
	matrix2_T() = default;
	explicit matrix2_T(const std::array<T, 4>& a);
	matrix2_T(const matrix2_T&) = default;
	matrix2_T(matrix2_T&&) = default;

	matrix2_T& operator=(const matrix2_T&) = default;
	matrix2_T& operator=(matrix2_T&&) = default;

	template <arithmetic U = T>
	std::array<U, 4> to_array() const
	{
		return {static_cast<U>(m00), static_cast<U>(m10), static_cast<U>(m01),
		        static_cast<U>(m11)};
	}

	template <arithmetic U>
	explicit operator matrix2_T<U>() const
	{
		return {U(m00), U(m10), U(m01), U(m11)};
	}

	static matrix2_T zero();
	static matrix2_T identity();
	static matrix2_T rows(vector2_T<T> row0, vector2_T<T> row1);
	static matrix2_T columns(vector2_T<T> column0, vector2_T<T> column1);
	static matrix2_T rotation(radian_T<T> angle);
	static matrix2_T rotation(normalized<complex_T<T>> angle);
	matrix2_T& transpose();
	matrix2_T get_transposed() const;
	T determinant() const;

	template <bool IsConstT>
	struct proxy
	{
	protected:
		using Type = std::conditional_t<IsConstT, const T, T>;

		Type& x;
		Type& y;

		constexpr Type& vector(int i)
		{
			return std::array<std::reference_wrapper<Type>, 2>{x, y}[i];
		}
		constexpr const Type& vector(int i) const
		{
			return std::array<std::reference_wrapper<const Type>, 2>{x, y}[i];
		}

	public:
		constexpr proxy(Type& e0, Type& e1) : x{e0}, y{e1} {}
		constexpr proxy(const proxy&) = default;
		constexpr proxy(proxy&&) = default;
		constexpr proxy& operator=(const proxy&) = default;
		constexpr proxy& operator=(proxy&&) = default;
		constexpr proxy& operator=(vector2_T<T> v)
		{
			x = v.x;
			y = v.y;
			return *this;
		}

		constexpr operator vector2_T<T>() const { return vector2_T<T>{x, y}; }
	};
	template <bool IsConstT>
	struct column_proxy : proxy<IsConstT>
	{
		using proxy<IsConstT>::proxy;
		using proxy<IsConstT>::operator=;
		using proxy<IsConstT>::x;
		using proxy<IsConstT>::y;
		constexpr proxy<IsConstT>::Type& row(int i)
		{
			return proxy<IsConstT>::vector(i);
		}
		constexpr const proxy<IsConstT>::Type& row(int i) const
		{
			return proxy<IsConstT>::vector(i);
		}
	};
	template <bool IsConstT>
	struct row_proxy : proxy<IsConstT>
	{
		using proxy<IsConstT>::proxy;
		using proxy<IsConstT>::operator=;
		using proxy<IsConstT>::x;
		using proxy<IsConstT>::y;
		constexpr proxy<IsConstT>::Type& column(int i)
		{
			return proxy<IsConstT>::vector(i);
		}
		constexpr const proxy<IsConstT>::Type& column(int i) const
		{
			return proxy<IsConstT>::vector(i);
		}
	};
	template <bool IsConstT>
	struct diag_proxy : proxy<IsConstT>
	{
		using proxy<IsConstT>::proxy;
		using proxy<IsConstT>::operator=;
		using proxy<IsConstT>::x;
		using proxy<IsConstT>::y;
	};

	constexpr column_proxy<false> column(int i)
	{
		return std::array{
		    column_proxy<false>{m00, m01},
		    column_proxy<false>{m10, m11},
		}[i];
	}
	constexpr column_proxy<true> column(int i) const
	{
		return std::array{
		    column_proxy<true>{m00, m01},
		    column_proxy<true>{m10, m11},
		}[i];
	}
	constexpr row_proxy<false> row(int i)
	{
		return std::array{
		    row_proxy<false>{m00, m10},
		    row_proxy<false>{m01, m11},
		}[i];
	}
	constexpr row_proxy<true> row(int i) const
	{
		return std::array{
		    row_proxy<true>{m00, m10},
		    row_proxy<true>{m01, m11},
		}[i];
	}
	constexpr diag_proxy<false> diag() { return diag_proxy<false>{m00, m11}; }
	constexpr diag_proxy<true> diag() const
	{
		return diag_proxy<true>{m00, m11};
	}
	constexpr T& diag(int i)
	{
		return std::array{
		    std::ref(m00),
		    std::ref(m11),
		}[i];
	}
	constexpr const T& diag(int i) const
	{
		return std::array{
		    std::cref(m00),
		    std::cref(m11),
		}[i];
	}

	vector2_T<T> operator*(vector2_T<T> v) const;
	matrix2_T operator*(matrix2_T m) const;

	using underlying_type = T;
};

template <arithmetic T>
matrix2_T<T> transpose(matrix2_T<T> m);
#pragma endregion

#pragma region declaration_matrix3_t
template <arithmetic T>
struct matrix3_T
{
	template <arithmetic U>
	friend class matrix2_T;
	template <arithmetic U>
	friend class matrix3_T;
	template <arithmetic U>
	friend class matrix4_T;

private:
	T m00, m01, m02, m10, m11, m12, m20, m21, m22;

	matrix3_T(T m00, T m10, T m20, T m01, T m11, T m21, T m02, T m12, T m22);

public:
	matrix3_T() = default;
	explicit matrix3_T(matrix2_T<T> m);
	matrix3_T(const std::array<T, 9>& a);
	matrix3_T(const matrix3_T&) = default;
	matrix3_T(matrix3_T&&) = default;

	matrix3_T& operator=(const matrix3_T&) = default;
	matrix3_T& operator=(matrix3_T&&) = default;

	template <arithmetic U>
	explicit operator matrix3_T<U>() const
	{
		return {U(m00), U(m10), U(m20), U(m01), U(m11),
		        U(m21), U(m02), U(m12), U(m22)};
	}

	operator std::array<T, 9>() const
	{
		return {
		    m00, m10, m20,  //
		    m01, m11, m21,  //
		    m02, m12, m22   //
		};
	}

	template <arithmetic U = T>
	std::array<U, 9> to_array() const
	{
		return {static_cast<U>(m00), static_cast<U>(m10), static_cast<U>(m20),
		        static_cast<U>(m01), static_cast<U>(m11), static_cast<U>(m21),
		        static_cast<U>(m02), static_cast<U>(m12), static_cast<U>(m22)};
	}

	static matrix3_T zero();
	static matrix3_T identity();
	static matrix3_T rows(vector3_T<T> row0,
	                      vector3_T<T> row1,
	                      vector3_T<T> row2);
	static matrix3_T columns(vector3_T<T> column0,
	                         vector3_T<T> column1,
	                         vector3_T<T> column2);

private:
	static matrix3_T rotation(direction3_T<T> axis, T sinAngle, T cosAngle);

public:
	static matrix3_T rotation(direction3_T<T> axis, radian_T<T> angle);
	static matrix3_T rotation(direction3_T<T> axis,
	                          normalized<complex_T<T>> angle);
	template <std::signed_integral U, U NumeratorT, U DenominatorT>
	static matrix3_T rotation(
	    direction3_T<T> axis,
	    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle);
	static matrix3_T rotation(euler_angles_T<T> r);
	static matrix3_T rotation(quaternion_T<T> q);

private:
	static matrix3_T rotation_around_x(T sinAngle, T cosAngle);
	static matrix3_T rotation_around_y(T sinAngle, T cosAngle);
	static matrix3_T rotation_around_z(T sinAngle, T cosAngle);

public:
	static matrix3_T rotation_around_x(radian_T<T> angle);
	static matrix3_T rotation_around_y(radian_T<T> angle);
	static matrix3_T rotation_around_z(radian_T<T> angle);
	static matrix3_T rotation_around_x(normalized<complex_T<T>> angle);
	static matrix3_T rotation_around_y(normalized<complex_T<T>> angle);
	static matrix3_T rotation_around_z(normalized<complex_T<T>> angle);
	template <std::signed_integral U, U NumeratorT, U DenominatorT>
	static matrix3_T rotation_around_x(
	    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle);
	template <std::signed_integral U, U NumeratorT, U DenominatorT>
	static matrix3_T rotation_around_y(
	    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle);
	template <std::signed_integral U, U NumeratorT, U DenominatorT>
	static matrix3_T rotation_around_z(
	    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle);
	static matrix3_T scale(vector3_T<T> t);
	static matrix3_T scale(T x, T y, T z);
	static matrix3_T scale(T s);
	matrix3_T& inverse();
	matrix3_T get_inversed() const;
	matrix3_T& transpose();
	matrix3_T get_transposed() const;
	T determinant() const;
	bool is_orthogonal() const;

	euler_angles_T<T> to_euler_angles() const;

	template <bool IsConstT>
	struct proxy
	{
	protected:
		using Type = std::conditional_t<IsConstT, const T, T>;

		Type& x;
		Type& y;
		Type& z;

		constexpr Type& vector(int i)
		{
			return std::array<std::reference_wrapper<Type>, 3>{x, y, z}[i];
		}
		constexpr const Type& vector(int i) const
		{
			return std::array<std::reference_wrapper<const Type>, 3>{x, y,
			                                                         z}[i];
		}

	public:
		constexpr proxy(Type& e0, Type& e1, Type& e2) : x{e0}, y{e1}, z{e2} {}
		constexpr proxy(const proxy&) = default;
		constexpr proxy(proxy&&) = default;
		constexpr proxy& operator=(const proxy&) = default;
		constexpr proxy& operator=(proxy&&) = default;
		constexpr proxy& operator=(vector3_T<T> v)
		{
			x = v.x;
			y = v.y;
			z = v.z;
			return *this;
		}

		constexpr vector2_T<T> xy() const { return vector2_T<T>{x, y}; }
		constexpr operator vector3_T<T>() const
		{
			return vector3_T<T>{x, y, z};
		}
	};
	template <bool IsConstT>
	struct column_proxy : proxy<IsConstT>
	{
		using proxy<IsConstT>::proxy;
		using proxy<IsConstT>::operator=;
		using proxy<IsConstT>::x;
		using proxy<IsConstT>::y;
		using proxy<IsConstT>::z;
		constexpr proxy<IsConstT>::Type& row(int i)
		{
			return proxy<IsConstT>::vector(i);
		}
		constexpr const proxy<IsConstT>::Type& row(int i) const
		{
			return proxy<IsConstT>::vector(i);
		}
	};
	template <bool IsConstT>
	struct row_proxy : proxy<IsConstT>
	{
		using proxy<IsConstT>::proxy;
		using proxy<IsConstT>::operator=;
		using proxy<IsConstT>::x;
		using proxy<IsConstT>::y;
		using proxy<IsConstT>::z;
		constexpr proxy<IsConstT>::Type& column(int i)
		{
			return proxy<IsConstT>::vector(i);
		}
		constexpr const proxy<IsConstT>::Type& column(int i) const
		{
			return proxy<IsConstT>::vector(i);
		}
	};
	template <bool IsConstT>
	struct diag_proxy : proxy<IsConstT>
	{
		using proxy<IsConstT>::proxy;
		using proxy<IsConstT>::operator=;
		using proxy<IsConstT>::x;
		using proxy<IsConstT>::y;
		using proxy<IsConstT>::z;
	};

	constexpr column_proxy<false> column(int i)
	{
		return std::array{
		    column_proxy<false>{m00, m01, m02},
		    column_proxy<false>{m10, m11, m12},
		    column_proxy<false>{m20, m21, m22},
		}[i];
	}
	constexpr column_proxy<true> column(int i) const
	{
		return std::array{
		    column_proxy<true>{m00, m01, m02},
		    column_proxy<true>{m10, m11, m12},
		    column_proxy<true>{m20, m21, m22},
		}[i];
	}
	constexpr row_proxy<false> row(int i)
	{
		return std::array{
		    row_proxy<false>{m00, m10, m20},
		    row_proxy<false>{m01, m11, m21},
		    row_proxy<false>{m02, m12, m22},
		}[i];
	}
	constexpr row_proxy<true> row(int i) const
	{
		return std::array{
		    row_proxy<true>{m00, m10, m20},
		    row_proxy<true>{m01, m11, m21},
		    row_proxy<true>{m02, m12, m22},
		}[i];
	}
	constexpr diag_proxy<false> diag()
	{
		return diag_proxy<false>{m00, m11, m22};
	}
	constexpr diag_proxy<true> diag() const
	{
		return diag_proxy<true>{m00, m11, m22};
	}
	constexpr T& diag(int i)
	{
		return std::array{
		    std::ref(m00),
		    std::ref(m11),
		    std::ref(m22),
		}[i];
	}
	constexpr const T& diag(int i) const
	{
		return std::array{
		    std::cref(m00),
		    std::cref(m11),
		    std::cref(m22),
		}[i];
	}

	matrix3_T operator*(T f) const;
	vector3_T<T> operator*(vector3_T<T> v) const;
	matrix3_T operator*(const matrix3_T& m) const;

	bool operator==(const matrix3_T& m) const
	{
		return m00 == m.m00 && m01 == m.m01 && m02 == m.m02 && m10 == m.m10 &&
		       m11 == m.m11 && m12 == m.m12 && m20 == m.m20 && m21 == m.m21 &&
		       m22 == m.m22;
	}

	using underlying_type = T;
};
#pragma endregion

#pragma region declaration_matrix4_t
template <arithmetic T>
struct matrix4_T
{
	template <arithmetic U>
	friend class matrix2_T;
	template <arithmetic U>
	friend class matrix3_T;
	template <arithmetic U>
	friend class matrix4_T;

private:
	T m00, m01, m02, m03,    //
	    m10, m11, m12, m13,  //
	    m20, m21, m22, m23,  //
	    m30, m31, m32, m33;

	matrix4_T(T m00,
	          T m10,
	          T m20,
	          T m30,
	          T m01,
	          T m11,
	          T m21,
	          T m31,
	          T m02,
	          T m12,
	          T m22,
	          T m32,
	          T m03,
	          T m13,
	          T m23,
	          T m33);

public:
	matrix4_T() = default;
	explicit matrix4_T(matrix2_T<T> m);
	explicit matrix4_T(const matrix3_T<T>& m);
	explicit matrix4_T(const perspective_T<T>& p);
	explicit matrix4_T(const transform3_T<T>& t);
	explicit matrix4_T(const uniform_transform3_T<T>& t);
	explicit matrix4_T(const unscaled_transform3_T<T>& t);
	explicit matrix4_T(const std::array<T, 16>& a);
	matrix4_T(const matrix4_T&) = default;
	matrix4_T(matrix4_T&&) = default;

	matrix4_T& operator=(const matrix4_T&) = default;
	matrix4_T& operator=(matrix4_T&&) = default;

	template <arithmetic U>
	explicit operator matrix4_T<U>() const
	{
		return {U(m00), U(m10), U(m20), U(m30), U(m01), U(m11),
		        U(m21), U(m31), U(m02), U(m12), U(m22), U(m32),
		        U(m03), U(m13), U(m23), U(m33)};
	}

	operator std::array<T, 16>() const
	{
		return {
		    m00, m10, m20, m30,  //
		    m01, m11, m21, m31,  //
		    m02, m12, m22, m32,  //
		    m03, m13, m23, m33   //
		};
	}

	template <arithmetic U = T>
	std::array<U, 16> to_array() const
	{
		return {
		    static_cast<U>(m00),  //
		    static_cast<U>(m10),  //
		    static_cast<U>(m20),  //
		    static_cast<U>(m30),  //
		    static_cast<U>(m01),  //
		    static_cast<U>(m11),  //
		    static_cast<U>(m21),  //
		    static_cast<U>(m31),  //
		    static_cast<U>(m02),  //
		    static_cast<U>(m12),  //
		    static_cast<U>(m22),  //
		    static_cast<U>(m32),  //
		    static_cast<U>(m03),  //
		    static_cast<U>(m13),  //
		    static_cast<U>(m23),  //
		    static_cast<U>(m33)   //
		};
	}

	static matrix4_T zero();
	static matrix4_T identity();
	static matrix4_T rows(vector4_T<T> row0,
	                      vector4_T<T> row1,
	                      vector4_T<T> row2,
	                      vector4_T<T> row3);
	static matrix4_T columns(vector4_T<T> column0,
	                         vector4_T<T> column1,
	                         vector4_T<T> column2,
	                         vector4_T<T> column3);
	static matrix4_T rotation(direction3_T<T> axis, radian_T<T> angle);
	static matrix4_T rotation(direction3_T<T> axis,
	                          normalized<complex_T<T>> angle);
	template <std::signed_integral U, U NumeratorT, U DenominatorT>
	static matrix4_T rotation(
	    direction3_T<T> axis,
	    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle);
	static matrix4_T rotation(euler_angles_T<T> r);
	static matrix4_T rotation(quaternion_T<T> q);
	static matrix4_T translation(vector3_T<T> t);
	static matrix4_T translation(T x, T y, T z);
	static matrix4_T scale(vector3_T<T> t);
	static matrix4_T scale(T x, T y, T z);
	static matrix4_T scale(T s);
	static matrix4_T look_at(
	    vector3_T<T> from,
	    vector3_T<T> to,
	    direction3_T<T> up = direction3_T<T>::already_normalized(
	        vector3_T<T>(T(0), T(1), T(0))));
	static matrix4_T orthographic(T left,
	                              T right,
	                              T top,
	                              T bottom,
	                              T near,
	                              T far);
	static matrix4_T rotation_around_x(radian_T<T> angle);
	static matrix4_T rotation_around_y(radian_T<T> angle);
	static matrix4_T rotation_around_z(radian_T<T> angle);
	static matrix4_T rotation_around_x(normalized<complex_T<T>> angle);
	static matrix4_T rotation_around_y(normalized<complex_T<T>> angle);
	static matrix4_T rotation_around_z(normalized<complex_T<T>> angle);
	template <std::signed_integral T, T NumeratorT, T DenominatorT>
	static matrix4_T rotation_around_x(
	    static_pi_fraction_T<T, NumeratorT, DenominatorT> angle);
	template <std::signed_integral T, T NumeratorT, T DenominatorT>
	static matrix4_T rotation_around_y(
	    static_pi_fraction_T<T, NumeratorT, DenominatorT> angle);
	template <std::signed_integral T, T NumeratorT, T DenominatorT>
	static matrix4_T rotation_around_z(
	    static_pi_fraction_T<T, NumeratorT, DenominatorT> angle);
	bool is_orthogonal() const;
	bool is_homogenous() const;
	matrix4_T& inverse();
	matrix4_T get_inversed() const;
	matrix4_T& transpose();
	matrix4_T get_transposed() const;
	T determinant() const;

	matrix3_T<T> extractRotationMatrix() const;

	vector3_T<T> transform_position(vector3_T<T> pos) const;
	vector3_T<T> transform_direction(vector3_T<T> dir) const;
	direction3_T<T> transform_direction(direction3_T<T> dir) const;
	vector3_T<T> transform_position_perspective(vector3_T<T> pos) const;

	template <bool IsConstT>
	struct proxy
	{
	protected:
		using Type = std::conditional_t<IsConstT, const T, T>;

		Type& x;
		Type& y;
		Type& z;
		Type& w;

		constexpr Type& vector(int i)
		{
			return std::array<std::reference_wrapper<Type>, 4>{x, y, z, w}[i];
		}
		constexpr const Type& vector(int i) const
		{
			return std::array<std::reference_wrapper<const Type>, 4>{x, y, z,
			                                                         w}[i];
		}

	public:
		constexpr proxy(Type& e0, Type& e1, Type& e2, Type& e3)
		    : x{e0}, y{e1}, z{e2}, w{e3}
		{
		}
		constexpr proxy(const proxy&) = default;
		constexpr proxy(proxy&&) = default;
		constexpr proxy& operator=(const proxy&) = default;
		constexpr proxy& operator=(proxy&&) = default;
		constexpr proxy& operator=(vector4_T<T> v)
		{
			x = v.x;
			y = v.y;
			z = v.z;
			w = v.w;
			return *this;
		}

		constexpr vector2_T<T> xy() const { return vector2_T<T>{x, y}; }
		constexpr vector3_T<T> xyz() const { return vector3_T<T>{x, y, z}; }
		constexpr operator vector4_T<T>() const
		{
			return vector4_T<T>{x, y, z, w};
		}
	};
	template <bool IsConstT>
	struct column_proxy : proxy<IsConstT>
	{
		using proxy<IsConstT>::proxy;
		using proxy<IsConstT>::operator=;
		using proxy<IsConstT>::x;
		using proxy<IsConstT>::y;
		using proxy<IsConstT>::z;
		using proxy<IsConstT>::w;
		constexpr proxy<IsConstT>::Type& row(int i)
		{
			return proxy<IsConstT>::vector(i);
		}
		constexpr const proxy<IsConstT>::Type& row(int i) const
		{
			return proxy<IsConstT>::vector(i);
		}
	};
	template <bool IsConstT>
	struct row_proxy : proxy<IsConstT>
	{
		using proxy<IsConstT>::proxy;
		using proxy<IsConstT>::operator=;
		using proxy<IsConstT>::x;
		using proxy<IsConstT>::y;
		using proxy<IsConstT>::z;
		using proxy<IsConstT>::w;
		constexpr proxy<IsConstT>::Type& column(int i)
		{
			return proxy<IsConstT>::vector(i);
		}
		constexpr const proxy<IsConstT>::Type& column(int i) const
		{
			return proxy<IsConstT>::vector(i);
		}
	};
	template <bool IsConstT>
	struct diag_proxy : proxy<IsConstT>
	{
		using proxy<IsConstT>::proxy;
		using proxy<IsConstT>::operator=;
		using proxy<IsConstT>::x;
		using proxy<IsConstT>::y;
		using proxy<IsConstT>::z;
		using proxy<IsConstT>::w;
	};

	constexpr column_proxy<false> column(int i)
	{
		return std::array{
		    column_proxy<false>{m00, m01, m02, m03},
		    column_proxy<false>{m10, m11, m12, m13},
		    column_proxy<false>{m20, m21, m22, m23},
		    column_proxy<false>{m30, m31, m32, m33},
		}[i];
	}
	constexpr column_proxy<true> column(int i) const
	{
		return std::array{
		    column_proxy<true>{m00, m01, m02, m03},
		    column_proxy<true>{m10, m11, m12, m13},
		    column_proxy<true>{m20, m21, m22, m23},
		    column_proxy<true>{m30, m31, m32, m33},
		}[i];
	}
	constexpr row_proxy<false> row(int i)
	{
		return std::array{
		    row_proxy<false>{m00, m10, m20, m30},
		    row_proxy<false>{m01, m11, m21, m31},
		    row_proxy<false>{m02, m12, m22, m32},
		    row_proxy<false>{m03, m13, m23, m33},
		}[i];
	}
	constexpr row_proxy<true> row(int i) const
	{
		return std::array{
		    row_proxy<true>{m00, m10, m20, m30},
		    row_proxy<true>{m01, m11, m21, m31},
		    row_proxy<true>{m02, m12, m22, m32},
		    row_proxy<true>{m03, m13, m23, m33},
		}[i];
	}
	constexpr diag_proxy<false> diag()
	{
		return diag_proxy<false>{m00, m11, m22, m33};
	}
	constexpr diag_proxy<true> diag() const
	{
		return diag_proxy<true>{m00, m11, m22, m33};
	}
	constexpr T& diag(int i)
	{
		return std::array{
		    std::ref(m00),
		    std::ref(m11),
		    std::ref(m22),
		    std::ref(m33),
		}[i];
	}
	constexpr const T& diag(int i) const
	{
		return std::array{
		    std::cref(m00),
		    std::cref(m11),
		    std::cref(m22),
		    std::cref(m33),
		}[i];
	}

	matrix4_T operator*(T f) const;
	matrix4_T operator/(T f) const;
	vector4_T<T> operator*(vector4_T<T> v) const;
	matrix4_T operator*(const matrix4_T& m) const;

	bool operator==(const matrix4_T& m) const
	{
		return m00 == m.m00 && m01 == m.m01 && m02 == m.m02 && m03 == m.m03 &&
		       m10 == m.m10 && m11 == m.m11 && m12 == m.m12 && m13 == m.m13 &&
		       m20 == m.m20 && m21 == m.m21 && m22 == m.m22 && m23 == m.m23 &&
		       m30 == m.m30 && m31 == m.m31 && m32 == m.m32 && m33 == m.m33;
	}

	using underlying_type = T;
};

template <arithmetic T>
matrix4_T<T> operator*(const transform3_T<T>& t, const matrix4_T<T>& m);
template <arithmetic T>
matrix4_T<T> operator*(const matrix4_T<T>& m, const transform3_T<T>& t);
template <arithmetic T>
matrix4_T<T> operator*(const uniform_transform3_T<T>& t,
                       const matrix4_T<T>& m);
template <arithmetic T>
matrix4_T<T> operator*(const matrix4_T<T>& m,
                       const uniform_transform3_T<T>& t);
template <arithmetic T>
matrix4_T<T> operator*(const unscaled_transform3_T<T>& t,
                       const matrix4_T<T>& m);
template <arithmetic T>
matrix4_T<T> operator*(const matrix4_T<T>& m,
                       const unscaled_transform3_T<T>& t);
#pragma endregion

#pragma region declaration_perspective_t
template <arithmetic T>
struct perspective_T
{
private:
	radian_T<T> m_angle;
	T m_ratio;
	T m_invRatio;
	T m_near;
	T m_far;
	T m_focalLength;

public:
	perspective_T(radian_T<T> angle, T ratio, T near, T far);
	template <std::signed_integral U, U NumeratorT, U DenominatorT>
	perspective_T(static_pi_fraction_T<U, NumeratorT, DenominatorT> angle,
	              T ratio,
	              T near,
	              T far);

	void set(radian_T<T> angle, T ratio, T near, T far);
	template <std::signed_integral U, U NumeratorT, U DenominatorT>
	void set(static_pi_fraction_T<U, NumeratorT, DenominatorT> angle,
	         T ratio,
	         T near,
	         T far);

	void set_angle(radian_T<T> angle);
	template <std::signed_integral U, U NumeratorT, U DenominatorT>
	void set_angle(static_pi_fraction_T<U, NumeratorT, DenominatorT> angle);
	radian_T<T> get_angle() const;

	void set_ratio(T ratio);
	T get_ratio() const;

	void set_near(T near_plane_t);
	T get_near() const;

	void set_far(T far_plane_t);
	T get_far() const;

	matrix4_T<T> to_matrix4() const;
	matrix4_T<T> to_inverse_matrix4() const;

	template <arithmetic U>
	friend matrix4_T<U> operator*(const matrix4_T<U>& m,
	                              const perspective_T<U>& p);
	template <arithmetic U>
	friend matrix4_T<U> operator*(const perspective_T<U>& p,
	                              const matrix4_T<U>& m);
	template <arithmetic U>
	friend vector4_T<U> operator*(const perspective_T<U>& p,
	                              const vector4_T<U>& v);
	template <arithmetic U>
	friend matrix4_T<U> operator*(const unscaled_transform3_T<U>& t,
	                              const perspective_T<U>& p);
	template <arithmetic U>
	friend matrix4_T<U> operator*(const perspective_T<U>& p,
	                              const unscaled_transform3_T<U>& t);

	using underlying_type = T;
};
#pragma endregion

#pragma region declaration_direction2_T
template <arithmetic T>
class direction2_T
{
public:
	using underlying_type = T;
	static constexpr size_t component_count = 2;

	constexpr direction2_T() = default;
	constexpr direction2_T(T x_, T y_);
	explicit constexpr direction2_T(const vector2_T<T>& v);

	constexpr const T& x() const;
	constexpr const T& y() const;

	constexpr std::array<T, 2> to_array() const;

	template <typename U>
	constexpr operator direction2_T<U>() const;

	constexpr direction2_T operator-() const;

	constexpr bool operator==(const direction2_T& v) const;
	constexpr bool operator!=(const direction2_T& v) const;

	constexpr bool nearly_equal(const direction2_T& v, T e = T(epsilon)) const;

	template <typename Index>
	constexpr const T& operator[](Index index) const;

	static constexpr direction2_T already_normalized(T x_, T y_);
	static constexpr direction2_T already_normalized(const vector2_T<T>& v);

	static constexpr direction2_T posX();
	static constexpr direction2_T posY();
	static constexpr direction2_T negX();
	static constexpr direction2_T negY();

private:
	T m_x;
	T m_y;
};

template <arithmetic T>
constexpr vector2_T<T> operator*(T t, const direction2_T<T>& d);
template <arithmetic T>
constexpr vector2_T<T> operator*(const direction2_T<T>& d, T t);

template <arithmetic T>
constexpr T dot(direction2_T<T> lhs, direction2_T<T> rhs);
template <arithmetic T>
constexpr T cross(direction2_T<T> lhs, direction2_T<T> rhs);

template <arithmetic T>
constexpr T dot(direction2_T<T> lhs, vector2_T<T> rhs);
template <arithmetic T>
constexpr T cross(direction2_T<T> lhs, vector2_T<T> rhs);

template <arithmetic T>
constexpr T dot(vector2_T<T> lhs, direction2_T<T> rhs);
template <arithmetic T>
constexpr T cross(vector2_T<T> lhs, direction2_T<T> rhs);
#pragma endregion

#pragma region declaration_direction3_T
template <arithmetic T>
class direction3_T
{
public:
	using underlying_type = T;
	static constexpr size_t component_count = 3;

	constexpr direction3_T() = default;
	constexpr direction3_T(T x_, T y_, T z_);
	explicit constexpr direction3_T(const vector3_T<T>& v);

	constexpr const T& x() const;
	constexpr const T& y() const;
	constexpr const T& z() const;

	constexpr std::array<T, 3> to_array() const;

	template <typename U>
	constexpr operator direction3_T<U>() const;

	constexpr direction3_T operator-() const;

	constexpr bool operator==(const direction3_T& v) const;
	constexpr bool operator!=(const direction3_T& v) const;

	constexpr bool nearly_equal(const direction3_T& v, T e = T(epsilon)) const;

	template <typename Index>
	constexpr const T& operator[](Index index) const;

	static constexpr direction3_T already_normalized(T x_, T y_, T z_);
	static constexpr direction3_T already_normalized(const vector3_T<T>& v);

	static constexpr direction3_T posX();
	static constexpr direction3_T posY();
	static constexpr direction3_T posZ();
	static constexpr direction3_T negX();
	static constexpr direction3_T negY();
	static constexpr direction3_T negZ();

private:
	T m_x;
	T m_y;
	T m_z;
};

template <arithmetic T>
constexpr vector3_T<T> operator*(T t, const direction3_T<T>& d);
template <arithmetic T>
constexpr vector3_T<T> operator*(const direction3_T<T>& d, T t);

template <arithmetic T>
constexpr T dot(direction3_T<T> lhs, direction3_T<T> rhs);
template <arithmetic T>
constexpr vector3_T<T> cross(direction3_T<T> lhs, direction3_T<T> rhs);

template <arithmetic T>
constexpr T dot(direction3_T<T> lhs, vector3_T<T> rhs);
template <arithmetic T>
constexpr vector3_T<T> cross(direction3_T<T> lhs, vector3_T<T> rhs);

template <arithmetic T>
constexpr T dot(vector3_T<T> lhs, direction3_T<T> rhs);
template <arithmetic T>
constexpr vector3_T<T> cross(vector3_T<T> lhs, direction3_T<T> rhs);
#pragma endregion

#pragma region declaration_euler_angles_t
template <arithmetic T>
struct euler_angles_T
{
	radian_T<T> x, y, z;

	using underlying_type = T;
};
#pragma endregion

#pragma region declaration_quaternion_t
template <arithmetic T>
struct quaternion_T
{
	T x, y, z, w;

	quaternion_T() = default;
	constexpr quaternion_T(T x_, T y_, T z_, T w_) : x{x_}, y{y_}, z{z_}, w{w_} {}
	constexpr quaternion_T(std::array<T, 4> a) : quaternion_T{a[0], a[1], a[2], a[3]} {}
	quaternion_T(direction3_T<T> axis, radian_T<T> angle);
	template <std::signed_integral U, U NumeratorT, U DenominatorT>
	quaternion_T(direction3_T<T> axis,
	             static_pi_fraction_T<U, NumeratorT, DenominatorT> angle);

	template <arithmetic U>
	explicit operator quaternion_T<U>() const
	{
		return quaternion_T<U>{static_cast<U>(x), static_cast<U>(y),
		                       static_cast<U>(z), static_cast<U>(w)};
	}

	operator std::array<T, 4>() const { return {x, y, z, w}; }

	template <arithmetic U = T>
	std::array<U, 4> to_array() const
	{
		return {static_cast<U>(x), static_cast<U>(y), static_cast<U>(z),
		        static_cast<U>(w)};
	}

	static quaternion_T identity();
	static quaternion_T from_rotation_matrix(const matrix3_T<T>& m);
	static quaternion_T look_at(direction3_T<T> forward, direction3_T<T> upward = direction3_T<T>::posY());

	T norm() const;
	T norm_squared() const;

	quaternion_T& inverse();
	quaternion_T& conjugate();
	quaternion_T& normalize();
	constexpr quaternion_T& clamp(quaternion_T min, quaternion_T max);

	quaternion_T get_inversed() const;
	quaternion_T get_conjugated() const;
	quaternion_T get_normalized() const;
	constexpr quaternion_T get_clamped(quaternion_T min, quaternion_T max) const;

	euler_angles_T<T> to_euler_angles() const;

	quaternion_T operator+(quaternion_T q) const;
	quaternion_T operator-(quaternion_T q) const;

	// TODO: implement operator* for normalized<>
	// if result is normalized too
	quaternion_T operator*(quaternion_T q) const;

	vector3_T<T> operator*(vector3_T<T> v) const;
	quaternion_T operator*(T f) const;
	quaternion_T operator/(T f) const;

	quaternion_T& operator+=(quaternion_T q);
	quaternion_T& operator-=(quaternion_T q);
	quaternion_T& operator*=(quaternion_T q);
	quaternion_T& operator*=(T f);
	quaternion_T& operator/=(T f);

	quaternion_T operator-() const;

	bool operator==(quaternion_T v) const;
	bool operator!=(quaternion_T v) const;

	using underlying_type = T;
};

template <arithmetic T>
constexpr T dot(quaternion_T<T> lhs, quaternion_T<T> rhs);
template <arithmetic T>
constexpr quaternion_T<T> cross(quaternion_T<T> lhs, quaternion_T<T> rhs);
template <arithmetic T>
quaternion_T<T> lerp(quaternion_T<T> begin, quaternion_T<T> end, T percent);
template <arithmetic T>
quaternion_T<T> nlerp(quaternion_T<T> begin, quaternion_T<T> end, T percent);
template <arithmetic T>
quaternion_T<T> slerp(quaternion_T<T> begin, quaternion_T<T> end, T percent);
#pragma endregion

#pragma region declaration_line_t
template <arithmetic T, line_representation representation>
struct line_T {};

#pragma region slope_intercept
template <arithmetic T>
struct line_T<T, line_representation::SlopeIntercept>
{
	T slope;
	T y_intercept;

	line_T() = default;
	line_T(T s, T yi);
	line_T(vector2_T<T> direction, T e = T(epsilon));
	line_T(vector2_T<T> point1, vector2_T<T> point2);

	vector2_T<T> resolve_for_x(T x) const;
	vector2_T<T> resolve_for_y(T y) const;

	using underlying_type = T;
};
#pragma endregion

#pragma region points
template <arithmetic T>
struct line_T<T, line_representation::Points>
{
	vector2_T<T> p1;
	vector2_T<T> p2;

	line_T() = default;
	line_T(vector2_T<T> point1, vector2_T<T> point2);

	//vector2_T<T> resolve_for_x(T x) const;
	//vector2_T<T> resolve_for_y(T y) const;

	using underlying_type = T;
};
#pragma endregion

//template <arithmetic T, line_representation representation1, line_representation representation2>
//bool are_parallel(line_T<T, representation1> l1, line_T<T, representation2> l2, T e = T(epsilon));
//
//template <arithmetic T, line_representation representation1, line_representation representation2>
//bool collides(line_T<T, representation1> l1, line_T<T, representation2> l2);
//template <arithmetic T, line_representation representation1, line_representation representation2>
//bool collides(line_T<T, representation1> l1, line_T<T, representation2> l2, vector2_T<T>& intersection);
#pragma endregion

#pragma region declaration_segment2_t
template <arithmetic T>
struct segment2_T
{
	vector2_T<T> origin;
	vector2_T<T> end;

	T length() const;
	T length_squared() const;

	segment2_T& invert();

	using underlying_type = T;
};

template <arithmetic T>
int collides(const segment2_T<T>& s0,
             const segment2_T<T>& s1,
             vector2_T<T>& outPoint0,
             vector2_T<T>& outPoint1,
             T e = T(epsilon));
#pragma endregion

#pragma region declaration_segment3_t
template <arithmetic T>
struct segment3_T
{
	vector3_T<T> origin;
	vector3_T<T> end;

	using underlying_type = T;
};

template <arithmetic T>
constexpr bool collides(const segment3_T<T>& s,
                        const plane_T<T>& p,
                        T e = T(epsilon)) noexcept;
template <arithmetic T>
constexpr bool collides(const segment3_T<T>& s,
                        const plane_T<T>& p,
                        vector3_T<T>& outPoint,
                        T e = T(epsilon)) noexcept;
#pragma endregion

#pragma region declaration_plane_t
template <arithmetic T>
struct plane_T
{
	vector3_T<T> origin;
	direction3_T<T> normal;

	template <arithmetic U>
	explicit operator plane_T<U>() const
	{
		return plane_T<U>{static_cast<vector3_T<U>>(origin),
		                  static_cast<direction3_T<U>>(normal)};
	}

	constexpr T evaluate(const vector3_T<T>& point) const;
	constexpr vector3_T<T> project3d(const vector3_T<T>& point) const;
	constexpr vector2_T<T> project2d(const vector3_T<T>& point,
	                       const direction3_T<T>& plane_tangent) const;
	constexpr vector3_T<T> unproject(const vector2_T<T>& point,
	                       const direction3_T<T>& plane_tangent) const;

	constexpr vector3_T<T> symmetry(const vector3_T<T>& point) const;

	constexpr direction3_T<T> computeTangent(T epsilon_ = T(epsilon)) const;

	using underlying_type = T;
};
#pragma endregion

#pragma region declaration_oriented_plane_t
template <arithmetic T>
struct oriented_plane_T
{
	vector3_T<T> origin;
	direction3_T<T> normal;
	direction3_T<T> tangent;

	oriented_plane_T() = default;
	oriented_plane_T(const oriented_plane_T&) = default;
	oriented_plane_T(oriented_plane_T&&) = default;
	oriented_plane_T& operator=(const oriented_plane_T&) = default;
	oriented_plane_T& operator=(oriented_plane_T&&) = default;

	constexpr oriented_plane_T(const vector3_T<T>& origin_, const direction3_T<T>& normal_, const direction3_T<T>& tangent_);
	constexpr oriented_plane_T(const plane_T<T>& plane, const direction3_T<T>& tangent_);

	template <arithmetic U>
	explicit operator oriented_plane_T<U>() const
	{
		return oriented_plane_T<U>{
		    static_cast<vector3_T<U>>(origin),
		    static_cast<direction3_T<U>>(normal),
		    static_cast<direction3_T<U>>(tangent),
		};
	}

	constexpr T evaluate(const vector3_T<T>& point) const;
	constexpr vector3_T<T> project3d(const vector3_T<T>& point) const;
	constexpr vector2_T<T> project2d(const vector3_T<T>& point) const;
	constexpr vector3_T<T> unproject(const vector2_T<T>& point) const;

	constexpr vector3_T<T> symmetry(const vector3_T<T>& point) const;

	using underlying_type = T;
};
#pragma endregion

#pragma region declaration_ray2_t
template <arithmetic T>
struct ray2_T
{
	vector2_T<T> origin;
	direction2_T<T> direction;

	using underlying_type = T;
};

template <arithmetic T>
constexpr bool collides(const ray2_T<T>& r, const segment2_T<T>& s);

template <arithmetic T>
constexpr std::optional<vector2_T<T>> intersection(const ray2_T<T>& r,
                                                  const segment2_T<T>& s);
#pragma endregion

#pragma region declaration_ray3_t
template <arithmetic T>
struct ray3_T
{
	vector3_T<T> origin;
	direction3_T<T> direction;

	using underlying_type = T;
};
#pragma endregion

#pragma region declaration_circle_t
template <arithmetic T>
struct circle_T
{
	vector2_T<T> origin;
	T radius;

	using underlying_type = T;
};

template <arithmetic T>
bool collides(circle_T<T> c, vector2_T<T> p);
template <arithmetic T>
bool collides(vector2_T<T> p, circle_T<T> c);
template <arithmetic T>
bool collides(ray3_T<T> r, plane_T<T> p);
template <arithmetic T>
bool collides(ray3_T<T> r, plane_T<T> p, vector3_T<T>& intersection);
template <arithmetic T>
bool collides(plane_T<T> p, ray3_T<T> r);
template <arithmetic T>
bool collides(plane_T<T> p, ray3_T<T> r, vector3_T<T>& intersection);
template <arithmetic T>
bool collides_bidirectional(const plane_T<T>& plane,
                            ray3_T<T> r,
                            vector3_T<T>& collision_point);
#pragma endregion

#pragma region declaration_aabb2_t
template <arithmetic T>
struct aabb2_T
{
	vector2_T<T> min;
	vector2_T<T> max;

	constexpr bool contains(vector2_T<T> p) const;

	constexpr vector2_T<T> get_center() const;
	constexpr std::array<vector2_T<T>, 4> get_points() const;

	constexpr vector2_T<T> clamp(vector2_T<T> p) const;

	constexpr aabb2_T& grow(const aabb2_T& box);
	constexpr aabb2_T& grow(vector2_T<T> point);

	constexpr aabb2_T operator+(const aabb2_T& rhs) const;
	constexpr aabb2_T operator+(vector2_T<T> rhs) const;

	constexpr aabb2_T& operator+=(const aabb2_T& rhs);
	constexpr aabb2_T& operator+=(vector2_T<T> rhs);

	static constexpr aabb2_T infinity();

	using underlying_type = T;
};

template <arithmetic T>
constexpr aabb2_T<T> operator+(vector2_T<T> lhs, const aabb2_T<T>& rhs);

template <arithmetic T>
constexpr bool collides(const aabb2_T<T>& b0, const aabb2_T<T>& b1);
#pragma endregion

#pragma region declaration_aabb3_t
template <arithmetic T>
struct aabb3_T
{
	vector3_T<T> min;
	vector3_T<T> max;

	constexpr bool contains(vector3_T<T> p) const;

	constexpr vector3_T<T> get_center() const;
	constexpr std::array<vector3_T<T>, 8> get_points() const;

	constexpr vector3_T<T> clamp(vector3_T<T> p) const;

	constexpr aabb3_T& grow(const aabb3_T& box);
	constexpr aabb3_T& grow(vector3_T<T> point);

	constexpr aabb3_T operator+(const aabb3_T& rhs) const;
	constexpr aabb3_T operator+(vector3_T<T> rhs) const;

	constexpr aabb3_T& operator+=(const aabb3_T& rhs);
	constexpr aabb3_T& operator+=(vector3_T<T> rhs);

	static constexpr aabb3_T infinity();

	using underlying_type = T;
};

template <arithmetic T>
constexpr aabb3_T<T> operator+(vector3_T<T> lhs, const aabb3_T<T>& rhs);

template <arithmetic T>
constexpr bool collides(const aabb3_T<T>& b, const ray3_T<T>& r);
template <arithmetic T>
constexpr bool collides(const aabb3_T<T>& b, const plane_T<T>& p);
template <arithmetic T>
constexpr bool collides(const aabb3_T<T>& b, const oriented_plane_T<T>& p);
template <arithmetic T>
constexpr bool collides(const aabb3_T<T>& b0, const aabb3_T<T>& b1);
#pragma endregion

#pragma region declaration_transform2_T
template <arithmetic T>
struct transform2_T
{
	vector2_T<T> position;
	normalized<complex_T<T>> rotation;
	vector2_T<T> scale;

	transform2_T operator*(const transform2_T& t) const;
	vector2_T<T> operator*(vector2_T<T> v) const;

	static transform2_T identity();

	using underlying_type = T;
};
#pragma endregion

#pragma region declaration_transform3_T
template <arithmetic T>
struct transform3_T
{
	vector3_T<T> position;
	quaternion_T<T> rotation;
	vector3_T<T> scale;

	transform3_T& translate_absolute(vector3_T<T> t);
	transform3_T& translate_relative(vector3_T<T> t);
	transform3_T& rotate(quaternion_T<T> r);

	matrix4_T<T> to_matrix4() const;

	transform3_T operator*(const transform3_T& t) const;
	vector3_T<T> operator*(vector3_T<T> v) const;
	quaternion_T<T> operator*(quaternion_T<T> q) const;

	static transform3_T identity();

	using underlying_type = T;
};
#pragma endregion

#pragma region declaration_uniform_transform2_T
template <arithmetic T>
struct uniform_transform2_T
{
	vector2_T<T> position;
	radian_T<T> rotation;
	T scale;

	uniform_transform2_T operator*(uniform_transform2_T t) const;
	vector2_T<T> operator*(vector2_T<T> v) const;

	static uniform_transform2_T identity();

	using underlying_type = T;
};
#pragma endregion

#pragma region declaration_uniform_transform3_T
template <arithmetic T>
struct uniform_transform3_T
{
	vector3_T<T> position;
	quaternion_T<T> rotation;
	T scale;

	matrix4_T<T> to_matrix4() const;

	uniform_transform3_T operator*(const uniform_transform3_T& t) const;
	vector3_T<T> operator*(vector3_T<T> v) const;
	quaternion_T<T> operator*(quaternion_T<T> q) const;

	static uniform_transform3_T identity();

	using underlying_type = T;
};
#pragma endregion

#pragma region declaration_unscaled_transform2_T
template <arithmetic T>
struct unscaled_transform2_T
{
	vector2_T<T> position;
	radian_T<T> rotation;

	unscaled_transform2_T operator*(unscaled_transform2_T t) const;
	vector2_T<T> operator*(vector2_T<T> v) const;

	static unscaled_transform2_T identity();

	using underlying_type = T;
};
#pragma endregion

#pragma region declaration_unscaled_transform3_T
template <arithmetic T>
struct unscaled_transform3_T
{
	vector3_T<T> position;
	quaternion_T<T> rotation;

	unscaled_transform3_T& translate_absolute(vector3_T<T> t);
	unscaled_transform3_T& translate_relative(vector3_T<T> t);
	unscaled_transform3_T& rotate(quaternion_T<T> r);

	unscaled_transform3_T& inverse();
	unscaled_transform3_T get_inversed() const;

	matrix4_T<T> to_matrix4() const;

	unscaled_transform3_T operator*(const unscaled_transform3_T& t) const;
	vector3_T<T> operator*(vector3_T<T> v) const;
	quaternion_T<T> operator*(quaternion_T<T> q) const;

	static unscaled_transform3_T identity();

	using underlying_type = T;
};

template <arithmetic T>
unscaled_transform3_T<T> inverse(unscaled_transform3_T<T> tr);
#pragma endregion

#pragma region declaration_value_domain
template<typename T> requires arithmetic<T> || vector<T>
class value_domain_T
{
	T m_lim0;
	T m_lim1;

public:
	value_domain_T() = default;
	constexpr value_domain_T(T lim0, T lim1) noexcept
	    : m_lim0{lim0}, m_lim1{lim1}
	{
	}
	value_domain_T(const value_domain_T&) = default;
	value_domain_T(value_domain_T&&) = default;
	~value_domain_T() = default;

	constexpr T lim0() const noexcept { return m_lim0; }
	constexpr T lim1() const noexcept { return m_lim1; }
	constexpr T range() const noexcept
	{
		return element_wise_abs(lim0() - lim1());
	}
	constexpr T min() const noexcept
	{
		return element_wise_min(lim0(), lim1());
	}
	constexpr T max() const noexcept
	{
		return element_wise_max(lim0(), lim1());
	}

	constexpr T clamp(T value) const noexcept
	{
		return cdm::clamp(value, min(), max());
	}
	constexpr T lerp(normalized_value_T<T> value) const noexcept;

	value_domain_T& operator=(const value_domain_T&) = default;
	value_domain_T& operator=(value_domain_T&&) = default;

	using underlying_type = T;
};
#pragma endregion

#pragma region declaration_unnormalized_value
template<typename T>
class unnormalized_value_T
{
	value_domain_T<T> m_domain;
	T m_value;

public:
	unnormalized_value_T() = default;
	constexpr explicit unnormalized_value_T(value_domain_T<T> domain,
	                                      T value) noexcept
	    : m_domain{domain}, m_value{m_domain.clamp(value)}
	{
	}
	constexpr explicit unnormalized_value_T(value_domain_T<T> domain) noexcept
	    : m_domain{domain}, m_value{m_domain.lim0()}
	{
	}
	constexpr explicit unnormalized_value_T(value_domain_T<T> domain,
	                                      normalized_value_T<T> value) noexcept;
	unnormalized_value_T(const unnormalized_value_T&) = default;
	unnormalized_value_T(unnormalized_value_T&&) = default;
	~unnormalized_value_T() = default;

	constexpr const value_domain_T<T>& domain() const noexcept
	{
		return m_domain;
	}
	constexpr const T& value() const noexcept { return m_value; }
	constexpr void setValue(T value) noexcept
	{
		m_value = m_domain.clamp(value);
	}

	constexpr normalized_value_T<T> to_normalized() const noexcept;

	constexpr unnormalized_value_T& operator=(
	    normalized_value_T<T> value) noexcept;
	unnormalized_value_T& operator=(const unnormalized_value_T&) = default;
	unnormalized_value_T& operator=(unnormalized_value_T&&) = default;

	using underlying_type = T;
};
#pragma endregion

#pragma region declaration_normalized_value
template<typename T>
class normalized_value_T
{
	T m_value;

public:
	normalized_value_T() = default;
	constexpr explicit normalized_value_T(T value) noexcept
	    : m_value{cdm::clamp<T>(value, zero<T>(), one<T>())}
	{
	}
	constexpr explicit normalized_value_T(unnormalized_value_T<T> value) noexcept
	    : m_value{(value.value() - value.domain().min()) /
	              value.domain().range()}
	{
		if (value.domain().min() != value.domain().lim0())
			inverse();
	}
	normalized_value_T(const normalized_value_T&) = default;
	normalized_value_T(normalized_value_T&&) = default;
	~normalized_value_T() = default;

	constexpr const T& value() const noexcept { return m_value; }

	constexpr normalized_value_T& inverse() noexcept
	{
		m_value = T(1) - m_value;
		return *this;
	}

	constexpr normalized_value_T get_inversed() const noexcept
	{
		auto res{*this};
		res.inverse();
		return res;
	}

	constexpr normalized_value_T& operator=(T value) noexcept
	{
		m_value = cdm::clamp<T>(value, T(0), T(1));
		return *this;
	}
	normalized_value_T& operator=(const normalized_value_T&) = default;
	normalized_value_T& operator=(normalized_value_T&&) = default;

	using underlying_type = T;
};
#pragma endregion

template <arithmetic T>
constexpr T domain_transfer(cdm::value_domain_T<T> from,
                            cdm::value_domain_T<T> to,
                            T value) noexcept
{
	const cdm::unnormalized_value_T<T> valueFrom{from, value};
	const cdm::normalized_value_T<T> norm{valueFrom};
	const cdm::unnormalized_value_T<T> valueTo{to, norm};
	return valueTo.value();
}
#pragma endregion

// ==========================================================================

#pragma region definitions
#pragma region definition_complex_T
template <arithmetic T>
complex_T<T>::complex_T(radian_T<T> angle) : r{cos(angle)}, i{sin(angle)}
{
}

template <arithmetic T>
T complex_T<T>::norm() const
{
	return std::sqrt(norm_squared());
}
template <arithmetic T>
T complex_T<T>::norm_squared() const
{
	return r * r + i * i;
}

template <arithmetic T>
complex_T<T>& complex_T<T>::normalize()
{
	const T n = norm();
	r /= n;
	i /= n;
	return *this;
}
template <arithmetic T>
complex_T<T> complex_T<T>::get_normalized() const
{
	complex_T<T> res{*this};
	res.normalize();
	return *this;
}

template <arithmetic T>
complex_T<T>& complex_T<T>::conjugate()
{
	i = -i;
	return *this;
}
template <arithmetic T>
complex_T<T> complex_T<T>::get_conjugated() const
{
	complex_T<T> res{*this};
	res.conjugate();
	return *this;
}

template <arithmetic T>
complex_T<T> complex_T<T>::operator+(complex_T<T> c) const
{
	return {r + c.r, i + c.i};
}
template <arithmetic T>
complex_T<T> complex_T<T>::operator-(complex_T<T> c) const
{
	return {r - c.r, i - c.i};
}

template <arithmetic T>
complex_T<T>& complex_T<T>::operator+=(complex_T<T> c)
{
	return *this = *this + c;
}
template <arithmetic T>
complex_T<T>& complex_T<T>::operator-=(complex_T<T> c)
{
	return *this = *this - c;
}
template <arithmetic T>
complex_T<T>& complex_T<T>::operator*=(complex_T<T> c)
{
	return *this = *this * c;
}

template <arithmetic T>
complex_T<T> complex_T<T>::identity() { return {T(1), T(0)}; }

template <arithmetic T>
complex_T<T> operator*(complex_T<T> c1, complex_T<T> c2)
{
	return {c1.r * c2.r - c1.i * c2.i, c1.r * c2.i + c1.i * c2.r};
}
template <arithmetic T>
complex_T<T> operator*(complex_T<T> c, T f)
{
	return {c.r * f, c.i * f};
}
template <arithmetic T>
complex_T<T> operator*(normalized<complex_T<T>> c1, complex_T<T> c2)
{
	return *c1 * c2;
}
template <arithmetic T>
complex_T<T> operator*(complex_T<T> c1, normalized<complex_T<T>> c2)
{
	return c1 * *c2;
}
template <arithmetic T>
normalized<complex_T<T>> operator*(normalized<complex_T<T>> c1,
                                   normalized<complex_T<T>> c2)
{
	return normalized<complex_T<T>>::already_normalized(complex_T<T>{
	    c1->r * c2->r - c1->i * c2->i, c1->r * c2->i + c1->i * c2->r});
}
template <arithmetic T>
vector2_T<T> operator*(normalized<complex_T<T>> c, vector2_T<T> v)
{
	return {c->r * v.x - c->i * v.y, c->r * v.y + c->i * v.x};
}

template <arithmetic T>
T norm(complex_T<T> c)
{
	return c.norm_squared();
}
template <arithmetic T>
T norm_squared(complex_T<T> c)
{
	return c.norm();
}
template <arithmetic T>
complex_T<T> normalize(complex_T<T> c)
{
	c.normalize();
	return c;
}
template <arithmetic T>
complex_T<T> conjugate(complex_T<T> c)
{
	c.conjugate();
	return c;
}
#pragma endregion

#pragma region definition_radian_T
template <arithmetic T>
constexpr radian_T<T>::radian_T(T f) : angle(f)
{
}
template <arithmetic T>
constexpr radian_T<T>::radian_T(degree_T<T> d)
    : angle(static_cast<T>(d) * T(deg_to_rad))
{
}

template <arithmetic T>
radian_T<T>::operator T() const
{
	return angle;
}

template <arithmetic T>
radian_T<T>& radian_T<T>::operator+=(radian_T<T> r)
{
	angle += r.angle;
	return *this;
}
template <arithmetic T>
radian_T<T>& radian_T<T>::operator-=(radian_T<T> r)
{
	angle -= r.angle;
	return *this;
}
template <arithmetic T>
radian_T<T>& radian_T<T>::operator*=(radian_T<T> r)
{
	angle *= r.angle;
	return *this;
}
template <arithmetic T>
radian_T<T>& radian_T<T>::operator/=(radian_T<T> r)
{
	angle /= r.angle;
	return *this;
}

template <arithmetic T>
radian_T<T> radian_T<T>::operator-() const
{
	return radian_T<T>{-angle};
}

template <arithmetic T>
radian_T<T>& radian_T<T>::operator=(degree_T<T> d)
{
	angle = d.angle * T(deg_to_rad);
	return *this;
}

template <arithmetic T>
radian_T<T> operator+(radian_T<T> r1, radian_T<T> r2)
{
	return radian_T<T>{static_cast<T>(r1) + static_cast<T>(r2)};
}
template <arithmetic T>
radian_T<T> operator-(radian_T<T> r1, radian_T<T> r2)
{
	return radian_T<T>{static_cast<T>(r1) - static_cast<T>(r2)};
}
template <arithmetic T>
radian_T<T> operator*(T f, radian_T<T> r)
{
	return radian_T<T>{static_cast<T>(r) * f};
}
template <arithmetic T>
radian_T<T> operator*(radian_T<T> r, T f)
{
	return radian_T<T>{f * static_cast<T>(r)};
}
template <arithmetic T>
radian_T<T> operator/(radian_T<T> r, T f)
{
	return radian_T<T>{static_cast<T>(r) / f};
}
template <arithmetic T>
radian_T<T>& operator+=(radian_T<T>& r, T f)
{
	r = r + f;
	return r;
}
template <arithmetic T>
radian_T<T>& operator-=(radian_T<T>& r, T f)
{
	r = r - f;
	return r;
}
template <arithmetic T>
radian_T<T>& operator*=(radian_T<T>& r, T f)
{
	r = r * f;
	return r;
}
template <arithmetic T>
radian_T<T>& operator/=(radian_T<T>& r, T f)
{
	r = r / f;
	return r;
}

template <arithmetic T>
bool operator<(radian_T<T> lhs, radian_T<T> rhs)
{
	return T(lhs) < T(rhs);
}
template <arithmetic T>
bool operator>(radian_T<T> lhs, radian_T<T> rhs)
{
	return T(lhs) > T(rhs);
}
template <arithmetic T>
bool operator==(radian_T<T> lhs, radian_T<T> rhs)
{
	return T(lhs) == T(rhs);
}
template <arithmetic T>
bool operator!=(radian_T<T> lhs, radian_T<T> rhs)
{
	return T(lhs) != T(rhs);
}
template <arithmetic T>
bool operator>=(radian_T<T> lhs, radian_T<T> rhs)
{
	return T(lhs) >= T(rhs);
}
template <arithmetic T>
bool operator<=(radian_T<T> lhs, radian_T<T> rhs)
{
	return T(lhs) <= T(rhs);
}

template <arithmetic T>
T sin(radian_T<T> r)
{
	return std::sin(static_cast<T>(r));
}
template <arithmetic T>
T cos(radian_T<T> r)
{
	return std::cos(static_cast<T>(r));
}
template <arithmetic T>
T tan(radian_T<T> r)
{
	return std::tan(static_cast<T>(r));
}
template <arithmetic T>
T asin(radian_T<T> r)
{
	return std::asin(static_cast<T>(r));
}
template <arithmetic T>
T acos(radian_T<T> r)
{
	return std::acos(static_cast<T>(r));
}
template <arithmetic T>
T atan(radian_T<T> r)
{
	return std::atan(static_cast<T>(r));
}
template <arithmetic T>
T sinh(radian_T<T> r)
{
	return std::sinh(static_cast<T>(r));
}
template <arithmetic T>
T cosh(radian_T<T> r)
{
	return std::cosh(static_cast<T>(r));
}
template <arithmetic T>
T tanh(radian_T<T> r)
{
	return std::tanh(static_cast<T>(r));
}
template <arithmetic T>
T asinh(radian_T<T> r)
{
	return std::asinh(static_cast<T>(r));
}
template <arithmetic T>
T acosh(radian_T<T> r)
{
	return std::acosh(static_cast<T>(r));
}
template <arithmetic T>
T atanh(radian_T<T> r)
{
	return std::atanh(static_cast<T>(r));
}
#pragma endregion

#pragma region definition_degree_T
template <arithmetic T>
constexpr degree_T<T>::degree_T(T f) : angle(f)
{
}
template <arithmetic T>
constexpr degree_T<T>::degree_T(radian_T<T> r)
    : angle(static_cast<T>(r) * T(rad_to_deg))
{
}

template <arithmetic T>
degree_T<T>::operator T() const
{
	return angle;
}

template <arithmetic T>
degree_T<T>& degree_T<T>::operator+=(degree_T d)
{
	angle += d.angle;
	return *this;
}
template <arithmetic T>
degree_T<T>& degree_T<T>::operator-=(degree_T d)
{
	angle -= d.angle;
	return *this;
}
template <arithmetic T>
degree_T<T>& degree_T<T>::operator*=(degree_T d)
{
	angle *= d.angle;
	return *this;
}
template <arithmetic T>
degree_T<T>& degree_T<T>::operator/=(degree_T d)
{
	angle /= d.angle;
	return *this;
}

template <arithmetic T>
degree_T<T> degree_T<T>::operator-() const
{
	return degree_T(-angle);
}

template <arithmetic T>
degree_T<T> operator+(degree_T<T> d1, degree_T<T> d2)
{
	return degree_T<T>{static_cast<T>(d1) + static_cast<T>(d2)};
}
template <arithmetic T>
degree_T<T> operator-(degree_T<T> d1, degree_T<T> d2)
{
	return degree_T<T>{static_cast<T>(d1) - static_cast<T>(d2)};
}
template <arithmetic T>
degree_T<T> operator*(T f, degree_T<T> d)
{
	return degree_T<T>{static_cast<T>(d) * f};
}
template <arithmetic T>
degree_T<T> operator*(degree_T<T> d, T f)
{
	return degree_T<T>{f * static_cast<T>(d)};
}
template <arithmetic T>
degree_T<T> operator/(degree_T<T> d, T f)
{
	return degree_T<T>{static_cast<T>(d) / f};
}
template <arithmetic T>
degree_T<T>& operator+=(degree_T<T>& d, T f)
{
	d = d + f;
	return d;
}
template <arithmetic T>
degree_T<T>& operator-=(degree_T<T>& d, T f)
{
	d = d - f;
	return d;
}
template <arithmetic T>
degree_T<T>& operator*=(degree_T<T>& d, T f)
{
	d = d * f;
	return d;
}
template <arithmetic T>
degree_T<T>& operator/=(degree_T<T>& d, T f)
{
	d = d / f;
	return d;
}

template <arithmetic T>
bool operator<(degree_T<T> lhs, degree_T<T> rhs)
{
	return T(lhs) < T(rhs);
}
template <arithmetic T>
bool operator>(degree_T<T> lhs, degree_T<T> rhs)
{
	return T(lhs) > T(rhs);
}
template <arithmetic T>
bool operator==(degree_T<T> lhs, degree_T<T> rhs)
{
	return T(lhs) == T(rhs);
}
template <arithmetic T>
bool operator!=(degree_T<T> lhs, degree_T<T> rhs)
{
	return T(lhs) != T(rhs);
}
template <arithmetic T>
bool operator>=(degree_T<T> lhs, degree_T<T> rhs)
{
	return T(lhs) >= T(rhs);
}
template <arithmetic T>
bool operator<=(degree_T<T> lhs, degree_T<T> rhs)
{
	return T(lhs) <= T(rhs);
}

template <arithmetic T>
T sin(degree_T<T> d)
{
	return sin(radian_T<T>(d));
}
template <arithmetic T>
T cos(degree_T<T> d)
{
	return cos(radian_T<T>(d));
}
template <arithmetic T>
T tan(degree_T<T> d)
{
	return tan(radian_T<T>(d));
}
template <arithmetic T>
T asin(degree_T<T> d)
{
	return asin(radian_T<T>(d));
}
template <arithmetic T>
T acos(degree_T<T> d)
{
	return acos(radian_T<T>(d));
}
template <arithmetic T>
T atan(degree_T<T> d)
{
	return atan(radian_T<T>(d));
}
template <arithmetic T>
T sinh(degree_T<T> d)
{
	return sinh(radian_T<T>(d));
}
template <arithmetic T>
T cosh(degree_T<T> d)
{
	return cosh(radian_T<T>(d));
}
template <arithmetic T>
T tanh(degree_T<T> d)
{
	return tanh(radian_T<T>(d));
}
template <arithmetic T>
T asinh(degree_T<T> d)
{
	return asinh(radian_T<T>(d));
}
template <arithmetic T>
T acosh(degree_T<T> d)
{
	return acosh(radian_T<T>(d));
}
template <arithmetic T>
T atanh(degree_T<T> d)
{
	return atanh(radian_T<T>(d));
}
#pragma endregion

#pragma region definition_vector2_t
template <arithmetic T>
constexpr vector2_T<T>::vector2_T(const direction2_T<T>& d) : vector2_T(d.x(), d.y()){}

template <arithmetic T>
vector2_T<T> vector2_T<T>::operator+(vector2_T<T> v) const
{
	return {x + v.x, y + v.y};
}
template <arithmetic T>
vector2_T<T> vector2_T<T>::operator-(vector2_T<T> v) const
{
	return {x - v.x, y - v.y};
}

template <arithmetic T>
vector2_T<T>& vector2_T<T>::operator+=(vector2_T<T> v)
{
	*this = *this + v;
	return *this;
}
template <arithmetic T>
vector2_T<T>& vector2_T<T>::operator-=(vector2_T<T> v)
{
	*this = *this - v;
	return *this;
}

template <arithmetic T>
vector2_T<T> vector2_T<T>::operator-() const
{
	return {-x, -y};
}

template <arithmetic T>
bool vector2_T<T>::operator==(vector2_T<T> v) const
{
	return x == v.x && y == v.y;
}
template <arithmetic T>
bool vector2_T<T>::operator!=(vector2_T<T> v) const
{
	return !operator==(v);
}

template <arithmetic T>
vector2_T<T> operator*(T f, vector2_T<T> v)
{
	return v * f;
}
template <arithmetic T>
vector2_T<T> operator*(vector2_T<T> v, T f)
{
	return {v.x * f, v.y * f};
}
template <arithmetic T>
vector2_T<T> operator/(vector2_T<T> v, T f)
{
	return {v.x / f, v.y / f};
}
template <arithmetic T>
vector2_T<T>& operator*=(vector2_T<T>& v, T f)
{
	v = v * f;
	return v;
}
template <arithmetic T>
vector2_T<T>& operator/=(vector2_T<T>& v, T f)
{
	v = v / f;
	return v;
}

template <arithmetic T>
T norm(vector2_T<T> v)
{
	return std::sqrt(norm_squared(v));
}
template <arithmetic T>
T norm_squared(vector2_T<T> v)
{
	return v.x * v.x + v.y * v.y;
}
template <arithmetic T>
vector2_T<T> normalize(vector2_T<T> v)
{
	T n = norm(v);
	v.x /= n;
	v.y /= n;
	return v;
}
//template <arithmetic T>
//constexpr vector2_T<T> clamp(vector2_T<T> v,
//                             vector2_T<T> min,
//                             vector2_T<T> max)
//{
//	v.clamp(min, max);
//	return v;
//}
template <arithmetic T>
constexpr T dot(vector2_T<T> lhs, vector2_T<T> rhs)
{
	return lhs.x * rhs.x + lhs.y * rhs.y;
}
template <arithmetic T>
constexpr T cross(vector2_T<T> lhs, vector2_T<T> rhs)
{
	return lhs.x * rhs.y - lhs.y * rhs.x;
}
template <arithmetic T>
constexpr vector2_T<T> lerp(vector2_T<T> begin, vector2_T<T> end, T percent)
{
	//return (end - begin) * percent + begin;
	return {
		std::lerp(begin.x, end.x, percent),
		std::lerp(begin.y, end.y, percent),
	};
}
template <arithmetic T>
vector2_T<T> nlerp(vector2_T<T> begin, vector2_T<T> end, T percent)
{
	return lerp(begin, end, percent).get_normalized();
}
template <arithmetic T>
vector2_T<T> slerp(vector2_T<T> begin, vector2_T<T> end, T percent)
{
	const radian_T angle = angle_between(begin, end) * percent;
	const T s = sin(angle);
	const T c = cos(angle);

	normalized<vector2_T<T>> res{
	    {c * begin.x - s * begin.y, s * begin.x + c * begin.y}};

	T f = cdm::lerp(norm(begin), norm(end), percent);
	return res * f;
}
template <arithmetic T>
T distance_between(vector2_T<T> v1, vector2_T<T> v2)
{
	return norm(v1 - v2);
}
template <arithmetic T>
T distance_squared_between(vector2_T<T> v1, vector2_T<T> v2)
{
	return norm_squared(v1 - v2);
}
template <arithmetic T>
vector2_T<T> from_to(vector2_T<T> from, vector2_T<T> to)
{
	return {to.x - from.x, to.y - from.y};
}
template <arithmetic T>
radian_T<T> angle_between(vector2_T<T> v1, vector2_T<T> v2)
{
	return radian_T{atan2(v2.y, v2.x) - atan2(v1.y, v1.x)};
}
template <arithmetic T>
bool nearly_equal(vector2_T<T> v1, vector2_T<T> v2, T e)
{
	return cdm::nearly_equal(v1.x, v2.x, e) &&
	       cdm::nearly_equal(v1.y, v2.y, e);
}
template <arithmetic T>
constexpr vector2_T<T> element_wise_min(vector2_T<T> v0, vector2_T<T> v1)
{
	return {
	    std::min(v0.x, v1.x),
	    std::min(v0.y, v1.y),
	};
}
template <arithmetic T>
constexpr vector2_T<T> element_wise_max(vector2_T<T> v0, vector2_T<T> v1)
{
	return {
	    std::max(v0.x, v1.x),
	    std::max(v0.y, v1.y),
	};
}
template <arithmetic T>
constexpr vector2_T<T> element_wise_abs(vector2_T<T> v)
{
	return {
	    std::abs(v.x),
	    std::abs(v.y),
	};
}
template <arithmetic T>
constexpr vector2_T<T> element_wise_lerp(vector2_T<T> begin, vector2_T<T> end, vector2_T<T> percent)
{
	return {
		std::lerp(begin.x, end.x, percent.x),
		std::lerp(begin.y, end.y, percent.y),
	};
}
#pragma endregion

#pragma region definition_vector3_t
template <arithmetic T>
constexpr vector3_T<T>::vector3_T(const direction3_T<T>& d) : vector3_T(d.x(), d.y(), d.z()){}

template <arithmetic T>
vector2_T<T> vector3_T<T>::xy() const
{
	return {x, y};
}

template <arithmetic T>
vector3_T<T> vector3_T<T>::operator+(vector3_T<T> v) const
{
	return {
	    x + v.x,
	    y + v.y,
	    z + v.z,
	};
}
template <arithmetic T>
vector3_T<T> vector3_T<T>::operator-(vector3_T<T> v) const
{
	return {
	    x - v.x,
	    y - v.y,
	    z - v.z,
	};
}

template <arithmetic T>
vector3_T<T>& vector3_T<T>::operator+=(vector3_T<T> v)
{
	*this = *this + v;
	return *this;
}
template <arithmetic T>
vector3_T<T>& vector3_T<T>::operator-=(vector3_T<T> v)
{
	*this = *this - v;
	return *this;
}

template <arithmetic T>
vector3_T<T> vector3_T<T>::operator-() const
{
	return {
	    -x,
	    -y,
	    -z,
	};
}

template <arithmetic T>
bool vector3_T<T>::operator==(vector3_T<T> v) const
{
	return x == v.x && y == v.y && z == v.z;
}
template <arithmetic T>
bool vector3_T<T>::operator!=(vector3_T<T> v) const
{
	return !operator==(v);
}

template <arithmetic T>
vector3_T<T> operator*(T f, vector3_T<T> v)
{
	return v * f;
}
template <arithmetic T>
vector3_T<T> operator*(vector3_T<T> v, T f)
{
	return {
	    v.x * f,
	    v.y * f,
	    v.z * f,
	};
}
template <arithmetic T>
vector3_T<T> operator/(vector3_T<T> v, T f)
{
	return {
	    v.x / f,
	    v.y / f,
	    v.z / f,
	};
}
template <arithmetic T>
vector3_T<T>& operator*=(vector3_T<T>& v, T f)
{
	v = v * f;
	return v;
}
template <arithmetic T>
vector3_T<T>& operator/=(vector3_T<T>& v, T f)
{
	v = v / f;
	return v;
}

template <arithmetic T>
T norm(vector3_T<T> v)
{
	return std::sqrt(norm_squared(v));
}
template <arithmetic T>
T norm_squared(vector3_T<T> v)
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}
template <arithmetic T>
vector3_T<T> normalize(vector3_T<T> v)
{
	const T n = norm(v);
	v.x /= n;
	v.y /= n;
	v.z /= n;
	return v;
}
//template <arithmetic T>
//constexpr vector3_T<T> clamp(vector3_T<T> v,
//                             vector3_T<T> min,
//                             vector3_T<T> max)
//{
//	v.clamp(min, max);
//	return v;
//}
template <arithmetic T>
constexpr T dot(vector3_T<T> lhs, vector3_T<T> rhs)
{
	return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}
template <arithmetic T>
constexpr vector3_T<T> cross(vector3_T<T> lhs, vector3_T<T> rhs)
{
	return {
	    lhs.y * rhs.z - lhs.z * rhs.y,
	    lhs.z * rhs.x - lhs.x * rhs.z,
	    lhs.x * rhs.y - lhs.y * rhs.x,
	};
}
template <arithmetic T>
vector3_T<T> lerp(vector3_T<T> begin, vector3_T<T> end, T percent)
{
	return (end - begin) * percent + begin;
}
template <arithmetic T>
vector3_T<T> nlerp(vector3_T<T> begin, vector3_T<T> end, T percent)
{
	return normalize(lerp(begin, end, percent));
}
template <arithmetic T>
T distance_between(vector3_T<T> v1, vector3_T<T> v2)
{
	return norm(v1 - v2);
}
template <arithmetic T>
T distance_squared_between(vector3_T<T> v1, vector3_T<T> v2)
{
	return norm_squared(v1 - v2);
}
template <arithmetic T>
vector3_T<T> from_to(vector3_T<T> from, vector3_T<T> to)
{
	return {to.x - from.x, to.y - from.y, to.z - from.z};
}
template <arithmetic T>
radian_T<T> angle_between(vector3_T<T> v1, vector3_T<T> v2)
{
	T divisor = std::sqrt(norm_squared(v1) * norm_squared(v2));
	if (divisor == T(0))
		return radian_T<T>{T(0)};

	T alpha = dot(v1, v2) / divisor;
	return radian_T<T>(std::acos(cdm::clamp(alpha, T(-1), T(1))));
}
template <arithmetic T>
bool nearly_equal(vector3_T<T> v1, vector3_T<T> v2, T e)
{
	return nearly_equal(v1.x, v2.x, e) && nearly_equal(v1.y, v2.y, e) &&
	       nearly_equal(v1.z, v2.z, e);
}
template <arithmetic T>
constexpr vector3_T<T> element_wise_min(vector3_T<T> v0, vector3_T<T> v1)
{
	return {
	    std::min(v0.x, v1.x),
	    std::min(v0.y, v1.y),
	    std::min(v0.z, v1.z),
	};
}
template <arithmetic T>
constexpr vector3_T<T> element_wise_max(vector3_T<T> v0, vector3_T<T> v1)
{
	return {
	    std::max(v0.x, v1.x),
	    std::max(v0.y, v1.y),
	    std::max(v0.z, v1.z),
	};
}
template <arithmetic T>
constexpr vector3_T<T> element_wise_abs(vector3_T<T> v)
{
	return {
	    std::abs(v.x),
	    std::abs(v.y),
	    std::abs(v.z),
	};
}
template <arithmetic T>
constexpr vector3_T<T> element_wise_lerp(vector3_T<T> begin, vector3_T<T> end, vector3_T<T> percent)
{
	return {
		std::lerp(begin.x, end.x, percent.x),
		std::lerp(begin.y, end.y, percent.y),
		std::lerp(begin.z, end.z, percent.z),
	};
}

template <arithmetic T>
radian_T<T> angle_around_axis(vector3_T<T> v0,
                              vector3_T<T> v1,
                              direction3_T<T> axis)
{
	vector3_T<T> c = cross(v0, v1);
	T angle = atan2(norm(c), dot(v0, v1));
	return radian_T<T>{dot(c, axis) < T(0) ? -angle : angle};
}
#pragma endregion

#pragma region definition_vector4_t
template <arithmetic T>
vector2_T<T> vector4_T<T>::xy() const
{
	return {x, y};
}
template <arithmetic T>
vector3_T<T> vector4_T<T>::xyz() const
{
	return {x, y, z};
}

template <arithmetic T>
vector4_T<T> vector4_T<T>::operator+(vector4_T<T> v) const
{
	return {x + v.x, y + v.y, z + v.z, w + v.w};
}
template <arithmetic T>
vector4_T<T> vector4_T<T>::operator-(vector4_T<T> v) const
{
	return {x - v.x, y - v.y, z - v.z, w - v.w};
}

template <arithmetic T>
vector4_T<T>& vector4_T<T>::operator+=(vector4_T<T> v)
{
	*this = *this + v;
	return *this;
}
template <arithmetic T>
vector4_T<T>& vector4_T<T>::operator-=(vector4_T<T> v)
{
	*this = *this - v;
	return *this;
}

template <arithmetic T>
vector4_T<T> vector4_T<T>::operator-() const
{
	return {-x, -y, -z, -w};
}

template <arithmetic T>
bool vector4_T<T>::operator==(vector4_T<T> v) const
{
	return x == v.x && y == v.y && z == v.z && w == v.w;
}
template <arithmetic T>
bool vector4_T<T>::operator!=(vector4_T<T> v) const
{
	return !operator==(v);
}

template <arithmetic T>
vector4_T<T> operator*(T f, vector4_T<T> v)
{
	return v * f;
}
template <arithmetic T>
vector4_T<T> operator*(vector4_T<T> v, T f)
{
	return {v.x * f, v.y * f, v.z * f, v.w * f};
}
template <arithmetic T>
vector4_T<T> operator/(vector4_T<T> v, T f)
{
	return {v.x / f, v.y / f, v.z / f, v.w / f};
}
template <arithmetic T>
vector4_T<T>& operator*=(vector4_T<T>& v, T f)
{
	v = v * f;
	return v;
}
template <arithmetic T>
vector4_T<T>& operator/=(vector4_T<T>& v, T f)
{
	v = v / f;
	return v;
}

template <arithmetic T>
T norm(vector4_T<T> v)
{
	return v.norm();
}
template <arithmetic T>
T norm_squared(vector4_T<T> v)
{
	return v.norm_squared();
}
template <arithmetic T>
vector4_T<T> normalize(vector4_T<T> v)
{
	v.normalize();
	return v;
}
//template <arithmetic T>
//constexpr vector4_T<T> clamp(vector4_T<T> v,
//                             vector4_T<T> min,
//                             vector4_T<T> max)
//{
//	v.clamp(min, max);
//	return v;
//}
template <arithmetic T>
T dot(vector4_T<T> lhs, vector4_T<T> rhs)
{
	return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w;
}
template <arithmetic T>
vector4_T<T> lerp(vector4_T<T> begin, vector4_T<T> end, T percent)
{
	return (end - begin) * percent + begin;
}
template <arithmetic T>
vector4_T<T> nlerp(vector4_T<T> begin, vector4_T<T> end, T percent)
{
	return normalize(lerp(begin, end, percent));
}
template <arithmetic T>
T distance_between(vector4_T<T> v1, vector4_T<T> v2)
{
	return norm(v1 - v2);
}
template <arithmetic T>
T distance_squared_between(vector4_T<T> v1, vector4_T<T> v2)
{
	return norm_squared(v1 - v2);
}
template <arithmetic T>
vector4_T<T> from_to(vector4_T<T> from, vector4_T<T> to)
{
	return to - from;
}
template <arithmetic T>
bool nearly_equal(vector4_T<T> v1, vector4_T<T> v2, T e)
{
	return nearly_equal(v1.x, v2.x, e) && nearly_equal(v1.y, v2.y, e) &&
	       nearly_equal(v1.z, v2.z, e) && nearly_equal(v1.w, v2.w, e);
}
template <arithmetic T>
constexpr vector4_T<T> element_wise_min(vector4_T<T> v0, vector4_T<T> v1)
{
	return {
	    std::min(v0.x, v1.x),
	    std::min(v0.y, v1.y),
	    std::min(v0.z, v1.z),
	    std::min(v0.w, v1.w),
	};
}
template <arithmetic T>
constexpr vector4_T<T> element_wise_max(vector4_T<T> v0, vector4_T<T> v1)
{
	return {
	    std::max(v0.x, v1.x),
	    std::max(v0.y, v1.y),
	    std::max(v0.z, v1.z),
	    std::max(v0.w, v1.w),
	};
}
template <arithmetic T>
constexpr vector4_T<T> element_wise_abs(vector4_T<T> v)
{
	return {
	    std::abs(v.x),
	    std::abs(v.y),
	    std::abs(v.z),
	    std::abs(v.w),
	};
}
template <arithmetic T>
constexpr vector4_T<T> element_wise_lerp(vector4_T<T> begin, vector4_T<T> end, vector4_T<T> percent)
{
	return {
		std::lerp(begin.x, end.x, percent.x),
		std::lerp(begin.y, end.y, percent.y),
		std::lerp(begin.z, end.z, percent.z),
		std::lerp(begin.w, end.w, percent.w),
	};
}
#pragma endregion

#pragma region definition_normalized
template <normalizable T>
normalized<T>::normalized(const T& t) : vector(normalize(t))
{
	assert(vector.norm_squared() != 0);
}

template <normalizable T>
normalized<T> normalized<T>::already_normalized(const T& t)
{
	assert(t.norm_squared() != 0);
	normalized<T> res;
	res.vector = t;
	return res;
}
template <normalizable T>
normalized<T> normalized<T>::already_normalized(T&& t) noexcept
{
	assert(t.norm_squared() != 0);
	normalized<T> res;
	res.vector = std::move(t);
	return res;
}

template <normalizable T>
const T& normalized<T>::operator*() const
{
	return vector;
}
template <normalizable T>
const T* normalized<T>::operator->() const
{
	return &vector;
}
template <normalizable T>
normalized<T>::operator const T&() const
{
	return vector;
}

template <normalizable T>
normalized<T> normalized<T>::operator-() const
{
	return -vector;
}

template <normalizable T>
bool normalized<T>::operator==(const normalized<T>& v) const
{
	return vector == v.vector;
}
template <normalizable T>
bool normalized<T>::operator!=(const normalized<T>& v) const
{
	return vector != v.vector;
}
template <normalizable T>
bool normalized<T>::operator==(const T& v) const
{
	return vector == v;
}
template <normalizable T>
bool normalized<T>::operator!=(const T& v) const
{
	return vector != v;
}

template <arithmetic T>
vector3_T<T> cross(normalized<vector3_T<T>> lhs, normalized<vector3_T<T>> rhs)
{
	return cross(*lhs, *rhs);
}
template <arithmetic T>
vector3_T<T> cross(vector3_T<T> lhs, normalized<vector3_T<T>> rhs)
{
	return cross(lhs, *rhs);
}
template <arithmetic T>
vector3_T<T> cross(normalized<vector3_T<T>> lhs, vector3_T<T> rhs)
{
	return cross(*lhs, rhs);
}

template <arithmetic T>
T dot(normalized<vector3_T<T>> lhs, normalized<vector3_T<T>> rhs)
{
	return dot(*lhs, *rhs);
}
template <arithmetic T>
T dot(vector3_T<T> lhs, normalized<vector3_T<T>> rhs)
{
	return dot(lhs, *rhs);
}
template <arithmetic T>
T dot(normalized<vector3_T<T>> lhs, vector3_T<T> rhs)
{
	return dot(*lhs, rhs);
}
#pragma endregion

#pragma region definition_matrix2_t
template <arithmetic T>
matrix2_T<T>::matrix2_T(T e00, T e10, T e01, T e11)
    : m00{e00}, m10{e10}, m01{e01}, m11{e11}
{
}
template <arithmetic T>
matrix2_T<T>::matrix2_T(const std::array<T, 4>& a)
    : matrix2_T(a[0], a[1], a[2], a[3])
{
}
template <arithmetic T>
matrix2_T<T> matrix2_T<T>::zero()
{
	return {T(0), T(0), T(0), T(0)};
}
template <arithmetic T>
matrix2_T<T> matrix2_T<T>::identity()
{
	return {T(1), T(0), T(0), T(1)};
}
template <arithmetic T>
matrix2_T<T> matrix2_T<T>::rows(vector2_T<T> row0, vector2_T<T> row1)
{
	return {row0.x, row0.y, row1.x, row1.y};
}
template <arithmetic T>
matrix2_T<T> matrix2_T<T>::columns(vector2_T<T> column0, vector2_T<T> column1)
{
	return {column0.x, column1.x, column0.y, column1.y};
}
template <arithmetic T>
matrix2_T<T> matrix2_T<T>::rotation(radian_T<T> angle)
{
	T c = cos(angle);
	T s = sin(angle);
	return {c, -s, s, c};
}
template <arithmetic T>
matrix2_T<T> matrix2_T<T>::rotation(normalized<complex_T<T>> angle)
{
	return {angle->r, -angle->i, angle->i, angle->r};
}

template <arithmetic T>
vector2_T<T> matrix2_T<T>::operator*(vector2_T<T> v) const
{
	return {m00 * v.x + m01 * v.y, m10 * v.x + m11 * v.y};
}
template <arithmetic T>
matrix2_T<T> matrix2_T<T>::operator*(matrix2_T<T> m) const
{
	return {m00 * m.m00 + m01 * m.m10, m00 * m.m01 + m01 * m.m11,

	        m10 * m.m00 + m11 * m.m10, m10 * m.m01 + m11 * m.m11};
}

template <arithmetic T>
matrix2_T<T> transpose(matrix2_T<T> m)
{
	std::swap(m.m01, m.m10);
	return m;
}
#pragma endregion

#pragma region definition_matrix3_t
template <arithmetic T>
matrix3_T<T>::matrix3_T(T e00,
                        T e10,
                        T e20,
                        T e01,
                        T e11,
                        T e21,
                        T e02,
                        T e12,
                        T e22)
    : m00{e00},
      m10{e10},
      m20{e20},
      m01{e01},
      m11{e11},
      m21{e21},
      m02{e02},
      m12{e12},
      m22{e22}
{
}
template <arithmetic T>
matrix3_T<T>::matrix3_T(matrix2_T<T> m)
    : matrix3_T(m.m00, m.m10, T(0), m.m10, m.m11, T(0), T(0), T(0), T(1))
{
}
template <arithmetic T>
matrix3_T<T>::matrix3_T(const std::array<T, 9>& a)
    : matrix3_T(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8])
{
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::zero()
{
	return {T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0)};
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::identity()
{
	return {T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1)};
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::rows(vector3_T<T> row0,
                                vector3_T<T> row1,
                                vector3_T<T> row2)
{
	return {row0.x, row0.y, row0.z, row1.x, row1.y,
	        row1.z, row2.x, row2.y, row2.z};
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::columns(vector3_T<T> column0,
                                   vector3_T<T> column1,
                                   vector3_T<T> column2)
{
	return {column0.x, column1.x, column2.x, column0.y, column1.y,
	        column2.y, column0.z, column1.z, column2.z};
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::rotation(direction3_T<T> axis,
                                    T sinAngle,
                                    T cosAngle)
{
	const T l{axis.x()};
	const T m{axis.y()};
	const T n{axis.z()};
	const T s{sinAngle};
	const T c{cosAngle};
	const T d{T(1) - cosAngle};
	const T ls{l * s};
	const T ms{m * s};
	const T ns{n * s};
	const T ld{l * d};
	const T md{m * d};
	const T nd{n * d};
	return {
	    (l * ld) + c,  (m * ld) - ns, (n * ld) + ms,  //
	    (l * md) + ns, (m * md) + c,  (n * md) - ls,  //
	    (l * nd) - ms, (m * nd) + ls, (n * nd) + c,   //
	};
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::rotation(direction3_T<T> axis, radian_T<T> angle)
{
	return matrix3_T<T>::rotation(axis, sin(angle), cos(angle));
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::rotation(direction3_T<T> axis,
                                    normalized<complex_T<T>> angle)
{
	return matrix3_T<T>::rotation(axis, sin(angle), cos(angle));
}
template <arithmetic T>
template <std::signed_integral U, U NumeratorT, U DenominatorT>
matrix3_T<T> matrix3_T<T>::rotation(
    direction3_T<T> axis,
    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle)
{
	return matrix3_T<T>::rotation(axis, sin<T>(angle), cos<T>(angle));
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::rotation(euler_angles_T<T> r)
{
	matrix3_T<T> RX = matrix3_T<T>::identity();
	RX.m11 = cos(-r.x);
	RX.m12 = -sin(-r.x);
	RX.m21 = sin(-r.x);
	RX.m22 = cos(-r.x);

	matrix3_T<T> RY = matrix3_T<T>::identity();
	RY.m00 = cos(-r.y);
	RY.m02 = sin(-r.y);
	RY.m20 = -sin(-r.y);
	RY.m22 = cos(-r.y);

	matrix3_T<T> RZ = matrix3_T<T>::identity();
	RZ.m00 = cos(-r.z);
	RZ.m01 = -sin(-r.z);
	RZ.m10 = sin(-r.z);
	RZ.m11 = cos(-r.z);

	return RY * RX * RZ;
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::rotation(quaternion_T<T> q)
{
	return {
	    detail::get_quaternion_t_matrix_element<0, 0>(q),
	    detail::get_quaternion_t_matrix_element<1, 0>(q),
	    detail::get_quaternion_t_matrix_element<2, 0>(q),
	    detail::get_quaternion_t_matrix_element<0, 1>(q),
	    detail::get_quaternion_t_matrix_element<1, 1>(q),
	    detail::get_quaternion_t_matrix_element<2, 1>(q),
	    detail::get_quaternion_t_matrix_element<0, 2>(q),
	    detail::get_quaternion_t_matrix_element<1, 2>(q),
	    detail::get_quaternion_t_matrix_element<2, 2>(q),
	};
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::rotation_around_x(T sinAngle, T cosAngle)
{
	return {
	    T(1), T(0),     T(0),       //
	    T(0), cosAngle, -sinAngle,  //
	    T(0), sinAngle, cosAngle    //
	};
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::rotation_around_y(T sinAngle, T cosAngle)
{
	return {
	    cosAngle,  T(0), sinAngle,  //
	    T(0),      T(1), T(0),      //
	    -sinAngle, T(0), cosAngle   //
	};
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::rotation_around_z(T sinAngle, T cosAngle)
{
	return {
	    cosAngle, -sinAngle, T(0),  //
	    sinAngle, cosAngle,  T(0),  //
	    T(0),     T(0),      T(1)   //
	};
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::rotation_around_x(radian_T<T> angle)
{
	return rotation_around_x(sin(angle), cos(angle));
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::rotation_around_y(radian_T<T> angle)
{
	return rotation_around_y(sin(angle), cos(angle));
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::rotation_around_z(radian_T<T> angle)
{
	return rotation_around_z(sin(angle), cos(angle));
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::rotation_around_x(normalized<complex_T<T>> angle)
{
	return rotation_around_x(angle->i, angle->r);
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::rotation_around_y(normalized<complex_T<T>> angle)
{
	return rotation_around_y(angle->i, angle->r);
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::rotation_around_z(normalized<complex_T<T>> angle)
{
	return rotation_around_z(angle->i, angle->r);
}
template <arithmetic T>
template <std::signed_integral U, U NumeratorT, U DenominatorT>
matrix3_T<T> matrix3_T<T>::rotation_around_x(
    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle)
{
	return rotation_around_x(sin<T>(angle), cos<T>(angle));
}
template <arithmetic T>
template <std::signed_integral U, U NumeratorT, U DenominatorT>
matrix3_T<T> matrix3_T<T>::rotation_around_y(
    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle)
{
	return rotation_around_y(sin<T>(angle), cos<T>(angle));
}
template <arithmetic T>
template <std::signed_integral U, U NumeratorT, U DenominatorT>
matrix3_T<T> matrix3_T<T>::rotation_around_z(
    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle)
{
	return rotation_around_z(sin<T>(angle), cos<T>(angle));
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::scale(vector3_T<T> s)
{
	return scale(s.x, s.y, s.z);
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::scale(T x, T y, T z)
{
	return {
	    x,    T(0), T(0),  //
	    T(0), y,    T(0),  //
	    T(0), T(0), z,     //
	};
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::scale(T s)
{
	return scale(s, s, s);
}
template <arithmetic T>
matrix3_T<T>& matrix3_T<T>::inverse()
{
	if (is_orthogonal())
	{
		return transpose();
	}

	const T det = determinant();
	assert(std::abs(det) > std::numeric_limits<T>::epsilon());

	const T invDet = T(1) / det;
	const T recM00 = m00;
	const T recM10 = m10;
	const T recM20 = m20;
	const T recM01 = m01;
	const T recM11 = m11;
	const T recM02 = m02;
	m00 = (recM11 * m22 - m21 * m12) * invDet;
	m10 = (m12 * recM20 - recM10 * m22) * invDet;
	m20 = (recM10 * m12 - recM20 * recM11) * invDet;
	m01 = (recM02 * m21 - recM01 * m22) * invDet;
	m11 = (recM00 * m22 - recM20 * recM02) * invDet;
	m21 = (recM01 * recM20 - recM00 * m21) * invDet;
	m02 = (recM01 * m12 - recM11 * recM02) * invDet;
	m12 = (recM02 * recM10 - recM00 * m12) * invDet;
	m22 = (recM00 * recM11 - recM10 * recM01) * invDet;

	return *this;
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::get_inversed() const
{
	matrix3_T<T> res = *this;
	res.inverse();
	return res;
}
template <arithmetic T>
matrix3_T<T>& matrix3_T<T>::transpose()
{
	std::swap(m01, m10);
	std::swap(m02, m20);
	std::swap(m12, m21);
	return *this;
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::get_transposed() const
{
	return {
	    m00, m01, m02,  //
	    m10, m11, m12,  //
	    m20, m21, m22   //
	};
}

template <arithmetic T>
T matrix3_T<T>::determinant() const
{
	return m00 * m11 * m22 + m01 * m12 * m20 + m02 * m10 * m21 -
	       m20 * m11 * m02 - m21 * m12 * m00 - m22 * m10 * m01;
}
template <arithmetic T>
bool matrix3_T<T>::is_orthogonal() const
{
	return nearly_equal(m00 * m01 + m10 * m11 + m20 * m21, T(0)) &&
	       nearly_equal(m00 * m02 + m10 * m12 + m20 * m22, T(0)) &&
	       nearly_equal(m01 * m02 + m11 * m12 + m21 * m22, T(0)) &&
	       nearly_equal(m00 * m00 + m10 * m10 + m20 * m20, T(1)) &&
	       nearly_equal(m01 * m01 + m11 * m11 + m21 * m21, T(1)) &&
	       nearly_equal(m02 * m02 + m12 * m12 + m22 * m22, T(1));
}

template <arithmetic T>
euler_angles_T<T> matrix3_T<T>::to_euler_angles() const
{
	const auto q = cdm::quaternion_T<T>::from_rotation_matrix(*this);
	return q.to_euler_angles();
}

template <arithmetic T>
matrix3_T<T> matrix3_T<T>::operator*(T f) const
{
	return {
	    m00 * f, m10 * f, m20 * f,  //
	    m01 * f, m11 * f, m21 * f,  //
	    m02 * f, m12 * f, m22 * f   //
	};
}
template <arithmetic T>
vector3_T<T> matrix3_T<T>::operator*(vector3_T<T> v) const
{
	return {m00 * v.x + m10 * v.y + m20 * v.z,
	        m01 * v.x + m11 * v.y + m21 * v.z,
	        m02 * v.x + m12 * v.y + m22 * v.z};
}
template <arithmetic T>
matrix3_T<T> matrix3_T<T>::operator*(const matrix3_T<T>& m) const
{
	return {m00 * m.m00 + m10 * m.m01 + m20 * m.m02,
	        m00 * m.m10 + m10 * m.m11 + m20 * m.m12,
	        m00 * m.m20 + m10 * m.m21 + m20 * m.m22,

	        m01 * m.m00 + m11 * m.m01 + m21 * m.m02,
	        m01 * m.m10 + m11 * m.m11 + m21 * m.m12,
	        m01 * m.m20 + m11 * m.m21 + m21 * m.m22,

	        m02 * m.m00 + m12 * m.m01 + m22 * m.m02,
	        m02 * m.m10 + m12 * m.m11 + m22 * m.m12,
	        m02 * m.m20 + m12 * m.m21 + m22 * m.m22};
}
#pragma endregion

#pragma region definition_matrix4_t
template <arithmetic T>
matrix4_T<T>::matrix4_T(T e00,
                        T e10,
                        T e20,
                        T e30,
                        T e01,
                        T e11,
                        T e21,
                        T e31,
                        T e02,
                        T e12,
                        T e22,
                        T e32,
                        T e03,
                        T e13,
                        T e23,
                        T e33)
    : m00{e00},
      m10{e10},
      m20{e20},
      m30{e30},
      m01{e01},
      m11{e11},
      m21{e21},
      m31{e31},
      m02{e02},
      m12{e12},
      m22{e22},
      m32{e32},
      m03{e03},
      m13{e13},
      m23{e23},
      m33{e33}
{
}
template <arithmetic T>
matrix4_T<T>::matrix4_T(matrix2_T<T> m)
    : matrix4_T(m.m00,
                m.m10,
                T(0),
                T(0),
                m.m10,
                m.m11,
                T(0),
                T(0),
                T(0),
                T(0),
                T(1),
                T(0),
                T(0),
                T(0),
                T(0),
                T(1))
{
}
template <arithmetic T>
matrix4_T<T>::matrix4_T(const matrix3_T<T>& m)
    : matrix4_T(m.m00,
                m.m10,
                m.m20,
                T(0),
                m.m01,
                m.m11,
                m.m21,
                T(0),
                m.m02,
                m.m12,
                m.m22,
                T(0),
                T(0),
                T(0),
                T(0),
                T(1))
{
}
template <arithmetic T>
matrix4_T<T>::matrix4_T(const perspective_T<T>& p) : matrix4_T(p.to_matrix4())
{
}
template <arithmetic T>
matrix4_T<T>::matrix4_T(const transform3_T<T>& t) : matrix4_T(t.to_matrix4())
{
}
template <arithmetic T>
matrix4_T<T>::matrix4_T(const uniform_transform3_T<T>& t)
    : matrix4_T(t.to_matrix4())
{
}
template <arithmetic T>
matrix4_T<T>::matrix4_T(const unscaled_transform3_T<T>& t)
    : matrix4_T(t.to_matrix4())
{
}
template <arithmetic T>
matrix4_T<T>::matrix4_T(const std::array<T, 16>& a)
    : matrix4_T(a[0],
                a[1],
                a[2],
                a[3],
                a[4],
                a[5],
                a[6],
                a[7],
                a[8],
                a[9],
                a[10],
                a[11],
                a[12],
                a[13],
                a[14],
                a[15])
{
}

template <arithmetic T>
matrix4_T<T> matrix4_T<T>::zero()
{
	return {T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0),
	        T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0)};
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::identity()
{
	return {T(1), T(0), T(0), T(0), T(0), T(1), T(0), T(0),
	        T(0), T(0), T(1), T(0), T(0), T(0), T(0), T(1)};
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::rows(vector4_T<T> row0,
                                vector4_T<T> row1,
                                vector4_T<T> row2,
                                vector4_T<T> row3)
{
	return {
	    row0.x, row0.y, row0.z, row0.w,  //
	    row1.x, row1.y, row1.z, row1.w,  //
	    row2.x, row2.y, row2.z, row2.w,  //
	    row3.x, row3.y, row3.z, row3.w   //
	};
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::columns(vector4_T<T> column0,
                                   vector4_T<T> column1,
                                   vector4_T<T> column2,
                                   vector4_T<T> column3)
{
	return {
	    column0.x, column1.x, column2.x, column3.x,  //
	    column0.y, column1.y, column2.y, column3.y,  //
	    column0.z, column1.z, column2.z, column3.z,  //
	    column0.w, column1.w, column2.w, column3.w   //
	};
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::rotation(direction3_T<T> axis, radian_T<T> angle)
{
	return matrix4_T<T>(matrix3_T<T>::rotation(axis, angle));
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::rotation(direction3_T<T> axis,
                                    normalized<complex_T<T>> angle)
{
	return matrix4_T<T>(matrix3_T<T>::rotation(axis, angle));
}
template <arithmetic T>
template <std::signed_integral U, U NumeratorT, U DenominatorT>
matrix4_T<T> matrix4_T<T>::rotation(
    direction3_T<T> axis,
    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle)
{
	return matrix4_T<T>(matrix3_T<T>::rotation(axis, angle));
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::rotation(euler_angles_T<T> r)
{
	return matrix4_T<T>(matrix3_T<T>::rotation(r));
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::rotation(quaternion_T<T> q)
{
	return matrix4_T<T>(matrix3_T<T>::rotation(q));
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::translation(vector3_T<T> t)
{
	return {
	    T(1), T(0), T(0), t.x,  //
	    T(0), T(1), T(0), t.y,  //
	    T(0), T(0), T(1), t.z,  //
	    T(0), T(0), T(0), T(1)  //
	};
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::translation(T x, T y, T z)
{
	return translation({x, y, z});
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::scale(vector3_T<T> s)
{
	return matrix4_T<T>(matrix3_T<T>::scale(s));
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::scale(T x, T y, T z)
{
	return matrix4_T<T>(matrix3_T<T>::scale(x, y, z));
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::scale(T s)
{
	return matrix4_T<T>(matrix3_T<T>::scale(s));
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::look_at(vector3_T<T> from,
                                   vector3_T<T> to,
                                   direction3_T<T> up)
{
	direction3_T<T> forward = (from - to);
	direction3_T<T> right = cross(up, forward);
	direction3_T<T> true_up = cross(forward, right);

	return {right.x(), true_up.x(), forward.x(), from.x(),
	        right.y(), true_up.y(), forward.y(), from.y(),
	        right.z(), true_up.z(), forward.z(), from.z(),
	        T(0),      T(0),        T(0),        T(1)};
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::orthographic(T left,
                                        T right,
                                        T top,
                                        T bottom,
                                        T near,
                                        T far)
{
	T a = T(1) / (right - left);
	T b = T(1) / (bottom - top);
	T c = T(1) / (near - far);

	return {
	    T(2) * a, T(0),     T(0), -(right + left) * a,
	    T(0),     T(2) * b, T(0), -(bottom + top) * b,
	    T(0),     T(0),     c,    near * c,
	    T(0),     T(0),     T(0), T(1),
	};
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::rotation_around_x(radian_T<T> angle)
{
	return matrix4_T<T>(matrix3_T<T>::rotation_around_x(angle));
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::rotation_around_y(radian_T<T> angle)
{
	return matrix4_T<T>(matrix3_T<T>::rotation_around_y(angle));
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::rotation_around_z(radian_T<T> angle)
{
	return matrix4_T<T>(matrix3_T<T>::rotation_around_z(angle));
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::rotation_around_x(normalized<complex_T<T>> angle)
{
	return matrix4_T<T>(matrix3_T<T>::rotation_around_x(angle));
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::rotation_around_y(normalized<complex_T<T>> angle)
{
	return matrix4_T<T>(matrix3_T<T>::rotation_around_y(angle));
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::rotation_around_z(normalized<complex_T<T>> angle)
{
	return matrix4_T<T>(matrix3_T<T>::rotation_around_z(angle));
}
template <arithmetic T>
template <std::signed_integral U, U NumeratorT, U DenominatorT>
matrix4_T<T> matrix4_T<T>::rotation_around_x(
    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle)
{
	return matrix4_T<T>(matrix3_T<T>::rotation_around_x(angle));
}
template <arithmetic T>
template <std::signed_integral U, U NumeratorT, U DenominatorT>
matrix4_T<T> matrix4_T<T>::rotation_around_y(
    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle)
{
	return matrix4_T<T>(matrix3_T<T>::rotation_around_y(angle));
}
template <arithmetic T>
template <std::signed_integral U, U NumeratorT, U DenominatorT>
matrix4_T<T> matrix4_T<T>::rotation_around_z(
    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle)
{
	return matrix4_T<T>(matrix3_T<T>::rotation_around_z(angle));
}
template <arithmetic T>
bool matrix4_T<T>::is_orthogonal() const
{
	return nearly_equal(m00 * m01 + m10 * m11 + m20 * m21 + m30 * m31, T(0)) &&
	       nearly_equal(m00 * m02 + m10 * m12 + m20 * m22 + m30 * m32, T(0)) &&
	       nearly_equal(m00 * m03 + m10 * m13 + m20 * m23 + m30 * m33, T(0)) &&
	       nearly_equal(m01 * m02 + m11 * m12 + m21 * m22 + m31 * m32, T(0)) &&
	       nearly_equal(m01 * m03 + m11 * m13 + m21 * m23 + m31 * m33, T(0)) &&
	       nearly_equal(m02 * m03 + m12 * m13 + m22 * m23 + m32 * m33, T(0)) &&
	       nearly_equal(m00 * m00 + m10 * m10 + m20 * m20 + m30 * m30, T(1)) &&
	       nearly_equal(m01 * m01 + m11 * m11 + m21 * m21 + m31 * m31, T(1)) &&
	       nearly_equal(m02 * m02 + m12 * m12 + m22 * m22 + m32 * m32, T(1)) &&
	       nearly_equal(m03 * m03 + m13 * m13 + m23 * m23 + m33 * m33, T(1));
}
template <arithmetic T>
bool matrix4_T<T>::is_homogenous() const
{
	return nearly_equal(m03, T(0)) &&  //
	       nearly_equal(m13, T(0)) &&  //
	       nearly_equal(m23, T(0)) &&  //
	       nearly_equal(m30, T(0)) &&  //
	       nearly_equal(m31, T(0)) &&  //
	       nearly_equal(m32, T(0)) &&  //
	       nearly_equal(m33, T(1));
}
template <arithmetic T>
T matrix4_T<T>::determinant() const
{
	if (is_homogenous())
	{
		return m00 * m11 * m22 + m01 * m12 * m20 + m02 * m10 * m21 -
		       m20 * m11 * m02 - m21 * m12 * m00 - m22 * m10 * m01;
	}

	const T det1 = m11 * (m22 * m33 - m32 * m23) -
	               m21 * (m12 * m33 - m32 * m13) +
	               m31 * (m12 * m23 - m22 * m13);
	const T det2 = m01 * (m22 * m33 - m32 * m23) -
	               m21 * (m02 * m33 - m32 * m03) +
	               m31 * (m02 * m23 - m22 * m03);
	const T det3 = m01 * (m12 * m33 - m32 * m13) -
	               m11 * (m02 * m33 - m32 * m03) +
	               m31 * (m02 * m13 - m12 * m03);
	const T det4 = m01 * (m12 * m23 - m22 * m13) -
	               m11 * (m02 * m23 - m22 * m03) +
	               m21 * (m02 * m13 - m12 * m03);
	return m00 * det1 - m10 * det2 + m20 * det3 - m30 * det4;
}

template <arithmetic T>
matrix3_T<T> matrix4_T<T>::extractRotationMatrix() const
{
	const T sx = column(0).xyz().norm();
	const T sy = column(1).xyz().norm();
	const T sz = column(2).xyz().norm();

	return matrix3_T<T>::rows(
		vector3_T<T>(row(0).column(0) / sx,
		             row(0).column(1) / sy,
		             row(0).column(2) / sz),
		vector3_T<T>(row(1).column(0) / sx,
                     row(1).column(1) / sy,
                     row(1).column(2) / sz),
		vector3_T<T>(row(2).column(0) / sx,
                     row(2).column(1) / sy,
                     row(2).column(2) / sz)
	);
}

template <arithmetic T>
vector3_T<T> matrix4_T<T>::transform_position(vector3_T<T> pos) const
{
	return (*this * vector4_T<T>(pos, T(1))).xyz();
}
template <arithmetic T>
vector3_T<T> matrix4_T<T>::transform_direction(vector3_T<T> dir) const
{
	return (*this * vector4_T<T>(dir, T(0))).xyz();
}
template <arithmetic T>
direction3_T<T> matrix4_T<T>::transform_direction(direction3_T<T> dir) const
{
	return direction3_T<T>((*this * vector4_T<T>(*dir, T(0))).xyz());
}
template <arithmetic T>
vector3_T<T> matrix4_T<T>::transform_position_perspective(vector3_T<T> pos) const
{
	const vector4_T<T> tmp = *this * vector4_T<T>(pos, T(1));
	return {
		tmp.x / tmp.w,
		tmp.y / tmp.w,
		tmp.z / tmp.w,
	};
}

template <arithmetic T>
matrix4_T<T> matrix4_T<T>::operator*(T f) const
{
	return {
	    m00 * f, m10 * f, m20 * f, m30 * f,  //
	    m00 * f, m10 * f, m20 * f, m30 * f,  //
	    m00 * f, m10 * f, m20 * f, m30 * f,  //
	    m00 * f, m10 * f, m20 * f, m30 * f   //
	};
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::operator/(T f) const
{
	return {
	    m00 / f, m10 / f, m20 / f, m30 / f,  //
	    m00 / f, m10 / f, m20 / f, m30 / f,  //
	    m00 / f, m10 / f, m20 / f, m30 / f,  //
	    m00 / f, m10 / f, m20 / f, m30 / f   //
	};
}
template <arithmetic T>
vector4_T<T> matrix4_T<T>::operator*(vector4_T<T> v) const
{
	return {
	    m00 * v.x + m10 * v.y + m20 * v.z + m30 * v.w,  //
	    m01 * v.x + m11 * v.y + m21 * v.z + m31 * v.w,  //
	    m02 * v.x + m12 * v.y + m22 * v.z + m32 * v.w,  //
	    m03 * v.x + m13 * v.y + m23 * v.z + m33 * v.w   //
	};
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::operator*(const matrix4_T<T>& m) const
{
	return {m.m00 * m00 + m.m01 * m10 + m.m02 * m20 + m.m03 * m30,
	        m.m10 * m00 + m.m11 * m10 + m.m12 * m20 + m.m13 * m30,
	        m.m20 * m00 + m.m21 * m10 + m.m22 * m20 + m.m23 * m30,
	        m.m30 * m00 + m.m31 * m10 + m.m32 * m20 + m.m33 * m30,

	        m.m00 * m01 + m.m01 * m11 + m.m02 * m21 + m.m03 * m31,
	        m.m10 * m01 + m.m11 * m11 + m.m12 * m21 + m.m13 * m31,
	        m.m20 * m01 + m.m21 * m11 + m.m22 * m21 + m.m23 * m31,
	        m.m30 * m01 + m.m31 * m11 + m.m32 * m21 + m.m33 * m31,

	        m.m00 * m02 + m.m01 * m12 + m.m02 * m22 + m.m03 * m32,
	        m.m10 * m02 + m.m11 * m12 + m.m12 * m22 + m.m13 * m32,
	        m.m20 * m02 + m.m21 * m12 + m.m22 * m22 + m.m23 * m32,
	        m.m30 * m02 + m.m31 * m12 + m.m32 * m22 + m.m33 * m32,

	        m.m00 * m03 + m.m01 * m13 + m.m02 * m23 + m.m03 * m33,
	        m.m10 * m03 + m.m11 * m13 + m.m12 * m23 + m.m13 * m33,
	        m.m20 * m03 + m.m21 * m13 + m.m22 * m23 + m.m23 * m33,
	        m.m30 * m03 + m.m31 * m13 + m.m32 * m23 + m.m33 * m33};
}

template <arithmetic T>
matrix4_T<T>& matrix4_T<T>::inverse()
{
	if (is_orthogonal())
	{
		return transpose();
	}

	const T det = determinant();
	assert(std::abs(det) > std::numeric_limits<T>::epsilon());

	const T invDet = T(1) / det;
	const T recM00 = m00;
	const T recM10 = m10;
	const T recM20 = m20;
	const T recM01 = m01;
	const T recM11 = m11;
	const T recM02 = m02;

	if (is_homogenous())
	{
		m00 = invDet * (recM11 * m22 - m21 * m12);
		m10 = invDet * (m12 * recM20 - recM10 * m22);
		m20 = invDet * (recM10 * m21 - recM20 * recM11);
		m01 = invDet * (recM02 * m21 - recM01 * m22);
		m11 = invDet * (recM00 * m22 - recM20 * recM02);
		m21 = invDet * (recM01 * recM20 - recM00 * m21);
		m02 = invDet * (recM01 * m12 - recM11 * recM02);
		m12 = invDet * (recM02 * recM10 - recM00 * m12);
		m22 = invDet * (recM00 * recM11 - recM10 * recM01);
		return *this;
	}

	const T recM30 = m30;
	const T recM21 = m21;
	const T recM31 = m31;
	const T recM12 = m12;
	const T recM22 = m22;
	const T recM03 = m03;
	const T recM13 = m13;

	m00 = invDet * (recM11 * recM22 * m33 - recM11 * m23 * m32 -
	                recM21 * recM12 * m33 + recM21 * recM13 * m32 +
	                recM31 * recM12 * m23 - recM31 * recM13 * recM22);

	m10 = invDet * (-recM10 * recM22 * m33 + recM10 * m23 * m32 +
	                recM20 * recM12 * m33 - recM20 * recM13 * m32 -
	                recM30 * recM12 * m23 + recM30 * recM13 * recM22);

	m20 = invDet * (recM10 * recM21 * m33 - recM10 * m23 * recM31 -
	                recM20 * recM11 * m33 + recM20 * recM13 * recM31 +
	                recM30 * recM11 * m23 - recM30 * recM13 * recM21);

	m30 = invDet * (-recM10 * recM21 * m32 + recM10 * recM22 * recM31 +
	                recM20 * recM11 * m32 - recM20 * recM12 * recM31 -
	                recM30 * recM11 * recM22 + recM30 * recM12 * recM21);

	m01 = invDet * (-recM01 * recM22 * m33 + recM01 * m23 * m32 +
	                recM21 * recM02 * m33 - recM21 * recM03 * m32 -
	                recM31 * recM02 * m23 + recM31 * recM03 * recM22);

	m11 = invDet * (recM00 * recM22 * m33 - recM00 * m23 * m32 -
	                recM20 * recM02 * m33 + recM20 * recM03 * m32 +
	                recM30 * recM02 * m23 - recM30 * recM03 * recM22);

	m21 = invDet * (-recM00 * recM21 * m33 + recM00 * m23 * recM31 +
	                recM20 * recM01 * m33 - recM20 * recM03 * recM31 -
	                recM30 * recM01 * m23 + recM30 * recM03 * recM21);

	m31 = invDet * (recM00 * recM21 * m32 - recM00 * recM22 * recM31 -
	                recM20 * recM01 * m32 + recM20 * recM02 * recM31 +
	                recM30 * recM01 * recM22 - recM30 * recM02 * recM21);

	m02 = invDet * (recM01 * recM12 * m33 - recM01 * recM13 * m32 -
	                recM11 * recM02 * m33 + recM11 * recM03 * m32 +
	                recM31 * recM02 * recM13 - recM31 * recM03 * recM12);

	m12 = invDet * (-recM00 * recM12 * m33 + recM00 * recM13 * m32 +
	                recM10 * recM02 * m33 - recM10 * recM03 * m32 -
	                recM30 * recM02 * recM13 + recM30 * recM03 * recM12);

	m22 = invDet * (recM00 * recM11 * m33 - recM00 * recM13 * recM31 -
	                recM10 * recM01 * m33 + recM10 * recM03 * recM31 +
	                recM30 * recM01 * recM13 - recM30 * recM03 * recM11);

	m32 = invDet * (-recM00 * recM11 * m32 + recM00 * recM12 * recM31 +
	                recM10 * recM01 * m32 - recM10 * recM02 * recM31 -
	                recM30 * recM01 * recM12 + recM30 * recM02 * recM11);

	m03 = invDet * (-recM01 * recM12 * m23 + recM01 * recM13 * recM22 +
	                recM11 * recM02 * m23 - recM11 * recM03 * recM22 -
	                recM21 * recM02 * recM13 + recM21 * recM03 * recM12);

	m13 = invDet * (recM00 * recM12 * m23 - recM00 * recM13 * recM22 -
	                recM10 * recM02 * m23 + recM10 * recM03 * recM22 +
	                recM20 * recM02 * recM13 - recM20 * recM03 * recM12);

	m23 = invDet * (-recM00 * recM11 * m23 + recM00 * recM13 * recM21 +
	                recM10 * recM01 * m23 - recM10 * recM03 * recM21 -
	                recM20 * recM01 * recM13 + recM20 * recM03 * recM11);

	m33 = invDet * (recM00 * recM11 * recM22 - recM00 * recM12 * recM21 -
	                recM10 * recM01 * recM22 + recM10 * recM02 * recM21 +
	                recM20 * recM01 * recM12 - recM20 * recM02 * recM11);

	return *this;
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::get_inversed() const
{
	matrix4_T<T> m{*this};
	m.inverse();
	return m;
}
template <arithmetic T>
matrix4_T<T>& matrix4_T<T>::transpose()
{
	std::swap(m01, m10);
	std::swap(m02, m20);
	std::swap(m03, m30);
	std::swap(m12, m21);
	std::swap(m13, m31);
	std::swap(m23, m32);
	return *this;
}
template <arithmetic T>
matrix4_T<T> matrix4_T<T>::get_transposed() const
{
	matrix4_T<T> m{*this};
	m.transpose();
	return m;
}

template <arithmetic T>
matrix4_T<T> operator*(const transform3_T<T>& t, const matrix4_T<T>& m)
{
	return t.to_matrix4() * m;
}
template <arithmetic T>
matrix4_T<T> operator*(const matrix4_T<T>& m, const transform3_T<T>& t)
{
	return m * t.to_matrix4();
}
template <arithmetic T>
matrix4_T<T> operator*(const uniform_transform3_T<T>& t, const matrix4_T<T>& m)
{
	return t.to_matrix4() * m;
}
template <arithmetic T>
matrix4_T<T> operator*(const matrix4_T<T>& m, const uniform_transform3_T<T>& t)
{
	return m * t.to_matrix4();
}
template <arithmetic T>
matrix4_T<T> operator*(const unscaled_transform3_T<T>& t,
                       const matrix4_T<T>& m)
{
	return t.to_matrix4() * m;
}
template <arithmetic T>
matrix4_T<T> operator*(const matrix4_T<T>& m,
                       const unscaled_transform3_T<T>& t)
{
	return m * t.to_matrix4();
}
#pragma endregion

#pragma region definition_perspective_t
template <arithmetic T>
perspective_T<T>::perspective_T(radian_T<T> angle, T ratio, T near, T far)
    : m_angle{angle},
      m_ratio{ratio},
      m_invRatio{T(1) / m_ratio},
      m_near{near},
      m_far{far},
      m_focalLength{T(1) / tan(m_angle / T(2))}
{
}
template <arithmetic T>
template <std::signed_integral U, U NumeratorT, U DenominatorT>
perspective_T<T>::perspective_T(
    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle,
    T ratio,
    T near,
    T far)
    : perspective_T(angle, ratio, near, far),
      m_focalLength{
          T(1) /
          tan<T>(static_pi_fraction_T<U, NumeratorT, DenominatorT * U(2)>{})}
{
}
template <arithmetic T>
void perspective_T<T>::set(radian_T<T> angle, T ratio, T near, T far)
{
	m_angle = angle;
	m_ratio = ratio;
	m_invRatio = T(1) / m_ratio;
	m_near = near;
	m_far = far;
	m_focalLength = T(1) / tan(m_angle / T(2));
}
template <arithmetic T>
template <std::signed_integral U, U NumeratorT, U DenominatorT>
void perspective_T<T>::set(
    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle,
    T ratio,
    T near,
    T far)
{
	set(angle, ratio, near, far);
	m_focalLength =
	    T(1) /
	    tan<T>(static_pi_fraction_T<U, NumeratorT, DenominatorT * U(2)>{});
}
template <arithmetic T>
void perspective_T<T>::set_angle(radian_T<T> angle)
{
	m_angle = angle;
	m_focalLength = T(1) / tan(angle / T(2));
}
template <arithmetic T>
template <std::signed_integral U, U NumeratorT, U DenominatorT>
void perspective_T<T>::set_angle(
    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle)
{
	m_angle = angle;
	m_focalLength =
	    T(1) /
	    tan<T>(static_pi_fraction_T<U, NumeratorT, DenominatorT * U(2)>{});
}
template <arithmetic T>
radian_T<T> perspective_T<T>::get_angle() const
{
	return m_angle;
}
template <arithmetic T>
void perspective_T<T>::set_ratio(T ratio)
{
	m_ratio = ratio;
	m_invRatio = T(1) / m_ratio;
}
template <arithmetic T>
T perspective_T<T>::get_ratio() const
{
	return m_ratio;
}
template <arithmetic T>
void perspective_T<T>::set_near(T near_plane_t)
{
	m_near = near_plane_t;
}
template <arithmetic T>
T perspective_T<T>::get_near() const
{
	return m_near;
}
template <arithmetic T>
void perspective_T<T>::set_far(T far_plane_t)
{
	m_far = far_plane_t;
}
template <arithmetic T>
T perspective_T<T>::get_far() const
{
	return m_far;
}
template <arithmetic T>
matrix4_T<T> perspective_T<T>::to_matrix4() const
{
	matrix4_T<T> res{matrix4_T<T>::zero()};
	const T a = m_focalLength * m_invRatio;
	const T b = m_focalLength;
	const T farMinusNear = m_far - m_near;
	const T invFarMinusNear = T(1) / farMinusNear;
	const T c = m_near * invFarMinusNear;
	const T d = m_far * c;

	res.column(0).row(0) = a;
	res.column(1).row(1) = b;
	res.column(2).row(2) = c;
	res.column(2).row(3) = T(-1);
	res.column(3).row(2) = d;

	//      0   1   2   3
	//   +-----------------+
	// 0 |  a   0   0   0  |
	// 1 |  0   b   0   0  |
	// 2 |  0   0   c   d  |
	// 3 |  0   0  -1   0  |
	//   +-----------------+

	return res;
}
template <arithmetic T>
matrix4_T<T> perspective_T<T>::to_inverse_matrix4() const
{
	matrix4_T<T> res{matrix4_T<T>::zero()};
	const T a = m_focalLength * m_invRatio;
	const T b = m_focalLength;
	const T farMinusNear = m_far - m_near;
	const T invFarMinusNear = T(1) / farMinusNear;
	const T c = m_near * invFarMinusNear;
	const T d = m_far * c;

	res.column(0).row(0) = T(1) / a;
	res.column(1).row(1) = T(1) / b;
	res.column(3).row(2) = T(-1);
	res.column(2).row(3) = T(1) / d;
	res.column(3).row(3) = c / d;

	//      0   1   2   3
	//   +-----------------+
	// 0 | 1/a  0   0   0  |
	// 1 |  0  1/b  0   0  |
	// 2 |  0   0   0  -1  |
	// 3 |  0   0  1/d c/d |
	//   +-----------------+

	return res;
}

template <arithmetic T>
matrix4_T<T> operator*(const matrix4_T<T>& m, const perspective_T<T>& p)
{
	const T a = p.m_focalLength * p.m_invRatio;
	const T b = p.m_focalLength;
	const T farMinusNear = p.m_far - p.m_near;
	const T invFarMinusNear = T(1) / farMinusNear;
	const T c = p.m_near * invFarMinusNear;
	const T d = p.m_far * c;

	matrix4_T<T> res;
	res.column(0).row(0) = a * m.column(0).row(0);
	res.column(1).row(0) = b * m.column(1).row(0);
	res.column(2).row(0) = c * m.column(2).row(0) - m.column(3).row(0);
	res.column(3).row(0) = d * m.column(2).row(0);
	res.column(0).row(1) = a * m.column(0).row(1);
	res.column(1).row(1) = b * m.column(1).row(1);
	res.column(2).row(1) = c * m.column(2).row(1) - m.column(3).row(1);
	res.column(3).row(1) = d * m.column(2).row(1);
	res.column(0).row(2) = a * m.column(0).row(2);
	res.column(1).row(2) = b * m.column(1).row(2);
	res.column(2).row(2) = c * m.column(2).row(2) - m.column(3).row(2);
	res.column(3).row(2) = d * m.column(2).row(2);
	res.column(0).row(3) = a * m.column(0).row(3);
	res.column(1).row(3) = b * m.column(1).row(3);
	res.column(2).row(3) = c * m.column(2).row(3) - m.column(3).row(3);
	res.column(3).row(3) = d * m.column(2).row(3);
	return res;
}
template <arithmetic T>
matrix4_T<T> operator*(const perspective_T<T>& p, const matrix4_T<T>& m)
{
	const T a = p.m_focalLength * p.m_invRatio;
	const T b = p.m_focalLength;
	const T farMinusNear = p.m_far - p.m_near;
	const T invFarMinusNear = T(1) / farMinusNear;
	const T c = p.m_near * invFarMinusNear;
	const T d = p.m_far * c;

	matrix4_T<T> res;
	res.column(0).row(0) = m.column(0).row(0) * a;
	res.column(1).row(0) = m.column(1).row(0) * a;
	res.column(2).row(0) = m.column(2).row(0) * a;
	res.column(3).row(0) = m.column(3).row(0) * a;
	res.column(0).row(1) = m.column(0).row(1) * b;
	res.column(1).row(1) = m.column(1).row(1) * b;
	res.column(2).row(1) = m.column(2).row(1) * b;
	res.column(3).row(1) = m.column(3).row(1) * b;
	res.column(0).row(2) = m.column(0).row(2) * c + m.column(0).row(3) * d;
	res.column(1).row(2) = m.column(1).row(2) * c + m.column(1).row(3) * d;
	res.column(2).row(2) = m.column(2).row(2) * c + m.column(2).row(3) * d;
	res.column(3).row(2) = m.column(3).row(2) * c + m.column(3).row(3) * d;
	res.column(0).row(3) = -m.column(0).row(2);
	res.column(1).row(3) = -m.column(1).row(2);
	res.column(2).row(3) = -m.column(2).row(2);
	res.column(3).row(3) = -m.column(3).row(2);
	return res;
}
template <arithmetic T>
vector4_T<T> operator*(const perspective_T<T>& p, const vector4_T<T>& v)
{
	const T a = p.m_focalLength * p.m_invRatio;
	const T b = p.m_focalLength;
	const T farMinusNear = p.m_far - p.m_near;
	const T invFarMinusNear = T(1) / farMinusNear;
	const T c = p.m_near * invFarMinusNear;
	const T d = p.m_far * c;

	//                       +---+
	//                       | x |
	//                       | y |
	//                       | z |
	//                       | w |
	//      0   1   2   3    +---+
	//   +-----------------+
	// 0 |  a   0   0   0  |
	// 1 |  0   b   0   0  |
	// 2 |  0   0   c   d  |
	// 3 |  0   0  -1   0  |
	//   +-----------------+

	return {
	    a * v.x,
	    b * v.y,
	    c * v.z + d * v.w,
	    -v.z,
	};
}
template <arithmetic T>
matrix4_T<T> operator*(const unscaled_transform3_T<T>& t,
                       const perspective_T<T>& p)
{
	const T a = p.m_focalLength * p.m_invRatio;
	const T b = p.m_focalLength;
	const T farMinusNear = p.m_far - p.m_near;
	const T invFarMinusNear = T(1) / farMinusNear;
	const T c = p.m_near * invFarMinusNear;
	const T d = p.m_far * c;
	const T xx = t.rotation.x * t.rotation.x;
	const T yy = t.rotation.y * t.rotation.y;
	const T zz = t.rotation.z * t.rotation.z;
	const T wx = t.rotation.w * t.rotation.x;
	const T wy = t.rotation.w * t.rotation.y;
	const T wz = t.rotation.w * t.rotation.z;
	const T xy = t.rotation.x * t.rotation.y;
	const T xz = t.rotation.x * t.rotation.z;
	const T yz = t.rotation.y * t.rotation.z;

	return {
	    a * (T(1) - T(2) * (yy + zz)),
	    b * (T(2) * (xy - wz)),
	    c * (T(2) * (xz + wy)) - t.position.x,
	    d * (T(2) * (xz + wy)),

	    a * (T(2) * (xy + wz)),
	    b * (T(1) - T(2) * (xx + zz)),
	    c * (T(2) * (yz - wx)) - t.position.y,
	    d * (T(2) * (yz - wx)),

	    a * (T(2) * (xz - wy)),
	    b * (T(2) * (yz + wx)),
	    c * (T(1) - T(2) * (xx + yy)) - t.position.z,
	    d * (T(1) - T(2) * (xx + yy)),

	    T(0),
	    T(0),
	    T(-1),
	    T(0),
	};
}
template <arithmetic T>
matrix4_T<T> operator*(const perspective_T<T>& p,
                       const unscaled_transform3_T<T>& t)
{
	const T a = p.m_focalLength * p.m_invRatio;
	const T b = p.m_focalLength;
	const T farMinusNear = p.m_far - p.m_near;
	const T invFarMinusNear = T(1) / farMinusNear;
	const T c = p.m_near * invFarMinusNear;
	const T d = p.m_far * c;
	const T xx = t.rotation.x * t.rotation.x;
	const T yy = t.rotation.y * t.rotation.y;
	const T zz = t.rotation.z * t.rotation.z;
	const T wx = t.rotation.w * t.rotation.x;
	const T wy = t.rotation.w * t.rotation.y;
	const T wz = t.rotation.w * t.rotation.z;
	const T xy = t.rotation.x * t.rotation.y;
	const T xz = t.rotation.x * t.rotation.z;
	const T yz = t.rotation.y * t.rotation.z;

	matrix4_T<T> res;

	res.row(0).column(0) = (T(1) - T(2) * (yy + zz)) * a;
	res.row(0).column(1) = (T(2) * (xy - wz)) * a;
	res.row(0).column(2) = (T(2) * (xz + wy)) * a;
	res.row(0).column(3) = t.position.x * a;

	res.row(1).column(0) = (T(2) * (xy + wz)) * b;
	res.row(1).column(1) = (T(1) - T(2) * (xx + zz)) * b;
	res.row(1).column(2) = (T(2) * (yz - wx)) * b;
	res.row(1).column(3) = t.position.y * b;

	res.row(2).column(0) = (T(2) * (xz - wy)) * c;
	res.row(2).column(1) = (T(2) * (yz + wx)) * c;
	res.row(2).column(2) = (T(1) - T(2) * (xx + yy)) * c;
	res.row(2).column(3) = t.position.z * c + d;

	res.row(3).column(0) = -(T(2) * (xz - wy));
	res.row(3).column(1) = -(T(2) * (yz + wx));
	res.row(3).column(2) = -(T(1) - T(2) * (xx + yy));
	res.row(3).column(3) = -t.position.z;

	return res;
}
#pragma endregion

#pragma region definition_direction2_T
template <arithmetic T>
constexpr direction2_T<T>::direction2_T(T x_, T y_)
    : direction2_T(vector2_T<T>(x_, y_))
{
}

template <arithmetic T>
constexpr direction2_T<T>::direction2_T(const vector2_T<T>& v)
{
	vector2_T<T> vn = v.get_normalized();
	m_x = vn.x;
	m_y = vn.y;
}

template <arithmetic T>
constexpr direction2_T<T> direction2_T<T>::operator-() const
{
	return already_normalized(-m_x, -m_y);
}

template <arithmetic T>
constexpr const T& direction2_T<T>::x() const
{
	return m_x;
}

template <arithmetic T>
constexpr const T& direction2_T<T>::y() const
{
	return m_y;
}

template <arithmetic T>
constexpr std::array<T, 2> direction2_T<T>::to_array() const
{
	return {
	    m_x,
	    m_y,
	};
}

template <arithmetic T>
template <typename U>
constexpr direction2_T<T>::operator direction2_T<U>() const
{
	return direction2_T<U>::already_normalized(static_cast<U>(m_x),
	                                           static_cast<U>(m_y));
}

template <arithmetic T>
constexpr bool direction2_T<T>::operator==(const direction2_T<T>& v) const
{
	return m_x == v.m_x && m_y == v.m_y;
}

template <arithmetic T>
constexpr bool direction2_T<T>::operator!=(const direction2_T<T>& v) const
{
	return !operator==(v);
}

template <arithmetic T>
constexpr bool direction2_T<T>::nearly_equal(const direction2_T<T>& v,
                                             T e) const
{
	return cdm::nearly_equal(*this, v, e);
}

template <arithmetic T>
template <typename Index>
constexpr const T& direction2_T<T>::operator[](Index index) const
{
	switch (index)
	{
	case 0:
		return m_x;
	case 1:
		return m_y;
	}

#ifndef NDEBUG
	throw std::out_of_range("attempting to index element " +
	                        std::to_string(index));
#else
	return x;
#endif
}

template <arithmetic T>
constexpr direction2_T<T> direction2_T<T>::already_normalized(T x_, T y_)
{
	direction2_T<T> res;
	res.m_x = x_;
	res.m_y = y_;
	return res;
}

template <arithmetic T>
constexpr direction2_T<T> direction2_T<T>::already_normalized(
    const vector2_T<T>& v)
{
	return already_normalized(v.x, v.y);
}

template <arithmetic T>
constexpr direction2_T<T> direction2_T<T>::posX()
{
	return already_normalized(T(1), T(0));
}

template <arithmetic T>
constexpr direction2_T<T> direction2_T<T>::posY()
{
	return already_normalized(T(0), T(1));
}

template <arithmetic T>
constexpr direction2_T<T> direction2_T<T>::negX()
{
	return already_normalized(T(-1), T(0));
}

template <arithmetic T>
constexpr direction2_T<T> direction2_T<T>::negY()
{
	return already_normalized(T(0), T(-1));
}

template <arithmetic T>
constexpr vector2_T<T> operator*(T t, const direction2_T<T>& d)
{
	return t * vector2_T<T>(d);
}

template <arithmetic T>
constexpr vector2_T<T> operator*(const direction2_T<T>& d, T t)
{
	return vector2_T<T>(d) * t;
}

template <arithmetic T>
constexpr T dot(direction2_T<T> lhs, direction2_T<T> rhs)
{
	return dot(vector2_T<T>(lhs), vector2_T<T>(rhs));
}
template <arithmetic T>
constexpr T cross(direction2_T<T> lhs, direction2_T<T> rhs)
{
	return cross(vector2_T<T>(lhs), vector2_T<T>(rhs));
}

template <arithmetic T>
constexpr T dot(direction2_T<T> lhs, vector2_T<T> rhs)
{
	return dot(vector2_T<T>(lhs), vector2_T<T>(rhs));
}
template <arithmetic T>
constexpr T cross(direction2_T<T> lhs, vector2_T<T> rhs)
{
	return cross(vector2_T<T>(lhs), vector2_T<T>(rhs));
}

template <arithmetic T>
constexpr T dot(vector2_T<T> lhs, direction2_T<T> rhs)
{
	return dot(vector2_T<T>(lhs), vector2_T<T>(rhs));
}
template <arithmetic T>
constexpr T cross(vector2_T<T> lhs, direction2_T<T> rhs)
{
	return cross(vector2_T<T>(lhs), vector2_T<T>(rhs));
}
#pragma endregion

#pragma region definition_direction3_T
template <arithmetic T>
constexpr direction3_T<T>::direction3_T(T x_, T y_, T z_)
    : direction3_T(vector3_T<T>(x_, y_, z_))
{
}

template <arithmetic T>
constexpr direction3_T<T>::direction3_T(const vector3_T<T>& v)
{
	vector3_T<T> vn = v.get_normalized();
	m_x = vn.x;
	m_y = vn.y;
	m_z = vn.z;
}

template <arithmetic T>
constexpr direction3_T<T> direction3_T<T>::operator-() const
{
	return already_normalized(-m_x, -m_y, -m_z);
}

template <arithmetic T>
constexpr const T& direction3_T<T>::x() const
{
	return m_x;
}

template <arithmetic T>
constexpr const T& direction3_T<T>::y() const
{
	return m_y;
}

template <arithmetic T>
constexpr const T& direction3_T<T>::z() const
{
	return m_z;
}

template <arithmetic T>
constexpr std::array<T, 3> direction3_T<T>::to_array() const
{
	return {
	    m_x,
	    m_y,
	    m_z,
	};
}

template <arithmetic T>
template <typename U>
constexpr direction3_T<T>::operator direction3_T<U>() const
{
	return direction3_T<U>::already_normalized(
	    static_cast<U>(m_x), static_cast<U>(m_y), static_cast<U>(m_z));
}

template <arithmetic T>
constexpr bool direction3_T<T>::operator==(const direction3_T<T>& v) const
{
	return m_x == v.m_x && m_y == v.m_y && m_z == v.m_z;
}

template <arithmetic T>
constexpr bool direction3_T<T>::operator!=(const direction3_T<T>& v) const
{
	return !operator==(v);
}

template <arithmetic T>
constexpr bool direction3_T<T>::nearly_equal(const direction3_T<T>& v,
                                             T e) const
{
	return cdm::nearly_equal(*this, v, e);
}

template <arithmetic T>
template <typename Index>
constexpr const T& direction3_T<T>::operator[](Index index) const
{
	switch (index)
	{
	case 0:
		return m_x;
	case 1:
		return m_y;
	case 2:
		return m_z;
	}

#ifndef NDEBUG
	throw std::out_of_range("attempting to index element " +
	                        std::to_string(index));
#else
	return x;
#endif
}

template <arithmetic T>
constexpr direction3_T<T> direction3_T<T>::already_normalized(T x_, T y_, T z_)
{
	direction3_T<T> res;
	res.m_x = x_;
	res.m_y = y_;
	res.m_z = z_;
	return res;
}

template <arithmetic T>
constexpr direction3_T<T> direction3_T<T>::already_normalized(
    const vector3_T<T>& v)
{
	return already_normalized(v.x, v.y, v.z);
}

template <arithmetic T>
constexpr direction3_T<T> direction3_T<T>::posX()
{
	return already_normalized(T(1), T(0), T(0));
}

template <arithmetic T>
constexpr direction3_T<T> direction3_T<T>::posY()
{
	return already_normalized(T(0), T(1), T(0));
}

template <arithmetic T>
constexpr direction3_T<T> direction3_T<T>::posZ()
{
	return already_normalized(T(0), T(0), T(1));
}

template <arithmetic T>
constexpr direction3_T<T> direction3_T<T>::negX()
{
	return already_normalized(T(-1), T(0), T(0));
}

template <arithmetic T>
constexpr direction3_T<T> direction3_T<T>::negY()
{
	return already_normalized(T(0), T(-1), T(0));
}

template <arithmetic T>
constexpr direction3_T<T> direction3_T<T>::negZ()
{
	return already_normalized(T(0), T(0), T(-1));
}

template <arithmetic T>
constexpr vector3_T<T> operator*(T t, const direction3_T<T>& d)
{
	return t * vector3_T<T>(d);
}

template <arithmetic T>
constexpr vector3_T<T> operator*(const direction3_T<T>& d, T t)
{
	return vector3_T<T>(d) * t;
}

template <arithmetic T>
constexpr T dot(direction3_T<T> lhs, direction3_T<T> rhs)
{
	return dot(vector3_T<T>(lhs), vector3_T<T>(rhs));
}
template <arithmetic T>
constexpr vector3_T<T> cross(direction3_T<T> lhs, direction3_T<T> rhs)
{
	return cross(vector3_T<T>(lhs), vector3_T<T>(rhs));
}

template <arithmetic T>
constexpr T dot(direction3_T<T> lhs, vector3_T<T> rhs)
{
	return dot(vector3_T<T>(lhs), vector3_T<T>(rhs));
}
template <arithmetic T>
constexpr vector3_T<T> cross(direction3_T<T> lhs, vector3_T<T> rhs)
{
	return cross(vector3_T<T>(lhs), vector3_T<T>(rhs));
}

template <arithmetic T>
constexpr T dot(vector3_T<T> lhs, direction3_T<T> rhs)
{
	return dot(vector3_T<T>(lhs), vector3_T<T>(rhs));
}
template <arithmetic T>
constexpr vector3_T<T> cross(vector3_T<T> lhs, direction3_T<T> rhs)
{
	return cross(vector3_T<T>(lhs), vector3_T<T>(rhs));
}
#pragma endregion

#pragma region definition_quaternion_t
template <arithmetic T>
quaternion_T<T>::quaternion_T(direction3_T<T> axis, radian_T<T> angle)
{
	vector3_T<T> tmpAxis = axis * sin(angle / T(2.0));
	*this = {tmpAxis.x, tmpAxis.y, tmpAxis.z, cos(angle / T(2.0))};
}
template <arithmetic T>
template <std::signed_integral U, U NumeratorT, U DenominatorT>
quaternion_T<T>::quaternion_T(
    direction3_T<T> axis,
    static_pi_fraction_T<U, NumeratorT, DenominatorT> angle)
{
	using HalfAngle = static_pi_fraction_T<U, NumeratorT, DenominatorT * U(2)>;

	vector3_T<T> tmpAxis = axis * sin<T>(HalfAngle{});
	*this = {tmpAxis.x, tmpAxis.y, tmpAxis.z, cos<T>(HalfAngle{})};
}

template <arithmetic T>
quaternion_T<T> quaternion_T<T>::identity()
{
	return {T(0), T(0), T(0), T(1)};
}

template <arithmetic T>
quaternion_T<T> quaternion_T<T>::from_rotation_matrix(const matrix3_T<T>& m)
{
	const auto diag = m.diag();

	const T fourXSquaredMinus1 = diag.x - diag.y - diag.z;
	const T fourYSquaredMinus1 = diag.y - diag.x - diag.z;
	const T fourZSquaredMinus1 = diag.z - diag.x - diag.y;
	const T fourWSquaredMinus1 = diag.x + diag.y + diag.z;

	int biggestIndex = 0;
	T fourBiggestSquaredMinus1 = fourWSquaredMinus1;
	if (fourXSquaredMinus1 > fourBiggestSquaredMinus1)
	{
		fourBiggestSquaredMinus1 = fourXSquaredMinus1;
		biggestIndex = 1;
	}
	if (fourYSquaredMinus1 > fourBiggestSquaredMinus1)
	{
		fourBiggestSquaredMinus1 = fourYSquaredMinus1;
		biggestIndex = 2;
	}
	if (fourZSquaredMinus1 > fourBiggestSquaredMinus1)
	{
		fourBiggestSquaredMinus1 = fourZSquaredMinus1;
		biggestIndex = 3;
	}

	const T biggestVal =
		std::sqrt(fourBiggestSquaredMinus1 + T(1)) *
		T(0.5);
	const T mult = T(0.25) / biggestVal;

	switch (biggestIndex)
	{
	case 0:
		return {
			(m.column(1).row(2) - m.column(2).row(1)) * mult,
			(m.column(2).row(0) - m.column(0).row(2)) * mult,
			(m.column(0).row(1) - m.column(1).row(0)) * mult,
			biggestVal,
		};
	case 1:
		return {
			biggestVal,
			(m.column(0).row(1) + m.column(1).row(0)) * mult,
			(m.column(2).row(0) + m.column(0).row(2)) * mult,
			(m.column(1).row(2) - m.column(2).row(1)) * mult,
		};
	case 2:
		return {
			(m.column(0).row(1) + m.column(1).row(0)) * mult,
			biggestVal,
			(m.column(1).row(2) + m.column(2).row(1)) * mult,
			(m.column(2).row(0) - m.column(0).row(2)) * mult,
		};
	case 3:
		return {
			(m.column(2).row(0) + m.column(0).row(2)) * mult,
			(m.column(1).row(2) + m.column(2).row(1)) * mult,
			biggestVal,
			(m.column(0).row(1) - m.column(1).row(0)) * mult,
		};

	default:
		assert(false);
		return {
			T(0),
			T(0),
			T(0),
			T(1),
		};
	}
}

template <arithmetic T>
quaternion_T<T> quaternion_T<T>::look_at(direction3_T<T> forward, direction3_T<T> upward)
{
	const direction3_T<T> backward = -forward;
	const direction3_T<T> rightward = direction3_T<T>(cross(upward, backward));
	upward = direction3_T<T>(cross(backward, rightward));

	const matrix3_T<T> rotMat = matrix3_T<T>::columns(rightward, upward, backward);

	return from_rotation_matrix(rotMat);
}

template <arithmetic T>
T quaternion_T<T>::norm() const
{
	return std::sqrt(norm_squared());
}
template <arithmetic T>
T quaternion_T<T>::norm_squared() const
{
	return x * x + y * y + z * z + w * w;
}

template <arithmetic T>
quaternion_T<T>& quaternion_T<T>::inverse()
{
	const T n = norm();
	if (nearly_equal(n, T(1)))
		return conjugate();

	x /= -n;
	y /= -n;
	z /= -n;
	w /= n;
	return *this;
}
template <arithmetic T>
quaternion_T<T>& quaternion_T<T>::conjugate()
{
	x = -x;
	y = -y;
	z = -z;
	return *this;
}
template <arithmetic T>
quaternion_T<T>& quaternion_T<T>::normalize()
{
	const T n = norm();
	x /= n;
	y /= n;
	z /= n;
	w /= n;
	return *this;
}
template <arithmetic T>
constexpr quaternion_T<T>& quaternion_T<T>::clamp(quaternion_T<T> min,
                                        quaternion_T<T> max)
{
	x = std::clamp(x, min.x, max.x);
	y = std::clamp(y, min.y, max.y);
	z = std::clamp(z, min.z, max.z);
	w = std::clamp(w, min.w, max.w);
	return *this;
}

template <arithmetic T>
quaternion_T<T> quaternion_T<T>::get_inversed() const
{
	quaternion_T<T> res = *this;
	res.inverse();
	return res;
}
template <arithmetic T>
quaternion_T<T> quaternion_T<T>::get_conjugated() const
{
	quaternion_T<T> res = *this;
	res.conjugate();
	return res;
}
template <arithmetic T>
quaternion_T<T> quaternion_T<T>::get_normalized() const
{
	quaternion_T<T> res = *this;
	res.normalize();
	return res;
}
template <arithmetic T>
constexpr quaternion_T<T> quaternion_T<T>::get_clamped(quaternion_T<T> min,
                                             quaternion_T<T> max) const
{
	quaternion_T<T> res = *this;
	res.clamp(min, max);
	return res;
}

template <arithmetic T>
euler_angles_T<T> quaternion_T<T>::to_euler_angles() const
{
	return cdm::euler_angles_T<T>{
	    cdm::radian_T<T>(
	        std::atan2(T(2) * (w * x + y * z), T(1) - T(2) * (x * x + y * y))),
	    cdm::radian_T<T>(std::asin(T(2) * (w * y - z * x))),
	    cdm::radian_T<T>(
	        std::atan2(T(2) * (w * z + x * y), T(1) - T(2) * (y * y + z * z))),
	};
}

template <arithmetic T>
quaternion_T<T> quaternion_T<T>::operator+(quaternion_T<T> q) const
{
	return {
	    x + q.x,
	    y + q.y,
	    z + q.z,
	    w + q.w,
	};
}
template <arithmetic T>
quaternion_T<T> quaternion_T<T>::operator-(quaternion_T<T> q) const
{
	return {
	    x - q.x,
	    y - q.y,
	    z - q.z,
	    w - q.w,
	};
}
template <arithmetic T>
vector3_T<T> quaternion_T<T>::operator*(vector3_T<T> v) const
{
	const vector3_T<T> qvec{x, y, z};
	vector3_T<T> uv{cross(qvec, v)};
	vector3_T<T> uuv{cross(qvec, uv)};
	uv = uv * (T(2) * w);
	uuv = uuv * T(2);
	return v + uv + uuv;
}
template <arithmetic T>
quaternion_T<T> quaternion_T<T>::operator*(quaternion_T<T> q) const
{
	return {
	    x * q.w + y * q.z - z * q.y + w * q.x,
	    -x * q.z + y * q.w + z * q.x + w * q.y,
	    x * q.y - y * q.x + z * q.w + w * q.z,
	    -x * q.x - y * q.y - z * q.z + w * q.w,
	};
}
template <arithmetic T>
quaternion_T<T> quaternion_T<T>::operator*(T f) const
{
	return {
	    x * f,
	    y * f,
	    z * f,
	    w * f,
	};
}
template <arithmetic T>
quaternion_T<T> quaternion_T<T>::operator/(T f) const
{
	return {
	    x / f,
	    y / f,
	    z / f,
	    w / f,
	};
}

template <arithmetic T>
quaternion_T<T>& quaternion_T<T>::operator+=(quaternion_T<T> q)
{
	*this = *this + q;
	return *this;
}
template <arithmetic T>
quaternion_T<T>& quaternion_T<T>::operator-=(quaternion_T<T> q)
{
	*this = *this - q;
	return *this;
}
template <arithmetic T>
quaternion_T<T>& quaternion_T<T>::operator*=(quaternion_T<T> q)
{
	*this = *this * q;
	return *this;
}
template <arithmetic T>
quaternion_T<T>& quaternion_T<T>::operator*=(T f)
{
	*this = *this * f;
	return *this;
}
template <arithmetic T>
quaternion_T<T>& quaternion_T<T>::operator/=(T f)
{
	*this = *this / f;
	return *this;
}

template <arithmetic T>
quaternion_T<T> quaternion_T<T>::operator-() const
{
	return {
	    -x,
	    -y,
	    -z,
	    -w,
	};
}

template <arithmetic T>
bool quaternion_T<T>::operator==(quaternion_T<T> q) const
{
	return x == q.x &&  //
	       y == q.y &&  //
	       z == q.z &&  //
	       w == q.w;    //
}
template <arithmetic T>
bool quaternion_T<T>::operator!=(quaternion_T<T> q) const
{
	return !operator==(q);
}

template <arithmetic T>
vector3_T<T> operator*(const normalized<quaternion_T<T>>& q, vector3_T<T> v)
{
	return *q * v;
}

template <arithmetic T>
constexpr T dot(quaternion_T<T> lhs, quaternion_T<T> rhs)
{
	return lhs.x * rhs.x +  //
	       lhs.y * rhs.y +  //
	       lhs.z * rhs.z +  //
	       lhs.w * rhs.w;   //
}
template <arithmetic T>
constexpr quaternion_T<T> cross(quaternion_T<T> lhs, quaternion_T<T> rhs)
{
	return {
	    lhs.y * rhs.z - lhs.z * rhs.y,  //
	    lhs.z * rhs.x - lhs.x * rhs.z,  //
	    lhs.x * rhs.y - lhs.y * rhs.x,  //
	    T(1)                            //
	};
}
template <arithmetic T>
quaternion_T<T> lerp(quaternion_T<T> begin, quaternion_T<T> end, T percent)
{
	return (end - begin) * percent + begin;
}
template <arithmetic T>
quaternion_T<T> nlerp(quaternion_T<T> begin, quaternion_T<T> end, T percent)
{
	return lerp(begin, end, percent).get_normalized();
}
template <arithmetic T>
quaternion_T<T> slerp(quaternion_T<T> begin, quaternion_T<T> end, T percent)
{
	quaternion_T<T> _end;
	T d = dot(begin, end);

	if (d > 0.9995f)
		return lerp(begin, end, percent);
	if (d < T(0))
	{
		_end = -end;
		d = -d;
	}
	else
		_end = end;

	percent = cdm::clamp(percent, T(0), T(1));
	d = cdm::clamp(d, T(-1), T(1));
	T theta = acos(d) * percent;

	quaternion_T<T> res = begin * d;
	res = _end - res;

	res = begin * cos(theta) + res * sin(theta);

	return res;
}
#pragma endregion

#pragma region definition_line_t
#pragma region slope_intercept
template <arithmetic T>
line_T<T, line_representation::SlopeIntercept>::line_T(T s, T yi) : slope(s), y_intercept(yi) {}
template <arithmetic T>
line_T<T, line_representation::SlopeIntercept>::line_T(vector2_T<T> direction, T e) : y_intercept(T(0))
{
	direction.normalize();

	if (nearly_equal(direction.y, T(1), e) || nearly_equal(direction.y, T(-1), e))
		throw std::runtime_error("this line can not be represented with this structure");

	// project the direction on the X axis
	const T d = dot(direction, vector2_T<T>(T(1), T(0)));

	// clang-format off
	//
	//   |
	//   |             +
	//   |             |
	//   |           /
	//   |             |
	//   |         /
	//   |             |
	//   |       /
	//   |             |
	//   |   v /
	//   |    +        |
	//   |   /|
	//   |  /          |
	//   | /  |
	//   |/            |
	// --+----+--------+----
	//   |0    d        1
	// / |
	//   |
	//
	// clang-format on

	// adjust the vector's length so that x == 1
	direction = direction / d;
	// y is now the coefficient

	slope = direction.y;
}
template <arithmetic T>
line_T<T, line_representation::SlopeIntercept>::line_T(vector2_T<T> point1,
                                                       vector2_T<T> point2)
{
	slope = (point2.y - point1.y) / (point2.x - point1.x);
	y_intercept = slope * (-point1.x) + point1.y;
}

template <arithmetic T>
vector2_T<T> line_T<T, line_representation::SlopeIntercept>::resolve_for_x(
    T x) const
{
	return {x, slope * x + y_intercept};
}
template <arithmetic T>
vector2_T<T> line_T<T, line_representation::SlopeIntercept>::resolve_for_y(
    T y) const
{
	return {(y - y_intercept) / slope, y};
}

template <arithmetic T>
bool are_parallel(line_T<T, line_representation::SlopeIntercept> l1,
                  line_T<T, line_representation::SlopeIntercept> l2,
                  T e)
{
	return nearly_equal(l1.coefficient, l2.coefficient, e);
}

template <arithmetic T>
bool collides(line_T<T, line_representation::SlopeIntercept> l1,
              line_T<T, line_representation::SlopeIntercept> l2)
{
	return !are_parallel(l1, l2);
}

template <arithmetic T>
bool collides(line_T<T, line_representation::SlopeIntercept> l1,
              line_T<T, line_representation::SlopeIntercept> l2,
              vector2_T<T>& intersection)
{
	bool collision = collides(l1, l2);

	if (collision)
	{
		T a = l1.coefficient - l2.coefficient;
		T b = l2.offset - l1.offset;
		intersection.x = b / a;
		intersection.y = l1.coefficient * intersection.x + l1.offset;

		assert(intersection.y ==
		       (l2.coefficient * intersection.x + l2.offset));
	}

	return collision;
}
#pragma endregion

#pragma region points
template <arithmetic T>
line_T<T, line_representation::Points>::line_T(vector2_T<T> point1,
                                               vector2_T<T> point2)
    : p1(point1), p2(point2)
{
}

// template <arithmetic T>
// vector2_T<T> line_T<T, line_representation::Points>::resolve_for_x(T x)
// const
//{
//	const auto direction = from_to(p1, p2).get_normalized();
//
//
// }
// template <arithmetic T>
// vector2_T<T> line_T<T, line_representation::Points>::resolve_for_y(T y)
// const
//{
// }

template <arithmetic T>
void collides(line_T<T, line_representation::Points> l1,
              line_T<T, line_representation::Points> l2,
              vector2_T<T>& intersection)
{
	const T x1 = l1.p1.x;
	const T x2 = l1.p2.x;
	const T x3 = l2.p1.x;
	const T x4 = l2.p2.x;
	const T y1 = l1.p1.y;
	const T y2 = l1.p2.y;
	const T y3 = l2.p1.y;
	const T y4 = l2.p2.y;

	return {
	    ((x1 * y1 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * x4 - y3 * x4)) /
	        ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)),
	    ((x1 * y1 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * x4 - y3 * x4)) /
	        ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))};
}
#pragma endregion

template <arithmetic T>
T distance_between(plane_T<T> p, vector3_T<T> point)
{
	// clang-format off
	return -point.x * p.normal.x()
	       -point.y * p.normal.y()
	       -point.z * p.normal.z();
	// clang-format on
}
template <arithmetic T>
T distance_between(vector3_T<T> v, plane_T<T> p)
{
	return distance_between(p, v);
}
#pragma endregion

#pragma region definition_segment2_t
template <arithmetic T>
T segment2_T<T>::length() const
{
	return distance_between(origin, end);
}
template <arithmetic T>
T segment2_T<T>::length_squared() const
{
	return distance_squared_between(origin, end);
}

template <arithmetic T>
segment2_T<T>& segment2_T<T>::invert()
{
	std::swap(origin, end);
	return *this;
	;
}

template <arithmetic T>
int collides(const segment2_T<T>& s0,
             const segment2_T<T>& s1,
             vector2_T<T>& outPoint0,
             vector2_T<T>& outPoint1,
             T e)
{
	auto det = [](T a, T b, T c, T d) -> T { return a * d - b * c; };

	auto betw = [&](T l, T r, T x) -> T
	{ return std::min(l, r) <= x + e && x <= std::max(l, r) + e; };

	auto intersect_1d = [&](T a, T b, T c, T d) -> bool
	{
		if (a > b)
			std::swap(a, b);
		if (c > d)
			std::swap(c, d);
		return std::max(a, c) <= std::min(b, d) + e;
	};

	auto compareP = [&](const cdm::vector2_T<T>& l,
	                    const cdm::vector2_T<T>& r) -> bool
	{ return l.x < r.x - e || (std::abs(l.x - r.x) < e && l.y < r.y - e); };

	vector2_T<T> a = s0.origin;
	vector2_T<T> b = s0.end;
	vector2_T<T> c = s1.origin;
	vector2_T<T> d = s1.end;
	if (!intersect_1d(a.x, b.x, c.x, d.x) || !intersect_1d(a.y, b.y, c.y, d.y))
		return 0;

	struct line
	{
		T a, b, c;

		line() = default;
		line(const cdm::vector2_T<T>& p, const cdm::vector2_T<T>& q)
		{
			a = p.y - q.y;
			b = q.x - p.x;
			c = -a * p.x - b * p.y;
			norm();
		}

		void norm()
		{
			T z = sqrt(a * a + b * b);
			if (nearly_equal(std::abs(z), T(0)))
			{
				a /= z;
				b /= z;
				c /= z;
			}
		}

		T dist(const cdm::vector2_T<T>& p) const
		{
			return a * p.x + b * p.y + c;
		}
	};

	line m(a, b);
	line n(c, d);
	T zn = det(m.a, m.b, n.a, n.b);
	if (abs(zn) < e)
	{
		if (abs(m.dist(c)) > e || abs(n.dist(a)) > e)
			return 0;
		if (compareP(b, a))
			std::swap(a, b);
		if (compareP(d, c))
			std::swap(c, d);

		if (compareP(a, c))
			outPoint0 = c;
		else
			outPoint0 = a;

		if (compareP(b, d))
			outPoint1 = b;
		else
			outPoint1 = d;

		return outPoint0 == outPoint1 ? 1 : 2;
	}
	else
	{
		outPoint0.x = outPoint1.x = -det(m.c, m.b, n.c, n.b) / zn;
		outPoint0.y = outPoint1.y = -det(m.a, m.c, n.a, n.c) / zn;
		return int(betw(a.x, b.x, outPoint0.x) &&
		           betw(a.y, b.y, outPoint0.y) &&
		           betw(c.x, d.x, outPoint0.x) && betw(c.y, d.y, outPoint0.y));
	}
}
#pragma endregion

#pragma region definition_segment3_t
template <arithmetic T>
constexpr bool collides(const segment3_T<T>& seg,
                        const plane_T<T>& plane,
                        T e) noexcept
{
	// distances to point
	const T p0 = plane.evaluate(seg.origin);
	const T p1 = plane.evaluate(seg.end);

	// are considered on the plane
	if (nearly_equal(p0, T(0), e) || nearly_equal(p1, T(0), e))
		return true;

	// if both points are on the same side of the plane
	if ((p0 > T(0) && p1 > T(0)) || (p0 < T(0) && p1 < T(0)))
		return false;  // no collision

	return true;
}

// I feel like there's a simpler and more efficient way to do this...
template <arithmetic T>
constexpr bool collides(const segment3_T<T>& seg,
                        const plane_T<T>& plane,
                        vector3_T<T>& outPoint,
                        T e) noexcept
{
	// distances to point
	const T p0 = plane.evaluate(seg.origin);
	const T p1 = plane.evaluate(seg.end);

	// are considered on the plane
	if (nearly_equal(p0, T(0), e))
	{
		outPoint = seg.origin;
		return true;
	}
	if (nearly_equal(p1, T(0), e))
	{
		outPoint = seg.end;
		return true;
	}

	// if both points are on the same side of the plane
	if ((p0 > T(0) && p1 > T(0)) || (p0 < T(0) && p1 < T(0)))
		return false;  // no collision

	const value_domain_T<T> d{p0, p1};
	const unnormalized_value_T<T> v{d, T(0)};
	const normalized_value_T<T> n{v};

	outPoint = lerp(seg.origin, seg.end, n.value());

	return true;
}
#pragma endregion

#pragma region definition_plane_t
template <arithmetic T>
constexpr T plane_T<T>::evaluate(const vector3_T<T>& point) const
{
	return normal.x() * (point.x - origin.x) +  //
	       normal.y() * (point.y - origin.y) +  //
	       normal.z() * (point.z - origin.z);   //
}
template <arithmetic T>
constexpr vector3_T<T> plane_T<T>::project3d(const vector3_T<T>& point) const
{
	T distance = evaluate(point);
	return point - normal * distance;
}
template <arithmetic T>
constexpr vector2_T<T> plane_T<T>::project2d(
    const vector3_T<T>& point,
    const direction3_T<T>& plane_tangent) const
{
	const auto bitangent = direction3_T<T>(cross(normal, plane_tangent));
	const auto tangent = direction3_T<T>(cross(bitangent, normal));

	const auto TBN = matrix3_T<T>::rows(tangent, bitangent, normal);

	const vector3_T<T> plane_space_point = point - origin;

	return (TBN * plane_space_point).xy();
}
template <arithmetic T>
constexpr vector3_T<T> plane_T<T>::unproject(
    const vector2_T<T>& point,
    const direction3_T<T>& plane_tangent) const
{
	const auto bitangent = direction3_T<T>(cross(normal, plane_tangent));
	const auto tangent = direction3_T<T>(cross(bitangent, normal));

	const auto invTBN =
	    matrix3_T<T>::rows(tangent, bitangent, normal).get_inversed();

	const auto plane_space_point = vector3_T<T>{point, T(0)};

	return (invTBN * plane_space_point) + origin;
}

template <arithmetic T>
constexpr vector3_T<T> plane_T<T>::symmetry(const vector3_T<T>& point) const
{
	auto projP = project3d(point);
	auto projP2P = cdm::from_to(projP, point);
	auto symP = point - T(2) * projP2P;
	return symP;
}

template <arithmetic T>
constexpr direction3_T<T> plane_T<T>::computeTangent(T epsilon_) const
{
	const auto tangent =
	    project3d(vector3_T{T(1), T(0), T(0)} + origin) - origin;
	if (norm_squared(tangent) < epsilon_)
		return direction3_T<T>(
		    project3d(vector3_T{T(0), T(1), T(0)} + origin) - origin);
	return direction3_T<T>(tangent);
}

template <arithmetic T>
bool collides(const plane_T<T>& plane,
              ray3_T<T> r,
              vector3_T<T>& collision_point)
{
	const T denom = dot(*plane.normal, r.direction);
	if (abs(denom) > T(0.0001))
	{
		const T t = -dot(plane.origin - r.origin, plane.normal) / denom;
		if (t >= T(0))
		{
			collision_point = r.origin + r.direction * t;
			return true;
		}
	}
	return false;
}
template <arithmetic T>
bool collides_bidirectional(const plane_T<T>& plane,
                            ray3_T<T> r,
                            vector3_T<T>& collision_point)
{
	const T denom = dot(plane.normal, r.direction);
	if (abs(denom) > T(0.0001))
	{
		const T t = -dot(plane.origin - r.origin, plane.normal) / denom;

		collision_point = r.origin + r.direction * t;
		return true;
	}
	return false;
}
#pragma endregion

#pragma region definition_oriented_plane_t
template <arithmetic T>
constexpr oriented_plane_T<T>::oriented_plane_T(
    const vector3_T<T>& origin_,
    const direction3_T<T>& normal_,
    const direction3_T<T>& tangent_)
    : origin(origin_), normal(normal_)
{
	const auto bitangent = direction3_T<T>(cross(normal, tangent_));
	tangent = direction3_T<T>(cross(bitangent, normal));
}
template <arithmetic T>
constexpr oriented_plane_T<T>::oriented_plane_T(
    const plane_T<T>& plane,
    const direction3_T<T>& tangent_)
    : oriented_plane_T(plane.origin, plane.normal, tangent_)
{
}

template <arithmetic T>
constexpr T oriented_plane_T<T>::evaluate(const vector3_T<T>& point) const
{
	return normal.x() * (point.x - origin.x) +  //
	       normal.y() * (point.y - origin.y) +  //
	       normal.z() * (point.z - origin.z);   //
}
template <arithmetic T>
constexpr vector3_T<T> oriented_plane_T<T>::project3d(
    const vector3_T<T>& point) const
{
	T distance = evaluate(point);
	return point - normal * distance;
}
template <arithmetic T>
constexpr vector2_T<T> oriented_plane_T<T>::project2d(
    const vector3_T<T>& point) const
{
	const auto bitangent = direction3_T<T>(cross(normal, tangent));

	const auto TBN = matrix3_T<T>::rows(tangent, bitangent, normal);

	const vector3_T<T> pointInPlaneSpace = point - origin;

	return (TBN * pointInPlaneSpace).xy();
}
template <arithmetic T>
constexpr vector3_T<T> oriented_plane_T<T>::unproject(
    const vector2_T<T>& point) const
{
	const auto bitangent = direction3_T<T>(cross(normal, tangent));

	const auto invTBN =
	    matrix3_T<T>::rows(tangent, bitangent, normal).get_inversed();

	const auto pointInPlaneSpace = vector3_T<T>{point, T(0)};

	return (invTBN * pointInPlaneSpace) + origin;
}

template <arithmetic T>
constexpr vector3_T<T> oriented_plane_T<T>::symmetry(
    const vector3_T<T>& point) const
{
	auto projP = project3d(point);
	auto projP2P = cdm::from_to(projP, point);
	auto symP = point - T(2) * projP2P;
	return symP;
}

template <arithmetic T>
bool collides(const oriented_plane_T<T>& plane,
              ray3_T<T> r,
              vector3_T<T>& collision_point)
{
	const T denom = dot(plane.normal, r.direction);
	if (abs(denom) > T(0.0001))
	{
		const T t = -dot(plane.origin - r.origin, plane.normal) / denom;
		if (t >= T(0))
		{
			collision_point = r.origin + r.direction * t;
			return true;
		}
	}
	return false;
}
template <arithmetic T>
bool collides_bidirectional(const oriented_plane_T<T>& plane,
                            ray3_T<T> r,
                            vector3_T<T>& collision_point)
{
	const T denom = dot(plane.normal, r.direction);
	if (abs(denom) > T(0.0001))
	{
		const T t = -dot(plane.origin - r.origin, plane.normal) / denom;

		collision_point = r.origin + r.direction * t;
		return true;
	}
	return false;
}
#pragma endregion

#pragma region definition_ray2_t
// https://rootllama.wordpress.com/2014/06/20/ray-line-segment-intersection-test-in-2d/
template <arithmetic T>
constexpr bool collides(const ray2_T<T>& r, const segment2_T<T>& s)
{
	const vector2_T<T> v1 = r.origin - s.origin;
	const vector2_T<T> v2 = s.end - s.origin;
	const vector2_T<T> v3 = {-r.direction.y(), r.direction.x()};

	const float V2DotV3 = dot(v2, v3);
	const float t1 = cross(v2, v1) / V2DotV3;
	const float t2 = dot(v1, v3) / V2DotV3;

	return t1 >= T(0) && t2 >= T(0) && t2 <= T(1);
}

// https://rootllama.wordpress.com/2014/06/20/ray-line-segment-intersection-test-in-2d/
template <arithmetic T>
constexpr std::optional<vector2_T<T>> intersection(const ray2_T<T>& r,
                                                   const segment2_T<T>& s)
{
	const vector2_T<T> v1 = r.origin - s.origin;
	const vector2_T<T> v2 = s.end - s.origin;
	const vector2_T<T> v3 = {-r.direction.y(), r.direction.x()};

	const float V2DotV3 = dot(v2, v3);
	const float t1 = cross(v2, v1) / V2DotV3;
	const float t2 = dot(v1, v3) / V2DotV3;

	if (t1 >= T(0) && t2 >= T(0) && t2 <= T(1))
		return r.origin + r.direction * t1;
	else
		return std::nullopt;
}
#pragma endregion

#pragma region definition_circle_t
template <arithmetic T>
bool collides(circle_T<T> c, vector2_T<T> p)
{
	return distance_between(c.origin, p) <= c.radius;
}
template <arithmetic T>
bool collides(vector2_T<T> p, circle_T<T> c)
{
	return collides(c, p);
}
#pragma endregion

#pragma region definition_aabb2_t
template <arithmetic T>
constexpr bool aabb2_T<T>::contains(vector2_T<T> p) const
{
	return p.x >= min.x && p.x <= max.x &&  //
	       p.y >= min.y && p.y <= max.y;    //
}

template <arithmetic T>
constexpr vector2_T<T> aabb2_T<T>::get_center() const
{
	return (max + min) / T(2);
}

template <arithmetic T>
constexpr std::array<vector2_T<T>, 4> aabb2_T<T>::get_points() const
{
	return std::array<vector2_T<T>, 4>{
	    vector2_T<T>{min.x, min.y},
	    vector2_T<T>{max.x, min.y},
	    vector2_T<T>{min.x, max.y},
	    vector2_T<T>{max.x, max.y},
	};
}

template <arithmetic T>
constexpr vector2_T<T> aabb2_T<T>::clamp(vector2_T<T> p) const
{
	p = element_wise_max(p, min);
	p = element_wise_min(p, max);
	return p;
}

template <arithmetic T>
constexpr aabb2_T<T>& aabb2_T<T>::grow(const aabb2_T<T>& box)
{
	grow(box.min);
	grow(box.max);
	return *this;
}

template <arithmetic T>
constexpr aabb2_T<T>& aabb2_T<T>::grow(vector2_T<T> point)
{
	min = cdm::element_wise_min(min, point);
	max = cdm::element_wise_max(max, point);
	return *this;
}

template <arithmetic T>
constexpr aabb2_T<T> aabb2_T<T>::operator+(const aabb2_T<T>& rhs) const
{
	aabb2_T<T> res = *this;
	res.grow(rhs);
	return res;
}

template <arithmetic T>
constexpr aabb2_T<T> aabb2_T<T>::operator+(vector2_T<T> rhs) const
{
	aabb2_T<T> res = *this;
	res.grow(rhs);
	return res;
}

template <arithmetic T>
constexpr aabb2_T<T>& aabb2_T<T>::operator+=(const aabb2_T<T>& rhs)
{
	return grow(rhs);
}

template <arithmetic T>
constexpr aabb2_T<T>& aabb2_T<T>::operator+=(vector2_T<T> rhs)
{
	return grow(rhs);
}

template <arithmetic T>
constexpr aabb2_T<T> aabb2_T<T>::infinity()
{
	return aabb2_T<T>{
	    .min{
	        std::numeric_limits<T>::infinity(),
	        std::numeric_limits<T>::infinity(),
	    },
	    .max{
	        -std::numeric_limits<T>::infinity(),
	        -std::numeric_limits<T>::infinity(),
	    },
	};
}

template <arithmetic T>
constexpr aabb2_T<T> operator+(vector2_T<T> lhs, const aabb2_T<T>& rhs)
{
	return rhs + lhs;
}

template <arithmetic T>
constexpr bool collides(const aabb2_T<T>& b0, const aabb2_T<T>& b1)
{
	return collides(b0, b1.origin) ||
	       collides(b0, b1.origin + vector2_T<T>(b1.dimension.x, T(0))) ||
	       collides(b0, b1.origin + vector2_T<T>(T(0), b1.dimension.y)) ||
	       collides(b0, b1.origin + b1.dimension);
}
#pragma endregion

#pragma region definition_aabb3_t
template <arithmetic T>
constexpr bool aabb3_T<T>::contains(vector3_T<T> p) const
{
	return p.x >= min.x && p.x <= max.x &&  //
	       p.y >= min.y && p.y <= max.y &&  //
	       p.z >= min.z && p.z <= max.z;    //
}

template <arithmetic T>
constexpr vector3_T<T> aabb3_T<T>::get_center() const
{
	return (max + min) / T(2);
}

template <arithmetic T>
constexpr std::array<vector3_T<T>, 8> aabb3_T<T>::get_points() const
{
	return std::array<vector3_T<T>, 8>{
	    vector3_T<T>{min.x, min.y, min.z}, vector3_T<T>{max.x, min.y, min.z},
	    vector3_T<T>{min.x, max.y, min.z}, vector3_T<T>{max.x, max.y, min.z},
	    vector3_T<T>{min.x, min.y, max.z}, vector3_T<T>{max.x, min.y, max.z},
	    vector3_T<T>{min.x, max.y, max.z}, vector3_T<T>{max.x, max.y, max.z},
	};
}

template <arithmetic T>
constexpr vector3_T<T> aabb3_T<T>::clamp(vector3_T<T> p) const
{
	p = element_wise_max(p, min);
	p = element_wise_min(p, max);
	return p;
}

template <arithmetic T>
constexpr aabb3_T<T>& aabb3_T<T>::grow(const aabb3_T<T>& box)
{
	grow(box.min);
	grow(box.max);
	return *this;
}

template <arithmetic T>
constexpr aabb3_T<T>& aabb3_T<T>::grow(vector3_T<T> point)
{
	min = cdm::element_wise_min(min, point);
	max = cdm::element_wise_max(max, point);
	return *this;
}

template <arithmetic T>
constexpr aabb3_T<T> aabb3_T<T>::operator+(const aabb3_T<T>& rhs) const
{
	aabb3_T<T> res = *this;
	res.grow(rhs);
	return res;
}

template <arithmetic T>
constexpr aabb3_T<T> aabb3_T<T>::operator+(vector3_T<T> rhs) const
{
	aabb3_T<T> res = *this;
	res.grow(rhs);
	return res;
}

template <arithmetic T>
constexpr aabb3_T<T>& aabb3_T<T>::operator+=(const aabb3_T<T>& rhs)
{
	return grow(rhs);
}

template <arithmetic T>
constexpr aabb3_T<T>& aabb3_T<T>::operator+=(vector3_T<T> rhs)
{
	return grow(rhs);
}

template <arithmetic T>
constexpr aabb3_T<T> aabb3_T<T>::infinity()
{
	return aabb3_T<T>{
	    .min{
	        std::numeric_limits<T>::infinity(),
	        std::numeric_limits<T>::infinity(),
	        std::numeric_limits<T>::infinity(),
	    },
	    .max{
	        -std::numeric_limits<T>::infinity(),
	        -std::numeric_limits<T>::infinity(),
	        -std::numeric_limits<T>::infinity(),
	    },
	};
}

template <arithmetic T>
constexpr aabb3_T<T> operator+(vector3_T<T> lhs, const aabb3_T<T>& rhs)
{
	return rhs + lhs;
}

template <arithmetic T>
constexpr bool collides(const aabb3_T<T>& b, const ray3_T<T>& r)
{
	constexpr vector3_T<T> inv{
	    T(1) / r.direction.x(),  //
	    T(1) / r.direction.y(),  //
	    T(1) / r.direction.z()   //
	};

	T t1 = (b.min.x - r.origin.x) * inv.x;
	T t2 = (b.max.x - r.origin.x) * inv.x;

	T tmin = std::min(t1, t2);
	T tmax = std::max(t1, t2);

	t1 = (b.min.y - r.origin.y) * inv.y;
	t2 = (b.max.y - r.origin.y) * inv.y;

	tmin = std::max(tmin, std::min(t1, t2));
	tmax = std::min(tmax, std::max(t1, t2));

	t1 = (b.min.z - r.origin.z) * inv.z;
	t2 = (b.max.z - r.origin.z) * inv.z;

	tmin = std::max(tmin, std::min(t1, t2));
	tmax = std::min(tmax, std::max(t1, t2));

	return tmax >= tmin;
}

template <arithmetic T>
constexpr bool collides(const aabb3_T<T>& b, const plane_T<T>& p)
{
	int evals = 0;
	for (const auto& point : b.get_points())
	{
		const T eval = p.evaluate(point);
		evals += -int(eval < T(0)) + int(eval > T(0));
	}

	return std::abs(evals) != 8;
}

template <arithmetic T>
constexpr bool collides(const aabb3_T<T>& b, const oriented_plane_T<T>& p)
{
	int evals = 0;
	for (const auto& point : b.get_points())
	{
		const T eval = p.evaluate(point);
		evals += -int(eval < T(0)) + int(eval > T(0));
	}

	return std::abs(evals) != 8;
}

template <arithmetic T>
constexpr bool collides(const aabb3_T<T>& b0, const aabb3_T<T>& b1)
{
	/// TODO this is wrong!
	///
	///         +--------------+
	///         |              |
	///         |              |
	///     +---|--------------|---+
	///     |   |              |   |
	///     |   |              |   |
	///     |   |              |   |
	///     |   |              |   |
	///     +---|--------------|---+
	///         |              |
	///         |              |
	///         +--------------+
	///
	/// ^ this return false
	///

	std::array<vector3_T<T>, 8> points = b0.get_points();

	for (auto& p : points)
		if (b1.contains(p))
			return true;

	points = b1.get_points();

	for (auto& p : points)
		if (b0.contains(p))
			return true;

	return false;
}
#pragma endregion

#pragma region definition_transform2_T
template <arithmetic T>
transform2_T<T> transform2_T<T>::operator*(const transform2_T<T>& t) const
{
	transform2_T<T> res;
	res.position =
	    rotation * vector2_T(scale.x * t.position.x, scale.y * t.position.y) +
	    position;
	res.rotation = rotation + t.rotation;
	res.scale = vector2_T(scale.x * t.scale.x, scale.y * t.scale.y);
	return res;
}
template <arithmetic T>
vector2_T<T> transform2_T<T>::operator*(vector2_T<T> v) const
{
	return rotation * vector2_T(scale.x * v.x, scale.y * v.y) + position;
}

template <arithmetic T>
transform2_T<T> transform2_T<T>::identity()
{
	return {
	    .position = vector2_T<T>::zero(),
	    .rotation = normalized<complex_T<T>>::already_normalized(
	        complex_T<T>::identity()),
	    .scale = vector2_T<T>::one(),
	};
}
#pragma endregion

#pragma region definition_transform3_T
template <arithmetic T>
matrix4_T<T> transform3_T<T>::to_matrix4() const
{
	matrix4_T<T> res{matrix4_T<T>::rotation(rotation)};

	res.column(0).row(0) = res.column(0).row(0) * scale.x;
	res.column(0).row(1) = res.column(0).row(1) * scale.y;
	res.column(0).row(2) = res.column(0).row(2) * scale.z;

	res.column(1).row(0) = res.column(1).row(0) * scale.x;
	res.column(1).row(1) = res.column(1).row(1) * scale.y;
	res.column(1).row(2) = res.column(1).row(2) * scale.z;

	res.column(2).row(0) = res.column(2).row(0) * scale.x;
	res.column(2).row(1) = res.column(2).row(1) * scale.y;
	res.column(2).row(2) = res.column(2).row(2) * scale.z;

	res.column(3).row(0) = position.x;
	res.column(3).row(1) = position.y;
	res.column(3).row(2) = position.z;

	return res;
}
template <arithmetic T>
transform3_T<T>& transform3_T<T>::translate_absolute(vector3_T<T> t)
{
	position += t;
	return *this;
}
template <arithmetic T>
transform3_T<T>& transform3_T<T>::translate_relative(vector3_T<T> t)
{
	position += rotation * t;
	return *this;
}
template <arithmetic T>
transform3_T<T>& transform3_T<T>::rotate(quaternion_T<T> r)
{
	rotation = r * rotation;
	return *this;
}

template <arithmetic T>
transform3_T<T> transform3_T<T>::operator*(const transform3_T<T>& t) const
{
	transform3_T<T> res{*this};
	// res.position = rotation *
	//                   vector3_T{
	//                       scale.x * t.position.x,  //
	//                       scale.y * t.position.y,  //
	//                       scale.z * t.position.z,  //
	//                   } +
	//               position;
	// res.rotation = rotation * t.rotation;
	// res.scale = vector3_T{
	//    scale.x * t.scale.x,  //
	//    scale.y * t.scale.y,  //
	//    scale.z * t.scale.z,  //
	//};

	// res.position = res.position * t.rotation;
	// res.rotation = res.rotation * t.rotation;
	////res.scale = t.rotation * res.scale;

	// res.position.x *= t.scale.x;
	// res.position.y *= t.scale.y;
	// res.position.z *= t.scale.z;
	// res.scale.x *= t.scale.x;
	// res.scale.y *= t.scale.y;
	// res.scale.z *= t.scale.z;

	res.position += t.position;

	return res;
}
template <arithmetic T>
vector3_T<T> transform3_T<T>::operator*(vector3_T<T> v) const
{
	vector3_T<T> res = rotation * v;
	res.x *= scale.x;
	res.y *= scale.y;
	res.z *= scale.z;
	return res + position;
}
template <arithmetic T>
quaternion_T<T> transform3_T<T>::operator*(quaternion_T<T> q) const
{
	return rotation * q;
}

template <arithmetic T>
transform3_T<T> transform3_T<T>::identity()
{
	return {
	    .position = vector3_T<T>::zero(),
	    .rotation = quaternion_T<T>::identity(),
	    .scale = vector3_T<T>::one(),
	};
}
#pragma endregion

#pragma region definition_uniform_transform2_T
template <arithmetic T>
uniform_transform2_T<T> uniform_transform2_T<T>::operator*(
    uniform_transform2_T<T> t) const
{
	uniform_transform2_T<T> res;
	matrix2_T r = matrix2_T<T>::rotation(rotation);
	res.position = r * (scale * t.position) + position;
	res.rotation = rotation + t.rotation;
	res.scale = scale * t.scale;
	return res;
}
template <arithmetic T>
vector2_T<T> uniform_transform2_T<T>::operator*(vector2_T<T> v) const
{
	matrix2_T r = matrix2_T::rotation(rotation);
	return r * (scale * v) + position;
}

template <arithmetic T>
uniform_transform2_T<T> uniform_transform2_T<T>::identity()
{
	return {
	    .position = vector2_T<T>::zero(),
	    .rotation = normalized<complex_T<T>>::already_normalized(
	        complex_T<T>::identity()),
	    .scale = T(1),
	};
}
#pragma endregion

#pragma region definition_uniform_transform3_T
template <arithmetic T>
matrix4_T<T> uniform_transform3_T<T>::to_matrix4() const
{
	matrix4_T<T> res{matrix4_T<T>::rotation(rotation)};

	res.column(0).row(0) = res.column(0).row(0) * scale;
	res.column(0).row(1) = res.column(0).row(1) * scale;
	res.column(0).row(2) = res.column(0).row(2) * scale;

	res.column(1).row(0) = res.column(1).row(0) * scale;
	res.column(1).row(1) = res.column(1).row(1) * scale;
	res.column(1).row(2) = res.column(1).row(2) * scale;

	res.column(2).row(0) = res.column(2).row(0) * scale;
	res.column(2).row(1) = res.column(2).row(1) * scale;
	res.column(2).row(2) = res.column(2).row(2) * scale;

	res.column(3).row(0) = position.x;
	res.column(3).row(1) = position.y;
	res.column(3).row(2) = position.z;

	return res;
}
template <arithmetic T>
uniform_transform3_T<T> uniform_transform3_T<T>::operator*(
    const uniform_transform3_T<T>& t) const
{
	uniform_transform3_T<T> res;
	res.position = rotation * (scale * t.position) + position;
	res.rotation = rotation * t.rotation;
	res.scale = scale * t.scale;
	return res;
}
template <arithmetic T>
vector3_T<T> uniform_transform3_T<T>::operator*(vector3_T<T> v) const
{
	return rotation * (scale * v) + position;
}
template <arithmetic T>
quaternion_T<T> uniform_transform3_T<T>::operator*(quaternion_T<T> q) const
{
	return rotation * q;
}

template <arithmetic T>
uniform_transform3_T<T> uniform_transform3_T<T>::identity()
{
	return {
	    .position = vector3_T<T>::zero(),
	    .rotation = quaternion_T<T>::identity(),
	    .scale = T(1),
	};
}
#pragma endregion

#pragma region definition_unscaled_transform2_T
template <arithmetic T>
unscaled_transform2_T<T> unscaled_transform2_T<T>::operator*(
    unscaled_transform2_T<T> t) const
{
	unscaled_transform2_T<T> res;
	matrix2_T<T> r = matrix2_T<T>::rotation(rotation);
	res.position = r * t.position + position;
	res.rotation = rotation + t.rotation;
	return res;
}
template <arithmetic T>
vector2_T<T> unscaled_transform2_T<T>::operator*(vector2_T<T> v) const
{
	matrix2_T<T> r = matrix2_T::rotation(rotation);
	return r * v + position;
}

template <arithmetic T>
unscaled_transform2_T<T> unscaled_transform2_T<T>::identity()
{
	return {
	    .position = vector2_T<T>::zero(),
	    .rotation = normalized<complex_T<T>>::already_normalized(
	        complex_T<T>::identity()),
	};
}
#pragma endregion

#pragma region definition_unscaled_transform3_T
template <arithmetic T>
unscaled_transform3_T<T>& unscaled_transform3_T<T>::translate_absolute(
    vector3_T<T> t)
{
	position += t;
	return *this;
}
template <arithmetic T>
unscaled_transform3_T<T>& unscaled_transform3_T<T>::translate_relative(
    vector3_T<T> t)
{
	position += rotation * t;
	return *this;
}
template <arithmetic T>
unscaled_transform3_T<T>& unscaled_transform3_T<T>::rotate(quaternion_T<T> r)
{
	rotation = r * rotation;
	return *this;
}

template <arithmetic T>
unscaled_transform3_T<T>& unscaled_transform3_T<T>::inverse()
{
	rotation.inverse();
	position = rotation * -position;
	return *this;
}
template <arithmetic T>
unscaled_transform3_T<T> unscaled_transform3_T<T>::get_inversed() const
{
	unscaled_transform3_T<T> res{*this};
	res.inverse();
	return res;
}

template <arithmetic T>
matrix4_T<T> unscaled_transform3_T<T>::to_matrix4() const
{
	matrix4_T<T> res{matrix4_T<T>::rotation(rotation)};

	res.column(3).row(0) = position.x;
	res.column(3).row(1) = position.y;
	res.column(3).row(2) = position.z;

	return res;
}

template <arithmetic T>
unscaled_transform3_T<T> inverse(unscaled_transform3_T<T> t)
{
	t.rotation = inverse(t.rotation);
	t.position = t.rotation * -t.position;
	return t;
}

template <arithmetic T>
unscaled_transform3_T<T> unscaled_transform3_T<T>::operator*(
    const unscaled_transform3_T<T>& t) const
{
	unscaled_transform3_T<T> res;
	res.position = rotation * t.position + position;
	res.rotation = rotation * t.rotation;
	return res;
}
template <arithmetic T>
vector3_T<T> unscaled_transform3_T<T>::operator*(vector3_T<T> v) const
{
	return rotation * v + position;
}
template <arithmetic T>
quaternion_T<T> unscaled_transform3_T<T>::operator*(quaternion_T<T> q) const
{
	return rotation * q;
}

template <arithmetic T>
unscaled_transform3_T<T> unscaled_transform3_T<T>::identity()
{
	return {
	    .position = vector3_T<T>::zero(),
	    .rotation = quaternion_T<T>::identity(),
	};
}
#pragma endregion

#pragma region definition_value_domain
template <typename T>
    requires arithmetic<T> || vector<T>
constexpr T value_domain_T<T>::lerp(normalized_value_T<T> value) const noexcept
{
	return element_wise_lerp(lim0(), lim1(), value.value());
}
#pragma endregion

#pragma region definition_unnormalized_value
template <typename T>
constexpr unnormalized_value_T<T>::unnormalized_value_T(
    value_domain_T<T> domain,
    normalized_value_T<T> value) noexcept
    : m_domain{domain}, m_value{m_domain.lerp(value)}
{
}

template <typename T>
constexpr normalized_value_T<T> unnormalized_value_T<T>::to_normalized()
    const noexcept
{
	return {*this};
}

template <typename T>
constexpr unnormalized_value_T<T>& unnormalized_value_T<T>::operator=(
    normalized_value_T<T> value) noexcept
{
	m_value = m_domain.lerp(value);
	return *this;
}
#pragma endregion

#pragma region definition_normalized_value
#pragma endregion

#pragma region definition_streams
template <arithmetic T>
std::ostream& operator<<(std::ostream& o, const radian_T<T>& a)
{
	return o << "radian_T(" << static_cast<T>(a) << ")";
}
template <arithmetic T>
std::ostream& operator<<(std::ostream& o, const degree_T<T>& a)
{
	return o << "degree_T(" << static_cast<T>(a) << ")";
}
template <arithmetic T>
std::ostream& operator<<(std::ostream& o, const pi_fraction_T<T>& f)
{
	return o << "pi_fraction_T(" << f.numerator << "pi / " << f.denominator
	         << ")";
}
template <arithmetic T, T NumT, T DenT>
std::ostream& operator<<(std::ostream& o,
                         const static_pi_fraction_T<T, NumT, DenT>& f)
{
	return o << "static_pi_fraction_T(" << f.numerator << "pi / "
	         << f.denominator << ")";
}
template <arithmetic T>
std::ostream& operator<<(std::ostream& os, vector2_T<T> v)
{
	return os << "vector2_T(" << v.x << ", " << v.y << ")";
}
template <arithmetic T>
std::ostream& operator<<(std::ostream& os, vector3_T<T> v)
{
	return os << "vector3_T(" << v.x << ", " << v.y << ", " << v.z << ")";
}
template <arithmetic T>
std::ostream& operator<<(std::ostream& os, vector4_T<T> v)
{
	return os << "vector4_T(" << v.x << ", " << v.y << ", " << v.z << ", "
	          << v.w << ")";
}
template <arithmetic T>
std::ostream& operator<<(std::ostream& os, quaternion_T<T> q)
{
	const auto flags = os.flags();
	const auto w = std::setw(13);
	os << std::right;
	// clang-format off
	os << "quaternion_T({" << w << q.x << ", " << w << q.y << ", " << w << q.z << "}, " << w << q.w << ")";
	// clang-format on
	os.setf(flags);
	return os;
}
template <arithmetic T>
std::ostream& operator<<(std::ostream& os, plane_T<T> p)
{
	return os << "plane_T(origin = " << p.origin << ", normal = " << *p.normal
	          << ")";
}
template <arithmetic T>
std::ostream& operator<<(std::ostream& os, oriented_plane_T<T> p)
{
	// clang-format off
	return os << "oriented_plane_T(origin = " << p.origin
	                         << ", normal = " << *p.normal
	                         << ", tangent = " << *p.tangent
	                         << ")";
	// clang-format on
}
template <arithmetic T>
std::ostream& operator<<(std::ostream& os, const cdm::matrix2_T<T>& m)
{
	const auto flags = os.flags();
	const auto w = std::setw(13);
	os << std::right;
	// clang-format off
	os << "matrix2_T(" << w << m.column(0).row(0) << " " << w << m.column(1).row(0) << "\n"
	      "          " << w << m.column(0).row(1) << " " << w << m.column(1).row(1) << ")";
	// clang-format on
	os.setf(flags);
	return os;
}
template <arithmetic T>
std::ostream& operator<<(std::ostream& os, const cdm::matrix3_T<T>& m)
{
	const auto flags = os.flags();
	const auto w = std::setw(13);
	os << std::right;
	// clang-format off
	os << "matrix3_T(" << w << m.column(0).row(0) << " " << w << m.column(1).row(0) << " " << w << m.column(2).row(0) << "\n"
	      "          " << w << m.column(0).row(1) << " " << w << m.column(1).row(1) << " " << w << m.column(2).row(1) << "\n"
	      "          " << w << m.column(0).row(2) << " " << w << m.column(1).row(2) << " " << w << m.column(2).row(2) << ")";
	// clang-format on
	os.setf(flags);
	return os;
}
template <arithmetic T>
std::ostream& operator<<(std::ostream& os, const cdm::matrix4_T<T>& m)
{
	const auto flags = os.flags();
	const auto w = std::setw(13);
	os << std::right;
	// clang-format off
	os << "matrix4_T(" << w << m.column(0).row(0) << " " << w << m.column(1).row(0) << " " << w << m.column(2).row(0) << " " << w << m.column(3).row(0) << "\n"
	      "          " << w << m.column(0).row(1) << " " << w << m.column(1).row(1) << " " << w << m.column(2).row(1) << " " << w << m.column(3).row(1) << "\n"
	      "          " << w << m.column(0).row(2) << " " << w << m.column(1).row(2) << " " << w << m.column(2).row(2) << " " << w << m.column(3).row(2) << "\n"
	      "          " << w << m.column(0).row(3) << " " << w << m.column(1).row(3) << " " << w << m.column(2).row(3) << " " << w << m.column(3).row(3) << ")";
	// clang-format on
	os.setf(flags);
	return os;
}
template <arithmetic T>
std::ostream& operator<<(std::ostream& os, const cdm::normalized<T>& n)
{
	return os << *n;
}
template <arithmetic T>
std::ostream& operator<<(std::ostream& os, transform3_T<T> t)
{
	// clang-format off
	return os << "transform3_T(position = " << t.position << ",\n"
	          << "             rotation = " << t.rotation << ",\n"
	          << "             scale =    " << t.scale << ")";
	// clang-format on
}
template <arithmetic T>
std::ostream& operator<<(std::ostream& os, segment2_T<T> t)
{
	// clang-format off
	return os << "segment2_T(origin = " << t.origin << ",\n"
	          << "           end =    " << t.end << ")";
	// clang-format on
}
template <arithmetic T>
std::ostream& operator<<(std::ostream& os, segment3_T<T> t)
{
	// clang-format off
	return os << "segment3_T(origin = " << t.origin << ",\n"
	          << "           end =    " << t.end << ")";
	// clang-format on
}
#pragma endregion
#pragma endregion

namespace literals
{
#pragma region definition_literals
inline cdm::radian_T<float> operator""_rad(long double d)
{
	return cdm::radian_T<float>(static_cast<float>(d));
}
inline cdm::radian_T<float> operator""_rad(unsigned long long int i)
{
	return cdm::radian_T<float>(static_cast<float>(i));
}
inline cdm::radian_T<float> operator""_pi(long double d)
{
	return cdm::radian_T<float>(static_cast<float>(d) * cdm::pi);
}
inline cdm::radian_T<float> operator""_pi(unsigned long long int i)
{
	return cdm::radian_T<float>(static_cast<float>(i) * cdm::pi);
}
inline cdm::degree_T<float> operator""_deg(long double d)
{
	return cdm::degree_T<float>(static_cast<float>(d));
}
inline cdm::degree_T<float> operator""_deg(unsigned long long int i)
{
	return cdm::degree_T<float>(static_cast<float>(i));
}
inline cdm::radian_T<double> operator""_radd(long double d)
{
	return cdm::radian_T<double>(static_cast<double>(d));
}
inline cdm::radian_T<double> operator""_radd(unsigned long long int i)
{
	return cdm::radian_T<double>(static_cast<double>(i));
}
inline cdm::radian_T<double> operator""_pid(long double d)
{
	return cdm::radian_T<double>(static_cast<double>(d) * cdm::pi);
}
inline cdm::radian_T<double> operator""_pid(unsigned long long int i)
{
	return cdm::radian_T<double>(static_cast<double>(i) * cdm::pi);
}
inline cdm::degree_T<double> operator""_degd(long double d)
{
	return cdm::degree_T<double>(static_cast<double>(d));
}
inline cdm::degree_T<double> operator""_degd(unsigned long long int i)
{
	return cdm::degree_T<double>(static_cast<double>(i));
}
#pragma endregion
}  // namespace literals

using namespace literals;

using complex = complex_T<float>;
using complexd = complex_T<double>;
using radian = radian_T<float>;
using radiand = radian_T<double>;
using degree = degree_T<float>;
using degreed = degree_T<double>;
using pi_fraction = pi_fraction_T<int32_t>;
template <int32_t NumT, int32_t DenT>
using static_pi_fraction = static_pi_fraction_T<int32_t, NumT, DenT>;
using vector2 = vector2_T<float>;
using vector2d = vector2_T<double>;
using vector3 = vector3_T<float>;
using vector3d = vector3_T<double>;
using vector4 = vector4_T<float>;
using vector4d = vector4_T<double>;
using matrix2 = matrix2_T<float>;
using matrix2d = matrix2_T<double>;
using matrix3 = matrix3_T<float>;
using matrix3d = matrix3_T<double>;
using matrix4 = matrix4_T<float>;
using matrix4d = matrix4_T<double>;
using perspective = perspective_T<float>;
using perspectived = perspective_T<double>;
using direction2 = direction2_T<float>;
using direction2d = direction2_T<double>;
using direction3 = direction3_T<float>;
using direction3d = direction3_T<double>;
using euler_angles = euler_angles_T<float>;
using euler_anglesd = euler_angles_T<double>;
using quaternion = quaternion_T<float>;
using quaterniond = quaternion_T<double>;
template <line_representation representation>
using line = line_T<float, representation>;
template <line_representation representation>
using lined = line_T<double, representation>;
using segment2 = segment2_T<float>;
using segment2d = segment2_T<double>;
using segment3 = segment3_T<float>;
using segment3d = segment3_T<double>;
using plane = plane_T<float>;
using planed = plane_T<double>;
using oriented_plane = oriented_plane_T<float>;
using oriented_planed = oriented_plane_T<double>;
using circle = circle_T<float>;
using circled = circle_T<double>;
using ray2 = ray2_T<float>;
using ray2d = ray2_T<double>;
using ray3 = ray3_T<float>;
using ray3d = ray3_T<double>;
using aabb2 = aabb2_T<float>;
using aabb2d = aabb2_T<double>;
using aabb3 = aabb3_T<float>;
using aabb3d = aabb3_T<double>;
using transform2 = transform2_T<float>;
using transform2d = transform2_T<double>;
using transform3 = transform3_T<float>;
using transform3d = transform3_T<double>;
using uniform_transform2 = uniform_transform2_T<float>;
using uniform_transform2d = uniform_transform2_T<double>;
using uniform_transform3 = uniform_transform3_T<float>;
using uniform_transform3d = uniform_transform3_T<double>;
using unscaled_transform2 = unscaled_transform2_T<float>;
using unscaled_transform2d = unscaled_transform2_T<double>;
using unscaled_transform3 = unscaled_transform3_T<float>;
using unscaled_transform3d = unscaled_transform3_T<double>;
using unnormalized_float = unnormalized_value_T<float>;
using unnormalized_double = unnormalized_value_T<double>;
using normalized_float = normalized_value_T<float>;
using normalized_double = normalized_value_T<double>;
}  // namespace cdm

namespace cdm_literals
{
using namespace cdm::literals;
}

#endif  // CDM_MATHS_HPP
