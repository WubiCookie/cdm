/* cdm_maths - v2.0.0 - geometric library - https://github.com/WubiCookie/cdm
   no warranty implied; use at your own risk

LICENSE

       DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
                   Version 2, December 2004

Copyright (C) 2021 Charles Seizilles de Mazancourt <charles DOT de DOT
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
#define CDM_MATHS_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <ostream>
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
#pragma region forward declarations
template <typename T>
class normalized;
template <typename T>
struct complex_t;
template <typename T>
struct radian_t;
template <typename T>
struct degree_t;
template <typename T>
struct pi_fraction_t;
template <typename T, T NumeratorT, T DenominatorT>
struct static_pi_fraction_t;
template <typename T>
struct vector2_t;
template <typename T>
struct vector3_t;
template <typename T>
struct vector4_t;
template <typename T>
struct matrix2_t;
template <typename T>
struct matrix3_t;
template <typename T>
struct matrix4_t;
template <typename T>
struct matrix4_t;
template <typename T>
struct perspective_t;
template <typename T>
struct euler_angles_t;
template <typename T>
struct quaternion_t;
template <typename T>
struct cartesian_direction2d_t;
template <typename T>
struct polar_direction_t;
template <typename T>
struct line_t;
template <typename T>
struct segment2d_t;
template <typename T>
struct plane_t;
template <typename T>
struct aa_rect_t;
template <typename T>
struct circle_t;
template <typename T>
struct ray2d_t;
template <typename T>
struct ray3d_t;
template <typename T>
struct aabb_t;
template <typename T>
struct transform2d_t;
template <typename T>
struct transform3d_t;
template <typename T>
struct uniform_transform2d_t;
template <typename T>
struct uniform_transform3d_t;
template <typename T>
struct unscaled_transform2d_t;
template <typename T>
struct unscaled_transform3d_t;

template <typename T>
using direction_t = normalized<vector3_t<T>>;

namespace detail
{
template <uint8_t x, uint8_t y, typename T>
T get_quaternion_t_matrix_element(quaternion_t<T> q);
}  // namespace detail

template <typename T>
constexpr T lerp(T begin, T end, T percent);
template <typename T>
constexpr T clamp(T f, T min, T max);

template <typename T>
bool nearly_equal(T f1, T f2, T e = T(epsilon));

template <typename T>
constexpr int sign(T val);
template <typename T>
constexpr int signnum(T val);
#pragma endregion

#pragma region constants declarations
constexpr double pi =
    3.1415926535897932384626433832795028841971693993751058209749445923;
constexpr double deg_to_rad = pi / 180.0;
constexpr double rad_to_deg = 180.0 / pi;
constexpr double sqrt2 =
    1.4142135623730950488016887242096980785696718753769480731766797379;
constexpr double sqrt3 =
    1.7320508075688772935274463415058723669428052538103806280558069794;
constexpr double inv_sqrt2 =
    0.7071067811865475244008443621048490392848359376884740365883398689;
constexpr double inv_sqrt3 =
    0.5773502691896257645091487805019574556476017512701268760186023264;
constexpr double sqrt3_over2 =
    0.8660254037844386467637231707529361834714026269051903140279034897;
constexpr double epsilon = 1.0e-05;
#pragma endregion

#pragma region declaration complex_t
template <typename T>
struct complex_t
{
	T r;
	T i;

	complex_t() = default;
	complex_t(T r_, T i_) : r{r_}, i{i_} {}
	explicit complex_t(radian_t<T> angle);
	complex_t(const complex_t&) = default;
	complex_t(complex_t&&) = default;
	complex_t& operator=(const complex_t&) = default;
	complex_t& operator=(complex_t&&) = default;

	complex_t operator+(complex_t c) const;
	complex_t operator-(complex_t c) const;

	complex_t& operator+=(complex_t c);
	complex_t& operator-=(complex_t c);
	complex_t& operator*=(complex_t c);
};

template <typename T>
complex_t<T> operator*(complex_t<T> c1, complex_t<T> c2);
template <typename T>
complex_t<T> operator*(complex_t<T> c, T f);
template <typename T>
complex_t<T> operator*(normalized<complex_t<T>> c1, complex_t<T> c2);
template <typename T>
complex_t<T> operator*(complex_t<T> c1, normalized<complex_t<T>> c2);
template <typename T>
normalized<complex_t<T>> operator*(normalized<complex_t<T>> c1,
                                   normalized<complex_t<T>> c2);
template <typename T>
vector2_t<T> operator*(normalized<complex_t<T>> c, vector2_t<T> v);

template <typename T>
T norm(complex_t<T> c);
template <typename T>
T norm_squared(complex_t<T> c);
template <typename T>
complex_t<T> normalize(complex_t<T> c);
template <typename T>
complex_t<T> conjugate(complex_t<T> c);
#pragma endregion

#pragma region declaration radian_t
template <typename T>
struct radian_t
{
private:
	T angle;

public:
	radian_t() = default;
	constexpr explicit radian_t(T f);
	radian_t(const radian_t& r) = default;
	radian_t(radian_t&& r) = default;
	constexpr radian_t(degree_t<T> d);

	explicit operator T() const;

	radian_t& operator+=(T f);
	radian_t& operator+=(radian_t r);
	radian_t& operator-=(T f);
	radian_t& operator-=(radian_t r);
	radian_t& operator*=(T f);
	radian_t& operator*=(radian_t r);
	radian_t& operator/=(T f);
	radian_t& operator/=(radian_t r);

	radian_t operator-() const;

	radian_t& operator=(const radian_t& r) = default;
	radian_t& operator=(degree_t<T> d);
};

template <typename T>
radian_t<T> operator+(radian_t<T>, radian_t<T>);
template <typename T>
radian_t<T> operator-(radian_t<T>, radian_t<T>);
template <typename T>
radian_t<T> operator*(radian_t<T>, T);
template <typename T>
radian_t<T> operator*(T, radian_t<T>);
template <typename T>
radian_t<T> operator/(radian_t<T>, T);
template <typename T>
bool operator<(radian_t<T>, radian_t<T>);
template <typename T>
bool operator>(radian_t<T>, radian_t<T>);
template <typename T>
bool operator==(radian_t<T>, radian_t<T>);
template <typename T>
bool operator!=(radian_t<T>, radian_t<T>);
template <typename T>
bool operator>=(radian_t<T>, radian_t<T>);
template <typename T>
bool operator<=(radian_t<T>, radian_t<T>);

template <typename T>
T sin(radian_t<T> r);
template <typename T>
T cos(radian_t<T> r);
template <typename T>
T tan(radian_t<T> r);
template <typename T>
T asin(radian_t<T> r);
template <typename T>
T acos(radian_t<T> r);
template <typename T>
T atan(radian_t<T> r);
template <typename T>
T sinh(radian_t<T> r);
template <typename T>
T cosh(radian_t<T> r);
template <typename T>
T tanh(radian_t<T> r);
template <typename T>
T asinh(radian_t<T> r);
template <typename T>
T acosh(radian_t<T> r);
template <typename T>
T atanh(radian_t<T> r);
#pragma endregion

#pragma region declaration degree_t
template <typename T>
struct degree_t
{
private:
	T angle;

public:
	degree_t() = default;
	constexpr explicit degree_t(T f);
	degree_t(const degree_t& d) = default;
	degree_t(degree_t&& d) = default;
	constexpr degree_t(radian_t<T> r);

	explicit operator T() const;

	degree_t& operator+=(T f);
	degree_t& operator+=(degree_t d);
	degree_t& operator-=(T f);
	degree_t& operator-=(degree_t d);
	degree_t& operator*=(T f);
	degree_t& operator*=(degree_t d);
	degree_t& operator/=(T f);
	degree_t& operator/=(degree_t d);

	degree_t operator-() const;

	degree_t& operator=(const degree_t& d) = default;
	degree_t& operator=(radian_t<T> r);
};

template <typename T>
degree_t<T> operator+(degree_t<T>, degree_t<T>);
template <typename T>
degree_t<T> operator-(degree_t<T>, degree_t<T>);
template <typename T>
degree_t<T> operator*(degree_t<T>, T);
template <typename T>
degree_t<T> operator*(T, degree_t<T>);
template <typename T>
degree_t<T> operator/(degree_t<T>, T);
template <typename T>
bool operator<(degree_t<T>, degree_t<T>);
template <typename T>
bool operator>(degree_t<T>, degree_t<T>);
template <typename T>
bool operator==(degree_t<T>, degree_t<T>);
template <typename T>
bool operator!=(degree_t<T>, degree_t<T>);
template <typename T>
bool operator>=(degree_t<T>, degree_t<T>);
template <typename T>
bool operator<=(degree_t<T>, degree_t<T>);

template <typename T>
T sin(degree_t<T> d);
template <typename T>
T cos(degree_t<T> d);
template <typename T>
T tan(degree_t<T> d);
template <typename T>
T asin(degree_t<T> d);
template <typename T>
T acos(degree_t<T> d);
template <typename T>
T atan(degree_t<T> d);
template <typename T>
T sinh(degree_t<T> d);
template <typename T>
T cosh(degree_t<T> d);
template <typename T>
T tanh(degree_t<T> d);
template <typename T>
T asinh(degree_t<T> d);
template <typename T>
T acosh(degree_t<T> d);
template <typename T>
T atanh(degree_t<T> d);
#pragma endregion

#pragma region declaration pi_fraction_t
template <typename T>
struct pi_fraction_t
{
	static_assert(std::is_arithmetic_v<T>, "T must be arithmetic");
	static_assert(std::is_signed_v<T>, "T must be signed");
	static_assert(std::is_integral_v<T>, "T must be integral");

	T numerator;
	T denominator;

	pi_fraction_t() = default;
	pi_fraction_t(T num, T den) : numerator{num}, denominator{den} {}
	pi_fraction_t(const pi_fraction_t& d) = default;
	pi_fraction_t(pi_fraction_t&& d) = default;

	pi_fraction_t& operator=(const pi_fraction_t&) = default;
	pi_fraction_t& operator=(pi_fraction_t&&) = default;

	template <typename U>
	operator radian_t<U>() const
	{
		return radian_t<U>((U(numerator) * U(pi)) / U(denominator));
	}
	template <typename U>
	operator degree_t<U>() const
	{
		return operator radian_t<U>();
	}

	pi_fraction_t operator-() const
	{
		return pi_fraction_t{-numerator, denominator};
	}
};

template <typename U, typename T>
U sin(pi_fraction_t<T> d)
{
	if constexpr (d.numerator == T(0))
		return U(0);
	else if constexpr (d == static_pi_fraction_t<T, T(1), T(1)>{})
		return U(0);
	else if constexpr (d == static_pi_fraction_t<T, T(1), T(2)>{})
		return U(1);
	else if constexpr (d == static_pi_fraction_t<T, T(1), T(3)>{})
		return U(sqrt3_over2);
	else if constexpr (d == static_pi_fraction_t<T, T(1), T(4)>{})
		return U(inv_sqrt2);
	else if constexpr (d == static_pi_fraction_t<T, T(1), T(6)>{})
		return U(0.5);

	else if constexpr (d == static_pi_fraction_t<T, T(-1), T(1)>{})
		return U(-0);
	else if constexpr (d == static_pi_fraction_t<T, T(-1), T(2)>{})
		return U(-1);
	else if constexpr (d == static_pi_fraction_t<T, T(-1), T(3)>{})
		return U(-sqrt3_over2);
	else if constexpr (d == static_pi_fraction_t<T, T(-1), T(4)>{})
		return U(-inv_sqrt2);
	else if constexpr (d == static_pi_fraction_t<T, T(-1), T(6)>{})
		return U(-0.5);

	else if constexpr (d == static_pi_fraction_t<T, T(2), T(1)>{})
		return U(0);
	else if constexpr (d == static_pi_fraction_t<T, T(2), T(3)>{})
		return U(sqrt3_over2);

	else if constexpr (d == static_pi_fraction_t<T, T(-2), T(1)>{})
		return U(-0);
	else if constexpr (d == static_pi_fraction_t<T, T(-2), T(3)>{})
		return U(-sqrt3_over2);

	else if constexpr (d == static_pi_fraction_t<T, T(3), T(2)>{})
		return U(-1);
	else if constexpr (d == static_pi_fraction_t<T, T(3), T(4)>{})
		return U(inv_sqrt2);

	else if constexpr (d == static_pi_fraction_t<T, T(-3), T(2)>{})
		return U(1);
	else if constexpr (d == static_pi_fraction_t<T, T(-3), T(4)>{})
		return U(-inv_sqrt2);

	constexpr radian_t<U> r = d;
	return sin(r);
}
template <typename U, typename T>
U cos(pi_fraction_t<T> d)
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
			return U(sqrt3_over2);
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

	constexpr radian_t<U> r = d;
	return cos(r);
}
template <typename U, typename T>
U tan(pi_fraction_t<T> d)
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

	constexpr radian_t<U> r = d;
	return tan(r);
}
template <typename U, typename T>
U asin(pi_fraction_t<T> d)
{
	/// TODO: implement
	constexpr radian_t<U> r = d;
	return asin(r);
}
template <typename U, typename T>
U acos(pi_fraction_t<T> d)
{
	/// TODO: implement
	return acos(radian_t<U>(d));
}
template <typename U, typename T>
U atan(pi_fraction_t<T> d)
{
	/// TODO: implement
	return atan(radian_t<U>(d));
}
template <typename U, typename T>
U sinh(pi_fraction_t<T> d)
{
	/// TODO: implement
	constexpr radian_t<U> r = d;
	return sinh(r);
}
template <typename U, typename T>
U cosh(pi_fraction_t<T> d)
{
	/// TODO: implement
	constexpr radian_t<U> r = d;
	return cosh(r);
}
template <typename U, typename T>
U tanh(pi_fraction_t<T> d)
{
	/// TODO: implement
	constexpr radian_t<U> r = d;
	return tanh(r);
}
template <typename U, typename T>
U asinh(pi_fraction_t<T> d)
{
	/// TODO: implement
	constexpr radian_t<U> r = d;
	return asinh(r);
}
template <typename U, typename T>
U acosh(pi_fraction_t<T> d)
{
	/// TODO: implement
	constexpr radian_t<U> r = d;
	return acosh(r);
}
template <typename U, typename T>
U atanh(pi_fraction_t<T> d)
{
	/// TODO: implement
	constexpr radian_t<U> r = d;
	return atanh(r);
}
#pragma endregion

#pragma region declaration static_pi_fraction_t
template <typename T, T NumeratorT, T DenominatorT>
struct static_pi_fraction_t
{
	static_assert(std::is_arithmetic_v<T>, "T must be arithmetic");
	static_assert(std::is_signed_v<T>, "T must be signed");
	static_assert(std::is_integral_v<T>, "T must be integral");
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

	template <typename U>
	constexpr operator radian_t<U>() const
	{
		if constexpr (numerator == T(0))
			return radian_t<U>(U(0));
		else
			return radian_t<U>((U(numerator) * U(pi)) / U(denominator));
	}
	template <typename U>
	constexpr operator degree_t<U>() const
	{
		if constexpr (numerator == T(0))
			return degree_t<U>(U(0));
		else
			return degree_t<U>(
			    radian_t<U>((U(numerator) * U(pi)) / U(denominator)));
	}
	constexpr operator pi_fraction_t<T>() const
	{
		return {numerator, denominator};
	}

	static_pi_fraction_t<T, -numerator, denominator> operator-() const
	{
		return {};
	}
};

template <typename T,
          T NumeratorTL,
          T DenominatorTL,
          T NumeratorTR,
          T DenominatorTR>
constexpr bool operator==(
    const static_pi_fraction_t<T, NumeratorTL, DenominatorTL>& lhs,
    const static_pi_fraction_t<T, NumeratorTR, DenominatorTR>& rhs)
{
	return lhs.numerator == rhs.numerator &&
	       lhs.denominator == rhs.denominator;
}
template <typename T,
          T NumeratorTL,
          T DenominatorTL,
          T NumeratorTR,
          T DenominatorTR>
constexpr bool operator!=(
    const static_pi_fraction_t<T, NumeratorTL, DenominatorTL>& lhs,
    const static_pi_fraction_t<T, NumeratorTR, DenominatorTR>& rhs)
{
	return !(lhs == rhs);
}

template <typename U, typename T, T NumeratorT, T DenominatorT>
U sin(static_pi_fraction_t<T, NumeratorT, DenominatorT> d)
{
	if constexpr (d.numerator == T(0))
		return U(0);
	else if constexpr (d == static_pi_fraction_t<T, T(1), T(1)>{})
		return U(0);
	else if constexpr (d == static_pi_fraction_t<T, T(1), T(2)>{})
		return U(1);
	else if constexpr (d == static_pi_fraction_t<T, T(1), T(3)>{})
		return U(sqrt3_over2);
	else if constexpr (d == static_pi_fraction_t<T, T(1), T(4)>{})
		return U(inv_sqrt2);
	else if constexpr (d == static_pi_fraction_t<T, T(1), T(6)>{})
		return U(0.5);

	else if constexpr (d == static_pi_fraction_t<T, T(-1), T(1)>{})
		return U(-0);
	else if constexpr (d == static_pi_fraction_t<T, T(-1), T(2)>{})
		return U(-1);
	else if constexpr (d == static_pi_fraction_t<T, T(-1), T(3)>{})
		return U(-sqrt3_over2);
	else if constexpr (d == static_pi_fraction_t<T, T(-1), T(4)>{})
		return U(-inv_sqrt2);
	else if constexpr (d == static_pi_fraction_t<T, T(-1), T(6)>{})
		return U(-0.5);

	else if constexpr (d == static_pi_fraction_t<T, T(2), T(1)>{})
		return U(0);
	else if constexpr (d == static_pi_fraction_t<T, T(2), T(3)>{})
		return U(sqrt3_over2);

	else if constexpr (d == static_pi_fraction_t<T, T(-2), T(1)>{})
		return U(-0);
	else if constexpr (d == static_pi_fraction_t<T, T(-2), T(3)>{})
		return U(-sqrt3_over2);

	else if constexpr (d == static_pi_fraction_t<T, T(3), T(2)>{})
		return U(-1);
	else if constexpr (d == static_pi_fraction_t<T, T(3), T(4)>{})
		return U(inv_sqrt2);

	else if constexpr (d == static_pi_fraction_t<T, T(-3), T(2)>{})
		return U(1);
	else if constexpr (d == static_pi_fraction_t<T, T(-3), T(4)>{})
		return U(-inv_sqrt2);

	constexpr radian_t<U> r = d;
	return sin(r);
}
template <typename U, typename T, T NumeratorT, T DenominatorT>
U cos(static_pi_fraction_t<T, NumeratorT, DenominatorT> d)
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
			return U(sqrt3_over2);
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

	constexpr radian_t<U> r = d;
	return cos(r);
}
template <typename U, typename T, T NumeratorT, T DenominatorT>
U tan(static_pi_fraction_t<T, NumeratorT, DenominatorT> d)
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

	constexpr radian_t<U> r = d;
	return tan(r);
}
template <typename U, typename T, T NumeratorT, T DenominatorT>
U asin(static_pi_fraction_t<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	constexpr radian_t<U> r = d;
	return asin(r);
}
template <typename U, typename T, T NumeratorT, T DenominatorT>
U acos(static_pi_fraction_t<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	return acos(radian_t<U>(d));
}
template <typename U, typename T, T NumeratorT, T DenominatorT>
U atan(static_pi_fraction_t<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	return atan(radian_t<U>(d));
}
template <typename U, typename T, T NumeratorT, T DenominatorT>
U sinh(static_pi_fraction_t<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	constexpr radian_t<U> r = d;
	return sinh(r);
}
template <typename U, typename T, T NumeratorT, T DenominatorT>
U cosh(static_pi_fraction_t<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	constexpr radian_t<U> r = d;
	return cosh(r);
}
template <typename U, typename T, T NumeratorT, T DenominatorT>
U tanh(static_pi_fraction_t<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	constexpr radian_t<U> r = d;
	return tanh(r);
}
template <typename U, typename T, T NumeratorT, T DenominatorT>
U asinh(static_pi_fraction_t<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	constexpr radian_t<U> r = d;
	return asinh(r);
}
template <typename U, typename T, T NumeratorT, T DenominatorT>
U acosh(static_pi_fraction_t<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	constexpr radian_t<U> r = d;
	return acosh(r);
}
template <typename U, typename T, T NumeratorT, T DenominatorT>
U atanh(static_pi_fraction_t<T, NumeratorT, DenominatorT> d)
{
	/// TODO: implement
	constexpr radian_t<U> r = d;
	return atanh(r);
}
#pragma endregion

#pragma region declaration vector2_t
template <typename T>
struct vector2_t
{
	T x, y;

	vector2_t() = default;
	vector2_t(T x_, T y_) : x{x_}, y{y_} {}
	vector2_t(std::array<T, 2> a) : vector2_t(a[0], a[1]) {}

	template <typename U>
	explicit operator vector2_t<U>() const
	{
		return vector2_t<U>{static_cast<U>(x), static_cast<U>(y)};
	}

	operator std::array<T, 2>() const { return {x, y}; }

	template <typename U = T>
	std::array<U, 2> to_array() const
	{
		return {static_cast<U>(x), static_cast<U>(y)};
	}

	T norm() const { return std::sqrt(norm_squared()); }
	T norm_squared() const { return x * x + y * y; }
	vector2_t& normalize()
	{
		T n{norm()};
		x /= n;
		y /= n;
		return *this;
	}
	vector2_t get_normalized() const
	{
		vector2_t res{*this};
		res.normalize();
		return res;
	}
	vector2_t& clamp(vector2_t min, vector2_t max)
	{
		x = cdm::clamp(x, min.x, max.x);
		y = cdm::clamp(y, min.y, max.y);
		return *this;
	}
	vector2_t get_clamped(vector2_t min, vector2_t max) const
	{
		vector2_t res{*this};
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

	vector2_t operator+(vector2_t v) const;
	vector2_t operator-(vector2_t v) const;
	vector2_t operator*(T f) const;
	vector2_t operator/(T f) const;

	vector2_t& operator+=(vector2_t v);
	vector2_t& operator-=(vector2_t v);
	vector2_t& operator*=(T f);
	vector2_t& operator/=(T f);

	vector2_t operator-() const;

	bool operator==(vector2_t v) const;
	bool operator!=(vector2_t v) const;
};

template <typename T>
vector2_t<T> operator*(T f, vector2_t<T> v);

template <typename T>
T dot(vector2_t<T> lhs, vector2_t<T> rhs);
template <typename T>
T cross(vector2_t<T> lhs, vector2_t<T> rhs);
template <typename T>
vector2_t<T> lerp(vector2_t<T> begin, vector2_t<T> end, T percent);
template <typename T>
vector2_t<T> nlerp(vector2_t<T> begin, vector2_t<T> end, T percent);
template <typename T>
vector2_t<T> slerp(vector2_t<T> begin, vector2_t<T> end, T percent);
template <typename T>
T distance_between(vector2_t<T> v1, vector2_t<T> v2);
template <typename T>
T distance_squared_between(vector2_t<T> v1, vector2_t<T> v2);
template <typename T>
vector2_t<T> from_to(vector2_t<T> from, vector2_t<T> to);
template <typename T>
radian_t<T> angle_between(vector2_t<T> v1, vector2_t<T> v2);
template <typename T>
bool nearly_equal(vector2_t<T> v1, vector2_t<T> v2, T e = T(epsilon));
template <typename T>
vector2_t<T> element_wise_min(vector2_t<T> v0, vector2_t<T> v1);
template <typename T>
vector2_t<T> element_wise_max(vector2_t<T> v0, vector2_t<T> v1);
#pragma endregion

#pragma region declaration vector3_t
template <typename T>
struct vector3_t
{
	T x, y, z;

	vector3_t() = default;
	vector3_t(T x_, T y_, T z_) : x{x_}, y{y_}, z{z_} {}
	vector3_t(vector2_t<T> v, T z_) : vector3_t(v.x, v.y, z_) {}
	vector3_t(std::array<T, 3> a) : vector3_t(a[0], a[1], a[2]) {}

	template <typename U>
	explicit operator vector3_t<U>() const
	{
		return vector3_t<U>{static_cast<U>(x), static_cast<U>(y),
		                    static_cast<U>(z)};
	}

	operator std::array<T, 3>() const { return {x, y, z}; }

	template <typename U = T>
	std::array<U, 3> to_array() const
	{
		return {static_cast<U>(x), static_cast<U>(y), static_cast<U>(z)};
	}

	vector2_t<T> xy() const;

	T norm() const { return std::sqrt(norm_squared()); }
	T norm_squared() const { return x * x + y * y + z * z; }
	vector3_t& normalize()
	{
		T n = norm();
		x /= n;
		y /= n;
		z /= n;
		return *this;
	}
	vector3_t get_normalized() const
	{
		vector3_t res{*this};
		res.normalize();
		return res;
	}
	vector3_t& clamp(vector3_t min, vector3_t max)
	{
		x = cdm::clamp(x, min.x, max.x);
		y = cdm::clamp(y, min.y, max.y);
		z = cdm::clamp(z, min.z, max.z);
		return *this;
	}
	vector3_t get_clamped(vector3_t min, vector3_t max) const
	{
		vector3_t res{*this};
		res.clamp(min, max);
		return res;
	}

	vector3_t operator+(vector3_t v) const;
	vector3_t operator-(vector3_t v) const;
	vector3_t operator*(T f) const;
	vector3_t operator/(T f) const;

	vector3_t& operator+=(vector3_t v);
	vector3_t& operator-=(vector3_t v);
	vector3_t& operator*=(T f);
	vector3_t& operator/=(T f);

	vector3_t operator-() const;

	bool operator==(vector3_t v) const;
	bool operator!=(vector3_t v) const;
};

template <typename T>
vector3_t<T> operator*(T f, vector3_t<T> v);

template <typename T>
T norm(vector3_t<T> v);
template <typename T>
T norm_squared(vector3_t<T> v);
template <typename T>
vector3_t<T> normalize(vector3_t<T> v);
template <typename T>
vector3_t<T> clamp(vector3_t<T> v, vector3_t<T> min, vector3_t<T> max);
template <typename T>
T dot(vector3_t<T> lhs, vector3_t<T> rhs);
template <typename T>
vector3_t<T> cross(vector3_t<T> lhs, vector3_t<T> rhs);
template <typename T>
vector3_t<T> lerp(vector3_t<T> begin, vector3_t<T> end, T percent);
template <typename T>
vector3_t<T> nlerp(vector3_t<T> begin, vector3_t<T> end, T percent);
template <typename T>
T distance_between(vector3_t<T> v1, vector3_t<T> v2);
template <typename T>
T distance_squared_between(vector3_t<T> v1, vector3_t<T> v2);
template <typename T>
vector3_t<T> from_to(vector3_t<T> from, vector3_t<T> to);
template <typename T>
radian_t<T> angle_between(vector3_t<T> v1, vector3_t<T> v2);
template <typename T>
bool nearly_equal(vector3_t<T> v1, vector3_t<T> v2, T e = epsilon);
template <typename T>
vector3_t<T> element_wise_min(vector3_t<T> v0, vector3_t<T> v1);
template <typename T>
vector3_t<T> element_wise_max(vector3_t<T> v0, vector3_t<T> v1);
template <typename T>
radian_t<T> angle_around_axis(vector3_t<T> v0,
                              vector3_t<T> v1,
                              direction_t<T> axis);
#pragma endregion

#pragma region declaration vector4_t
template <typename T>
struct vector4_t
{
	T x, y, z, w;

	vector4_t() = default;
	vector4_t(T x_, T y_, T z_, T w_) : x{x_}, y{y_}, z{z_}, w{w_} {}
	vector4_t(vector2_t<T> v, T z_, T w_) : vector4_t(v.x, v.y, z_, w_) {}
	vector4_t(vector3_t<T> v, T w_) : vector4_t(v.x, v.y, v.z, w_) {}
	vector4_t(std::array<T, 4>&) : vector4_t(a[0], a[1], a[2], a[3]) {}

	template <typename U>
	explicit operator vector4_t<U>() const
	{
		return vector4_t<U>{static_cast<U>(x), static_cast<U>(y),
		                    static_cast<U>(z), static_cast<U>(w)};
	}

	operator std::array<T, 4>() const { return {x, y, z, w}; }

	template <typename U = T>
	std::array<U, 4> to_array() const
	{
		return {static_cast<U>(x), static_cast<U>(y), static_cast<U>(z),
		        static_cast<U>(w)};
	}

	vector2_t<T> xy() const;
	vector3_t<T> xyz() const;

	vector4_t operator+(vector4_t v) const;
	vector4_t operator-(vector4_t v) const;
	vector4_t operator*(T f) const;
	vector4_t operator/(T f) const;

	vector4_t& operator+=(vector4_t v);
	vector4_t& operator-=(vector4_t v);
	vector4_t& operator*=(T f);
	vector4_t& operator/=(T f);

	vector4_t operator-() const;

	bool operator==(vector4_t v) const;
	bool operator!=(vector4_t v) const;
};

template <typename T>
T dot(vector4_t<T> lhs, vector4_t<T> rhs);

template <typename T>
T norm(vector4_t<T> v);
template <typename T>
T norm_squared(vector4_t<T> v);
template <typename T>
vector4_t<T> normalize(vector4_t<T> v);
template <typename T>
vector4_t<T> clamp(vector4_t<T> v, vector4_t<T> min, vector4_t<T> max);
template <typename T>
T dot(vector4_t<T> lhs, vector4_t<T> rhs);
template <typename T>
vector4_t<T> lerp(vector4_t<T> begin, vector4_t<T> end, T percent);
template <typename T>
vector4_t<T> nlerp(vector4_t<T> begin, vector4_t<T> end, T percent);
template <typename T>
T distance_between(vector4_t<T> v1, vector4_t<T> v2);
template <typename T>
T distance_squared_between(vector4_t<T> v1, vector4_t<T> v2);
template <typename T>
vector4_t<T> from_to(vector4_t<T> from, vector4_t<T> to);
template <typename T>
bool nearly_equal(vector4_t<T> v1, vector4_t<T> v2, T e = epsilon);
template <typename T>
vector4_t<T> element_wise_min(vector4_t<T> v0, vector4_t<T> v1);
template <typename T>
vector4_t<T> element_wise_max(vector4_t<T> v0, vector4_t<T> v1);
#pragma endregion

#pragma region declaration normalized
template <typename T>
class normalized
{
	T vector;

public:
	normalized() = default;
	normalized(const normalized&) = default;
	normalized(normalized&&) = default;
	normalized(const T& t);

	static normalized already_normalized(const T& t);
	static normalized already_normalized(T&& t) noexcept;

	const T& operator*() const;
	const T* operator->() const;
	operator const T&() const;

	normalized operator+(const T& v) const;
	normalized operator-(const T& v) const;
	T operator*(T f) const;
	T operator/(T f) const;

	normalized& operator+=(const T& v);
	normalized& operator-=(const T& v);

	normalized operator-() const;

	bool operator==(const T& v) const;
	bool operator!=(const T& v) const;

	normalized& operator=(const normalized&) = default;
	normalized& operator=(normalized&&) = default;
	normalized& operator=(const T& t);
	normalized& operator=(T&& t) noexcept;
};

template <typename VecT>
VecT cross(normalized<VecT> lhs, VecT rhs)
{
	return cross(*lhs, rhs);
}
template <typename VecT>
VecT cross(VecT lhs, normalized<VecT> rhs)
{
	return cross(lhs, *rhs);
}
template <typename VecT>
normalized<VecT> cross(normalized<VecT> lhs, normalized<VecT> rhs)
{
	return normalized<VecT>::already_normalized(cross(*lhs, *rhs));
}
#pragma endregion

#pragma region declaration matrix2_t
template <typename T>
struct matrix2_t
{
	template <typename U>
	friend class matrix2_t;
	template <typename U>
	friend class matrix3_t;
	template <typename U>
	friend class matrix4_t;

private:
	T m00, m01, m10, m11;

	matrix2_t(T m00, T m10, T m01, T m11);

public:
	matrix2_t() = default;
	explicit matrix2_t(const std::array<T, 4>& a);
	matrix2_t(const matrix2_t&) = default;
	matrix2_t(matrix2_t&&) = default;

	matrix2_t& operator=(const matrix2_t&) = default;
	matrix2_t& operator=(matrix2_t&&) = default;

	template <typename U = T>
	std::array<U, 4> to_array() const
	{
		return {static_cast<U>(m00), static_cast<U>(m10), static_cast<U>(m01),
		        static_cast<U>(m11)};
	}

	template <typename U>
	explicit operator matrix2_t<U>() const
	{
		return {U(m00), U(m10), U(m01), U(m11)};
	}

	static matrix2_t zero();
	static matrix2_t identity();
	static matrix2_t rows(vector2_t<T> row0, vector2_t<T> row1);
	static matrix2_t columns(vector2_t<T> column0, vector2_t<T> column1);
	static matrix2_t rotation(radian_t<T> angle);
	static matrix2_t rotation(normalized<complex_t<T>> angle);
	matrix2_t& transpose();
	matrix2_t get_transposed() const;
	T determinant() const;

	struct proxy
	{
	protected:
		std::reference_wrapper<T> x;
		std::reference_wrapper<T> y;

		constexpr T& vector(int i)
		{
			return std::array<std::reference_wrapper<T>, 2>{x, y}[i];
		}
		constexpr const T& vector(int i) const
		{
			return std::array<std::reference_wrapper<const T>, 2>{x, y}[i];
		}

	public:
		constexpr proxy(T& e0, T& e1) : x{e0}, y{e1} {}
		constexpr proxy(const proxy&) = default;
		constexpr proxy(proxy&&) = default;
		constexpr proxy& operator=(const proxy&) = default;
		constexpr proxy& operator=(proxy&&) = default;
		constexpr proxy& operator=(vector2_t<T> v)
		{
			x = v.x;
			y = v.y;
			return *this;
		}

		constexpr operator vector2_t<T>() const { return vector2_t<T>{x, y}; }
	};
	struct column_proxy : proxy
	{
		using proxy::proxy;
		constexpr T& row(int i) { return proxy::vector(i); }
		constexpr const T& row(int i) const { return proxy::vector(i); }
	};
	struct row_proxy : proxy
	{
		using proxy::proxy;
		constexpr T& column(int i) { return proxy::vector(i); }
		constexpr const T& column(int i) const { return proxy::vector(i); }
	};

	constexpr column_proxy column(int i)
	{
		return std::array{column_proxy(m00, m01), column_proxy(m10, m11)}[i];
	}
	constexpr row_proxy row(int i)
	{
		return std::array{row_proxy(m00, m10), row_proxy(m01, m11)}[i];
	}

	vector2_t<T> operator*(vector2_t<T> v) const;
	matrix2_t operator*(matrix2_t m) const;
};

template <typename T>
matrix2_t<T> transpose(matrix2_t<T> m);
#pragma endregion

#pragma region declaration matrix3_t
template <typename T>
struct matrix3_t
{
	template <typename U>
	friend class matrix2_t;
	template <typename U>
	friend class matrix3_t;
	template <typename U>
	friend class matrix4_t;

private:
	T m00, m01, m02, m10, m11, m12, m20, m21, m22;

	matrix3_t(T m00, T m10, T m20, T m01, T m11, T m21, T m02, T m12, T m22);

public:
	matrix3_t() = default;
	explicit matrix3_t(matrix2_t<T> m);
	matrix3_t(const std::array<T, 9>& a);
	matrix3_t(const matrix3_t&) = default;
	matrix3_t(matrix3_t&&) = default;

	matrix3_t& operator=(const matrix3_t&) = default;
	matrix3_t& operator=(matrix3_t&&) = default;

	template <typename U>
	explicit operator matrix3_t<U>() const
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

	template <typename U = T>
	std::array<U, 9> to_array() const
	{
		return {static_cast<U>(m00), static_cast<U>(m10), static_cast<U>(m20),
		        static_cast<U>(m01), static_cast<U>(m11), static_cast<U>(m21),
		        static_cast<U>(m02), static_cast<U>(m12), static_cast<U>(m22)};
	}

	static matrix3_t zero();
	static matrix3_t identity();
	static matrix3_t rows(vector3_t<T> row0,
	                      vector3_t<T> row1,
	                      vector3_t<T> row2);
	static matrix3_t columns(vector3_t<T> column0,
	                         vector3_t<T> column1,
	                         vector3_t<T> column2);

private:
	static matrix3_t rotation(direction_t<T> axis, T sinAngle, T cosAngle);

public:
	static matrix3_t rotation(direction_t<T> axis, radian_t<T> angle);
	static matrix3_t rotation(direction_t<T> axis,
	                          normalized<complex_t<T>> angle);
	template <typename U, U NumeratorT, U DenominatorT>
	static matrix3_t rotation(
	    direction_t<T> axis,
	    static_pi_fraction_t<U, NumeratorT, DenominatorT> angle);
	static matrix3_t rotation(euler_angles_t<T> r);
	static matrix3_t rotation(quaternion_t<T> q);

private:
	static matrix3_t rotation_around_x(T sinAngle, T cosAngle);
	static matrix3_t rotation_around_y(T sinAngle, T cosAngle);
	static matrix3_t rotation_around_z(T sinAngle, T cosAngle);

public:
	static matrix3_t rotation_around_x(radian_t<T> angle);
	static matrix3_t rotation_around_y(radian_t<T> angle);
	static matrix3_t rotation_around_z(radian_t<T> angle);
	static matrix3_t rotation_around_x(normalized<complex_t<T>> angle);
	static matrix3_t rotation_around_y(normalized<complex_t<T>> angle);
	static matrix3_t rotation_around_z(normalized<complex_t<T>> angle);
	template <typename U, U NumeratorT, U DenominatorT>
	static matrix3_t rotation_around_x(
	    static_pi_fraction_t<U, NumeratorT, DenominatorT> angle);
	template <typename U, U NumeratorT, U DenominatorT>
	static matrix3_t rotation_around_y(
	    static_pi_fraction_t<U, NumeratorT, DenominatorT> angle);
	template <typename U, U NumeratorT, U DenominatorT>
	static matrix3_t rotation_around_z(
	    static_pi_fraction_t<U, NumeratorT, DenominatorT> angle);
	static matrix3_t scale(vector3_t<T> t);
	static matrix3_t scale(T x, T y, T z);
	static matrix3_t scale(T s);
	matrix3_t& inverse();
	matrix3_t get_inversed() const;
	matrix3_t& transpose();
	matrix3_t get_transposed() const;
	T determinant() const;
	bool is_orthogonal() const;

	template <bool IsConstT>
	struct proxy
	{
	protected:
		using Type = std::conditional_t<IsConstT, const T, T>;

		std::reference_wrapper<Type> x;
		std::reference_wrapper<Type> y;
		std::reference_wrapper<Type> z;

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
		constexpr proxy& operator=(vector4_t<T> v)
		{
			x = v.x;
			y = v.y;
			z = v.z;
			return *this;
		}

		constexpr vector2_t<T> xy() const { return vector2_t<T>{x, y}; }
		constexpr operator vector3_t<T>() const
		{
			return vector3_t<T>{x, y, z};
		}
	};
	template <bool IsConstT>
	struct column_proxy : proxy<IsConstT>
	{
		using proxy::proxy;
		constexpr proxy<IsConstT>::Type& row(int i)
		{
			return proxy::vector(i);
		}
		constexpr const proxy<IsConstT>::Type& row(int i) const
		{
			return proxy::vector(i);
		}
	};
	template <bool IsConstT>
	struct row_proxy : proxy<IsConstT>
	{
		using proxy::proxy;
		constexpr proxy<IsConstT>::Type& column(int i)
		{
			return proxy::vector(i);
		}
		constexpr const proxy<IsConstT>::Type& column(int i) const
		{
			return proxy::vector(i);
		}
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

	matrix3_t operator*(T f) const;
	vector3_t<T> operator*(vector3_t<T> v) const;
	matrix3_t operator*(const matrix3_t& m) const;

	bool operator==(const matrix3_t& m) const
	{
		return m00 == m.m00 && m01 == m.m01 && m02 == m.m02 && m10 == m.m10 &&
		       m11 == m.m11 && m12 == m.m12 && m20 == m.m20 && m21 == m.m21 &&
		       m22 == m.m22;
	}
};
#pragma endregion

#pragma region declaration matrix4_t
template <typename T>
struct matrix4_t
{
	template <typename U>
	friend class matrix2_t;
	template <typename U>
	friend class matrix3_t;
	template <typename U>
	friend class matrix4_t;

private:
	T m00, m01, m02, m03,    //
	    m10, m11, m12, m13,  //
	    m20, m21, m22, m23,  //
	    m30, m31, m32, m33;

	matrix4_t(T m00,
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
	matrix4_t() = default;
	explicit matrix4_t(matrix2_t<T> m);
	explicit matrix4_t(const matrix3_t<T>& m);
	explicit matrix4_t(const perspective_t<T>& p);
	explicit matrix4_t(const transform3d_t<T>& t);
	explicit matrix4_t(const uniform_transform3d_t<T>& t);
	explicit matrix4_t(const unscaled_transform3d_t<T>& t);
	explicit matrix4_t(const std::array<T, 16>& a);
	matrix4_t(const matrix4_t&) = default;
	matrix4_t(matrix4_t&&) = default;

	matrix4_t& operator=(const matrix4_t&) = default;
	matrix4_t& operator=(matrix4_t&&) = default;

	template <typename U>
	explicit operator matrix4_t<U>() const
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

	template <typename U = T>
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

	static matrix4_t zero();
	static matrix4_t identity();
	static matrix4_t rows(vector4_t<T> row0,
	                      vector4_t<T> row1,
	                      vector4_t<T> row2,
	                      vector4_t<T> row3);
	static matrix4_t columns(vector4_t<T> column0,
	                         vector4_t<T> column1,
	                         vector4_t<T> column2,
	                         vector4_t<T> column3);
	static matrix4_t rotation(direction_t<T> axis, radian_t<T> angle);
	static matrix4_t rotation(direction_t<T> axis,
	                          normalized<complex_t<T>> angle);
	template <typename U, U NumeratorT, U DenominatorT>
	static matrix4_t rotation(
	    direction_t<T> axis,
	    static_pi_fraction_t<U, NumeratorT, DenominatorT> angle);
	static matrix4_t rotation(euler_angles_t<T> r);
	static matrix4_t rotation(quaternion_t<T> q);
	static matrix4_t translation(vector3_t<T> t);
	static matrix4_t translation(T x, T y, T z);
	static matrix4_t scale(vector3_t<T> t);
	static matrix4_t scale(T x, T y, T z);
	static matrix4_t scale(T s);
	static matrix4_t look_at(
	    vector3_t<T> from,
	    vector3_t<T> to,
	    direction_t<T> up = direction_t<T>::already_normalized(
	        vector3_t<T>(T(0), T(1), T(0))));
	static matrix4_t orthographic(T left,
	                              T right,
	                              T top,
	                              T bottom,
	                              T near,
	                              T far);
	static matrix4_t rotation_around_x(radian_t<T> angle);
	static matrix4_t rotation_around_y(radian_t<T> angle);
	static matrix4_t rotation_around_z(radian_t<T> angle);
	static matrix4_t rotation_around_x(normalized<complex_t<T>> angle);
	static matrix4_t rotation_around_y(normalized<complex_t<T>> angle);
	static matrix4_t rotation_around_z(normalized<complex_t<T>> angle);
	template <typename T, T NumeratorT, T DenominatorT>
	static matrix4_t rotation_around_x(
	    static_pi_fraction_t<T, NumeratorT, DenominatorT> angle);
	template <typename T, T NumeratorT, T DenominatorT>
	static matrix4_t rotation_around_y(
	    static_pi_fraction_t<T, NumeratorT, DenominatorT> angle);
	template <typename T, T NumeratorT, T DenominatorT>
	static matrix4_t rotation_around_z(
	    static_pi_fraction_t<T, NumeratorT, DenominatorT> angle);
	bool is_orthogonal() const;
	bool is_homogenous() const;
	matrix4_t& inverse();
	matrix4_t get_inversed() const;
	matrix4_t& transpose();
	matrix4_t get_transposed() const;
	T determinant() const;

	template <bool IsConstT>
	struct proxy
	{
	protected:
		using Type = std::conditional_t<IsConstT, const T, T>;

		std::reference_wrapper<Type> x;
		std::reference_wrapper<Type> y;
		std::reference_wrapper<Type> z;
		std::reference_wrapper<Type> w;

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
		constexpr proxy& operator=(vector4_t<T> v)
		{
			x = v.x;
			y = v.y;
			z = v.z;
			w = v.w;
			return *this;
		}

		constexpr vector2_t<T> xy() const { return vector2_t<T>{x, y}; }
		constexpr vector3_t<T> xyz() const { return vector3_t<T>{x, y, z}; }
		constexpr operator vector4_t<T>() const
		{
			return vector4_t<T>{x, y, z, w};
		}
	};
	template <bool IsConstT>
	struct column_proxy : proxy<IsConstT>
	{
		using proxy::proxy;
		constexpr proxy<IsConstT>::Type& row(int i)
		{
			return proxy::vector(i);
		}
		constexpr const proxy<IsConstT>::Type& row(int i) const
		{
			return proxy::vector(i);
		}
	};
	template <bool IsConstT>
	struct row_proxy : proxy<IsConstT>
	{
		using proxy::proxy;
		constexpr proxy<IsConstT>::Type& column(int i)
		{
			return proxy::vector(i);
		}
		constexpr const proxy<IsConstT>::Type& column(int i) const
		{
			return proxy::vector(i);
		}
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

	matrix4_t operator*(T f) const;
	matrix4_t operator/(T f) const;
	vector4_t<T> operator*(vector4_t<T> v) const;
	matrix4_t operator*(const matrix4_t& m) const;

	bool operator==(const matrix4_t& m) const
	{
		return m00 == m.m00 && m01 == m.m01 && m02 == m.m02 && m03 == m.m03 &&
		       m10 == m.m10 && m11 == m.m11 && m12 == m.m12 && m13 == m.m13 &&
		       m20 == m.m20 && m21 == m.m21 && m22 == m.m22 && m23 == m.m23 &&
		       m30 == m.m30 && m31 == m.m31 && m32 == m.m32 && m33 == m.m33;
	}
};
#pragma endregion

#pragma region declaration perspective_t
template <typename T>
struct perspective_t
{
private:
	radian_t<T> m_angle;
	T m_ratio;
	T m_invRatio;
	T m_near;
	T m_far;
	T m_invTanHalfFovy;

public:
	perspective_t(radian_t<T> angle, T ratio, T near, T far);

	void set(radian_t<T> angle, T ratio, T near, T far);

	void set_angle(radian_t<T> angle);
	radian_t<T> get_angle() const;

	void set_ratio(T ratio);
	T get_ratio() const;

	void set_near(T near_plane_t);
	T get_near() const;

	void set_far(T far_plane_t);
	T get_far() const;

	matrix4_t<T> to_matrix4() const;

	template <typename U>
	friend matrix4_t<U> operator*(const matrix4_t<U>& m,
	                              const perspective_t<U>& p);
	template <typename U>
	friend matrix4_t<U> operator*(const perspective_t<U>& p,
	                              const matrix4_t<U>& m);

	template <typename U>
	friend matrix4_t<U> operator*(const unscaled_transform3d_t<U>& t,
	                              const perspective_t<U>& p);
	template <typename U>
	friend matrix4_t<U> operator*(const perspective_t<U>& p,
	                              const unscaled_transform3d_t<U>& t);
};
#pragma endregion

#pragma region declaration euler_angles_t
template <typename T>
struct euler_angles_t
{
	radian_t<T> x, y, z;
};
#pragma endregion

#pragma region declaration quaternion_t
template <typename T>
struct quaternion_t
{
	T x, y, z, w;

	quaternion_t() = default;
	quaternion_t(T x_, T y_, T z_, T w_) : x{x_}, y{y_}, z{z_}, w{w_} {}
	quaternion_t(direction_t<T> axis, radian_t<T> angle);
	template <typename U, U NumeratorT, U DenominatorT>
	quaternion_t(direction_t<T> axis,
	             static_pi_fraction_t<U, NumeratorT, DenominatorT> angle);

	template <typename U>
	explicit operator quaternion_t<U>() const
	{
		return quaternion_t<U>{static_cast<U>(x), static_cast<U>(y),
		                       static_cast<U>(z), static_cast<U>(w)};
	}

	template <typename U = T>
	std::array<U, 4> to_array() const
	{
		return {static_cast<U>(x), static_cast<U>(y), static_cast<U>(z),
		        static_cast<U>(w)};
	}

	static quaternion_t zero();
	static quaternion_t identity();

	T norm() const;
	T norm_squared() const;

	quaternion_t& inverse();
	quaternion_t& conjugate();
	quaternion_t& normalize();
	quaternion_t& clamp(quaternion_t min, quaternion_t max);

	quaternion_t get_inversed() const;
	quaternion_t get_conjugated() const;
	quaternion_t get_normalized() const;
	quaternion_t get_clamped(quaternion_t min, quaternion_t max) const;

	quaternion_t operator+(quaternion_t q) const;
	quaternion_t operator-(quaternion_t q) const;
	quaternion_t operator*(
	    quaternion_t q) const;  // TODO: implement operator* for normalized<>
	                            // if result is normalized too
	vector3_t<T> operator*(vector3_t<T> v) const;
	quaternion_t operator*(T f) const;
	quaternion_t operator/(T f) const;

	quaternion_t& operator+=(quaternion_t q);
	quaternion_t& operator-=(quaternion_t q);
	quaternion_t& operator*=(quaternion_t q);
	quaternion_t& operator*=(T f);
	quaternion_t& operator/=(T f);

	quaternion_t operator-() const;

	bool operator==(quaternion_t v) const;
	bool operator!=(quaternion_t v) const;
};

template <typename T>
T dot(quaternion_t<T> lhs, quaternion_t<T> rhs);
template <typename T>
quaternion_t<T> cross(quaternion_t<T> lhs, quaternion_t<T> rhs);
template <typename T>
quaternion_t<T> lerp(quaternion_t<T> begin, quaternion_t<T> end, T percent);
template <typename T>
quaternion_t<T> nlerp(quaternion_t<T> begin, quaternion_t<T> end, T percent);
template <typename T>
quaternion_t<T> slerp(quaternion_t<T> begin, quaternion_t<T> end, T percent);
#pragma endregion

#pragma region declaration cartesian_direction2d_t
template <typename T>
struct cartesian_direction2d_t
{
	T x, y;

	static cartesian_direction2d_t from(
	    const normalized<vector2_t<T>>& direction);
	static cartesian_direction2d_t from(polar_direction_t<T> direction);
	static cartesian_direction2d_t from(radian_t<T> angle);

	cartesian_direction2d_t& rotate(radian_t<T> angle);
};
#pragma endregion

#pragma region declaration polar_direction_t
template <typename T>
struct polar_direction_t
{
	radian_t<T> angle;

	static polar_direction_t from(const normalized<vector2_t<T>>& direction);
	static polar_direction_t from(cartesian_direction2d_t<T> direction);
	static polar_direction_t from(radian_t<T> angle);

	polar_direction_t& rotate(radian_t<T> angle);

	polar_direction_t& wrap();
};
#pragma endregion

#pragma region declaration line_t
template <typename T>
struct line_t
{
	T coefficient;
	T offset;

	static line_t from(vector2_t<T> v);
};

template <typename T>
bool are_parallel(line_t<T> l1, line_t<T> l2);

template <typename T>
bool collides(line_t<T> l1, line_t<T> l2);
template <typename T>
bool collides(line_t<T> l1, line_t<T> l2, vector2_t<T>& intersection);
#pragma endregion

#pragma region declaration segment2d_t
template <typename T>
struct segment2d_t
{
	vector2_t<T> origin;
	vector2_t<T> end;
};

template <typename T>
bool collides(const segment2d_t<T>& s0,
              const segment2d_t<T>& s1,
              vector2_t<T>& outPoint,
              T e = T(epsilon));
#pragma endregion

#pragma region declaration plane_t
template <typename T>
struct plane_t
{
	vector3_t<T> origin;
	direction_t<T> normal;

	template <typename U>
	explicit operator plane_t<U>() const
	{
		return plane_t<U>{static_cast<vector3_t<U>>(origin),
		                  static_cast<vector3_t<U>>(*normal)};
	}

	T evaluate(vector3_t<T> point) const;
	vector3_t<T> project3d(vector3_t<T> point) const;
	vector2_t<T> project2d(vector3_t<T> point,
	                       direction_t<T> plane_t_tangent) const;
	vector3_t<T> unproject(vector2_t<T> point,
	                       direction_t<T> plane_t_tangent) const;

	direction_t<T> computeTangent(T e = T(epsilon)) const;
};
#pragma endregion

#pragma region declaration ray2d_t
template <typename T>
struct ray2d_t
{
	vector2_t<T> origin;
	normalized<vector2_t<T>> direction;
};
#pragma endregion

#pragma region declaration ray3d_t
template <typename T>
struct ray3d_t
{
	vector3_t<T> origin;
	direction_t<T> direction;
};
#pragma endregion

#pragma region declaration aa_rect_t
template <typename T>
struct aa_rect_t
{
	vector2_t<T> origin;
	vector2_t<T> dimention;

	bool contains(vector2_t<T> v) const;
};

template <typename T>
std::size_t collides(aa_rect_t<T> r1,
                     aa_rect_t<T> r2,
                     vector2_t<T>* intersection1,
                     vector2_t<T>* insersection2);
template <typename T>
std::size_t collides(aa_rect_t<T> r,
                     line_t<T> l,
                     vector2_t<T>* intersection1,
                     vector2_t<T>* insersection2);
template <typename T>
std::size_t collides(line_t<T> l,
                     aa_rect_t<T> r,
                     vector2_t<T>* intersection1,
                     vector2_t<T>* insersection2);
#pragma endregion

#pragma region declaration circle_t
template <typename T>
struct circle_t
{
	vector2_t<T> origin;
	T radius;
};

template <typename T>
std::size_t collides(circle_t<T> c1,
                     circle_t<T> c2,
                     vector2_t<T>* intersection1,
                     vector2_t<T>* insersection2);
template <typename T>
std::size_t collides(circle_t<T> c,
                     aa_rect_t<T> r,
                     vector2_t<T>* intersection1,
                     vector2_t<T>* insersection2);
template <typename T>
std::size_t collides(aa_rect_t<T> r,
                     circle_t<T> c,
                     vector2_t<T>* intersection1,
                     vector2_t<T>* insersection2);
template <typename T>
std::size_t collides(circle_t<T> c,
                     line_t<T> l,
                     vector2_t<T>* intersection1,
                     vector2_t<T>* insersection2);
template <typename T>
std::size_t collides(line_t<T> l,
                     circle_t<T> c,
                     vector2_t<T>* intersection1,
                     vector2_t<T>* insersection2);

template <typename T>
bool collides(ray3d_t<T> r, plane_t<T> p);
template <typename T>
bool collides(ray3d_t<T> r, plane_t<T> p, vector3_t<T>& intersection);
template <typename T>
bool collides(plane_t<T> p, ray3d_t<T> r);
template <typename T>
bool collides(plane_t<T> p, ray3d_t<T> r, vector3_t<T>& intersection);
template <typename T>
bool collides_bidirectional(const plane_t<T>& plane,
                            ray3d_t<T> r,
                            vector3_t<T>& collision_point);
#pragma endregion

#pragma region declaration aabb_t
template <typename T>
struct aabb_t
{
	vector3_t<T> min;
	vector3_t<T> max;

	bool contains(vector3_t<T> p) const;
	vector3_t<T> get_center() const;

	aabb_t operator+(aabb_t rhs) const;
};

template <typename T>
bool collides(aabb_t<T> b, ray3d_t<T> r);
template <typename T>
bool collides(aabb_t<T> b1, aabb_t<T> b2);
#pragma endregion

#pragma region declaration transform2d_t
template <typename T>
struct transform2d_t
{
	vector2_t<T> position;
	normalized<complex_t<T>> rotation;
	vector2_t<T> scale;

	transform2d_t operator*(const transform2d_t& t) const;
	vector2_t<T> operator*(vector2_t<T> v) const;
};
#pragma endregion

#pragma region declaration transform3d_t
template <typename T>
struct transform3d_t
{
	vector3_t<T> position;
	quaternion_t<T> rotation;
	vector3_t<T> scale;

	matrix4_t<T> to_matrix4() const
	{
		matrix4_t<T> Tr{matrix4_t<T>::translation(position)};
		matrix4_t<T> R{matrix4_t<T>::rotation(rotation)};
		matrix4_t<T> S{matrix4_t<T>::scale(scale)};
		return Tr * S * R;
	}

	transform3d_t operator*(const transform3d_t& t) const;
	vector3_t<T> operator*(vector3_t<T> v) const;
	quaternion_t<T> operator*(quaternion_t<T> q) const;
};
#pragma endregion

#pragma region declaration uniform_transform2d_t
template <typename T>
struct uniform_transform2d_t
{
	vector2_t<T> position;
	radian_t<T> rotation;
	T scale;

	uniform_transform2d_t operator*(uniform_transform2d_t t) const;
	vector2_t<T> operator*(vector2_t<T> v) const;
};
#pragma endregion

#pragma region declaration uniform_transform3d_t
template <typename T>
struct uniform_transform3d_t
{
	vector3_t<T> position;
	quaternion_t<T> rotation;
	T scale;

	matrix4_t<T> to_matrix4() const
	{
		matrix4_t<T> Tr{matrix4_t<T>::translation(position)};
		matrix4_t<T> R{matrix4_t<T>::rotation(rotation)};
		matrix4_t<T> S{matrix4_t<T>::scale(scale)};
		return Tr * S * R;
	}

	uniform_transform3d_t operator*(const uniform_transform3d_t& t) const;
	vector3_t<T> operator*(vector3_t<T> v) const;
	quaternion_t<T> operator*(quaternion_t<T> q) const;
};
#pragma endregion

#pragma region declaration unscaled_transform2d_t
template <typename T>
struct unscaled_transform2d_t
{
	vector2_t<T> position;
	radian_t<T> rotation;

	unscaled_transform2d_t operator*(unscaled_transform2d_t t) const;
	vector2_t<T> operator*(vector2_t<T> v) const;
};
#pragma endregion

#pragma region declaration unscaled_transform3d_t
template <typename T>
struct unscaled_transform3d_t
{
	vector3_t<T> position;
	quaternion_t<T> rotation;

	unscaled_transform3d_t& translate_absolute(vector3_t<T> t);
	unscaled_transform3d_t& translate_relative(vector3_t<T> t);
	unscaled_transform3d_t& rotate(quaternion_t<T> r);

	unscaled_transform3d_t& inverse()
	{
		rotation.inverse();
		position = rotation * -position;
		return *this;
	}
	unscaled_transform3d_t get_inversed() const
	{
		unscaled_transform3d_t res = *this;
		res.inverse();
		return res;
	}

	matrix4_t<T> to_matrix4() const
	{
		matrix4_t<T> Tr{matrix4_t<T>::translation(position)};
		matrix4_t<T> R{matrix4_t<T>::rotation(rotation)};
		return Tr * R;
	}

	unscaled_transform3d_t operator*(const unscaled_transform3d_t& t) const;
	vector3_t<T> operator*(vector3_t<T> v) const;
	quaternion_t<T> operator*(quaternion_t<T> q) const;
};
#pragma endregion

#pragma region definition misc
template <typename T>
unscaled_transform3d_t<T> inverse(unscaled_transform3d_t<T> tr);

template <typename T>
constexpr T lerp(T begin, T end, T percent)
{
	return (end - begin) * percent + begin;
}

template <typename T>
constexpr T clamp(T f, T min, T max)
{
	return std::min(std::max(f, min), max);
}

template <typename T>
bool nearly_equal(T f1, T f2, T e)
{
	return std::abs(f1 - f2) < e;
}

template <typename T>
constexpr int sign(T val)
{
	return (T(0) <= val) - (val < T(0));
}
template <typename T>
constexpr int signnum(T val)
{
	return (T(0) < val) - (val < T(0));
}
#pragma endregion

#pragma region definition complex_t
template <typename T>
complex_t<T>::complex_t(radian_t<T> angle) : r{cos(angle)}, i{sin(angle)}
{
}

template <typename T>
complex_t<T> complex_t<T>::operator+(complex_t<T> c) const
{
	return {r + c.r, i + c.i};
}
template <typename T>
complex_t<T> complex_t<T>::operator-(complex_t<T> c) const
{
	return {r - c.r, i - c.i};
}

template <typename T>
complex_t<T>& complex_t<T>::operator+=(complex_t<T> c)
{
	return *this = *this + c;
}
template <typename T>
complex_t<T>& complex_t<T>::operator-=(complex_t<T> c)
{
	return *this = *this - c;
}
template <typename T>
complex_t<T>& complex_t<T>::operator*=(complex_t<T> c)
{
	return *this = *this * c;
}

template <typename T>
complex_t<T> operator*(complex_t<T> c1, complex_t<T> c2)
{
	return {c1.r * c2.r - c1.i * c2.i, c1.r * c2.i + c1.i * c2.r};
}
template <typename T>
complex_t<T> operator*(complex_t<T> c, T f)
{
	return {c.r * f, c.i * f};
}
template <typename T>
complex_t<T> operator*(normalized<complex_t<T>> c1, complex_t<T> c2)
{
	return *c1 * c2;
}
template <typename T>
complex_t<T> operator*(complex_t<T> c1, normalized<complex_t<T>> c2)
{
	return c1 * *c2;
}
template <typename T>
normalized<complex_t<T>> operator*(normalized<complex_t<T>> c1,
                                   normalized<complex_t<T>> c2)
{
	return normalized<complex_t<T>>::already_normalized(complex_t<T>{
	    c1->r * c2->r - c1->i * c2->i, c1->r * c2->i + c1->i * c2->r});
}
template <typename T>
vector2_t<T> operator*(normalized<complex_t<T>> c, vector2_t<T> v)
{
	return {c->r * v.x - c->i * v.y, c->r * v.y + c->i * v.x};
}

template <typename T>
T norm(complex_t<T> c)
{
	return std::sqrt(norm_squared(c));
}
template <typename T>
T norm_squared(complex_t<T> c)
{
	return c.r * c.r + c.i * c.i;
}
template <typename T>
complex_t<T> normalize(complex_t<T> c)
{
	T n = norm(c);
	c.r /= n;
	c.i /= n;
	return c;
}
template <typename T>
complex_t<T> conjugate(complex_t<T> c)
{
	c.i = -c.i;
	return c;
}
#pragma endregion

#pragma region definition radian_t
template <typename T>
constexpr radian_t<T>::radian_t(T f) : angle(f)
{
}
template <typename T>
constexpr radian_t<T>::radian_t(degree_t<T> d)
    : angle(static_cast<T>(d) * T(deg_to_rad))
{
}

template <typename T>
radian_t<T>::operator T() const
{
	return angle;
}

template <typename T>
radian_t<T>& radian_t<T>::operator+=(T f)
{
	angle += f;
	return *this;
}
template <typename T>
radian_t<T>& radian_t<T>::operator+=(radian_t<T> r)
{
	angle += r.angle;
	return *this;
}
template <typename T>
radian_t<T>& radian_t<T>::operator-=(T f)
{
	angle -= f;
	return *this;
}
template <typename T>
radian_t<T>& radian_t<T>::operator-=(radian_t<T> r)
{
	angle -= r.angle;
	return *this;
}
template <typename T>
radian_t<T>& radian_t<T>::operator*=(T f)
{
	angle *= f;
	return *this;
}
template <typename T>
radian_t<T>& radian_t<T>::operator*=(radian_t<T> r)
{
	angle *= r.angle;
	return *this;
}
template <typename T>
radian_t<T>& radian_t<T>::operator/=(T f)
{
	angle /= f;
	return *this;
}
template <typename T>
radian_t<T>& radian_t<T>::operator/=(radian_t<T> r)
{
	angle /= r.angle;
	return *this;
}

template <typename T>
radian_t<T> radian_t<T>::operator-() const
{
	return radian_t<T>{-angle};
}

template <typename T>
radian_t<T>& radian_t<T>::operator=(degree_t<T> d)
{
	angle = d.angle * T(deg_to_rad);
	return *this;
}

template <typename T>
radian_t<T> operator+(radian_t<T> r1, radian_t<T> r2)
{
	return radian_t<T>{static_cast<T>(r1) + static_cast<T>(r2)};
}
template <typename T>
radian_t<T> operator-(radian_t<T> r1, radian_t<T> r2)
{
	return radian_t<T>{static_cast<T>(r1) - static_cast<T>(r2)};
}
template <typename T>
radian_t<T> operator*(radian_t<T> r, T f)
{
	return radian_t<T>{static_cast<T>(r) * f};
}
template <typename T>
radian_t<T> operator*(T f, radian_t<T> r)
{
	return radian_t<T>{f * static_cast<T>(r)};
}
template <typename T>
radian_t<T> operator/(radian_t<T> r, T f)
{
	return radian_t<T>{static_cast<T>(r) / f};
}
template <typename T>
bool operator<(radian_t<T> lhs, radian_t<T> rhs)
{
	return T(lhs) < T(rhs);
}
template <typename T>
bool operator>(radian_t<T> lhs, radian_t<T> rhs)
{
	return T(lhs) > T(rhs);
}
template <typename T>
bool operator==(radian_t<T> lhs, radian_t<T> rhs)
{
	return T(lhs) == T(rhs);
}
template <typename T>
bool operator!=(radian_t<T> lhs, radian_t<T> rhs)
{
	return T(lhs) != T(rhs);
}
template <typename T>
bool operator>=(radian_t<T> lhs, radian_t<T> rhs)
{
	return T(lhs) >= T(rhs);
}
template <typename T>
bool operator<=(radian_t<T> lhs, radian_t<T> rhs)
{
	return T(lhs) <= T(rhs);
}

template <typename T>
T sin(radian_t<T> r)
{
	return std::sin(static_cast<T>(r));
}
template <typename T>
T cos(radian_t<T> r)
{
	return std::cos(static_cast<T>(r));
}
template <typename T>
T tan(radian_t<T> r)
{
	return std::tan(static_cast<T>(r));
}
template <typename T>
T asin(radian_t<T> r)
{
	return std::asin(static_cast<T>(r));
}
template <typename T>
T acos(radian_t<T> r)
{
	return std::acos(static_cast<T>(r));
}
template <typename T>
T atan(radian_t<T> r)
{
	return std::atan(static_cast<T>(r));
}
template <typename T>
T sinh(radian_t<T> r)
{
	return std::sinh(static_cast<T>(r));
}
template <typename T>
T cosh(radian_t<T> r)
{
	return std::cosh(static_cast<T>(r));
}
template <typename T>
T tanh(radian_t<T> r)
{
	return std::tanh(static_cast<T>(r));
}
template <typename T>
T asinh(radian_t<T> r)
{
	return std::asinh(static_cast<T>(r));
}
template <typename T>
T acosh(radian_t<T> r)
{
	return std::acosh(static_cast<T>(r));
}
template <typename T>
T atanh(radian_t<T> r)
{
	return std::atanh(static_cast<T>(r));
}
#pragma endregion

#pragma region definition degree_t
template <typename T>
constexpr degree_t<T>::degree_t(T f) : angle(f)
{
}
template <typename T>
constexpr degree_t<T>::degree_t(radian_t<T> r)
    : angle(static_cast<T>(r) * T(rad_to_deg))
{
}

template <typename T>
degree_t<T>::operator T() const
{
	return angle;
}

template <typename T>
degree_t<T>& degree_t<T>::operator+=(T f)
{
	angle += f;
	return *this;
}
template <typename T>
degree_t<T>& degree_t<T>::operator+=(degree_t d)
{
	angle += d.angle;
	return *this;
}
template <typename T>
degree_t<T>& degree_t<T>::operator-=(T f)
{
	angle -= f;
	return *this;
}
template <typename T>
degree_t<T>& degree_t<T>::operator-=(degree_t d)
{
	angle -= d.angle;
	return *this;
}
template <typename T>
degree_t<T>& degree_t<T>::operator*=(T f)
{
	angle *= f;
	return *this;
}
template <typename T>
degree_t<T>& degree_t<T>::operator*=(degree_t d)
{
	angle *= d.angle;
	return *this;
}
template <typename T>
degree_t<T>& degree_t<T>::operator/=(T f)
{
	angle /= f;
	return *this;
}
template <typename T>
degree_t<T>& degree_t<T>::operator/=(degree_t d)
{
	angle /= d.angle;
	return *this;
}

template <typename T>
degree_t<T> degree_t<T>::operator-() const
{
	return degree_t(-angle);
}

template <typename T>
degree_t<T> operator+(degree_t<T> d1, degree_t<T> d2)
{
	return degree_t<T>{static_cast<T>(d1) + static_cast<T>(d2)};
}
template <typename T>
degree_t<T> operator-(degree_t<T> d1, degree_t<T> d2)
{
	return degree_t<T>{static_cast<T>(d1) - static_cast<T>(d2)};
}
template <typename T>
degree_t<T> operator*(degree_t<T> d, T f)
{
	return degree_t<T>{static_cast<T>(d) * f};
}
template <typename T>
degree_t<T> operator*(T f, degree_t<T> d)
{
	return degree_t<T>{f * static_cast<T>(d)};
}
template <typename T>
degree_t<T> operator/(degree_t<T> d, T f)
{
	return degree_t<T>{static_cast<T>(d) / f};
}
template <typename T>
bool operator<(degree_t<T> lhs, degree_t<T> rhs)
{
	return T(lhs) < T(rhs);
}
template <typename T>
bool operator>(degree_t<T> lhs, degree_t<T> rhs)
{
	return T(lhs) > T(rhs);
}
template <typename T>
bool operator==(degree_t<T> lhs, degree_t<T> rhs)
{
	return T(lhs) == T(rhs);
}
template <typename T>
bool operator!=(degree_t<T> lhs, degree_t<T> rhs)
{
	return T(lhs) != T(rhs);
}
template <typename T>
bool operator>=(degree_t<T> lhs, degree_t<T> rhs)
{
	return T(lhs) >= T(rhs);
}
template <typename T>
bool operator<=(degree_t<T> lhs, degree_t<T> rhs)
{
	return T(lhs) <= T(rhs);
}

template <typename T>
T sin(degree_t<T> d)
{
	return sin(radian_t<T>(d));
}
template <typename T>
T cos(degree_t<T> d)
{
	return cos(radian_t<T>(d));
}
template <typename T>
T tan(degree_t<T> d)
{
	return tan(radian_t<T>(d));
}
template <typename T>
T asin(degree_t<T> d)
{
	return asin(radian_t<T>(d));
}
template <typename T>
T acos(degree_t<T> d)
{
	return acos(radian_t<T>(d));
}
template <typename T>
T atan(degree_t<T> d)
{
	return atan(radian_t<T>(d));
}
template <typename T>
T sinh(degree_t<T> d)
{
	return sinh(radian_t<T>(d));
}
template <typename T>
T cosh(degree_t<T> d)
{
	return cosh(radian_t<T>(d));
}
template <typename T>
T tanh(degree_t<T> d)
{
	return tanh(radian_t<T>(d));
}
template <typename T>
T asinh(degree_t<T> d)
{
	return asinh(radian_t<T>(d));
}
template <typename T>
T acosh(degree_t<T> d)
{
	return acosh(radian_t<T>(d));
}
template <typename T>
T atanh(degree_t<T> d)
{
	return atanh(radian_t<T>(d));
}
#pragma endregion

#pragma region definition vector2_t
template <typename T>
vector2_t<T> vector2_t<T>::operator+(vector2_t<T> v) const
{
	return {x + v.x, y + v.y};
}
template <typename T>
vector2_t<T> vector2_t<T>::operator-(vector2_t<T> v) const
{
	return {x - v.x, y - v.y};
}
template <typename T>
vector2_t<T> vector2_t<T>::operator*(T f) const
{
	return {x * f, y * f};
}
template <typename T>
vector2_t<T> vector2_t<T>::operator/(T f) const
{
	return {x / f, y / f};
}

template <typename T>
vector2_t<T>& vector2_t<T>::operator+=(vector2_t<T> v)
{
	*this = *this + v;
	return *this;
}
template <typename T>
vector2_t<T>& vector2_t<T>::operator-=(vector2_t<T> v)
{
	*this = *this - v;
	return *this;
}
template <typename T>
vector2_t<T>& vector2_t<T>::operator*=(T f)
{
	*this = *this * f;
	return *this;
}
template <typename T>
vector2_t<T>& vector2_t<T>::operator/=(T f)
{
	*this = *this / f;
	return *this;
}

template <typename T>
vector2_t<T> vector2_t<T>::operator-() const
{
	return {-x, -y};
}

template <typename T>
bool vector2_t<T>::operator==(vector2_t<T> v) const
{
	return x == v.x && y == v.y;
}
template <typename T>
bool vector2_t<T>::operator!=(vector2_t<T> v) const
{
	return !operator==(v);
}

template <typename T>
vector2_t<T> operator*(T f, vector2_t<T> v)
{
	return v * f;
}

template <typename T>
T norm(vector2_t<T> v)
{
	return std::sqrt(norm_squared(v));
}
template <typename T>
T norm_squared(vector2_t<T> v)
{
	return v.x * v.x + v.y * v.y;
}
template <typename T>
vector2_t<T>& normalize(vector2_t<T>& v)
{
	T n = norm(v);
	v.x /= n;
	v.y /= n;
	return v;
}
template <typename T>
vector2_t<T> clamp(vector2_t<T> v, vector2_t<T> min, vector2_t<T> max)
{
	v.x = cdm::clamp(v.x, min.x, max.x);
	v.y = cdm::clamp(v.y, min.y, max.y);
	return v;
}
template <typename T>
T dot(vector2_t<T> lhs, vector2_t<T> rhs)
{
	return lhs.x * rhs.x + lhs.y * rhs.y;
}
template <typename T>
T cross(vector2_t<T> lhs, vector2_t<T> rhs)
{
	return lhs.x * rhs.y - lhs.y * rhs.x;
}
template <typename T>
vector2_t<T> lerp(vector2_t<T> begin, vector2_t<T> end, T percent)
{
	return (end - begin) * percent + begin;
}
template <typename T>
vector2_t<T> nlerp(vector2_t<T> begin, vector2_t<T> end, T percent)
{
	return lerp(begin, end, percent).get_normalized();
}
template <typename T>
vector2_t<T> slerp(vector2_t<T> begin, vector2_t<T> end, T percent)
{
	const radian_t angle = angle_between(begin, end) * percent;
	const T s = sin(angle);
	const T c = cos(angle);

	normalized<vector2_t<T>> res{
	    {c * begin.x - s * begin.y, s * begin.x + c * begin.y}};

	T f = cdm::lerp(norm(begin), norm(end), percent);
	return res * f;
}
template <typename T>
T distance_between(vector2_t<T> v1, vector2_t<T> v2)
{
	return norm(v1 - v2);
}
template <typename T>
T distance_squared_between(vector2_t<T> v1, vector2_t<T> v2)
{
	return norm_squared(v1 - v2);
}
template <typename T>
vector2_t<T> from_to(vector2_t<T> from, vector2_t<T> to)
{
	return {to.x - from.x, to.y - from.y};
}
template <typename T>
radian_t<T> angle_between(vector2_t<T> v1, vector2_t<T> v2)
{
	return radian_t{atan2f(v2.y, v2.x) - atan2f(v1.y, v1.x)};
}
template <typename T>
bool nearly_equal(vector2_t<T> v1, vector2_t<T> v2, T e)
{
	return cdm::nearly_equal(v1.x, v2.x, e) &&
	       cdm::nearly_equal(v1.y, v2.y, e);
}
template <typename T>
vector2_t<T> element_wise_min(vector2_t<T> v0, vector2_t<T> v1)
{
	return {v0.x < v1.x ? v0.x : v1.x, v0.y < v1.y ? v0.y : v1.y};
};
template <typename T>
vector2_t<T> element_wise_max(vector2_t<T> v0, vector2_t<T> v1)
{
	return {v0.x > v1.x ? v0.x : v1.x, v0.y > v1.y ? v0.y : v1.y};
};
#pragma endregion

#pragma region definition vector3_t
template <typename T>
vector2_t<T> vector3_t<T>::xy() const
{
	return {x, y};
}

template <typename T>
vector3_t<T> vector3_t<T>::operator+(vector3_t<T> v) const
{
	return {x + v.x, y + v.y, z + v.z};
}
template <typename T>
vector3_t<T> vector3_t<T>::operator-(vector3_t<T> v) const
{
	return {x - v.x, y - v.y, z - v.z};
}
template <typename T>
vector3_t<T> vector3_t<T>::operator*(T f) const
{
	return {x * f, y * f, z * f};
}
template <typename T>
vector3_t<T> vector3_t<T>::operator/(T f) const
{
	return {x / f, y / f, z / f};
}

template <typename T>
vector3_t<T>& vector3_t<T>::operator+=(vector3_t<T> v)
{
	*this = *this + v;
	return *this;
}
template <typename T>
vector3_t<T>& vector3_t<T>::operator-=(vector3_t<T> v)
{
	*this = *this - v;
	return *this;
}
template <typename T>
vector3_t<T>& vector3_t<T>::operator*=(T f)
{
	*this = *this * f;
	return *this;
}
template <typename T>
vector3_t<T>& vector3_t<T>::operator/=(T f)
{
	*this = *this / f;
	return *this;
}

template <typename T>
vector3_t<T> vector3_t<T>::operator-() const
{
	return {-x, -y, -z};
}

template <typename T>
bool vector3_t<T>::operator==(vector3_t<T> v) const
{
	return x == v.x && y == v.y && z == v.z;
}
template <typename T>
bool vector3_t<T>::operator!=(vector3_t<T> v) const
{
	return !operator==(v);
}

template <typename T>
vector3_t<T> operator*(T f, vector3_t<T> v)
{
	return v * f;
}

template <typename T>
T norm(vector3_t<T> v)
{
	return std::sqrt(norm_squared(v));
}
template <typename T>
T norm_squared(vector3_t<T> v)
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}
template <typename T>
vector3_t<T> normalize(vector3_t<T> v)
{
	T n = norm(v);
	v.x /= n;
	v.y /= n;
	v.z /= n;
	return v;
}
template <typename T>
vector3_t<T> clamp(vector3_t<T> v, vector3_t<T> min, vector3_t<T> max)
{
	v.x = cdm::clamp(v.x, min.x, max.x);
	v.y = cdm::clamp(v.y, min.y, max.y);
	v.z = cdm::clamp(v.z, min.z, max.z);
	return v;
}
template <typename T>
T dot(vector3_t<T> lhs, vector3_t<T> rhs)
{
	return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}
template <typename T>
vector3_t<T> cross(vector3_t<T> lhs, vector3_t<T> rhs)
{
	return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z,
	        lhs.x * rhs.y - lhs.y * rhs.x};
}
template <typename T>
vector3_t<T> lerp(vector3_t<T> begin, vector3_t<T> end, T percent)
{
	return (end - begin) * percent + begin;
}
template <typename T>
vector3_t<T> nlerp(vector3_t<T> begin, vector3_t<T> end, T percent)
{
	return normalize(lerp(begin, end, percent));
}
template <typename T>
T distance_between(vector3_t<T> v1, vector3_t<T> v2)
{
	return norm(v1 - v2);
}
template <typename T>
T distance_squared_between(vector3_t<T> v1, vector3_t<T> v2)
{
	return norm_squared(v1 - v2);
}
template <typename T>
vector3_t<T> from_to(vector3_t<T> from, vector3_t<T> to)
{
	return {to.x - from.x, to.y - from.y, to.z - from.z};
}
template <typename T>
radian_t<T> angle_between(vector3_t<T> v1, vector3_t<T> v2)
{
	T divisor = std::sqrt(norm_squared(v1) * norm_squared(v2));
	T alpha = dot(v1, v2) / divisor;
	return radian_t<T>(std::acosf(cdm::clamp(alpha, T(-1), T(1))));
}
template <typename T>
bool nearly_equal(vector3_t<T> v1, vector3_t<T> v2, T e)
{
	return nearly_equal(v1.x, v2.x, e) && nearly_equal(v1.y, v2.y, e) &&
	       nearly_equal(v1.z, v2.z, e);
}
template <typename T>
vector3_t<T> element_wise_min(vector3_t<T> v0, vector3_t<T> v1)
{
	return {v0.x < v1.x ? v0.x : v1.x, v0.y < v1.y ? v0.y : v1.y,
	        v0.z < v1.z ? v0.z : v1.z};
};
template <typename T>
vector3_t<T> element_wise_max(vector3_t<T> v0, vector3_t<T> v1)
{
	return {v0.x > v1.x ? v0.x : v1.x, v0.y > v1.y ? v0.y : v1.y,
	        v0.z > v1.z ? v0.z : v1.z};
};
template <typename T>
radian_t<T> angle_around_axis(vector3_t<T> v0,
                              vector3_t<T> v1,
                              direction_t<T> axis)
{
	vector3_t<T> c = cross(v0, v1);
	T angle = atan2f(norm(c), dot(v0, v1));
	return radian_t<T>{dot(c, *axis) < T(0) ? -angle : angle};
}
#pragma endregion

#pragma region definition vector4_t
template <typename T>
vector2_t<T> vector4_t<T>::xy() const
{
	return {x, y};
}
template <typename T>
vector3_t<T> vector4_t<T>::xyz() const
{
	return {x, y, z};
}

template <typename T>
vector4_t<T> vector4_t<T>::operator+(vector4_t<T> v) const
{
	return {x + v.x, y + v.y, z + v.z, w + v.w};
}
template <typename T>
vector4_t<T> vector4_t<T>::operator-(vector4_t<T> v) const
{
	return {x - v.x, y - v.y, z - v.z, w - v.w};
}
template <typename T>
vector4_t<T> vector4_t<T>::operator*(T f) const
{
	return {x * f, y * f, z * f, w * f};
}
template <typename T>
vector4_t<T> vector4_t<T>::operator/(T f) const
{
	return {x / f, y / f, z / f, w / f};
}

template <typename T>
vector4_t<T>& vector4_t<T>::operator+=(vector4_t<T> v)
{
	*this = *this + v;
	return *this;
}
template <typename T>
vector4_t<T>& vector4_t<T>::operator-=(vector4_t<T> v)
{
	*this = *this - v;
	return *this;
}
template <typename T>
vector4_t<T>& vector4_t<T>::operator*=(T f)
{
	*this = *this * f;
	return *this;
}
template <typename T>
vector4_t<T>& vector4_t<T>::operator/=(T f)
{
	*this = *this / f;
	return *this;
}

template <typename T>
vector4_t<T> vector4_t<T>::operator-() const
{
	return {-x, -y, -z, -w};
}

template <typename T>
bool vector4_t<T>::operator==(vector4_t<T> v) const
{
	return x == v.x && y == v.y && z == v.z && w == v.w;
}
template <typename T>
bool vector4_t<T>::operator!=(vector4_t<T> v) const
{
	return !operator==(v);
}

template <typename T>
T norm(vector4_t<T> v)
{
	return std::sqrt(norm_squared(v));
}
template <typename T>
T norm_squared(vector4_t<T> v)
{
	return v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w;
}
template <typename T>
vector4_t<T> normalize(vector4_t<T> v)
{
	T n = norm(v);
	v.x /= n;
	v.y /= n;
	v.z /= n;
	v.w /= n;
	return v;
}
template <typename T>
vector4_t<T> clamp(vector4_t<T> v, vector4_t<T> min, vector4_t<T> max)
{
	v.x = cdm::clamp(v.x, min.x, max.x);
	v.y = cdm::clamp(v.y, min.y, max.y);
	v.z = cdm::clamp(v.z, min.z, max.z);
	v.w = cdm::clamp(v.w, min.w, max.w);
	return v;
}
template <typename T>
vector4_t<T> negate(vector4_t<T> v)
{
	v.x = -v.x;
	v.y = -v.y;
	v.z = -v.z;
	v.w = -v.w;
	return v;
}
template <typename T>
T dot(vector4_t<T> lhs, vector4_t<T> rhs)
{
	return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w;
}
template <typename T>
vector4_t<T> lerp(vector4_t<T> begin, vector4_t<T> end, T percent)
{
	return (end - begin) * percent + begin;
}
template <typename T>
vector4_t<T> nlerp(vector4_t<T> begin, vector4_t<T> end, T percent)
{
	return normalized(lerp(begin, end, percent));
}
template <typename T>
T distance_between(vector4_t<T> v1, vector4_t<T> v2)
{
	return norm(v1 - v2);
}
template <typename T>
T distance_squared_between(vector4_t<T> v1, vector4_t<T> v2)
{
	return norm_squared(v1 - v2);
}
template <typename T>
vector4_t<T> from_to(vector4_t<T> from, vector4_t<T> to)
{
	return {to.x - from.x, to.y - from.y, to.z - from.z, to.w - from.w};
}
template <typename T>
bool nearly_equal(vector4_t<T> v1, vector4_t<T> v2, T e)
{
	return nearly_equal(v1.x, v2.x, e) && nearly_equal(v1.y, v2.y, e) &&
	       nearly_equal(v1.z, v2.z, e) && nearly_equal(v1.w, v2.w, e);
}
template <typename T>
vector4_t<T> element_wise_min(vector4_t<T> v0, vector4_t<T> v1)
{
	return {v0.x < v1.x ? v0.x : v1.x, v0.y < v1.y ? v0.y : v1.y,
	        v0.z < v1.z ? v0.z : v1.z, v0.w < v1.w ? v0.w : v1.w};
};
template <typename T>
vector4_t<T> element_wise_max(vector4_t<T> v0, vector4_t<T> v1)
{
	return {v0.x > v1.x ? v0.x : v1.x, v0.y > v1.y ? v0.y : v1.y,
	        v0.z > v1.z ? v0.z : v1.z, v0.w > v1.w ? v0.w : v1.w};
};
#pragma endregion

#pragma region definition normalized
template <typename T>
normalized<T>::normalized(const T& t) : vector(normalize(t))
{
}

template <typename T>
normalized<T> normalized<T>::already_normalized(const T& t)
{
	normalized res;
	res.vector = t;
	return res;
}
template <typename T>
normalized<T> normalized<T>::already_normalized(T&& t) noexcept
{
	normalized res;
	res.vector = std::move(t);
	return res;
}

template <typename T>
const T& normalized<T>::operator*() const
{
	return vector;
}
template <typename T>
const T* normalized<T>::operator->() const
{
	return &vector;
}
template <typename T>
normalized<T>::operator const T&() const
{
	return vector;
}

template <typename T>
normalized<T> normalized<T>::operator+(const T& v) const
{
	return normalize(vector + v);
}
template <typename T>
normalized<T> normalized<T>::operator-(const T& v) const
{
	return normalize(vector - v);
}
template <typename T>
T normalized<T>::operator*(T f) const
{
	return vector * f;
}
template <typename T>
T normalized<T>::operator/(T f) const
{
	return vector / f;
}

template <typename T>
normalized<T>& normalized<T>::operator+=(const T& v)
{
	vector = normalize(vector + v);
	return *this;
}
template <typename T>
normalized<T>& normalized<T>::operator-=(const T& v)
{
	vector = normalize(vector - v);
	return *this;
}

template <typename T>
normalized<T> normalized<T>::operator-() const
{
	return -vector;
}

template <typename T>
bool normalized<T>::operator==(const T& v) const
{
	return vector == v;
}
template <typename T>
bool normalized<T>::operator!=(const T& v) const
{
	return vector != v;
}

template <typename T>
normalized<T>& normalized<T>::operator=(const T& t)
{
	vector = normalize(t);
	return *this;
}
template <typename T>
normalized<T>& normalized<T>::operator=(T&& t) noexcept
{
	vector = normalize(t);
	return *this;
}
#pragma endregion

#pragma region definition matrix2_t
template <typename T>
matrix2_t<T>::matrix2_t(T e00, T e10, T e01, T e11)
    : m00{e00}, m10{e10}, m01{e01}, m11{e11}
{
}
template <typename T>
matrix2_t<T>::matrix2_t(const std::array<T, 4>& a)
    : matrix2_t(a[0], a[1], a[2], a[3])
{
}
template <typename T>
matrix2_t<T> matrix2_t<T>::zero()
{
	return {T(0), T(0), T(0), T(0)};
}
template <typename T>
matrix2_t<T> matrix2_t<T>::identity()
{
	return {T(1), T(0), T(0), T(1)};
}
template <typename T>
matrix2_t<T> matrix2_t<T>::rows(vector2_t<T> row0, vector2_t<T> row1)
{
	return {row0.x, row0.y, row1.x, row1.y};
}
template <typename T>
matrix2_t<T> matrix2_t<T>::columns(vector2_t<T> column0, vector2_t<T> column1)
{
	return {column0.x, column1.x, column0.y, column1.y};
}
template <typename T>
matrix2_t<T> matrix2_t<T>::rotation(radian_t<T> angle)
{
	T c = cos(angle);
	T s = sin(angle);
	return {c, -s, s, c};
}
template <typename T>
matrix2_t<T> matrix2_t<T>::rotation(normalized<complex_t<T>> angle)
{
	return {angle->r, -angle->i, angle->i, angle->r};
}

template <typename T>
vector2_t<T> matrix2_t<T>::operator*(vector2_t<T> v) const
{
	return {m00 * v.x + m01 * v.y, m10 * v.x + m11 * v.y};
}
template <typename T>
matrix2_t<T> matrix2_t<T>::operator*(matrix2_t<T> m) const
{
	return {m00 * m.m00 + m01 * m.m10, m00 * m.m01 + m01 * m.m11,

	        m10 * m.m00 + m11 * m.m10, m10 * m.m01 + m11 * m.m11};
}

template <typename T>
matrix2_t<T> transpose(matrix2_t<T> m)
{
	std::swap(m.m01, m.m10);
	return m;
}
#pragma endregion

#pragma region definition matrix3_t
template <typename T>
matrix3_t<T>::matrix3_t(T e00,
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
template <typename T>
matrix3_t<T>::matrix3_t(matrix2_t<T> m)
    : matrix3_t(m.m00, m.m10, T(0), m.m10, m.m11, T(0), T(0), T(0), T(1))
{
}
template <typename T>
matrix3_t<T>::matrix3_t(const std::array<T, 9>& a)
    : matrix3_t(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8])
{
}
template <typename T>
matrix3_t<T> matrix3_t<T>::zero()
{
	return {T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0)};
}
template <typename T>
matrix3_t<T> matrix3_t<T>::identity()
{
	return {T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1)};
}
template <typename T>
matrix3_t<T> matrix3_t<T>::rows(vector3_t<T> row0,
                                vector3_t<T> row1,
                                vector3_t<T> row2)
{
	return {row0.x, row0.y, row0.z, row1.x, row1.y,
	        row1.z, row2.x, row2.y, row2.z};
}
template <typename T>
matrix3_t<T> matrix3_t<T>::columns(vector3_t<T> column0,
                                   vector3_t<T> column1,
                                   vector3_t<T> column2)
{
	return {column0.x, column1.x, column2.x, column0.y, column1.y,
	        column2.y, column0.z, column1.z, column2.z};
}
template <typename T>
matrix3_t<T> matrix3_t<T>::rotation(direction_t<T> axis,
                                    T sinAngle,
                                    T cosAngle)
{
	const T l{axis->x};
	const T m{axis->y};
	const T n{axis->z};
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
template <typename T>
matrix3_t<T> matrix3_t<T>::rotation(direction_t<T> axis, radian_t<T> angle)
{
	return matrix3_t<T>::rotation(axis, sin(angle), cos(angle));
}
template <typename T>
matrix3_t<T> matrix3_t<T>::rotation(direction_t<T> axis,
                                    normalized<complex_t<T>> angle)
{
	return matrix3_t<T>::rotation(axis, sin(angle), cos(angle));
}
template <typename T>
template <typename U, U NumeratorT, U DenominatorT>
matrix3_t<T> matrix3_t<T>::rotation(
    direction_t<T> axis,
    static_pi_fraction_t<U, NumeratorT, DenominatorT> angle)
{
	return matrix3_t<T>::rotation(axis, sin<T>(angle), cos<T>(angle));
}
template <typename T>
matrix3_t<T> matrix3_t<T>::rotation(euler_angles_t<T> r)
{
	matrix3_t<T> RX = matrix3_t<T>::identity();
	RX.m11 = cos(r.x);
	RX.m12 = -sin(r.x);
	RX.m21 = sin(r.x);
	RX.m22 = cos(r.x);

	matrix3_t<T> RY = matrix3_t<T>::identity();
	RY.m00 = cos(r.y);
	RY.m02 = sin(r.y);
	RY.m20 = -sin(r.y);
	RY.m22 = cos(r.y);

	matrix3_t<T> RZ = matrix3_t<T>::identity();
	RZ.m00 = cos(r.z);
	RZ.m01 = -sin(r.z);
	RZ.m10 = sin(r.z);
	RZ.m11 = cos(r.z);

	return RY * RX * RZ;
}
template <typename T>
matrix3_t<T> matrix3_t<T>::rotation(quaternion_t<T> q)
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
template <typename T>
matrix3_t<T> matrix3_t<T>::rotation_around_x(T sinAngle, T cosAngle)
{
	return {
	    T(1), T(0),     T(0),       //
	    T(0), cosAngle, -sinAngle,  //
	    T(0), sinAngle, cosAngle    //
	};
}
template <typename T>
matrix3_t<T> matrix3_t<T>::rotation_around_y(T sinAngle, T cosAngle)
{
	return {
	    cosAngle,  T(0), sinAngle,  //
	    T(0),      T(1), T(0),      //
	    -sinAngle, T(0), cosAngle   //
	};
}
template <typename T>
matrix3_t<T> matrix3_t<T>::rotation_around_z(T sinAngle, T cosAngle)
{
	return {
	    cosAngle, -sinAngle, T(0),  //
	    sinAngle, cosAngle,  T(0),  //
	    T(0),     T(0),      T(1)   //
	};
}
template <typename T>
matrix3_t<T> matrix3_t<T>::rotation_around_x(radian_t<T> angle)
{
	return rotation_around_x(sin(angle), cos(angle));
}
template <typename T>
matrix3_t<T> matrix3_t<T>::rotation_around_y(radian_t<T> angle)
{
	return rotation_around_y(sin(angle), cos(angle));
}
template <typename T>
matrix3_t<T> matrix3_t<T>::rotation_around_z(radian_t<T> angle)
{
	return rotation_around_z(sin(angle), cos(angle));
}
template <typename T>
matrix3_t<T> matrix3_t<T>::rotation_around_x(normalized<complex_t<T>> angle)
{
	return rotation_around_x(angle->i, angle->r);
}
template <typename T>
matrix3_t<T> matrix3_t<T>::rotation_around_y(normalized<complex_t<T>> angle)
{
	return rotation_around_y(angle->i, angle->r);
}
template <typename T>
matrix3_t<T> matrix3_t<T>::rotation_around_z(normalized<complex_t<T>> angle)
{
	return rotation_around_z(angle->i, angle->r);
}
template <typename T>
template <typename U, U NumeratorT, U DenominatorT>
matrix3_t<T> matrix3_t<T>::rotation_around_x(
    static_pi_fraction_t<U, NumeratorT, DenominatorT> angle)
{
	return rotation_around_x(sin<T>(angle), cos<T>(angle));
}
template <typename T>
template <typename U, U NumeratorT, U DenominatorT>
matrix3_t<T> matrix3_t<T>::rotation_around_y(
    static_pi_fraction_t<U, NumeratorT, DenominatorT> angle)
{
	return rotation_around_y(sin<T>(angle), cos<T>(angle));
}
template <typename T>
template <typename U, U NumeratorT, U DenominatorT>
matrix3_t<T> matrix3_t<T>::rotation_around_z(
    static_pi_fraction_t<U, NumeratorT, DenominatorT> angle)
{
	return rotation_around_z(sin<T>(angle), cos<T>(angle));
}
template <typename T>
matrix3_t<T> matrix3_t<T>::scale(vector3_t<T> s)
{
	return scale(s.x, s.y, s.z);
}
template <typename T>
matrix3_t<T> matrix3_t<T>::scale(T x, T y, T z)
{
	return {
	    x,    T(0), T(0),  //
	    T(0), y,    T(0),  //
	    T(0), T(0), z,     //
	};
}
template <typename T>
matrix3_t<T> matrix3_t<T>::scale(T s)
{
	return scale(s, s, s);
}
template <typename T>
matrix3_t<T>& matrix3_t<T>::inverse()
{
	if (is_orthogonal())
	{
		return transpose();
	}

	T det = determinant();
	T recM00 = m00;
	T recM10 = m10;
	T recM20 = m20;
	T recM01 = m01;
	T recM11 = m11;
	T recM02 = m02;
	m00 = (recM11 * m22 - m21 * m12) / det;
	m10 = (m12 * recM20 - recM10 * m22) / det;
	m20 = (recM10 * m12 - recM20 * recM11) / det;
	m01 = (recM02 * m21 - recM01 * m22) / det;
	m11 = (recM00 * m22 - recM20 * recM02) / det;
	m21 = (recM01 * recM20 - recM00 * m21) / det;
	m02 = (recM01 * m12 - recM11 * recM02) / det;
	m12 = (recM02 * recM10 - recM00 * m12) / det;
	m22 = (recM00 * recM11 - recM10 * recM01) / det;

	return *this;
}
template <typename T>
matrix3_t<T> matrix3_t<T>::get_inversed() const
{
	matrix3_t<T> res = *this;
	res.inverse();
	return res;
}
template <typename T>
matrix3_t<T>& matrix3_t<T>::transpose()
{
	std::swap(m01, m10);
	std::swap(m02, m20);
	std::swap(m12, m21);
	return *this;
}
template <typename T>
matrix3_t<T> matrix3_t<T>::get_transposed() const
{
	return {
	    m00, m01, m02,  //
	    m10, m11, m12,  //
	    m20, m21, m22   //
	};
}

template <typename T>
T matrix3_t<T>::determinant() const
{
	return m00 * m11 * m22 + m01 * m12 * m20 + m02 * m10 * m21 -
	       m20 * m11 * m02 - m21 * m12 * m00 - m22 * m10 * m01;
}
template <typename T>
bool matrix3_t<T>::is_orthogonal() const
{
	return nearly_equal(m00 * m01 + m10 * m11 + m20 * m21, T(0)) &&
	       nearly_equal(m00 * m02 + m10 * m12 + m20 * m22, T(0)) &&
	       nearly_equal(m01 * m02 + m11 * m12 + m21 * m22, T(0)) &&
	       nearly_equal(m00 * m00 + m10 * m10 + m20 * m20, T(1)) &&
	       nearly_equal(m01 * m01 + m11 * m11 + m21 * m21, T(1)) &&
	       nearly_equal(m02 * m02 + m12 * m12 + m22 * m22, T(1));
}

template <typename T>
matrix3_t<T> matrix3_t<T>::operator*(T f) const
{
	return {
	    m00 * f, m10 * f, m20 * f,  //
	    m01 * f, m11 * f, m21 * f,  //
	    m02 * f, m12 * f, m22 * f   //
	};
}
template <typename T>
vector3_t<T> matrix3_t<T>::operator*(vector3_t<T> v) const
{
	return {m00 * v.x + m10 * v.y + m20 * v.z,
	        m01 * v.x + m11 * v.y + m21 * v.z,
	        m02 * v.x + m12 * v.y + m22 * v.z};
}
template <typename T>
matrix3_t<T> matrix3_t<T>::operator*(const matrix3_t<T>& m) const
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

#pragma region definition matrix4_t
template <typename T>
matrix4_t<T>::matrix4_t(T e00,
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
template <typename T>
matrix4_t<T>::matrix4_t(matrix2_t<T> m)
    : matrix4_t(m.m00,
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
template <typename T>
matrix4_t<T>::matrix4_t(const matrix3_t<T>& m)
    : matrix4_t(m.m00,
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
template <typename T>
matrix4_t<T>::matrix4_t(const perspective_t<T>& p) : matrix4_t(p.to_matrix4())
{
}
template <typename T>
matrix4_t<T>::matrix4_t(const transform3d_t<T>& t) : matrix4_t(t.to_matrix4())
{
}
template <typename T>
matrix4_t<T>::matrix4_t(const uniform_transform3d_t<T>& t)
    : matrix4_t(t.to_matrix4())
{
}
template <typename T>
matrix4_t<T>::matrix4_t(const unscaled_transform3d_t<T>& t)
    : matrix4_t(t.to_matrix4())
{
}
template <typename T>
matrix4_t<T>::matrix4_t(const std::array<T, 16>& a)
    : matrix4_t(a[0],
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

template <typename T>
matrix4_t<T> matrix4_t<T>::zero()
{
	return {T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0),
	        T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0)};
}
template <typename T>
matrix4_t<T> matrix4_t<T>::identity()
{
	return {T(1), T(0), T(0), T(0), T(0), T(1), T(0), T(0),
	        T(0), T(0), T(1), T(0), T(0), T(0), T(0), T(1)};
}
template <typename T>
matrix4_t<T> matrix4_t<T>::rows(vector4_t<T> row0,
                                vector4_t<T> row1,
                                vector4_t<T> row2,
                                vector4_t<T> row3)
{
	return {row0.x, row0.y, row0.z, row0.w, row1.x, row1.y, row1.z, row1.w,
	        row2.x, row2.y, row2.z, row2.w, row3.x, row3.y, row3.z, row3.w};
}
template <typename T>
matrix4_t<T> matrix4_t<T>::columns(vector4_t<T> column0,
                                   vector4_t<T> column1,
                                   vector4_t<T> column2,
                                   vector4_t<T> column3)
{
	return {column0.x, column1.x, column2.x, column3.x, column0.y, column1.y,
	        column2.y, column3.y, column0.z, column1.z, column2.z, column3.z,
	        column0.w, column1.w, column2.w, column3.w};
}
template <typename T>
matrix4_t<T> matrix4_t<T>::rotation(direction_t<T> axis, radian_t<T> angle)
{
	return matrix4_t<T>(matrix3_t<T>::rotation(axis, angle));
}
template <typename T>
matrix4_t<T> matrix4_t<T>::rotation(direction_t<T> axis,
                                    normalized<complex_t<T>> angle)
{
	return matrix4_t<T>(matrix3_t<T>::rotation(axis, angle));
}
template <typename T>
template <typename U, U NumeratorT, U DenominatorT>
matrix4_t<T> matrix4_t<T>::rotation(
    direction_t<T> axis,
    static_pi_fraction_t<U, NumeratorT, DenominatorT> angle)
{
	return matrix4_t<T>(matrix3_t<T>::rotation(axis, angle));
}
template <typename T>
matrix4_t<T> matrix4_t<T>::rotation(euler_angles_t<T> r)
{
	return matrix4_t<T>(matrix3_t<T>::rotation(r));
}
template <typename T>
matrix4_t<T> matrix4_t<T>::rotation(quaternion_t<T> q)
{
	return matrix4_t<T>(matrix3_t<T>::rotation(q));
}
template <typename T>
matrix4_t<T> matrix4_t<T>::translation(vector3_t<T> t)
{
	return {T(1), T(0), T(0), t.x, T(0), T(1), T(0), t.y,
	        T(0), T(0), T(1), t.z, T(0), T(0), T(0), T(1)};
}
template <typename T>
matrix4_t<T> matrix4_t<T>::translation(T x, T y, T z)
{
	return translation({x, y, z});
}
template <typename T>
matrix4_t<T> matrix4_t<T>::scale(vector3_t<T> s)
{
	return matrix4_t<T>(matrix3_t<T>::scale(s));
}
template <typename T>
matrix4_t<T> matrix4_t<T>::scale(T x, T y, T z)
{
	return matrix4_t<T>(matrix3_t<T>::scale(x, y, z));
}
template <typename T>
matrix4_t<T> matrix4_t<T>::scale(T s)
{
	return matrix4_t<T>(matrix3_t<T>::scale(s));
}
template <typename T>
matrix4_t<T> matrix4_t<T>::look_at(vector3_t<T> from,
                                   vector3_t<T> to,
                                   direction_t<T> up)
{
	direction_t<T> forward = (from - to);
	direction_t<T> right = cross(up, forward);
	direction_t<T> true_up = cross(forward, right);

	return {right->x, true_up->x, forward->x, from->x,
	        right->y, true_up->y, forward->y, from->y,
	        right->z, true_up->z, forward->z, from->z,
	        T(0),     T(0),       T(0),       T(1)};
}
template <typename T>
matrix4_t<T> matrix4_t<T>::orthographic(T left,
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
template <typename T>
matrix4_t<T> matrix4_t<T>::rotation_around_x(radian_t<T> angle)
{
	return matrix4_t<T>(matrix3_t<T>::rotation_around_x(angle));
}
template <typename T>
matrix4_t<T> matrix4_t<T>::rotation_around_y(radian_t<T> angle)
{
	return matrix4_t<T>(matrix3_t<T>::rotation_around_y(angle));
}
template <typename T>
matrix4_t<T> matrix4_t<T>::rotation_around_z(radian_t<T> angle)
{
	return matrix4_t<T>(matrix3_t<T>::rotation_around_z(angle));
}
template <typename T>
matrix4_t<T> matrix4_t<T>::rotation_around_x(normalized<complex_t<T>> angle)
{
	return matrix4_t<T>(matrix3_t<T>::rotation_around_x(angle));
}
template <typename T>
matrix4_t<T> matrix4_t<T>::rotation_around_y(normalized<complex_t<T>> angle)
{
	return matrix4_t<T>(matrix3_t<T>::rotation_around_y(angle));
}
template <typename T>
matrix4_t<T> matrix4_t<T>::rotation_around_z(normalized<complex_t<T>> angle)
{
	return matrix4_t<T>(matrix3_t<T>::rotation_around_z(angle));
}
template <typename T>
template <typename U, U NumeratorT, U DenominatorT>
matrix4_t<T> matrix4_t<T>::rotation_around_x(
    static_pi_fraction_t<U, NumeratorT, DenominatorT> angle)
{
	return matrix4_t<T>(matrix3_t<T>::rotation_around_x(angle));
}
template <typename T>
template <typename U, U NumeratorT, U DenominatorT>
matrix4_t<T> matrix4_t<T>::rotation_around_y(
    static_pi_fraction_t<U, NumeratorT, DenominatorT> angle)
{
	return matrix4_t<T>(matrix3_t<T>::rotation_around_y(angle));
}
template <typename T>
template <typename U, U NumeratorT, U DenominatorT>
matrix4_t<T> matrix4_t<T>::rotation_around_z(
    static_pi_fraction_t<U, NumeratorT, DenominatorT> angle)
{
	return matrix4_t<T>(matrix3_t<T>::rotation_around_z(angle));
}
template <typename T>
bool matrix4_t<T>::is_orthogonal() const
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
template <typename T>
bool matrix4_t<T>::is_homogenous() const
{
	return nearly_equal(m03, T(0)) && nearly_equal(m13, T(0)) &&
	       nearly_equal(m23, T(0)) && nearly_equal(m30, T(0)) &&
	       nearly_equal(m31, T(0)) && nearly_equal(m32, T(0)) &&
	       nearly_equal(m33, T(1));
}
template <typename T>
T matrix4_t<T>::determinant() const
{
	if (is_homogenous())
	{
		return m00 * m11 * m22 + m01 * m12 * m20 + m02 * m10 * m21 -
		       m20 * m11 * m02 - m21 * m12 * m00 - m22 * m10 * m01;
	}

	T det1 = m11 * (m22 * m33 - m32 * m23) - m21 * (m12 * m33 - m32 * m13) +
	         m31 * (m12 * m23 - m22 * m13);
	T det2 = m01 * (m22 * m33 - m32 * m23) - m21 * (m02 * m33 - m32 * m03) +
	         m31 * (m02 * m23 - m22 * m03);
	T det3 = m01 * (m12 * m33 - m32 * m13) - m11 * (m02 * m33 - m32 * m03) +
	         m31 * (m02 * m13 - m12 * m03);
	T det4 = m01 * (m12 * m23 - m22 * m13) - m11 * (m02 * m23 - m22 * m03) +
	         m21 * (m02 * m13 - m12 * m03);
	return m00 * det1 - m10 * det2 + m20 * det3 - m30 * det4;
}

template <typename T>
matrix4_t<T> matrix4_t<T>::operator*(T f) const
{
	return {m00 * f, m10 * f, m20 * f, m30 * f, m00 * f, m10 * f,
	        m20 * f, m30 * f, m00 * f, m10 * f, m20 * f, m30 * f,
	        m00 * f, m10 * f, m20 * f, m30 * f};
}
template <typename T>
matrix4_t<T> matrix4_t<T>::operator/(T f) const
{
	return {m00 / f, m10 / f, m20 / f, m30 / f, m00 / f, m10 / f,
	        m20 / f, m30 / f, m00 / f, m10 / f, m20 / f, m30 / f,
	        m00 / f, m10 / f, m20 / f, m30 / f};
}
template <typename T>
vector4_t<T> matrix4_t<T>::operator*(vector4_t<T> v) const
{
	return {m00 * v.x + m10 * v.y + m20 * v.z + m30 * v.w,
	        m01 * v.x + m11 * v.y + m21 * v.z + m31 * v.w,
	        m02 * v.x + m12 * v.y + m22 * v.z + m32 * v.w,
	        m03 * v.x + m13 * v.y + m23 * v.z + m33 * v.w};
}
template <typename T>
matrix4_t<T> matrix4_t<T>::operator*(const matrix4_t<T>& m) const
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

template <typename T>
matrix4_t<T>& matrix4_t<T>::inverse()
{
	if (is_orthogonal())
	{
		return transpose();
	}

	T det = determinant();
	T invDet = T(1) / det;
	T recM00 = m00;
	T recM10 = m10;
	T recM20 = m20;
	T recM01 = m01;
	T recM11 = m11;
	T recM02 = m02;

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

	T recM30 = m30;
	T recM21 = m21;
	T recM31 = m31;
	T recM12 = m12;
	T recM22 = m22;
	T recM03 = m03;
	T recM13 = m13;

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
template <typename T>
matrix4_t<T> matrix4_t<T>::get_inversed() const
{
	matrix4_t<T> m{*this};
	m.inverse();
	return m;
}
template <typename T>
matrix4_t<T>& matrix4_t<T>::transpose()
{
	std::swap(m01, m10);
	std::swap(m02, m20);
	std::swap(m03, m30);
	std::swap(m12, m21);
	std::swap(m13, m31);
	std::swap(m23, m32);
	return *this;
}
template <typename T>
matrix4_t<T> matrix4_t<T>::get_transposed() const
{
	matrix4_t<T> m{*this};
	m.transpose();
	return m;
}

template <typename T>
matrix4_t<T> inverse(matrix4_t<T> m)
{
	if (m.is_orthogonal())
	{
		return transpose(m);
	}

	T det = m.determinant();
	T invDet = T(1) / det;
	T recM00 = m.m00;
	T recM10 = m.m10;
	T recM20 = m.m20;
	T recM01 = m.m01;
	T recM11 = m.m11;
	T recM02 = m.m02;

	if (m.is_homogenous())
	{
		m.m00 = invDet * (recM11 * m.m22 - m.m21 * m.m12);
		m.m10 = invDet * (m.m12 * recM20 - recM10 * m.m22);
		m.m20 = invDet * (recM10 * m.m21 - recM20 * recM11);
		m.m01 = invDet * (recM02 * m.m21 - recM01 * m.m22);
		m.m11 = invDet * (recM00 * m.m22 - recM20 * recM02);
		m.m21 = invDet * (recM01 * recM20 - recM00 * m.m21);
		m.m02 = invDet * (recM01 * m.m12 - recM11 * recM02);
		m.m12 = invDet * (recM02 * recM10 - recM00 * m.m12);
		m.m22 = invDet * (recM00 * recM11 - recM10 * recM01);
		return m;
	}

	T recM30 = m.m30;
	T recM21 = m.m21;
	T recM31 = m.m31;
	T recM12 = m.m12;
	T recM22 = m.m22;
	T recM03 = m.m03;
	T recM13 = m.m13;

	m.m00 = invDet * (recM11 * recM22 * m.m33 - recM11 * m.m23 * m.m32 -
	                  recM21 * recM12 * m.m33 + recM21 * recM13 * m.m32 +
	                  recM31 * recM12 * m.m23 - recM31 * recM13 * recM22);

	m.m10 = invDet * (-recM10 * recM22 * m.m33 + recM10 * m.m23 * m.m32 +
	                  recM20 * recM12 * m.m33 - recM20 * recM13 * m.m32 -
	                  recM30 * recM12 * m.m23 + recM30 * recM13 * recM22);

	m.m20 = invDet * (recM10 * recM21 * m.m33 - recM10 * m.m23 * recM31 -
	                  recM20 * recM11 * m.m33 + recM20 * recM13 * recM31 +
	                  recM30 * recM11 * m.m23 - recM30 * recM13 * recM21);

	m.m30 = invDet * (-recM10 * recM21 * m.m32 + recM10 * recM22 * recM31 +
	                  recM20 * recM11 * m.m32 - recM20 * recM12 * recM31 -
	                  recM30 * recM11 * recM22 + recM30 * recM12 * recM21);

	m.m01 = invDet * (-recM01 * recM22 * m.m33 + recM01 * m.m23 * m.m32 +
	                  recM21 * recM02 * m.m33 - recM21 * recM03 * m.m32 -
	                  recM31 * recM02 * m.m23 + recM31 * recM03 * recM22);

	m.m11 = invDet * (recM00 * recM22 * m.m33 - recM00 * m.m23 * m.m32 -
	                  recM20 * recM02 * m.m33 + recM20 * recM03 * m.m32 +
	                  recM30 * recM02 * m.m23 - recM30 * recM03 * recM22);

	m.m21 = invDet * (-recM00 * recM21 * m.m33 + recM00 * m.m23 * recM31 +
	                  recM20 * recM01 * m.m33 - recM20 * recM03 * recM31 -
	                  recM30 * recM01 * m.m23 + recM30 * recM03 * recM21);

	m.m31 = invDet * (recM00 * recM21 * m.m32 - recM00 * recM22 * recM31 -
	                  recM20 * recM01 * m.m32 + recM20 * recM02 * recM31 +
	                  recM30 * recM01 * recM22 - recM30 * recM02 * recM21);

	m.m02 = invDet * (recM01 * recM12 * m.m33 - recM01 * recM13 * m.m32 -
	                  recM11 * recM02 * m.m33 + recM11 * recM03 * m.m32 +
	                  recM31 * recM02 * recM13 - recM31 * recM03 * recM12);

	m.m12 = invDet * (-recM00 * recM12 * m.m33 + recM00 * recM13 * m.m32 +
	                  recM10 * recM02 * m.m33 - recM10 * recM03 * m.m32 -
	                  recM30 * recM02 * recM13 + recM30 * recM03 * recM12);

	m.m22 = invDet * (recM00 * recM11 * m.m33 - recM00 * recM13 * recM31 -
	                  recM10 * recM01 * m.m33 + recM10 * recM03 * recM31 +
	                  recM30 * recM01 * recM13 - recM30 * recM03 * recM11);

	m.m32 = invDet * (-recM00 * recM11 * m.m32 + recM00 * recM12 * recM31 +
	                  recM10 * recM01 * m.m32 - recM10 * recM02 * recM31 -
	                  recM30 * recM01 * recM12 + recM30 * recM02 * recM11);

	m.m03 = invDet * (-recM01 * recM12 * m.m23 + recM01 * recM13 * recM22 +
	                  recM11 * recM02 * m.m23 - recM11 * recM03 * recM22 -
	                  recM21 * recM02 * recM13 + recM21 * recM03 * recM12);

	m.m13 = invDet * (recM00 * recM12 * m.m23 - recM00 * recM13 * recM22 -
	                  recM10 * recM02 * m.m23 + recM10 * recM03 * recM22 +
	                  recM20 * recM02 * recM13 - recM20 * recM03 * recM12);

	m.m23 = invDet * (-recM00 * recM11 * m.m23 + recM00 * recM13 * recM21 +
	                  recM10 * recM01 * m.m23 - recM10 * recM03 * recM21 -
	                  recM20 * recM01 * recM13 + recM20 * recM03 * recM11);

	m.m33 = invDet * (recM00 * recM11 * recM22 - recM00 * recM12 * recM21 -
	                  recM10 * recM01 * recM22 + recM10 * recM02 * recM21 +
	                  recM20 * recM01 * recM12 - recM20 * recM02 * recM11);

	return *this;
}
template <typename T>
matrix4_t<T> transpose(matrix4_t<T> m)
{
	std::swap(m.m01, m.m10);
	std::swap(m.m02, m.m20);
	std::swap(m.m03, m.m30);
	std::swap(m.m12, m.m21);
	std::swap(m.m13, m.m31);
	std::swap(m.m23, m.m32);
	return m;
}
#pragma endregion

#pragma region definition perspective_t
template <typename T>
perspective_t<T>::perspective_t(radian_t<T> angle, T ratio, T near, T far)
    : m_angle{angle},
      m_ratio{ratio},
      m_invRatio{T(1) / m_ratio},
      m_near{near},
      m_far{far},
      m_invTanHalfFovy{T(1) / tan(m_angle / T(2))}
{
}
template <typename T>
void perspective_t<T>::set(radian_t<T> angle, T ratio, T near, T far)
{
	m_angle = angle;
	m_ratio = ratio;
	m_invRatio = T(1) / m_ratio;
	m_near = near;
	m_far = far;
	m_invTanHalfFovy = T(1) / tan(m_angle / T(2));
}
template <typename T>
void perspective_t<T>::set_angle(radian_t<T> angle)
{
	m_angle = angle;
	m_invTanHalfFovy = T(1) / tan(angle / T(2));
}
template <typename T>
radian_t<T> perspective_t<T>::get_angle() const
{
	return m_angle;
}
template <typename T>
void perspective_t<T>::set_ratio(T ratio)
{
	m_ratio = ratio;
	m_invRatio = T(1) / m_ratio;
}
template <typename T>
T perspective_t<T>::get_ratio() const
{
	return m_ratio;
}
template <typename T>
void perspective_t<T>::set_near(T near_plane_t)
{
	m_near = near_plane_t;
}
template <typename T>
T perspective_t<T>::get_near() const
{
	return m_near;
}
template <typename T>
void perspective_t<T>::set_far(T far_plane_t)
{
	m_far = far_plane_t;
}
template <typename T>
T perspective_t<T>::get_far() const
{
	return m_far;
}
template <typename T>
matrix4_t<T> perspective_t<T>::to_matrix4() const
{
	matrix4_t<T> res{matrix4_t<T>::zero()};
	res.column(0).row(0) = m_invRatio * m_invTanHalfFovy;
	res.column(1).row(1) = -m_invTanHalfFovy;
	res.column(2).row(2) = m_far / (m_near - m_far);
	res.column(2).row(3) = -T(1);
	res.column(3).row(2) = -(m_far * m_near) / (m_far - m_near);

	return res;
}

template <typename T>
matrix4_t<T> operator*(const matrix4_t<T>& m, const perspective_t<T>& p)
{
	T a = p.m_invRatio * p.m_invTanHalfFovy;
	T b = -p.m_invTanHalfFovy;
	T c = p.m_far / (p.m_near - p.m_far);
	T d = -(p.m_far * p.m_near) / (p.m_far - p.m_near);
	return {a * m.m00, b * m.m10, c * m.m20 - m.m30, d * m.m20,

	        a * m.m01, b * m.m11, c * m.m21 - m.m31, d * m.m21,

	        a * m.m02, b * m.m12, c * m.m22 - m.m32, d * m.m22,

	        a * m.m03, b * m.m13, c * m.m23 - m.m33, d * m.m23};
}
template <typename T>
matrix4_t<T> operator*(const perspective_t<T>& p, const matrix4_t<T>& m)
{
	T a = p.m_invRatio * p.m_invTanHalfFovy;
	T b = -p.m_invTanHalfFovy;
	T c = p.m_far / (p.m_near - p.m_far);
	T d = -(p.m_far * p.m_near) / (p.m_far - p.m_near);
	return {m.m00 * a,
	        m.m10 * a,
	        m.m20 * a,
	        m.m30 * a,

	        m.m01 * b,
	        m.m11 * b,
	        m.m21 * b,
	        m.m31 * b,

	        m.m02 * c + m.m03 * d,
	        m.m12 * c + m.m13 * d,
	        m.m22 * c + m.m23 * d,
	        m.m32 * c + m.m33 * d,

	        -m.m02,
	        -m.m12,
	        -m.m22,
	        -m.m32};
}

template <typename T>
matrix4_t<T> operator*(const unscaled_transform3d_t<T>& t,
                       const perspective_t<T>& p)
{
	T a = p.m_invRatio * p.m_invTanHalfFovy;
	T b = -p.m_invTanHalfFovy;
	T c = p.m_far / (p.m_near - p.m_far);
	T d = -(p.m_far * p.m_near) / (p.m_far - p.m_near);
	T xx = t.rotation.x * t.rotation.x;
	T yy = t.rotation.y * t.rotation.y;
	T zz = t.rotation.z * t.rotation.z;
	T wx = t.rotation.w * t.rotation.x;
	T wy = t.rotation.w * t.rotation.y;
	T wz = t.rotation.w * t.rotation.z;
	T xy = t.rotation.x * t.rotation.y;
	T xz = t.rotation.x * t.rotation.z;
	T yz = t.rotation.y * t.rotation.z;
	return {a * (T(1) - T(2) * (yy + zz)),
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
	        -T(1),
	        T(0)};
}
template <typename T>
matrix4_t<T> operator*(const perspective_t<T>& p,
                       const unscaled_transform3d_t<T>& t)
{
	T a = p.m_invRatio * p.m_invTanHalfFovy;
	T b = -p.m_invTanHalfFovy;
	T c = p.m_far / (p.m_near - p.m_far);
	T d = -(p.m_far * p.m_near) / (p.m_far - p.m_near);
	T xx = t.rotation.x * t.rotation.x;
	T yy = t.rotation.y * t.rotation.y;
	T zz = t.rotation.z * t.rotation.z;
	T wx = t.rotation.w * t.rotation.x;
	T wy = t.rotation.w * t.rotation.y;
	T wz = t.rotation.w * t.rotation.z;
	T xy = t.rotation.x * t.rotation.y;
	T xz = t.rotation.x * t.rotation.z;
	T yz = t.rotation.y * t.rotation.z;

	matrix4_t<T> res;
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

#pragma region definition quaternion_t
template <typename T>
quaternion_t<T>::quaternion_t(direction_t<T> axis, radian_t<T> angle)
{
	vector3_t<T> tmpAxis = *axis * sin(angle / 2.0f);
	*this = {tmpAxis.x, tmpAxis.y, tmpAxis.z, cos(angle / 2.0f)};
}
template <typename T>
template <typename U, U NumeratorT, U DenominatorT>
quaternion_t<T>::quaternion_t(
    direction_t<T> axis,
    static_pi_fraction_t<U, NumeratorT, DenominatorT> angle)
{
	vector3_t<T> tmpAxis =
	    *axis * sin<T>(static_pi_fraction_t<U, angle.numerator,
	                                        angle.denominator * U(2)>{});
	*this = {tmpAxis.x, tmpAxis.y, tmpAxis.z,
	         cos<T>(static_pi_fraction_t<U, angle.numerator,
	                                     angle.denominator * U(2)>{})};
}

template <typename T>
quaternion_t<T> quaternion_t<T>::zero()
{
	return {T(0), T(0), T(0), T(0)};
}
template <typename T>
quaternion_t<T> quaternion_t<T>::identity()
{
	return {T(0), T(0), T(0), T(1)};
}

template <typename T>
T quaternion_t<T>::norm() const
{
	return std::sqrt(norm_squared());
}
template <typename T>
T quaternion_t<T>::norm_squared() const
{
	return x * x + y * y + z * z + w * w;
}

template <typename T>
quaternion_t<T>& quaternion_t<T>::inverse()
{
	T n = norm();
	if (nearly_equal(n, T(1)))
		return conjugate();

	x /= -n;
	y /= -n;
	z /= -n;
	w /= n;
	return *this;
}
template <typename T>
quaternion_t<T>& quaternion_t<T>::conjugate()
{
	x = -x;
	y = -y;
	z = -z;
	return *this;
}
template <typename T>
quaternion_t<T>& quaternion_t<T>::normalize()
{
	T n = norm();
	x /= n;
	y /= n;
	z /= n;
	w /= n;
	return *this;
}
template <typename T>
quaternion_t<T>& quaternion_t<T>::clamp(quaternion_t<T> min,
                                        quaternion_t<T> max)
{
	x = cdm::clamp(x, min.x, max.x);
	y = cdm::clamp(y, min.y, max.y);
	z = cdm::clamp(z, min.z, max.z);
	w = cdm::clamp(w, min.w, max.w);
	return *this;
}

template <typename T>
quaternion_t<T> quaternion_t<T>::get_inversed() const
{
	quaternion_t<T> res = *this;
	res.inverse();
	return res;
}
template <typename T>
quaternion_t<T> quaternion_t<T>::get_conjugated() const
{
	quaternion_t<T> res = *this;
	res.conjugate();
	return res;
}
template <typename T>
quaternion_t<T> quaternion_t<T>::get_normalized() const
{
	quaternion_t<T> res = *this;
	res.normalize();
	return res;
}
template <typename T>
quaternion_t<T> quaternion_t<T>::get_clamped(quaternion_t<T> min,
                                             quaternion_t<T> max) const
{
	quaternion_t<T> res = *this;
	clamp(res, min, max);
	return res;
}

template <typename T>
quaternion_t<T> quaternion_t<T>::operator+(quaternion_t<T> q) const
{
	return {x + q.x, y + q.y, z + q.z, w + q.w};
}
template <typename T>
quaternion_t<T> quaternion_t<T>::operator-(quaternion_t<T> q) const
{
	return {x - q.x, y - q.y, z - q.z, w - q.w};
}
template <typename T>
vector3_t<T> quaternion_t<T>::operator*(vector3_t<T> v) const
{
	vector3_t<T> qvec{x, y, z};
	vector3_t<T> uv{cross(qvec, v)};
	vector3_t<T> uuv{cross(qvec, uv)};
	uv = uv * (T(2) * w);
	uuv = uuv * T(2);
	return v + uv + uuv;
}
template <typename T>
quaternion_t<T> quaternion_t<T>::operator*(quaternion_t<T> q) const
{
	quaternion_t<T> res;
	res.x = x * q.w + y * q.z - z * q.y + w * q.x;
	res.y = -x * q.z + y * q.w + z * q.x + w * q.y;
	res.z = x * q.y - y * q.x + z * q.w + w * q.z;
	res.w = -x * q.x - y * q.y - z * q.z + w * q.w;
	return res;
}
template <typename T>
quaternion_t<T> quaternion_t<T>::operator*(T f) const
{
	return {x * f, y * f, z * f, w * f};
}
template <typename T>
quaternion_t<T> quaternion_t<T>::operator/(T f) const
{
	return {x / f, y / f, z / f, w / f};
}

template <typename T>
quaternion_t<T>& quaternion_t<T>::operator+=(quaternion_t<T> q)
{
	*this = *this + q;
	return *this;
}
template <typename T>
quaternion_t<T>& quaternion_t<T>::operator-=(quaternion_t<T> q)
{
	*this = *this - q;
	return *this;
}
template <typename T>
quaternion_t<T>& quaternion_t<T>::operator*=(quaternion_t<T> q)
{
	*this = *this * q;
	return *this;
}
template <typename T>
quaternion_t<T>& quaternion_t<T>::operator*=(T f)
{
	*this = *this * f;
	return *this;
}
template <typename T>
quaternion_t<T>& quaternion_t<T>::operator/=(T f)
{
	*this = *this / f;
	return *this;
}

template <typename T>
quaternion_t<T> quaternion_t<T>::operator-() const
{
	return {-x, -y, -z, -w};
}

template <typename T>
bool quaternion_t<T>::operator==(quaternion_t<T> q) const
{
	return x == q.x && y == q.y && z == q.z && w == q.w;
}
template <typename T>
bool quaternion_t<T>::operator!=(quaternion_t<T> q) const
{
	return !operator==(q);
}

template <typename T>
vector3_t<T> operator*(const normalized<quaternion_t<T>>& q, vector3_t<T> v)
{
	vector3_t<T> qvec = {q->x, q->y, q->z};
	vector3_t<T> uv = cross(qvec, v);
	vector3_t<T> uuv = cross(qvec, uv);
	uv = uv * (T(2) * q->w);
	uuv = uuv * T(2);
	return v + uv + uuv;
}

template <typename T>
T dot(quaternion_t<T> lhs, quaternion_t<T> rhs)
{
	return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w;
}
template <typename T>
quaternion_t<T> cross(quaternion_t<T> lhs, quaternion_t<T> rhs)
{
	return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z,
	        lhs.x * rhs.y - lhs.y * rhs.x, T(1)};
}
template <typename T>
quaternion_t<T> lerp(quaternion_t<T> begin, quaternion_t<T> end, T percent)
{
	return (end - begin) * percent + begin;
}
template <typename T>
quaternion_t<T> nlerp(quaternion_t<T> begin, quaternion_t<T> end, T percent)
{
	return lerp(begin, end, percent).get_normalized();
}
template <typename T>
quaternion_t<T> slerp(quaternion_t<T> begin, quaternion_t<T> end, T percent)
{
	quaternion_t<T> _end;
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
	d = cdm::clamp(d, -T(1), T(1));
	T theta = acosf(d) * percent;

	quaternion_t<T> res = begin * d;
	res = _end - res;

	res = begin * cosf(theta) + res * sinf(theta);

	return res;
}
#pragma endregion

#pragma region definition cartesian_direction2d_t
template <typename T>
cartesian_direction2d_t<T> cartesian_direction2d_t<T>::from(
    const normalized<vector2_t<T>>& direction)
{
	cartesian_direction2d_t<T> d;
	d.x = direction->x;
	d.y = direction->y;
	return d;
}
template <typename T>
cartesian_direction2d_t<T> cartesian_direction2d_t<T>::from(
    polar_direction_t<T> direction)
{
	return from(direction.angle);
}
template <typename T>
cartesian_direction2d_t<T> cartesian_direction2d_t<T>::from(radian_t<T> angle)
{
	cartesian_direction2d_t<T> d;
	d.x = cos(angle);
	d.y = sin(angle);
	return d;
}

template <typename T>
cartesian_direction2d_t<T>& cartesian_direction2d_t<T>::rotate(
    radian_t<T> angle)
{
	const T c = cos(angle);
	const T s = sin(angle);
	x = x * c - y * s;
	y = x * s + y * c;

	return *this;
}
#pragma endregion

#pragma region definition polar_direction_t
template <typename T>
polar_direction_t<T> polar_direction_t<T>::from(
    const normalized<vector2_t<T>>& direction)
{
	angle = radian_t<T>(atan(direction->x));
}
template <typename T>
polar_direction_t<T> polar_direction_t<T>::from(
    cartesian_direction2d_t<T> direction)
{
	angle = radian_t<T>(atan(direction.x));
}
template <typename T>
polar_direction_t<T> polar_direction_t<T>::from(radian_t<T> angle_)
{
	angle = angle_;
}

template <typename T>
polar_direction_t<T>& polar_direction_t<T>::rotate(radian_t<T> angle_)
{
	angle += angle_;
	return *this;
}
template <typename T>
polar_direction_t<T>& polar_direction_t<T>::wrap()
{
	int s = sign(T(angle));
	while (angle > radian_t(pi) || angle <= -radian_t(pi))
		angle -= radian_t(pi) * T(s);
	return *this;
}
#pragma endregion

#pragma region definition line_t
template <typename T>
bool are_parallel(line_t<T> l1, line_t<T> l2)
{
	return nearly_equal(l1.coefficient, l2.coefficient) ||
	       nearly_equal(l1.coefficient, -l2.coefficient);
}

template <typename T>
bool collides(line_t<T> l1, line_t<T> l2)
{
	return !line_t<T>::are_parallel(l1, l2);
}
template <typename T>
bool collides(line_t<T> l1, line_t<T> l2, vector2_t<T>& intersection)
{
	bool collision = collides(l1, l2);

	if (collision)
	{
		T a = l1.coefficient - l2.coefficient;
		T b = l2.offset - l1.offset;
		intersection.x = b / a;
		intersection.y = l1.coefficient * intersection.x + l1.offset;

		// assert(intersection.y == (l2.coefficient * intersection.x +
		// l2.offset));
	}

	return collision;
}
template <typename T>
bool collides(line_t<T> l, vector2_t<T> v)
{
	return nearly_equal(l.coefficient * v.x + l.offset, v.y);
}
template <typename T>
bool collides(vector2_t<T> v, line_t<T> l)
{
	return collides(l, v);
}

template <typename T>
T distance_between(plane_t<T> p, vector3_t<T> v)
{
	return -v.x * p.normal->x - v.y * p.normal->y - v.z * p.normal->z;
}
template <typename T>
T distance_between(vector3_t<T> v, plane_t<T> p)
{
	return distance_between(p, v);
}
#pragma endregion

#pragma region definition aa_rect_t
template <typename T>
bool aa_rect_t<T>::contains(vector2_t<T> v) const
{
	return (v.x >= origin.x) && (v.x <= origin.x + dimention.x) &&
	       (v.y >= origin.y) && (v.y <= origin.y + dimention.y);
}

template <typename T>
bool collides(aa_rect_t<T> r1, aa_rect_t<T> r2)
{
	return collides(r1, r2.origin) ||
	       collides(r1, r2.origin + vector2_t(r2.dimention.x, 0)) ||
	       collides(r1, r2.origin + vector2_t(0, r2.dimention.y)) ||
	       collides(r1, r2.origin + r2.dimention);
}

// https://stackoverflow.com/a/565282
template <typename T>
bool intersects(segment2d_t<T> s0,
                segment2d_t<T> s1,
                vector2_t<T>& outPoint,
                T e = T(epsilon))
{
	vector2_t<T> p = s0.origin;
	vector2_t<T> q = s1.origin;
	vector2_t<T> r = from_to(p, s0.end);
	vector2_t<T> s = from_to(q, s1.end);

	vector2_t<T> qmp = q - p;

	T rcs = cross(r, s);

	if (nearly_equal(rcs, T(0), e))  // coline_tar
	{
		if (nearly_equal(cross(qmp, r), T(0),
		                 e))  // same line_t, different segments
		{
			T invrdr = T(1) / dot(r, r);
			T t0 = dot(qmp, r) * invrdr;
			T t1 = dot(qmp + s, r) * invrdr;

			if (dot(s, r) < T(0))
				std::swap(t0, t1);

			if (t0 <= T(1) && t0 >= T(0) || t1 <= T(1) && t1 >= T(0))
			{
				outPoint = s0.origin;  /// TODO: better than that
				return true;
			}
			else
				return false;
		}
		else  // parallel, non-intersecting
		{
			return false;
		}
	}
	else
	{
		T invrcs = T(1) / rcs;

		T t = cross(qmp, s) * invrcs;
		T u = cross(qmp, r) * invrcs;

		if (t >= T(0) && t <= T(1) && u >= T(0) && u <= T(1))  // intersects
		{
			outPoint = p + r * t;

			return true;
		}
		else
		{
			return false;
		}
	}
}
#pragma endregion

#pragma region definition plane_t
template <typename T>
T plane_t<T>::evaluate(vector3_t<T> point) const
{
	return normal->x * (point.x - origin.x) +
	       normal->y * (point.y - origin.y) + normal->z * (point.z - origin.z);
}
template <typename T>
vector3_t<T> plane_t<T>::project3d(vector3_t<T> point) const
{
	T distance = evaluate(point);
	return point - *normal * distance;
}
template <typename T>
vector2_t<T> plane_t<T>::project2d(vector3_t<T> point,
                                   direction_t<T> plane_t_tangent) const
{
	direction_t<T> bitangent = cross(normal, plane_t_tangent);
	direction_t<T> tangent = cross(bitangent, normal);

	matrix3_t<T> TBN{matrix3_t<T>::rows(*tangent, *bitangent, *normal)};

	vector3_t<T> plane_t_space_point = point - origin;

	return (TBN * plane_t_space_point).xy();
}
template <typename T>
vector3_t<T> plane_t<T>::unproject(vector2_t<T> point,
                                   direction_t<T> plane_t_tangent) const
{
	direction_t<T> bitangent = cross(normal, plane_t_tangent);
	direction_t<T> tangent = cross(bitangent, normal);

	matrix3_t<T> invTBN{
	    matrix3_t<T>::rows(*tangent, *bitangent, *normal).get_inversed()};

	vector3_t<T> plane_t_space_point = vector3_t<T>{point.x, point.y, T(0)};

	return (invTBN * plane_t_space_point) + origin;
}
template <typename T>
direction_t<T> plane_t<T>::computeTangent(T e) const
{
	vector3_t<T> tangent =
	    project3d(vector3_t{T(1), T(0), T(0)} + origin) - origin;
	if (norm_squared(tangent) < e)
		tangent = project3d(vector3_t{T(0), T(1), T(0)} + origin) - origin;
	return direction_t<T>(tangent);
}

template <typename T>
bool collides(const plane_t<T>& plane_t,
              ray3d_t<T> r,
              vector3_t<T>& collision_point)
{
	T denom = dot(*plane_t.normal, r.direction);
	if (abs(denom) > T(0.0001))
	{
		T t = -dot(plane_t.origin - r.origin, plane_t.normal) / denom;
		if (t >= T(0))
		{
			collision_point = r.origin + r.direction * t;
			return true;
		}
	}
	return false;
}
template <typename T>
bool collides_bidirectional(const plane_t<T>& plane,
                            ray3d_t<T> r,
                            vector3_t<T>& collision_point)
{
	T denom = dot(*plane.normal, *r.direction);
	if (abs(denom) > T(0.0001))
	{
		T t = -dot(plane.origin - r.origin, *plane.normal) / denom;

		collision_point = r.origin + *r.direction * t;
		return true;
	}
	return false;
}
#pragma endregion

#pragma region definition aabb_t
template <typename T>
bool aabb_t<T>::contains(vector3_t<T> p) const
{
	return p.x >= min.x && p.x <= max.x && p.y >= min.y && p.y <= max.y &&
	       p.z >= min.z && p.z <= max.z;
}
template <typename T>
vector3_t<T> aabb_t<T>::get_center() const
{
	return (max + min) / T(2);
}

template <typename T>
aabb_t<T> aabb_t<T>::operator+(aabb_t<T> rhs) const
{
	aabb_t<T> res;
	res.min.x = std::min(min.x, rhs.min.x);
	res.min.y = std::min(min.y, rhs.min.y);
	res.min.z = std::min(min.z, rhs.min.z);
	res.max.x = std::max(max.x, rhs.max.x);
	res.max.y = std::max(max.y, rhs.max.y);
	res.max.z = std::max(max.z, rhs.max.z);

	return res;
}

template <typename T>
bool collides(aabb_t<T> b, ray3d_t<T> r)
{
	vector3_t<T> inv{T(1) / r.direction->x, T(1) / r.direction->y,
	                 T(1) / r.direction->z};

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
template <typename T>
bool collides(aabb_t<T> b1, aabb_t<T> b2)
{
	if (b1.contains(b2.min))
		return true;
	if (b1.contains(b2.max))
		return true;
	if (b2.contains(b1.min))
		return true;
	if (b2.contains(b1.max))
		return true;

	std::array<vector3_t<T>, 6> otherPoints;
	otherPoints[0] = {b1.min.x, b1.min.y, b1.max.z};
	otherPoints[1] = {b1.min.x, b1.max.y, b1.min.z};
	otherPoints[2] = {b1.max.x, b1.min.y, b1.min.z};
	otherPoints[3] = {b1.min.x, b1.max.y, b1.max.z};
	otherPoints[4] = {b1.max.x, b1.min.y, b1.max.z};
	otherPoints[5] = {b1.max.x, b1.max.y, b1.min.z};

	for (auto& p : otherPoints)
		if (b2.contains(p))
			return true;

	otherPoints[0] = {b2.min.x, b2.min.y, b2.max.z};
	otherPoints[1] = {b2.min.x, b2.max.y, b2.min.z};
	otherPoints[2] = {b2.max.x, b2.min.y, b2.min.z};
	otherPoints[3] = {b2.min.x, b2.max.y, b2.max.z};
	otherPoints[4] = {b2.max.x, b2.min.y, b2.max.z};
	otherPoints[5] = {b2.max.x, b2.max.y, b2.min.z};

	for (auto& p : otherPoints)
		if (b1.contains(p))
			return true;

	return false;
}
#pragma endregion

#pragma region definition transform2d_t
template <typename T>
transform2d_t<T> transform2d_t<T>::operator*(const transform2d_t<T>& t) const
{
	transform2d_t<T> res;
	res.position =
	    rotation * vector2_t(scale.x * t.position.x, scale.y * t.position.y) +
	    position;
	res.rotation = rotation + t.rotation;
	res.scale = vector2_t(scale.x * t.scale.x, scale.y * t.scale.y);
	return res;
}
template <typename T>
vector2_t<T> transform2d_t<T>::operator*(vector2_t<T> v) const
{
	return rotation * vector2_t(scale.x * v.x, scale.y * v.y) + position;
}
#pragma endregion

#pragma region definition transform3d_t
template <typename T>
transform3d_t<T> transform3d_t<T>::operator*(const transform3d_t<T>& t) const
{
	transform3d_t<T> res;
	res.position =
	    rotation * vector3_t(scale.x * t.position.x, scale.y * t.position.y,
	                         scale.z * t.position.z) +
	    position;
	res.rotation = rotation * t.rotation;
	res.scale = vector3_t(scale.x * t.scale.x, scale.y * t.scale.y,
	                      scale.z * t.scale.z);
	return res;
}
template <typename T>
vector3_t<T> transform3d_t<T>::operator*(vector3_t<T> v) const
{
	return rotation * vector3_t(scale.x * v.x, scale.y * v.y, scale.z * v.z) +
	       position;
}
template <typename T>
quaternion_t<T> transform3d_t<T>::operator*(quaternion_t<T> q) const
{
	return rotation * q;
}
#pragma endregion

#pragma region definition uniform_transform2d_t
template <typename T>
uniform_transform2d_t<T> uniform_transform2d_t<T>::operator*(
    uniform_transform2d_t<T> t) const
{
	uniform_transform2d_t<T> res;
	matrix2_t r = matrix2_t<T>::rotation(rotation);
	res.position = r * (scale * t.position) + position;
	res.rotation = rotation + t.rotation;
	res.scale = scale * t.scale;
	return res;
}
template <typename T>
vector2_t<T> uniform_transform2d_t<T>::operator*(vector2_t<T> v) const
{
	matrix2_t r = matrix2_t::rotation(rotation);
	return r * (scale * v) + position;
}
#pragma endregion

#pragma region definition uniform_transform3d_t
template <typename T>
uniform_transform3d_t<T> uniform_transform3d_t<T>::operator*(
    const uniform_transform3d_t<T>& t) const
{
	uniform_transform3d_t<T> res;
	res.position = rotation * (scale * t.position) + position;
	res.rotation = rotation * t.rotation;
	res.scale = scale * t.scale;
	return res;
}
template <typename T>
vector3_t<T> uniform_transform3d_t<T>::operator*(vector3_t<T> v) const
{
	return rotation * (scale * v) + position;
}
template <typename T>
quaternion_t<T> uniform_transform3d_t<T>::operator*(quaternion_t<T> q) const
{
	return rotation * q;
}
#pragma endregion

#pragma region definition unscaled_transform2d_t
template <typename T>
unscaled_transform2d_t<T> unscaled_transform2d_t<T>::operator*(
    unscaled_transform2d_t<T> t) const
{
	unscaled_transform2d_t<T> res;
	matrix2_t r = matrix2_t<T>::rotation(rotation);
	res.position = r * t.position + position;
	res.rotation = rotation + t.rotation;
	return res;
}
template <typename T>
vector2_t<T> unscaled_transform2d_t<T>::operator*(vector2_t<T> v) const
{
	matrix2_t<T> r = matrix2_t::rotation(rotation);
	return r * v + position;
}
#pragma endregion

#pragma region definition unscaled_transform3d_t
template <typename T>
unscaled_transform3d_t<T>& unscaled_transform3d_t<T>::translate_absolute(
    vector3_t<T> t)
{
	position += t;
	return *this;
}

template <typename T>
unscaled_transform3d_t<T>& unscaled_transform3d_t<T>::translate_relative(
    vector3_t<T> t)
{
	position += rotation * t;
	return *this;
}

template <typename T>
unscaled_transform3d_t<T>& unscaled_transform3d_t<T>::rotate(quaternion_t<T> r)
{
	rotation = r * rotation;
	return *this;
}

template <typename T>
unscaled_transform3d_t<T> inverse(unscaled_transform3d_t<T> t)
{
	t.rotation = inverse(t.rotation);
	t.position = t.rotation * -t.position;
	return t;
}

template <typename T>
unscaled_transform3d_t<T> unscaled_transform3d_t<T>::operator*(
    const unscaled_transform3d_t<T>& t) const
{
	unscaled_transform3d_t<T> res;
	res.position = rotation * t.position + position;
	res.rotation = rotation * t.rotation;
	return res;
}
template <typename T>
vector3_t<T> unscaled_transform3d_t<T>::operator*(vector3_t<T> v) const
{
	return rotation * v + position;
}
template <typename T>
quaternion_t<T> unscaled_transform3d_t<T>::operator*(quaternion_t<T> q) const
{
	return rotation * q;
}
#pragma endregion

#pragma region definition streams
template <typename T>
std::ostream& operator<<(std::ostream& o, const radian_t<T>& a)
{
	return o << "radian_t(" << static_cast<T>(a) << ")";
}
template <typename T>
std::ostream& operator<<(std::ostream& o, const degree_t<T>& a)
{
	return o << "degree_t(" << static_cast<T>(a) << ")";
}
template <typename T>
std::ostream& operator<<(std::ostream& o, const pi_fraction_t<T>& f)
{
	return o << "pi_fraction_t(" << f.numerator << "pi / " << f.denominator
	         << ")";
}
template <typename T, T NumT, T DenT>
std::ostream& operator<<(std::ostream& o,
                         const static_pi_fraction_t<T, NumT, DenT>& f)
{
	return o << "static_pi_fraction_t(" << f.numerator << "pi / "
	         << f.denominator << ")";
}
template <typename T>
std::ostream& operator<<(std::ostream& os, vector2_t<T> v)
{
	return os << "vector2_t(" << v.x << ", " << v.y << ")";
}
template <typename T>
std::ostream& operator<<(std::ostream& os, vector3_t<T> v)
{
	return os << "vector3_t(" << v.x << ", " << v.y << ", " << v.z << ")";
}
template <typename T>
std::ostream& operator<<(std::ostream& os, vector4_t<T> v)
{
	return os << "vector4_t(" << v.x << ", " << v.y << ", " << v.z << ", "
	          << v.w << ")";
}
template <typename T>
std::ostream& operator<<(std::ostream& os, quaternion_t<T> q)
{
	return os << "quaternion_t({" << q.x << ", " << q.y << ", " << q.z << "}, "
	          << q.w << ")";
}
template <typename T>
std::ostream& operator<<(std::ostream& os, plane_t<T> p)
{
	return os << "plane_t(origin = " << p.origin << ", normal = " << p.normal
	          << ")";
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const cdm::matrix2_t<T>& m)
{
	return os << "matrix2_t(" << m.m00 << "\t" << m.m10 << "\n        "
	          << m.m01 << "\t" << m.m11 << ")";
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const cdm::matrix3_t<T>& m)
{
	return os << "matrix3_t(" << m.column(0).row(0) << "\t"
	          << m.column(1).row(0) << "\t" << m.column(2).row(0)
	          << "\n        " << m.column(0).row(1) << "\t"
	          << m.column(1).row(1) << "\t" << m.column(2).row(1)
	          << "\n        " << m.column(0).row(2) << "\t"
	          << m.column(1).row(2) << "\t" << m.column(2).row(2) << ")";
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const cdm::matrix4_t<T>& m)
{
	return os << "matrix4_t(" << m.column(0).row(0) << "\t"
	          << m.column(1).row(0) << "\t" << m.column(2).row(0) << "\t"
	          << m.column(3).row(0) << "\n        " << m.column(0).row(1)
	          << "\t" << m.column(1).row(1) << "\t" << m.column(2).row(1)
	          << "\t" << m.column(3).row(1) << "\n        "
	          << m.column(0).row(2) << "\t" << m.column(1).row(2) << "\t"
	          << m.column(2).row(2) << "\t" << m.column(3).row(2)
	          << "\n        " << m.column(0).row(3) << "\t"
	          << m.column(1).row(3) << "\t" << m.column(2).row(3) << "\t"
	          << m.column(3).row(3) << ")";
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const cdm::normalized<T>& n)
{
	return os << *n;
}
#pragma endregion

#pragma region definition misc
namespace detail
{
template <unsigned char x, unsigned char y, typename T>
T get_quaternion_t_matrix_element(quaternion_t<T> q)
{
	static_assert(x < 3, "x must be [0;2]");
	static_assert(y < 3, "y must be [0;2]");
	if constexpr (y == 0)
	{
		if constexpr (x == 0)
			return T(1) - T(2) * (q.y * q.y + q.z * q.z);
		if constexpr (x == 1)
			return T(2) * (q.x * q.y - q.z * q.w);
		else
			return T(2) * (q.x * q.z + q.y * q.w);
	}
	if constexpr (y == 1)
	{
		if constexpr (x == 0)
			return T(2) * (q.x * q.y + q.z * q.w);
		if constexpr (x == 1)
			return T(1) - T(2) * (q.x * q.x + q.z * q.z);
		else
			return T(2) * (q.y * q.z - q.x * q.w);
	}
	else
	{
		if constexpr (x == 0)
			return T(2) * (q.x * q.z - q.y * q.w);
		if constexpr (x == 1)
			return T(2) * (q.y * q.z + q.x * q.w);
		else
			return T(1) - T(2) * (q.x * q.x + q.y * q.y);
	}
}
}  // namespace detail

template <typename Functor, typename T>
std::vector<vector3_t<T>> function2D_sampler(const Functor& functor,
                                             T min,
                                             T max,
                                             T step)
{
	std::vector<vector3_t<T>> res;

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

template <typename Functor, typename T>
std::vector<std::vector<vector3_t<T>>> function3D_sampler(
    const Functor& functor,
    vector2_t<T> min,
    vector2_t<T> max,
    vector2_t<T> step)
{
	std::vector<std::vector<vector3_t<T>>> res;

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
		std::vector<vector3_t<T>> row;
		row.reserve(xCount);
		for (T x = min.x; x < max.x; x += step.x)
		{
			row.push_back(functor(x, y));
		}
		res.push_back(std::move(row));
	}

	return res;
}
#pragma endregion

namespace literals
{
#pragma region definition literals
inline cdm::radian_t<float> operator""_rad(long double d)
{
	return cdm::radian_t<float>(static_cast<float>(d));
}
inline cdm::radian_t<float> operator""_rad(unsigned long long int i)
{
	return cdm::radian_t<float>(static_cast<float>(i));
}
inline cdm::radian_t<float> operator""_pi(long double d)
{
	return cdm::radian_t<float>(static_cast<float>(d) * cdm::pi);
}
inline cdm::radian_t<float> operator""_pi(unsigned long long int i)
{
	return cdm::radian_t<float>(static_cast<float>(i) * cdm::pi);
}
inline cdm::degree_t<float> operator""_deg(long double d)
{
	return cdm::degree_t<float>(static_cast<float>(d));
}
inline cdm::degree_t<float> operator""_deg(unsigned long long int i)
{
	return cdm::degree_t<float>(static_cast<float>(i));
}
inline cdm::radian_t<double> operator""_radd(long double d)
{
	return cdm::radian_t<double>(static_cast<double>(d));
}
inline cdm::radian_t<double> operator""_radd(unsigned long long int i)
{
	return cdm::radian_t<double>(static_cast<double>(i));
}
inline cdm::radian_t<double> operator""_pid(long double d)
{
	return cdm::radian_t<double>(static_cast<double>(d) * cdm::pi);
}
inline cdm::radian_t<double> operator""_pid(unsigned long long int i)
{
	return cdm::radian_t<double>(static_cast<double>(i) * cdm::pi);
}
inline cdm::degree_t<double> operator""_degd(long double d)
{
	return cdm::degree_t<double>(static_cast<double>(d));
}
inline cdm::degree_t<double> operator""_degd(unsigned long long int i)
{
	return cdm::degree_t<double>(static_cast<double>(i));
}
#pragma endregion
}  // namespace literals

using namespace literals;

using complex = complex_t<float>;
using complexd = complex_t<double>;
using radian = radian_t<float>;
using radiand = radian_t<double>;
using degree = degree_t<float>;
using degreed = degree_t<double>;
using pi_fraction = pi_fraction_t<int>;
template <int NumT, int DenT>
using static_pi_fraction = static_pi_fraction_t<int, NumT, DenT>;
using vector2 = vector2_t<float>;
using vector2d = vector2_t<double>;
using vector3 = vector3_t<float>;
using vector3d = vector3_t<double>;
using vector4 = vector4_t<float>;
using vector4d = vector4_t<double>;
using matrix2 = matrix2_t<float>;
using matrix2d = matrix2_t<double>;
using matrix3 = matrix3_t<float>;
using matrix3d = matrix3_t<double>;
using matrix4 = matrix4_t<float>;
using matrix4d = matrix4_t<double>;
using perspective = perspective_t<float>;
using perspectived = perspective_t<double>;
using euler_angles = euler_angles_t<float>;
using euler_anglesd = euler_angles_t<double>;
using quaternion = quaternion_t<float>;
using quaterniond = quaternion_t<double>;
using cartesian_direction2d = cartesian_direction2d_t<float>;
using cartesian_direction2dd = cartesian_direction2d_t<double>;
using polar_direction = polar_direction_t<float>;
using polar_directiond = polar_direction_t<double>;
using line = line_t<float>;
using lined = line_t<double>;
using segment2d = segment2d_t<float>;
using segment2dd = segment2d_t<double>;
using plane = plane_t<float>;
using planed = plane_t<double>;
using aa_rect = aa_rect_t<float>;
using aa_rectd = aa_rect_t<double>;
using circle = circle_t<float>;
using circled = circle_t<double>;
using ray2d = ray2d_t<float>;
using ray2dd = ray2d_t<double>;
using ray3d = ray3d_t<float>;
using ray3dd = ray3d_t<double>;
using aabb = aabb_t<float>;
using aabbd = aabb_t<double>;
using transform2d = transform2d_t<float>;
using transform2dd = transform2d_t<double>;
using transform3d = transform3d_t<float>;
using transform3dd = transform3d_t<double>;
using uniform_transform2d = uniform_transform2d_t<float>;
using uniform_transform2dd = uniform_transform2d_t<double>;
using uniform_transform3d = uniform_transform3d_t<float>;
using uniform_transform3dd = uniform_transform3d_t<double>;
using unscaled_transform2d = unscaled_transform2d_t<float>;
using unscaled_transform2dd = unscaled_transform2d_t<double>;
using unscaled_transform3d = unscaled_transform3d_t<float>;
using unscaled_transform3dd = unscaled_transform3d_t<double>;
using direction = direction_t<float>;
using directiond = direction_t<double>;
}  // namespace cdm

#endif  // CDM_MATHS_HPP
