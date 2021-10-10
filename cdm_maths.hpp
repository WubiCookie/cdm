/* cdm_maths - v2.0.0 - geometric library - https://github.com/WubiCookie/cdm
   no warranty implied; use at your own risk

LICENSE

       DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
                   Version 2, December 2004

Copyright (C) 2021 Charles Seizilles de Mazancourt <charles DOT de DOT mazancourt AT hotmail DOT fr>

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
template<typename T>
class normalized;
template<typename T>
struct complex;
template<typename T>
struct radian;
template<typename T>
struct degree;
template<typename T>
struct vector2;
template<typename T>
struct vector3;
template<typename T>
struct vector4;
template<typename T>
struct matrix2;
template<typename T>
struct matrix3;
template<typename T>
struct matrix4;
template<typename T>
struct perspective;
template<typename T>
struct euler_angles;
template<typename T>
struct quaternion;
template<typename T>
struct cartesian_direction2d;
template<typename T>
struct polar_direction;
template<typename T>
struct line;
template<typename T>
struct segment2d;
template<typename T>
struct plane;
//struct rect;
template<typename T>
struct aa_rect;
template<typename T>
struct circle;
//struct ellipse;
//struct aa_ellipse;
template<typename T>
struct ray2d;
template<typename T>
struct ray3d;
template<typename T>
struct aabb;
template<typename T>
struct transform2d;
template<typename T>
struct transform3d;
template<typename T>
struct uniform_transform2d;
template<typename T>
struct uniform_transform3d;
template<typename T>
struct unscaled_transform2d;
template<typename T>
struct unscaled_transform3d;

namespace detail
{
template<uint8_t x, uint8_t y, typename T>
T get_quaternion_matrix_element(quaternion<T> q);
} // namespace detail

template<typename T>
constexpr T lerp(T begin, T end, T percent);
template<typename T>
constexpr T clamp(T f, T min, T max);

template<typename T>
bool nearly_equal(T f1, T f2, T e = T(epsilon));

template<typename T>
constexpr int sign(T val);
template<typename T>
constexpr int signnum(T val);
#pragma endregion

#pragma region constants declarations
constexpr double pi = 3.141592653589793238462643;
constexpr double deg_to_rad = pi / 180.0;
constexpr double rad_to_deg = 180.0 / pi;
constexpr double epsilon = 1.0e-05;
#pragma endregion

#pragma region declaration complex
template<typename T>
struct complex
{
	T r;
	T i;

	static complex from(radian<T> angle);

	complex operator+(complex c) const;
	complex operator-(complex c) const;

	complex& operator+=(complex c);
	complex& operator-=(complex c);
	complex& operator*=(complex c);
};

template<typename T>
complex<T> operator*(complex<T> c1, complex<T> c2);
template<typename T>
complex<T> operator*(complex<T> c, T f);
template<typename T>
complex<T> operator*(normalized<complex<T>> c1, complex<T> c2);
template<typename T>
complex<T> operator*(complex<T> c1, normalized<complex<T>> c2);
template<typename T>
normalized<complex<T>> operator*(normalized<complex<T>> c1, normalized<complex<T>> c2);
template<typename T>
vector2<T> operator*(normalized<complex<T>> c, vector2<T> v);

template<typename T>
T norm(complex<T> c);
template<typename T>
T norm_squared(complex<T> c);
template<typename T>
complex<T> normalize(complex<T> c);
template<typename T>
complex<T> conjugate(complex<T> c);
#pragma endregion

#pragma region declaration radian
template<typename T>
struct radian
{
	T angle;

	radian() = default;
	explicit radian(T f);
	radian(const radian& r) = default;
	radian(radian&& r) = default;
	radian(degree<T> d);

	explicit operator T() const;

	radian& operator+=(T f);
	radian& operator+=(radian r);
	radian& operator-=(T f);
	radian& operator-=(radian r);
	radian& operator*=(T f);
	radian& operator*=(radian r);
	radian& operator/=(T f);
	radian& operator/=(radian r);

	radian operator-() const;

	radian& operator=(const radian& r) = default;
	radian& operator=(degree<T> d);
};

template<typename T>
radian<T> operator+(radian<T>, radian<T>);
template<typename T>
radian<T> operator-(radian<T>, radian<T>);
template<typename T>
radian<T> operator*(radian<T>, radian<T>);
template<typename T>
radian<T> operator/(radian<T>, radian<T>);
template<typename T>
radian<T> operator*(radian<T>, T);
template<typename T>
radian<T> operator/(radian<T>, T);
template<typename T>
radian<T> operator*(T, radian<T>);
template<typename T>
radian<T> operator/(T, radian<T>);
template<typename T>
bool operator<(T, radian<T>);
template<typename T>
bool operator>(T, radian<T>);
template<typename T>
bool operator==(T, radian<T>);
template<typename T>
bool operator!=(T, radian<T>);
template<typename T>
bool operator>=(T, radian<T>);
template<typename T>
bool operator<=(T, radian<T>);
template<typename T>
bool operator<(radian<T>, T);
template<typename T>
bool operator>(radian<T>, T);
template<typename T>
bool operator==(radian<T>, T);
template<typename T>
bool operator!=(radian<T>, T);
template<typename T>
bool operator>=(radian<T>, T);
template<typename T>
bool operator<=(radian<T>, T);
template<typename T>
bool operator<(radian<T>, radian<T>);
template<typename T>
bool operator>(radian<T>, radian<T>);
template<typename T>
bool operator==(radian<T>, radian<T>);
template<typename T>
bool operator!=(radian<T>, radian<T>);
template<typename T>
bool operator>=(radian<T>, radian<T>);
template<typename T>
bool operator<=(radian<T>, radian<T>);

template<typename T>
T sin(radian<T> r);
template<typename T>
T cos(radian<T> r);
template<typename T>
T tan(radian<T> r);
template<typename T>
T asin(radian<T> r);
template<typename T>
T acos(radian<T> r);
template<typename T>
T atan(radian<T> r);
template<typename T>
T sinh(radian<T> r);
template<typename T>
T cosh(radian<T> r);
template<typename T>
T tanh(radian<T> r);
template<typename T>
T asinh(radian<T> r);
template<typename T>
T acosh(radian<T> r);
template<typename T>
T atanh(radian<T> r);
#pragma endregion

#pragma region declaration degree
template<typename T>
struct degree
{
	T angle;

	degree() = default;
	explicit degree(T f);
	degree(const degree& d) = default;
	degree(degree&& d) = default;
	degree(radian<T> r);

	explicit operator T() const;

	degree& operator+=(T f);
	degree& operator+=(degree d);
	degree& operator-=(T f);
	degree& operator-=(degree d);
	degree& operator*=(T f);
	degree& operator*=(degree d);
	degree& operator/=(T f);
	degree& operator/=(degree d);

	degree operator-() const;

	degree& operator=(const degree& d) = default;
	degree& operator=(radian<T> r);
};

template<typename T>
degree<T> operator+(degree<T>, degree<T>);
template<typename T>
degree<T> operator-(degree<T>, degree<T>);
template<typename T>
degree<T> operator*(degree<T>, degree<T>);
template<typename T>
degree<T> operator/(degree<T>, degree<T>);
template<typename T>
degree<T> operator*(degree<T>, T);
template<typename T>
degree<T> operator/(degree<T>, T);
template<typename T>
degree<T> operator*(T, degree<T>);
template<typename T>
degree<T> operator/(T, degree<T>);
template<typename T>
bool operator<(T, degree<T>);
template<typename T>
bool operator>(T, degree<T>);
template<typename T>
bool operator==(T, degree<T>);
template<typename T>
bool operator!=(T, degree<T>);
template<typename T>
bool operator>=(T, degree<T>);
template<typename T>
bool operator<=(T, degree<T>);
template<typename T>
bool operator<(degree<T>, T);
template<typename T>
bool operator>(degree<T>, T);
template<typename T>
bool operator==(degree<T>, T);
template<typename T>
bool operator!=(degree<T>, T);
template<typename T>
bool operator>=(degree<T>, T);
template<typename T>
bool operator<=(degree<T>, T);
template<typename T>
bool operator<(degree<T>, degree<T>);
template<typename T>
bool operator>(degree<T>, degree<T>);
template<typename T>
bool operator==(degree<T>, degree<T>);
template<typename T>
bool operator!=(degree<T>, degree<T>);
template<typename T>
bool operator>=(degree<T>, degree<T>);
template<typename T>
bool operator<=(degree<T>, degree<T>);

template<typename T>
T sin(degree<T> d);
template<typename T>
T cos(degree<T> d);
template<typename T>
T tan(degree<T> d);
template<typename T>
T asin(degree<T> d);
template<typename T>
T acos(degree<T> d);
template<typename T>
T atan(degree<T> d);
template<typename T>
T sinh(degree<T> d);
template<typename T>
T cosh(degree<T> d);
template<typename T>
T tanh(degree<T> d);
template<typename T>
T asinh(degree<T> d);
template<typename T>
T acosh(degree<T> d);
template<typename T>
T atanh(degree<T> d);
#pragma endregion

#pragma region declaration vector2
template<typename T>
struct vector2
{
	T x, y;

	template<typename U = T>
	std::array<U, 2> to_array() const;

	bool lies_on(line<T> l);
	//bool lies_on(const rect& r);
	bool lies_on(aa_rect<T> r);
	bool lies_on(circle<T> c);
	//bool lies_on(const ellipe& e);
	//bool lies_on(const aa_ellipe& e);

	vector2 operator+(vector2 v) const;
	vector2 operator-(vector2 v) const;
	vector2 operator*(T f) const;
	vector2 operator/(T f) const;

	vector2& operator+=(vector2 v);
	vector2& operator-=(vector2 v);
	vector2& operator*=(T f);
	vector2& operator/=(T f);

	vector2 operator-() const;

	bool operator==(vector2 v) const;
	bool operator!=(vector2 v) const;
};

template<typename T>
vector2<T> operator*(T f, vector2<T> v);

template<typename T>
T norm(vector2<T> v);
template<typename T>
T norm_squared(vector2<T> v);
template<typename T>
vector2<T> normalize(vector2<T> v);
template<typename T>
vector2<T> clamp(vector2<T> v, vector2<T> min, vector2<T> max);
template<typename T>
T dot(vector2<T> lhs, vector2<T> rhs);
template<typename T>
T cross(vector2<T> lhs, vector2<T> rhs);
template<typename T>
vector2<T> lerp(vector2<T> begin, vector2<T> end, T percent);
template<typename T>
vector2<T> nlerp(vector2<T> begin, vector2<T> end, T percent);
template<typename T>
vector2<T> slerp(vector2<T> begin, vector2<T> end, T percent);
template<typename T>
T distance_between(vector2<T> v1, vector2<T> v2);
template<typename T>
T distance_squared_between(vector2<T> v1, vector2<T> v2);
template<typename T>
vector2<T> from_to(vector2<T> from, vector2<T> to);
template<typename T>
radian<T> angle_between(vector2<T> v1, vector2<T> v2);
template<typename T>
bool nearly_equal(vector2<T> v1, vector2<T> v2, T e = T(epsilon));
#pragma endregion

#pragma region declaration vector3
template<typename T>
struct vector3
{
	T x, y, z;

	template<typename U = T>
	std::array<U, 3> to_array() const;

	radian<T> angle_around_axis(vector3 v, vector3 axis);

	vector2<T> xy() const;

	vector3 operator+(vector3 v) const;
	vector3 operator-(vector3 v) const;
	vector3 operator*(T f) const;
	vector3 operator/(T f) const;

	vector3& operator+=(vector3 v);
	vector3& operator-=(vector3 v);
	vector3& operator*=(T f);
	vector3& operator/=(T f);

	vector3 operator-() const;

	bool operator==(vector3 v) const;
	bool operator!=(vector3 v) const;
};

template<typename T>
vector3<T> operator*(T f, vector3<T> v);

template<typename T>
T norm(vector3<T> v);
template<typename T>
T norm_squared(vector3<T> v);
template<typename T>
vector3<T> normalize(vector3<T> v);
template<typename T>
vector3<T> clamp(vector3<T> v, vector3<T> min, vector3<T> max);
template<typename T>
T dot(vector3<T> lhs, vector3<T> rhs);
template<typename T>
vector3<T> cross(vector3<T> lhs, vector3<T> rhs);
template<typename T>
vector3<T> lerp(vector3<T> begin, vector3<T> end, T percent);
template<typename T>
vector3<T> nlerp(vector3<T> begin, vector3<T> end, T percent);
template<typename T>
T distance_between(vector3<T> v1, vector3<T> v2);
template<typename T>
T distance_squared_between(vector3<T> v1, vector3<T> v2);
template<typename T>
vector3<T> from_to(vector3<T> from, vector3<T> to);
template<typename T>
radian<T> angle_between(vector3<T> v1, vector3<T> v2);
template<typename T>
bool nearly_equal(vector3<T> v1, vector3<T> v2, T e = epsilon);
#pragma endregion

#pragma region declaration vector4
template<typename T>
struct vector4
{
	T x, y, z, w;

	template<typename U = T>
	std::array<U, 4> to_array() const;

	vector2<T> xy() const;
	vector3<T> xyz() const;

	vector4 operator+(vector4 v) const;
	vector4 operator-(vector4 v) const;
	vector4 operator*(T f) const;
	vector4 operator/(T f) const;

	vector4& operator+=(vector4 v);
	vector4& operator-=(vector4 v);
	vector4& operator*=(T f);
	vector4& operator/=(T f);

	vector4 operator-() const;

	bool operator==(vector4 v) const;
	bool operator!=(vector4 v) const;
};

template<typename T>
T dot(vector4<T> lhs, vector4<T> rhs);

template<typename T>
T norm(vector4<T> v);
template<typename T>
T norm_squared(vector4<T> v);
template<typename T>
vector4<T> normalize(vector4<T> v);
template<typename T>
vector4<T> clamp(vector4<T> v, vector4<T> min, vector4<T> max);
template<typename T>
T dot(vector4<T> lhs, vector4<T> rhs);
// vector4<T> cross(vector4<T> lhs, vector4<T> rhs);
template<typename T>
vector4<T> lerp(vector4<T> begin, vector4<T> end, T percent);
template<typename T>
vector4<T> nlerp(vector4<T> begin, vector4<T> end, T percent);
// vector4<T> slerp(vector4<T> begin, vector4<T> end, T percent);
template<typename T>
T distance_between(vector4<T> v1, vector4<T> v2);
template<typename T>
T distance_squared_between(vector4<T> v1, vector4<T> v2);
template<typename T>
vector4<T> from_to(vector4<T> from, vector4<T> to);
// radian angle_between(vector4<T> lhs, vector4<T> rhs);
template<typename T>
bool nearly_equal(vector4<T> v1, vector4<T> v2, T e = epsilon);
#pragma endregion

#pragma region declaration normalized
template<typename T>
class normalized
{
	T vector;

public:
	normalized() = default;
	normalized(const normalized&) = default;
	normalized(normalized&&) = default;
	normalized(const T& t);
	normalized(T&& t) noexcept;

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
#pragma endregion

#pragma region declaration matrix2
template<typename T>
struct matrix2
{
	T m00, m10,
	  m01, m11;

	template<typename U = T>
	std::array<U, 4> to_array() const;

	static matrix2 zero();
	static matrix2 identity();
	static matrix2 rotation(radian<T> angle);
	static matrix2 rotation(normalized<complex<T>> angle);
	matrix2 transposed() const;
	T determinant() const;

	T& at(uint8_t x, uint8_t y);
	const T& at(uint8_t x, uint8_t y) const;

	vector2<T> operator*(vector2<T> v) const;
	matrix2 operator*(matrix2 m) const;
};

template<typename T>
matrix2<T> transpose(matrix2<T> m);
#pragma endregion

#pragma region declaration matrix3
template<typename T>
struct matrix3
{
	T m00, m10, m20,
	  m01, m11, m21,
	  m02, m12, m22;

	template<typename U = T>
	std::array<U, 9> to_array() const;

	static matrix3 zero();
	static matrix3 identity();
	static matrix3 rotation(euler_angles<T> r);
	static matrix3 rotation(quaternion<T> q);
	static matrix3 rotation_around_x(radian<T> angle);
	static matrix3 rotation_around_y(radian<T> angle);
	static matrix3 rotation_around_z(radian<T> angle);
	static matrix3 rotation_around_x(normalized<complex<T>> angle);
	static matrix3 rotation_around_y(normalized<complex<T>> angle);
	static matrix3 rotation_around_z(normalized<complex<T>> angle);
	matrix3 inversed() const;
	matrix3 transposed() const;
	T determinant() const;
	bool is_orthogonal() const;

	T& at(uint8_t x, uint8_t y);
	const T& at(uint8_t x, uint8_t y) const;

	matrix3 operator*(T f) const;
	vector3<T> operator*(vector3<T> v) const;
	matrix3 operator*(const matrix3& m) const;
};

template<typename T>
matrix3<T> inverse(matrix3<T> m);
template<typename T>
matrix3<T> transpose(matrix3<T> m);
#pragma endregion

#pragma region declaration matrix4
template<typename T>
struct matrix4
{
	T m00, m10, m20, m30,
	  m01, m11, m21, m31,
	  m02, m12, m22, m32,
	  m03, m13, m23, m33;

	template<typename U = T>
	std::array<U, 16> to_array() const;

	static matrix4 zero();
	static matrix4 identity();
	static matrix4 rotation(euler_angles<T> r);
	static matrix4 rotation(quaternion<T> q);
	static matrix4 translation(vector3<T> t);
	static matrix4 translation(T x, T y, T z);
	static matrix4 scale(vector3<T> t);
	static matrix4 look_at(vector3<T> from, vector3<T> to, vector3<T> up = { T(0), T(1), T(0) });
	static matrix4 orthographic(T left, T right, T top, T bottom, T near, T far);
	static matrix4 rotation_around_x(radian<T> angle);
	static matrix4 rotation_around_y(radian<T> angle);
	static matrix4 rotation_around_z(radian<T> angle);
	static matrix4 rotation_around_x(normalized<complex<T>> angle);
	static matrix4 rotation_around_y(normalized<complex<T>> angle);
	static matrix4 rotation_around_z(normalized<complex<T>> angle);
	bool is_orthogonal() const;
	bool is_homogenous() const;
	matrix4 inversed() const;
	matrix4 transposed() const;
	T determinant() const;

	T& at(uint8_t x, uint8_t y);
	const T& at(uint8_t x, uint8_t y) const;

	matrix4 operator*(T f) const;
	matrix4 operator/(T f) const;
	vector4<T> operator*(vector4<T> v) const;
	matrix4 operator*(const matrix4& m) const;
};
#pragma endregion

#pragma region declaration perspective
template<typename T>
struct perspective
{
private:
	radian<T> m_angle;
	T m_ratio;
	T m_inv_ratio;
	T m_near;
	T m_far;
	T m_invTanHalfFovy;

public:
	void set(radian<T> angle, T ratio, T near, T far);

	void set_angle(radian<T> angle);
	radian<T> get_angle() const;

	void set_ratio(T ratio);
	T get_ratio() const;

	void set_near(T near_plane);
	T get_near() const;

	void set_far(T far_plane);
	T get_far() const;

	matrix4<T> to_matrix4() const;

	friend matrix4<T> operator*(const matrix4<T>& m, const perspective& p);
	friend matrix4<T> operator*(const perspective& p, const matrix4<T>& m);

	friend matrix4<T> operator*(const unscaled_transform3d<T>& t, const perspective& p);
	friend matrix4<T> operator*(const perspective& p, const unscaled_transform3d<T>& t);
};

template<typename T>
matrix4<T> operator*(const matrix4<T>& m, const perspective<T>& p);
template<typename T>
matrix4<T> operator*(const perspective<T>& p, const matrix4<T>& m);

template<typename T>
matrix4<T> operator*(const unscaled_transform3d<T>& t, const perspective<T>& p);
template<typename T>
matrix4<T> operator*(const perspective<T>& p, const unscaled_transform3d<T>& t);
#pragma endregion

#pragma region declaration euler_angles
template<typename T>
struct euler_angles
{
	radian<T> x, y, z;
};
#pragma endregion

#pragma region declaration quaternion
template<typename T>
struct quaternion
{
	T x, y, z, w;

	template<typename U = T>
	std::array<U, 4> to_array() const;

	static quaternion zero();
	static quaternion identity();
	static quaternion from(const normalized<vector3<T>>& axis, radian<T> angle);

	quaternion inversed() const;
	quaternion conjugated() const;
	quaternion normalized() const;
	quaternion clamped(quaternion min, quaternion max) const;

	quaternion operator+(quaternion q) const;
	quaternion operator-(quaternion q) const;
	quaternion operator*(quaternion q) const; // TODO: implement operator* for normalized<> if result is normalized too
	quaternion operator*(T f) const;
	quaternion operator/(T f) const;

	quaternion& operator+=(quaternion q);
	quaternion& operator-=(quaternion q);
	quaternion& operator*=(quaternion q);
	quaternion& operator*=(T f);
	quaternion& operator/=(T f);

	quaternion operator-() const;

	bool operator==(quaternion v) const;
	bool operator!=(quaternion v) const;
};

template<typename T>
vector3<T> operator*(const normalized<quaternion<T>>& q, vector3<T> v);

template<typename T>
T norm(quaternion<T> q);
template<typename T>
T norm_squared(quaternion<T> q);
template<typename T>
quaternion<T>& normalize(quaternion<T>& q);
template<typename T>
quaternion<T>& clamp(quaternion<T>& q, quaternion<T> min, quaternion<T> max);
template<typename T>
quaternion<T>& negate(quaternion<T>& q);
template<typename T>
quaternion<T>& inverse(quaternion<T>& q);
template<typename T>
quaternion<T>& conjugate(quaternion<T>& q);
template<typename T>
T dot(quaternion<T> lhs, quaternion<T> rhs);
template<typename T>
quaternion<T> cross(quaternion<T> lhs, quaternion<T> rhs);
template<typename T>
quaternion<T> lerp(quaternion<T> begin, quaternion<T> end, T percent);
template<typename T>
quaternion<T> nlerp(quaternion<T> begin, quaternion<T> end, T percent);
template<typename T>
quaternion<T> slerp(quaternion<T> begin, quaternion<T> end, T percent);
#pragma endregion

#pragma region declaration cartesian_direction2d
template<typename T>
struct cartesian_direction2d
{
	T x, y;

	static cartesian_direction2d from(const normalized<vector2<T>>& direction);
	static cartesian_direction2d from(polar_direction<T> direction);
	static cartesian_direction2d from(radian<T> angle);

	cartesian_direction2d& rotate(radian<T> angle);
};
#pragma endregion

#pragma region declaration polar_direction
template<typename T>
struct polar_direction
{
	radian<T> angle;

	static polar_direction from(const normalized<vector2<T>>& direction);
	static polar_direction from(cartesian_direction2d<T> direction);
	static polar_direction from(radian<T> angle);

	polar_direction& rotate(radian<T> angle);

	polar_direction& wrap();
};
#pragma endregion

#pragma region declaration line
template<typename T>
struct line
{
	T coefficient;
	T offset;

	static line from(vector2<T> v);
};

template<typename T>
bool are_parallel(line<T> l1, line<T> l2);

template<typename T>
bool collides(line<T> l1, line<T> l2);
template<typename T>
bool collides(line<T> l1, line<T> l2, vector2<T>& intersection);
#pragma endregion

#pragma region declaration segment2d
template<typename T>
struct segment2d
{
	vector2<T> origin;
	vector2<T> end;
};

template<typename T>
bool collides(const segment2d<T>& s0,
              const segment2d<T>& s1,
              vector2<T>& outPoint,
              T e = T(epsilon));
#pragma endregion

#pragma region declaration plane
template<typename T>
struct plane
{
	vector3<T> origin;
	normalized<vector3<T>> normal;

	T evaluate(vector3<T> point) const;
	vector3<T> project3d(vector3<T> point) const;
	vector2<T> project2d(vector3<T> point, normalized<vector3<T>> plane_tangent) const;
	vector3<T> unproject(vector2<T> point, normalized<vector3<T>> plane_tangent) const;

	normalized<vector3<T>> computeTangent(T e = T(epsilon)) const;
};
#pragma endregion

#pragma region declaration ray2d
template<typename T>
struct ray2d
{
	vector2<T> origin;
	normalized<vector2<T>> direction;
};
#pragma endregion

#pragma region declaration ray3d
template<typename T>
struct ray3d
{
	vector3<T> origin;
	normalized<vector3<T>> direction;
};
#pragma endregion

//struct rect
//{
//	vector2 origin;
//	vector2 dimention;
//	complex angle;
//};
//
//bool collides(const rect& r1, const rect& r2);
//std::size_t collides(const rect& r1, const rect& r2, vector2* intersection1, vector2* insersection2);
//bool collides(const rect& r, const line& l);
//std::size_t collides(const rect& r, const line& l, vector2* intersection1, vector2* insersection2);
//bool collides(const line& l, const rect& r);
//std::size_t collides(const line& l, const rect& r, vector2* intersection1, vector2* insersection2);
//bool collides(const rect& r, vector2 v);
//bool collides(vector2 v, const rect& r);

#pragma region declaration aa_rect
template<typename T>
struct aa_rect
{
	vector2<T> origin;
	vector2<T> dimention;

	bool contains(vector2<T> v) const;
};

template<typename T>
std::size_t collides(aa_rect<T> r1, aa_rect<T> r2, vector2<T>* intersection1, vector2<T>* insersection2);
//std::size_t collides(const aa_rect& r1, const rect& r2, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const rect& r1, const aa_rect& r2, vector2* intersection1, vector2* insersection2);
template<typename T>
std::size_t collides(aa_rect<T> r, line<T> l, vector2<T>* intersection1, vector2<T>* insersection2);
template<typename T>
std::size_t collides(line<T> l, aa_rect<T> r, vector2<T>* intersection1, vector2<T>* insersection2);
#pragma endregion

#pragma region declaration circle
template<typename T>
struct circle
{
	vector2<T> origin;
	T radius;
};

template<typename T>
std::size_t collides(circle<T> c1, circle<T> c2, vector2<T>* intersection1, vector2<T>* insersection2);
template<typename T>
std::size_t collides(circle<T> c, aa_rect<T> r, vector2<T>* intersection1, vector2<T>* insersection2);
template<typename T>
std::size_t collides(aa_rect<T> r, circle<T> c, vector2<T>* intersection1, vector2<T>* insersection2);
//std::size_t collides(const circle& c, const rect& r, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const rect& r, const circle& c, vector2* intersection1, vector2* insersection2);
template<typename T>
std::size_t collides(circle<T> c, line<T> l, vector2<T>* intersection1, vector2<T>* insersection2);
template<typename T>
std::size_t collides(line<T> l, circle<T> c, vector2<T>* intersection1, vector2<T>* insersection2);

//struct ellipse
//{
//	vector2 origin;
//	T semi_minor = T(1);
//	T semi_major = T(1);
//	complex angle;
//};
//
//std::size_t collides(const ellipse& e1, const ellipse& e2, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const ellipse& e, const circle& c, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const circle& c, const ellipse& e, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const ellipse& e, const aa_rect& r, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const aa_rect& r, const ellipse& e, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const ellipse& e, const rect& r, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const rect& r, const ellipse& e, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const ellipse& e, const line& l, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const line& l, const ellipse& e, vector2* intersection1, vector2* insersection2);
//bool collides(const ellipse& e, vector2 v);
//bool collides(vector2 v, const ellipse& e);
//
//struct aa_ellipse
//{
//	vector2 origin;
//	T semi_minor = T(1);
//	T semi_major = T(1);
//};
//
//std::size_t collides(const aa_ellipse& e1, const aa_ellipse& e2, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const aa_ellipse& e1, const ellipse& e2, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const ellipse& e1, const aa_ellipse& e2, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const aa_ellipse& e, const circle& c, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const circle& c, const aa_ellipse& e, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const aa_ellipse& e, const aa_rect& r, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const aa_rect& r, const aa_ellipse& e, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const aa_ellipse& e, const rect& r, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const rect& r, const aa_ellipse& e, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const aa_ellipse& e, const line& l, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const line& l, const aa_ellipse& e, vector2* intersection1, vector2* insersection2);
//bool collides(const aa_ellipse& e, vector2 v);
//bool collides(vector2 v, const aa_ellipse& e);

template<typename T>
bool collides(ray3d<T> r, plane<T> p);
template<typename T>
bool collides(ray3d<T> r, plane<T> p, vector3<T>& intersection);
template<typename T>
bool collides(plane<T> p, ray3d<T> r);
template<typename T>
bool collides(plane<T> p, ray3d<T> r, vector3<T>& intersection);
template<typename T>
bool collides_bidirectional(plane<T> p, ray3d<T> r, vector3<T>& intersection);
#pragma endregion

#pragma region declaration aabb
template<typename T>
struct aabb
{
	vector3<T> min;
	vector3<T> max;

	bool contains(vector3<T> p) const;
	vector3<T> get_center() const;

	aabb operator+(aabb rhs) const;
};

template<typename T>
bool collides(aabb<T> b, ray3d<T> r);
template<typename T>
bool collides(aabb<T> b1, aabb<T> b2);
#pragma endregion

#pragma region declaration transform2d
template<typename T>
struct transform2d
{
	vector2<T> position;
	normalized<complex<T>> rotation;
	vector2<T> scale;

	transform2d operator*(const transform2d& t) const;
	vector2<T> operator*(vector2<T> v) const;
};
#pragma endregion

#pragma region declaration transform3d
template<typename T>
struct transform3d
{
	vector3<T> position;
	quaternion<T> rotation;
	vector3<T> scale;

	transform3d operator*(const transform3d& t) const;
	vector3<T> operator*(vector3<T> v) const;
	quaternion<T> operator*(quaternion<T> q) const;
};
#pragma endregion

#pragma region declaration uniform_transform2d
template<typename T>
struct uniform_transform2d
{
	vector2<T> position;
	radian<T> rotation;
	T scale;

	uniform_transform2d operator*(uniform_transform2d t) const;
	vector2<T> operator*(vector2<T> v) const;
};
#pragma endregion

#pragma region declaration uniform_transform3d
template<typename T>
struct uniform_transform3d
{
	vector3<T> position;
	quaternion<T> rotation;
	T scale;

	uniform_transform3d operator*(const uniform_transform3d& t) const;
	vector3<T> operator*(vector3<T> v) const;
	quaternion<T> operator*(quaternion<T> q) const;
};
#pragma endregion

#pragma region declaration unscaled_transform2d
template<typename T>
struct unscaled_transform2d
{
	vector2<T> position;
	radian<T> rotation;

	unscaled_transform2d operator*(unscaled_transform2d t) const;
	vector2<T> operator*(vector2<T> v) const;
};
#pragma endregion

#pragma region declaration unscaled_transform3d
template<typename T>
struct unscaled_transform3d
{
	vector3<T> position;
	quaternion<T> rotation;

	unscaled_transform3d& translate_absolute(vector3<T> t);
	unscaled_transform3d& translate_relative(vector3<T> t);
	unscaled_transform3d& rotate(quaternion<T> r);

	matrix4<T> to_matrix() const;

	unscaled_transform3d operator*(const unscaled_transform3d& t) const;
	vector3<T> operator*(vector3<T> v) const;
	quaternion<T> operator*(quaternion<T> q) const;
};
#pragma endregion

#pragma region definition misc
template<typename T>
unscaled_transform3d<T> inverse(unscaled_transform3d<T> tr);

template<typename T>
constexpr T lerp(T begin, T end, T percent)
{
	return (end - begin) * percent + begin;
}

template<typename T>
constexpr T clamp(T f, T min, T max)
{
	return std::min(std::max(f, min), max);
}

template<typename T>
bool nearly_equal(T f1, T f2, T e)
{
	return std::abs(f1 - f2) < e;
}

template<typename T>
constexpr int sign(T val)
{
    return (T(0) <= val) - (val < T(0));
}
template<typename T>
constexpr int signnum(T val)
{
    return (T(0) < val) - (val < T(0));
}
#pragma endregion

#pragma region definition complex
template<typename T>
complex<T> complex<T>::from(radian<T> angle)
{
	complex<T> c;
	c.r = cos(angle);
	c.i = sin(angle);
	return c;
}

template<typename T>
complex<T> complex<T>::operator+(complex<T> c) const
{
	return {
		r + c.r,
		i + c.i
	};
}
template<typename T>
complex<T> complex<T>::operator-(complex<T> c) const
{
	return {
		r - c.r,
		i - c.i
	};
}

template<typename T>
complex<T>& complex<T>::operator+=(complex<T> c) { return *this = *this + c; }
template<typename T>
complex<T>& complex<T>::operator-=(complex<T> c) { return *this = *this - c; }
template<typename T>
complex<T>& complex<T>::operator*=(complex<T> c) { return *this = *this * c; }

template<typename T>
complex<T> operator*(complex<T> c1, complex<T> c2)
{
	return {
		c1.r * c2.r - c1.i * c2.i,
		c1.r * c2.i + c1.i * c2.r
	};
}
template<typename T>
complex<T> operator*(complex<T> c, T f) { return {c.r * f, c.i * f}; }
template<typename T>
complex<T> operator*(normalized<complex<T>> c1, complex<T> c2) { return *c1 * c2; }
template<typename T>
complex<T> operator*(complex<T> c1, normalized<complex<T>> c2) { return c1 * *c2; }
template<typename T>
normalized<complex<T>> operator*(normalized<complex<T>> c1, normalized<complex<T>> c2)
{
	return normalized<complex<T>>::already_normalized(complex<T>{
		c1->r * c2->r - c1->i * c2->i,
		c1->r * c2->i + c1->i * c2->r
	});
}
template<typename T>
vector2<T> operator*(normalized<complex<T>> c, vector2<T> v)
{
	return {
		c->r * v.x - c->i * v.y,
		c->r * v.y + c->i * v.x
	};
}

template<typename T>
T norm(complex<T> c) { return std::sqrt(norm_squared(c)); }
template<typename T>
T norm_squared(complex<T> c) { return c.r * c.r + c.i * c.i; }
template<typename T>
complex<T> normalize(complex<T> c)
{
	T n = norm(c);
	c.r /= n;
	c.i /= n;
	return c;
}
template<typename T>
complex<T> conjugate(complex<T> c) { c.i = -c.i; return c; }
#pragma endregion

#pragma region definition radian
template<typename T>
radian<T>::radian(T f) : angle(f) {}
template<typename T>
radian<T>::radian(degree<T> d) : angle(d.angle * T(deg_to_rad)) {}

template<typename T>
radian<T>::operator T() const { return angle; }

template<typename T>
radian<T>& radian<T>::operator+=(T f) { angle += f; return *this; }
template<typename T>
radian<T>& radian<T>::operator+=(radian<T> r) { angle += r.angle; return *this; }
template<typename T>
radian<T>& radian<T>::operator-=(T f) { angle -= f; return *this; }
template<typename T>
radian<T>& radian<T>::operator-=(radian<T> r) { angle -= r.angle; return *this; }
template<typename T>
radian<T>& radian<T>::operator*=(T f) { angle *= f; return *this; }
template<typename T>
radian<T>& radian<T>::operator*=(radian<T> r) { angle *= r.angle; return *this; }
template<typename T>
radian<T>& radian<T>::operator/=(T f) { angle /= f; return *this; }
template<typename T>
radian<T>& radian<T>::operator/=(radian<T> r) { angle /= r.angle; return *this; }

template<typename T>
radian<T> radian<T>::operator-() const { return radian<T>{ -angle }; }

template<typename T>
radian<T>& radian<T>::operator=(degree<T> d) { angle = d.angle * T(deg_to_rad); return *this; }

template<typename T>
radian<T> operator+(radian<T> r1, radian<T> r2) { return radian<T>{ r1.angle + r2.angle }; }
template<typename T>
radian<T> operator-(radian<T> r1, radian<T> r2) { return radian<T>{ r1.angle - r2.angle }; }
template<typename T>
radian<T> operator*(radian<T> r1, radian<T> r2) { return radian<T>{ r1.angle * r2.angle }; }
template<typename T>
radian<T> operator/(radian<T> r1, radian<T> r2) { return radian<T>{ r1.angle / r2.angle }; }
template<typename T>
radian<T> operator*(radian<T> r, T f) { return radian<T>{ r.angle * f }; }
template<typename T>
radian<T> operator/(radian<T> r, T f) { return radian<T>{ r.angle / f }; }
template<typename T>
radian<T> operator*(T f, radian<T> r) { return radian<T>{ f * r.angle }; }
template<typename T>
radian<T> operator/(T f, radian<T> r) { return radian<T>{ f / r.angle }; }
template<typename T>
bool operator<(T lhs, radian<T> rhs) { return T(lhs) < T(rhs); }
template<typename T>
bool operator>(T lhs, radian<T> rhs) { return T(lhs) > T(rhs); }
template<typename T>
bool operator==(T lhs, radian<T> rhs) { return T(lhs) == T(rhs); }
template<typename T>
bool operator!=(T lhs, radian<T> rhs) { return T(lhs) != T(rhs); }
template<typename T>
bool operator>=(T lhs, radian<T> rhs) { return T(lhs) >= T(rhs); }
template<typename T>
bool operator<=(T lhs, radian<T> rhs) { return T(lhs) <= T(rhs); }
template<typename T>
bool operator<(radian<T> lhs, T rhs) { return T(lhs) < T(rhs); }
template<typename T>
bool operator>(radian<T> lhs, T rhs) { return T(lhs) > T(rhs); }
template<typename T>
bool operator==(radian<T> lhs, T rhs) { return T(lhs) == T(rhs); }
template<typename T>
bool operator!=(radian<T> lhs, T rhs) { return T(lhs) != T(rhs); }
template<typename T>
bool operator>=(radian<T> lhs, T rhs) { return T(lhs) >= T(rhs); }
template<typename T>
bool operator<=(radian<T> lhs, T rhs) { return T(lhs) <= T(rhs); }
template<typename T>
bool operator<(radian<T> lhs, radian<T> rhs) { return T(lhs) < T(rhs); }
template<typename T>
bool operator>(radian<T> lhs, radian<T> rhs) { return T(lhs) > T(rhs); }
template<typename T>
bool operator==(radian<T> lhs, radian<T> rhs) { return T(lhs) == T(rhs); }
template<typename T>
bool operator!=(radian<T> lhs, radian<T> rhs) { return T(lhs) != T(rhs); }
template<typename T>
bool operator>=(radian<T> lhs, radian<T> rhs) { return T(lhs) >= T(rhs); }
template<typename T>
bool operator<=(radian<T> lhs, radian<T> rhs) { return T(lhs) <= T(rhs); }

template<typename T>
T sin(radian<T> r) { return std::sin(r.angle); }
template<typename T>
T cos(radian<T> r) { return std::cos(r.angle); }
template<typename T>
T tan(radian<T> r) { return std::tan(r.angle); }
template<typename T>
T asin(radian<T> r) { return std::asin(r.angle); }
template<typename T>
T acos(radian<T> r) { return std::acos(r.angle); }
template<typename T>
T atan(radian<T> r) { return std::atan(r.angle); }
template<typename T>
T sinh(radian<T> r) { return std::sinh(r.angle); }
template<typename T>
T cosh(radian<T> r) { return std::cosh(r.angle); }
template<typename T>
T tanh(radian<T> r) { return std::tanh(r.angle); }
template<typename T>
T asinh(radian<T> r) { return std::asinh(r.angle); }
template<typename T>
T acosh(radian<T> r) { return std::acosh(r.angle); }
template<typename T>
T atanh(radian<T> r) { return std::atanh(r.angle); }
#pragma endregion

#pragma region definition degree
template<typename T>
degree<T>::degree(T f) : angle(f) {}
template<typename T>
degree<T>::degree(radian<T> r) : angle(r.angle * T(rad_to_deg)) {}

template<typename T>
degree<T>::operator T() const { return angle; }

template<typename T>
degree<T>& degree<T>::operator+=(T f) { angle += f; return *this; }
template<typename T>
degree<T>& degree<T>::operator+=(degree d) { angle += d.angle; return *this; }
template<typename T>
degree<T>& degree<T>::operator-=(T f) { angle -= f; return *this; }
template<typename T>
degree<T>& degree<T>::operator-=(degree d) { angle -= d.angle; return *this; }
template<typename T>
degree<T>& degree<T>::operator*=(T f) { angle *= f; return *this; }
template<typename T>
degree<T>& degree<T>::operator*=(degree d) { angle *= d.angle; return *this; }
template<typename T>
degree<T>& degree<T>::operator/=(T f) { angle /= f; return *this; }
template<typename T>
degree<T>& degree<T>::operator/=(degree d) { angle /= d.angle; return *this; }

template<typename T>
degree<T> degree<T>::operator-() const { return radian(-angle); }

template<typename T>
degree<T> operator+(degree<T> d1, degree<T> d2) { return degree<T>{ d1.angle + d2.angle }; }
template<typename T>
degree<T> operator-(degree<T> d1, degree<T> d2) { return degree<T>{ d1.angle - d2.angle }; }
template<typename T>
degree<T> operator*(degree<T> d1, degree<T> d2) { return degree<T>{ d1.angle * d2.angle }; }
template<typename T>
degree<T> operator/(degree<T> d1, degree<T> d2) { return degree<T>{ d1.angle / d2.angle }; }
template<typename T>
degree<T> operator*(degree<T> d, T f) { return degree<T>{ d.angle * f }; }
template<typename T>
degree<T> operator/(degree<T> d, T f) { return degree<T>{ d.angle / f }; }
template<typename T>
degree<T> operator*(T f, degree<T> d) { return degree<T>{ f * d.angle }; }
template<typename T>
degree<T> operator/(T f, degree<T> d) { return degree<T>{ f / d.angle }; }
template<typename T>
bool operator<(T lhs, degree<T> rhs) { return T(lhs) < T(rhs); }
template<typename T>
bool operator>(T lhs, degree<T> rhs) { return T(lhs) > T(rhs); }
template<typename T>
bool operator==(T lhs, degree<T> rhs) { return T(lhs) == T(rhs); }
template<typename T>
bool operator!=(T lhs, degree<T> rhs) { return T(lhs) != T(rhs); }
template<typename T>
bool operator>=(T lhs, degree<T> rhs) { return T(lhs) >= T(rhs); }
template<typename T>
bool operator<=(T lhs, degree<T> rhs) { return T(lhs) <= T(rhs); }
template<typename T>
bool operator<(degree<T> lhs, T rhs) { return T(lhs) < T(rhs); }
template<typename T>
bool operator>(degree<T> lhs, T rhs) { return T(lhs) > T(rhs); }
template<typename T>
bool operator==(degree<T> lhs, T rhs) { return T(lhs) == T(rhs); }
template<typename T>
bool operator!=(degree<T> lhs, T rhs) { return T(lhs) != T(rhs); }
template<typename T>
bool operator>=(degree<T> lhs, T rhs) { return T(lhs) >= T(rhs); }
template<typename T>
bool operator<=(degree<T> lhs, T rhs) { return T(lhs) <= T(rhs); }
template<typename T>
bool operator<(degree<T> lhs, degree<T> rhs) { return T(lhs) < T(rhs); }
template<typename T>
bool operator>(degree<T> lhs, degree<T> rhs) { return T(lhs) > T(rhs); }
template<typename T>
bool operator==(degree<T> lhs, degree<T> rhs) { return T(lhs) == T(rhs); }
template<typename T>
bool operator!=(degree<T> lhs, degree<T> rhs) { return T(lhs) != T(rhs); }
template<typename T>
bool operator>=(degree<T> lhs, degree<T> rhs) { return T(lhs) >= T(rhs); }
template<typename T>
bool operator<=(degree<T> lhs, degree<T> rhs) { return T(lhs) <= T(rhs); }

template<typename T>
T sin(degree<T> d) { return sin(radian<T>::from(d)); }
template<typename T>
T cos(degree<T> d) { return cos(radian<T>::from(d)); }
template<typename T>
T tan(degree<T> d) { return tan(radian<T>::from(d)); }
template<typename T>
T asin(degree<T> d) { return asin(radian<T>::from(d)); }
template<typename T>
T acos(degree<T> d) { return acos(radian<T>::from(d)); }
template<typename T>
T atan(degree<T> d) { return atan(radian<T>::from(d)); }
template<typename T>
T sinh(degree<T> d) { return sinh(radian<T>::from(d)); }
template<typename T>
T cosh(degree<T> d) { return cosh(radian<T>::from(d)); }
template<typename T>
T tanh(degree<T> d) { return tanh(radian<T>::from(d)); }
template<typename T>
T asinh(degree<T> d) { return asinh(radian<T>::from(d)); }
template<typename T>
T acosh(degree<T> d) { return acosh(radian<T>::from(d)); }
template<typename T>
T atanh(degree<T> d) { return atanh(radian<T>::from(d)); }
#pragma endregion

#pragma region definition vector2
template<typename T>
bool vector2<T>::lies_on(line<T> l) { return cdm::nearly_equal(l.coefficient * x + l.offset, y); }
//bool vector2::lies_on(const rect& r)
template<typename T>
bool vector2<T>::lies_on(aa_rect<T> r)
{
	return ((x >= r.origin.x) && (x <= r.origin.x + r.dimention.x) && (lies_on(line{ 0, r.origin.y }) || lies_on(line{ 0, r.origin.y + r.dimention.y }))) ||
		((y >= r.origin.y) && (y <= r.origin.y + r.dimention.y) && (lies_on(line{ r.origin.x, 0 }) || lies_on(line{ r.origin.x + r.dimention.x, 0 })));
}
//bool vector2::lies_on(const circle& c)
//bool vector2::lies_on(const ellipe& e)
//bool vector2::lies_on(const aa_ellipe& e)

template<typename T>
vector2<T> vector2<T>::operator+(vector2<T> v) const { return {x + v.x, y + v.y}; }
template<typename T>
vector2<T> vector2<T>::operator-(vector2<T> v) const { return {x - v.x, y - v.y}; }
template<typename T>
vector2<T> vector2<T>::operator*(T f) const { return {x * f, y * f}; }
template<typename T>
vector2<T> vector2<T>::operator/(T f) const { return {x / f, y / f}; }

template<typename T>
vector2<T>& vector2<T>::operator+=(vector2<T> v) { *this = *this + v; return *this; }
template<typename T>
vector2<T>& vector2<T>::operator-=(vector2<T> v) { *this = *this - v; return *this; }
template<typename T>
vector2<T>& vector2<T>::operator*=(T f) { *this = *this * f; return *this; }
template<typename T>
vector2<T>& vector2<T>::operator/=(T f) { *this = *this / f; return *this; }

template<typename T>
vector2<T> vector2<T>::operator-() const { return {-x, -y}; }

template<typename T>
bool vector2<T>::operator==(vector2<T> v) const { return x == v.x && y == v.y; }
template<typename T>
bool vector2<T>::operator!=(vector2<T> v) const { return !operator==(v); }

template<typename T>
vector2<T> operator*(T f, vector2<T> v) { return v * f; }

template<typename T>
T norm(vector2<T> v) { return std::sqrt(norm_squared(v)); }
template<typename T>
T norm_squared(vector2<T> v) { return v.x * v.x + v.y * v.y; }
template<typename T>
vector2<T>& normalize(vector2<T>& v)
{
	T n = norm(v);
	v.x /= n;
	v.y /= n;
	return v;
}
template<typename T>
vector2<T> clamp(vector2<T> v, vector2<T> min, vector2<T> max)
{
	v.x = cdm::clamp(v.x, min.x, max.x);
	v.y = cdm::clamp(v.y, min.y, max.y);
	return v;
}
template<typename T>
T dot(vector2<T> lhs, vector2<T> rhs) { return lhs.x * rhs.x + lhs.y * rhs.y; }
template<typename T>
T cross(vector2<T> lhs, vector2<T> rhs) { return lhs.x * rhs.y - lhs.y * rhs.x; }
template<typename T>
vector2<T> lerp(vector2<T> begin, vector2<T> end, T percent) { return (end - begin) * percent + begin; }
template<typename T>
vector2<T> nlerp(vector2<T> begin, vector2<T> end, T percent) { return lerp(begin, end, percent).get_normalized(); }
template<typename T>
vector2<T> slerp(vector2<T> begin, vector2<T> end, T percent)
{
	const radian angle = angle_between(begin, end) * percent;
	const T s = sin(angle);
	const T c = cos(angle);

	normalized<vector2<T>> res {
		{
			c * begin.x - s * begin.y,
			s * begin.x + c * begin.y
		}
	};

	T f = cdm::lerp(norm(begin), norm(end), percent);
	return res * f;
}
template<typename T>
T distance_between(vector2<T> v1, vector2<T> v2) { return norm(v1 - v2); }
template<typename T>
T distance_squared_between(vector2<T> v1, vector2<T> v2) { return norm_squared(v1 - v2); }
template<typename T>
vector2<T> from_to(vector2<T> from, vector2<T> to) { return {to.x - from.x, to.y - from.y}; }
template<typename T>
radian<T> angle_between(vector2<T> v1, vector2<T> v2) { return radian{ atan2f(v2.y, v2.x) - atan2f(v1.y, v1.x) }; }
template<typename T>
bool nearly_equal(vector2<T> v1, vector2<T> v2, T e) { return cdm::nearly_equal(v1.x, v2.x, e) && cdm::nearly_equal(v1.y, v2.y, e); }
#pragma endregion

#pragma region definition vector3
template<typename T>
radian<T> vector3<T>::angle_around_axis(vector3<T> v, vector3<T> axis)
{
	vector3<T> c = cross(*this, v);
	T angle = atan2f(norm(c), dot(*this, v));
	return radian<T>{ dot(c, axis) < T(0) ? -angle : angle };
}

template<typename T>
vector2<T> vector3<T>::xy() const { return {x, y}; }

template<typename T>
vector3<T> vector3<T>::operator+(vector3<T> v) const { return {x + v.x, y + v.y, z + v.z}; }
template<typename T>
vector3<T> vector3<T>::operator-(vector3<T> v) const { return {x - v.x, y - v.y, z - v.z}; }
template<typename T>
vector3<T> vector3<T>::operator*(T f) const { return {x * f, y * f, z * f}; }
template<typename T>
vector3<T> vector3<T>::operator/(T f) const { return {x / f, y / f, z / f}; }

template<typename T>
vector3<T>& vector3<T>::operator+=(vector3<T> v) { *this = *this + v; return *this; }
template<typename T>
vector3<T>& vector3<T>::operator-=(vector3<T> v) { *this = *this - v; return *this; }
template<typename T>
vector3<T>& vector3<T>::operator*=(T f) { *this = *this * f; return *this; }
template<typename T>
vector3<T>& vector3<T>::operator/=(T f) { *this = *this / f; return *this; }

template<typename T>
vector3<T> vector3<T>::operator-() const { return {-x, -y, -z}; }

template<typename T>
bool vector3<T>::operator==(vector3<T> v) const { return x == v.x && y == v.y && z == v.z; }
template<typename T>
bool vector3<T>::operator!=(vector3<T> v) const { return !operator==(v); }

template<typename T>
vector3<T> operator*(T f, vector3<T> v) { return v * f; }

template<typename T>
T norm(vector3<T> v) { return std::sqrt(norm_squared(v)); }
template<typename T>
T norm_squared(vector3<T> v) { return v.x * v.x + v.y * v.y + v.z * v.z; }
template<typename T>
vector3<T> normalize(vector3<T> v)
{
	T n = norm(v);
	v.x /= n;
	v.y /= n;
	v.z /= n;
	return v;
}
template<typename T>
vector3<T> clamp(vector3<T> v, vector3<T> min, vector3<T> max)
{
	v.x = cdm::clamp(v.x, min.x, max.x);
	v.y = cdm::clamp(v.y, min.y, max.y);
	v.z = cdm::clamp(v.z, min.z, max.z);
	return v;
}
template<typename T>
T dot(vector3<T> lhs, vector3<T> rhs) { return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z; }
template<typename T>
vector3<T> cross(vector3<T> lhs, vector3<T> rhs)
{
	return {
		lhs.y * rhs.z - lhs.z * rhs.y,
		lhs.z * rhs.x - lhs.x * rhs.z,
		lhs.x * rhs.y - lhs.y * rhs.x
	};
}
template<typename T>
vector3<T> lerp(vector3<T> begin, vector3<T> end, T percent) { return (end - begin) * percent + begin; }
template<typename T>
vector3<T> nlerp(vector3<T> begin, vector3<T> end, T percent) { return normalize(lerp(begin, end, percent)); }
template<typename T>
T distance_between(vector3<T> v1, vector3<T> v2) { return norm(v1 - v2); }
template<typename T>
T distance_squared_between(vector3<T> v1, vector3<T> v2) { return norm_squared(v1 - v2); }
template<typename T>
vector3<T> from_to(vector3<T> from, vector3<T> to) { return {to.x - from.x, to.y - from.y, to.z - from.z}; }
template<typename T>
radian<T> angle_between(vector3<T> v1, vector3<T> v2)
{
	T divisor = std::sqrt(norm_squared(v1) * norm_squared(v2));
	T alpha = dot(v1, v2) / divisor;
	return radian<T>(std::acosf(cdm::clamp(alpha, T(-1), T(1))));
}
template<typename T>
bool nearly_equal(vector3<T> v1, vector3<T> v2, T e)
{
	return nearly_equal(v1.x, v2.x, e) &&
	       nearly_equal(v1.y, v2.y, e) &&
	       nearly_equal(v1.z, v2.z, e);
}
#pragma endregion

#pragma region definition vector4
template<typename T>
vector2<T> vector4<T>::xy() const { return { x, y }; }
template<typename T>
vector3<T> vector4<T>::xyz() const { return { x, y, z }; }

template<typename T>
vector4<T> vector4<T>::operator+(vector4<T> v) const { return {x + v.x, y + v.y, z + v.z, w + v.w}; }
template<typename T>
vector4<T> vector4<T>::operator-(vector4<T> v) const { return {x - v.x, y - v.y, z - v.z, w - v.w}; }
template<typename T>
vector4<T> vector4<T>::operator*(T f) const { return {x * f, y * f, z * f, w * f}; }
template<typename T>
vector4<T> vector4<T>::operator/(T f) const { return {x / f, y / f, z / f, w / f}; }

template<typename T>
vector4<T>& vector4<T>::operator+=(vector4<T> v) { *this = *this + v; return *this; }
template<typename T>
vector4<T>& vector4<T>::operator-=(vector4<T> v) { *this = *this - v; return *this; }
template<typename T>
vector4<T>& vector4<T>::operator*=(T f) { *this = *this * f; return *this; }
template<typename T>
vector4<T>& vector4<T>::operator/=(T f) { *this = *this / f; return *this; }

template<typename T>
vector4<T> vector4<T>::operator-() const { return {-x, -y, -z, -w}; }

template<typename T>
bool vector4<T>::operator==(vector4<T> v) const { return x == v.x && y == v.y && z == v.z && w == v.w; }
template<typename T>
bool vector4<T>::operator!=(vector4<T> v) const { return !operator==(v); }

template<typename T>
T norm(vector4<T> v) { return std::sqrt(norm_squared(v)); }
template<typename T>
T norm_squared(vector4<T> v) { return v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w; }
template<typename T>
vector4<T> normalize(vector4<T> v)
{
	T n = norm(v);
	v.x /= n;
	v.y /= n;
	v.z /= n;
	v.w /= n;
	return v;
}
template<typename T>
vector4<T> clamp(vector4<T> v, vector4<T> min, vector4<T> max)
{
	v.x = cdm::clamp(v.x, min.x, max.x);
	v.y = cdm::clamp(v.y, min.y, max.y);
	v.z = cdm::clamp(v.z, min.z, max.z);
	v.w = cdm::clamp(v.w, min.w, max.w);
	return v;
}
template<typename T>
vector4<T> negate(vector4<T> v)
{
	v.x = -v.x;
	v.y = -v.y;
	v.z = -v.z;
	v.w = -v.w;
	return v;
}
template<typename T>
T dot(vector4<T> lhs, vector4<T> rhs) { return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w; }
template<typename T>
vector4<T> lerp(vector4<T> begin, vector4<T> end, T percent) { return (end - begin) * percent + begin; }
template<typename T>
vector4<T> nlerp(vector4<T> begin, vector4<T> end, T percent) { return normalized(lerp(begin, end, percent)); }
template<typename T>
T distance_between(vector4<T> v1, vector4<T> v2) { return norm(v1 - v2); }
template<typename T>
T distance_squared_between(vector4<T> v1, vector4<T> v2) { return norm_squared(v1 - v2); }
template<typename T>
vector4<T> from_to(vector4<T> from, vector4<T> to) { return {to.x - from.x, to.y - from.y, to.z - from.z, to.w - from.w}; }
template<typename T>
bool nearly_equal(vector4<T> v1, vector4<T> v2, T e)
{
	return nearly_equal(v1.x, v2.x, e) &&
	       nearly_equal(v1.y, v2.y, e) &&
	       nearly_equal(v1.z, v2.z, e) &&
	       nearly_equal(v1.w, v2.w, e);
}
#pragma endregion

#pragma region definition normalized
template<typename T>
normalized<T>::normalized(const T& t) : vector(normalize(t)) {}
//template<typename T>
//normalized<T>::normalized(T&& t) noexcept : vector(normalize(std::move(t)) {}

template<typename T>
normalized<T> normalized<T>::already_normalized(const T& t)
{
	normalized res;
	res.vector = t;
	return res;
}
template<typename T>
normalized<T> normalized<T>::already_normalized(T&& t) noexcept
{
	normalized res;
	res.vector = std::move(t);
	return res;
}

template<typename T>
const T& normalized<T>::operator*() const { return vector; }
template<typename T>
const T* normalized<T>::operator->() const { return &vector; }
template<typename T>
normalized<T>::operator const T&() const { return vector; }

template<typename T>
normalized<T> normalized<T>::operator+(const T& v) const { return normalize(vector + v); }
template<typename T>
normalized<T> normalized<T>::operator-(const T& v) const { return normalize(vector - v); }
template<typename T>
T normalized<T>::operator*(T f) const { return vector * f; }
template<typename T>
T normalized<T>::operator/(T f) const { return vector / f; }

template<typename T>
normalized<T>& normalized<T>::operator+=(const T& v) { vector = normalize(vector + v); return *this; }
template<typename T>
normalized<T>& normalized<T>::operator-=(const T& v) { vector = normalize(vector - v); return *this; }

template<typename T>
normalized<T> normalized<T>::operator-() const { return -vector; }

template<typename T>
bool normalized<T>::operator==(const T& v) const { return vector == v; }
template<typename T>
bool normalized<T>::operator!=(const T& v) const { return vector != v; }

template<typename T>
normalized<T>& normalized<T>::operator=(const T& t) { vector = normalize(t); return *this; }
template<typename T>
normalized<T>& normalized<T>::operator=(T&& t) noexcept { vector = normalize(t); return *this; }
#pragma endregion

#pragma region definition matrix2
template<typename T>
matrix2<T> matrix2<T>::zero() { return {T(0), T(0), T(0), T(0)}; }
template<typename T>
matrix2<T> matrix2<T>::identity() { return {T(1), T(0), T(0), T(1)}; }
template<typename T>
matrix2<T> matrix2<T>::rotation(radian<T> angle)
{
	T c = cos(angle);
	T s = sin(angle);
	return {
		c, -s,
		s, c
	};
}
template<typename T>
matrix2<T> matrix2<T>::rotation(normalized<complex<T>> angle)
{
	return {
		angle->r, -angle->i,
		angle->i, angle->r
	};
}

template<typename T>
T& matrix2<T>::at(uint8_t x, uint8_t y) { return reinterpret_cast<T*>(this)[x + 2 * y]; }
template<typename T>
const T& matrix2<T>::at(uint8_t x, uint8_t y) const { return reinterpret_cast<const T*>(this)[x + 2 * y]; }

template<typename T>
vector2<T> matrix2<T>::operator*(vector2<T> v) const
{
	return {
		m00 * v.x + m01 * v.y,
		m10 * v.x + m11 * v.y
	};
}
template<typename T>
matrix2<T> matrix2<T>::operator*(matrix2<T> m) const
{
	return {
		m00 * m.m00 + m01 * m.m10,
		m00 * m.m01 + m01 * m.m11,

		m10 * m.m00 + m11 * m.m10,
		m10 * m.m01 + m11 * m.m11
	};
}

template<typename T>
matrix2<T> transpose(matrix2<T> m)
{
	std::swap(m.m01, m.m10);
	return m;
}
#pragma endregion

#pragma region definition matrix3
template<typename T>
matrix3<T> matrix3<T>::zero() { return {T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0)}; }
template<typename T>
matrix3<T> matrix3<T>::identity() { return {T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1)}; }
template<typename T>
matrix3<T> matrix3<T>::rotation(euler_angles<T> r)
{
	matrix3<T> RX = matrix3<T>::identity();
	RX.m11 = cos(r.x);
	RX.m12 = -sin(r.x);
	RX.m21 = sin(r.x);
	RX.m22 = cos(r.x);

	matrix3<T> RY = matrix3<T>::identity();
	RY.m00 = cos(r.y);
	RY.m02 = sin(r.y);
	RY.m20 = -sin(r.y);
	RY.m22 = cos(r.y);

	matrix3<T> RZ = matrix3<T>::identity();
	RZ.m00 = cos(r.z);
	RZ.m01 = -sin(r.z);
	RZ.m10 = sin(r.z);
	RZ.m11 = cos(r.z);

	//return RZ * RX * RY;
	//return RZ * RY * RX;
	return RY * RX * RZ;
	//return RY * RZ * RX;
	//return RY * RX * RZ;
	//return RX * RZ * RY;
}
template<typename T>
matrix3<T> matrix3<T>::rotation(quaternion<T> q)
{
	return {
		detail::get_quaternion_matrix_element<0, 0>(q), detail::get_quaternion_matrix_element<1, 0>(q), detail::get_quaternion_matrix_element<2, 0>(q),
		detail::get_quaternion_matrix_element<0, 1>(q), detail::get_quaternion_matrix_element<1, 1>(q), detail::get_quaternion_matrix_element<2, 1>(q),
		detail::get_quaternion_matrix_element<0, 2>(q), detail::get_quaternion_matrix_element<1, 2>(q), detail::get_quaternion_matrix_element<2, 2>(q),
	};
}
template<typename T>
matrix3<T> matrix3<T>::rotation_around_x(radian<T> angle)
{
	T c = cos(angle);
	T s = sin(angle);
	return {
		T(1), T(0), T(0),
		T(0), c,    s,
		T(0), -s,   c
	};
}
template<typename T>
matrix3<T> matrix3<T>::rotation_around_y(radian<T> angle)
{
	T c = cos(angle);
	T s = sin(angle);
	return {
		c,    T(0), s,
		T(0), T(1), T(0),
		-s,   T(0), c
	};
}
template<typename T>
matrix3<T> matrix3<T>::rotation_around_z(radian<T> angle)
{
	T c = cos(angle);
	T s = sin(angle);
	return {
		c,    -s,   T(0),
		s,    c,    T(0),
		T(0), T(0), T(1)
	};
}
template<typename T>
matrix3<T> matrix3<T>::rotation_around_x(normalized<complex<T>> angle)
{
	T c = angle->r;
	T s = angle->i;
	return {
		T(1), T(0), T(0),
		T(0), c,    s,
		T(0), -s,   c
	};
}
template<typename T>
matrix3<T> matrix3<T>::rotation_around_y(normalized<complex<T>> angle)
{
	T c = angle->r;
	T s = angle->i;
	return {
		c,    T(0), s,
		T(0), T(1), T(0),
		-s,   T(0), c
	};
}
template<typename T>
matrix3<T> matrix3<T>::rotation_around_z(normalized<complex<T>> angle)
{
	T c = angle->r;
	T s = angle->i;
	return {
		c,    -s,   T(0),
		s,    c,    T(0),
		T(0), T(0), T(1)
	};
}
template<typename T>
matrix3<T> matrix3<T>::inversed() const
{
	matrix3<T> res = *this;
	res.inverse();
	return res;
}
template<typename T>
matrix3<T> matrix3<T>::transposed() const
{
	return {
		m00, m10, m20,
		m01, m11, m21,
		m02, m12, m22
	};
}
template<typename T>
T matrix3<T>::determinant() const
{
	return
		m00 * m11 * m22 +
		m01 * m12 * m20 +
		m02 * m10 * m21 -
		m20 * m11 * m02 -
		m21 * m12 * m00 -
		m22 * m10 * m01;
}
template<typename T>
bool matrix3<T>::is_orthogonal() const
{
	return
		nearly_equal(m00 * m01 + m10 * m11 + m20 * m21, T(0)) &&
		nearly_equal(m00 * m02 + m10 * m12 + m20 * m22, T(0)) &&
		nearly_equal(m01 * m02 + m11 * m12 + m21 * m22, T(0)) &&
		nearly_equal(m00 * m00 + m10 * m10 + m20 * m20, T(1)) &&
		nearly_equal(m01 * m01 + m11 * m11 + m21 * m21, T(1)) &&
		nearly_equal(m02 * m02 + m12 * m12 + m22 * m22, T(1));
}

template<typename T>
T& matrix3<T>::at(uint8_t x, uint8_t y) { return reinterpret_cast<T*>(this)[x + 3 * y]; }
template<typename T>
const T& matrix3<T>::at(uint8_t x, uint8_t y) const { return reinterpret_cast<const T*>(this)[x + 3 * y]; }

template<typename T>
matrix3<T> matrix3<T>::operator*(T f) const
{
	return {
		m00 * f, m10 * f, m20 * f,
		m01 * f, m11 * f, m21 * f,
		m02 * f, m12 * f, m22 * f
	};
}
template<typename T>
vector3<T> matrix3<T>::operator*(vector3<T> v) const
{
	return {
		m00 * v.x + m10 * v.y + m20 * v.z,
		m01 * v.x + m11 * v.y + m21 * v.z,
		m02 * v.x + m12 * v.y + m22 * v.z
	};
}
template<typename T>
matrix3<T> matrix3<T>::operator*(const matrix3<T>& m) const
{
	return {
		m00 * m.m00 + m10 * m.m01 + m20 * m.m02,
		m00 * m.m10 + m10 * m.m11 + m20 * m.m12,
		m00 * m.m20 + m10 * m.m21 + m20 * m.m22,

		m01 * m.m00 + m11 * m.m01 + m21 * m.m02,
		m01 * m.m10 + m11 * m.m11 + m21 * m.m12,
		m01 * m.m20 + m11 * m.m21 + m21 * m.m22,

		m02 * m.m00 + m12 * m.m01 + m22 * m.m02,
		m02 * m.m10 + m12 * m.m11 + m22 * m.m12,
		m02 * m.m20 + m12 * m.m21 + m22 * m.m22
	};
}

template<typename T>
matrix3<T> inverse(matrix3<T> m)
{
	if (m.is_orthogonal())
	{
		return transpose(m);
	}

	T det = m.determinant();
	T recM00 = m.m00;
	T recM10 = m.m10;
	T recM20 = m.m20;
	T recM01 = m.m01;
	T recM11 = m.m11;
	T recM02 = m.m02;
	m.m00 = (recM11 * m.m22 - m.m21 * m.m12) / det;
	m.m10 = (m.m12 * recM20 - recM10 * m.m22) / det;
	m.m20 = (recM10 * m.m12 - recM20 * recM11) / det;
	m.m01 = (recM02 * m.m21 - recM01 * m.m22) / det;
	m.m11 = (recM00 * m.m22 - recM20 * recM02) / det;
	m.m21 = (recM01 * recM20 - recM00 * m.m21) / det;
	m.m02 = (recM01 * m.m12 - recM11 * recM02) / det;
	m.m12 = (recM02 * recM10 - recM00 * m.m12) / det;
	m.m22 = (recM00 * recM11 - recM10 * recM01) / det;
	return m;
}
template<typename T>
matrix3<T> transpose(matrix3<T> m)
{
	std::swap(m.m01, m.m10);
	std::swap(m.m02, m.m20);
	std::swap(m.m12, m.m21);
	return m;
}
#pragma endregion

#pragma region definition matrix4
template<typename T>
matrix4<T> matrix4<T>::zero() { return {T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0)}; }
template<typename T>
matrix4<T> matrix4<T>::identity() { return {T(1), T(0), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(0), T(1)}; }
template<typename T>
matrix4<T> matrix4<T>::rotation(euler_angles<T> r) { return matrix4(matrix3::rotation(r)); }
template<typename T>
matrix4<T> matrix4<T>::rotation(quaternion<T> q) { return matrix4(matrix3::rotation(q)); }
template<typename T>
matrix4<T> matrix4<T>::translation(vector3<T> t) { return {T(1), T(0), T(0), t.x, T(0), T(1), T(0), t.y, T(0), T(0), T(1), t.z, T(0), T(0), T(0), T(1)}; }
template<typename T>
matrix4<T> matrix4<T>::translation(T x, T y, T z) { return translation({x, y, z}); }
template<typename T>
matrix4<T> matrix4<T>::scale(vector3<T> s) { return {s.x, T(0), T(0), T(0), T(0), s.y, T(0), T(0), T(0), T(0), s.z, T(0), T(0), T(0), T(0), T(1)}; }
template<typename T>
matrix4<T> matrix4<T>::look_at(vector3<T> from, vector3<T> to, vector3<T> up)
{
	vector3<T> forward = (from - to).get_normalized();
	vector3<T> right = cross(up.get_normalized(), forward);
	vector3<T> true_up = cross(forward, right);

	return {
		right.x, true_up.x, forward.x, from.x,
		right.y, true_up.y, forward.y, from.y,
		right.z, true_up.z, forward.z, from.z,
		T(0),    T(0),      T(0),      T(1)
	};
}
template<typename T>
matrix4<T> matrix4<T>::orthographic(T left, T right, T top, T bottom, T near, T far)
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
template<typename T>
matrix4<T> matrix4<T>::rotation_around_x(radian<T> angle) { return matrix4<T>(matrix3<T>::rotation_around_x(angle)); }
template<typename T>
matrix4<T> matrix4<T>::rotation_around_y(radian<T> angle) { return matrix4<T>(matrix3<T>::rotation_around_y(angle)); }
template<typename T>
matrix4<T> matrix4<T>::rotation_around_z(radian<T> angle) { return matrix4<T>(matrix3<T>::rotation_around_z(angle)); }
template<typename T>
matrix4<T> matrix4<T>::rotation_around_x(normalized<complex<T>> angle) { return matrix4<T>(matrix3<T>::rotation_around_x(angle)); }
template<typename T>
matrix4<T> matrix4<T>::rotation_around_y(normalized<complex<T>> angle) { return matrix4<T>(matrix3<T>::rotation_around_y(angle)); }
template<typename T>
matrix4<T> matrix4<T>::rotation_around_z(normalized<complex<T>> angle) { return matrix4<T>(matrix3<T>::rotation_around_z(angle)); }
template<typename T>
bool matrix4<T>::is_orthogonal() const
{
	return
		nearly_equal(m00 * m01 + m10 * m11 + m20 * m21 + m30 * m31, T(0)) &&
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
template<typename T>
bool matrix4<T>::is_homogenous() const
{
	return
		nearly_equal(m03, T(0)) &&
		nearly_equal(m13, T(0)) &&
		nearly_equal(m23, T(0)) &&
		nearly_equal(m30, T(0)) &&
		nearly_equal(m31, T(0)) &&
		nearly_equal(m32, T(0)) &&
		nearly_equal(m33, T(1));
}
template<typename T>
matrix4<T> matrix4<T>::inversed() const
{
	matrix4<T> res = *this;
	res.inverse();
	return res;
}
template<typename T>
matrix4<T> matrix4<T>::transposed() const
{
	matrix4<T> res = *this;
	res.transpose();
	return res;
}
template<typename T>
T matrix4<T>::determinant() const
{
	if (is_homogenous())
	{
		return
			m00 * m11 * m22 +
			m01 * m12 * m20 +
			m02 * m10 * m21 -
			m20 * m11 * m02 -
			m21 * m12 * m00 -
			m22 * m10 * m01;
	}

	T det1 = m11 * (m22 * m33 - m32 * m23) - m21 * (m12 * m33 - m32 * m13) + m31 * (m12 * m23 - m22 * m13);
	T det2 = m01 * (m22 * m33 - m32 * m23) - m21 * (m02 * m33 - m32 * m03) + m31 * (m02 * m23 - m22 * m03);
	T det3 = m01 * (m12 * m33 - m32 * m13) - m11 * (m02 * m33 - m32 * m03) + m31 * (m02 * m13 - m12 * m03);
	T det4 = m01 * (m12 * m23 - m22 * m13) - m11 * (m02 * m23 - m22 * m03) + m21 * (m02 * m13 - m12 * m03);
	return m00 * det1 - m10 * det2 + m20 * det3 - m30 * det4;
}

template<typename T>
T& matrix4<T>::at(uint8_t x, uint8_t y) { return reinterpret_cast<T*>(this)[x + 4 * y]; }
template<typename T>
const T& matrix4<T>::at(uint8_t x, uint8_t y) const { return reinterpret_cast<const T*>(this)[x + 4 * y]; }

template<typename T>
matrix4<T> matrix4<T>::operator*(T f) const
{
	return {
		m00 * f, m10 * f, m20 * f, m30 * f,
		m00 * f, m10 * f, m20 * f, m30 * f,
		m00 * f, m10 * f, m20 * f, m30 * f,
		m00 * f, m10 * f, m20 * f, m30 * f
	};
}
template<typename T>
matrix4<T> matrix4<T>::operator/(T f) const
{
	return {
		m00 / f, m10 / f, m20 / f, m30 / f,
		m00 / f, m10 / f, m20 / f, m30 / f,
		m00 / f, m10 / f, m20 / f, m30 / f,
		m00 / f, m10 / f, m20 / f, m30 / f
	};
}
template<typename T>
vector4<T> matrix4<T>::operator*(vector4<T> v) const
{
	return {
		m00 * v.x + m10 * v.y + m20 * v.z + m30 * v.w,
		m01 * v.x + m11 * v.y + m21 * v.z + m31 * v.w,
		m02 * v.x + m12 * v.y + m22 * v.z + m32 * v.w,
		m03 * v.x + m13 * v.y + m23 * v.z + m33 * v.w
	};
}
template<typename T>
matrix4<T> matrix4<T>::operator*(const matrix4<T>& m) const
{
	return {
		m.m00 * m00 + m.m01 * m10 + m.m02 * m20 + m.m03 * m30,
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
		m.m30 * m03 + m.m31 * m13 + m.m32 * m23 + m.m33 * m33
	};
}

template<typename T>
template<typename U>
std::array<U, 16> matrix4<T>::to_array() const
{
	std::array<U, 16> res {
		U(m00), U(m10), U(m20), U(m30),
		U(m01), U(m11), U(m21), U(m31),
		U(m02), U(m12), U(m22), U(m32),
		U(m03), U(m13), U(m23), U(m33)
	};

	return res;
}

template<typename T>
matrix4<T> inverse(matrix4<T> m)
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

	m.m00 = invDet * (
		recM11 * recM22 * m.m33 -
		recM11 * m.m23 * m.m32 -
		recM21 * recM12 * m.m33 +
		recM21 * recM13 * m.m32 +
		recM31 * recM12 * m.m23 -
		recM31 * recM13 * recM22);

	m.m10 = invDet * (
		-recM10 * recM22 * m.m33 +
		recM10 * m.m23 * m.m32 +
		recM20 * recM12 * m.m33 -
		recM20 * recM13 * m.m32 -
		recM30 * recM12 * m.m23 +
		recM30 * recM13 * recM22
		);

	m.m20 = invDet * (
		recM10 * recM21 * m.m33 -
		recM10 * m.m23 * recM31 -
		recM20 * recM11 * m.m33 +
		recM20 * recM13 * recM31 +
		recM30 * recM11 * m.m23 -
		recM30 * recM13 * recM21
		);

	m.m30 = invDet * (
		-recM10 * recM21 * m.m32 +
		recM10 * recM22 * recM31 +
		recM20 * recM11 * m.m32 -
		recM20 * recM12 * recM31 -
		recM30 * recM11 * recM22 +
		recM30 * recM12 * recM21
		);

	m.m01 = invDet * (
		-recM01 * recM22 * m.m33 +
		recM01 * m.m23 * m.m32 +
		recM21 * recM02 * m.m33 -
		recM21 * recM03 * m.m32 -
		recM31 * recM02 * m.m23 +
		recM31 * recM03 * recM22
		);

	m.m11 = invDet * (
		recM00 * recM22 * m.m33 -
		recM00 * m.m23 * m.m32 -
		recM20 * recM02 * m.m33 +
		recM20 * recM03 * m.m32 +
		recM30 * recM02 * m.m23 -
		recM30 * recM03 * recM22
		);

	m.m21 = invDet * (
		-recM00 * recM21 * m.m33 +
		recM00 * m.m23 * recM31 +
		recM20 * recM01 * m.m33 -
		recM20 * recM03 * recM31 -
		recM30 * recM01 * m.m23 +
		recM30 * recM03 * recM21
		);

	m.m31 = invDet * (
		recM00 * recM21 * m.m32 -
		recM00 * recM22 * recM31 -
		recM20 * recM01 * m.m32 +
		recM20 * recM02 * recM31 +
		recM30 * recM01 * recM22 -
		recM30 * recM02 * recM21
		);

	m.m02 = invDet * (
		recM01 * recM12 * m.m33 -
		recM01 * recM13 * m.m32 -
		recM11 * recM02 * m.m33 +
		recM11 * recM03 * m.m32 +
		recM31 * recM02 * recM13 -
		recM31 * recM03 * recM12
		);

	m.m12 = invDet * (
		-recM00 * recM12 * m.m33 +
		recM00 * recM13 * m.m32 +
		recM10 * recM02 * m.m33 -
		recM10 * recM03 * m.m32 -
		recM30 * recM02 * recM13 +
		recM30 * recM03 * recM12
		);

	m.m22 = invDet * (
		recM00 * recM11 * m.m33 -
		recM00 * recM13 * recM31 -
		recM10 * recM01 * m.m33 +
		recM10 * recM03 * recM31 +
		recM30 * recM01 * recM13 -
		recM30 * recM03 * recM11
		);

	m.m32 = invDet * (
		-recM00 * recM11 * m.m32 +
		recM00 * recM12 * recM31 +
		recM10 * recM01 * m.m32 -
		recM10 * recM02 * recM31 -
		recM30 * recM01 * recM12 +
		recM30 * recM02 * recM11
		);

	m.m03 = invDet * (
		-recM01 * recM12 * m.m23 +
		recM01 * recM13 * recM22 +
		recM11 * recM02 * m.m23 -
		recM11 * recM03 * recM22 -
		recM21 * recM02 * recM13 +
		recM21 * recM03 * recM12
		);

	m.m13 = invDet * (
		recM00 * recM12 * m.m23 -
		recM00 * recM13 * recM22 -
		recM10 * recM02 * m.m23 +
		recM10 * recM03 * recM22 +
		recM20 * recM02 * recM13 -
		recM20 * recM03 * recM12
		);

	m.m23 = invDet * (
		-recM00 * recM11 * m.m23 +
		recM00 * recM13 * recM21 +
		recM10 * recM01 * m.m23 -
		recM10 * recM03 * recM21 -
		recM20 * recM01 * recM13 +
		recM20 * recM03 * recM11
		);

	m.m33 = invDet * (
		recM00 * recM11 * recM22 -
		recM00 * recM12 * recM21 -
		recM10 * recM01 * recM22 +
		recM10 * recM02 * recM21 +
		recM20 * recM01 * recM12 -
		recM20 * recM02 * recM11
		);

	return *this;
}
template<typename T>
matrix4<T> transpose(matrix4<T> m)
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

#pragma region definition perspective
template<typename T>
void perspective<T>::set(radian<T> angle, T ratio, T near, T far)
{
	m_angle = angle;
	m_ratio = ratio;
	m_inv_ratio = T(1) / m_ratio;
	m_near = near;
	m_far = far;
	m_invTanHalfFovy = T(1) / tan(m_angle / T(2));
}
template<typename T>
void perspective<T>::set_angle(radian<T> angle)
{
	m_angle = angle;
	m_invTanHalfFovy = T(1) / tan(angle / T(2));
}
template<typename T>
radian<T> perspective<T>::get_angle() const { return m_angle; }
template<typename T>
void perspective<T>::set_ratio(T ratio)
{
	m_ratio = ratio;
	m_inv_ratio = T(1) / m_ratio;
}
template<typename T>
T perspective<T>::get_ratio() const { return m_ratio; }
template<typename T>
void perspective<T>::set_near(T near_plane) { m_near = near_plane; }
template<typename T>
T perspective<T>::get_near() const { return m_near; }
template<typename T>
void perspective<T>::set_far(T far_plane) { m_far = far_plane; }
template<typename T>
T perspective<T>::get_far() const { return m_far; }
template<typename T>
matrix4<T> perspective<T>::to_matrix4() const
{
	return {
		m_inv_ratio * m_invTanHalfFovy, T(0),              T(0),                     T(0),
		T(0),                           -m_invTanHalfFovy, T(0),                     T(0),
		T(0),                           T(0),              m_far / (m_near - m_far), -(m_far * m_near) / (m_far - m_near),
		T(0),                           T(0),              -T(1),                    T(0)
	};
}

template<typename T>
matrix4<T> operator*(const matrix4<T>& m, const perspective<T>& p)
{
	T a = p.m_inv_ratio * p.m_invTanHalfFovy;
	T b = -p.m_invTanHalfFovy;
	T c = p.m_far / (p.m_near - p.m_far);
	T d = -(p.m_far * p.m_near) / (p.m_far - p.m_near);
	return {
		a * m.m00,
		b * m.m10,
		c * m.m20 - m.m30,
		d * m.m20,

		a * m.m01,
		b * m.m11,
		c * m.m21 - m.m31,
		d * m.m21,

		a * m.m02,
		b * m.m12,
		c * m.m22 - m.m32,
		d * m.m22,

		a * m.m03,
		b * m.m13,
		c * m.m23 - m.m33,
		d * m.m23
	};
}
template<typename T>
matrix4<T> operator*(const perspective<T>& p, const matrix4<T>& m)
{
	T a = p.m_inv_ratio * p.m_invTanHalfFovy;
	T b = -p.m_invTanHalfFovy;
	T c = p.m_far / (p.m_near - p.m_far);
	T d = -(p.m_far * p.m_near) / (p.m_far - p.m_near);
	return {
		m.m00 * a,
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
		-m.m32
	};
}

template<typename T>
matrix4<T> operator*(const unscaled_transform3d<T>& t, const perspective<T>& p)
{
	T a = p.m_inv_ratio * p.m_invTanHalfFovy;
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
		-T(1),
		T(0)
	};
}
template<typename T>
matrix4<T> operator*(const perspective<T>& p, const unscaled_transform3d<T>& t)
{
	T a = p.m_inv_ratio * p.m_invTanHalfFovy;
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
	return {
		(T(1) - T(2) * (yy + zz)) * a,
		(T(2) * (xy - wz)) * a,
		(T(2) * (xz + wy)) * a,
		t.position.x * a,

		(T(2) * (xy + wz)) * b,
		(T(1) - T(2) * (xx + zz)) * b,
		(T(2) * (yz - wx)) * b,
		t.position.y * b,

		(T(2) * (xz - wy)) * c,
		(T(2) * (yz + wx)) * c,
		(T(1) - T(2) * (xx + yy)) * c,
		t.position.z * c + d,

		-(T(2) * (xz - wy)),
		-(T(2) * (yz + wx)),
		-(T(1) - T(2) * (xx + yy)),
		-t.position.z
	};
}
#pragma endregion

#pragma region definition quaternion
template<typename T>
quaternion<T> quaternion<T>::zero() { return {T(0), T(0), T(0), T(0)}; }
template<typename T>
quaternion<T> quaternion<T>::identity() { return {T(0), T(0), T(0), T(1)}; }

template<typename T>
quaternion<T> quaternion<T>::inversed() const
{
	quaternion<T> res = *this;
	inverse(res);
	return res;
}
template<typename T>
quaternion<T> quaternion<T>::conjugated() const
{
	quaternion<T> res = *this;
	conjugate(res);
	return res;
}
template<typename T>
quaternion<T> quaternion<T>::normalized() const
{
	quaternion<T> res = *this;
	normalize(res);
	return res;
}
template<typename T>
quaternion<T> quaternion<T>::clamped(quaternion<T> min, quaternion<T> max) const
{
	quaternion<T> res = *this;
	clamp(res, min, max);
	return res;
}

template<typename T>
quaternion<T> quaternion<T>::operator+(quaternion<T> q) const { return {x + q.x, y + q.y, z + q.z, w + q.w}; }
template<typename T>
quaternion<T> quaternion<T>::operator-(quaternion<T> q) const { return {x - q.x, y - q.y, z - q.z, w - q.w}; }
template<typename T>
quaternion<T> quaternion<T>::operator*(quaternion<T> q) const
{
	quaternion<T> res;
	res.x =  x * q.w + y * q.z - z * q.y + w * q.x;
	res.y = -x * q.z + y * q.w + z * q.x + w * q.y;
	res.z =  x * q.y - y * q.x + z * q.w + w * q.z;
	res.w = -x * q.x - y * q.y - z * q.z + w * q.w;
	return res;
}
template<typename T>
quaternion<T> quaternion<T>::operator*(T f) const { return {x * f, y * f, z * f, w * f}; }
template<typename T>
quaternion<T> quaternion<T>::operator/(T f) const { return {x / f, y / f, z / f, w / f}; }

template<typename T>
quaternion<T>& quaternion<T>::operator+=(quaternion<T> q) { *this = *this + q; return *this; }
template<typename T>
quaternion<T>& quaternion<T>::operator-=(quaternion<T> q) { *this = *this - q; return *this; }
template<typename T>
quaternion<T>& quaternion<T>::operator*=(quaternion<T> q) { *this = *this * q; return *this; }
template<typename T>
quaternion<T>& quaternion<T>::operator*=(T f) { *this = *this * f; return *this; }
template<typename T>
quaternion<T>& quaternion<T>::operator/=(T f) { *this = *this / f; return *this; }

template<typename T>
quaternion<T> quaternion<T>::operator-() const { return get_negated(); }

template<typename T>
bool quaternion<T>::operator==(quaternion<T> q) const { return x == q.x && y == q.y && z == q.z && w == q.w; }
template<typename T>
bool quaternion<T>::operator!=(quaternion<T> q) const { return !operator==(q); }

template<typename T>
vector3<T> operator*(const normalized<quaternion<T>>& q, vector3<T> v)
{
	vector3 qvec = {q->x, q->y, q->z};
	vector3 uv = cross(qvec, v);
	vector3 uuv = cross(qvec, uv);
	uv = uv * (T(2) * q->w);
	uuv = uuv * T(2);
	return v + uv + uuv;
}

template<typename T>
T norm(quaternion<T> q) { return std::sqrt(norm_squared(q)); }
template<typename T>
T norm_squared(quaternion<T> q) { return q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w; }
template<typename T>
quaternion<T> normalize(quaternion<T> q)
{
	T n = norm(q);
	q.x /= n;
	q.y /= n;
	q.z /= n;
	q.w /= n;
	return q;
}
template<typename T>
quaternion<T> clamp(quaternion<T> q, quaternion<T> min, quaternion<T> max)
{
	q.x = cdm::clamp(q.x, min.x, max.x);
	q.y = cdm::clamp(q.y, min.y, max.y);
	q.z = cdm::clamp(q.z, min.z, max.z);
	q.w = cdm::clamp(q.w, min.w, max.w);
	return q;
}
template<typename T>
quaternion<T> inverse(quaternion<T> q)
{
	T n = norm(q);
	if (nearly_equal(n, T(1)))
		return conjugate(q);

	q.x /= -n;
	q.y /= -n;
	q.z /= -n;
	q.w /= n;
	return q;
}
template<typename T>
quaternion<T> conjugate(quaternion<T> q)
{
	q.x = -q.x;
	q.y = -q.y;
	q.z = -q.z;
	return q;
}
template<typename T>
T dot(quaternion<T> lhs, quaternion<T> rhs) { return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z + lhs.w * rhs.w; }
template<typename T>
quaternion<T> cross(quaternion<T> lhs, quaternion<T> rhs)
{
	return {
		lhs.y * rhs.z - lhs.z * rhs.y,
		lhs.z * rhs.x - lhs.x * rhs.z,
		lhs.x * rhs.y - lhs.y * rhs.x,
		T(1)
	};
}
template<typename T>
quaternion<T> lerp(quaternion<T> begin, quaternion<T> end, T percent) { return (end - begin) * percent + begin; }
template<typename T>
quaternion<T> nlerp(quaternion<T> begin, quaternion<T> end, T percent) { return lerp(begin, end, percent).get_normalized(); }
template<typename T>
quaternion<T> slerp(quaternion<T> begin, quaternion<T> end, T percent)
{
	quaternion<T> _end;
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

	quaternion<T> res = begin * d;
	res = _end - res;

	res = begin * cosf(theta) + res * sinf(theta);

	return res;
}
#pragma endregion

#pragma region definition cartesian_direction2d
template<typename T>
cartesian_direction2d<T> cartesian_direction2d<T>::from(const normalized<vector2<T>>& direction)
{
	cartesian_direction2d<T> d;
	d.x = direction->x;
	d.y = direction->y;
	return d;
}
template<typename T>
cartesian_direction2d<T> cartesian_direction2d<T>::from(polar_direction<T> direction) { return from(direction.angle); }
template<typename T>
cartesian_direction2d<T> cartesian_direction2d<T>::from(radian<T> angle)
{
	cartesian_direction2d<T> d;
	d.x = cos(angle);
	d.y = sin(angle);
	return d;
}

template<typename T>
cartesian_direction2d<T>& cartesian_direction2d<T>::rotate(radian<T> angle)
{
	const T c = cos(angle);
	const T s = sin(angle);
    x = x * c - y * s;
    y = x * s + y * c;

	return *this;
}
#pragma endregion

#pragma region definition polar_direction
template<typename T>
polar_direction<T> polar_direction<T>::from(const normalized<vector2<T>>& direction)
{
	angle = radian<T>(atan(direction->x));
}
template<typename T>
polar_direction<T> polar_direction<T>::from(cartesian_direction2d<T> direction)
{
	angle = radian<T>(atan(direction.x));
}
template<typename T>
polar_direction<T> polar_direction<T>::from(radian<T> angle_) { angle = angle_; }

template<typename T>
polar_direction<T>& polar_direction<T>::rotate(radian<T> angle_)
{
	angle += angle_;
	return *this;
}
template<typename T>
polar_direction<T>& polar_direction<T>::wrap()
{
	int s = sign(T(angle));
	while (angle > radian(pi) || angle <= -radian(pi))
		angle -= radian(pi) * T(s);
	return *this;
}
#pragma endregion

#pragma region definition line
template<typename T>
bool are_parallel(line<T> l1, line<T> l2)
{
	return nearly_equal(l1.coefficient, l2.coefficient) || nearly_equal(l1.coefficient, -l2.coefficient);
}

template<typename T>
bool collides(line<T> l1, line<T> l2)
{
	return !line<T>::are_parallel(l1, l2);
}
template<typename T>
bool collides(line<T> l1, line<T> l2, vector2<T>& intersection)
{
	bool collision = collides(l1, l2);

	if (collision)
	{
		T a = l1.coefficient - l2.coefficient;
		T b = l2.offset - l1.offset;
		intersection.x = b / a;
		intersection.y = l1.coefficient * intersection.x + l1.offset;

		//assert(intersection.y == (l2.coefficient * intersection.x + l2.offset));
	}

	return collision;
}
template<typename T>
bool collides(line<T> l, vector2<T> v)
{
	return nearly_equal(l.coefficient * v.x + l.offset, v.y);
}
template<typename T>
bool collides(vector2<T> v, line<T> l)
{
	return collides(l, v);
}

template<typename T>
T distance_between(plane<T> p, vector3<T> v)
{
	return -v.x * p.normal->x - v.y * p.normal->y - v.z * p.normal->z;
}
template<typename T>
T distance_between(vector3<T> v, plane<T> p)
{
	return distance_between(p, v);
}
#pragma endregion

#pragma region definition aa_rect
template<typename T>
bool aa_rect<T>::contains(vector2<T> v) const
{
	return (v.x >= origin.x) && (v.x <= origin.x + dimention.x) && (v.y >= origin.y) && (v.y <= origin.y + dimention.y);
}

template<typename T>
bool collides(aa_rect<T> r1, aa_rect<T> r2)
{
	return collides(r1, r2.origin) ||
		collides(r1, r2.origin + vector2(r2.dimention.x, 0)) ||
		collides(r1, r2.origin + vector2(0, r2.dimention.y)) ||
		collides(r1, r2.origin + r2.dimention);
}

//bool collides(const ray3d& r, const plane& p)
//{
//	T DdotN = r.direction->dot(p.normal);
//	if (std::abs(DdotN) > epsilon)
//	{
//		T t = -(r.origin.dot(p.normal) + p.distance) / DdotN;
//		return t >= 0;
//	}
//	return false;
//}
//bool collides(const ray3d& r, const plane& p, vector3& intersection)
//{
//	T DdotN = r.direction->dot(p.normal);
//	if (std::abs(DdotN) > epsilon)
//	{
//		T t = -(r.origin.dot(p.normal) + p.distance) / DdotN;
//		intersection = r.origin + t * r.direction;
//		return t >= 0;
//	}
//	return false;
//}
//bool collides(const plane& p, const ray3d& r) { return collides(r, p); }
//bool collides(const plane& p, const ray3d& r, vector3& intersection) { return collides(r, p, intersection); }

// https://stackoverflow.com/a/565282
template<typename T>
bool intersects(segment2d<T> s0,
                segment2d<T> s1,
                vector2<T>& outPoint,
                T e)
{
	vector2<T> p = s0.origin;
	vector2<T> q = s1.origin;
	vector2<T> r = from_to(p, s0.end);
	vector2<T> s = from_to(q, s1.end);

	vector2<T> qmp = q - p;

	T rcs = cross(r, s);

	if (nearly_equal(rcs, T(0), e)) // colinear
	{
		if (nearly_equal(cross(qmp, r), T(0), e)) // same line, different segments
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

#pragma region definition plane
template<typename T>
T plane<T>::evaluate(vector3<T> point) const
{
	return normal->x * (point.x - origin.x) +
	       normal->y * (point.y - origin.y) +
	       normal->z * (point.z - origin.z);
}
template<typename T>
vector3<T> plane<T>::project3d(vector3<T> point) const
{
	T distance = evaluate(point);
	return point - normal * distance;
}
template<typename T>
vector2<T> plane<T>::project2d(vector3<T> point, normalized<vector3<T>> plane_tangent) const
{
	normalized<vector3<T>> bitangent = cross(*normal, plane_tangent);
	normalized<vector3<T>> tangent = cross(*bitangent, normal);

	matrix3<T> TBN(tangent->x,   tangent->y,   tangent->z,
	            bitangent->x, bitangent->y, bitangent->z,
	            normal->x,    normal->y,    normal->z
	);

	vector3<T> plane_space_point = point - origin;

	return (TBN * plane_space_point).xy();
}
template<typename T>
vector3<T> plane<T>::unproject(vector2<T> point, normalized<vector3<T>> plane_tangent) const
{
	normalized<vector3<T>> bitangent = cross(*normal, plane_tangent);
	normalized<vector3<T>> tangent = cross(*bitangent, normal);

	matrix3<T> invTBN = matrix3(tangent->x,   tangent->y,   tangent->z,
	                         bitangent->x, bitangent->y, bitangent->z,
	                         normal->x,    normal->y,    normal->z
	).get_inversed();

	vector3<T> plane_space_point = vector3<T>{ point.x, point.y, T(0) };

	return (invTBN * plane_space_point) + origin;
}
template<typename T>
normalized<vector3<T>> plane<T>::computeTangent(T e) const
{
	vector3<T> tangent = project3d(vector3{T(1), T(0), T(0)} + origin) - origin;
	if (norm_squared(tangent) < e)
		tangent = project3d(vector3{T(0), T(1), T(0)} + origin) - origin;
	return normalized<vector3<T>>(tangent);
}

template<typename T>
bool collides(const plane<T>& plane, ray3d<T> r, vector3<T>& collision_point)
{
	T denom = dot(*plane.normal, r.direction);
	if (abs(denom) > T(0.0001))
	{
		T t = -dot(plane.origin - r.origin, plane.normal) / denom;
		if (t >= T(0))
		{
			collision_point = r.origin + r.direction * t;
			return true;
		}
	}
	return false;
}
template<typename T>
bool collides_bidirectional(const plane<T>& plane, ray3d<T> r, vector3<T>& collision_point)
{
	T denom = dot(*plane.normal, r.direction);
	if (abs(denom) > T(0.0001))
	{
		T t = -dot(plane.origin - r.origin, plane.normal) / denom;

		collision_point = r.origin + r.direction * t;
		return true;
	}
	return false;
}
//bool collides(const plane& p, const ray3d& r, vector3& collision_point)
//{
//	p.
//}
//bool collides(const ray3d& r, const plane& p, vector3& collision_point)
//{
//	return collides(p, r, collision_point);
//}
#pragma endregion

#pragma region definition aabb
template<typename T>
bool aabb<T>::contains(vector3<T> p) const
{
	return p.x >= min.x && p.x <= max.x &&
		   p.y >= min.y && p.y <= max.y &&
	       p.z >= min.z && p.z <= max.z;
}template<typename T>
vector3<T> aabb<T>::get_center() const { return (max + min) / T(2); }

template<typename T>
aabb<T> aabb<T>::operator+(aabb<T> rhs) const
{
	aabb<T> res;
	res.min.x = std::min(min.x, rhs.min.x);
	res.min.y = std::min(min.y, rhs.min.y);
	res.min.z = std::min(min.z, rhs.min.z);
	res.max.x = std::max(max.x, rhs.max.x);
	res.max.y = std::max(max.y, rhs.max.y);
	res.max.z = std::max(max.z, rhs.max.z);

	return res;
}

template<typename T>
bool collides(aabb<T> b, ray3d<T> r)
{
	vector3<T> inv
	{
		T(1) / r.direction->x,
		T(1) / r.direction->y,
		T(1) / r.direction->z
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
template<typename T>
bool collides(aabb<T> b1, aabb<T> b2)
{
	if (b1.contains(b2.min))
		return true;
	if (b1.contains(b2.max))
		return true;
	if (b2.contains(b1.min))
		return true;
	if (b2.contains(b1.max))
		return true;

	std::array<vector3<T>, 6> otherPoints;
	otherPoints[0] = { b1.min.x, b1.min.y, b1.max.z };
	otherPoints[1] = { b1.min.x, b1.max.y, b1.min.z };
	otherPoints[2] = { b1.max.x, b1.min.y, b1.min.z };
	otherPoints[3] = { b1.min.x, b1.max.y, b1.max.z };
	otherPoints[4] = { b1.max.x, b1.min.y, b1.max.z };
	otherPoints[5] = { b1.max.x, b1.max.y, b1.min.z };

	for (auto& p : otherPoints)
		if (b2.contains(p))
			return true;

	otherPoints[0] = { b2.min.x, b2.min.y, b2.max.z };
	otherPoints[1] = { b2.min.x, b2.max.y, b2.min.z };
	otherPoints[2] = { b2.max.x, b2.min.y, b2.min.z };
	otherPoints[3] = { b2.min.x, b2.max.y, b2.max.z };
	otherPoints[4] = { b2.max.x, b2.min.y, b2.max.z };
	otherPoints[5] = { b2.max.x, b2.max.y, b2.min.z };

	for (auto& p : otherPoints)
		if (b1.contains(p))
			return true;

	return false;
}
#pragma endregion

#pragma region definition transform2d
template<typename T>
transform2d<T> transform2d<T>::operator*(const transform2d<T>& t) const
{
	transform2d<T> res;
	res.position = rotation * vector2(scale.x * t.position.x, scale.y * t.position.y) + position;
	res.rotation = rotation + t.rotation;
	res.scale = vector2(scale.x * t.scale.x, scale.y * t.scale.y);
	return res;
}
template<typename T>
vector2<T> transform2d<T>::operator*(vector2<T> v) const
{
	return rotation * vector2(scale.x * v.x, scale.y * v.y) + position;
}
#pragma endregion

#pragma region definition transform3d
template<typename T>
transform3d<T> transform3d<T>::operator*(const transform3d<T>& t) const
{
	transform3d<T> res;
	res.position = rotation * vector3(scale.x * t.position.x, scale.y * t.position.y, scale.z * t.position.z) + position;
	res.rotation = rotation * t.rotation;
	res.scale = vector3(scale.x * t.scale.x, scale.y * t.scale.y, scale.z * t.scale.z);
	return res;
}
template<typename T>
vector3<T> transform3d<T>::operator*(vector3<T> v) const
{
	return rotation * vector3(scale.x * v.x, scale.y * v.y, scale.z * v.z) + position;
}
template<typename T>
quaternion<T> transform3d<T>::operator*(quaternion<T> q) const { return rotation * q; }
#pragma endregion

#pragma region definition uniform_transform2d
template<typename T>
uniform_transform2d<T> uniform_transform2d<T>::operator*(uniform_transform2d<T> t) const
{
	uniform_transform2d<T> res;
	matrix2 r = matrix2<T>::rotation(rotation);
	res.position = r * (scale * t.position) + position;
	res.rotation = rotation + t.rotation;
	res.scale = scale * t.scale;
	return res;
}
template<typename T>
vector2<T> uniform_transform2d<T>::operator*(vector2<T> v) const
{
	matrix2 r = matrix2::rotation(rotation);
	return r * (scale * v) + position;
}
#pragma endregion

#pragma region definition uniform_transform3d
template<typename T>
uniform_transform3d<T> uniform_transform3d<T>::operator*(const uniform_transform3d<T>& t) const
{
	uniform_transform3d<T> res;
	res.position = rotation * (scale * t.position) + position;
	res.rotation = rotation * t.rotation;
	res.scale = scale * t.scale;
	return res;
}
template<typename T>
vector3<T> uniform_transform3d<T>::operator*(vector3<T> v) const { return rotation * (scale * v) + position; }
template<typename T>
quaternion<T> uniform_transform3d<T>::operator*(quaternion<T> q) const { return rotation * q; }
#pragma endregion

#pragma region definition unscaled_transform2d
template<typename T>
unscaled_transform2d<T> unscaled_transform2d<T>::operator*(unscaled_transform2d<T> t) const
{
	unscaled_transform2d<T> res;
	matrix2 r = matrix2<T>::rotation(rotation);
	res.position = r * t.position + position;
	res.rotation = rotation + t.rotation;
	return res;
}
template<typename T>
vector2<T> unscaled_transform2d<T>::operator*(vector2<T> v) const
{
	matrix2<T> r = matrix2::rotation(rotation);
	return r * v + position;
}
#pragma endregion

#pragma region definition unscaled_transform3d
template<typename T>
unscaled_transform3d<T>& unscaled_transform3d<T>::translate_absolute(vector3<T> t)
{
	position += t;
	return *this;
}

template<typename T>
unscaled_transform3d<T>& unscaled_transform3d<T>::translate_relative(vector3<T> t)
{
	position += rotation * t;
	return *this;
}

template<typename T>
unscaled_transform3d<T>& unscaled_transform3d<T>::rotate(quaternion<T> r)
{
	rotation = r * rotation;
	return *this;
}

template<typename T>
unscaled_transform3d<T> inverse(unscaled_transform3d<T> t)
{
	t.rotation = inverse(t.rotation);
	t.position = t.rotation * -t.position;
	return t;
}

template<typename T>
matrix4<T> unscaled_transform3d<T>::to_matrix() const
{
	return matrix4<T>::from(*this);
}

template<typename T>
unscaled_transform3d<T> unscaled_transform3d<T>::operator*(const unscaled_transform3d<T>& t) const
{
	unscaled_transform3d<T> res;
	res.position = rotation * t.position + position;
	res.rotation = rotation * t.rotation;
	return res;
}
template<typename T>
vector3<T> unscaled_transform3d<T>::operator*(vector3<T> v) const { return rotation * v + position; }
template<typename T>
quaternion<T> unscaled_transform3d<T>::operator*(quaternion<T> q) const { return rotation * q; }
#pragma endregion

#pragma region definition streams
template<typename T>
std::ostream& operator<<(std::ostream& os, vector2<T> v)
{
	return os << "vector2(" << v.x << ", " << v.y << ")";
}
template<typename T>
std::ostream& operator<<(std::ostream& os, vector3<T> v)
{
	return os << "vector3(" << v.x << ", " << v.y << ", " << v.z << ")";
}
template<typename T>
std::ostream& operator<<(std::ostream& os, vector4<T> v)
{
	return os << "vector4(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")";
}
template<typename T>
std::ostream& operator<<(std::ostream& os, quaternion<T> q)
{
	return os << "quaternion({" << q.x << ", " << q.y << ", " << q.z << "}, " << q.w << ")";
}
template<typename T>
std::ostream& operator<<(std::ostream& os, plane<T> p)
{
	return os << "plane(origin = " << p.origin << ", normal = " << p.normal << ")";
}
template<typename T>
std::ostream& operator<<(std::ostream& os, const cdm::matrix2<T>& m)
{
	return os
		<<   "matrix2(" << m.m00 << "\t" << m.m10
		<< "\n        " << m.m01 << "\t" << m.m11
		<< ")";
}
template<typename T>
std::ostream& operator<<(std::ostream& os, const cdm::matrix3<T>& m)
{
	return os
		<<   "matrix3(" << m.m00 << "\t" << m.m10 << "\t" << m.m20
		<< "\n        " << m.m01 << "\t" << m.m11 << "\t" << m.m21
		<< "\n        " << m.m02 << "\t" << m.m12 << "\t" << m.m22
		<< ")";
}
template<typename T>
std::ostream& operator<<(std::ostream& os, const cdm::matrix4<T>& m)
{
	return os
		<<   "matrix4(" << m.m00 << "\t" << m.m10 << "\t" << m.m20 << "\t" << m.m30
		<< "\n        " << m.m01 << "\t" << m.m11 << "\t" << m.m21 << "\t" << m.m31
		<< "\n        " << m.m02 << "\t" << m.m12 << "\t" << m.m22 << "\t" << m.m32
		<< "\n        " << m.m03 << "\t" << m.m13 << "\t" << m.m23 << "\t" << m.m33
		<< ")";
}
#pragma endregion

#pragma region definition misc
namespace detail
{
template<unsigned char x, unsigned char y, typename T>
T get_quaternion_matrix_element(quaternion<T> q)
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
} // namespace detail

template<typename Functor, typename T>
std::vector<vector3<T>> function2D_sampler(const Functor& functor, T min, T max, T step)
{
	std::vector<vector3<T>> res;

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

template<typename Functor, typename T>
std::vector<std::vector<vector3<T>>> function3D_sampler(const Functor& functor, vector2<T> min, vector2<T> max, vector2<T> step)
{
	std::vector<std::vector<vector3<T>>> res;

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
		std::vector<vector3<T>> row;
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
cdm::radian<float> operator""_rad(long double d) { return cdm::radian<float>(static_cast<float>(d)); }
cdm::radian<float> operator""_rad(unsigned long long int i) { return cdm::radian<float>(static_cast<float>(i)); }
cdm::radian<float> operator""_pi(long double d) { return cdm::radian<float>(static_cast<float>(d) * cdm::pi); }
cdm::radian<float> operator""_pi(unsigned long long int i) { return cdm::radian<float>(static_cast<float>(i) * cdm::pi); }
cdm::degree<float> operator""_deg(long double d) { return cdm::degree<float>(static_cast<float>(d)); }
cdm::degree<float> operator""_deg(unsigned long long int i) { return cdm::degree<float>(static_cast<float>(i)); }
cdm::radian<double> operator""_radd(long double d) { return cdm::radian<double>(static_cast<double>(d)); }
cdm::radian<double> operator""_radd(unsigned long long int i) { return cdm::radian<double>(static_cast<double>(i)); }
cdm::radian<double> operator""_pid(long double d) { return cdm::radian<double>(static_cast<double>(d) * cdm::pi); }
cdm::radian<double> operator""_pid(unsigned long long int i) { return cdm::radian<double>(static_cast<double>(i) * cdm::pi); }
cdm::degree<double> operator""_degd(long double d) { return cdm::degree<double>(static_cast<double>(d)); }
cdm::degree<double> operator""_degd(unsigned long long int i) { return cdm::degree<double>(static_cast<double>(i)); }
#pragma endregion
} // namespace literals
} // namespace cdm

#endif // CDM_MATHS_HPP
