/* cdm_maths - v0.2.1 - geometric library - https://github.com/WubiCookie/cdm
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
template<typename T>
class normalized;
struct complex;
struct radian;
struct degree;
struct vector2;
struct vector3;
struct vector4;
struct matrix2;
struct matrix3;
struct matrix3x4;
struct matrix4x3;
struct matrix4;
struct perspective;
struct euler_angles;
struct quaternion;
struct cartesian_direction2d;
struct polar_direction2d;
struct line;
struct plane;
//struct rect;
struct aa_rect;
struct circle;
//struct ellipse;
//struct aa_ellipse;
struct ray2d;
struct ray3d;
struct aabb;
struct transform2d;
struct transform3d;
struct uniform_transform2d;
struct uniform_transform3d;
struct unscaled_transform2d;
struct unscaled_transform3d;

namespace detail
{
template<unsigned char x, unsigned char y>
float get_quaternion_matrix_element(quaternion q);
} // namespace detail

inline constexpr float pi = 3.141592653589793238462643f;
inline constexpr float deg_to_rad = pi / 180.0f;
inline constexpr float rad_to_deg = 180.0f / pi;

template<typename T>
T clamp(T f, T min_, T max_) { return std::min(max_, std::max(min_, f)); }

template<typename T, typename U>
T lerp(T begin, T end, U percent) { return (end - begin) * percent + begin; }

template<typename T>
bool nearly_equal_epsilon(T f1, T f2, T epsilon) { return std::abs(f1 - f2) < epsilon; }

template<typename T>
bool nearly_equal(T f1, T f2) { return nearly_equal_epsilon(f1, f2, std::numeric_limits<float>::epsilon()); }

template<typename T>
constexpr int sign(T val);
template<typename T>
constexpr int signnum(T val);

struct complex
{
	float r = 0.0f;
	float i = 0.0f;

	complex() = default;
	complex(const complex&) = default;
	complex(complex&&) = default;
	complex(float real, float imaginary);
	explicit complex(radian angle);

	float norm() const;
	float norm_squared() const;
	complex& normalize();
	complex get_normalized() const;
	complex& conjugate();
	complex get_conjugated() const;

	complex operator+(complex c) const;
	complex operator-(complex c) const;

	complex& operator+=(complex c);
	complex& operator-=(complex c);
	complex& operator*=(complex c);

	complex& operator=(const complex&) = default;
	complex& operator=(complex&&) = default;
};

complex operator*(complex c1, complex c2);
complex operator*(complex c, float f);
complex operator*(normalized<complex> c1, complex c2);
complex operator*(complex c1, normalized<complex> c2);
normalized<complex> operator*(normalized<complex> c1, normalized<complex> c2);
vector2 operator*(normalized<complex> c, vector2 v);

struct radian
{
	float angle = 0.0f;

	radian() = default;
	explicit radian(float f);
	radian(const radian& r);
	radian(degree d);

	explicit operator float() const;

	float sin() const;
	float cos() const;
	float tan() const;
	float asin() const;
	float acos() const;
	float atan() const;
	float sinh() const;
	float cosh() const;
	float tanh() const;
	float asinh() const;
	float acosh() const;
	float atanh() const;

	radian& operator+=(float f);
	radian& operator+=(radian r);
	radian& operator-=(float f);
	radian& operator-=(radian r);
	radian& operator*=(float f);
	radian& operator*=(radian r);
	radian& operator/=(float f);
	radian& operator/=(radian r);

	radian operator-() const;

	radian& operator=(radian r);
	radian& operator=(degree d);
};

radian operator+(radian, radian);
radian operator-(radian, radian);
radian operator*(radian, radian);
radian operator/(radian, radian);
radian operator*(radian, float);
radian operator/(radian, float);
radian operator*(float, radian);
radian operator/(float, radian);
bool operator<(float, radian);
bool operator>(float, radian);
bool operator==(float, radian);
bool operator!=(float, radian);
bool operator>=(float, radian);
bool operator<=(float, radian);
bool operator<(radian, float);
bool operator>(radian, float);
bool operator==(radian, float);
bool operator!=(radian, float);
bool operator>=(radian, float);
bool operator<=(radian, float);
bool operator<(radian, radian);
bool operator>(radian, radian);
bool operator==(radian, radian);
bool operator!=(radian, radian);
bool operator>=(radian, radian);
bool operator<=(radian, radian);

struct degree
{
	float angle = 0.0f;

	degree() = default;
	explicit degree(float f);
	degree(const degree& d);
	degree(radian r);

	explicit operator float() const;

	float sin() const;
	float cos() const;
	float tan() const;
	float asin() const;
	float acos() const;
	float atan() const;
	float sinh() const;
	float cosh() const;
	float tanh() const;
	float asinh() const;
	float acosh() const;
	float atanh() const;

	degree& operator+=(float f);
	degree& operator+=(degree d);
	degree& operator-=(float f);
	degree& operator-=(degree d);
	degree& operator*=(float f);
	degree& operator*=(degree d);
	degree& operator/=(float f);
	degree& operator/=(degree d);

	degree operator-() const;

	degree& operator=(degree d);
	degree& operator=(radian r);
};

degree operator+(degree, degree);
degree operator-(degree, degree);
degree operator*(degree, degree);
degree operator/(degree, degree);
degree operator*(degree, float);
degree operator/(degree, float);
degree operator*(float, degree);
degree operator/(float, degree);
bool operator<(float, degree);
bool operator>(float, degree);
bool operator==(float, degree);
bool operator!=(float, degree);
bool operator>=(float, degree);
bool operator<=(float, degree);
bool operator<(degree, float);
bool operator>(degree, float);
bool operator==(degree, float);
bool operator!=(degree, float);
bool operator>=(degree, float);
bool operator<=(degree, float);
bool operator<(degree, degree);
bool operator>(degree, degree);
bool operator==(degree, degree);
bool operator!=(degree, degree);
bool operator>=(degree, degree);
bool operator<=(degree, degree);

struct vector2
{
	float x = 0.0f, y = 0.0f;

	vector2() = default;
	vector2(float x, float y);

	float norm() const;
	float norm_squared() const;
	vector2& normalize();
	vector2 get_normalized() const;
	vector2& clamp(vector2 min, vector2 max);
	vector2 get_clamped(vector2 min, vector2 max) const;
	vector2& negate();
	vector2 get_negated() const;
	float dot(vector2 v) const;
	static vector2 lerp(vector2 begin, vector2 end, float percent);
	static vector2 nlerp(vector2 begin, vector2 end, float percent);
	static vector2 slerp(vector2 begin, vector2 end, float percent);
	float distance_from(vector2 v) const;
	static float distance_between(vector2 v1, vector2 v2);
	static vector2 from_to(vector2 from, vector2 to);
	vector2 to(vector2 v) const;
	static radian angle(vector2 v1, vector2 v2);
	static bool nearly_equal(vector2 v1, vector2 v2);

	bool lies_on(const line& l);
	//bool lies_on(const rect& r);
	bool lies_on(const aa_rect& r);
	bool lies_on(const circle& c);
	//bool lies_on(const ellipe& e);
	//bool lies_on(const aa_ellipe& e);

	float& operator[](size_t i);
	float operator[](size_t i) const;

	vector2 operator+(vector2 v) const;
	vector2 operator-(vector2 v) const;
	vector2 operator*(float f) const;
	vector2 operator/(float f) const;

	vector2& operator+=(vector2 v);
	vector2& operator-=(vector2 v);
	vector2& operator*=(float f);
	vector2& operator/=(float f);

	vector2 operator-() const;

	bool operator==(vector2 v) const;
	bool operator!=(vector2 v) const;
};

vector2 operator*(float f, vector2 v);
float dot(vector2 lhs, vector2 rhs);

struct vector3
{
	float x = 0.0f, y = 0.0f, z = 0.0f;

	vector3() = default;
	vector3(float x, float y, float z);
	vector3(vector2 v, float z);

	float norm() const;
	float norm_squared() const;
	vector3& normalize();
	vector3 get_normalized() const;
	vector3& clamp(vector3 min, vector3 max);
	vector3 get_clamped(vector3 min, vector3 max) const;
	vector3& negate();
	vector3 get_negated() const;
	float dot(vector3 v) const;
	vector3 cross(vector3 v) const;
	static vector3 lerp(vector3 begin, vector3 end, float percent);
	static vector3 nlerp(vector3 begin, vector3 end, float percent);
	float distance_from(vector3 v) const;
	static float distance_between(vector3 v1, vector3 v2);
	static float distance_squared_between(vector3 v1, vector3 v2);
	static vector3 from_to(vector3 from, vector3 to);
	vector3 to(vector3 v) const;
	radian angle(vector3 v) const;
	radian angle_around_axis(vector3 v, vector3 axis);
	static bool nearly_equal(vector3 v1, vector3 v2, float epsilon = std::numeric_limits<float>::epsilon());

	vector2 xy() const;

	float& operator[](size_t i);
	float operator[](size_t i) const;

	vector3 operator+(vector3 v) const;
	vector3 operator-(vector3 v) const;
	vector3 operator*(float f) const;
	vector3 operator/(float f) const;

	vector3& operator+=(vector3 v);
	vector3& operator-=(vector3 v);
	vector3& operator*=(float f);
	vector3& operator/=(float f);

	vector3 operator-() const;

	bool operator==(vector3 v) const;
	bool operator!=(vector3 v) const;
};

vector3 operator*(float f, vector3 v);
float dot(vector3 lhs, vector3 rhs);

struct vector4
{
	float x = 0.0f, y = 0.0f, z = 0.0f, w = 0.0f;

	vector4() = default;
	vector4(float x, float y, float z, float w);
	vector4(vector2 v, float z, float w);
	vector4(vector3 v, float w);

	float norm() const;
	float norm_squared() const;
	vector4& normalize();
	vector4 get_normalized() const;
	vector4& clamp(vector4 min, vector4 max);
	vector4 get_clamped(vector4 min, vector4 max) const;
	vector4& negate();
	vector4 get_negated() const;
	float dot(vector4 v) const;
	static vector4 lerp(vector4 begin, vector4 end, float percent);
	static vector4 nlerp(vector4 begin, vector4 end, float percent);
	float distance_from(vector4 v) const;
	static float distance_between(vector4 v1, vector4 v2);
	static vector4 from_to(vector4 from, vector4 to);
	vector4 to(vector4 v) const;
	static radian angle(vector4 v1, vector4 v2);

	vector2 xy() const;
	vector3 xyz() const;

	vector4 operator+(vector4 v) const;
	vector4 operator-(vector4 v) const;
	vector4 operator*(float f) const;
	vector4 operator/(float f) const;

	vector4& operator+=(vector4 v);
	vector4& operator-=(vector4 v);
	vector4& operator*=(float f);
	vector4& operator/=(float f);

	vector4 operator-() const;

	bool operator==(vector4 v) const;
	bool operator!=(vector4 v) const;
};

float dot(vector4 lhs, vector4 rhs);

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
	T operator*(float f) const;
	T operator/(float f) const;

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

struct matrix2
{
	float
		m00 = 0.0f, m10 = 0.0f,
		m01 = 0.0f, m11 = 0.0f;

	matrix2() = default;
	matrix2(float m00, float m10, float m01, float m11);
	explicit matrix2(const matrix3& m);
	explicit matrix2(const matrix4& m);

	static matrix2 zero(void);
	static matrix2 identity(void);
	static matrix2 rotation(radian angle);
	static matrix2 rotation(normalized<complex> angle);
	matrix2& transpose();
	matrix2 get_transposed() const;
	float get_determinant() const;

	float& at(uint8_t x, uint8_t y);
	const float& at(uint8_t x, uint8_t y) const;

	vector2 operator*(vector2 v) const;
	matrix2 operator*(matrix2 m) const;
};

struct matrix3
{
	float
		m00 = 0.0f, m10 = 0.0f, m20 = 0.0f,
		m01 = 0.0f, m11 = 0.0f, m21 = 0.0f,
		m02 = 0.0f, m12 = 0.0f, m22 = 0.0f;

	matrix3() = default;
	matrix3(float m00, float m10, float m20,
		float m01, float m11, float m21,
		float m02, float m12, float m22);

	explicit matrix3(matrix2 m);
	explicit matrix3(const matrix4& m);
	explicit matrix3(const transform2d& t);
	explicit matrix3(uniform_transform2d t);
	explicit matrix3(unscaled_transform2d t);

	static matrix3 zero();
	static matrix3 identity();
	static matrix3 rotation(euler_angles r);
	static matrix3 rotation(quaternion q);
	static matrix3 rotation_around_x(radian angle);
	static matrix3 rotation_around_y(radian angle);
	static matrix3 rotation_around_z(radian angle);
	static matrix3 rotation_around_x(normalized<complex> angle);
	static matrix3 rotation_around_y(normalized<complex> angle);
	static matrix3 rotation_around_z(normalized<complex> angle);
	matrix3& inverse();
	matrix3 get_inversed() const;
	matrix3& transpose();
	matrix3 get_transposed() const;
	float get_determinant() const;
	bool is_orthogonal() const;

	float& at(uint8_t x, uint8_t y);
	const float& at(uint8_t x, uint8_t y) const;

	matrix3 operator*(float f) const;
	vector3 operator*(vector3 v) const;
	matrix3 operator*(const matrix3& m) const;
	matrix4x3 operator*(const matrix4x3& m) const;
};

struct matrix3x4
{
	float
		m00 = 0.0f, m10 = 0.0f, m20 = 0.0f,
		m01 = 0.0f, m11 = 0.0f, m21 = 0.0f,
		m02 = 0.0f, m12 = 0.0f, m22 = 0.0f,
		m03 = 0.0f, m13 = 0.0f, m23 = 0.0f;

	matrix3x4() = default;
	matrix3x4(float m00, float m10, float m20,
		float m01, float m11, float m21,
		float m02, float m12, float m22,
		float m03, float m13, float m23);

	explicit matrix3x4(matrix2 m);
	explicit matrix3x4(const matrix3& m);
	explicit matrix3x4(const matrix4& m);
	explicit matrix3x4(const transform3d& t);
	explicit matrix3x4(const uniform_transform3d& t);
	explicit matrix3x4(const unscaled_transform3d& t);

	static matrix3x4 zero();
	static matrix3x4 identity();
	static matrix3x4 rotation(euler_angles r);
	static matrix3x4 rotation(quaternion q);
	static matrix3x4 rotation_around_x(radian angle);
	static matrix3x4 rotation_around_y(radian angle);
	static matrix3x4 rotation_around_z(radian angle);
	static matrix3x4 rotation_around_x(normalized<complex> angle);
	static matrix3x4 rotation_around_y(normalized<complex> angle);
	static matrix3x4 rotation_around_z(normalized<complex> angle);
	matrix4x3 get_transposed() const;

	float& at(uint8_t x, uint8_t y);
	const float& at(uint8_t x, uint8_t y) const;

	matrix3x4 operator*(float f) const;
	vector4 operator*(vector3 v) const;
	matrix3x4 operator*(const matrix3& m) const;
	matrix4 operator*(const matrix4x3& m) const;
};

struct matrix4x3
{
	float
		m00 = 0.0f, m10 = 0.0f, m20 = 0.0f, m30,
		m01 = 0.0f, m11 = 0.0f, m21 = 0.0f, m31,
		m02 = 0.0f, m12 = 0.0f, m22 = 0.0f, m32;

	matrix4x3() = default;
	matrix4x3(float m00, float m10, float m20, float m30,
		float m01, float m11, float m21, float m31,
		float m02, float m12, float m22, float m32);

	explicit matrix4x3(matrix2 m);
	explicit matrix4x3(const matrix3& m);
	explicit matrix4x3(const matrix4& m);
	explicit matrix4x3(const transform3d& t);
	explicit matrix4x3(const uniform_transform3d& t);
	explicit matrix4x3(const unscaled_transform3d& t);

	static matrix4x3 zero();
	static matrix4x3 identity();
	static matrix4x3 rotation(euler_angles r);
	static matrix4x3 rotation(quaternion q);
	static matrix4x3 translation(vector3 t);
	static matrix4x3 rotation_around_x(radian angle);
	static matrix4x3 rotation_around_y(radian angle);
	static matrix4x3 rotation_around_z(radian angle);
	static matrix4x3 rotation_around_x(normalized<complex> angle);
	static matrix4x3 rotation_around_y(normalized<complex> angle);
	static matrix4x3 rotation_around_z(normalized<complex> angle);
	matrix3x4 get_transposed() const;

	float& at(uint8_t x, uint8_t y);
	const float& at(uint8_t x, uint8_t y) const;

	matrix4x3 operator*(float f) const;
	vector3 operator*(vector4 v) const;
	matrix3 operator*(const matrix3x4& m) const;
	matrix4x3 operator*(const matrix4& m) const;
};

struct matrix4
{
	float
		m00 = 0.0f, m10 = 0.0f, m20 = 0.0f, m30 = 0.0f,
		m01 = 0.0f, m11 = 0.0f, m21 = 0.0f, m31 = 0.0f,
		m02 = 0.0f, m12 = 0.0f, m22 = 0.0f, m32 = 0.0f,
		m03 = 0.0f, m13 = 0.0f, m23 = 0.0f, m33 = 0.0f;

	matrix4() = default;
	matrix4(float m00, float m10, float m20, float m30,
		float m01, float m11, float m21, float m31,
		float m02, float m12, float m22, float m32,
		float m03, float m13, float m23, float m33);
	explicit matrix4(matrix2 m);
	explicit matrix4(const matrix3& m);
	matrix4(const matrix4& m) = default;
	explicit matrix4(const perspective& p);
	explicit matrix4(const transform3d& t);
	explicit matrix4(const uniform_transform3d& t);
	explicit matrix4(const unscaled_transform3d& t);
	explicit matrix4(const std::array<float, 16>& a);
	explicit matrix4(const std::array<double, 16>& a);

	static matrix4 zero();
	static matrix4 identity();
	static matrix4 rotation(euler_angles r);
	static matrix4 rotation(quaternion q);
	static matrix4 translation(vector3 t);
	static matrix4 translation(float x, float y, float z);
	static matrix4 scale(vector3 t);
	static matrix4 look_at(vector3 from, vector3 to, vector3 up = { 0.0f, 1.0f, 0.0f });
	static matrix4 orthographic(float left, float right, float top, float bottom, float near, float far);
	static matrix4 rotation_around_x(radian angle);
	static matrix4 rotation_around_y(radian angle);
	static matrix4 rotation_around_z(radian angle);
	static matrix4 rotation_around_x(normalized<complex> angle);
	static matrix4 rotation_around_y(normalized<complex> angle);
	static matrix4 rotation_around_z(normalized<complex> angle);
	bool is_orthogonal() const;
	bool is_homogenous() const;
	matrix4& inverse();
	matrix4 get_inversed() const;
	matrix4& transpose();
	matrix4 get_transposed() const;
	float get_determinant() const;

	float& at(uint8_t x, uint8_t y);
	const float& at(uint8_t x, uint8_t y) const;

	matrix4 operator*(float f) const;
	matrix4 operator/(float f) const;
	vector4 operator*(vector4 v) const;
	matrix3x4 operator*(const matrix3x4& m) const;
	matrix4 operator*(const matrix4& m) const;

	operator std::array<float, 16>() const;
	operator std::array<double, 16>() const;
};

struct perspective
{
private:
	radian m_angle = degree(90.0f);
	float m_ratio = 1.0f;
	float m_inv_ratio = 1.0f;
	float m_near = 0.01f;
	float m_far = 10000.0f;
	float m_invTanHalfFovy = 1.0f / (m_angle / 2.0f).tan();

public:
	perspective() = default;
	perspective(radian angle, float ratio, float near, float far);
	perspective(const perspective&) = default;
	perspective(perspective&&) = default;

	perspective& operator=(const perspective&) = default;
	perspective& operator=(perspective&&) = default;

	void set_angle(radian angle);
	radian get_angle() const;

	void set_ratio(float ratio);
	float get_ratio() const;

	void set_near(float near_plane);
	float get_near() const;

	void set_far(float far_plane);
	float get_far() const;

	matrix4 to_matrix4() const;

	friend matrix4 operator*(const matrix4& m, const perspective& p);
	friend matrix4 operator*(const perspective& p, const matrix4& m);

	friend matrix4 operator*(const unscaled_transform3d& t, const perspective& p);
	friend matrix4 operator*(const perspective& p, const unscaled_transform3d& t);
};

matrix4 operator*(const matrix4& m, const perspective& p);
matrix4 operator*(const perspective& p, const matrix4& m);

matrix4 operator*(const unscaled_transform3d& t, const perspective& p);
matrix4 operator*(const perspective& p, const unscaled_transform3d& t);

struct euler_angles
{
	radian x, y, z;

	euler_angles() = default;
	euler_angles(radian x, radian y, radian z);
};

struct quaternion
{
	float x = 0.0f, y = 0.0f, z = 0.0f, w = 1.0f;

	quaternion() = default;
	quaternion(float x, float y, float z, float w);
	quaternion(const normalized<vector3>& axis, radian angle);

	static quaternion zero();
	static quaternion identity();
	quaternion& inverse();
	quaternion get_inversed() const;
	quaternion& conjugate();
	quaternion get_conjugated() const;
	quaternion& normalize();
	quaternion get_normalized() const;
	quaternion& negate();
	quaternion get_negated() const;
	quaternion& clamp(quaternion min, quaternion max);
	quaternion get_clamped(quaternion min, quaternion max) const;
	float norm() const;
	float norm_squared() const;
	static quaternion lerp(quaternion begin, quaternion end, float percent);
	static quaternion nlerp(quaternion begin, quaternion end, float percent);
	static quaternion slerp(quaternion begin, quaternion end, float percent);
	float dot(quaternion) const;
	quaternion cross(quaternion) const;

	quaternion operator+(quaternion q) const;
	quaternion operator-(quaternion q) const;
	quaternion operator*(quaternion q) const; // TODO: implement operator* for normalized<> if result is normalized too
	quaternion operator*(float f) const;
	quaternion operator/(float f) const;

	quaternion& operator+=(quaternion q);
	quaternion& operator-=(quaternion q);
	quaternion& operator*=(quaternion q);
	quaternion& operator*=(float f);
	quaternion& operator/=(float f);

	quaternion operator-() const;

	bool operator==(quaternion v) const;
	bool operator!=(quaternion v) const;
};

vector3 operator*(const normalized<quaternion>& q, vector3 v);

struct cartesian_direction2d
{
	float x, y;

	cartesian_direction2d() = default;
	cartesian_direction2d(float x, float y);
	cartesian_direction2d(const normalized<vector2>& direction);
	cartesian_direction2d(polar_direction2d direction);
	cartesian_direction2d(radian angle);
	cartesian_direction2d(const cartesian_direction2d&) = default;
	cartesian_direction2d(cartesian_direction2d&&) noexcept = default;
	~cartesian_direction2d() = default;

	cartesian_direction2d& operator=(const cartesian_direction2d&) = default;
	cartesian_direction2d& operator=(cartesian_direction2d&&) noexcept = default;

	cartesian_direction2d& rotate(radian angle);
	cartesian_direction2d get_rotated(radian angle);
};

struct polar_direction2d
{
	radian angle;

	polar_direction2d() = default;
	polar_direction2d(float x, float y);
	polar_direction2d(const normalized<vector2>& direction);
	polar_direction2d(cartesian_direction2d direction);
	polar_direction2d(radian angle);
	polar_direction2d(const polar_direction2d&) = default;
	polar_direction2d(polar_direction2d&&) noexcept = default;
	~polar_direction2d() = default;

	polar_direction2d& operator=(const polar_direction2d&) = default;
	polar_direction2d& operator=(polar_direction2d&&) noexcept = default;

	polar_direction2d& rotate(radian angle);
	polar_direction2d get_rotated(radian angle);

	polar_direction2d& wrap();
	polar_direction2d get_wrapped();
};

struct line
{
	float coefficient;
	float offset;

	static bool is_parallel(line l1, line l2);
};

inline bool collides(line l1, line l2);
inline bool collides(line l1, line l2, vector2& intersection);
inline bool collides(line l, vector2 v);
inline bool collides(vector2 v, line l);

struct plane
{
	vector3 origin;
	normalized<vector3> normal;

	float evaluate(vector3 point) const;
	vector3 project3d(vector3 point) const;
	vector2 project2d(vector3 point, normalized<vector3> plane_tangent) const;
	vector3 unproject(vector2 point, normalized<vector3> plane_tangent) const;
};

struct ray2d
{
	vector2 origin;
	normalized<vector2> direction;
};

struct ray3d
{
	vector3 origin;
	normalized<vector3> direction;
};

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

struct aa_rect
{
	vector2 origin;
	vector2 dimention;

	bool contains(vector2 v) const;
};

inline std::size_t collides(const aa_rect& r1, const aa_rect& r2, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const aa_rect& r1, const rect& r2, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const rect& r1, const aa_rect& r2, vector2* intersection1, vector2* insersection2);
inline std::size_t collides(const aa_rect& r, const line& l, vector2* intersection1, vector2* insersection2);
inline std::size_t collides(const line& l, const aa_rect& r, vector2* intersection1, vector2* insersection2);
inline bool collides(const aa_rect& r, vector2 v);
inline bool collides(vector2 v, const aa_rect& r);

struct circle
{
	vector2 origin;
	float radius;
};

inline std::size_t collides(const circle& c1, const circle& c2, vector2* intersection1, vector2* insersection2);
inline std::size_t collides(const circle& c, const aa_rect& r, vector2* intersection1, vector2* insersection2);
inline std::size_t collides(const aa_rect& r, const circle& c, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const circle& c, const rect& r, vector2* intersection1, vector2* insersection2);
//std::size_t collides(const rect& r, const circle& c, vector2* intersection1, vector2* insersection2);
inline std::size_t collides(const circle& c, const line& l, vector2* intersection1, vector2* insersection2);
inline std::size_t collides(const line& l, const circle& c, vector2* intersection1, vector2* insersection2);
inline bool collides(const circle& c, vector2 v);
inline bool collides(vector2 v, const circle& c);

//struct ellipse
//{
//	vector2 origin;
//	float semi_minor = 1.0f;
//	float semi_major = 1.0f;
//	complex angle;
//};
//
//inline std::size_t collides(const ellipse& e1, const ellipse& e2, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const ellipse& e, const circle& c, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const circle& c, const ellipse& e, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const ellipse& e, const aa_rect& r, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const aa_rect& r, const ellipse& e, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const ellipse& e, const rect& r, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const rect& r, const ellipse& e, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const ellipse& e, const line& l, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const line& l, const ellipse& e, vector2* intersection1, vector2* insersection2);
//inline bool collides(const ellipse& e, vector2 v);
//inline bool collides(vector2 v, const ellipse& e);
//
//struct aa_ellipse
//{
//	vector2 origin;
//	float semi_minor = 1.0f;
//	float semi_major = 1.0f;
//};
//
//inline std::size_t collides(const aa_ellipse& e1, const aa_ellipse& e2, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const aa_ellipse& e1, const ellipse& e2, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const ellipse& e1, const aa_ellipse& e2, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const aa_ellipse& e, const circle& c, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const circle& c, const aa_ellipse& e, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const aa_ellipse& e, const aa_rect& r, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const aa_rect& r, const aa_ellipse& e, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const aa_ellipse& e, const rect& r, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const rect& r, const aa_ellipse& e, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const aa_ellipse& e, const line& l, vector2* intersection1, vector2* insersection2);
//inline std::size_t collides(const line& l, const aa_ellipse& e, vector2* intersection1, vector2* insersection2);
//inline bool collides(const aa_ellipse& e, vector2 v);
//inline bool collides(vector2 v, const aa_ellipse& e);

inline bool collides(const ray3d& r, const plane& p);
inline bool collides(const ray3d& r, const plane& p, vector3& intersection);
inline bool collides(const plane& p, const ray3d& r);
inline bool collides(const plane& p, const ray3d& r, vector3& intersection);
inline bool collides_bidirectional(const plane& p, const ray3d& r, vector3& intersection);

struct aabb
{
	vector3 min;
	vector3 max;

	bool contains(vector3 p) const;
	vector3 get_center() const;

	aabb operator+(aabb rhs) const;
};

inline bool collides(aabb b, ray3d r);
inline bool collides(const aabb& b1, const aabb& b2);

struct transform2d
{
	vector2 position;
	normalized<complex> rotation;
	vector2 scale = {1.0f, 1.0f};

	transform2d() = default;
	transform2d(vector2 position, normalized<complex> rotation, vector2 scale);

	transform2d operator*(const transform2d& t) const;
	vector2 operator*(vector2 v) const;
};

struct transform3d
{
	vector3 position;
	quaternion rotation;
	vector3 scale = {1.0f, 1.0f, 1.0f};

	transform3d() = default;
	transform3d(vector3 position, quaternion rotation, vector3 scale);

	transform3d operator*(const transform3d& t) const;
	vector3 operator*(vector3 v) const;
	quaternion operator*(quaternion q) const;
};

struct uniform_transform2d
{
	vector2 position;
	radian rotation;
	float scale = 1.0f;

	uniform_transform2d() = default;
	uniform_transform2d(vector2 position, radian rotation, float scale);

	uniform_transform2d operator*(uniform_transform2d t) const;
	vector2 operator*(vector2 v) const;
};

struct uniform_transform3d
{
	vector3 position;
	quaternion rotation;
	float scale = 1.0f;

	uniform_transform3d() = default;
	uniform_transform3d(vector3 position, quaternion rotation, float scale);

	uniform_transform3d operator*(const uniform_transform3d& t) const;
	vector3 operator*(vector3 v) const;
	quaternion operator*(quaternion q) const;
};

struct unscaled_transform2d
{
	vector2 position;
	radian rotation;

	unscaled_transform2d() = default;
	unscaled_transform2d(vector2 position, radian rotation);

	unscaled_transform2d operator*(unscaled_transform2d t) const;
	vector2 operator*(vector2 v) const;
};

struct unscaled_transform3d
{
	vector3 position;
	quaternion rotation;

	unscaled_transform3d() = default;
	unscaled_transform3d(vector3 position, quaternion rotation);

	unscaled_transform3d& translate_absolute(vector3 t);
	unscaled_transform3d& translate_relative(vector3 t);
	unscaled_transform3d& rotate(quaternion r);

	unscaled_transform3d& inverse();
	unscaled_transform3d get_inversed() const;

	matrix4 to_matrix() const;

	unscaled_transform3d operator*(const unscaled_transform3d& t) const;
	vector3 operator*(vector3 v) const;
	quaternion operator*(quaternion q) const;
};

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

inline complex::complex(float real, float imaginary) : r(real), i(imaginary) {}
inline complex::complex(radian angle) : r(angle.cos()), i(angle.sin()) {}

inline float complex::norm() const { return sqrtf(norm_squared()); }
inline float complex::norm_squared() const { return r * r + i * i; }
inline complex& complex::normalize()
{
	float n = norm();
	r /= n;
	i /= n;
	return *this;
}
inline complex complex::get_normalized() const
{
	complex res = *this;
	res.normalize();
	return res;
}
inline complex& complex::conjugate() { i = -i; return *this; }
inline complex complex::get_conjugated() const
{
	complex res = *this;
	res.conjugate();
	return res;
}

inline complex complex::operator+(complex c) const
{
	return {
		r + c.r,
		i + c.i
	};
}
inline complex complex::operator-(complex c) const
{
	return {
		r - c.r,
		i - c.i
	};
}

inline complex& complex::operator+=(complex c) { return *this = *this + c; }
inline complex& complex::operator-=(complex c) { return *this = *this - c; }
inline complex& complex::operator*=(complex c) { return *this = *this * c; }

inline complex operator*(complex c1, complex c2)
{
	return {
		c1.r * c2.r - c1.i * c2.i,
		c1.r * c2.i + c1.i * c2.r
	};
}
inline complex operator*(complex c, float f) { return {c.r * f, c.i * f}; }
inline complex operator*(normalized<complex> c1, complex c2) { return *c1 * c2; }
inline complex operator*(complex c1, normalized<complex> c2) { return c1 * *c2; }
inline normalized<complex> operator*(normalized<complex> c1, normalized<complex> c2)
{
	return normalized<complex>::already_normalized(complex(
		c1->r * c2->r - c1->i * c2->i,
		c1->r * c2->i + c1->i * c2->r
	));
}
inline vector2 operator*(normalized<complex> c, vector2 v)
{
	return {
		c->r * v.x - c->i * v.y,
		c->r * v.y + c->i * v.x
	};
}

inline radian::radian(float f) : angle(f) {}
inline radian::radian(const radian& r) : angle(r.angle) {}
inline radian::radian(degree d) : angle(d.angle * deg_to_rad) {}

inline radian::operator float() const { return angle; }

inline float radian::sin() const { return std::sinf(angle); }
inline float radian::cos() const { return std::cosf(angle); }
inline float radian::tan() const { return std::tanf(angle); }
inline float radian::asin() const { return std::asinf(angle); }
inline float radian::acos() const { return std::acosf(angle); }
inline float radian::atan() const { return std::atanf(angle); }
inline float radian::sinh() const { return std::sinhf(angle); }
inline float radian::cosh() const { return std::coshf(angle); }
inline float radian::tanh() const { return std::tanhf(angle); }
inline float radian::asinh() const { return std::asinhf(angle); }
inline float radian::acosh() const { return std::acoshf(angle); }
inline float radian::atanh() const { return std::atanhf(angle); }

inline radian& radian::operator+=(float f) { angle += f; return *this; }
inline radian& radian::operator+=(radian r) { angle += r.angle; return *this; }
inline radian& radian::operator-=(float f) { angle -= f; return *this; }
inline radian& radian::operator-=(radian r) { angle -= r.angle; return *this; }
inline radian& radian::operator*=(float f) { angle *= f; return *this; }
inline radian& radian::operator*=(radian r) { angle *= r.angle; return *this; }
inline radian& radian::operator/=(float f) { angle /= f; return *this; }
inline radian& radian::operator/=(radian r) { angle /= r.angle; return *this; }

inline radian radian::operator-() const { return radian(-angle); }

inline radian& radian::operator=(radian r) { angle = r.angle; return *this; }
inline radian& radian::operator=(degree d) { angle = d.angle * deg_to_rad; return *this; }

inline degree::degree(float f) : angle(f) {}
inline degree::degree(const degree& d) : angle(d.angle) {}
inline degree::degree(radian r) : angle(r.angle * rad_to_deg) {}

inline degree::operator float() const { return angle; }

inline float degree::sin() const { return radian(*this).sin(); }
inline float degree::cos() const { return radian(*this).cos(); }
inline float degree::tan() const { return radian(*this).tan(); }
inline float degree::asin() const { return radian(*this).asin(); }
inline float degree::acos() const { return radian(*this).acos(); }
inline float degree::atan() const { return radian(*this).atan(); }
inline float degree::sinh() const { return radian(*this).sinh(); }
inline float degree::cosh() const { return radian(*this).cosh(); }
inline float degree::tanh() const { return radian(*this).tanh(); }
inline float degree::asinh() const { return radian(*this).asinh(); }
inline float degree::acosh() const { return radian(*this).acosh(); }
inline float degree::atanh() const { return radian(*this).atanh(); }

inline degree& degree::operator+=(float f) { angle += f; return *this; }
inline degree& degree::operator+=(degree d) { angle += d.angle; return *this; }
inline degree& degree::operator-=(float f) { angle -= f; return *this; }
inline degree& degree::operator-=(degree d) { angle -= d.angle; return *this; }
inline degree& degree::operator*=(float f) { angle *= f; return *this; }
inline degree& degree::operator*=(degree d) { angle *= d.angle; return *this; }
inline degree& degree::operator/=(float f) { angle /= f; return *this; }
inline degree& degree::operator/=(degree d) { angle /= d.angle; return *this; }

inline degree degree::operator-() const { return radian(-angle); }

inline degree& degree::operator=(degree d) { angle = d.angle; return *this; }
inline degree& degree::operator=(radian r) { angle = r.angle * rad_to_deg; return *this; }

inline radian operator+(radian r1, radian r2) { return radian(r1.angle + r2.angle); }
inline radian operator-(radian r1, radian r2) { return radian(r1.angle - r2.angle); }
inline radian operator*(radian r1, radian r2) { return radian(r1.angle * r2.angle); }
inline radian operator/(radian r1, radian r2) { return radian(r1.angle / r2.angle); }
inline radian operator*(radian r, float f) { return radian(r.angle * f); }
inline radian operator/(radian r, float f) { return radian(r.angle / f); }
inline radian operator*(float f, radian r) { return radian(f * r.angle); }
inline radian operator/(float f, radian r) { return radian(f / r.angle); }
inline bool operator<(float lhs, radian rhs) { return float(lhs) < float(rhs); }
inline bool operator>(float lhs, radian rhs) { return float(lhs) > float(rhs); }
inline bool operator==(float lhs, radian rhs) { return float(lhs) == float(rhs); }
inline bool operator!=(float lhs, radian rhs) { return float(lhs) != float(rhs); }
inline bool operator>=(float lhs, radian rhs) { return float(lhs) >= float(rhs); }
inline bool operator<=(float lhs, radian rhs) { return float(lhs) <= float(rhs); }
inline bool operator<(radian lhs, float rhs) { return float(lhs) < float(rhs); }
inline bool operator>(radian lhs, float rhs) { return float(lhs) > float(rhs); }
inline bool operator==(radian lhs, float rhs) { return float(lhs) == float(rhs); }
inline bool operator!=(radian lhs, float rhs) { return float(lhs) != float(rhs); }
inline bool operator>=(radian lhs, float rhs) { return float(lhs) >= float(rhs); }
inline bool operator<=(radian lhs, float rhs) { return float(lhs) <= float(rhs); }
inline bool operator<(radian lhs, radian rhs) { return float(lhs) < float(rhs); }
inline bool operator>(radian lhs, radian rhs) { return float(lhs) > float(rhs); }
inline bool operator==(radian lhs, radian rhs) { return float(lhs) == float(rhs); }
inline bool operator!=(radian lhs, radian rhs) { return float(lhs) != float(rhs); }
inline bool operator>=(radian lhs, radian rhs) { return float(lhs) >= float(rhs); }
inline bool operator<=(radian lhs, radian rhs) { return float(lhs) <= float(rhs); }

inline degree operator+(degree d1, degree d2) { return degree(d1.angle + d2.angle); }
inline degree operator-(degree d1, degree d2) { return degree(d1.angle - d2.angle); }
inline degree operator*(degree d1, degree d2) { return degree(d1.angle * d2.angle); }
inline degree operator/(degree d1, degree d2) { return degree(d1.angle / d2.angle); }
inline degree operator*(degree d, float f) { return degree(d.angle * f); }
inline degree operator/(degree d, float f) { return degree(d.angle / f); }
inline degree operator*(float f, degree d) { return degree(f * d.angle); }
inline degree operator/(float f, degree d) { return degree(f / d.angle); }
inline bool operator<(float lhs, degree rhs) { return float(lhs) < float(rhs); }
inline bool operator>(float lhs, degree rhs) { return float(lhs) > float(rhs); }
inline bool operator==(float lhs, degree rhs) { return float(lhs) == float(rhs); }
inline bool operator!=(float lhs, degree rhs) { return float(lhs) != float(rhs); }
inline bool operator>=(float lhs, degree rhs) { return float(lhs) >= float(rhs); }
inline bool operator<=(float lhs, degree rhs) { return float(lhs) <= float(rhs); }
inline bool operator<(degree lhs, float rhs) { return float(lhs) < float(rhs); }
inline bool operator>(degree lhs, float rhs) { return float(lhs) > float(rhs); }
inline bool operator==(degree lhs, float rhs) { return float(lhs) == float(rhs); }
inline bool operator!=(degree lhs, float rhs) { return float(lhs) != float(rhs); }
inline bool operator>=(degree lhs, float rhs) { return float(lhs) >= float(rhs); }
inline bool operator<=(degree lhs, float rhs) { return float(lhs) <= float(rhs); }
inline bool operator<(degree lhs, degree rhs) { return float(lhs) < float(rhs); }
inline bool operator>(degree lhs, degree rhs) { return float(lhs) > float(rhs); }
inline bool operator==(degree lhs, degree rhs) { return float(lhs) == float(rhs); }
inline bool operator!=(degree lhs, degree rhs) { return float(lhs) != float(rhs); }
inline bool operator>=(degree lhs, degree rhs) { return float(lhs) >= float(rhs); }
inline bool operator<=(degree lhs, degree rhs) { return float(lhs) <= float(rhs); }

inline vector2::vector2(float x_, float y_) : x(x_), y(y_) {}

inline float vector2::norm() const { return sqrtf(norm_squared()); }
inline float vector2::norm_squared() const { return x * x + y * y; }
inline vector2& vector2::normalize()
{
	float n = norm();
	x /= n;
	y /= n;
	return *this;
}
inline vector2 vector2::get_normalized() const
{
	vector2 res = *this;
	res.normalize();
	return res;
}
inline vector2& vector2::clamp(vector2 min, vector2 max)
{
	x = cdm::clamp(x, min.x, max.x);
	y = cdm::clamp(y, min.y, max.y);
	return *this;
}
inline vector2 vector2::get_clamped(vector2 min, vector2 max) const
{
	vector2 res = *this;
	res.clamp(min, max);
	return res;
}
inline vector2& vector2::negate()
{
	x = -x;
	y = -y;
	return *this;
}
inline vector2 vector2::get_negated() const
{
	vector2 res = *this;
	res.negate();
	return res;
}
inline float vector2::dot(vector2 v) const { return x * v.x + y * v.y; }
inline vector2 vector2::lerp(vector2 begin, vector2 end, float percent) {	return (end - begin) * percent + begin; }
inline vector2 vector2::nlerp(vector2 begin, vector2 end, float percent) { return lerp(begin, end, percent).get_normalized(); }
inline vector2 vector2::slerp(vector2 begin, vector2 end, float percent)
{
	const radian angle = vector2::angle(begin, end) * percent;
	const float s = angle.sin();
	const float c = angle.cos();

	/*vector2 res{
		c * begin.x - s * begin.y,
		s * begin.x + c * begin.y
	};
	res.normalize();*/

	normalized<vector2> res {
		{
			c * begin.x - s * begin.y,
			s * begin.x + c * begin.y
		}
	};

	float f = cdm::lerp(begin.norm(), end.norm(), percent);
	return res * f;
}
inline float vector2::distance_from(vector2 v) const { return distance_between(*this, v); }
inline float vector2::distance_between(vector2 v1, vector2 v2) { return (v1 - v2).norm(); }
inline vector2 vector2::from_to(vector2 from, vector2 to) { return {to.x - from.x, to.y - from.y}; }
inline vector2 vector2::to(vector2 v) const { return from_to(*this, v); }
inline radian vector2::angle(vector2 v1, vector2 v2) { return radian(atan2f(v2.y, v2.x) - atan2f(v1.y, v1.x)); }
inline bool vector2::nearly_equal(vector2 v1, vector2 v2) { return cdm::nearly_equal(v1.x, v2.x) && cdm::nearly_equal(v1.y, v2.y); }

inline bool vector2::lies_on(const line& l) { return cdm::nearly_equal(l.coefficient * x + l.offset, y); }
//bool vector2::lies_on(const rect& r)
inline bool vector2::lies_on(const aa_rect& r)
{
	return ((x >= r.origin.x) && (x <= r.origin.x + r.dimention.x) && (lies_on(line{ 0, r.origin.y }) || lies_on(line{ 0, r.origin.y + r.dimention.y }))) ||
		((y >= r.origin.y) && (y <= r.origin.y + r.dimention.y) && (lies_on(line{ r.origin.x, 0 }) || lies_on(line{ r.origin.x + r.dimention.x, 0 })));
}
//bool vector2::lies_on(const circle& c)
//bool vector2::lies_on(const ellipe& e)
//bool vector2::lies_on(const aa_ellipe& e)

inline float& vector2::operator[](size_t i) { return reinterpret_cast<float*>(this)[i]; }
inline float vector2::operator[](size_t i) const { return reinterpret_cast<const float*>(this)[i]; }

inline vector2 vector2::operator+(vector2 v) const { return {x + v.x, y + v.y}; }
inline vector2 vector2::operator-(vector2 v) const { return {x - v.x, y - v.y}; }
inline vector2 vector2::operator*(float f) const { return {x * f, y * f}; }
inline vector2 vector2::operator/(float f) const { return {x / f, y / f}; }

inline vector2& vector2::operator+=(vector2 v) { *this = *this + v; return *this; }
inline vector2& vector2::operator-=(vector2 v) { *this = *this - v; return *this; }
inline vector2& vector2::operator*=(float f) { *this = *this * f; return *this; }
inline vector2& vector2::operator/=(float f) { *this = *this / f; return *this; }

inline vector2 vector2::operator-() const { return get_negated(); }

inline bool vector2::operator==(vector2 v) const { return x == v.x && y == v.y; }
inline bool vector2::operator!=(vector2 v) const { return !operator==(v); }

inline vector2 operator*(float f, vector2 v) { return v * f; }

inline float dot(vector2 lhs, vector2 rhs) { return lhs.dot(rhs); }

inline vector3::vector3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}
inline vector3::vector3(vector2 v, float z_) : vector3(v.x, v.y, z_) {}

inline float vector3::norm() const { return sqrtf(norm_squared()); }
inline float vector3::norm_squared() const { return x * x + y * y + z * z; }
inline vector3& vector3::normalize()
{
	float n = norm();
	x /= n;
	y /= n;
	z /= n;
	return *this;
}
inline vector3 vector3::get_normalized() const
{
	vector3 res = *this;
	res.normalize();
	return res;
}
inline vector3& vector3::clamp(vector3 min, vector3 max)
{
	x = cdm::clamp(x, min.x, max.x);
	y = cdm::clamp(y, min.y, max.y);
	z = cdm::clamp(z, min.z, max.z);
	return *this;
}
inline vector3 vector3::get_clamped(vector3 min, vector3 max) const
{
	vector3 res = *this;
	res.clamp(min, max);
	return res;
}
inline vector3& vector3::negate()
{
	x = -x;
	y = -y;
	z = -z;
	return *this;
}
inline vector3 vector3::get_negated() const
{
	vector3 res = *this;
	res.negate();
	return res;
}
inline float vector3::dot(vector3 v) const { return x * v.x + y * v.y + z * v.z; }
inline vector3 vector3::cross(vector3 v) const
{
	return {
		y * v.z - z * v.y,
		z * v.x - x * v.z,
		x * v.y - y * v.x
	};
}
inline vector3 vector3::lerp(vector3 begin, vector3 end, float percent) {	return (end - begin) * percent + begin; }
inline vector3 vector3::nlerp(vector3 begin, vector3 end, float percent) { return lerp(begin, end, percent).get_normalized(); }
inline float vector3::distance_from(vector3 v) const { return distance_between(*this, v); }
inline float vector3::distance_between(vector3 v1, vector3 v2) { return (v1 - v2).norm(); }
inline float vector3::distance_squared_between(vector3 v1, vector3 v2) { return (v1 - v2).norm_squared(); }
inline vector3 vector3::from_to(vector3 from, vector3 to) { return {to.x - from.x, to.y - from.y, to.z - from.z}; }
inline vector3 vector3::to(vector3 v) const { return from_to(*this, v); }
inline radian vector3::angle(vector3 v) const
{
	float divisor = sqrtf(norm_squared() * v.norm_squared());
	float alpha = dot(v) / divisor;
	return radian(acosf(cdm::clamp(alpha, -1.0f, 1.0f)));
}
inline radian vector3::angle_around_axis(vector3 v, vector3 axis)
{
	vector3 c = cross(v);
	float angle = atan2f(c.norm(), dot(v));
	return radian(c.dot(axis) < 0.0f ? -angle : angle);
}
inline bool vector3::nearly_equal(vector3 v1, vector3 v2, float epsilon)
{
	return nearly_equal_epsilon(v1.x, v2.x, epsilon) &&
	       nearly_equal_epsilon(v1.y, v2.y, epsilon) &&
	       nearly_equal_epsilon(v1.z, v2.z, epsilon);
}

inline vector2 vector3::xy() const { return {x, y}; }

inline float& vector3::operator[](size_t i) { return reinterpret_cast<float*>(this)[i]; }
inline float vector3::operator[](size_t i) const { return reinterpret_cast<const float*>(this)[i]; }

inline vector3 vector3::operator+(vector3 v) const { return {x + v.x, y + v.y, z + v.z}; }
inline vector3 vector3::operator-(vector3 v) const { return {x - v.x, y - v.y, z - v.z}; }
inline vector3 vector3::operator*(float f) const { return {x * f, y * f, z * f}; }
inline vector3 vector3::operator/(float f) const { return {x / f, y / f, z / f}; }

inline vector3& vector3::operator+=(vector3 v) { *this = *this + v; return *this; }
inline vector3& vector3::operator-=(vector3 v) { *this = *this - v; return *this; }
inline vector3& vector3::operator*=(float f) { *this = *this * f; return *this; }
inline vector3& vector3::operator/=(float f) { *this = *this / f; return *this; }

inline vector3 vector3::operator-() const { return get_negated(); }

inline bool vector3::operator==(vector3 v) const { return x == v.x && y == v.y && z == v.z; }
inline bool vector3::operator!=(vector3 v) const { return !operator==(v); }

inline vector3 operator*(float f, vector3 v) { return v * f; }

inline float dot(vector3 lhs, vector3 rhs) { return lhs.dot(rhs); }

inline vector4::vector4(float x_, float y_, float z_, float w_) : x(x_), y(y_), z(z_), w(w_) {}
inline vector4::vector4(vector2 v, float z_, float w_) : vector4(v.x, v.y, z_, w_) {}
inline vector4::vector4(vector3 v, float w_) : vector4(v.x, v.y, v.z, w_) {}

inline float vector4::norm() const { return sqrtf(norm_squared()); }
inline float vector4::norm_squared() const { return x * x + y * y + z * z + w * w; }
inline vector4& vector4::normalize()
{
	float n = norm();
	x /= n;
	y /= n;
	z /= n;
	w /= n;
	return *this;
}
inline vector4 vector4::get_normalized() const
{
	vector4 res = *this;
	res.normalize();
	return res;
}
inline vector4& vector4::clamp(vector4 min, vector4 max)
{
	x = cdm::clamp(x, min.x, max.x);
	y = cdm::clamp(y, min.y, max.y);
	z = cdm::clamp(z, min.z, max.z);
	w = cdm::clamp(w, min.w, max.w);
	return *this;
}
inline vector4 vector4::get_clamped(vector4 min, vector4 max) const
{
	vector4 res = *this;
	res.clamp(min, max);
	return res;
}
inline vector4& vector4::negate()
{
	x = -x;
	y = -y;
	z = -z;
	w = -w;
	return *this;
}
inline vector4 vector4::get_negated() const
{
	vector4 res = *this;
	res.negate();
	return res;
}
inline float vector4::dot(vector4 v) const { return x * v.x + y * v.y + z * v.z + w * v.w; }
inline vector4 vector4::lerp(vector4 begin, vector4 end, float percent) {	return (end - begin) * percent + begin; }
inline vector4 vector4::nlerp(vector4 begin, vector4 end, float percent) { return lerp(begin, end, percent).get_normalized(); }
inline float vector4::distance_from(vector4 v) const { return distance_between(*this, v); }
inline float vector4::distance_between(vector4 v1, vector4 v2) { return (v1 - v2).norm(); }
inline vector4 vector4::from_to(vector4 from, vector4 to) { return {to.x - from.x, to.y - from.y, to.z - from.z, to.w - from.w}; }
inline vector4 vector4::to(vector4 v) const { return from_to(*this, v); }

inline vector2 vector4::xy() const { return { x, y }; }
inline vector3 vector4::xyz() const { return { x, y, z }; }

inline vector4 vector4::operator+(vector4 v) const { return {x + v.x, y + v.y, z + v.z, w + v.w}; }
inline vector4 vector4::operator-(vector4 v) const { return {x - v.x, y - v.y, z - v.z, w - v.w}; }
inline vector4 vector4::operator*(float f) const { return {x * f, y * f, z * f, w * f}; }
inline vector4 vector4::operator/(float f) const { return {x / f, y / f, z / f, w / f}; }

inline vector4& vector4::operator+=(vector4 v) { *this = *this + v; return *this; }
inline vector4& vector4::operator-=(vector4 v) { *this = *this - v; return *this; }
inline vector4& vector4::operator*=(float f) { *this = *this * f; return *this; }
inline vector4& vector4::operator/=(float f) { *this = *this / f; return *this; }

inline vector4 vector4::operator-() const { return get_negated(); }

inline bool vector4::operator==(vector4 v) const { return x == v.x && y == v.y && z == v.z && w == v.w; }
inline bool vector4::operator!=(vector4 v) const { return !operator==(v); }

inline float dot(vector4 lhs, vector4 rhs) { return lhs.dot(rhs); }

template<typename T>
normalized<T>::normalized(const T& t) : vector(t) { vector.normalize(); }
template<typename T>
normalized<T>::normalized(T&& t) noexcept : vector(std::move(t)) { vector.normalize(); }

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
normalized<T> normalized<T>::operator+(const T& v) const { return (vector + v).get_normalized(); }
template<typename T>
normalized<T> normalized<T>::operator-(const T& v) const { return (vector - v).get_normalized(); }
template<typename T>
T normalized<T>::operator*(float f) const { return vector * f; }
template<typename T>
T normalized<T>::operator/(float f) const { return vector / f; }

template<typename T>
normalized<T>& normalized<T>::operator+=(const T& v) { vector = (vector + v).get_normalized(); return *this; }
template<typename T>
normalized<T>& normalized<T>::operator-=(const T& v) { vector = (vector - v).get_normalized(); return *this; }

template<typename T>
normalized<T> normalized<T>::operator-() const { return vector.get_negated(); }

template<typename T>
bool normalized<T>::operator==(const T& v) const { return vector == v; }
template<typename T>
bool normalized<T>::operator!=(const T& v) const { return vector != v; }

template<typename T>
normalized<T>& normalized<T>::operator=(const T& t) { vector = t.get_normalized(); return *this; }
template<typename T>
normalized<T>& normalized<T>::operator=(T&& t) noexcept { vector = t.get_normalized(); return *this; }

inline matrix2::matrix2(float m00_, float m10_,
	float m01_, float m11_) :
	m00(m00_), m10(m10_),
	m01(m01_), m11(m11_)
{}
inline matrix2::matrix2(const matrix3& m) :
	matrix2(
		m.m00, m.m10,
		m.m01, m.m11
	)
{}
inline matrix2::matrix2(const matrix4& m) :
	matrix2(
		m.m00, m.m10,
		m.m01, m.m11
	)
{}

inline matrix2 matrix2::zero() { return {0.0f, 0.0f, 0.0f, 0.0f}; }
inline matrix2 matrix2::identity() { return {1.0f, 0.0f, 0.0f, 1.0f}; }
inline matrix2 matrix2::rotation(radian angle)
{
	float c = angle.cos();
	float s = angle.sin();
	return {
		c, -s,
		s, c
	};
}
inline matrix2 matrix2::rotation(normalized<complex> angle)
{
	return {
		angle->r, -angle->i,
		angle->i, angle->r
	};
}
inline matrix2& matrix2::transpose() { std::swap(m01, m10); return *this; }
inline matrix2 matrix2::get_transposed() const
{
	return {
		m00, m01,
		m10, m11
	};
}

inline float& matrix2::at(uint8_t x, uint8_t y) { return reinterpret_cast<float*>(this)[x + 2 * y]; }
inline const float& matrix2::at(uint8_t x, uint8_t y) const { return reinterpret_cast<const float*>(this)[x + 2 * y]; }

inline vector2 matrix2::operator*(vector2 v) const
{
	return {
		m00 * v.x + m01 * v.y,
		m10 * v.x + m11 * v.y
	};
}
inline matrix2 matrix2::operator*(matrix2 m) const
{
	return {
		m00 * m.m00 + m01 * m.m10,
		m00 * m.m01 + m01 * m.m11,

		m10 * m.m00 + m11 * m.m10,
		m10 * m.m01 + m11 * m.m11
	};
}

inline matrix3::matrix3(float m00_, float m10_, float m20_,
	float m01_, float m11_, float m21_,
	float m02_, float m12_, float m22_) :
	m00(m00_), m10(m10_), m20(m20_),
	m01(m01_), m11(m11_), m21(m21_),
	m02(m02_), m12(m12_), m22(m22_)
{}
inline matrix3::matrix3(matrix2 m) :
	matrix3(
		m.m00, m.m10, 0.0f,
		m.m01, m.m11, 0.0f,
		0.0f,  0.0f,  1.0f
	)
{}
inline matrix3::matrix3(const matrix4& m) :
	matrix3(
		m.m00, m.m10, m.m20,
		m.m01, m.m11, m.m21,
		m.m02, m.m12, m.m22
	)
{}
inline matrix3::matrix3(const transform2d& t)
{
	matrix2 r = matrix2::rotation(t.rotation);
	*this = {
		t.scale.x * r.m00, t.scale.y * r.m10, t.position.x,
		t.scale.x * r.m01, t.scale.y * r.m11, t.position.y,
		0.0f,              0.0f,              1.0f
	};
}
inline matrix3::matrix3(uniform_transform2d t)
{
	matrix2 r = matrix2::rotation(t.rotation);
	*this = {
		t.scale * r.m00, t.scale * r.m10, t.position.x,
		t.scale * r.m01, t.scale * r.m11, t.position.y,
		0.0f,            0.0f,            1.0f
	};
}
inline matrix3::matrix3(unscaled_transform2d t)
{
	matrix2 r = matrix2::rotation(t.rotation);
	*this = {
		r.m00, r.m10, t.position.x,
		r.m01, r.m11, t.position.y,
		0.0f, 0.0f, 1.0f
	};
}

inline matrix3 matrix3::zero() { return {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}; }
inline matrix3 matrix3::identity() { return {1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f}; }
inline matrix3 matrix3::rotation(euler_angles r)
{
	matrix3 RX = matrix3::identity();
	RX.m11 = r.x.cos();
	RX.m12 = -r.x.sin();
	RX.m21 = r.x.sin();
	RX.m22 = r.x.cos();

	matrix3 RY = matrix3::identity();
	RY.m00 = r.y.cos();
	RY.m02 = r.y.sin();
	RY.m20 = -r.y.sin();
	RY.m22 = r.y.cos();

	matrix3 RZ = matrix3::identity();
	RZ.m00 = r.z.cos();
	RZ.m01 = -r.z.sin();
	RZ.m10 = r.z.sin();
	RZ.m11 = r.z.cos();

	//return RZ * RX * RY;
	//return RZ * RY * RX;
	return RY * RX * RZ;
	//return RY * RZ * RX;
	//return RY * RX * RZ;
	//return RX * RZ * RY;
}
inline matrix3 matrix3::rotation(quaternion q)
{
	return {
		detail::get_quaternion_matrix_element<0, 0>(q), detail::get_quaternion_matrix_element<1, 0>(q), detail::get_quaternion_matrix_element<2, 0>(q),
		detail::get_quaternion_matrix_element<0, 1>(q), detail::get_quaternion_matrix_element<1, 1>(q), detail::get_quaternion_matrix_element<2, 1>(q),
		detail::get_quaternion_matrix_element<0, 2>(q), detail::get_quaternion_matrix_element<1, 2>(q), detail::get_quaternion_matrix_element<2, 2>(q),
	};
}
inline matrix3 matrix3::rotation_around_x(radian angle)
{
	float c = angle.cos();
	float s = angle.sin();
	return {
		1.0f, 0.0f, 0.0f,
		0.0f, c,    s,
		0.0f, -s,   c
	};
}
inline matrix3 matrix3::rotation_around_y(radian angle)
{
	float c = angle.cos();
	float s = angle.sin();
	return {
		c,    0.0f, s,
		0.0f, 1.0f, 0.0f,
		-s,   0.0f, c
	};
}
inline matrix3 matrix3::rotation_around_z(radian angle)
{
	float c = angle.cos();
	float s = angle.sin();
	return {
		c,    -s,   0.0f,
		s,    c,    0.0f,
		0.0f, 0.0f, 1.0f
	};
}
inline matrix3 matrix3::rotation_around_x(normalized<complex> angle)
{
	float c = angle->r;
	float s = angle->i;
	return {
		1.0f, 0.0f, 0.0f,
		0.0f, c,    s,
		0.0f, -s,   c
	};
}
inline matrix3 matrix3::rotation_around_y(normalized<complex> angle)
{
	float c = angle->r;
	float s = angle->i;
	return {
		c,    0.0f, s,
		0.0f, 1.0f, 0.0f,
		-s,   0.0f, c
	};
}
inline matrix3 matrix3::rotation_around_z(normalized<complex> angle)
{
	float c = angle->r;
	float s = angle->i;
	return {
		c,    -s,   0.0f,
		s,    c,    0.0f,
		0.0f, 0.0f, 1.0f
	};
}
inline matrix3& matrix3::inverse()
{
	if (is_orthogonal())
	{
		transpose();
		return *this;
	}

	float det = get_determinant();
	float recM00 = m00;
	float recM10 = m10;
	float recM20 = m20;
	float recM01 = m01;
	float recM11 = m11;
	float recM02 = m02;
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
inline matrix3 matrix3::get_inversed() const
{
	matrix3 res = *this;
	res.inverse();
	return res;
}
inline matrix3& matrix3::transpose()
{
	std::swap(m01, m10);
	std::swap(m02, m20);
	std::swap(m12, m21);
	return *this;
}
inline matrix3 matrix3::get_transposed() const
{
	return {
		m00, m10, m20,
		m01, m11, m21,
		m02, m12, m22
	};
}
inline float matrix3::get_determinant() const
{
	return
		m00 * m11 * m22 +
		m01 * m12 * m20 +
		m02 * m10 * m21 -
		m20 * m11 * m02 -
		m21 * m12 * m00 -
		m22 * m10 * m01;
}
inline bool matrix3::is_orthogonal() const
{
	return
		nearly_equal(m00 * m01 + m10 * m11 + m20 * m21, 0.0f) &&
		nearly_equal(m00 * m02 + m10 * m12 + m20 * m22, 0.0f) &&
		nearly_equal(m01 * m02 + m11 * m12 + m21 * m22, 0.0f) &&
		nearly_equal(m00 * m00 + m10 * m10 + m20 * m20, 1.0f) &&
		nearly_equal(m01 * m01 + m11 * m11 + m21 * m21, 1.0f) &&
		nearly_equal(m02 * m02 + m12 * m12 + m22 * m22, 1.0f);
}

inline float& matrix3::at(uint8_t x, uint8_t y) { return reinterpret_cast<float*>(this)[x + 3 * y]; }
inline const float& matrix3::at(uint8_t x, uint8_t y) const { return reinterpret_cast<const float*>(this)[x + 3 * y]; }

inline matrix3 matrix3::operator*(float f) const
{
	return {
		m00 * f, m10 * f, m20 * f,
		m01 * f, m11 * f, m21 * f,
		m02 * f, m12 * f, m22 * f
	};
}
inline vector3 matrix3::operator*(vector3 v) const
{
	return {
		m00 * v.x + m10 * v.y + m20 * v.z,
		m01 * v.x + m11 * v.y + m21 * v.z,
		m02 * v.x + m12 * v.y + m22 * v.z
	};
}
inline matrix3 matrix3::operator*(const matrix3& m) const
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
inline matrix4x3 matrix3::operator*(const matrix4x3& m) const
{
	return {
		m.m00 * m00 + m.m01 * m10 + m.m02 * m20,
		m.m10 * m00 + m.m11 * m10 + m.m12 * m20,
		m.m20 * m00 + m.m21 * m10 + m.m22 * m20,
		m.m30 * m00 + m.m31 * m10 + m.m32 * m20,

		m.m00 * m01 + m.m01 * m11 + m.m02 * m21,
		m.m10 * m01 + m.m11 * m11 + m.m12 * m21,
		m.m20 * m01 + m.m21 * m11 + m.m22 * m21,
		m.m30 * m01 + m.m31 * m11 + m.m32 * m21,

		m.m00 * m02 + m.m01 * m12 + m.m02 * m22,
		m.m10 * m02 + m.m11 * m12 + m.m12 * m22,
		m.m20 * m02 + m.m21 * m12 + m.m22 * m22,
		m.m30 * m02 + m.m31 * m12 + m.m32 * m22
	};
}

inline matrix3x4::matrix3x4(float m00_, float m10_, float m20_,
	float m01_, float m11_, float m21_,
	float m02_, float m12_, float m22_,
	float m03_, float m13_, float m23_) :
	m00(m00_), m10(m10_), m20(m20_),
	m01(m01_), m11(m11_), m21(m21_),
	m02(m02_), m12(m12_), m22(m22_),
	m03(m03_), m13(m13_), m23(m23_)
{}
inline matrix3x4::matrix3x4(matrix2 m) :
	matrix3x4(
		m.m00, m.m10, 0.0f,
		m.m01, m.m11, 0.0f,
		0.0f, 0.0f, 1.0f,
		0.0f, 0.0f, 0.0f
	)
{}
inline matrix3x4::matrix3x4(const matrix3& m) :
	matrix3x4(
		m.m00, m.m10, m.m20,
		m.m01, m.m11, m.m21,
		m.m02, m.m12, m.m22,
		0.0f, 0.0f, 0.0f
	)
{}
inline matrix3x4::matrix3x4(const matrix4& m) :
	matrix3x4(
		m.m00, m.m10, m.m20,
		m.m01, m.m11, m.m21,
		m.m02, m.m12, m.m22,
		m.m03, m.m13, m.m23
	)
{}
inline matrix3x4::matrix3x4(const transform3d& t)
{
	float xx = t.rotation.x * t.rotation.x;
	float yy = t.rotation.y * t.rotation.y;
	float zz = t.rotation.z * t.rotation.z;
	float wx = t.rotation.w * t.rotation.x;
	float wy = t.rotation.w * t.rotation.y;
	float wz = t.rotation.w * t.rotation.z;
	float xy = t.rotation.x * t.rotation.y;
	float xz = t.rotation.x * t.rotation.z;
	float yz = t.rotation.y * t.rotation.z;
	*this = {
		t.scale.x * (1.0f - 2.0f * (yy + zz)), t.scale.y * (2.0f * (xy - wz)),        t.scale.z * (2.0f * (xz + wy)),
		t.scale.x * (2.0f * (xy + wz)),        t.scale.y * (1.0f - 2.0f * (xx + zz)), t.scale.z * (2.0f * (yz - wx)),
		t.scale.x * (2.0f * (xz - wy)),        t.scale.y * (2.0f * (yz + wx)),        t.scale.z * (1.0f - 2.0f * (xx + yy)),
		0.0f,                                  0.0f,                                  0.0f
	};
}
inline matrix3x4::matrix3x4(const uniform_transform3d& t)
{
	float xx = t.rotation.x * t.rotation.x;
	float yy = t.rotation.y * t.rotation.y;
	float zz = t.rotation.z * t.rotation.z;
	float wx = t.rotation.w * t.rotation.x;
	float wy = t.rotation.w * t.rotation.y;
	float wz = t.rotation.w * t.rotation.z;
	float xy = t.rotation.x * t.rotation.y;
	float xz = t.rotation.x * t.rotation.z;
	float yz = t.rotation.y * t.rotation.z;
	*this = {
		t.scale * (1.0f - 2.0f * (yy + zz)), t.scale * (2.0f * (xy - wz)),        t.scale * (2.0f * (xz + wy)),
		t.scale * (2.0f * (xy + wz)),        t.scale * (1.0f - 2.0f * (xx + zz)), t.scale * (2.0f * (yz - wx)),
		t.scale * (2.0f * (xz - wy)),        t.scale * (2.0f * (yz + wx)),        t.scale * (1.0f - 2.0f * (xx + yy)),
		0.0f,                                0.0f,                                0.0f
	};
}
inline matrix3x4::matrix3x4(const unscaled_transform3d& t)
{
	float xx = t.rotation.x * t.rotation.x;
	float yy = t.rotation.y * t.rotation.y;
	float zz = t.rotation.z * t.rotation.z;
	float wx = t.rotation.w * t.rotation.x;
	float wy = t.rotation.w * t.rotation.y;
	float wz = t.rotation.w * t.rotation.z;
	float xy = t.rotation.x * t.rotation.y;
	float xz = t.rotation.x * t.rotation.z;
	float yz = t.rotation.y * t.rotation.z;
	*this = {
		(1.0f - 2.0f * (yy + zz)), (2.0f * (xy - wz)),        (2.0f * (xz + wy)),
		(2.0f * (xy + wz)),        (1.0f - 2.0f * (xx + zz)), (2.0f * (yz - wx)),
		(2.0f * (xz - wy)),        (2.0f * (yz + wx)),        (1.0f - 2.0f * (xx + yy)),
		0.0f,                      0.0f,                      0.0f
	};
}

inline matrix3x4 matrix3x4::zero() { return { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f }; }
inline matrix3x4 matrix3x4::identity() { return { 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f }; }
inline matrix3x4 matrix3x4::rotation(euler_angles r) { return matrix3x4(matrix3::rotation(r)); }
inline matrix3x4 matrix3x4::rotation(quaternion q) { return matrix3x4(matrix3::rotation(q)); }
inline matrix3x4 matrix3x4::rotation_around_x(radian angle) { return matrix3x4(matrix3::rotation_around_x(angle)); }
inline matrix3x4 matrix3x4::rotation_around_y(radian angle) { return matrix3x4(matrix3::rotation_around_y(angle)); }
inline matrix3x4 matrix3x4::rotation_around_z(radian angle) { return matrix3x4(matrix3::rotation_around_z(angle)); }
inline matrix3x4 matrix3x4::rotation_around_x(normalized<complex> angle) { return matrix3x4(matrix3::rotation_around_x(angle)); }
inline matrix3x4 matrix3x4::rotation_around_y(normalized<complex> angle) { return matrix3x4(matrix3::rotation_around_y(angle)); }
inline matrix3x4 matrix3x4::rotation_around_z(normalized<complex> angle) { return matrix3x4(matrix3::rotation_around_z(angle)); }
inline matrix4x3 matrix3x4::get_transposed() const
{
	return {
		m00, m01, m02, m03,
		m10, m11, m12, m13,
		m20, m21, m22, m23
	};
}

inline float& matrix3x4::at(uint8_t x, uint8_t y) { return reinterpret_cast<float*>(this)[x + 3 * y]; }
inline const float& matrix3x4::at(uint8_t x, uint8_t y) const { return reinterpret_cast<const float*>(this)[x + 3 * y]; }

inline matrix3x4 matrix3x4::operator*(float f) const
{
	return {
		m00 * f, m10 * f, m20 * f,
		m01 * f, m11 * f, m21 * f,
		m02 * f, m12 * f, m22 * f,
		m03 * f, m13 * f, m23 * f
	};
}
inline vector4 matrix3x4::operator*(vector3 v) const
{
	return {
		m00 * v.x + m10 * v.y + m20 * v.z,
		m01 * v.x + m11 * v.y + m21 * v.z,
		m02 * v.x + m12 * v.y + m22 * v.z,
		m03 * v.x + m13 * v.y + m23 * v.z
	};
}
inline matrix3x4 matrix3x4::operator*(const matrix3& m) const
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
		   m02 * m.m20 + m12 * m.m21 + m22 * m.m22,

		   m03 * m.m00 + m13 * m.m01 + m23 * m.m02,
		   m03 * m.m10 + m13 * m.m11 + m23 * m.m12,
		   m03 * m.m20 + m13 * m.m21 + m23 * m.m22
	};
}
inline matrix4 matrix3x4::operator*(const matrix4x3& m) const
{
	return {
		m.m00 * m00 + m.m01 * m10 + m.m02 * m20,
		m.m10 * m00 + m.m11 * m10 + m.m12 * m20,
		m.m20 * m00 + m.m21 * m10 + m.m22 * m20,
		m.m30 * m00 + m.m31 * m10 + m.m32 * m20,

		m.m00 * m01 + m.m01 * m11 + m.m02 * m21,
		m.m10 * m01 + m.m11 * m11 + m.m12 * m21,
		m.m20 * m01 + m.m21 * m11 + m.m22 * m21,
		m.m30 * m01 + m.m31 * m11 + m.m32 * m21,

		m.m00 * m02 + m.m01 * m12 + m.m02 * m22,
		m.m10 * m02 + m.m11 * m12 + m.m12 * m22,
		m.m20 * m02 + m.m21 * m12 + m.m22 * m22,
		m.m30 * m02 + m.m31 * m12 + m.m32 * m22,

		m.m00 * m03 + m.m01 * m13 + m.m02 * m23,
		m.m10 * m03 + m.m11 * m13 + m.m12 * m23,
		m.m20 * m03 + m.m21 * m13 + m.m22 * m23,
		m.m30 * m03 + m.m31 * m13 + m.m32 * m23
	};
}

inline matrix4x3::matrix4x3(float m00_, float m10_, float m20_, float m30_,
	float m01_, float m11_, float m21_, float m31_,
	float m02_, float m12_, float m22_, float m32_) :
	m00(m00_), m10(m10_), m20(m20_), m30(m30_),
	m01(m01_), m11(m11_), m21(m21_), m31(m31_),
	m02(m02_), m12(m12_), m22(m22_), m32(m32_)
{}
inline matrix4x3::matrix4x3(matrix2 m) :
	matrix4x3(
		m.m00, m.m10, 0.0f, 0.0f,
		m.m01, m.m11, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f
	)
{}
inline matrix4x3::matrix4x3(const matrix3& m) :
	matrix4x3(
		m.m00, m.m10, m.m20, 0.0f,
		m.m01, m.m11, m.m21, 0.0f,
		m.m02, m.m12, m.m22, 0.0f
	)
{}
inline matrix4x3::matrix4x3(const matrix4& m) :
	matrix4x3(
		m.m00, m.m10, m.m20, m.m30,
		m.m01, m.m11, m.m21, m.m31,
		m.m02, m.m12, m.m22, m.m32
	)
{}
inline matrix4x3::matrix4x3(const transform3d& t)
{
	float xx = t.rotation.x * t.rotation.x;
	float yy = t.rotation.y * t.rotation.y;
	float zz = t.rotation.z * t.rotation.z;
	float wx = t.rotation.w * t.rotation.x;
	float wy = t.rotation.w * t.rotation.y;
	float wz = t.rotation.w * t.rotation.z;
	float xy = t.rotation.x * t.rotation.y;
	float xz = t.rotation.x * t.rotation.z;
	float yz = t.rotation.y * t.rotation.z;
	*this = {
		t.scale.x * (1.0f - 2.0f * (yy + zz)), t.scale.y * (2.0f * (xy - wz)),        t.scale.z * (2.0f * (xz + wy)),        t.position.x,
		t.scale.x * (2.0f * (xy + wz)),        t.scale.y * (1.0f - 2.0f * (xx + zz)), t.scale.z * (2.0f * (yz - wx)),        t.position.y,
		t.scale.x * (2.0f * (xz - wy)),        t.scale.y * (2.0f * (yz + wx)),        t.scale.z * (1.0f - 2.0f * (xx + yy)), t.position.z
	};
}
inline matrix4x3::matrix4x3(const uniform_transform3d& t)
{
	float xx = t.rotation.x * t.rotation.x;
	float yy = t.rotation.y * t.rotation.y;
	float zz = t.rotation.z * t.rotation.z;
	float wx = t.rotation.w * t.rotation.x;
	float wy = t.rotation.w * t.rotation.y;
	float wz = t.rotation.w * t.rotation.z;
	float xy = t.rotation.x * t.rotation.y;
	float xz = t.rotation.x * t.rotation.z;
	float yz = t.rotation.y * t.rotation.z;
	*this = {
		t.scale * (1.0f - 2.0f * (yy + zz)), t.scale * (2.0f * (xy - wz)),        t.scale * (2.0f * (xz + wy)),        t.position.x,
		t.scale * (2.0f * (xy + wz)),        t.scale * (1.0f - 2.0f * (xx + zz)), t.scale * (2.0f * (yz - wx)),        t.position.y,
		t.scale * (2.0f * (xz - wy)),        t.scale * (2.0f * (yz + wx)),        t.scale * (1.0f - 2.0f * (xx + yy)), t.position.z
	};
}
inline matrix4x3::matrix4x3(const unscaled_transform3d& t)
{
	float xx = t.rotation.x * t.rotation.x;
	float yy = t.rotation.y * t.rotation.y;
	float zz = t.rotation.z * t.rotation.z;
	float wx = t.rotation.w * t.rotation.x;
	float wy = t.rotation.w * t.rotation.y;
	float wz = t.rotation.w * t.rotation.z;
	float xy = t.rotation.x * t.rotation.y;
	float xz = t.rotation.x * t.rotation.z;
	float yz = t.rotation.y * t.rotation.z;
	*this = {
		(1.0f - 2.0f * (yy + zz)), (2.0f * (xy - wz)),        (2.0f * (xz + wy)),        t.position.x,
		(2.0f * (xy + wz)),        (1.0f - 2.0f * (xx + zz)), (2.0f * (yz - wx)),        t.position.y,
		(2.0f * (xz - wy)),        (2.0f * (yz + wx)),        (1.0f - 2.0f * (xx + yy)), t.position.z
	};
}

inline matrix4x3 matrix4x3::zero() { return { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f }; }
inline matrix4x3 matrix4x3::identity() { return { 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f }; }
inline matrix4x3 matrix4x3::rotation(euler_angles r) { return matrix4x3(matrix3::rotation(r)); }
inline matrix4x3 matrix4x3::rotation(quaternion q) { return matrix4x3(matrix3::rotation(q)); }
inline matrix4x3 matrix4x3::translation(vector3 t) { return { 1.0f, 0.0f, 0.0f, t.x, 0.0f, 1.0f, 0.0f, t.y, 0.0f, 0.0f, 1.0f, t.z }; }
inline matrix4x3 matrix4x3::rotation_around_x(radian angle) { return matrix4x3(matrix3::rotation_around_x(angle)); }
inline matrix4x3 matrix4x3::rotation_around_y(radian angle) { return matrix4x3(matrix3::rotation_around_y(angle)); }
inline matrix4x3 matrix4x3::rotation_around_z(radian angle) { return matrix4x3(matrix3::rotation_around_z(angle)); }
inline matrix4x3 matrix4x3::rotation_around_x(normalized<complex> angle) { return matrix4x3(matrix3::rotation_around_x(angle)); }
inline matrix4x3 matrix4x3::rotation_around_y(normalized<complex> angle) { return matrix4x3(matrix3::rotation_around_y(angle)); }
inline matrix4x3 matrix4x3::rotation_around_z(normalized<complex> angle) { return matrix4x3(matrix3::rotation_around_z(angle)); }
inline matrix3x4 matrix4x3::get_transposed() const
{
	return {
		m00, m01, m02,
		m10, m11, m12,
		m20, m21, m22,
		m30, m31, m32
	};
}

inline float& matrix4x3::at(uint8_t x, uint8_t y) { return reinterpret_cast<float*>(this)[x + 4 * y]; }
inline const float& matrix4x3::at(uint8_t x, uint8_t y) const { return reinterpret_cast<const float*>(this)[x + 4 * y]; }

inline matrix4x3 matrix4x3::operator*(float f) const
{
	return {
		m00 * f, m10 * f, m20 * f, m20 * f,
		m01 * f, m11 * f, m21 * f, m21 * f,
		m02 * f, m12 * f, m22 * f, m22 * f
	};
}
inline vector3 matrix4x3::operator*(vector4 v) const
{
	return {
		m00 * v.x + m10 * v.y + m20 * v.z + m30 * v.w,
		m01 * v.x + m11 * v.y + m21 * v.z + m31 * v.w,
		m02 * v.x + m12 * v.y + m22 * v.z + m32 * v.w
	};
}
inline matrix3 matrix4x3::operator*(const matrix3x4& m) const
{
	return {
		m.m00 * m00 + m.m01 * m10 + m.m02 * m20 + m.m03 * m30,
		m.m10 * m00 + m.m11 * m10 + m.m12 * m20 + m.m13 * m30,
		m.m20 * m00 + m.m21 * m10 + m.m22 * m20 + m.m23 * m30,

		m.m00 * m01 + m.m01 * m11 + m.m02 * m21 + m.m03 * m31,
		m.m10 * m01 + m.m11 * m11 + m.m12 * m21 + m.m13 * m31,
		m.m20 * m01 + m.m21 * m11 + m.m22 * m21 + m.m23 * m31,

		m.m00 * m02 + m.m01 * m12 + m.m02 * m22 + m.m03 * m32,
		m.m10 * m02 + m.m11 * m12 + m.m12 * m22 + m.m13 * m32,
		m.m20 * m02 + m.m21 * m12 + m.m22 * m22 + m.m23 * m32
	};
}
inline matrix4x3 matrix4x3::operator*(const matrix4& m) const
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
		m.m30 * m02 + m.m31 * m12 + m.m32 * m22 + m.m33 * m32
	};
}

inline matrix4::matrix4(float m00_, float m10_, float m20_, float m30_,
	float m01_, float m11_, float m21_, float m31_,
	float m02_, float m12_, float m22_, float m32_,
	float m03_, float m13_, float m23_, float m33_) :
	m00(m00_), m10(m10_), m20(m20_), m30(m30_),
	m01(m01_), m11(m11_), m21(m21_), m31(m31_),
	m02(m02_), m12(m12_), m22(m22_), m32(m32_),
	m03(m03_), m13(m13_), m23(m23_), m33(m33_)
{}
inline matrix4::matrix4(matrix2 m) :
	matrix4(
		m.m00, m.m10, 0.0f, 0.0f,
		m.m01, m.m11, 0.0f, 0.0f,
		0.0f,  0.0f,  1.0f, 0.0f,
		0.0f,  0.0f,  0.0f, 1.0f
	)
{}
inline matrix4::matrix4(const matrix3& m) :
	matrix4(
		m.m00, m.m10, m.m20, 0.0f,
		m.m01, m.m11, m.m21, 0.0f,
		m.m02, m.m12, m.m22, 0.0f,
		0.0f,  0.0f,  0.0f,  1.0f
	)
{}
inline matrix4::matrix4(const perspective& p) { *this = p.to_matrix4(); }
inline matrix4::matrix4(const transform3d& t)
{
	float xx = t.rotation.x * t.rotation.x;
	float yy = t.rotation.y * t.rotation.y;
	float zz = t.rotation.z * t.rotation.z;
	float wx = t.rotation.w * t.rotation.x;
	float wy = t.rotation.w * t.rotation.y;
	float wz = t.rotation.w * t.rotation.z;
	float xy = t.rotation.x * t.rotation.y;
	float xz = t.rotation.x * t.rotation.z;
	float yz = t.rotation.y * t.rotation.z;
	*this = {
		t.scale.x * (1.0f - 2.0f * (yy + zz)), t.scale.y * (2.0f * (xy - wz)),        t.scale.z * (2.0f * (xz + wy)),        t.position.x,
		t.scale.x * (2.0f * (xy + wz)),        t.scale.y * (1.0f - 2.0f * (xx + zz)), t.scale.z * (2.0f * (yz - wx)),        t.position.y,
		t.scale.x * (2.0f * (xz - wy)),        t.scale.y * (2.0f * (yz + wx)),        t.scale.z * (1.0f - 2.0f * (xx + yy)), t.position.z,
		0.0f,                                  0.0f,                                  0.0f,                                  1.0f
	};
}
inline matrix4::matrix4(const uniform_transform3d& t)
{
	float xx = t.rotation.x * t.rotation.x;
	float yy = t.rotation.y * t.rotation.y;
	float zz = t.rotation.z * t.rotation.z;
	float wx = t.rotation.w * t.rotation.x;
	float wy = t.rotation.w * t.rotation.y;
	float wz = t.rotation.w * t.rotation.z;
	float xy = t.rotation.x * t.rotation.y;
	float xz = t.rotation.x * t.rotation.z;
	float yz = t.rotation.y * t.rotation.z;
	*this = {
		t.scale * (1.0f - 2.0f * (yy + zz)), t.scale * (2.0f * (xy - wz)),        t.scale * (2.0f * (xz + wy)),        t.position.x,
		t.scale * (2.0f * (xy + wz)),        t.scale * (1.0f - 2.0f * (xx + zz)), t.scale * (2.0f * (yz - wx)),        t.position.y,
		t.scale * (2.0f * (xz - wy)),        t.scale * (2.0f * (yz + wx)),        t.scale * (1.0f - 2.0f * (xx + yy)), t.position.z,
		0.0f,                                0.0f,                                0.0f,                                1.0f
	};
}
inline matrix4::matrix4(const unscaled_transform3d& t)
{
	float xx = t.rotation.x * t.rotation.x;
	float yy = t.rotation.y * t.rotation.y;
	float zz = t.rotation.z * t.rotation.z;
	float wx = t.rotation.w * t.rotation.x;
	float wy = t.rotation.w * t.rotation.y;
	float wz = t.rotation.w * t.rotation.z;
	float xy = t.rotation.x * t.rotation.y;
	float xz = t.rotation.x * t.rotation.z;
	float yz = t.rotation.y * t.rotation.z;
	*this = {
		(1.0f - 2.0f * (yy + zz)), (2.0f * (xy - wz)),        (2.0f * (xz + wy)),        t.position.x,
		(2.0f * (xy + wz)),        (1.0f - 2.0f * (xx + zz)), (2.0f * (yz - wx)),        t.position.y,
		(2.0f * (xz - wy)),        (2.0f * (yz + wx)),        (1.0f - 2.0f * (xx + yy)), t.position.z,
		0.0f,                      0.0f,                      0.0f,                      1.0f
	};
}
inline matrix4::matrix4(const std::array<float, 16>& a) : matrix4(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15]) {}
inline matrix4::matrix4(const std::array<double, 16>& a) : matrix4(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15]) {}

inline matrix4 matrix4::zero() { return {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}; }
inline matrix4 matrix4::identity() { return {1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f}; }
inline matrix4 matrix4::rotation(euler_angles r) { return matrix4(matrix3::rotation(r)); }
inline matrix4 matrix4::rotation(quaternion q) { return matrix4(matrix3::rotation(q)); }
inline matrix4 matrix4::translation(vector3 t) { return {1.0f, 0.0f, 0.0f, t.x, 0.0f, 1.0f, 0.0f, t.y, 0.0f, 0.0f, 1.0f, t.z, 0.0f, 0.0f, 0.0f, 1.0f}; }
inline matrix4 matrix4::translation(float x, float y, float z) { return translation({x, y, z}); }
inline matrix4 matrix4::scale(vector3 s) { return {s.x, 0.0f, 0.0f, 0.0f, 0.0f, s.y, 0.0f, 0.0f, 0.0f, 0.0f, s.z, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f}; }
inline matrix4 matrix4::look_at(vector3 from, vector3 to, vector3 up)
{
	vector3 forward = (from - to).get_normalized();
    vector3 right = up.get_normalized().cross(forward);
    vector3 true_up = forward.cross(right);

	return {
		right.x, true_up.x, forward.x, from.x,
		right.y, true_up.y, forward.y, from.y,
		right.z, true_up.z, forward.z, from.z,
		0.0f,    0.0f,      0.0f,      1.0f
	};
}
inline matrix4 matrix4::orthographic(float left, float right, float top, float bottom, float near, float far)
{
	float a = 1.0f / (right - left);
	float b = 1.0f / (bottom - top);
	float c = 1.0f / (near - far);

	return {
		2.0f * a, 0.0f,     0.0f, -(right + left) * a,
		0.0f,     2.0f * b, 0.0f, -(bottom + top) * b,
		0.0f,     0.0f,     c,    near * c,
		0.0f,     0.0f,     0.0f, 1.0f,
	};
}
inline matrix4 matrix4::rotation_around_x(radian angle) { return matrix4(matrix3::rotation_around_x(angle)); }
inline matrix4 matrix4::rotation_around_y(radian angle) { return matrix4(matrix3::rotation_around_y(angle)); }
inline matrix4 matrix4::rotation_around_z(radian angle) { return matrix4(matrix3::rotation_around_z(angle)); }
inline matrix4 matrix4::rotation_around_x(normalized<complex> angle) { return matrix4(matrix3::rotation_around_x(angle)); }
inline matrix4 matrix4::rotation_around_y(normalized<complex> angle) { return matrix4(matrix3::rotation_around_y(angle)); }
inline matrix4 matrix4::rotation_around_z(normalized<complex> angle) { return matrix4(matrix3::rotation_around_z(angle)); }
inline bool matrix4::is_orthogonal() const
{
	return
		nearly_equal(m00 * m01 + m10 * m11 + m20 * m21 + m30 * m31, 0.0f) &&
		nearly_equal(m00 * m02 + m10 * m12 + m20 * m22 + m30 * m32, 0.0f) &&
		nearly_equal(m00 * m03 + m10 * m13 + m20 * m23 + m30 * m33, 0.0f) &&
		nearly_equal(m01 * m02 + m11 * m12 + m21 * m22 + m31 * m32, 0.0f) &&
		nearly_equal(m01 * m03 + m11 * m13 + m21 * m23 + m31 * m33, 0.0f) &&
		nearly_equal(m02 * m03 + m12 * m13 + m22 * m23 + m32 * m33, 0.0f) &&
		nearly_equal(m00 * m00 + m10 * m10 + m20 * m20 + m30 * m30, 1.0f) &&
		nearly_equal(m01 * m01 + m11 * m11 + m21 * m21 + m31 * m31, 1.0f) &&
		nearly_equal(m02 * m02 + m12 * m12 + m22 * m22 + m32 * m32, 1.0f) &&
		nearly_equal(m03 * m03 + m13 * m13 + m23 * m23 + m33 * m33, 1.0f);
}
inline bool matrix4::is_homogenous() const
{
	return
		nearly_equal(m03, 0.0f) &&
		nearly_equal(m13, 0.0f) &&
		nearly_equal(m23, 0.0f) &&
		nearly_equal(m30, 0.0f) &&
		nearly_equal(m31, 0.0f) &&
		nearly_equal(m32, 0.0f) &&
		nearly_equal(m33, 1.0f);
}
inline matrix4& matrix4::inverse()
{
	if (is_orthogonal())
	{
		transpose();
		return *this;
	}

	float det = get_determinant();
	float invDet = (1.0f / det);
	float recM00 = m00;
	float recM10 = m10;
	float recM20 = m20;
	float recM01 = m01;
	float recM11 = m11;
	float recM02 = m02;

	if (is_homogenous())
	{
		m00 = invDet * (recM11 * m22    - m21    * m12);
		m10 = invDet * (m12    * recM20	- recM10 * m22);
		m20 = invDet * (recM10 * m21    - recM20 * recM11);
		m01 = invDet * (recM02 * m21    - recM01 * m22);
		m11 = invDet * (recM00 * m22    - recM20 * recM02);
		m21 = invDet * (recM01 * recM20 - recM00 * m21);
		m02 = invDet * (recM01 * m12    - recM11 * recM02);
		m12 = invDet * (recM02 * recM10 - recM00 * m12);
		m22 = invDet * (recM00 * recM11 - recM10 * recM01);
		return *this;
	}

	float recM30 = m30;
	float recM21 = m21;
	float recM31 = m31;
	float recM12 = m12;
	float recM22 = m22;
	float recM03 = m03;
	float recM13 = m13;

	m00 = invDet * (
		recM11 * recM22 * m33 -
		recM11 * m23 * m32 -
		recM21 * recM12 * m33 +
		recM21 * recM13 * m32 +
		recM31 * recM12 * m23 -
		recM31 * recM13 * recM22);

	m10 = invDet * (
		-recM10 * recM22 * m33 +
		recM10 * m23 * m32 +
		recM20 * recM12 * m33 -
		recM20 * recM13 * m32 -
		recM30 * recM12 * m23 +
		recM30 * recM13 * recM22
	);

	m20 = invDet * (
		recM10 * recM21 * m33 -
		recM10 * m23 * recM31 -
		recM20 * recM11 * m33 +
		recM20 * recM13 * recM31 +
		recM30 * recM11 * m23 -
		recM30 * recM13 * recM21
	);

	m30 = invDet * (
		-recM10 * recM21 * m32 +
		recM10 * recM22 * recM31 +
		recM20 * recM11 * m32 -
		recM20 * recM12 * recM31 -
		recM30 * recM11 * recM22 +
		recM30 * recM12 * recM21
	);

	m01 = invDet * (
		-recM01 * recM22 * m33 +
		recM01 * m23 * m32 +
		recM21 * recM02 * m33 -
		recM21 * recM03 * m32 -
		recM31 * recM02 * m23 +
		recM31 * recM03 * recM22
	);

	m11 = invDet * (
		recM00 * recM22 * m33 -
		recM00 * m23 * m32 -
		recM20 * recM02 * m33 +
		recM20 * recM03 * m32 +
		recM30 * recM02 * m23 -
		recM30 * recM03 * recM22
	);

	m21 = invDet * (
		-recM00 * recM21 * m33 +
		recM00 * m23 * recM31 +
		recM20 * recM01 * m33 -
		recM20 * recM03 * recM31 -
		recM30 * recM01 * m23 +
		recM30 * recM03 * recM21
	);

	m31 = invDet * (
		recM00 * recM21 * m32 -
		recM00 * recM22 * recM31 -
		recM20 * recM01 * m32 +
		recM20 * recM02 * recM31 +
		recM30 * recM01 * recM22 -
		recM30 * recM02 * recM21
	);

	m02 = invDet * (
		recM01 * recM12 * m33 -
		recM01 * recM13 * m32 -
		recM11 * recM02 * m33 +
		recM11 * recM03 * m32 +
		recM31 * recM02 * recM13 -
		recM31 * recM03 * recM12
	);

	m12 = invDet * (
		-recM00 * recM12 * m33 +
		recM00 * recM13 * m32 +
		recM10 * recM02 * m33 -
		recM10 * recM03 * m32 -
		recM30 * recM02 * recM13 +
		recM30 * recM03 * recM12
	);

	m22 = invDet * (
		recM00 * recM11 * m33 -
		recM00 * recM13 * recM31 -
		recM10 * recM01 * m33 +
		recM10 * recM03 * recM31 +
		recM30 * recM01 * recM13 -
		recM30 * recM03 * recM11
	);

	m32 = invDet * (
		-recM00 * recM11 * m32 +
		recM00 * recM12 * recM31 +
		recM10 * recM01 * m32 -
		recM10 * recM02 * recM31 -
		recM30 * recM01 * recM12 +
		recM30 * recM02 * recM11
	);

	m03 = invDet * (
		-recM01 * recM12 * m23 +
		recM01 * recM13 * recM22 +
		recM11 * recM02 * m23 -
		recM11 * recM03 * recM22 -
		recM21 * recM02 * recM13 +
		recM21 * recM03 * recM12
	);

	m13 = invDet * (
		recM00 * recM12 * m23 -
		recM00 * recM13 * recM22 -
		recM10 * recM02 * m23 +
		recM10 * recM03 * recM22 +
		recM20 * recM02 * recM13 -
		recM20 * recM03 * recM12
	);

	m23 = invDet * (
		-recM00 * recM11 * m23 +
		recM00 * recM13 * recM21 +
		recM10 * recM01 * m23 -
		recM10 * recM03 * recM21 -
		recM20 * recM01 * recM13 +
		recM20 * recM03 * recM11
	);

	m33 = invDet * (
		recM00 * recM11 * recM22 -
		recM00 * recM12 * recM21 -
		recM10 * recM01 * recM22 +
		recM10 * recM02 * recM21 +
		recM20 * recM01 * recM12 -
		recM20 * recM02 * recM11
	);

	return *this;
}
inline matrix4 matrix4::get_inversed() const
{
	matrix4 res = *this;
	res.inverse();
	return res;
}
inline matrix4& matrix4::transpose()
{
	std::swap(m01, m10);
	std::swap(m02, m20);
	std::swap(m03, m30);
	std::swap(m12, m21);
	std::swap(m13, m31);
	std::swap(m23, m32);
	return *this;
}
inline matrix4 matrix4::get_transposed() const
{
	matrix4 res = *this;
	res.transpose();
	return res;
}
inline float matrix4::get_determinant() const
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

	float det1 = m11 * (m22 * m33 - m32 * m23) - m21 * (m12 * m33 - m32 * m13) + m31 * (m12 * m23 - m22 * m13);
	float det2 = m01 * (m22 * m33 - m32 * m23) - m21 * (m02 * m33 - m32 * m03) + m31 * (m02 * m23 - m22 * m03);
	float det3 = m01 * (m12 * m33 - m32 * m13) - m11 * (m02 * m33 - m32 * m03) + m31 * (m02 * m13 - m12 * m03);
	float det4 = m01 * (m12 * m23 - m22 * m13) - m11 * (m02 * m23 - m22 * m03) + m21 * (m02 * m13 - m12 * m03);
	return m00 * det1 - m10 * det2 + m20 * det3 - m30 * det4;
}

inline float& matrix4::at(uint8_t x, uint8_t y) { return reinterpret_cast<float*>(this)[x + 4 * y]; }
inline const float& matrix4::at(uint8_t x, uint8_t y) const { return reinterpret_cast<const float*>(this)[x + 4 * y]; }

inline matrix4 matrix4::operator*(float f) const
{
	return {
		m00 * f, m10 * f, m20 * f, m30 * f,
		m00 * f, m10 * f, m20 * f, m30 * f,
		m00 * f, m10 * f, m20 * f, m30 * f,
		m00 * f, m10 * f, m20 * f, m30 * f
	};
}
inline matrix4 matrix4::operator/(float f) const
{
	return {
		m00 / f, m10 / f, m20 / f, m30 / f,
		m00 / f, m10 / f, m20 / f, m30 / f,
		m00 / f, m10 / f, m20 / f, m30 / f,
		m00 / f, m10 / f, m20 / f, m30 / f
	};
}
inline vector4 matrix4::operator*(vector4 v) const
{
	return {
		m00 * v.x + m10 * v.y + m20 * v.z + m30 * v.w,
		m01 * v.x + m11 * v.y + m21 * v.z + m31 * v.w,
		m02 * v.x + m12 * v.y + m22 * v.z + m32 * v.w,
		m03 * v.x + m13 * v.y + m23 * v.z + m33 * v.w
	};
}
inline matrix3x4 matrix4::operator*(const matrix3x4& m) const
{
	return {
		m.m00 * m00 + m.m01 * m10 + m.m02 * m20 + m.m03 * m30,
		m.m10 * m00 + m.m11 * m10 + m.m12 * m20 + m.m13 * m30,
		m.m20 * m00 + m.m21 * m10 + m.m22 * m20 + m.m23 * m30,

		m.m00 * m01 + m.m01 * m11 + m.m02 * m21 + m.m03 * m31,
		m.m10 * m01 + m.m11 * m11 + m.m12 * m21 + m.m13 * m31,
		m.m20 * m01 + m.m21 * m11 + m.m22 * m21 + m.m23 * m31,

		m.m00 * m02 + m.m01 * m12 + m.m02 * m22 + m.m03 * m32,
		m.m10 * m02 + m.m11 * m12 + m.m12 * m22 + m.m13 * m32,
		m.m20 * m02 + m.m21 * m12 + m.m22 * m22 + m.m23 * m32,

		m.m00 * m03 + m.m01 * m13 + m.m02 * m23 + m.m03 * m33,
		m.m10 * m03 + m.m11 * m13 + m.m12 * m23 + m.m13 * m33,
		m.m20 * m03 + m.m21 * m13 + m.m22 * m23 + m.m23 * m33
	};
}
inline matrix4 matrix4::operator*(const matrix4& m) const
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

inline matrix4::operator std::array<float, 16>() const
{
	std::array<float, 16> res {
		m00, m10, m20, m30,
		m01, m11, m21, m31,
		m02, m12, m22, m32,
		m03, m13, m23, m33
	};

	return res;
}
inline matrix4::operator std::array<double, 16>() const
{
	std::array<double, 16> res {
		m00, m10, m20, m30,
		m01, m11, m21, m31,
		m02, m12, m22, m32,
		m03, m13, m23, m33
	};

	return res;
}

inline perspective::perspective(radian angle, float ratio, float near, float far) :
	m_angle(angle),
    m_ratio(ratio),
    m_inv_ratio(1.0f / m_ratio),
    m_near(near),
    m_far(far),
    m_invTanHalfFovy(1.0f / (m_angle / 2.0f).tan())
{}
inline void perspective::set_angle(radian angle)
	{
	m_angle = angle;
	m_invTanHalfFovy = 1.0f / (angle / 2.0f).tan();
}
inline radian perspective::get_angle() const { return m_angle; }
inline void perspective::set_ratio(float ratio)
{
	m_ratio = ratio;
	m_inv_ratio = 1.0f / m_ratio;
}
inline float perspective::get_ratio() const { return m_ratio; }
inline void perspective::set_near(float near_plane) { m_near = near_plane; }
inline float perspective::get_near() const { return m_near; }
inline void perspective::set_far(float far_plane) { m_far = far_plane; }
inline float perspective::get_far() const { return m_far; }
inline matrix4 perspective::to_matrix4() const
{
	return {
		m_inv_ratio * m_invTanHalfFovy, 0.0f,              0.0f,                     0.0f,
		0.0f,                           -m_invTanHalfFovy, 0.0f,                     0.0f,
		0.0f,                           0.0f,              m_far / (m_near - m_far), -(m_far * m_near) / (m_far - m_near),
		0.0f,                           0.0f,              -1.0f,                    0.0f
	};
}

inline matrix4 operator*(const matrix4& m, const perspective& p)
{
	float a = p.m_inv_ratio * p.m_invTanHalfFovy;
	float b = -p.m_invTanHalfFovy;
	float c = p.m_far / (p.m_near - p.m_far);
	float d = -(p.m_far * p.m_near) / (p.m_far - p.m_near);
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
inline matrix4 operator*(const perspective& p, const matrix4& m)
{
	float a = p.m_inv_ratio * p.m_invTanHalfFovy;
	float b = -p.m_invTanHalfFovy;
	float c = p.m_far / (p.m_near - p.m_far);
	float d = -(p.m_far * p.m_near) / (p.m_far - p.m_near);
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

inline matrix4 operator*(const unscaled_transform3d& t, const perspective& p)
{
	float a = p.m_inv_ratio * p.m_invTanHalfFovy;
	float b = -p.m_invTanHalfFovy;
	float c = p.m_far / (p.m_near - p.m_far);
	float d = -(p.m_far * p.m_near) / (p.m_far - p.m_near);
	float xx = t.rotation.x * t.rotation.x;
	float yy = t.rotation.y * t.rotation.y;
	float zz = t.rotation.z * t.rotation.z;
	float wx = t.rotation.w * t.rotation.x;
	float wy = t.rotation.w * t.rotation.y;
	float wz = t.rotation.w * t.rotation.z;
	float xy = t.rotation.x * t.rotation.y;
	float xz = t.rotation.x * t.rotation.z;
	float yz = t.rotation.y * t.rotation.z;
	return {
		a * (1.0f - 2.0f * (yy + zz)),
		b * (2.0f * (xy - wz)),
		c * (2.0f * (xz + wy)) - t.position.x,
		d * (2.0f * (xz + wy)),

		a * (2.0f * (xy + wz)),
		b * (1.0f - 2.0f * (xx + zz)),
		c * (2.0f * (yz - wx)) - t.position.y,
		d * (2.0f * (yz - wx)),

		a * (2.0f * (xz - wy)),
		b * (2.0f * (yz + wx)),
		c * (1.0f - 2.0f * (xx + yy)) - t.position.z,
		d * (1.0f - 2.0f * (xx + yy)),

		0.0f,
		0.0f,
		-1.0f,
		0.0f
	};
}
inline matrix4 operator*(const perspective& p, const unscaled_transform3d& t)
{
	float a = p.m_inv_ratio * p.m_invTanHalfFovy;
	float b = -p.m_invTanHalfFovy;
	float c = p.m_far / (p.m_near - p.m_far);
	float d = -(p.m_far * p.m_near) / (p.m_far - p.m_near);
	float xx = t.rotation.x * t.rotation.x;
	float yy = t.rotation.y * t.rotation.y;
	float zz = t.rotation.z * t.rotation.z;
	float wx = t.rotation.w * t.rotation.x;
	float wy = t.rotation.w * t.rotation.y;
	float wz = t.rotation.w * t.rotation.z;
	float xy = t.rotation.x * t.rotation.y;
	float xz = t.rotation.x * t.rotation.z;
	float yz = t.rotation.y * t.rotation.z;
	return {
		(1.0f - 2.0f * (yy + zz)) * a,
		(2.0f * (xy - wz)) * a,
		(2.0f * (xz + wy)) * a,
		t.position.x * a,

		(2.0f * (xy + wz)) * b,
		(1.0f - 2.0f * (xx + zz)) * b,
		(2.0f * (yz - wx)) * b,
		t.position.y * b,

		(2.0f * (xz - wy)) * c,
		(2.0f * (yz + wx)) * c,
		(1.0f - 2.0f * (xx + yy)) * c,
		t.position.z * c + d,

		-(2.0f * (xz - wy)),
		-(2.0f * (yz + wx)),
		-(1.0f - 2.0f * (xx + yy)),
		-t.position.z
	};
}

inline euler_angles::euler_angles(radian x_, radian y_, radian z_) : x(x_), y(y_), z(z_) {}

inline quaternion::quaternion(float x_, float y_, float z_, float w_) : x(x_), y(y_), z(z_), w(w_) {}
inline quaternion::quaternion(const normalized<vector3>& axis, radian angle)
{
	/*vector3 tmpAxis = vector3_get_normalized(axis);
	tmpAxis = vector3_times_float(&tmpAxis, sinf(angle / 2.0f));*/

	vector3 tmpAxis = axis * (angle / 2.0f).sin();
	*this = {tmpAxis.x, tmpAxis.y, tmpAxis.z, (angle / 2.0f).cos()};
}

inline quaternion quaternion::zero() { return {0.0f, 0.0f, 0.0f, 0.0f}; }
inline quaternion quaternion::identity() { return {0.0f, 0.0f, 0.0f, 1.0f}; }
inline quaternion& quaternion::inverse()
{
	float n = norm();
	if (nearly_equal(n, 1.0f))
		return conjugate();

	x /= -n;
	y /= -n;
	z /= -n;
	w /= n;
	return *this;
}
inline quaternion quaternion::get_inversed() const
{
	quaternion res = *this;
	res.inverse();
	return res;
}
inline quaternion& quaternion::conjugate()
{
	x = -x;
	y = -y;
	z = -z;
	return *this;
}
inline quaternion quaternion::get_conjugated() const
{
	quaternion res = *this;
	res.conjugate();
	return res;
}
inline quaternion& quaternion::normalize()
{
	float n = norm();
	x /= n;
	y /= n;
	z /= n;
	w /= n;
	return *this;
}
inline quaternion quaternion::get_normalized() const
{
	quaternion res = *this;
	res.normalize();
	return res;
}
inline quaternion& quaternion::negate()
{
	x = -x;
	y = -y;
	z = -z;
	w = -w;
	return *this;
}
inline quaternion quaternion::get_negated() const
{
	quaternion res = *this;
	res.negate();
	return res;
}
inline quaternion& quaternion::clamp(quaternion min, quaternion max)
{
	x = cdm::clamp(x, min.x, max.x);
	y = cdm::clamp(y, min.y, max.y);
	z = cdm::clamp(z, min.z, max.z);
	w = cdm::clamp(w, min.w, max.w);
	return *this;
}
inline quaternion quaternion::get_clamped(quaternion min, quaternion max) const
{
	quaternion res = *this;
	res.clamp(min, max);
	return res;
}
inline float quaternion::norm() const { return sqrtf(norm_squared()); }
inline float quaternion::norm_squared() const { return x * x + y * y + z * z + w * w; }
inline quaternion quaternion::lerp(quaternion begin, quaternion end, float percent) { return (end - begin) * percent + begin; }
inline quaternion quaternion::nlerp(quaternion begin, quaternion end, float percent) { return lerp(begin, end, percent).get_normalized(); }
inline quaternion quaternion::slerp(quaternion begin, quaternion end, float percent)
{
	quaternion _end;
	float dot = begin.dot(end);

	if (dot > 0.9995f)
		return quaternion::lerp(begin, end, percent);
	if (dot < 0.0f)
	{
		_end = -end;
		dot = -dot;
	}
	else
		_end = end;

	percent = cdm::clamp(percent, 0.0f, 1.0f);
	dot = cdm::clamp(dot, -1.0f, 1.0f);
	float theta = acosf(dot) * percent;

	quaternion res = begin * dot;
	res = _end - res;

	res = begin * cosf(theta) + res * sinf(theta);

	return res;
}
inline float quaternion::dot(quaternion q) const { return x * q.x + y * q.y + z * q.z + w * q.w; }
inline quaternion quaternion::cross(quaternion q) const
{
	return {
		y * q.z - z * q.y,
		z * q.x - x * q.z,
		x * q.y - y * q.x,
		1.0f
	};
}

inline quaternion quaternion::operator+(quaternion q) const { return {x + q.x, y + q.y, z + q.z, w + q.w}; }
inline quaternion quaternion::operator-(quaternion q) const { return {x - q.x, y - q.y, z - q.z, w - q.w}; }
inline quaternion quaternion::operator*(quaternion q) const
{
	quaternion res;
	res.x =  x * q.w + y * q.z - z * q.y + w * q.x;
	res.y = -x * q.z + y * q.w + z * q.x + w * q.y;
	res.z =  x * q.y - y * q.x + z * q.w + w * q.z;
	res.w = -x * q.x - y * q.y - z * q.z + w * q.w;
	return res;
}
inline quaternion quaternion::operator*(float f) const { return {x * f, y * f, z * f, w * f}; }
inline quaternion quaternion::operator/(float f) const { return {x / f, y / f, z / f, w / f}; }
inline vector3 operator*(const normalized<quaternion>& q, vector3 v)
{
	vector3 qvec = {q->x, q->y, q->z};
	vector3 uv = qvec.cross(v);
	vector3 uuv = qvec.cross(uv);
	uv = uv * (2.0f * q->w);
	uuv = uuv * 2.0f;
	return v + uv + uuv;
}

inline quaternion& quaternion::operator+=(quaternion q) { *this = *this + q; return *this; }
inline quaternion& quaternion::operator-=(quaternion q) { *this = *this - q; return *this; }
inline quaternion& quaternion::operator*=(quaternion q) { *this = *this * q; return *this; }
inline quaternion& quaternion::operator*=(float f) { *this = *this * f; return *this; }
inline quaternion& quaternion::operator/=(float f) { *this = *this / f; return *this; }

inline quaternion quaternion::operator-() const { return get_negated(); }

inline bool quaternion::operator==(quaternion q) const { return x == q.x && y == q.y && z == q.z && w == q.w; }
inline bool quaternion::operator!=(quaternion q) const { return !operator==(q); }

inline cartesian_direction2d::cartesian_direction2d(float x_, float y_) : cartesian_direction2d(normalized<vector2>({ x_, y_ })) {}
inline cartesian_direction2d::cartesian_direction2d(const normalized<vector2>& direction) : x(direction->x), y(direction->y) {}
inline cartesian_direction2d::cartesian_direction2d(polar_direction2d direction) : cartesian_direction2d(direction.angle) {}
inline cartesian_direction2d::cartesian_direction2d(radian angle) : cartesian_direction2d(angle.cos(), angle.sin()) {}

inline cartesian_direction2d& cartesian_direction2d::rotate(radian angle)
{
	const float c = angle.cos();
	const float s = angle.sin();
    x = x * c - y * s;
    y = x * s + y * c;

	return *this;
}
inline cartesian_direction2d cartesian_direction2d::get_rotated(radian angle)
{
	cartesian_direction2d res = *this;
	res.rotate(angle);
	return res;
}

inline polar_direction2d::polar_direction2d(float x, float y) : polar_direction2d(radian(atan2(y, x))) {}
inline polar_direction2d::polar_direction2d(const normalized<vector2>& direction) : polar_direction2d(direction->x, direction->y) {}
inline polar_direction2d::polar_direction2d(cartesian_direction2d direction) : polar_direction2d(direction.y, direction.x) {}
inline polar_direction2d::polar_direction2d(radian angle_) : angle(angle_) {}

inline polar_direction2d& polar_direction2d::rotate(radian angle_)
{
	angle += angle_;
	return *this;
}
inline polar_direction2d polar_direction2d::get_rotated(radian angle)
{
	polar_direction2d res = *this;
	res.rotate(angle);
	return res;
}
inline polar_direction2d& polar_direction2d::wrap()
{
	int s = sign(float(angle));
	while (angle > radian(pi) || angle <= -radian(pi))
		angle -= radian(pi) * float(s);
	return *this;
}
inline polar_direction2d polar_direction2d::get_wrapped()
{
	polar_direction2d res = *this;
	res.wrap();
	return res;
}

inline bool line::is_parallel(line l1, line l2)
{
	return nearly_equal(l1.coefficient, l2.coefficient) || nearly_equal(l1.coefficient, -l2.coefficient);
}

inline bool collides(line l1, line l2)
{
	return !line::is_parallel(l1, l2);
}
inline bool collides(line l1, line l2, vector2& intersection)
{
	bool collision = collides(l1, l2);

	if (collision)
	{
		float a = l1.coefficient - l2.coefficient;
		float b = l2.offset - l1.offset;
		intersection.x = b / a;
		intersection.y = l1.coefficient * intersection.x + l1.offset;

		//assert(intersection.y == (l2.coefficient * intersection.x + l2.offset));
	}

	return collision;
}
inline bool collides(line l, vector2 v)
{
	return nearly_equal(l.coefficient * v.x + l.offset, v.y);
}
inline bool collides(vector2 v, const line& l)
{
	return collides(l, v);
}

inline float distance_between(const plane& p, vector3 v)
{
	return -v.x * p.normal->x - v.y * p.normal->y - v.z * p.normal->z;
}
inline float distance_between(vector3 v, const plane& p)
{
	return distance_between(p, v);
}

inline bool aa_rect::contains(vector2 v) const
{
	return (v.x >= origin.x) && (v.x <= origin.x + dimention.x) && (v.y >= origin.y) && (v.y <= origin.y + dimention.y);
}

inline bool collides(const aa_rect& r1, const aa_rect& r2)
{
	return collides(r1, r2.origin) ||
		collides(r1, r2.origin + vector2(r2.dimention.x, 0)) ||
		collides(r1, r2.origin + vector2(0, r2.dimention.y)) ||
		collides(r1, r2.origin + r2.dimention);
}
inline bool collides(const aa_rect&, vector2)
{
	return false;
}
inline bool collides(vector2 v, const aa_rect& r)
{
	return collides(r, v);
}

//inline bool collides(const ray3d& r, const plane& p)
//{
//	float DdotN = r.direction->dot(p.normal);
//	if (std::abs(DdotN) > std::numeric_limits<float>::epsilon())
//	{
//		float t = -(r.origin.dot(p.normal) + p.distance) / DdotN;
//		return t >= 0;
//	}
//	return false;
//}
//inline bool collides(const ray3d& r, const plane& p, vector3& intersection)
//{
//	float DdotN = r.direction->dot(p.normal);
//	if (std::abs(DdotN) > std::numeric_limits<float>::epsilon())
//	{
//		float t = -(r.origin.dot(p.normal) + p.distance) / DdotN;
//		intersection = r.origin + t * r.direction;
//		return t >= 0;
//	}
//	return false;
//}
//inline bool collides(const plane& p, const ray3d& r) { return collides(r, p); }
//inline bool collides(const plane& p, const ray3d& r, vector3& intersection) { return collides(r, p, intersection); }

inline float plane::evaluate(vector3 point) const
{
	return normal->x * (point.x - origin.x) +
	       normal->y * (point.y - origin.y) +
	       normal->z * (point.z - origin.z);
}
inline vector3 plane::project3d(vector3 point) const
{
	float distance = evaluate(point);
	return point - normal * distance;
}
inline vector2 plane::project2d(vector3 point, normalized<vector3> plane_tangent) const
{
	normalized<vector3> bitangent = normal->cross(plane_tangent);
	normalized<vector3> tangent = bitangent->cross(normal);

	matrix3 TBN(tangent->x,   tangent->y,   tangent->z,
	            bitangent->x, bitangent->y, bitangent->z,
	            normal->x,    normal->y,    normal->z
	);

	vector3 plane_space_point = point - origin;

	return (TBN * plane_space_point).xy();
}
inline vector3 plane::unproject(vector2 point, normalized<vector3> plane_tangent) const
{
	normalized<vector3> bitangent = normal->cross(plane_tangent);
	normalized<vector3> tangent = bitangent->cross(normal);

	matrix3 invTBN = matrix3(tangent->x,   tangent->y,   tangent->z,
	                         bitangent->x, bitangent->y, bitangent->z,
	                         normal->x,    normal->y,    normal->z
	).get_inversed();

	vector3 plane_space_point = vector3(point, 0.0f);

	return (invTBN * plane_space_point) + origin;
}

inline bool collides(const plane& plane, const ray3d& r, vector3& collision_point)
{
	float denom = plane.normal->dot(r.direction);
	if (abs(denom) > 0.0001f)
	{
		float t = -(plane.origin - r.origin).dot(plane.normal) / denom;
		if (t >= 0.0f)
		{
			collision_point = r.origin + r.direction * t;
			return true;
		}
	}
	return false;
}
inline bool collides_bidirectional(const plane& plane, const ray3d& r, vector3& collision_point)
{
	float denom = plane.normal->dot(r.direction);
	if (abs(denom) > 0.0001f)
	{
		float t = -(plane.origin - r.origin).dot(plane.normal) / denom;
		//if (t >= 0.0f)
		//{
			collision_point = r.origin + r.direction * t;
			return true;
		//}
		//else
		//{
		//	collision_point = r.origin - r.direction * t;
		//	return true;
		//}
	}
	return false;
}
//inline bool collides(const plane& p, const ray3d& r, vector3& collision_point)
//{
//	p.
//}
//inline bool collides(const ray3d& r, const plane& p, vector3& collision_point)
//{
//	return collides(p, r, collision_point);
//}

inline bool aabb::contains(vector3 p) const
{
	return p.x >= min.x && p.x <= max.x &&
		   p.y >= min.y && p.y <= max.y &&
	       p.z >= min.z && p.z <= max.z;
}
inline vector3 aabb::get_center() const { return (max + min) / 2.0f; }

inline aabb aabb::operator+(aabb rhs) const
{
	aabb res;
	res.min.x = std::min(min.x, rhs.min.x);
	res.min.y = std::min(min.y, rhs.min.y);
	res.min.z = std::min(min.z, rhs.min.z);
	res.max.x = std::max(max.x, rhs.max.x);
	res.max.y = std::max(max.y, rhs.max.y);
	res.max.z = std::max(max.z, rhs.max.z);

	return res;
}

inline bool collides(aabb b, ray3d r)
{
	vector3 inv
	{
		1.0f / r.direction->x,
		1.0f / r.direction->y,
		1.0f / r.direction->z
	};

	float t1 = (b.min.x - r.origin.x) * inv.x;
	float t2 = (b.max.x - r.origin.x) * inv.x;

	float tmin = std::min(t1, t2);
	float tmax = std::max(t1, t2);

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
inline bool collides(const aabb& b1, const aabb& b2)
{
	if (b1.contains(b2.min))
		return true;
	if (b1.contains(b2.max))
		return true;
	if (b2.contains(b1.min))
		return true;
	if (b2.contains(b1.max))
		return true;

	std::array<vector3, 6> otherPoints;
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

inline transform2d::transform2d(vector2 position_, normalized<complex> rotation_, vector2 scale_) :
	position(position_),
	rotation(rotation_),
	scale(scale_)
{}

inline transform2d transform2d::operator*(const transform2d& t) const
{
	transform2d res;
	res.position = rotation * vector2(scale.x * t.position.x, scale.y * t.position.y) + position;
	res.rotation = rotation + t.rotation;
	res.scale = vector2(scale.x * t.scale.x, scale.y * t.scale.y);
	return res;
}
inline vector2 transform2d::operator*(vector2 v) const
{
	return rotation * vector2(scale.x * v.x, scale.y * v.y) + position;
}

inline transform3d::transform3d(vector3 position_, quaternion rotation_, vector3 scale_) :
	position(position_),
	rotation(rotation_),
	scale(scale_)
{}

inline transform3d transform3d::operator*(const transform3d& t) const
{
	transform3d res;
	res.position = rotation * vector3(scale.x * t.position.x, scale.y * t.position.y, scale.z * t.position.z) + position;
	res.rotation = rotation * t.rotation;
	res.scale = vector3(scale.x * t.scale.x, scale.y * t.scale.y, scale.z * t.scale.z);
	return res;
}
inline vector3 transform3d::operator*(vector3 v) const
{
	return rotation * vector3(scale.x * v.x, scale.y * v.y, scale.z * v.z) + position;
}
inline quaternion transform3d::operator*(quaternion q) const { return rotation * q; }

inline uniform_transform2d::uniform_transform2d(vector2 position_, radian rotation_, float scale_) :
	position(position_),
	rotation(rotation_),
	scale(scale_)
{}

inline uniform_transform2d uniform_transform2d::operator*(uniform_transform2d t) const
{
	uniform_transform2d res;
	matrix2 r = matrix2::rotation(rotation);
	res.position = r * (scale * t.position) + position;
	res.rotation = rotation + t.rotation;
	res.scale = scale * t.scale;
	return res;
}
inline vector2 uniform_transform2d::operator*(vector2 v) const
{
	matrix2 r = matrix2::rotation(rotation);
	return r * (scale * v) + position;
}

inline uniform_transform3d::uniform_transform3d(vector3 position_, quaternion rotation_, float scale_) :
	position(position_),
	rotation(rotation_),
	scale(scale_)
{}

inline uniform_transform3d uniform_transform3d::operator*(const uniform_transform3d& t) const
{
	uniform_transform3d res;
	res.position = rotation * (scale * t.position) + position;
	res.rotation = rotation * t.rotation;
	res.scale = scale * t.scale;
	return res;
}
inline vector3 uniform_transform3d::operator*(vector3 v) const { return rotation * (scale * v) + position; }
inline quaternion uniform_transform3d::operator*(quaternion q) const { return rotation * q; }

inline unscaled_transform2d::unscaled_transform2d(vector2 position_, radian rotation_) :
	position(position_),
	rotation(rotation_)
{}

inline unscaled_transform2d unscaled_transform2d::operator*(unscaled_transform2d t) const
{
	unscaled_transform2d res;
	matrix2 r = matrix2::rotation(rotation);
	res.position = r * t.position + position;
	res.rotation = rotation + t.rotation;
	return res;
}
inline vector2 unscaled_transform2d::operator*(vector2 v) const
{
	matrix2 r = matrix2::rotation(rotation);
	return r * v + position;
}

inline unscaled_transform3d::unscaled_transform3d(vector3 position_, quaternion rotation_) :
	position(position_),
	rotation(rotation_)
{}

inline unscaled_transform3d& unscaled_transform3d::translate_absolute(vector3 t)
{
	position += t;
	return *this;
}

inline unscaled_transform3d& unscaled_transform3d::translate_relative(vector3 t)
{
	position += rotation * t;
	return *this;
}

inline unscaled_transform3d& unscaled_transform3d::rotate(quaternion r)
{
	rotation = r * rotation;
	return *this;
}

inline unscaled_transform3d& unscaled_transform3d::inverse()
{
	rotation.inverse();
	position = rotation * -position;
	return *this;
}
inline unscaled_transform3d unscaled_transform3d::get_inversed() const
{
	unscaled_transform3d res = *this;
	res.inverse();
	return res;
}

inline matrix4 unscaled_transform3d::to_matrix() const
{
	return matrix4(*this);
}

inline unscaled_transform3d unscaled_transform3d::operator*(const unscaled_transform3d& t) const
{
	unscaled_transform3d res;
	res.position = rotation * t.position + position;
	res.rotation = rotation * t.rotation;
	return res;
}
inline vector3 unscaled_transform3d::operator*(vector3 v) const { return rotation * v + position; }
inline quaternion unscaled_transform3d::operator*(quaternion q) const { return rotation * q; }

inline std::ostream& operator<<(std::ostream& os, vector2 v)
{
	return os << "vector2(" << v.x << ", " << v.y << ")";
}
inline std::ostream& operator<<(std::ostream& os, vector3 v)
{
	return os << "vector3(" << v.x << ", " << v.y << ", " << v.z << ")";
}
inline std::ostream& operator<<(std::ostream& os, vector4 v)
{
	return os << "vector4(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")";
}
inline std::ostream& operator<<(std::ostream& os, quaternion q)
{
	return os << "quaternion({" << q.x << ", " << q.y << ", " << q.z << "}, " << q.w << ")";
}
inline std::ostream& operator<<(std::ostream& os, plane p)
{
	return os << "plane(origin = " << p.origin << ", normal = " << p.normal << ")";
}

inline std::ostream& operator<<(std::ostream& os, const cdm::matrix4& m)
{
	return os << "matrix4(" << m.m00 << "\t" << m.m10 << "\t" << m.m20 << "\t" << m.m30
		    << "\n        " << m.m01 << "\t" << m.m11 << "\t" << m.m21 << "\t" << m.m31
			<< "\n        " << m.m02 << "\t" << m.m12 << "\t" << m.m22 << "\t" << m.m32
			<< "\n        " << m.m03 << "\t" << m.m13 << "\t" << m.m23 << "\t" << m.m33
			<< ")";
}

namespace detail
{
template<unsigned char x, unsigned char y>
float get_quaternion_matrix_element(quaternion q)
{
	static_assert(x < 3, "x must be [0;2]");
	static_assert(y < 3, "y must be [0;2]");
	if constexpr (y == 0)
	{
		if constexpr (x == 0)
			return 1.0f - 2.0f * (q.y * q.y + q.z * q.z);
		if constexpr (x == 1)
			return 2.0f * (q.x * q.y - q.z * q.w);
		else
			return 2.0f * (q.x * q.z + q.y * q.w);
	}
	if constexpr (y == 1)
	{
		if constexpr (x == 0)
			return 2.0f * (q.x * q.y + q.z * q.w);
		if constexpr (x == 1)
			return 1.0f - 2.0f * (q.x * q.x + q.z * q.z);
		else
			return 2.0f * (q.y * q.z - q.x * q.w);
	}
	else
	{
		if constexpr (x == 0)
			return 2.0f * (q.x * q.z - q.y * q.w);
		if constexpr (x == 1)
			return 2.0f * (q.y * q.z + q.x * q.w);
		else
			return 1.0f - 2.0f * (q.x * q.x + q.y * q.y);
	}
}
} // namespace detail

template<typename Functor>
std::vector<vector3> function2D_sampler(const Functor& functor, float min, float max, float step)
{
	std::vector<vector3> res;

	if (min > max)
		std::swap(min, max);

	if (step < std::numeric_limits<float>::epsilon())
		step = std::numeric_limits<float>::epsilon();

	for (float f = min; f < max; f += step)
	{
		res.push_back(functor(f));
	}

	return res;
}

template<typename Functor>
std::vector<std::vector<vector3>> function3D_sampler(const Functor& functor, vector2 min, vector2 max, vector2 step)
{
	std::vector<std::vector<vector3>> res;

	if (min.x > max.x)
		std::swap(min.x, max.x);
	if (min.y > max.y)
		std::swap(min.y, max.y);

	if (step.x < std::numeric_limits<float>::epsilon())
		step.x = std::numeric_limits<float>::epsilon();
	if (step.y < std::numeric_limits<float>::epsilon())
		step.y = std::numeric_limits<float>::epsilon();

	size_t xCount = 0;
	size_t yCount = 0;

	for (float x = min.x; x < max.x; x += step.x)
		xCount++;
	for (float y = min.y; y < max.y; y += step.y)
		yCount++;
	res.reserve(yCount);

	for (float y = min.y; y < max.y; y += step.y)
	{
		std::vector<vector3> row;
		row.reserve(xCount);
		for (float x = min.x; x < max.x; x += step.x)
		{
			row.push_back(functor(x, y));
		}
		res.push_back(std::move(row));
	}

	return res;
}

namespace literals
{
inline cdm::radian operator""_rad(long double d) { return cdm::radian(static_cast<float>(d)); }
inline cdm::radian operator""_rad(unsigned long long int i) { return cdm::radian(static_cast<float>(i)); }
inline cdm::radian operator""_pi(long double d) { return cdm::radian(static_cast<float>(d) * cdm::pi); }
inline cdm::radian operator""_pi(unsigned long long int i) { return cdm::radian(static_cast<float>(i) * cdm::pi); }
inline cdm::degree operator""_deg(long double d) { return cdm::degree(static_cast<float>(d)); }
inline cdm::degree operator""_deg(unsigned long long int i) { return cdm::degree(static_cast<float>(i)); }
} // namespace literals
} // namespace cdm

#endif // CDM_MATHS_HPP
