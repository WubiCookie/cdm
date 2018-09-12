/* cdm_maths - v0.0 - geometric library - https://github.com/WubiCookie/cdm
   no warranty implied; use at your own risk

LICENSE

       DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
                   Version 2, December 2004

Copyright (C) 2018 Charles Seizilles de Mazancourt <charles DOT de DOT mazancourt AT hotmail DOT fr>

Everyone is permitted to copy and distribute verbatim or modified
copies of this license document, and changing it is allowed as long
as the name is changed.

           DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
  TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

 0. You just DO WHAT THE FUCK YOU WANT TO.

CREDITS

Written by Charles Seizilles de Mazancourt
*/

#ifndef CDM_MATH_HPP
#define CDM_MATH_HPP

#include <algorithm>
#include <cmath>
#include <utility>

#ifndef CDM_PI
#define CDM_PI (3.14159265f)
#endif

#ifndef CDM_FLT_EPSILON
#define CDM_FLT_EPSILON (1.192092896E-7f)
#endif

#ifndef CDM_DEG_TO_RAD
#define CDM_DEG_TO_RAD (0.01745329f)
#endif

#ifndef CDM_RAD_TO_DEG
#define CDM_RAD_TO_DEG (57.29577951f)
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
struct matrix4;
struct euler_angles;
struct quaternion;
struct line;
struct plane;
struct ray2d;
struct ray3d;
struct transform2d;
struct transform3d;
struct uniform_transform2d;
struct uniform_transform3d;
struct unscaled_transform2d;
struct unscaled_transform3d;

namespace detail
{
template<unsigned char x, unsigned char y>
float get_quaternion_matrix_element(const quaternion& q);
} // namespace detail

template<typename T>
T clamp(T f, T min_, T max_) { return std::min(max_, std::max(min_, f)); }

template<typename T>
T lerp(T begin, T end, T percent) { return (end - begin) * percent + begin; }

template<typename T>
bool nearly_equal_epsilon(T f1, T f2, T epsilon) { return fabs(f1 - f2) < epsilon; }

template<typename T>
bool nearly_equal(T f1, T f2) { return nearly_equal_epsilon(f1, f2, CDM_FLT_EPSILON); }

struct complex
{
	float r = 0.0f;
	float i = 0.0f;

	complex() = default;
	complex(const complex&) = default;
	complex(complex&&) = default;
	complex(float real, float imaginary);
	explicit complex(const radian& angle);

	float norm() const;
	float norm_squared() const;
	complex& normalize();
	complex get_normalized() const;
	complex& conjugate();
	complex get_conjugated() const;

	complex operator+(const complex& c) const;
	complex operator-(const complex& c) const;

	complex& operator+=(const complex& c);
	complex& operator-=(const complex& c);
	complex& operator*=(const complex& c);

	complex& operator=(const complex&) = default;
	complex& operator=(complex&&) = default;
};

complex operator*(const complex& c1, const complex& c2);
complex operator*(const complex& c, float f);
complex operator*(const normalized<complex>& c1, const complex& c2);
complex operator*(const complex& c1, const normalized<complex>& c2);
normalized<complex> operator*(const normalized<complex>& c1, const normalized<complex>& c2);
vector2 operator*(const normalized<complex>& c, const vector2& v);

struct radian
{
	float angle = 0.0f;

	radian() = default;
	explicit radian(float f);
	radian(const radian& r);
	radian(const degree& d);

	operator float();
	operator float() const;

	radian& operator+=(float f);
	radian& operator-=(float f);
	radian& operator*=(float f);
	radian& operator/=(float f);

	radian operator-();
};

radian operator+(const radian&, const radian&);
radian operator-(const radian&, const radian&);
radian operator*(const radian&, const radian&);
radian operator/(const radian&, const radian&);
radian operator*(const radian&, float);
radian operator*(float, const radian&);

struct degree
{
	float angle = 0.0f;

	degree() = default;
	explicit degree(float f);
	degree(const degree& d);
	degree(const radian& r);

	operator float();
	operator float() const;

	degree& operator+=(float f);
	degree& operator-=(float f);
	degree& operator*=(float f);
	degree& operator/=(float f);

	degree operator-();
};

degree operator+(const degree&, const degree&);
degree operator-(const degree&, const degree&);
degree operator*(const degree&, const degree&);
degree operator/(const degree&, const degree&);
degree operator*(const degree&, float);
degree operator*(float, const degree&);

radian operator""_rad(long double);
radian operator""_rad(unsigned long long int);
radian operator""_pi(long double);
radian operator""_pi(unsigned long long int);
degree operator""_deg(long double);
degree operator""_deg(unsigned long long int);

struct vector2
{
	float x = 0.0f, y = 0.0f;

	vector2() = default;
	vector2(float x, float y);

	float norm() const;
	float norm_squared() const;
	vector2& normalize();
	vector2 get_normalized() const;
	vector2& clamp(const vector2& min, const vector2& max);
	vector2 get_clamped(const vector2& min, const vector2& max) const;
	vector2& negate();
	vector2 get_negated() const;
	float dot(const vector2& v) const;
	static vector2 lerp(const vector2& begin, const vector2& end, float percent);
	static vector2 nlerp(const vector2& begin, const vector2& end, float percent);
	static vector2 slerp(const vector2& begin, const vector2& end, float percent);
	float distance_from(const vector2& v) const;
	static float distance_between(const vector2& v1, const vector2& v2);
	static vector2 from_to(const vector2& from, const vector2& to);
	vector2 to(const vector2& v) const;
	static radian angle(const vector2& v1, const vector2& v2);
	static bool nearly_equal(const vector2& v1, const vector2& v2);

	vector2 operator+(const vector2& v) const;
	vector2 operator-(const vector2& v) const;
	vector2 operator*(float f) const;
	vector2 operator/(float f) const;

	vector2& operator+=(const vector2& v);
	vector2& operator-=(const vector2& v);
	vector2& operator*=(float f);
	vector2& operator/=(float f);

	vector2 operator-() const;

	bool operator==(const vector2& v) const;
	bool operator!=(const vector2& v) const;
};

vector2 operator*(float f, const vector2& v);

struct vector3
{
	float x = 0.0f, y = 0.0f, z = 0.0f;

	vector3() = default;
	vector3(float x, float y, float z);

	float norm() const;
	float norm_squared() const;
	vector3& normalize();
	vector3 get_normalized() const;
	vector3& clamp(const vector3& min, const vector3& max);
	vector3 get_clamped(const vector3& min, const vector3& max) const;
	vector3& negate();
	vector3 get_negated() const;
	float dot(const vector3& v) const;
	vector3 cross(const vector3& v) const;
	static vector3 lerp(const vector3& begin, const vector3& end, float percent);
	static vector3 nlerp(const vector3& begin, const vector3& end, float percent);
	float distance_from(const vector3& v) const;
	static float distance_between(const vector3& v1, const vector3& v2);
	static vector3 from_to(const vector3& from, const vector3& to);
	vector3 to(const vector3& v) const;
	radian angle(const vector3& v) const;
	radian angle_around_axis(const vector3& v, const vector3& axis);

	vector3 operator+(const vector3& v) const;
	vector3 operator-(const vector3& v) const;
	vector3 operator*(float f) const;
	vector3 operator/(float f) const;

	vector3& operator+=(const vector3& v);
	vector3& operator-=(const vector3& v);
	vector3& operator*=(float f);
	vector3& operator/=(float f);

	vector3 operator-() const;

	bool operator==(const vector3& v) const;
	bool operator!=(const vector3& v) const;
};

vector3 operator*(float f, const vector3& v);

struct vector4
{
	float x = 0.0f, y = 0.0f, z = 0.0f, w = 0.0f;

	vector4() = default;
	vector4(float x, float y, float z, float w);

	float norm() const;
	float norm_squared() const;
	vector4& normalize();
	vector4 get_normalized() const;
	vector4& clamp(const vector4& min, const vector4& max);
	vector4 get_clamped(const vector4& min, const vector4& max) const;
	vector4& negate();
	vector4 get_negated() const;
	float dot(const vector4& v) const;
	static vector4 lerp(const vector4& begin, const vector4& end, float percent);
	static vector4 nlerp(const vector4& begin, const vector4& end, float percent);
	float distance_from(const vector4& v) const;
	static float distance_between(const vector4& v1, const vector4& v2);
	static vector4 from_to(const vector4& from, const vector4& to);
	vector4 to(const vector4& v) const;
	static radian angle(const vector4& v1, const vector4& v2);

	vector4 operator+(const vector4& v) const;
	vector4 operator-(const vector4& v) const;
	vector4 operator*(float f) const;
	vector4 operator/(float f) const;

	vector4& operator+=(const vector4& v);
	vector4& operator-=(const vector4& v);
	vector4& operator*=(float f);
	vector4& operator/=(float f);

	vector4 operator-() const;

	bool operator==(const vector4& v) const;
	bool operator!=(const vector4& v) const;
};

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
	matrix2(const matrix3& m);
	matrix2(const matrix4& m);

	static matrix2 zero(void);
	static matrix2 identity(void);
	static matrix2 rotation(const radian& angle);
	static matrix2 rotation(const normalized<complex>& angle);
	matrix2& transpose();
	matrix2 get_transposed() const;
	float get_determinant() const;

	float& at(uint8_t x, uint8_t y);
	const float& at(uint8_t x, uint8_t y) const;

	vector2 operator*(const vector2& v) const;
	matrix2 operator*(const matrix2& m) const;
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

	matrix3(const matrix2& m);
	matrix3(const matrix4& m);
	matrix3(const transform2d& t);
	matrix3(const uniform_transform2d& t);
	matrix3(const unscaled_transform2d& t);

	static matrix3 zero();
	static matrix3 identity();
	static matrix3 rotation(const euler_angles& r);
	static matrix3 rotation(const quaternion& q);
	static matrix3 rotation_around_x(const radian& angle);
	static matrix3 rotation_around_y(const radian& angle);
	static matrix3 rotation_around_z(const radian& angle);
	static matrix3 rotation_around_x(const normalized<complex>& angle);
	static matrix3 rotation_around_y(const normalized<complex>& angle);
	static matrix3 rotation_around_z(const normalized<complex>& angle);
	matrix3& inverse();
	matrix3 get_inversed() const;
	matrix3& transpose();
	matrix3 get_transposed() const;
	float get_determinant() const;
	bool is_orthogonal() const;

	float& at(uint8_t x, uint8_t y);
	const float& at(uint8_t x, uint8_t y) const;

	matrix3 operator*(float f) const;
	vector3 operator*(const vector3& v) const;
	matrix3 operator*(const matrix3& m) const;
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
	matrix4(const matrix2& m);
	matrix4(const matrix3& m);
	matrix4(const transform3d& t);
	matrix4(const uniform_transform3d& t);
	matrix4(const unscaled_transform3d& t);

	static matrix4 zero();
	static matrix4 identity();
	static matrix4 rotation(const euler_angles& r);
	static matrix4 rotation(const quaternion& q);
	static matrix4 translation(const vector3& t);
	static matrix4 perspective(const radian& angle, float ratio, float near_plane, float far_plane);
	bool is_orthogonal() const;
	bool is_homogenous() const;
	matrix4& inverse();
	matrix4 get_inversed() const;
	matrix4& transpose();
	matrix4 get_transposed() const;
	float get_determinant() const;

	float& at(uint8_t x, uint8_t y);
	const float& at(uint8_t x, uint8_t y) const;

	vector4 operator*(const vector4& v) const;
	matrix4 operator*(const matrix4& m) const;
};

struct euler_angles
{
	radian x, y, z;

	euler_angles() = default;
	euler_angles(const radian& x, const radian& y, const radian& z);
};

struct quaternion
{
	float x = 0.0f, y = 0.0f, z = 0.0f, w = 1.0f;

	quaternion() = default;
	quaternion(float x, float y, float z, float w);
	quaternion(const normalized<vector3>& axis, const radian& angle);

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
	quaternion& clamp(const quaternion& min, const quaternion& max);
	quaternion get_clamped(const quaternion& min, const quaternion& max) const;
	float norm() const;
	float norm_squared() const;
	static quaternion lerp(const quaternion& begin, const quaternion& end, float percent);
	static quaternion nlerp(const quaternion& begin, const quaternion& end, float percent);
	static quaternion slerp(const quaternion& begin, const quaternion& end, float percent);
	float dot(const quaternion&) const;
	quaternion cross(const quaternion&) const;

	quaternion operator+(const quaternion& q) const;
	quaternion operator-(const quaternion& q) const;
	quaternion operator*(const quaternion& q) const; // TODO: implement operator* for normalized<> if result is normalized too
	quaternion operator*(float f) const;
	quaternion operator/(float f) const;

	quaternion& operator+=(const quaternion& q);
	quaternion& operator-=(const quaternion& q);
	quaternion& operator*=(float f);
	quaternion& operator/=(float f);

	quaternion operator-() const;

	bool operator==(const quaternion& v) const;
	bool operator!=(const quaternion& v) const;
};

vector3 operator*(const normalized<quaternion>& q, const vector3& v);

struct line
{
	vector2 normal;
	float distance = 0.0f;
};

struct plane
{
	vector3 normal;
	float distance = 0.0f;
};

struct ray2d
{
	vector2 origin;
	vector2 direction;
};

struct ray3d
{
	vector3 origin;
	vector3 direction;
};

struct transform2d
{
	vector2 position;
	normalized<complex> rotation;
	vector2 scale = {1.0f, 1.0f};

	transform2d() = default;
	transform2d(const vector2& position, const normalized<complex>& rotation, const vector2& scale);

	transform2d operator*(const transform2d& t) const;
	vector2 operator*(const vector2& v) const;
};

struct transform3d
{
	vector3 position;
	quaternion rotation;
	vector3 scale = {1.0f, 1.0f, 1.0f};

	transform3d() = default;
	transform3d(const vector3& position, const quaternion& rotation, const vector3& scale);

	transform3d operator*(const transform3d& t) const;
	vector3 operator*(const vector3& v) const;
	quaternion operator*(const quaternion& q) const;
};

struct uniform_transform2d
{
	vector2 position;
	radian rotation;
	float scale = 1.0f;

	uniform_transform2d() = default;
	uniform_transform2d(const vector2& position, const radian& rotation, float scale);

	uniform_transform2d operator*(const uniform_transform2d& t) const;
	vector2 operator*(const vector2& v) const;
};

struct uniform_transform3d
{
	vector3 position;
	quaternion rotation;
	float scale = 1.0f;

	uniform_transform3d() = default;
	uniform_transform3d(const vector3& position, const quaternion& rotation, float scale);

	uniform_transform3d operator*(const uniform_transform3d& t) const;
	vector3 operator*(const vector3& v) const;
	quaternion operator*(const quaternion& q) const;
};

struct unscaled_transform2d
{
	vector2 position;
	radian rotation;

	unscaled_transform2d() = default;
	unscaled_transform2d(const vector2& position, const radian& rotation);

	unscaled_transform2d operator*(const unscaled_transform2d& t) const;
	vector2 operator*(const vector2& v) const;
};

struct unscaled_transform3d
{
	vector3 position;
	quaternion rotation;

	unscaled_transform3d() = default;
	unscaled_transform3d(const vector3& position, const quaternion& rotation);

	unscaled_transform3d operator*(const unscaled_transform3d& t) const;
	vector3 operator*(const vector3& v) const;
	quaternion operator*(const quaternion& q) const;
};

inline complex::complex(float real, float imaginary) : r(real), i(imaginary) {}
inline complex::complex(const radian& angle) : r(cosf(angle)), i(sinf(angle)) {}

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

inline complex complex::operator+(const complex& c) const
{
	return {
		r + c.r,
		i + c.i
	};
}
inline complex complex::operator-(const complex& c) const
{
	return {
		r - c.r,
		i - c.i
	};
}

inline complex& complex::operator+=(const complex& c) { return *this = *this + c; }
inline complex& complex::operator-=(const complex& c) { return *this = *this - c; }
inline complex& complex::operator*=(const complex& c) { return *this = *this * c; }

inline complex operator*(const complex& c1, const complex& c2)
{
	return {
		c1.r * c2.r - c1.i * c2.i,
		c1.r * c2.i + c1.i * c2.r
	};
}
inline complex operator*(const complex& c, float f) { return {c.r * f, c.i * f}; }
inline complex operator*(const normalized<complex>& c1, const complex& c2) { return *c1 * c2; }
inline complex operator*(const complex& c1, const normalized<complex>& c2) { return c1 * *c2; }
inline normalized<complex> operator*(const normalized<complex>& c1, const normalized<complex>& c2)
{
	return normalized<complex>::already_normalized(complex(
		c1->r * c2->r - c1->i * c2->i,
		c1->r * c2->i + c1->i * c2->r
	));
}
inline vector2 operator*(const normalized<complex>& c, const vector2& v)
{
	return {
		c->r * v.x - c->i * v.y,
		c->r * v.y + c->i * v.x
	};
}

inline radian::radian(float f) : angle(f) {}
inline radian::radian(const radian& r) : angle(r.angle) {}
inline radian::radian(const degree& d) : angle(d.angle * CDM_DEG_TO_RAD) {}

inline radian::operator float() { return angle; }
inline radian::operator float() const { return angle; }

inline radian& radian::operator+=(float f) { angle += f; return *this; }
inline radian& radian::operator-=(float f) { angle -= f; return *this; }
inline radian& radian::operator*=(float f) { angle *= f; return *this; }
inline radian& radian::operator/=(float f) { angle /= f; return *this; }

inline radian radian::operator-() { return radian(-angle); }

inline degree::degree(float f) : angle(f) {}
inline degree::degree(const degree& d) : angle(d.angle) {}
inline degree::degree(const radian& r) : angle(r.angle * CDM_RAD_TO_DEG) {}

inline degree::operator float() { return angle; }
inline degree::operator float() const { return angle; }

inline degree& degree::operator+=(float f) { angle += f; return *this; }
inline degree& degree::operator-=(float f) { angle -= f; return *this; }
inline degree& degree::operator*=(float f) { angle *= f; return *this; }
inline degree& degree::operator/=(float f) { angle /= f; return *this; }

inline degree degree::operator-() { return radian(-angle); }

inline radian operator+(const radian& r1, const radian& r2) { return radian(r1.angle + r2.angle); }
inline radian operator-(const radian& r1, const radian& r2) { return radian(r1.angle - r2.angle); }
inline radian operator*(const radian& r1, const radian& r2) { return radian(r1.angle * r2.angle); }
inline radian operator/(const radian& r1, const radian& r2) { return radian(r1.angle / r2.angle); }
inline radian operator*(const radian& r, float f) { return radian(r.angle * f); }
inline radian operator*(float f, const radian& r) { return radian(r.angle * f); }

inline degree operator+(const degree& d1, const degree& d2) { return degree(d1.angle + d2.angle); }
inline degree operator-(const degree& d1, const degree& d2) { return degree(d1.angle - d2.angle); }
inline degree operator*(const degree& d1, const degree& d2) { return degree(d1.angle * d2.angle); }
inline degree operator/(const degree& d1, const degree& d2) { return degree(d1.angle / d2.angle); }
inline degree operator*(const degree& d, float f) { return degree(d.angle * f); }
inline degree operator*(float f, const degree& d) { return degree(d.angle * f); }

inline radian operator""_rad(long double d) { return radian(d); }
inline radian operator""_rad(unsigned long long int i) { return radian(i); }
inline radian operator""_pi(long double d) { return radian(d * CDM_PI); }
inline radian operator""_pi(unsigned long long int i) { return radian(i * CDM_PI); }
inline degree operator""_deg(long double d) { return degree(d); }
inline degree operator""_deg(unsigned long long int i) { return degree(i); }

inline vector2::vector2(float x, float y) : x(x), y(y) {}

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
inline vector2& vector2::clamp(const vector2& min, const vector2& max)
{
	x = cdm::clamp(x, min.x, max.x);
	y = cdm::clamp(y, min.y, max.y);
	return *this;
}
inline vector2 vector2::get_clamped(const vector2& min, const vector2& max) const
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
inline float vector2::dot(const vector2& v) const { return x * v.x + y * v.y; }
inline vector2 vector2::lerp(const vector2& begin, const vector2& end, float percent) {	return (end - begin) * percent + begin; }
inline vector2 vector2::nlerp(const vector2& begin, const vector2& end, float percent) { return lerp(begin, end, percent).get_normalized(); }
inline vector2 vector2::slerp(const vector2& begin, const vector2& end, float percent)
{
	float angle = vector2::angle(begin, end) * percent;
	float s = sinf(angle);
	float c = cosf(angle);

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
inline float vector2::distance_from(const vector2& v) const { return distance_between(*this, v); }
inline float vector2::distance_between(const vector2& v1, const vector2& v2) { return (v1 - v2).norm(); }
inline vector2 vector2::from_to(const vector2& from, const vector2& to) { return {to.x - from.x, to.y - from.y}; }
inline vector2 vector2::to(const vector2& v) const { return from_to(*this, v); }
inline radian vector2::angle(const vector2& v1, const vector2& v2) { return radian(atan2f(v2.y, v2.x) - atan2f(v1.y, v1.x)); }
inline bool vector2::nearly_equal(const vector2& v1, const vector2& v2) { return cdm::nearly_equal(v1.x, v2.x) && cdm::nearly_equal(v1.y, v2.y); }

inline vector2 vector2::operator+(const vector2& v) const { return {x + v.x, y + v.y}; }
inline vector2 vector2::operator-(const vector2& v) const { return {x - v.x, y - v.y}; }
inline vector2 vector2::operator*(float f) const { return {x * f, y * f}; }
inline vector2 vector2::operator/(float f) const { return {x / f, y / f}; }

inline vector2& vector2::operator+=(const vector2& v) { *this = *this + v; return *this; }
inline vector2& vector2::operator-=(const vector2& v) { *this = *this - v; return *this; }
inline vector2& vector2::operator*=(float f) { *this = *this * f; return *this; }
inline vector2& vector2::operator/=(float f) { *this = *this / f; return *this; }

inline vector2 vector2::operator-() const { return get_negated(); }

inline bool vector2::operator==(const vector2& v) const { return x == v.x && y == v.y; }
inline bool vector2::operator!=(const vector2& v) const { return !operator==(v); }

inline vector2 operator*(float f, const vector2& v) { return v * f; }

inline vector3::vector3(float x, float y, float z) : x(x), y(y), z(z) {}

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
inline vector3& vector3::clamp(const vector3& min, const vector3& max)
{
	x = cdm::clamp(x, min.x, max.x);
	y = cdm::clamp(y, min.y, max.y);
	z = cdm::clamp(z, min.z, max.z);
	return *this;
}
inline vector3 vector3::get_clamped(const vector3& min, const vector3& max) const
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
inline float vector3::dot(const vector3& v) const { return x * v.x + y * v.y + z * v.z; }
inline vector3 vector3::cross(const vector3& v) const
{
	return {
		y * v.z - z * v.y,
		z * v.x - x * v.z,
		x * v.y - y * v.x
	};
}
inline vector3 vector3::lerp(const vector3& begin, const vector3& end, float percent) {	return (end - begin) * percent + begin; }
inline vector3 vector3::nlerp(const vector3& begin, const vector3& end, float percent) { return lerp(begin, end, percent).get_normalized(); }
inline float vector3::distance_from(const vector3& v) const { return distance_between(*this, v); }
inline float vector3::distance_between(const vector3& v1, const vector3& v2) { return (v1 - v2).norm(); }
inline vector3 vector3::from_to(const vector3& from, const vector3& to) { return {to.x - from.x, to.y - from.y, to.z - from.z}; }
inline vector3 vector3::to(const vector3& v) const { return from_to(*this, v); }
inline radian vector3::angle(const vector3& v) const
{
	float divisor = sqrtf(norm_squared() * v.norm_squared());
	float alpha = dot(v) / divisor;
	return radian(acosf(cdm::clamp(alpha, -1.0f, 1.0f)));
}
inline radian vector3::angle_around_axis(const vector3& v, const vector3& axis)
{
	vector3 c = cross(v);
	float angle = atan2f(c.norm(), dot(v));
	return radian(c.dot(axis) < 0.0f ? -angle : angle);
}

inline vector3 vector3::operator+(const vector3& v) const { return {x + v.x, y + v.y, z + v.z}; }
inline vector3 vector3::operator-(const vector3& v) const { return {x - v.x, y - v.y, z - v.z}; }
inline vector3 vector3::operator*(float f) const { return {x * f, y * f, z * f}; }
inline vector3 vector3::operator/(float f) const { return {x / f, y / f, z / f}; }

inline vector3& vector3::operator+=(const vector3& v) { *this = *this + v; return *this; }
inline vector3& vector3::operator-=(const vector3& v) { *this = *this - v; return *this; }
inline vector3& vector3::operator*=(float f) { *this = *this * f; return *this; }
inline vector3& vector3::operator/=(float f) { *this = *this / f; return *this; }

inline vector3 vector3::operator-() const { return get_negated(); }

inline bool vector3::operator==(const vector3& v) const { return x == v.x && y == v.y && z == v.z; }
inline bool vector3::operator!=(const vector3& v) const { return !operator==(v); }

inline vector3 operator*(float f, const vector3& v) { return v * f; }

inline vector4::vector4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}

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
inline vector4& vector4::clamp(const vector4& min, const vector4& max)
{
	x = cdm::clamp(x, min.x, max.x);
	y = cdm::clamp(y, min.y, max.y);
	z = cdm::clamp(z, min.z, max.z);
	w = cdm::clamp(w, min.w, max.w);
	return *this;
}
inline vector4 vector4::get_clamped(const vector4& min, const vector4& max) const
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
inline float vector4::dot(const vector4& v) const { return x * v.x + y * v.y + z * v.z + w * v.w; }
inline vector4 vector4::lerp(const vector4& begin, const vector4& end, float percent) {	return (end - begin) * percent + begin; }
inline vector4 vector4::nlerp(const vector4& begin, const vector4& end, float percent) { return lerp(begin, end, percent).get_normalized(); }
inline float vector4::distance_from(const vector4& v) const { return distance_between(*this, v); }
inline float vector4::distance_between(const vector4& v1, const vector4& v2) { return (v1 - v2).norm(); }
inline vector4 vector4::from_to(const vector4& from, const vector4& to) { return {to.x - from.x, to.y - from.y, to.z - from.z, to.w - from.w}; }
inline vector4 vector4::to(const vector4& v) const { return from_to(*this, v); }

inline vector4 vector4::operator+(const vector4& v) const { return {x + v.x, y + v.y, z + v.z, w + v.w}; }
inline vector4 vector4::operator-(const vector4& v) const { return {x - v.x, y - v.y, z - v.z, w - v.w}; }
inline vector4 vector4::operator*(float f) const { return {x * f, y * f, z * f, w * f}; }
inline vector4 vector4::operator/(float f) const { return {x / f, y / f, z / f, w / f}; }

inline vector4& vector4::operator+=(const vector4& v) { *this = *this + v; return *this; }
inline vector4& vector4::operator-=(const vector4& v) { *this = *this - v; return *this; }
inline vector4& vector4::operator*=(float f) { *this = *this * f; return *this; }
inline vector4& vector4::operator/=(float f) { *this = *this / f; return *this; }

inline vector4 vector4::operator-() const { return get_negated(); }

inline bool vector4::operator==(const vector4& v) const { return x == v.x && y == v.y && z == v.z && w == v.w; }
inline bool vector4::operator!=(const vector4& v) const { return !operator==(v); }

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

inline matrix2::matrix2(float m00, float m10,
	float m01, float m11) :
	m00(m00), m10(m10),
	m01(m01), m11(m11)
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
inline matrix2 matrix2::rotation(const radian& angle)
{
	float c = cosf(angle);
	float s = sinf(angle);
	return {
		c, -s,
		s, c
	};
}
inline matrix2 matrix2::rotation(const normalized<complex>& angle)
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

inline vector2 matrix2::operator*(const vector2& v) const
{
	return {
		m00 * v.x + m01 * v.y,
		m10 * v.x + m11 * v.y
	};
}
inline matrix2 matrix2::operator*(const matrix2& m) const
{
	return {
		m00 * m.m00 + m01 * m.m10,
		m00 * m.m01 + m01 * m.m11,

		m10 * m.m00 + m11 * m.m10,
		m10 * m.m01 + m11 * m.m11
	};
}

inline matrix3::matrix3(float m00, float m10, float m20,
	float m01, float m11, float m21,
	float m02, float m12, float m22) :
	m00(m00), m10(m10), m20(m20),
	m01(m01), m11(m11), m21(m21),
	m02(m02), m12(m12), m22(m22)
{}
inline matrix3::matrix3(const matrix2& m) :
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
inline matrix3::matrix3(const uniform_transform2d& t)
{
	matrix2 r = matrix2::rotation(t.rotation);
	*this = {
		t.scale * r.m00, t.scale * r.m10, t.position.x,
		t.scale * r.m01, t.scale * r.m11, t.position.y,
		0.0f,            0.0f,            1.0f
	};
}
inline matrix3::matrix3(const unscaled_transform2d& t)
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
inline matrix3 matrix3::rotation(const euler_angles& r)
{
	matrix3 RX = matrix3::identity();
	RX.m11 = cosf(r.x);
	RX.m12 = -sinf(r.x);
	RX.m21 = sinf(r.x);
	RX.m22 = cosf(r.x);

	matrix3 RY = matrix3::identity();
	RY.m00 = cosf(r.y);
	RY.m02 = sinf(r.y);
	RY.m20 = -sinf(r.y);
	RY.m22 = cosf(r.y);

	matrix3 RZ = matrix3::identity();
	RZ.m00 = cosf(r.z);
	RZ.m01 = -sinf(r.z);
	RZ.m10 = sinf(r.z);
	RZ.m11 = cosf(r.z);

	//return RZ * RX * RY;
	//return RZ * RY * RX;
	return RY * RX * RZ;
	//return RY * RZ * RX;
	//return RY * RX * RZ;
	//return RX * RZ * RY;
}
inline matrix3 matrix3::rotation(const quaternion& q)
{
	return {
		detail::get_quaternion_matrix_element<0, 0>(q), detail::get_quaternion_matrix_element<1, 0>(q), detail::get_quaternion_matrix_element<2, 0>(q),
		detail::get_quaternion_matrix_element<0, 1>(q), detail::get_quaternion_matrix_element<1, 1>(q), detail::get_quaternion_matrix_element<2, 1>(q),
		detail::get_quaternion_matrix_element<0, 2>(q), detail::get_quaternion_matrix_element<1, 2>(q), detail::get_quaternion_matrix_element<2, 2>(q),
	};
}
inline matrix3 matrix3::rotation_around_x(const radian& angle)
{
	float c = cosf(angle);
	float s = sinf(angle);
	return {
		1.0f, 0.0f, 0.0f,
		0.0f, c,    s,
		0.0f, -s,   c
	};
}
inline matrix3 matrix3::rotation_around_y(const radian& angle)
{
	float c = cosf(angle);
	float s = sinf(angle);
	return {
		c,    0.0f, s,
		0.0f, 1.0f, 0.0f,
		-s,   0.0f, c
	};
}
inline matrix3 matrix3::rotation_around_z(const radian& angle)
{
	float c = cosf(angle);
	float s = sinf(angle);
	return {
		c,    -s,   0.0f,
		s,    c,    0.0f,
		0.0f, 0.0f, 1.0f
	};
}
inline matrix3 matrix3::rotation_around_x(const normalized<complex>& angle)
{
	float c = angle->r;
	float s = angle->i;
	return {
		1.0f, 0.0f, 0.0f,
		0.0f, c,    s,
		0.0f, -s,   c
	};
}
inline matrix3 matrix3::rotation_around_y(const normalized<complex>& angle)
{
	float c = angle->r;
	float s = angle->i;
	return {
		c,    0.0f, s,
		0.0f, 1.0f, 0.0f,
		-s,   0.0f, c
	};
}
inline matrix3 matrix3::rotation_around_z(const normalized<complex>& angle)
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
inline vector3 matrix3::operator*(const vector3& v) const
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

inline matrix4::matrix4(float m00, float m10, float m20, float m30,
	float m01, float m11, float m21, float m31,
	float m02, float m12, float m22, float m32,
	float m03, float m13, float m23, float m33) :
	m00(m00), m10(m10), m20(m20), m30(m30),
	m01(m01), m11(m11), m21(m21), m31(m31),
	m02(m02), m12(m12), m22(m22), m32(m32),
	m03(m03), m13(m13), m23(m23), m33(m33)
{}
inline matrix4::matrix4(const matrix2& m) :
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

inline matrix4 matrix4::zero() { return {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}; }
inline matrix4 matrix4::identity() { return {1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f}; }
inline matrix4 matrix4::rotation(const euler_angles& r) { return matrix3::rotation(r); }
inline matrix4 matrix4::rotation(const quaternion& q) { return matrix3::rotation(q); }
inline matrix4 matrix4::translation(const vector3& t) { return {1.0f, 0.0f, 0.0f, t.x, 0.0f, 1.0f, 0.0f, t.y, 0.0f, 0.0f, 1.0f, t.z, 0.0f, 0.0f, 0.0f, 1.0f}; }
inline matrix4 matrix4::perspective(const radian& angle, float ratio, float near_plane, float far_plane)
{
	const float invTanHalfFovy = 1.0f / tanf(angle * 0.5f);
	const float k = far_plane / (2.0f * (near_plane - far_plane));
	return {
		invTanHalfFovy / ratio, 0.0f,           0.0f,  0.0f,
		0.0f,                   invTanHalfFovy, 0.0f,  0.0f,
		0.0f,                   0.0f,           k,     k * (1.0f + 2.0f * near_plane),
		0.0f,                   0.0f,           -1.0f, 0.0f
	};
}
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

inline vector4 matrix4::operator*(const vector4& v) const
{
	return {
		m00 * v.x + m10 * v.y + m20 * v.z + m20 * v.w,
		m01 * v.x + m11 * v.y + m21 * v.z + m21 * v.w,
		m02 * v.x + m12 * v.y + m22 * v.z + m22 * v.w,
		m03 * v.x + m13 * v.y + m23 * v.z + m23 * v.w
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

inline euler_angles::euler_angles(const radian& x, const radian& y, const radian& z) : x(x), y(y), z(z) {}

inline quaternion::quaternion(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}
inline quaternion::quaternion(const normalized<vector3>& axis, const radian& angle)
{
	/*vector3 tmpAxis = vector3_get_normalized(axis);
	tmpAxis = vector3_times_float(&tmpAxis, sinf(angle / 2.0f));*/

	vector3 tmpAxis = axis * sinf(angle / 2.0f);
	*this = {tmpAxis.x, tmpAxis.y, tmpAxis.z, cosf(angle / 2.0f)};
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
inline quaternion& quaternion::clamp(const quaternion& min, const quaternion& max)
{
	x = cdm::clamp(x, min.x, max.x);
	y = cdm::clamp(y, min.y, max.y);
	z = cdm::clamp(z, min.z, max.z);
	w = cdm::clamp(w, min.w, max.w);
	return *this;
}
inline quaternion quaternion::get_clamped(const quaternion& min, const quaternion& max) const
{
	quaternion res = *this;
	res.clamp(min, max);
	return res;
}
inline float quaternion::norm() const { return sqrtf(norm_squared()); }
inline float quaternion::norm_squared() const { return x * x + y * y + z * z + w * w; }
inline quaternion quaternion::lerp(const quaternion& begin, const quaternion& end, float percent) { return (end - begin) * percent + begin; }
inline quaternion quaternion::nlerp(const quaternion& begin, const quaternion& end, float percent) { return lerp(begin, end, percent).get_normalized(); }
inline quaternion quaternion::slerp(const quaternion& begin, const quaternion& end, float percent)
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
inline float quaternion::dot(const quaternion& q) const { return x * q.x + y * q.y + z * q.z + w * q.w; }
inline quaternion quaternion::cross(const quaternion& q) const
{
	return {
		y * q.z - z * q.y,
		z * q.x - x * q.z,
		x * q.y - y * q.x,
		1.0f
	};
}

inline quaternion quaternion::operator+(const quaternion& q) const { return {x + q.x, y + q.y, z + q.z, w + q.w}; }
inline quaternion quaternion::operator-(const quaternion& q) const { return {x - q.x, y - q.y, z - q.z, w - q.w}; }
inline quaternion quaternion::operator*(const quaternion& q) const
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
inline vector3 operator*(const normalized<quaternion>& q, const vector3& v)
{
	vector3 qvec = {q->x, q->y, q->z};
	vector3 uv = qvec.cross(v);
	vector3 uuv = qvec.cross(uv);
	uv = uv * (2.0f * q->w);
	uuv = uuv * 2.0f;
	return v + uv + uuv;
}

inline quaternion& quaternion::operator+=(const quaternion& q) { *this = *this + q; return *this; }
inline quaternion& quaternion::operator-=(const quaternion& q) { *this = *this - q; return *this; }
inline quaternion& quaternion::operator*=(float f) { *this = *this * f; return *this; }
inline quaternion& quaternion::operator/=(float f) { *this = *this / f; return *this; }

inline quaternion quaternion::operator-() const { return get_negated(); }

inline bool quaternion::operator==(const quaternion& q) const { return x == q.x && y == q.y && z == q.z && w == q.w; }
inline bool quaternion::operator!=(const quaternion& q) const { return !operator==(q); }

inline transform2d::transform2d(const vector2& position, const normalized<complex>& rotation, const vector2& scale) :
	position(position),
	rotation(rotation),
	scale(scale)
{}

inline transform2d transform2d::operator*(const transform2d& t) const
{
	transform2d res;
	res.position = rotation * vector2(scale.x * t.position.x, scale.y * t.position.y) + position;
	res.rotation = rotation + t.rotation;
	res.scale = vector2(scale.x * t.scale.x, scale.y * t.scale.y);
	return res;
}
inline vector2 transform2d::operator*(const vector2& v) const
{
	return rotation * vector2(scale.x * v.x, scale.y * v.y) + position;
}

inline transform3d::transform3d(const vector3& position, const quaternion& rotation, const vector3& scale) :
	position(position),
	rotation(rotation),
	scale(scale)
{}

inline transform3d transform3d::operator*(const transform3d& t) const
{
	transform3d res;
	// throw "transform3d transform3d::operator*(const transform3d& t) const not implemented";
	res.position = rotation * vector3(scale.x * t.position.x, scale.y * t.position.y, scale.z * t.position.z) + position;
	res.rotation = rotation * t.rotation;
	res.scale = vector3(scale.x * t.scale.x, scale.y * t.scale.y, scale.z * t.scale.z);
	return res;
}
inline vector3 transform3d::operator*(const vector3& v) const
{
	// throw "vector3 transform3d::operator*(const vector3& v) const not implemented";
	return rotation * vector3(scale.x * v.x, scale.y * v.y, scale.z * v.z) + position;
}
inline quaternion transform3d::operator*(const quaternion& q) const { return rotation * q; }

inline uniform_transform2d::uniform_transform2d(const vector2& position, const radian& rotation, float scale) :
	position(position),
	rotation(rotation),
	scale(scale)
{}

inline uniform_transform2d uniform_transform2d::operator*(const uniform_transform2d& t) const
{
	uniform_transform2d res;
	// throw "uniform_transform2d uniform_transform2d::operator*(const uniform_transform2d& t) const not implemented";
	matrix2 r = matrix2::rotation(rotation);
	res.position = r * (scale * t.position) + position;
	res.rotation = rotation + t.rotation;
	res.scale = scale * t.scale;
	return res;
}
inline vector2 uniform_transform2d::operator*(const vector2& v) const
{
	// throw "vector2 uniform_transform2d::operator*(const vector2& v) const not implemented";
	matrix2 r = matrix2::rotation(rotation);
	return r * (scale * v) + position;
}

inline uniform_transform3d::uniform_transform3d(const vector3& position, const quaternion& rotation, float scale) :
	position(position),
	rotation(rotation),
	scale(scale)
{}

inline uniform_transform3d uniform_transform3d::operator*(const uniform_transform3d& t) const
{
	uniform_transform3d res;
	res.position = rotation * (scale * t.position) + position;
	res.rotation = rotation * t.rotation;
	res.scale = scale * t.scale;
	return res;
}
inline vector3 uniform_transform3d::operator*(const vector3& v) const { return rotation * (scale * v) + position; }
inline quaternion uniform_transform3d::operator*(const quaternion& q) const { return rotation * q; }

inline unscaled_transform2d::unscaled_transform2d(const vector2& position, const radian& rotation) :
	position(position),
	rotation(rotation)
{}

inline unscaled_transform2d unscaled_transform2d::operator*(const unscaled_transform2d& t) const
{
	unscaled_transform2d res;
	matrix2 r = matrix2::rotation(rotation);
	res.position = r * t.position + position;
	res.rotation = rotation + t.rotation;
	return res;
}
inline vector2 unscaled_transform2d::operator*(const vector2& v) const
{
	matrix2 r = matrix2::rotation(rotation);
	return r * v + position;
}

inline unscaled_transform3d::unscaled_transform3d(const vector3& position, const quaternion& rotation) :
	position(position),
	rotation(rotation)
{}

inline unscaled_transform3d unscaled_transform3d::operator*(const unscaled_transform3d& t) const
{
	unscaled_transform3d res;
	res.position = rotation * t.position + position;
	res.rotation = rotation * t.rotation;
	return res;
}
inline vector3 unscaled_transform3d::operator*(const vector3& v) const { return rotation * v + position; }
inline quaternion unscaled_transform3d::operator*(const quaternion& q) const { return rotation * q; }

namespace detail
{
template<unsigned char x, unsigned char y>
float get_quaternion_matrix_element(const quaternion& q)
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
} // namespace cdm

#endif // CDM_MATH_HPP