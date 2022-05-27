/* cdm_curve v0.0
   C++ splines and curves library
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

#ifndef CDM_CURVE_HPP
#define CDM_CURVE_HPP 1

namespace cdm
{
template <typename point>
struct bezier_curve
{
	point p0, p1, p2, p3;

	point interpolate(float t)
	{
		return bezier_interpolator(p0, p1, p2, p3, t);
	}
};

template <typename point>
inline point bezier_interpolator(point p0,
                                 point p1,
                                 point p2,
                                 point p3,
                                 float t)
{
	const float t2 = t * t;
	const float t3 = t2 * t;

	return (1 - t) * (1 - t) * (1 - t) * p0 + 3 * t * (1 - t) * (1 - t) * p1 +
	       3 * t2 * (1 - t) * p2 + t3 * p3;
}

template <typename point>
inline point bezier_interpolator(bezier_curve<point> b0,
                                 bezier_curve<point> b1,
                                 bezier_curve<point> b2,
                                 bezier_curve<point> b3,
                                 float t,
                                 float u)
{
	bezier_curve<point> b{b0.interpolate(t), b1.interpolate(t),
	                      b2.interpolate(t), b3.interpolate(t)};

	return b.interpolate(u);
}

template <typename point>
struct hermitian_curve
{
	point p0, p1, r0, r1;
};

template <typename point>
inline point hermitian_interpolator(point p0,
                                    point p1,
                                    point r0,
                                    point r1,
                                    float t)
{
	const float t2 = t * t;
	const float t3 = t2 * t;

	return (2 * t3 - 3 * t2 + 1) * p0 + (-2 * t3 + 3 * t2) * p1 +
	       (t3 - 2 * t2 + t) * r0 + (t3 - t2) * r1;
}

template <typename point>
struct bspline_curve
{
	point p0, p1, p2, p3;
};

template <typename point>
inline point bspline_interpolator(point p0,
                                  point p1,
                                  point p2,
                                  point p3,
                                  float t)
{
	const float t2 = t * t;
	const float t3 = t2 * t;

	return ((1 - t) * (1 - t) * (1 - t) * p0 + (3 * t3 - 6 * t2 + 4) * p1 +
	        (-3 * t3 + 3 * t2 + 3 * t + 1) * p2 + t3 * p3) /
	       6.0f;
}
}  // namespace cdm

#endif  // CDM_CURVE_HPP
