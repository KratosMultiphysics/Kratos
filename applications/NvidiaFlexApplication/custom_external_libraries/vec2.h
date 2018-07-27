// This code contains NVIDIA Confidential Information and is disclosed to you
// under a form of NVIDIA software license agreement provided separately to you.
//
// Notice
// NVIDIA Corporation and its licensors retain all intellectual property and
// proprietary rights in and to this software and related documentation and
// any modifications thereto. Any use, reproduction, disclosure, or
// distribution of this software and related documentation without an express
// license agreement from NVIDIA Corporation is strictly prohibited.
//
// ALL NVIDIA DESIGN SPECIFICATIONS, CODE ARE PROVIDED "AS IS.". NVIDIA MAKES
// NO WARRANTIES, EXPRESSED, IMPLIED, STATUTORY, OR OTHERWISE WITH RESPECT TO
// THE MATERIALS, AND EXPRESSLY DISCLAIMS ALL IMPLIED WARRANTIES OF NONINFRINGEMENT,
// MERCHANTABILITY, AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Information and code furnished is believed to be accurate and reliable.
// However, NVIDIA Corporation assumes no responsibility for the consequences of use of such
// information or for any infringement of patents or other rights of third parties that may
// result from its use. No license is granted by implication or otherwise under any patent
// or patent rights of NVIDIA Corporation. Details are subject to change without notice.
// This code supersedes and replaces all information previously supplied.
// NVIDIA Corporation products are not authorized for use as critical
// components in life support devices or systems without express written approval of
// NVIDIA Corporation.
//
// Copyright (c) 2013-2016 NVIDIA Corporation. All rights reserved.

#pragma once

#if defined(_WIN32) && !defined(__CUDACC__)
#if defined(_DEBUG)

#define VEC2_VALIDATE() {	assert(_finite(x));\
	assert(!_isnan(x));\
	\
	assert(_finite(y));\
	assert(!_isnan(y));\
						 }
#else

#define VEC2_VALIDATE() {\
	assert(isfinite(x));\
	assert(isfinite(y)); }\

#endif // _WIN32

#else
#define VEC2_VALIDATE()
#endif

#ifdef _DEBUG
#define FLOAT_VALIDATE(f) { assert(_finite(f)); assert(!_isnan(f)); }
#else
#define FLOAT_VALIDATE(f)
#endif


// vec2
template <typename T>
class XVector2
{
public:

	typedef T value_type;

	CUDA_CALLABLE XVector2() : x(0.0f), y(0.0f) { VEC2_VALIDATE(); }
	CUDA_CALLABLE XVector2(T _x) : x(_x), y(_x) { VEC2_VALIDATE(); }
	CUDA_CALLABLE XVector2(T _x, T _y) : x(_x), y(_y) { VEC2_VALIDATE(); }
	CUDA_CALLABLE XVector2(const T* p) : x(p[0]), y(p[1]) {}

	template <typename U>
	CUDA_CALLABLE explicit XVector2(const XVector2<U>& v) : x(v.x), y(v.y) {}

	CUDA_CALLABLE operator T* () { return &x; }
	CUDA_CALLABLE operator const T* () const { return &x; };

	CUDA_CALLABLE void Set(T x_, T y_) { VEC2_VALIDATE(); x = x_; y = y_; }

	CUDA_CALLABLE XVector2<T> operator * (T scale) const { XVector2<T> r(*this); r *= scale; VEC2_VALIDATE(); return r; }
	CUDA_CALLABLE XVector2<T> operator / (T scale) const { XVector2<T> r(*this); r /= scale; VEC2_VALIDATE(); return r; }
	CUDA_CALLABLE XVector2<T> operator + (const XVector2<T>& v) const { XVector2<T> r(*this); r += v; VEC2_VALIDATE(); return r; }
	CUDA_CALLABLE XVector2<T> operator - (const XVector2<T>& v) const { XVector2<T> r(*this); r -= v; VEC2_VALIDATE(); return r; }

	CUDA_CALLABLE XVector2<T>& operator *=(T scale) {x *= scale; y *= scale; VEC2_VALIDATE(); return *this;}
	CUDA_CALLABLE XVector2<T>& operator /=(T scale) {T s(1.0f/scale); x *= s; y *= s; VEC2_VALIDATE(); return *this;}
	CUDA_CALLABLE XVector2<T>& operator +=(const XVector2<T>& v) {x += v.x; y += v.y; VEC2_VALIDATE(); return *this;}
	CUDA_CALLABLE XVector2<T>& operator -=(const XVector2<T>& v) {x -= v.x; y -= v.y; VEC2_VALIDATE(); return *this;}

	CUDA_CALLABLE XVector2<T>& operator *=(const XVector2<T>& scale) {x *= scale.x; y *= scale.y; VEC2_VALIDATE(); return *this;}

	// negate
	CUDA_CALLABLE XVector2<T> operator -() const { VEC2_VALIDATE(); return XVector2<T>(-x, -y); }

	// returns this vector
	CUDA_CALLABLE void Normalize() { *this /= Length(*this); }
	CUDA_CALLABLE void SafeNormalize(const XVector2<T>& v=XVector2<T>(0.0f,0.0f))
	{
		T length = Length(*this);
		*this = (length==0.00001f)?v:(*this /= length);
	}

	T x;
	T y;
};

typedef XVector2<float> Vec2;
typedef XVector2<float> Vector2;

// lhs scalar scale
template <typename T>
CUDA_CALLABLE XVector2<T> operator *(T lhs, const XVector2<T>& rhs)
{
	XVector2<T> r(rhs);
	r *= lhs;
	return r;
}

template <typename T>
CUDA_CALLABLE XVector2<T> operator*(const XVector2<T>& lhs, const XVector2<T>& rhs)
{
	XVector2<T> r(lhs);
	r *= rhs;
	return r;
}

template <typename T>
CUDA_CALLABLE bool operator==(const XVector2<T>& lhs, const XVector2<T>& rhs)
{
	return (lhs.x == rhs.x && lhs.y == rhs.y);
}


template <typename T>
CUDA_CALLABLE T Dot(const XVector2<T>& v1, const XVector2<T>& v2)
{
	return v1.x * v2.x + v1.y * v2.y; 
}

// returns the ccw perpendicular vector 
template <typename T>
CUDA_CALLABLE XVector2<T> PerpCCW(const XVector2<T>& v)
{
	return XVector2<T>(-v.y, v.x);
}

template <typename T>
CUDA_CALLABLE XVector2<T> PerpCW(const XVector2<T>& v)
{
	return XVector2<T>(v.y, -v.x);
}

// component wise min max functions
template <typename T>
CUDA_CALLABLE XVector2<T> Max(const XVector2<T>& a, const XVector2<T>& b)
{
	return XVector2<T>(Max(a.x, b.x), Max(a.y, b.y));
}

template <typename T>
CUDA_CALLABLE XVector2<T> Min(const XVector2<T>& a, const XVector2<T>& b)
{
	return XVector2<T>(Min(a.x, b.x), Min(a.y, b.y));
}

// 2d cross product, treat as if a and b are in the xy plane and return magnitude of z
template <typename T>
CUDA_CALLABLE T Cross(const XVector2<T>& a, const XVector2<T>& b)
{
	return (a.x*b.y - a.y*b.x);
}




