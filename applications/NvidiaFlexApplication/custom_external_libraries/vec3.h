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

#if 0 //_DEBUG
#define VEC3_VALIDATE() {	\
	assert(_finite(x));\
	assert(!_isnan(x));\
	\
	assert(_finite(y));\
	assert(!_isnan(y));\
	\
	assert(_finite(z));\
	assert(!_isnan(z));\
						 }
#else
#define VEC3_VALIDATE()
#endif

template <typename T=float>
class XVector3
{
public:

	typedef T value_type;

	CUDA_CALLABLE inline XVector3() : x(0.0f), y(0.0f), z(0.0f) {}
	CUDA_CALLABLE inline XVector3(T a) : x(a), y(a), z(a) {}
	CUDA_CALLABLE inline XVector3(const T* p) : x(p[0]), y(p[1]), z(p[2]) {}
	CUDA_CALLABLE inline XVector3(T x_, T y_, T z_) : x(x_), y(y_), z(z_) 
	{
		VEC3_VALIDATE();
	}

	CUDA_CALLABLE inline operator T* () { return &x; }
	CUDA_CALLABLE inline operator const T* () const { return &x; };

	CUDA_CALLABLE inline void Set(T x_, T y_, T z_) { VEC3_VALIDATE(); x = x_; y = y_; z = z_;}

	CUDA_CALLABLE inline XVector3<T> operator * (T scale) const { XVector3<T> r(*this); r *= scale; return r; VEC3_VALIDATE();}
	CUDA_CALLABLE inline XVector3<T> operator / (T scale) const { XVector3<T> r(*this); r /= scale; return r; VEC3_VALIDATE();}
	CUDA_CALLABLE inline XVector3<T> operator + (const XVector3<T>& v) const { XVector3<T> r(*this); r += v; return r; VEC3_VALIDATE();}
	CUDA_CALLABLE inline XVector3<T> operator - (const XVector3<T>& v) const { XVector3<T> r(*this); r -= v; return r; VEC3_VALIDATE();}
    CUDA_CALLABLE inline XVector3<T> operator /(const XVector3<T>& v) const { XVector3<T> r(*this); r /= v; return r; VEC3_VALIDATE();}
    CUDA_CALLABLE inline XVector3<T> operator *(const XVector3<T>& v) const { XVector3<T> r(*this); r *= v; return r; VEC3_VALIDATE();}

	CUDA_CALLABLE inline XVector3<T>& operator *=(T scale) {x *= scale; y *= scale; z*= scale; VEC3_VALIDATE(); return *this;}
	CUDA_CALLABLE inline XVector3<T>& operator /=(T scale) {T s(1.0f/scale); x *= s; y *= s; z *= s; VEC3_VALIDATE(); return *this;}
	CUDA_CALLABLE inline XVector3<T>& operator +=(const XVector3<T>& v) {x += v.x; y += v.y; z += v.z; VEC3_VALIDATE(); return *this;}
	CUDA_CALLABLE inline XVector3<T>& operator -=(const XVector3<T>& v) {x -= v.x; y -= v.y; z -= v.z; VEC3_VALIDATE(); return *this;}
    CUDA_CALLABLE inline XVector3<T>& operator /=(const XVector3<T>& v) {x /= v.x; y /= v.y; z /= v.z; VEC3_VALIDATE(); return *this; }
    CUDA_CALLABLE inline XVector3<T>& operator *=(const XVector3<T>& v) {x *= v.x; y *= v.y; z *= v.z; VEC3_VALIDATE(); return *this; }

	CUDA_CALLABLE inline bool operator != (const XVector3<T>& v) const { return (x != v.x || y != v.y || z != v.z); }

	// negate
	CUDA_CALLABLE inline XVector3<T> operator -() const { VEC3_VALIDATE(); return XVector3<T>(-x, -y, -z); }

	CUDA_CALLABLE void Validate()
	{
		VEC3_VALIDATE();
	}

	T x,y,z;
};

typedef XVector3<float> Vec3;
typedef XVector3<float> Vector3;

// lhs scalar scale
template <typename T>
CUDA_CALLABLE XVector3<T> operator *(T lhs, const XVector3<T>& rhs)
{
	XVector3<T> r(rhs);
	r *= lhs;
	return r;
}

template <typename T>
CUDA_CALLABLE bool operator==(const XVector3<T>& lhs, const XVector3<T>& rhs)
{
	return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z);
}

template <typename T>
CUDA_CALLABLE typename T::value_type Dot3(const T& v1, const T& v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z*v2.z; 
}

CUDA_CALLABLE inline float Dot3(const float* v1, const float * v2)
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; 
}


template <typename T>
CUDA_CALLABLE inline T Dot(const XVector3<T>& v1, const XVector3<T>& v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

CUDA_CALLABLE inline Vec3 Cross(const Vec3& b, const Vec3& c)
{
	return Vec3(b.y*c.z - b.z*c.y,
			    b.z*c.x - b.x*c.z,
				b.x*c.y - b.y*c.x);
}

// component wise min max functions
template <typename T>
CUDA_CALLABLE inline XVector3<T> Max(const XVector3<T>& a, const XVector3<T>& b)
{
	return XVector3<T>(Max(a.x, b.x), Max(a.y, b.y), Max(a.z, b.z));
}

template <typename T>
CUDA_CALLABLE inline XVector3<T> Min(const XVector3<T>& a, const XVector3<T>& b)
{
	return XVector3<T>(Min(a.x, b.x), Min(a.y, b.y), Min(a.z, b.z));
}

