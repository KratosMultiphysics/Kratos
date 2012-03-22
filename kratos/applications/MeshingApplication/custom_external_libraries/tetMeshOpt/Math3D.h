#pragma once

#include <iostream>
#include <algorithm>
#include <math.h>
#include <windows.h>
#include "stdafx.h"

double Min ( double a, double b );
double Min ( float a, float b );
double Max ( double a, double b );
double Max ( float a, float b );
double MaxValue(const double *Numbers, const int Count);
double MinValue(const double *Numbers, const int Count);


 struct float4
{
	float x,y,z,w ;
	
	float4 operator-(float4 v1) ;
	   
	float4 operator+(float4 v1) ;
	   
	float4 operator*(float v1) ;
	   
	float4 operator*(float4 v1) ;
	   
	float4 operator=(float v1) ;
	   
   float4 Float4(float v);
   
};

 float4 Float4(float x,float y, float z);
 float4 Float4(double x,double y, double z);
 float4 Float4(float x);

 float4 max4(float4 a,float4 b);

  typedef  struct 
{	
	float4 min;
	float4 max;
	float4 size;
	float4 center;
	
} BoundBox;

  BoundBox calcBound(float4 p,float4 s);
  float dot(float4 v1,float4 v2);
    float4 cross(float4 vVector1,float4 vVector2);
	float magnitude(float4 v1);
	 float4 normalize(float4 vNormal);
	 float4 Normal(float4 v0,float4 v1,float4 v2);
	 float distance(float4 Point1, float4 Point2);
double Area(float4 Point1,float4 Point2,float4 Point3);
double tetraVolume(float4 v0, float4 v1, float4 v2, float4 v3);