

#include "Math3D.h"

float4 float4::operator-(float4 v1) 
	   { float4 r; 
	       r.x = x-v1.x;   
		   r.y = y-v1.y;   
		   r.z = z-v1.z;   
		   r.w = w-v1.w;   
		   return r; }
	float4 float4::operator+(float4 v1) 
	   { float4 r; 
	       r.x = x+v1.x;   
		   r.y = y+v1.y;   
		   r.z = z+v1.z;   
		   r.w = w+v1.w;   
		   return r; }
	float4 float4::operator*(float v1) 
	   { float4 r; 
	       r.x = x*v1;   
		   r.y = y*v1;   
		   r.z = z*v1;   
		   r.w = w*v1;   
		   return r; }
	float4 float4::operator*(float4 v1) 
	   { float4 r; 
	       r.x = x*v1.x;   
		   r.y = y*v1.y;   
		   r.z = z*v1.z;   
		   r.w = w*v1.w;   
		   return r; }
	float4 float4::operator=(float v1) 
	   { float4 r; 
	       r.x = v1;   
		   r.y = v1;   
		   r.z = v1;   
		   r.w = v1;   
		   return r; }
  /* float4 operator=(float4 v1) 
	   { float4 r; 
	       r.x = v1.x;   
		   r.y = v1.y;   
		   r.z = v1.z;   
		   r.w = v1.w;   
		   return r; }*/
   float4 float4::Float4(float v)
   {
	   float4 r; 
	       r.x = v;   
		   r.y = v;   
		   r.z = v;   
		   r.w = v;   
		   return r;
   }

int Min ( int a, int b )
{
  return (a<b)?a:b;     // or: return !comp(b,a)?a:b; for the comp version
}

double Min ( double a, double b )
{
  return (a<b)?a:b;     // or: return !comp(b,a)?a:b; for the comp version
}

float Min ( float a, float b )
{
  return (a<b)?a:b;     // or: return !comp(b,a)?a:b; for the comp version
}

double Max ( double a, double b )
{
  return (a>b)?a:b;     // or: return !comp(b,a)?a:b; for the comp version
}

float Max ( float a, float b )
{
  return (a>b)?a:b;     // or: return !comp(b,a)?a:b; for the comp version
}


template <class T> const T& Min ( const T& a, const T& b )
{
  return !(b<a)?a:b;     // or: return !comp(b,a)?a:b; for the comp version
}

double MaxValue(const double *Numbers, const int Count)
{
	double Maximum = Numbers[0];
	for(int i = 1; i < Count; i++)
		if( Maximum < Numbers[i] )
			Maximum = Numbers[i];
	return Maximum;
}

double MinValue(const double *Numbers, const int Count)
{
	double Minimun = Numbers[0];
	for(int i = 1; i < Count; i++)
		if( Minimun > Numbers[i] )
			Minimun = Numbers[i];
	return Minimun;
}


float4 Float4(float x,float y, float z)
 {
	 float4 f;
	 f.x = x;
	 f.y = y;
	 f.z = z;
	 return f;
 }

 float4 Float4(double x,double y, double z)
 {
	 float4 f;
	 f.x = (float)(x);
	 f.y = (float)(y);
	 f.z = (float)(z);
	 return f;
 }

 float4 Float4(float x)
 {
	 float4 f;
	 f.x = x;
	 f.y = x;
	 f.z = x;
	 return f;
 }

 float4 max4(float4 a,float4 b)
{
	 float4 f;
	 f.x = Max(a.x,b.x);
	 f.y = Max(a.y,b.y);
	 f.z = Max(a.z,b.z);
	 return f;
}

 BoundBox calcBound(float4 p,float4 s)
{
  BoundBox result ;
  result.min.x =Min(p.x,p.x+s.x);
  result.min.y =Min(p.y,p.y+s.y);
  result.min.z =Min(p.z,p.z+s.z);

  result.max.x =Max(p.x,p.x+s.x);
  result.max.y =Max(p.y,p.y+s.y);
  result.max.z =Max(p.z,p.z+s.z);

	result.size.x = result.max.x-result.min.x;
   result.size.y = result.max.y-result.min.y;
   result.size.z = result.max.z-result.min.z;

   result.center = (result.min + result.max)*0.5;

   return result;
}


float dot(float4 v1,float4 v2)
{
  return  v1.x * v2.x + v1.y*v2.y + v1.z * v2.z;
}




  float4 cross(float4 vVector1,float4 vVector2)
  {
    float4 vNormal;

 // The X value for the vector is:  (V1.y * V2.z) - (V1.z * V2.y)													// Get the X value
   vNormal.x= ((vVector1.y * vVector2.z) - (vVector1.z * vVector2.y));
 // The Y value for the vector is:  (V1.z * V2.x) - (V1.x * V2.z)
   vNormal.y= ((vVector1.z * vVector2.x) - (vVector1.x * vVector2.z));
 // The Z value for the vector is:  (V1.x * V2.y) - (V1.y * V2.x)
  vNormal.z= ((vVector1.x * vVector2.y) - (vVector1.y * vVector2.x));
  return vNormal;
  }

  float magnitude(float4 v1)
{
   return sqrt(v1.x*v1.x + v1.y * v1.y + v1.z * v1.z);
}
  float4 normalize(float4 vNormal)
{
  float mmagnitude;

  mmagnitude= magnitude(vNormal);				// Get the magnitude of our normal
  // Now that we have the magnitude, we can divide our normal by that magnitude.
  // That will make our normal a total length of 1.  This makes it easier to work with too.
  vNormal.x=vNormal.x / mmagnitude; // Divide the X value of our normal by it's magnitude
  vNormal.y=vNormal.y / mmagnitude; // Divide the Y value of our normal by it's magnitude
  vNormal.z=vNormal.z / mmagnitude; // Divide the Z value of our normal by it's magnitude
// Finally, return our normalized normal.
  return vNormal; // Return the new normal of length 1.
}

  float4 Normal(float4 v0,float4 v1,float4 v2)
{
  float4 vVector1,vVector2,vNormal;
// Get 2 vectors from the polygon (2 sides), Remember the order!
  vVector1.x= v2.x - v0.x;
  vVector1.y= v2.y - v0.y;
  vVector1.z= v2.z - v0.z;
  vVector2.x= v1.x - v0.x;
  vVector2.y= v1.y - v0.y;
  vVector2.z= v1.z - v0.z;
  vNormal = cross(vVector1, vVector2);
  return normalize(vNormal);
}

  float distance(float4 Point1, float4 Point2)
  {
   float dx,dy,dz;
  
  dx = Point2.x - Point1.x;
  dy = Point2.y - Point1.y;
  dz = Point2.z - Point1.z;
   return  sqrt(dx * dx + dy * dy + dz * dz);
  }


double Area(float4 Point1,float4 Point2,float4 Point3)
{
  float Dx1,  Dx2,   Dy1 ,  Dy2 ,   Dz1 ,  Dz2 ,   Cx  ,  Cy  ,   Cz  ;
  
  Dx1 = Point2.x - Point1.x;
  Dy1 = Point2.y - Point1.y;
  Dz1 = Point2.z - Point1.z;

  Dx2 = Point3.x - Point1.x;
  Dy2 = Point3.y - Point1.y;
  Dz2 = Point3.z - Point1.z;

  Cx  = Dy1 * Dz2 - Dy2 * Dz1;
  Cy  = Dx2 * Dz1 - Dx1 * Dz2;
  Cz  = Dx1 * Dy2 - Dx2 * Dy1;

  return (sqrt(Cx * Cx + Cy * Cy + Cz * Cz) * 0.5);
}

double tetraVolume(float4 v0, float4 v1, float4 v2, float4 v3)
{
    double x21,y21,z21, x31,y31,z31, x41,y41,z41 ,  x123,y123,z123;


    x21  = v1.x - v0.x;
    y21  = v1.y - v0.y;
    z21  = v1.z - v0.z;

    x31  = v2.x - v0.x;
    y31  = v2.y - v0.y;
    z31  = v2.z - v0.z;

    x41  = v3.x - v0.x;
    y41  = v3.y - v0.y;
    z41  = v3.z - v0.z;

    x123 = y21*z31 - z21*y31;
    y123 = z21*x31 - x21*z31;
    z123 = x21*y31 - y21*x31;

    return (x123*x41 + y123*y41 + z123*z41) / 6;
}


