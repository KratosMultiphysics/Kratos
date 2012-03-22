
#include <math.h>
#include "Math3D.h"
#include "u_Types.h"
#include "u_delphiClasses.h"
#include "u_TetraFunctions.h"

#define PI 3.1416
#define PRECISSION 1024*1024*8

void CalcAng(float4 coor0,float4 coor1,float4 coor2,float4 coor3 , double* dangle, double *cangle, double *vol)
{
	double   x21,y21,z21, x31,y31,z31, x41,y41,z41, x32,y32,z32, x42,y42,z42, x43,y43,z43  ,
		x123,y123,z123,o123, x124,y124,z124,o124, x134,y134,z134,o134, x234,y234,z234,o234,
		o21,o31,o41,o32,o42,o43;
	double d[6];
	double c[12];
	x21  = coor1.x - coor0.x;
	y21  = coor1.y - coor0.y;
	z21  = coor1.z - coor0.z;

	x31  = coor2.x - coor0.x;
	y31  = coor2.y - coor0.y;
	z31  = coor2.z - coor0.z;

	x41  = coor3.x - coor0.x;
	y41  = coor3.y - coor0.y;
	z41  = coor3.z - coor0.z;

	x32  = coor2.x - coor1.x;
	y32  = coor2.y - coor1.y;
	z32  = coor2.z - coor1.z;

	x42  = coor3.x - coor1.x;
	y42  = coor3.y - coor1.y;
	z42  = coor3.z - coor1.z;

	x43  = coor3.x - coor2.x;
	y43  = coor3.y - coor2.y;
	z43  = coor3.z - coor2.z;

	x123 = y21*z31 - z21*y31;
	y123 = z21*x31 - x21*z31;
	z123 = x21*y31 - y21*x31;
	o123 = sqrt(x123*x123 + y123*y123 + z123*z123);

	x124 = y41*z21 - z41*y21; 
	y124 = z41*x21 - x41*z21;
	z124 = x41*y21 - y41*x21;
	o124 = sqrt(x124*x124 + y124*y124 + z124*z124);

	x134 = y31*z41 - z31*y41;
	y134 = z31*x41 - x31*z41;
	z134 = x31*y41 - y31*x41;
	o134 = sqrt(x134*x134 + y134*y134 + z134*z134);

	x234 = y42*z32 - z42*y32;
	y234 = z42*x32 - x42*z32;
	z234 = x42*y32 - y42*x32;
	o234 = sqrt(x234*x234 + y234*y234 + z234*z234);

	d[0] = (x123*x124 + y123*y124 + z123*z124)/(o123*o124);
	d[1] = (x123*x134 + y123*y134 + z123*z134)/(o123*o134);
	d[2] = (x134*x124 + y134*y124 + z134*z124)/(o134*o124);
	d[3] = (x123*x234 + y123*y234 + z123*z234)/(o123*o234);
	d[4] = (x124*x234 + y124*y234 + z124*z234)/(o124*o234);
	d[5] = (x134*x234 + y134*y234 + z134*z234)/(o134*o234);
	if (d[0] < -1.0)  d[0] = -1.0;
	if (d[0] >  1.0)  d[0] =  1.0;
	if (d[1] < -1.0)  d[1] = -1.0;
	if (d[1] >  1.0)  d[1] =  1.0;
	if (d[2] < -1.0)  d[2] = -1.0;
	if (d[2] >  1.0)  d[2] =  1.0;
	if (d[3] < -1.0)  d[3] = -1.0;
	if (d[3] >  1.0)  d[3] =  1.0;
	if (d[4] < -1.0)  d[4] = -1.0;
	if (d[4] >  1.0)  d[4] =  1.0;
	if (d[5] < -1.0)  d[5] = -1.0;
	if (d[5] >  1.0)  d[5] =  1.0;

	d[0] = (PI-acos(d[0])) * 180 / PI;
	d[1] = (PI-acos(d[1])) * 180 / PI;
	d[2] = (PI-acos(d[2])) * 180 / PI;
	d[3] = (PI-acos(d[3])) * 180 / PI;
	d[4] = (PI-acos(d[4])) * 180 / PI;
	d[5] = (PI-acos(d[5])) * 180 / PI;

	o21 = sqrt(x21*x21 + y21*y21 + z21*z21);
	o31 = sqrt(x31*x31 + y31*y31 + z31*z31);
	o41 = sqrt(x41*x41 + y41*y41 + z41*z41);
	o32 = sqrt(x32*x32 + y32*y32 + z32*z32);
	o42 = sqrt(x42*x42 + y42*y42 + z42*z42);
	o43 = sqrt(x43*x43 + y43*y43 + z43*z43);
	c[0] = (x21*x31+y21*y31+z21*z31) / (o21*o31);
	c[1] =-(x21*x32+y21*y32+z21*z32) / (o21*o32);
	c[3] = (x21*x41+y21*y41+z21*z41) / (o21*o41);
	c[4] =-(x21*x42+y21*y42+z21*z42) / (o21*o42);
	c[6] = (x32*x42+y32*y42+z32*z42) / (o32*o42);
	c[7] =-(x32*x43+y32*y43+z32*z43) / (o32*o43);
	c[9] = (x31*x41+y31*y41+z31*z41) / (o31*o41);
	c[10]=-(x31*x43+y31*y43+z31*z43) / (o31*o43);
	if (c[0] < -1.0)  c[0] = -1.0;
	if (c[0] >  1.0)  c[0] =  1.0;
	if (c[1] < -1.0)  c[1] = -1.0;
	if (c[1] >  1.0)  c[1] =  1.0;
	if (c[3] < -1.0)  c[3] = -1.0;
	if (c[3] >  1.0)  c[3] =  1.0;
	if (c[4] < -1.0)  c[4] = -1.0;
	if (c[4] >  1.0)  c[4] =  1.0;
	if (c[6] < -1.0)  c[6] = -1.0;
	if (c[6] >  1.0)  c[6] =  1.0;
	if (c[7] < -1.0)  c[7] = -1.0;
	if (c[7] >  1.0)  c[7] =  1.0;
	if (c[9] < -1.0)  c[9] = -1.0;
	if (c[9] >  1.0)  c[9] =  1.0;
	if (c[10]< -1.0)  c[10]= -1.0;
	if (c[10]>  1.0)  c[10]=  1.0;

	c[0] = acos(c[0]);
	c[1] = acos(c[1]);
	c[2] = PI - c[0] - c[1];
	c[3] = acos(c[3]);
	c[4] = acos(c[4]);
	c[5] = PI - c[3] - c[4];
	c[6] = acos(c[6]);
	c[7] = acos(c[7]);
	c[8] = PI - c[6] - c[7];
	c[9] = acos(c[9]);
	c[10]= acos(c[10]);
	c[11]= PI - c[9] - c[10];

	c[0]  = c[0] * 180 / PI;
	c[1]  = c[1] * 180 / PI;
	c[2]  = c[2] * 180 / PI;
	c[3]  = c[3] * 180 / PI;
	c[4]  = c[4] * 180 / PI;
	c[5]  = c[5] * 180 / PI;
	c[6]  = c[6] * 180 / PI;
	c[7]  = c[7] * 180 / PI;
	c[8]  = c[8] * 180 / PI;
	c[9]  = c[9] * 180 / PI;
	c[10] = c[10]* 180 / PI;
	c[11] = c[11]* 180 / PI;

	*dangle = MinValue(d,6);
	*cangle = MinValue(c,12);
	*vol = (x123*x41 + y123*y41 + z123*z41) / 6;

};



double calidadxArea(float4 v1,float4 v2,float4 v3)
{
	float   surf, sumLados;

	sumLados = magnitude(v2-v1)+magnitude(v3-v2)+magnitude(v3-v1);
	surf = (float)Area(v1,v2,v3);
	if (sumLados==0)
		return 0;
	else
		return  (surf ) /(sumLados*sumLados);
}

double calidadxArea(TVertex* v1,TVertex* v2,TVertex* v3)
{
	float   surf, sumLados;

	sumLados = magnitude(v2->fPos-v1->fPos)+magnitude(v3->fPos-v2->fPos)+magnitude(v3->fPos-v1->fPos);
	surf = (float)Area(v1->fPos,v2->fPos,v3->fPos);
	if (sumLados==0)
		return 0;
	else
		return  (surf ) /(sumLados*sumLados);
}

double getmaxEdgeLength(TVertex* v0,TVertex* v1,TVertex* v2,TVertex* v3)
{
	double Nbrs[] = { distance(  v0->pos(),v1->pos()  ),
		distance(  v0->pos(),v2->pos()  ),
		distance(  v0->pos(),v3->pos()  ),
		distance(  v1->pos(),v2->pos()  ),
		distance(  v1->pos(),v3->pos()  ),
		distance(  v2->pos(),v3->pos()  ) };

	return  MaxValue(Nbrs,6 );
}

double vrelaxQuality(TVertex* v0, TVertex* v1,TVertex* v2,TVertex* v3)
{
	double qlty ,
		// Calculo el volumen
		x21,y21,z21,x31,y31,z31,x41,y41,z41 ,
		d1,d2,d3,vol ,
		// Calculo la calidad
		dx,dy,dz ,
		maxleng2, leng2 ;
	float4 p1,p2,p3, p4 ;

	p1 = v0->fPos;
	p2 = v1->fPos;
	p3 = v2->fPos;
	p4 = v3->fPos;

	x21 = p2.x-p1.x; y21 = p2.y-p1.y; z21 = p2.z-p1.z;
	x31 = p3.x-p1.x; y31 = p3.y-p1.y; z31 = p3.z-p1.z;
	x41 = p4.x-p1.x; y41 = p4.y-p1.y; z41 = p4.z-p1.z;
	d1 = y31*z41 - z31*y41;
	d2 = y41*z21 - z41*y21;
	d3 = y21*z31 - z21*y31;
	vol = (x21 * d1 + x31 * d2 + x41 * d3) / 6.0;
	qlty = vol;

	dx = p2.x-p1.x; 
	dy = p2.y-p1.y; 
	dz = p2.z-p1.z;
	maxleng2=dx*dx+dy*dy+dz*dz;
	dx = p3.x-p1.x; dy = p3.y-p1.y; dz = p3.z-p1.z;
	leng2=dx*dx+dy*dy+dz*dz;
	if(leng2 > maxleng2) maxleng2=leng2;
	dx = p4.x-p1.x; dy = p4.y-p1.y; dz = p4.z-p1.z;
	leng2=dx*dx+dy*dy+dz*dz;
	if(leng2 > maxleng2) maxleng2=leng2;
	dx = p3.x-p2.x; dy = p3.y-p2.y; dz = p3.z-p2.z;
	leng2=dx*dx+dy*dy+dz*dz;
	if(leng2 > maxleng2) maxleng2=leng2;
	dx = p4.x-p2.x; dy = p4.y-p2.y; dz = p4.z-p2.z;
	leng2=dx*dx+dy*dy+dz*dz;
	if(leng2 > maxleng2) maxleng2=leng2;
	dx = p4.x-p3.x; dy = p4.y-p3.y; dz = p4.z-p3.z;
	leng2=dx*dx+dy*dy+dz*dz;
	if(leng2 > maxleng2) maxleng2=leng2;

	if (vol>0) 
		qlty =(vol*(vol))/(maxleng2*maxleng2*maxleng2);
	else
		qlty =(vol*(-vol))/(maxleng2*maxleng2*maxleng2);


	int iqlty =dround( qlty * PRECISSION);    

	return 72.0 * ( (double)(iqlty)/(double)(PRECISSION));  //calcula fast quality
}


double vrelaxQuality(TVertex* vertexes[])
{
	return vrelaxQuality(vertexes[0], vertexes[1], vertexes[2], vertexes[3]);
}

float relaxQuality(TObject* o)
{
	if ( (TTetra*)o != NULL)
	{
		TTetra* t = (TTetra*)(o);
		return (float)vrelaxQuality( t->vertexes );
	}
	else
		return 0.0f;
}

double diedralAngle(float4 v0, float4 v1, float4 v2, float4 v3)
{
	double d[6];
	double x21,y21,z21, x31,y31,z31, x41,y41,z41, x32,y32,z32, x42,y42,z42, x43,y43,z43 ,
		x123,y123,z123,o123, x124,y124,z124,o124, x134,y134,z134,o134, x234,y234,z234,o234    ;

	x21  = v1.x - v0.x;
	y21  = v1.y - v0.y;
	z21  = v1.z - v0.z;

	x31  = v2.x - v0.x;
	y31  = v2.y - v0.y;
	z31  = v2.z - v0.z;

	x41  = v3.x - v0.x;
	y41  = v3.y - v0.y;
	z41  = v3.z - v0.z;

	x32  = v2.x - v1.x;
	y32  = v2.y - v1.y;
	z32  = v2.z - v1.z;

	x42  = v3.x - v1.x;
	y42  = v3.y - v1.y;
	z42  = v3.z - v1.z;

	x43  = v3.x - v2.x;
	y43  = v3.y - v2.y;
	z43  = v3.z - v2.z;

	x123 = y21*z31 - z21*y31;
	y123 = z21*x31 - x21*z31;
	z123 = x21*y31 - y21*x31;
	o123 = sqrt(x123*x123 + y123*y123 + z123*z123);

	x124 = y41*z21 - z41*y21;
	y124 = z41*x21 - x41*z21;
	z124 = x41*y21 - y41*x21;
	o124 = sqrt(x124*x124 + y124*y124 + z124*z124);

	x134 = y31*z41 - z31*y41;
	y134 = z31*x41 - x31*z41;
	z134 = x31*y41 - y31*x41;
	o134 = sqrt(x134*x134 + y134*y134 + z134*z134);

	x234 = y42*z32 - z42*y32;
	y234 = z42*x32 - x42*z32;
	z234 = x42*y32 - y42*x32;
	o234 = sqrt(x234*x234 + y234*y234 + z234*z234);

	d[0] = (x123*x124 + y123*y124 + z123*z124)/(o123*o124);
	d[1] = (x123*x134 + y123*y134 + z123*z134)/(o123*o134);
	d[2] = (x134*x124 + y134*y124 + z134*z124)/(o134*o124);
	d[3] = (x123*x234 + y123*y234 + z123*z234)/(o123*o234);
	d[4] = (x124*x234 + y124*y234 + z124*z234)/(o124*o234);
	d[5] = (x134*x234 + y134*y234 + z134*z234)/(o134*o234);
	if (d[0] < -1.0)  d[0] = -1.0;
	if (d[0] >  1.0)  d[0] =  1.0;
	if (d[1] < -1.0)  d[1] = -1.0;
	if (d[1] >  1.0)  d[1] =  1.0;
	if (d[2] < -1.0)  d[2] = -1.0;
	if (d[2] >  1.0)  d[2] =  1.0;
	if (d[3] < -1.0)  d[3] = -1.0;
	if (d[3] >  1.0)  d[3] =  1.0;
	if (d[4] < -1.0)  d[4] = -1.0;
	if (d[4] >  1.0)  d[4] =  1.0;
	if (d[5] < -1.0)  d[5] = -1.0;
	if (d[5] >  1.0)  d[5] =  1.0;

	d[0] = (PI-acos(d[0])) * 180 / PI;
	d[1] = (PI-acos(d[1])) * 180 / PI;
	d[2] = (PI-acos(d[2])) * 180 / PI;
	d[3] = (PI-acos(d[3])) * 180 / PI;
	d[4] = (PI-acos(d[4])) * 180 / PI;
	d[5] = (PI-acos(d[5])) * 180 / PI;

	return MinValue( d,6);
}
