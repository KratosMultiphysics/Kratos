#pragma OPENCL EXTENSION cl_amd_fp64: enable
#pragma OPENCL EXTENSION cl_amd_printf: enable

#ifndef cl_amd_fp64

	// Failed, probably we are not on an ATI platform, so try Khronos version
	#pragma OPENCL EXTENSION cl_khr_fp64: enable

#endif

#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics: enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics: enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_extended_atomics : enable
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics : enable

/* Geometry tools for BinsObject */

inline double4 cross3(double4 a, double4 b)
{
    double4 auxa, auxb;

    auxa = (double4)(a.x,a.y,a.z,0);
    auxb = (double4)(b.x,b.y,b.z,0);

    return cross(auxa,auxb);
}

inline double dot3(double4 a, double4 b)
{
    double4 auxa, auxb;

    auxa = (double4)(a.x,a.y,a.z,0);
    auxb = (double4)(b.x,b.y,b.z,0);

    return dot(auxa,auxb);
}

inline bool planeBoxOverlap(double4 normal, double d, double4 maxbox)
{
    int q;

    double4  vmin,vmax;

    if(normal.x>0.0f)
    {
	vmin.x=-maxbox.x;
	vmax.x= maxbox.x;
    }
    else
    {
	vmin.x=maxbox.x;
	vmax.x=-maxbox.x;
    }

    if(normal.y>0.0f)
    {
	vmin.y=-maxbox.y;
	vmax.y= maxbox.y;
    }
    else
    {
	vmin.y=maxbox.y;
	vmax.y=-maxbox.y;
    }


    if(normal.z>0.0f)
    {
	vmin.z=-maxbox.z;
	vmax.z= maxbox.z;
    }
    else
    {
	vmin.z=maxbox.z;
	vmax.z=-maxbox.z;
    }

    if(dot3(normal,vmin)+d>0.0f) 
	return false;
    if(dot3(normal,vmax)+d>=0.0f) 
	return true;
    
    return false;
}

inline void FindMinMax(double x0, 
		       double x1,
		       double x2,
		       double *min, //Estas han de ir asi.
		       double *max) 
{
    (*min) = (*max) = x0;   

    if(x1<(*min)) 
	(*min)=x1;

    if(x1>(*max)) 
	(*max)=x1;

    if(x2<(*min)) 
	(*min)=x2;

    if(x2>(*max)) 
	(*max)=x2;
}
 

/*======================== X-tests ========================*/
inline unsigned int  AxisTest_X01(double a,   double b, 
				  double fa,  double fb,
				  double p0,  double p2,
				  double min, double max, double rad, 
				  double4 v0,
				  double4 v2,
				  double4 boxhalfsize
				  )
{
    p0 = a*v0.y - b*v0.z;			       	   
    p2 = a*v2.y - b*v2.z;	
		       	   
    if(p0<p2) 
    {
	min=p0; 
	max=p2;
    } 
    else 
    {
	min=p2; 
	max=p0;
    } 

    rad = fa * boxhalfsize.y + fb * boxhalfsize.z;   

    if(min>rad || max<-rad) 
	return 0;
    else 
	return 1;
}

inline unsigned int  AxisTest_X2(double a,   double b, 
				 double fa,  double fb,
				 double p0,  double p1,
				 double min, double max, double rad, 
				 double4 v0,
				 double4 v1,
				 double4 boxhalfsize
				 )
{
    p0 = a*v0.y - b*v0.z;			           
    p1 = a*v1.y - b*v1.z;
			       	   
    if(p0<p1) 
    {
	min=p0; 
	max=p1;
    } 
    else 
    {
	min=p1; 
	max=p0;
    } 

    rad = fa * boxhalfsize.y + fb * boxhalfsize.z;   

    if(min>rad || max<-rad) 
	return 0;
    else 
	return 1;
      }

/*======================== Y-tests ========================*/
inline unsigned int  AxisTest_Y02(double a,   double b, 
				  double fa,  double fb,
				  double p0,  double p2,
				  double min, double max, double rad, 
				  double4 v0,
				  double4 v2,
				  double4 boxhalfsize
				  )
{
    p0 = -a*v0.x + b*v0.z;		      	   
    p2 = -a*v2.x + b*v2.z;	
       	       	   
    if(p0<p2) 
    {
	min=p0;
	max=p2;
    }
    else
    {
	min=p2;
	max=p0;
    } 

    rad = fa * boxhalfsize.x + fb * boxhalfsize.z;   

    if(min>rad || max<-rad) 
	return 0;
    else 
	return 1;
}
        
inline unsigned int  AxisTest_Y1(double a,   double b, 
				 double fa,  double fb,
				 double p0,  double p1,
				 double min, double max, double rad, 
				 double4 v0,
				 double4 v1,
				 double4 boxhalfsize
				 )			   
{
    p0 = -a*v0.x + b*v0.z;
    p1 = -a*v1.x + b*v1.z;

    if(p0<p1) 
    {
	min=p0; 
	max=p1;
    }
    else
    {
	min=p1; 
	max=p0;
    } 

    rad = fa * boxhalfsize.x + fb * boxhalfsize.z;  

    if(min>rad || max<-rad) 
	return 0;
    else 
	return 1;
}
	 
/*======================== Z-tests ========================*/

inline unsigned int  AxisTest_Z12(double a,   double b, 
				  double fa,  double fb,
				  double p1,  double p2,
				  double min, double max, double rad, 
				  double4 v1,
				  double4 v2,
				  double4 boxhalfsize
				  )
{
    p1 = a*v1.x - b*v1.y;			           
    p2 = a*v2.x - b*v2.y;	
			    
    if(p2<p1) 
    {
	min=p2; 
	max=p1;
    } 
    else
    {
	min=p1; 
	max=p2;
    } 

    rad = fa * boxhalfsize.x + fb * boxhalfsize.y; 

    if(min>rad || max<-rad) 
	return 0;
    else 
	return 1;
}  
			 
inline unsigned int  AxisTest_Z0(double a,   double b, 
				 double fa,  double fb,
				 double p0,  double p1,
				 double min, double max, double rad, 
				 double4 v0,
				 double4 v1,
				 double4 boxhalfsize
				 ) 
{
    p0 = a*v0.x - b*v0.y;
    p1 = a*v1.x - b*v1.y;	
			    
    if(p0<p1) 
    {
	min=p0; 
	max=p1;
    } 
    else 
    {
	min=p1; 
	max=p0;
    } 

    rad = fa * boxhalfsize.x + fb * boxhalfsize.y;   
    if(min>rad || max<-rad) 
	return 0;
    else 
	return 1;
} 


inline bool TriBoxOverlap(double4 boxcenter, 
			  double4 boxhalfsize, 
			  double4 * triverts)
{

    /*    use separating axis theorem to test overlap between triangle and box */
    /*    need to test for overlap in these directions: */
    /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
    /*       we do not even need to test these) */
    /*    2) normal of the triangle */
    /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
    /*       this gives 3x3=9 more tests */
    
    double min,max,d,p0,p1,p2,rad,fex,fey,fez;  
    double4 v0,v1,v2;
    double4 axis;
    double4 normal, e0, e1 ,e2;
// 
// 		  /* This is the fastest branch on Sun */
// 		  /* move everything so that the boxcenter is in (0,0,0) */
    v0 = triverts[0]- boxcenter;
    v1 = triverts[1]- boxcenter;
    v2 = triverts[2]- boxcenter;
// 
// 		  /* compute triangle edges */
    e0 = v1 - v0;      /* tri edge 0 */
    e1 = v2 - v1;      /* tri edge 1 */
    e2 = v0 - v2;      /* tri edge 2 */
// 
// 		  /* Bullet 3:  */
// 		  /*  test the 9 tests first (this was faster) */
    fex = fabs(e0.x);
    fey = fabs(e0.y);
    fez = fabs(e0.z);		  
    //AXISTEST_X01(e0[2], e0[1], fez, fey);
    if (AxisTest_X01(e0.z, e0.y, fez, fey, p0, p2, min,max, rad, v0, v2, boxhalfsize)==0) 
	return false; 
    //AXISTEST_Y02(e0[2], e0[0], fez, fex);
    if(AxisTest_Y02( e0.z, e0.x, fez, fex, p0, p2, min,max, rad, v0, v2, boxhalfsize)==0) 
	return false;
    //AXISTEST_Z12(e0[1], e0[0], fey, fex);
    if(AxisTest_Z12(e0.y, e0.x, fey, fex, p1, p2, min, max, rad, v1,v2, boxhalfsize )==0) 
	return false;
	      
    
    
    fex = fabs(e1.x);
    fey = fabs(e1.y);
    fez = fabs(e1.z);
    //AXISTEST_X01(e1[2], e1[1], fez, fey);
    if( AxisTest_X01(e1.z, e1.y, fez, fey, p0, p2, min, max, rad, v0, v2, boxhalfsize)==0) 
	return false; 
    //AXISTEST_Y02(e1[2], e1[0], fez, fex);
    if( AxisTest_Y02(e1.z, e1.x, fez, fex, p0, p2, min, max, rad, v0, v2, boxhalfsize)==0) 
	return false;
    //AXISTEST_Z0(e1[1], e1[0], fey, fex);
    if( AxisTest_Z0(e1.y, e1.x, fey, fex, p0,  p1, min, max, rad, v0, v1, boxhalfsize)==0) 
	return false; 
    

    fex = fabs(e2.x);
    fey = fabs(e2.y);
    fez = fabs(e2.z);
    //AXISTEST_X2(e2[2], e2[1], fez, fey);
    if (AxisTest_X2(e2.z, e2.y, fez, fey, p0, p1, min, max, rad, v0, v1, boxhalfsize )==0) 
	return false;
    //AXISTEST_Y1(e2[2], e2[0], fez, fex);
    if (AxisTest_Y1(e2.z, e2.x, fez, fex, p0, p1, min, max, rad, v0, v1, boxhalfsize )==0) 
	return false; 
    //AXISTEST_Z12(e2[1], e2[0], fey, fex);
    if(AxisTest_Z12(e2.y, e2.x, fey, fex, p1, p2, min, max, rad, v1, v2, boxhalfsize ) ==0) 
	return false;
    

    /* Bullet 1: */
    /*  first test overlap in the {x,y,z}-directions */
    /*  find min, max of the triangle each direction, and test for overlap in */
    /*  that direction -- this is equivalent to testing a minimal AABB around */
    /*  the triangle against the AABB */

    /* test in X-direction */
    FindMinMax(v0.x,v1.x,v2.x,&min,&max);
    if(min>boxhalfsize.x || max<-boxhalfsize.x) 
	return false;

    /* test in Y-direction */
    FindMinMax(v0.y,v1.y,v2.y,&min,&max);
    if(min>boxhalfsize.y || max<-boxhalfsize.y) 
	return false;

    /* test in Z-direction */
    FindMinMax(v0.z,v1.z,v2.z,&min,&max);
    if(min>boxhalfsize.z || max<-boxhalfsize.z) 
	return false;

    /* Bullet 2: */
    /*  test if the box intersects the plane of the triangle */
    /*  compute plane equation of triangle: normal*x+d=0 */
//     MathUtils<double>::CrossProduct(normal, e0, e1);

    normal = cross3(e0,e1); //Vigilar aui con la .w

//     d=-inner_prod(normal,v0);  /* plane eq: normal.x+d=0 */ Se puede usar dot?

    d = -dot3(normal,v0); //Vigilar aqui tambien con la .w

    if(!planeBoxOverlap(normal,d,boxhalfsize)) 
	return false;

    return true;   /* box and triangle overlaps */
}

bool HasIntersection(__global double4 * Points, __global int4 * Triangles, double4 rLowPoint, double4 rHighPoint) 
{
    int t = get_global_id(0);
    int dim = 3;

    double4 boxcenter;
    double4 boxhalfsize;

    boxcenter.x = 0.50 * (rLowPoint.x + rHighPoint.x); 
    boxcenter.y = 0.50 * (rLowPoint.y + rHighPoint.y);

    if(dim == 3) 
      boxcenter.z = 0.50 * (rLowPoint.z + rHighPoint.z);
    else 
      boxcenter.z = 0;

    boxhalfsize.x = 0.50 * (rHighPoint.x - rLowPoint.x); 
    boxhalfsize.y = 0.50 * (rHighPoint.y - rLowPoint.y);
    if(dim == 3) 
      boxhalfsize.z = 0.50 * (rHighPoint.z - rLowPoint.z);
    else 
      boxhalfsize.z = 0;

    double4 triverts[3];

    triverts[0] = Points[Triangles[t].x-1];
    triverts[1] = Points[Triangles[t].y-1];
    triverts[2] = Points[Triangles[t].z-1];

    return TriBoxOverlap(boxcenter, boxhalfsize, triverts);
}

bool HasIntersection3D(__global double4 * Points, __global int4 * Triangles, double4 rLowPoint, double4 rHighPoint) 
{
    int t = get_global_id(0);
    int dim = 3;

    double4 boxcenter;
    double4 boxhalfsize;

    boxcenter.x = 0.50 * (rLowPoint.x + rHighPoint.x); 
    boxcenter.y = 0.50 * (rLowPoint.y + rHighPoint.y);

    if(dim == 3) 
      boxcenter.z = 0.50 * (rLowPoint.z + rHighPoint.z);
    else 
      boxcenter.z = 0;

    boxhalfsize.x = 0.50 * (rHighPoint.x - rLowPoint.x); 
    boxhalfsize.y = 0.50 * (rHighPoint.y - rLowPoint.y);
    if(dim == 3) 
      boxhalfsize.z = 0.50 * (rHighPoint.z - rLowPoint.z);
    else 
      boxhalfsize.z = 0;

//??¿¿??¿¿??¿¿
    double4 trivertsFaceA[3];
    double4 trivertsFaceB[3];
    double4 trivertsFaceC[3];
    double4 trivertsFaceD[3];

    trivertsFaceA[0] = Points[Triangles[t].x-1];
    trivertsFaceA[1] = Points[Triangles[t].y-1];
    trivertsFaceA[2] = Points[Triangles[t].z-1];

    trivertsFaceB[0] = Points[Triangles[t].x-1];
    trivertsFaceB[1] = Points[Triangles[t].y-1];
    trivertsFaceB[2] = Points[Triangles[t].w-1];

    trivertsFaceC[0] = Points[Triangles[t].x-1];
    trivertsFaceC[1] = Points[Triangles[t].w-1];
    trivertsFaceC[2] = Points[Triangles[t].z-1];

    trivertsFaceD[0] = Points[Triangles[t].w-1];
    trivertsFaceD[1] = Points[Triangles[t].y-1];
    trivertsFaceD[2] = Points[Triangles[t].z-1];

    return TriBoxOverlap(boxcenter, boxhalfsize, trivertsFaceA) |
	   TriBoxOverlap(boxcenter, boxhalfsize, trivertsFaceB) |
	   TriBoxOverlap(boxcenter, boxhalfsize, trivertsFaceC) |
	   TriBoxOverlap(boxcenter, boxhalfsize, trivertsFaceD);
}

double CalculateVol3(double4 a, double4 b, double4 c)
{
    return 0.5 * ( (b.x-a.x)*(c.y-a.y)- (b.y-a.y)*(c.x-a.x) );
}
	
double CalculateVol4(double4 a, double4 b, double4 c, double4 d)
{
    double x10 = a.x - d.x;
    double y10 = a.y - d.y;
    double z10 = a.z - d.z;
	           
    double x20 = b.x - d.x;
    double y20 = b.y - d.y;
    double z20 = b.z - d.z;
	           
    double x30 = c.x - d.x;
    double y30 = c.y - d.y;
    double z30 = c.z - d.z;
           
    double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
    return  detJ * 0.1666666666666666666667;
}

inline double4 calculatePositionT3(double4 a, double4 b, double4 c, double4 p)
{

    double4 N = (double4)(-1,-1,-1,0);

//     double inv_vol = native_divide(1.0f,(float)CalculateVol3(a,b,c)); //Fast
    double inv_vol = 1.0 / CalculateVol3(a,b,c);
    N = (double4)(CalculateVol3(b,c,p),CalculateVol3(c,a,p),CalculateVol3(a,b,p),0);

    return N * inv_vol;
}
	
inline double4 calculatePositionT4(double4 a, double4 b, double4 c, double4 d, double4 p)
{
    double4 N = (double4)(-1,-1,-1,0);

    double inv_vol = native_divide(1.0f,(float)CalculateVol4(a,b,c,d)); //Fast
//     double inv_vol = 1.0 / CalculateVol4(a,b,c,d);
    N = (double4)(CalculateVol4(b,d,c,p),CalculateVol4(d,a,c,p),CalculateVol4(d,b,a,p),CalculateVol4(a,b,c,p));

    return N * inv_vol;
} 