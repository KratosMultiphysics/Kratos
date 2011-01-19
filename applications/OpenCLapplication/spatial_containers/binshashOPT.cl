#pragma OPENCL EXTENSION cl_amd_fp64: enable
#pragma OPENCL EXTENSION cl_amd_printf: enable
//#pragma OPENCL EXTENSION cl_khr_fp64: enable

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

#include "Scan.cl"

/////////////////////////////////////////////////////////
////                      AUX                        ////
/////////////////////////////////////////////////////////

double CalculateVol3(double x0, double y0,
                     double x1, double y1,
                     double x2, double y2
                   )
{
    return 0.5 * ( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
}
	
double CalculateVol4(double x0, double y0, double z0,
                     double x1, double y1, double z1,
                     double x2, double y2, double z2,
                     double x3, double y3, double z3
                   )
{
    double x10 = x0 - x3;
    double y10 = y0 - y3;
    double z10 = z0 - z3;
	           
    double x20 = x1 - x3;
    double y20 = y1 - y3;
    double z20 = z1 - z3;
	           
    double x30 = x2 - x3;
    double y30 = y2 - y3;
    double z30 = z2 - z3;
           
    double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
    return  detJ * 0.1666666666666666666667;
}

double4 calculatePositionT3(double4 a, double4 b, double4 c, double4 p)
{
    double vol = CalculateVol3(a.x,a.y,b.x,b.y,c.x,c.y);
    double4 N = (double4)(-1,-1,-1,0);

    if(vol == 0.0)
      return N;

    double inv_vol = 1.0 / vol;

    N.x = CalculateVol3(b.x,b.y,c.x,c.y,p.x,p.y) * inv_vol;
    N.y = CalculateVol3(c.x,c.y,a.x,a.y,p.x,p.y) * inv_vol;
    N.z = CalculateVol3(a.x,a.y,b.x,b.y,p.x,p.y) * inv_vol;

    return N;
}
	
double4 calculatePositionT4(double4 a, double4 b, double4 c, double4 d, double4 p)
{
    double4 N = (double4)(-1,-1,-1,0);
    double vol = CalculateVol4(a.x,a.y,a.z,b.x,b.y,b.z,c.x,c.y,c.z,d.x,d.y,d.z);
    double inv_vol = 1.0 / vol;

//     if(vol == 0.0)
//       return N;

    N.x = CalculateVol4(b.x,b.y,b.z,d.x,d.y,d.z,c.x,c.y,c.z,p.x,p.y,p.z);
    N.y = CalculateVol4(d.x,d.y,d.z,a.x,a.y,a.z,c.x,c.y,c.z,p.x,p.y,p.z);
    N.z = CalculateVol4(d.x,d.y,d.z,b.x,b.y,b.z,a.x,a.y,a.z,p.x,p.y,p.z);
    N.w = CalculateVol4(a.x,a.y,a.z,b.x,b.y,b.z,c.x,c.y,c.z,p.x,p.y,p.z);

    return N * inv_vol;
} 
    
double functionDistance(double4 a, double4 b) {
    double4 temp;
    temp = (a-b) * (a-b);
    return temp.x + temp.y + temp.z;
}

int calculatePosition(double ThisCoord, 
		      int ThisDimension, 
		      __global double * N, 
		      __global double4 * MinPoint, 
		      __global double * InvCellSize 
		      )
{
    int d_index;
    
    if(ThisDimension == 2)
      d_index = (ThisCoord - MinPoint[0].z) * InvCellSize[ThisDimension];
    else if(ThisDimension == 1)
      d_index = (ThisCoord - MinPoint[0].y) * InvCellSize[ThisDimension];
    else
      d_index = (ThisCoord - MinPoint[0].x) * InvCellSize[ThisDimension];

    int index = (int)( (d_index < 0.00) ? 0.00 : d_index );
    
    return  (index > N[ThisDimension]-1) ? N[ThisDimension]-1 : index;
}

int4 calculatePosition4(double4 ThisCoord, 
		      __global double * N, 
		      __global double4 * MinPoint, 
		      __global double * InvCellSize 
		      )
{
    int4 d_index;
    
    d_index.w = -1;
    d_index.z = (ThisCoord.z - MinPoint[0].z) * InvCellSize[2];
    d_index.y = (ThisCoord.y - MinPoint[0].y) * InvCellSize[1];
    d_index.x = (ThisCoord.x - MinPoint[0].x) * InvCellSize[0];

    int4 index;
    index.w = -1;
    index.z = (int)((d_index.z < 0.00) ? 0.00 : d_index.z);
    index.y = (int)((d_index.y < 0.00) ? 0.00 : d_index.y);
    index.x = (int)((d_index.x < 0.00) ? 0.00 : d_index.x);

    index.z = (index.z > N[2]-1) ? N[2]-1 : index.z;
    index.y = (index.y > N[1]-1) ? N[1]-1 : index.y;
    index.x = (index.x > N[0]-1) ? N[0]-1 : index.x;
    
    return  index;
}

int calculateIndex(double4 ThisPoint, 
		    __global double * N, 
		    __global double4 * MinPoint, 
		    __global double * InvCellSize
		  ) 
{
    int Index = 0;
    Index += calculatePosition(ThisPoint.z,2,N,MinPoint,InvCellSize);
    Index *= N[1];
    Index += calculatePosition(ThisPoint.y,1,N,MinPoint,InvCellSize);
    Index *= N[0];
    Index += calculatePosition(ThisPoint.x,0,N,MinPoint,InvCellSize);
    return Index;
}

int calculateIndexForCell( double4 ThisIndex, 
			   __global double * N 
			  )
{
    int Index = 0;
    Index += ThisIndex.z;
    Index *= N[1];
    Index += ThisIndex.y;
    Index *= N[0];
    Index += ThisIndex.x;
    return Index;
}

int4 calculateCell(double4 ThisPoint, 
		      double Radius, 
		      __global double * N, 
		      __global double4 * MinPoint, 
		      __global double * InvCellSize )
{
    int4 Cell = calculatePosition4(ThisPoint+Radius,N,MinPoint,InvCellSize);
    return Cell;
}


/////////////////////////////////////////////////////////
////                 GENERATE BINS                   ////
/////////////////////////////////////////////////////////

__kernel void GenerateBinsA(__global double4 * points, 
			    __global int * IndexCellReference,
			    __global double * InvCellSize,
			    __global double * N,
			    __global double4 * MinPoint
			  ) 
{
    //GenerateBins
    int t = get_global_id(0);

    if(t < POINT_SIZE) 
    {	
      //Update storage counter, storing ahead
      int index = calculateIndex(points[t],N,MinPoint,InvCellSize);
      atom_inc(&IndexCellReference[index]);
    }
}

__kernel void GenerateBinsC(__global double4 * points, 
			    __global int * IndexCellReference,
			    __global double * InvCellSize,
			    __global double * N,
			    __global double4 * MinPoint,
			    __global double4 * BinsContainer
			  ) 
{
    int t = get_global_id(0);

    if(t < POINT_SIZE) 
    {	
	int index = calculateIndex(points[t],N,MinPoint,InvCellSize);
	BinsContainer[atom_dec(&IndexCellReference[index])-1] = points[t];
// 	printf("BinsTriangle: %f\n",points[t].w);
    }
}


/////////////////////////////////////////////////////////
////               SEARCH FUNCTIONS                  ////
/////////////////////////////////////////////////////////


///////////////////////////////////SearchInRadius///////////////////////////////////
__kernel void SearchInRadiusMultiple(__global int * IndexCellReference,
			     __global double4 * BinsContainer,
			     __global double * InvCellSize,
			     __global double * N,
			     __global const double * radius,
			     __global const double * radius2,
			     __global double4 * point,
			     __global double4 * MinPoint,
			     __global int * outdata,
			     __global int * results,
			     __global int * maxResults
			    ) 
{
    int ip = get_global_id(0);

    __private double4 pointIp = point[ip];

    __private int4 cellBegin = calculateCell(pointIp,-radius[0],N,MinPoint,InvCellSize);
    __private int4 cellEnd   = calculateCell(pointIp, radius[0],N,MinPoint,InvCellSize);
    results[ip] = 0;

    __private int4 auxData; //x-> index for IndexCellReference
			    //y-> lowIndex of BinsContainer
			    //z-> highIndex of BinsContainer

    for(size_t i = cellBegin.x; i <= cellEnd.x; i++) 
    {
	for(size_t j = cellBegin.y; j <= cellEnd.y; j++) 
	{
	    for(size_t k = cellBegin.z; k <= cellEnd.z; k++) 
	    {
		auxData.x = calculateIndexForCell((double4)(i,j,k,-1),N);
		
		auxData.y = IndexCellReference[auxData.x];
		auxData.z = IndexCellReference[auxData.x+1];

		for(size_t l = auxData.y; l < auxData.z; l++) 
		{
		    if(functionDistance(pointIp,BinsContainer[l]) < radius2[0]) 
		    {
			if(results[ip] < maxResults[0])
			    outdata[mad_hi(ip,maxResults[0],results[ip])] = BinsContainer[l].w;
			results[ip]++;
		    }
		}
	    }
	}
    }
}

///////////////////////////////////SearchNearest///////////////////////////////////
__kernel void SearchNearestMultiple(__global int * IndexCellReference,
				    __global double4 * BinsContainer,
				    __global double * InvCellSize,
				    __global double * N,
				    __global const double * radius,
				    __global const double * radius2,
				    __global double4 * point,
				    __global double4 * MinPoint,
				    __global double * distance,
				    __global int * results
				    ) 
{
    int ip = get_global_id(0);
    int t = get_local_id(0);

    __private double4 pointIp = point[ip];

    __private int4 cellBegin = calculateCell(pointIp,-radius[0],N,MinPoint,InvCellSize);
    __private int4 cellEnd   = calculateCell(pointIp, radius[0],N,MinPoint,InvCellSize);

    distance[ip] = INFINITY;

    __private int4 auxData; //x-> index for IndexCellReference
			    //y-> lowIndex of BinsContainer
			    //z-> highIndex of BinsContainer

    for(size_t i = cellBegin.x; i <= cellEnd.x; i++) 
    {
	for(size_t j = cellBegin.y; j <= cellEnd.y; j++) 
	{
	    for(size_t k = cellBegin.z; k <= cellEnd.z; k++) 
	    {
		auxData.x = calculateIndexForCell((double4)(i,j,k,-1),N);
		
		auxData.y = IndexCellReference[auxData.x];
		auxData.z = IndexCellReference[auxData.x+1];
		
		for(size_t l = auxData.y; l < auxData.z; l++) 
		{
		    __private double thisDist = functionDistance(pointIp,BinsContainer[l]); 
		    if(isless(thisDist,distance[ip]))
		    {
			distance[ip] = thisDist;
	                results[ip] = BinsContainer[l].w;
		    }
		}
	    }
	}
    }
}

///////////////////////////////////SearchTriangle///////////////////////////////////
__kernel void SearchTriangle2D(__global int * IndexCellReference,
			     __global double4 * BinsContainer,
			     __global double4 * PointsTriangle,
			     __global int4 * Triangles,
			     __global double * InvCellSize,
			     __global double * N,
			     __global double * radius,
			     __global double * radius2,
			     __global double4 * point,
			     __global double4 * MinPoint,
			     __global double * uVector,
			     __global double4 * Nresults,
			     __global int * pSize
			    ) 
{
    int ip = get_global_id(0);

    if(ip < pSize[0])
    {
	__private bool F = false;;
	__private int4 Triangle;

	__private double4 pointIp = point[ip];
	__private double4 fN;

	int4 cellBegin = calculateCell(pointIp,-radius[0],N,MinPoint,InvCellSize);
	int4 cellEnd   = calculateCell(pointIp, radius[0],N,MinPoint,InvCellSize);

	for(size_t i = cellBegin.x; !F && (i <= cellEnd.x); i++) 
	{
	    for(size_t j = cellBegin.y; !F && (j <= cellEnd.y); j++) 
	    {
		for(size_t k = cellBegin.z; !F && (k <= cellEnd.z); k++) 
		{
		    int index = calculateIndexForCell((double4)(i,j,k,-1),N);
		    
		    int loIndex = IndexCellReference[index];
		    int hiIndex = IndexCellReference[index+1];


		    for(size_t l = loIndex; !F && (l < hiIndex); l++) 
		    {
			if(functionDistance(pointIp,BinsContainer[l]) < radius2[0])
			{
			    Triangle = Triangles[(int)BinsContainer[l].w];

			    fN = calculatePositionT3(PointsTriangle[Triangle.x-1],
						     PointsTriangle[Triangle.y-1],
						     PointsTriangle[Triangle.z-1],
					             pointIp);

			    F = fN.x > -0.00001 && 
				fN.y > -0.00001 && 
				fN.z > -0.00001 &&
				fN.x <  1.00001 && 
				fN.y <  1.00001 && 
				fN.z <  1.00001 ; 
			}
		    }
		}
	    }
	}

	Nresults[ip] = F ? (fN * uVector[ip]) : - 1;
    }
}

///////////////////////////////////SearchTriangle///////////////////////////////////
//Prefetch all index. si va se puede aprovechar con los demas kernels.
// __kernel void SearchTriangle3DA(__global int * IndexCellReference,
// 			        __global double * InvCellSize,
// 			        __global double * N,
// 			        __global double * radius,
// 			        __global double4 * point,
// 			        __global double4 * MinPoint,
// 			        __global int2 * searchIndex,
// 			        __global int * pSize
// 			    ) 
// {
//     size_t ip = get_global_id(0);
//     int results = 0;
// 
//     if(ip < pSize[0])
//     {
// 	int4 cellBegin = calculateCell(point[ip],-radius[0],N,MinPoint,InvCellSize);
// 	int4 cellEnd   = calculateCell(point[ip], radius[0],N,MinPoint,InvCellSize);
// 
// 	for(int i = cellBegin.x; i <= cellEnd.x; i++) 
// 	{
// 	    for(int j = cellBegin.y; j <= cellEnd.y; j++) 
// 	    {
// 		for(int k = cellBegin.z; k <= cellEnd.z; k++) 
// 		{
// 		    int index = calculateIndexForCell((double4)(i,j,k,-1),N);
// 
// 		    searchIndex[ip*300+results].x = IndexCellReference[index];
// 		    searchIndex[ip*300+results].y = IndexCellReference[index+1];
// 		    results++;
// 		}
// 	    }
// 	}
// 
//     searchIndex[ip*300+results].x = -1;
//     searchIndex[ip*300+results].y = -1;
// 
// // 	searchIndex[ip].x = IndexCellReference[calculateIndexForCell((double4)(cellBegin.x,cellBegin.y,cellBegin.z,-1),N)];
// // 	searchIndex[ip].y = IndexCellReference[calculateIndexForCell((double4)(cellEnd.x  ,cellEnd.y  ,cellEnd.z,  -1),N)];
//     }
// }

// __kernel void SearchTriangle3DB(__global double4 * BinsContainer,
// 			       __global double4 * PointsTriangle,
// 			       __global int4 * Triangles,
// 			       __global double * radius2,
// 			       __global double4 * point,
// 			       __global double * uVector,
// 			       __global double4 * Nresults,
// 			       __global int2 * searchIndex,
// 			       __global int * pSize
// 			    ) 
// {
//     size_t ip = get_global_id(0);
// 
//     if(ip < pSize[0])
//     {
// 	short F = 0;
// 	double4 fN;
// 
// int i = 0;
// int index = ip*300+i;
// while(searchIndex[index].x != -1)
// {
// 	for(size_t l = searchIndex[index].x; !F && (l < searchIndex[index].y); l++) 
// 	{
// 	    if(functionDistance(point[ip],BinsContainer[l]) < radius2[0])
// 	    {
// 		int4 triangleIndex = Triangles[(int)BinsContainer[l].w];
// 
// 		fN = calculatePositionT4(PointsTriangle[triangleIndex.x-1],
// 					 PointsTriangle[triangleIndex.y-1],
// 					 PointsTriangle[triangleIndex.z-1],
// 					 PointsTriangle[triangleIndex.w-1],
// 					 point[ip]);
// 
// 		F = fN.x > -0.00001f && 
// 		    fN.y > -0.00001f && 
// 		    fN.z > -0.00001f &&
// 		    fN.w > -0.00001f &&
// 		    fN.x <  1.00001f && 
// 		    fN.y <  1.00001f && 
// 		    fN.z <  1.00001f && 
// 		    fN.w <  1.00001f; 
// 	    }
// 	}
// i++;
// index = ip*300+i;
// }
// 
//  	Nresults[ip] = F ? (fN * uVector[ip]) : -1;
// // Nresults[ip].x = searchIndex[ip].x;
// // Nresults[ip].y = searchIndex[ip].y;
//     }
// }

/////////////////////////////////////////////////////////
////                      TEST                       ////
/////////////////////////////////////////////////////////

///////////////////////////////////Not sure if this code is correct///////////////////////////////////
__kernel void SearchNearestMultipleCubic(__global int * IndexCellReference,
				    __global double4 * BinsContainer,
				    __global double * InvCellSize,
				    __global double * N,
				    __global double * radius,
				    __global double * radius2,
				    __global double4 * point,
				    __global double4 * MinPoint,
				    __global double * distance,
				    __global int * results
				    ) 
{
    int ip = get_global_id(0);

    __private double4 pointIp = point[ip];

    __private int4 cellBegin = calculateCell(pointIp,-radius[0],N,MinPoint,InvCellSize);
    __private int4 cellEnd   = calculateCell(pointIp, radius[0],N,MinPoint,InvCellSize);
    __private int4 cellPoint = calculateCell(pointIp, 0,N,MinPoint,InvCellSize);

    distance[ip] = INFINITY;

    bool Found = false;
    bool partialFound = false;

    int cellcuberadius = 0;
    int maxradiusx = max(cellEnd.x - cellPoint.x, cellPoint.x - cellBegin.x);
    int maxradiusy = max(cellEnd.y - cellPoint.y, cellPoint.y - cellBegin.y);
    int maxradiusz = max(cellEnd.z - cellPoint.z, cellPoint.z - cellBegin.z);

    int maxradius = max(max(maxradiusx,maxradiusy),maxradiusz);

    while(!Found && cellcuberadius <= maxradius)
    {
      int kini = cellPoint.z - cellcuberadius;
      int kend = cellPoint.z + cellcuberadius;
      int k = kini;

      while( k <= kend )
      {
	int jini = cellPoint.y - cellcuberadius;
	int jend = cellPoint.y + cellcuberadius;
	int j = jini;
	while( j <= jend )
	{
          int iini = cellPoint.x - cellcuberadius;
	  int iend = cellPoint.x + cellcuberadius;
	  int i = iini;
	  while( i <= iend )
	  {
	    if( (
		  (k == kini) || 
		  (k == kend) || 
		  (j == 0)    || 
		  (j == jend) || 
		  (i == 0)    || 
		  (i == iend)                        ) &&
		( k >= cellBegin.z && k <= cellEnd.z ) &&
		( j >= cellBegin.y && j <= cellEnd.y ) &&
		( i >= cellBegin.x && i <= cellEnd.x ) )
	    {
	      __private int index = calculateIndexForCell((double4)(i,j,k,-1),N);

	      __private int loIndex = IndexCellReference[index];
	      __private int hiIndex = IndexCellReference[index+1];
	      
	      for(size_t l = loIndex; l < hiIndex; l++) 
	      {
		  __private double thisDist = functionDistance(pointIp,BinsContainer[l]); 
		  if(isless(thisDist,distance[ip]))
		  {
		      distance[ip] = thisDist;
		      results[ip] = BinsContainer[l].w;
		      partialFound = true;
		  }
	      }
	    } // end if
	    i++;
	  } // endl while i
	  j++;
	} // endl while j
	k++;
      } // endl while k
      if( partialFound && cellcuberadius!=0 )
	Found = true; //congratulations! you found it!!!
      cellcuberadius++;
    } //endl while cuberadius
}

//Triangle3D ORIGIANAL
__kernel void SearchTriangle3D2(__global int * IndexCellReference,
			       __global double4 * BinsContainer,
			       __global double4 * PointsTriangle,
			       __global int4 * Triangles,
			       __global double * InvCellSize,
			       __global double * N,
			       __global double * radius,
			       __global double * radius2,
			       __global double4 * point,
			       __global double4 * MinPoint,
			       __global double * uVector,
			       __global double4 * Nresults,
			       __global int * pSize
			    ) 
{
    size_t ip = get_global_id(0);

    if(ip < pSize[0])
    {
	short F = 0;
	double4 fN;

	int4 cellBegin = calculateCell(point[ip],-radius[0],N,MinPoint,InvCellSize);
	int4 cellEnd   = calculateCell(point[ip], radius[0],N,MinPoint,InvCellSize);

	for(int i = cellBegin.x; !F && (i <= cellEnd.x); i++) 
	{
	    for(int j = cellBegin.y; !F && (j <= cellEnd.y); j++) 
	    {
		for(int k = cellBegin.z; !F && (k <= cellEnd.z); k++) 
		{
		    int index = calculateIndexForCell((double4)(i,j,k,-1),N);
		    
		    int loIndex = IndexCellReference[index];
		    int hiIndex = IndexCellReference[index+1];

		    for(size_t l = loIndex; !F && (l < hiIndex); l++) 
		    {
			if(functionDistance(point[ip],BinsContainer[l]) < radius2[0])
			{
			    int4 triangleIndex = Triangles[(int)BinsContainer[l].w];

			    fN = calculatePositionT4(PointsTriangle[triangleIndex.x-1],
						     PointsTriangle[triangleIndex.y-1],
						     PointsTriangle[triangleIndex.z-1],
						     PointsTriangle[triangleIndex.w-1],
						     point[ip]);

			    F = fN.x > -0.00001f && 
				fN.y > -0.00001f && 
				fN.z > -0.00001f &&
				fN.w > -0.00001f &&
				fN.x <  1.00001f && 
				fN.y <  1.00001f && 
				fN.z <  1.00001f && 
				fN.w <  1.00001f; 
			}
		    }
		}
	    }
	}

	Nresults[ip] = F ? (fN * uVector[ip]) : -1;
    }
}

//Triangle3D 
__kernel void SearchTriangle3D(__global int * IndexCellReference,
			       __global double4 * BinsContainer,
			       __global double4 * PointsTriangle,
			       __global int4 * Triangles,
			       __global double * InvCellSize,
			       __global double * N,
			       __global double * radius,
			       __global double * radius2,
			       __global double4 * point,
			       __global double4 * MinPoint,
			       __global double * uVector,
			       __global double4 * Nresults,
			       __global int * pSize
			       __local 
			    ) 
{
    size_t ip = get_global_id(0);

    if(ip < pSize[0])
    {
	short F = 0;
	double4 fN;

	int4 cellBegin = calculateCell(point[ip],-radius[0],N,MinPoint,InvCellSize);
	int4 cellEnd   = calculateCell(point[ip], radius[0],N,MinPoint,InvCellSize);

	for(int i = cellBegin.x; !F && (i <= cellEnd.x); i++) 
	{
	    for(int j = cellBegin.y; !F && (j <= cellEnd.y); j++) 
	    {
		for(int k = cellBegin.z; !F && (k <= cellEnd.z); k++) 
		{
		    int index = calculateIndexForCell((double4)(i,j,k,-1),N);
		    
		    int loIndex = IndexCellReference[index];
		    int hiIndex = IndexCellReference[index+1];

		    //if(hiIndex-loIndex > 15) hiIndex = loIndex + 5;

		    for(size_t l = loIndex; !F && (l < hiIndex) ; l++) 
		    {
			if(functionDistance(point[ip],BinsContainer[l]) < radius2[0])
			{
			    int4 triangleIndex = Triangles[(int)BinsContainer[l].w]-1;

			    fN = calculatePositionT4(PointsTriangle[triangleIndex.x],
						     PointsTriangle[triangleIndex.y],
						     PointsTriangle[triangleIndex.z],
						     PointsTriangle[triangleIndex.w],
						     point[ip]);

			    F =  fN.x > -0.00001 && 
				 fN.y > -0.00001 && 
				 fN.z > -0.00001 &&
				 fN.w > -0.00001 &&
				 fN.x <  1.00001 && 
				 fN.y <  1.00001 && 
				 fN.z <  1.00001 && 
				 fN.w <  1.00001; 
			}
		    }
		}
	    }
	}

	Nresults[ip] = F ? (fN * uVector[ip]) : -1;
    }
}

