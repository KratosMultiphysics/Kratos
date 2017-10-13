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

    double inv_vol = native_divide(1.0f,(float)CalculateVol3(a,b,c)); //Fast
//     double inv_vol = 1.0 / CalculateVol3(a,b,c);
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
    
double functionDistance(double4 a, double4 b) 
{
    double4 temp;

    temp = (a-b) * (a-b);

    return temp.x + temp.y + temp.z;
}

int calculatePosition(double ThisCoord, 
		      int ThisDimension, 
		      __constant double * N, 
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
		      __constant double * N, 
		      __global double4 * MinPoint, 
		      __global double * InvCellSize 
		      )
{
    int4 d_index = convert_int4_rte((ThisCoord - MinPoint[0]) * (double4)(InvCellSize[0],InvCellSize[1],InvCellSize[2],0));
    int4 index = (d_index < (int4)(0,0,0,0)) ? (int4)(0,0,0,0) : d_index;
    return index = (index > (int4)(N[0]-1,N[1]-1,N[2]-1,0)) ? (int4)(N[0]-1,N[1]-1,N[2]-1,0) : index;
}

int calculateIndex(double4 ThisPoint, 
		    __constant double * N, 
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
			   __constant double * N 
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
		      __constant double * N, 
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
			    __constant double * N,
			    __global double4 * MinPoint,
			    __constant int * amount
			  ) 
{
    int t = get_global_id(0);

    if(t < amount[0]) 
    {	
      int index = calculateIndex(points[t],N,MinPoint,InvCellSize);
      atom_inc(&IndexCellReference[index]);
    }
}

__kernel void GenerateBinsC(__global double4 * points, 
			    __global int * IndexCellReference,
			    __global double * InvCellSize,
			    __constant double * N,
			    __global double4 * MinPoint,
			    __global double4 * BinsContainer,
			    __constant int * amount
			  ) 
{
    int t = get_global_id(0);

    if(t < amount[0]) 
    {	
	int index = calculateIndex(points[t],N,MinPoint,InvCellSize);
	BinsContainer[atom_dec(&IndexCellReference[index+1])-1] = points[t];
    }
}


/////////////////////////////////////////////////////////
////               SEARCH FUNCTIONS                  ////
/////////////////////////////////////////////////////////


///////////////////////////////////SearchInRadius///////////////////////////////////
__kernel void SearchInRadiusMultiple(__global int * IndexCellReference,
			     __global double4 * BinsContainer,
			     __global double * InvCellSize,
			     __constant double * N,
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
				    __constant double * N,
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
			     __constant double * N,
			     __constant double * radius,
			     __global double4 * point,
			     __global double4 * MinPoint,
			     __global double4 * Nresults,
			     __constant int * pSize,
			     __local int4 * cellBegin,
			     __local int4 * cellEnd
			    ) 
{
    int ip = get_global_id(0);
    int il = get_local_id(0);

    if(ip < pSize[0])
    {
	int Found = 0;;
	double4 fN;

	cellBegin[il] = calculateCell(point[ip],-radius[0],N,MinPoint,InvCellSize);
	cellEnd[il]   = calculateCell(point[ip], radius[0],N,MinPoint,InvCellSize);

	barrier(CLK_LOCAL_MEM_FENCE);

	for(int k = cellBegin[il].z; !Found && (k <= cellEnd[il].z); k++) 
	{
	    for(int j = cellBegin[il].y; !Found && (j <= cellEnd[il].y); j++) 
	    {
		for(int i = cellBegin[il].x; !Found && (i <= cellEnd[il].x); i++) 
		{
		    int index = calculateIndexForCell((double4)(i,j,k,-1),N);
		    
		    int loIndex = IndexCellReference[index];
		    int hiIndex = IndexCellReference[index+1];

		    for(int l = loIndex; !Found && (l < hiIndex); l++) 
		    {
			if(functionDistance(point[ip],BinsContainer[l]) < radius[1])
			{
			    int4 triangleIndex = Triangles[(int)BinsContainer[l].w];

			    fN = calculatePositionT3(PointsTriangle[triangleIndex.x-1],
						     PointsTriangle[triangleIndex.y-1],
						     PointsTriangle[triangleIndex.z-1],
					             point[ip]);

			    Found = fN.x >= 0.0 && fN.y >= 0.0 && fN.z >= 0.0 &&
				    fN.x <= 1.0 && fN.y <= 1.0 && fN.z <= 1.0; 
			}
		    }
		}
	    }
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	Nresults[ip] = Found ? (fN * point[ip].w) : - 1;
    }
}

__kernel void SearchTriangle3D(__global int * IndexCellReference,
			       __global double4 * BinsContainer,
			       __global double4 * PointsTriangle,
			       __global int4 * Triangles,
			       __global double * InvCellSize,
			       __constant double * N,
			       __constant double * radius,
			       __global double4 * point,
			       __global double4 * MinPoint,
			       __global double4 * Nresults,
			       __constant int * pSize,
			       __local int4 * cellBegin,
			       __local int4 * cellEnd
			    ) 
{
    int ip = get_global_id(0);
    int il = get_local_id(0);

    if(ip < pSize[0])
    {
	int Found = 0;
	double4 fN;

	cellBegin[il] = calculateCell(point[ip],-radius[0],N,MinPoint,InvCellSize);
	cellEnd[il]   = calculateCell(point[ip], radius[0],N,MinPoint,InvCellSize);

	barrier(CLK_LOCAL_MEM_FENCE);

	for(int k = cellBegin[il].z; !Found && (k <= cellEnd[il].z); k++) 
	{
	    for(int j = cellBegin[il].y; !Found && (j <= cellEnd[il].y); j++) 
	    {
		for(int i = cellBegin[il].x; !Found && (i <= cellEnd[il].x); i++) 
		{
		    int index = calculateIndexForCell((double4)(i,j,k,-1),N);
		    
		    int loIndex = IndexCellReference[index+1];
		    int hiIndex = IndexCellReference[index+2];

		    for(size_t l = loIndex; !Found && (l < hiIndex); l++) 
		    {
			if(functionDistance(point[ip],BinsContainer[l]) < radius[1])
			{
			    int4 triangleIndex = Triangles[(int)BinsContainer[l].w];

			    fN = calculatePositionT4(PointsTriangle[triangleIndex.x-1],
						     PointsTriangle[triangleIndex.y-1],
						     PointsTriangle[triangleIndex.z-1],
						     PointsTriangle[triangleIndex.w-1],
						     point[ip]);

			    Found =  fN.x >= 0.0 && fN.y >= 0.0 && fN.z >= 0.0 && fN.w >= 0.0 &&
				     fN.x <= 1.0 && fN.y <= 1.0 && fN.z <= 1.0 && fN.w <= 1.0; 
			}
		    }
		}
	    }
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	Nresults[ip] = Found ? (fN * point[ip].w) : -1;
    }
}


/////////////////////////////////////////////////////////
////                      TEST                       ////
/////////////////////////////////////////////////////////

///////////////////////////////////Not sure if this code is correct///////////////////////////////////
__kernel void SearchNearestMultipleCubic(__global int * IndexCellReference,
				    __global double4 * BinsContainer,
				    __global double * InvCellSize,
				    __constant double * N,
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

