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

// #ifndef NVIDIA

	#pragma OPENCL EXTENSION cl_nv_compiler_options : enable
	#pragma OPENCL EXTENSION cl_nv_pragma_unroll : enable

// #endif

#include "opencl_operations.cl"

/////////////////////////////////////////////////////////

__constant double one_third = 1.0 / 3.0;
__constant double one_sixt = 1.0 / 6.0;
__constant double two_third = 2.0 * 1.0 / 3.0;
__constant double one_mid = 1.0 / 2.0;

inline void updateV(double4 fN,
	     __global double4 * V) 
{
    int ip = get_global_id(0);
    V[ip] = V[ip] * fN;
}

inline void updateF(double4 fN,
	     __global double4 * F)
{
    int ip = get_global_id(0);
    F[ip] = F[ip] * fN; 
}

inline void updateR(double4 fN,
	     __global double4 * R)
{
    int ip = get_global_id(0);
    R[ip] = R[ip] * fN; 
}

void updateY(double dT,
	     __global double * T,
	     __global double * Y,
	     __global double4 * R)
{

}

void updatePos(double4 fN,
	       __global double4 * point,
	       __global double4 * V)
{
    int ip = get_global_id(0);
    double4 vp = V[ip] * fabs(fN);
    double veuler = vp.x + vp.y + vp.z + vp.w;
    point[ip].x += veuler*0.005;
    point[ip].w = 1;
}

/////////////////////////////////////////////////////////

void clearBoundingBox(double4 point, double4 *min, double4 *max)
{
    (*min) = (*max) = point;
}

void addToBoundingBox(double4 point, double4 *min, double4 *max)
{
    *min = point < *min ? point : *min;
    *max = point > *max ? point : *max;
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

int calculateIndexForCell( int4 ThisIndex, 
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
    int4 Cell;
    Cell.x = calculatePosition(ThisPoint.x+Radius,0,N,MinPoint,InvCellSize);
    Cell.y = calculatePosition(ThisPoint.y+Radius,1,N,MinPoint,InvCellSize);
    Cell.z = calculatePosition(ThisPoint.z+Radius,2,N,MinPoint,InvCellSize);
    return Cell;
}


/////////////////////////////////////////////////////////

__kernel void GenerateBinsObjectsA(__global double4 * Points,
				   __global int4 * Triangles,
				   __global int * IndexCellReference,
				   __global double * InvCellSize,
				   __constant double * N,
				   __global double4 * MinPoint,
				   int size
			  ) 
{
    int t = get_global_id(0);

    if(t < size && Triangles[t].x != -1) 
    {	
	double4 minBox, maxBox;
	double4 minCell, maxCell;
	double4 low, high;

	clearBoundingBox(Points[Triangles[t].x-1],&low,&high);
	addToBoundingBox(Points[Triangles[t].y-1],&low,&high);
	addToBoundingBox(Points[Triangles[t].z-1],&low,&high);
//  	addToBoundingBox(Points[Triangles[t].w-1],&low,&high);

	int4 cellBegin = calculateCell(low,0,N,MinPoint,InvCellSize);
	int4 cellEnd  = calculateCell(high,0,N,MinPoint,InvCellSize);

	if(HasIntersection(Points,Triangles,low,high))
	{
	    for(int i = cellBegin.x; i <= cellEnd.x; i++)
	    {
		for(int j = cellBegin.y; j <= cellEnd.y; j++)
		{
		    for(int k = cellBegin.z; k <= cellEnd.z; k++)
		    {
			atom_inc(&IndexCellReference[calculateIndexForCell((int4)(i,j,k,-1),N)]);
		    }
		}
	    }
	}
    }
}

__kernel void GenerateBinsC(__global double4 * Points,
			    __global int4 * Triangles,
			    __global int * IndexCellReference,
			    __global double * InvCellSize,
			    __constant double * N,
			    __global double4 * MinPoint,
			    __global int * BinsObjectContainer,
			    int size
			  ) 
{
    int t = get_global_id(0);

    if(t < size && Triangles[t].x != -1) 
    {	
	double4 minBox, maxBox;
	double4 minCell, maxCell;
	double4 low, high;

	clearBoundingBox(Points[Triangles[t].x-1],&low,&high);
	addToBoundingBox(Points[Triangles[t].y-1],&low,&high);
	addToBoundingBox(Points[Triangles[t].z-1],&low,&high);
//  	addToBoundingBox(Points[Triangles[t].w-1],&low,&high);

	int4 cellBegin = calculateCell(low,0,N,MinPoint,InvCellSize);
	int4 cellEnd  = calculateCell(high,0,N,MinPoint,InvCellSize);

	if(HasIntersection(Points,Triangles,low,high))
	{
	    for(int i = cellBegin.x; i <= cellEnd.x; i++)
	    {
		for(int j = cellBegin.y; j <= cellEnd.y; j++)
		{
		    for(int k = cellBegin.z; k <= cellEnd.z; k++)
		    {
			int index = calculateIndexForCell((int4)(i,j,k,-1),N);
			int pos = atom_dec(&IndexCellReference[index+1])-1;

			BinsObjectContainer[pos] = t;
		    }
		}
	    }
	}
    }
}

/////////////////////MoveParticles//////////////////////
__kernel void Move(__global int * IndexCellReference,
		   __global int * BinsObjectContainer,
		   __global double4 * PointsTriangle,
		   __global int4 * Triangles,
		   __global double * InvCellSize,
		   __constant double * N,
		   __global double4 * MinPoint,
		   __global double4 * point,
		   __global double4 * pointVelocity,
		   __global double4 * pointDisplace,
		   __global double4 * pointForce,
		   __global double4 * nodesV,
		   __global double4 * nodesF,
		   __global double4 * nodesP,
		   __global double4 * body_force,
		   __global double4 * pointVelocityOld,
		   __global double4 * pointDisplaceOld,
		   double idensity,
		   double small_dt,
		   double subSteps,
		   int use_eulerian
		  )
{

      int gid = get_global_id(0);

      int4 triangleIndex = 0;

      double4 eulerian_vel = (double4)(0.0,0.0,0.0,0.0);

      pointVelocityOld[gid] = pointVelocity[gid];
      pointDisplaceOld[gid] = pointDisplace[gid];
      point[gid].xyz += pointDisplace[gid].xyz;

      for(int subStep = 0; subStep < subSteps; subStep++)
      {
	int index = calculateIndex(point[gid],N,MinPoint,InvCellSize);
		    
	int loIndex = IndexCellReference[index+1];
	int hiIndex = IndexCellReference[index+2];

	int Found = 0;

	double4 fN = (double4)(0.0,0.0,0.0,0.0);

	for(int l = loIndex; !Found && (l < hiIndex ); l++) 
	{
	    triangleIndex = Triangles[BinsObjectContainer[l]];

	    fN = calculatePositionT3(PointsTriangle[triangleIndex.x-1],
				     PointsTriangle[triangleIndex.y-1],
				     PointsTriangle[triangleIndex.z-1],
// 				     PointsTriangle[triangleIndex.w-1],
				     point[gid]);

	    Found = fN.x >= 0.0 && fN.y >= 0.0 && fN.z >= 0.0 /*&& fN.w >= 0.0*/ &&
		    fN.x <= 1.0 && fN.y <= 1.0 && fN.z <= 1.0 /*&& fN.w <= 1.0*/;
	}

	if (Found) {
	    eulerian_vel = ((nodesV[triangleIndex.x-1] * fN.x) +
			    (nodesV[triangleIndex.y-1] * fN.y) +
			    (nodesV[triangleIndex.z-1] * fN.z));

	    fN *= idensity;

	    pointForce[gid] = (double4)(0.0,0.0,0.0,0.0);
	    pointForce[gid] += fN.x * nodesF[triangleIndex.x-1];
	    pointForce[gid] += fN.y * nodesF[triangleIndex.y-1];
	    pointForce[gid] += fN.z * nodesF[triangleIndex.z-1];

	    pointVelocity[gid].xyz = (use_eulerian && !subStep) ? eulerian_vel.xyz : pointVelocity[gid].xyz;
	    pointVelocityOld[gid].xyz = (use_eulerian && !subStep) ? eulerian_vel.xyz : pointVelocityOld[gid].xyz;

	    pointVelocity[gid].xyz += ((body_force[0].xyz - fN.x * nodesP[triangleIndex.x-1].xyz - fN.y * nodesP[triangleIndex.y-1].xyz - fN.z * nodesP[triangleIndex.z-1].xyz + pointForce[gid].xyz) * small_dt); 

	    point[gid].xyz 	   += eulerian_vel.xyz * small_dt;
	    pointDisplace[gid].xyz += eulerian_vel.xyz * small_dt;
	} else {
	    point[gid].w = min(-point[gid].w,point[gid].w);
	}
    }
}

__kernel void Move2(__global double4 * point,
		    __global double4 * pointVelocity,
		    __global double4 * pointDisplace,
		    __global double4 * pointVelocityOld,
		    __global double4 * pointDisplaceOld,
		    double dt
		   )
{
    int gid = get_global_id(0);

    double4 stp_dis = pointDisplace[gid] - pointDisplaceOld[gid];

    double norm_d = sqrt((stp_dis.x * stp_dis.x) + (stp_dis.y * stp_dis.y));
    double norm_v = sqrt((pointVelocityOld[gid].x * pointVelocityOld[gid].x) + (pointVelocityOld[gid].y * pointVelocityOld[gid].y));

    if(norm_d*10.0 < norm_v*dt) {
	point[gid].w = min(-point[gid].w,point[gid].w);
    }
}

__kernel void resetCounter(__global int * gElementCounter) {
      gElementCounter[get_global_id(0)] = 0;
}

__kernel void calculateField(__global int * IndexCellReference,
		   __global int * BinsObjectContainer,
		   __global double4 * PointsTriangle,
		   __global int4 * Triangles,
		   __global double * InvCellSize,
		   __constant double * N,
		   __constant double * radius,
		   __global double4 * MinPoint,
		   __global double4 * point,
		   __global double4 * pointVelocity,
		   __global double4 * nodesV,
		   __global int * gIndex,
		   __global double4 * gFn,
		   int size,
		   __local int * l
		  )
{
    int gid = get_global_id(0);

    if(gid < size) {

	gFn[gid]    = 0;
	gIndex[gid] = -1;

	if(point[gid].w >= 0) 
	{
	    double4 fN;
	    int4 triangleIndex;

	    int index = calculateIndex(point[gid],N,MinPoint,InvCellSize);
			
	    int loIndex = IndexCellReference[index+1];
	    int hiIndex = IndexCellReference[index+2];

	    int Found = 0;

	    for(l[get_local_id(0)] = loIndex; !Found && (l[get_local_id(0)] < hiIndex); l[get_local_id(0)]++) 
	    {
		triangleIndex = Triangles[BinsObjectContainer[l[get_local_id(0)]]];

		fN = calculatePositionT3(PointsTriangle[triangleIndex.x-1],
					 PointsTriangle[triangleIndex.y-1],
					 PointsTriangle[triangleIndex.z-1],
    //   			 	 PointsTriangle[triangleIndex.w-1],
					 point[gid]);

		Found = fN.x >= 0.0 && fN.y >= 0.0 && fN.z >= 0.0 /*&& fN.w >= 0.0*/ &&
			fN.x <= 1.0 && fN.y <= 1.0 && fN.z <= 1.0 /*&& fN.w <= 1.0*/;
	    }

	    gFn[gid]    = fN * Found;
	    int lindex = BinsObjectContainer[l[get_local_id(0)]-1] * Found + -1 * !Found;
	    gIndex[gid] = lindex;
	}
    }
}

/////////////////////////////////////////////////////////
////               PARALLEL SCAN                     ////
/////////////////////////////////////////////////////////


/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */


//Passed down with -D option on clBuildProgram
//Must be a power of two
//#define WORKGROUP_SIZE 256



////////////////////////////////////////////////////////////////////////////////
// Scan codelets
////////////////////////////////////////////////////////////////////////////////
#if(1)
    //Naive inclusive scan: O(N * log2(N)) operations
    //Allocate 2 * 'size' local memory, initialize the first half
    //with 'size' zeros avoiding if(pos >= offset) condition evaluation
    //and saving instructions
    inline uint scan1Inclusive(uint idata, __local uint *l_Data, uint size){
        uint pos = 2 * get_local_id(0) - (get_local_id(0) & (size - 1));
        l_Data[pos] = 0;
        pos += size;
        l_Data[pos] = idata;

        for(uint offset = 1; offset < size; offset <<= 1){
            barrier(CLK_LOCAL_MEM_FENCE);
            uint t = l_Data[pos] + l_Data[pos - offset];
            barrier(CLK_LOCAL_MEM_FENCE);
            l_Data[pos] = t;
        }

        return l_Data[pos];
    }

    inline uint scan1Exclusive(uint idata, __local uint *l_Data, uint size){
        return scan1Inclusive(idata, l_Data, size) - idata;
    }

#else
    #define LOG2_WARP_SIZE 5U
    #define      WARP_SIZE (1U << LOG2_WARP_SIZE)

    //Almost the same as naive scan1Inclusive but doesn't need barriers
    //and works only for size <= WARP_SIZE
    inline uint warpScanInclusive(uint idata, __local uint *l_Data, uint size){
        uint pos = 2 * get_local_id(0) - (get_local_id(0) & (size - 1));
        l_Data[pos] = 0;
        pos += size;
        l_Data[pos] = idata;

        if(size >=  2) l_Data[pos] += l_Data[pos -  1];
        if(size >=  4) l_Data[pos] += l_Data[pos -  2];
        if(size >=  8) l_Data[pos] += l_Data[pos -  4];
        if(size >= 16) l_Data[pos] += l_Data[pos -  8];
        if(size >= 32) l_Data[pos] += l_Data[pos - 16];

        return l_Data[pos];
    }

    inline uint warpScanExclusive(uint idata, __local uint *l_Data, uint size){
        return warpScanInclusive(idata, l_Data, size) - idata;
    }

    inline uint scan1Inclusive(uint idata, __local uint *l_Data, uint size){
        if(size > WARP_SIZE){
            //Bottom-level inclusive warp scan
            uint warpResult = warpScanInclusive(idata, l_Data, WARP_SIZE);

            //Save top elements of each warp for exclusive warp scan
            //sync to wait for warp scans to complete (because l_Data is being overwritten)
            barrier(CLK_LOCAL_MEM_FENCE);
            if( (get_local_id(0) & (WARP_SIZE - 1)) == (WARP_SIZE - 1) )
                l_Data[get_local_id(0) >> LOG2_WARP_SIZE] = warpResult;

            //wait for warp scans to complete
            barrier(CLK_LOCAL_MEM_FENCE);
            if( get_local_id(0) < (WORKGROUP_SIZE / WARP_SIZE) ){
                //grab top warp elements
                uint val = l_Data[get_local_id(0)];
                //calculate exclsive scan and write back to shared memory
                l_Data[get_local_id(0)] = warpScanExclusive(val, l_Data, size >> LOG2_WARP_SIZE);
            }

            //return updated warp scans with exclusive scan results
            barrier(CLK_LOCAL_MEM_FENCE);
            return warpResult + l_Data[get_local_id(0) >> LOG2_WARP_SIZE];
        }else{
            return warpScanInclusive(idata, l_Data, size);
        }
    }

    inline uint scan1Exclusive(uint idata, __local uint *l_Data, uint size){
        return scan1Inclusive(idata, l_Data, size) - idata;
    }
#endif


//Vector scan: the array to be scanned is stored
//in work-item private memory as uint4
inline uint4 scan4Inclusive(uint4 data4, __local uint *l_Data, uint size){
    //Level-0 inclusive scan
    data4.y += data4.x;
    data4.z += data4.y;
    data4.w += data4.z;

    //Level-1 exclusive scan
    uint val = scan1Inclusive(data4.w, l_Data, size / 4) - data4.w;

    return (data4 + (uint4)val);
}

inline uint4 scan4Exclusive(uint4 data4, __local uint *l_Data, uint size){
    return scan4Inclusive(data4, l_Data, size) - data4;
}

////////////////////////////////////////////////////////////////////////////////
// Scan kernels
////////////////////////////////////////////////////////////////////////////////
__kernel __attribute__((reqd_work_group_size(WORKGROUP_SIZE, 1, 1)))
void scanExclusiveLocal1(
    __global uint4 *d_Dst,
    __global uint4 *d_Src,
    __local uint *l_Data,
    uint size
){
    //Load data
    uint4 idata4 = d_Src[get_global_id(0)];

    //Calculate exclusive scan
    uint4 odata4  = scan4Exclusive(idata4, l_Data, size);

    //Write back
    d_Dst[get_global_id(0)] = odata4;
}

//Exclusive scan of top elements of bottom-level scans (4 * THREADBLOCK_SIZE)
__kernel __attribute__((reqd_work_group_size(WORKGROUP_SIZE, 1, 1)))
void scanExclusiveLocal2(
    __global uint *d_Buf,
    __global uint *d_Dst,
    __global uint *d_Src,
    __local uint *l_Data,
    uint N,
    uint arrayLength
){
    //Load top elements
    //Convert results of bottom-level scan back to inclusive
    //Skip loads and stores for inactive work-items of the work-group with highest index(pos >= N)
    uint data = 0;
    if(get_global_id(0) < N)
    data =
        d_Dst[(4 * WORKGROUP_SIZE - 1) + (4 * WORKGROUP_SIZE) * get_global_id(0)] + 
        d_Src[(4 * WORKGROUP_SIZE - 1) + (4 * WORKGROUP_SIZE) * get_global_id(0)];

    //Compute
    uint odata = scan1Exclusive(data, l_Data, arrayLength);

    //Avoid out-of-bound access
    if(get_global_id(0) < N)
        d_Buf[get_global_id(0)] = odata;
}

//Final step of large-array scan: combine basic inclusive scan with exclusive scan of top elements of input arrays
__kernel __attribute__((reqd_work_group_size(WORKGROUP_SIZE, 1, 1)))
void uniformUpdate(
    __global uint4 *d_Data,
    __global uint *d_Buf
){
    __local uint buf[1];

    uint4 data4 = d_Data[get_global_id(0)];

    if(get_local_id(0) == 0)
        buf[0] = d_Buf[get_group_id(0)];

    barrier(CLK_LOCAL_MEM_FENCE);
    data4 += (uint4)buf[0];
    d_Data[get_global_id(0)] = data4;
}



