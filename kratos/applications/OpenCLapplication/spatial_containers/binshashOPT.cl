#pragma OPENCL EXTENSION cl_amd_fp64: enable
//#pragma OPENCL EXTENSION cl_khr_fp64: enable

#ifndef cl_amd_fp64

	// Failed, probably we are not on an ATI platform, so try Khronos version
	#pragma OPENCL EXTENSION cl_khr_fp64: enable

#endif


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
    int4 Cell;

    Cell.x = calculatePosition(ThisPoint.x+Radius,0,N,MinPoint,InvCellSize);
    Cell.y = calculatePosition(ThisPoint.y+Radius,1,N,MinPoint,InvCellSize);
    Cell.z = calculatePosition(ThisPoint.z+Radius,2,N,MinPoint,InvCellSize);

    return Cell;
}

__kernel void GenerateBins(__global double4 * points, 
			   __global int * IndexCellReference,
			   __global double * InvCellSize,
			   __global double * N,
			   __global double4 * MinPoint,
			   __global double4 * BinsContainer
			  ) 
{
    //GenerateBins

    int t = get_global_id(0);

    if (t < 1) {
	//Init vectors
	for(int i = 0; i < CELL_SIZE; i++)
	{
	    IndexCellReference[i] = 0;
	}
      
	//Update storage counter, storing ahead
	for(int i = 0; i < POINT_SIZE; i++) 
	{
	    size_t index = calculateIndex(points[i],N,MinPoint,InvCellSize);
	    IndexCellReference[index+1]++; 
	}

	for(int i = 1; i < CELL_SIZE; i++) 
	{
	    IndexCellReference[i] += IndexCellReference[i-1];
	}

	//Storing in bins
	size_t k = 0;
	size_t pos,cellIndex;

	while(k < POINT_SIZE) 
	{
	    cellIndex = calculateIndex(points[k],N,MinPoint,InvCellSize);
	    pos = IndexCellReference[cellIndex]++;
	    double4 point = points[k];
	    BinsContainer[pos] = point;
	    k++;
	}

	for(int i = CELL_SIZE-1; i > 0; i--) 
	{
	    IndexCellReference[i] = IndexCellReference[i-1];
	}
	
	IndexCellReference[0] = 0;

    }

}

__kernel void SearchInRadiusMultiple(__global int * IndexCellReference,
			     __global double4 * BinsContainer,
			     __global double * InvCellSize,
			     __global double * N,
			     __global double * radius,
			     __global double * radius2,
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

    for(size_t i = cellBegin.x; i <= cellEnd.x; i++) 
    {
	for(size_t j = cellBegin.y; j <= cellEnd.y; j++) 
	{
	    for(size_t k = cellBegin.z; k <= cellEnd.z; k++) 
	    {
		__private int index = calculateIndexForCell((double4)(i,j,k,-1),N);
		
		__private int loIndex = IndexCellReference[index];
		__private int hiIndex = IndexCellReference[index+1];
		
		for(size_t l = loIndex; l < hiIndex; l++) 
		{
		    if(isless(functionDistance(pointIp,BinsContainer[l]),radius2[0])) 
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

__kernel void SearchNearestMultiple(__global int * IndexCellReference,
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
    distance[ip] = INFINITY;

    for(size_t i = cellBegin.x; i <= cellEnd.x; i++) 
    {
	for(size_t j = cellBegin.y; j <= cellEnd.y; j++) 
	{
	    for(size_t k = cellBegin.z; k <= cellEnd.z; k++) 
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
		    }
		}
	    }
	}
    }
}

//Not sure if this code is correct
__kernel void SearchNearestMultipleSpiral(__global int * IndexCellReference,
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

      while( !Found && k <= kend )
      {
	int jini = cellPoint.y - cellcuberadius;
	int jend = cellPoint.y + cellcuberadius;
	int j = jini;
	while( !Found && j <= jend )
	{
          int iini = cellPoint.x - cellcuberadius;
	  int iend = cellPoint.x + cellcuberadius;
	  int i = iini;
	  while( !Found && i <= iend )
	  {
	    if( (
		  (k == kini) || 
		  (k == kend) || 
		  (j == 0)    || 
		  (j == jend) || 
		  (i == 0)    || 
		  (i == iend)               ) &&
		( k >= cellBegin.z && k <= cellEnd.z ) &&
		( j >= cellBegin.y && j <= cellEnd.y ) &&
		( i >= cellBegin.x && i <= cellEnd.z ) )
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
	if(!(k < kend) && partialFound)
	  Found = true; //congratulations! you found it!!!
      } // endl while k
      cellcuberadius++;
    } //endl while cuberadius
}
