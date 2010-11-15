#pragma OPENCL EXTENSION cl_khr_fp64: enable


inline double functionDistance(double4 a, double4 b) {
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

    switch(ThisDimension) {
      case 2: d_index = (ThisCoord - MinPoint[0].z) * InvCellSize[ThisDimension]; break;
      case 1: d_index = (ThisCoord - MinPoint[0].y) * InvCellSize[ThisDimension]; break;
      case 0: d_index = (ThisCoord - MinPoint[0].x) * InvCellSize[ThisDimension]; break;
    }
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

double4 calculateCell(double4 ThisPoint, 
		      double Radius, 
		      __global double * N, 
		      __global double4 * MinPoint, 
		      __global double * InvCellSize )
{
    double4 Cell;
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
			     __global int * w_size,
			     __global int * maxResults
			    ) 
{
    int index;
    int tope;
    double4 cellBegin;
    double4 cellEnd;

    int ip = get_global_id(0);

    if(ip < w_size[0]) {
	cellBegin = calculateCell(point[ip],-radius[0],N,MinPoint,InvCellSize);
	cellEnd   = calculateCell(point[ip], radius[0],N,MinPoint,InvCellSize);
	results[ip] = 0;

	for(size_t i = cellBegin.x; i <= cellEnd.x; i++) 
	{
	    for(size_t j = cellBegin.y; j <= cellEnd.y; j++) 
	    {
		for(size_t k = cellBegin.z; k <= cellEnd.z; k++) 
		{
		    index = calculateIndexForCell((double4)(i,j,k,-1),N);
		    
		    for(size_t l = 0; (IndexCellReference[index]+l) < IndexCellReference[index+1]; l++) 
		    {
			if(functionDistance(point[ip],BinsContainer[(int)IndexCellReference[index]+l]) < radius2[0]) 
			{
			    if(results[ip] <= maxResults[0])
				outdata[ip*maxResults[0]+results[ip]] = BinsContainer[(int)IndexCellReference[index]+l].w;
			    results[ip]++;
			}
		    }
		}
	    }
	}
    }
}