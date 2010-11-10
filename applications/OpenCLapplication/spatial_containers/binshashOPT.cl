#pragma OPENCL EXTENSION cl_khr_fp64: enable


inline double functionDistance(double4 a, double4 b) {
    double4 temp;
    temp = (a-b) * (a-b);
    return temp.x + temp.y + temp.z;
}

int calculatePosition(double ThisCoord, 
		      int ThisDimension, 
		      __global double * mN, 
		      __global double4 * mMinPoint, 
		      __global double * mInvCellSize 
		      )
{
    int d_index;

    switch(ThisDimension) {
      case 2: d_index = (ThisCoord - mMinPoint[0].z) * mInvCellSize[ThisDimension]; break;
      case 1: d_index = (ThisCoord - mMinPoint[0].y) * mInvCellSize[ThisDimension]; break;
      case 0: d_index = (ThisCoord - mMinPoint[0].x) * mInvCellSize[ThisDimension]; break;
    }
    int index = (int)( (d_index < 0.00) ? 0.00 : d_index );
    
    return  (index > mN[ThisDimension]-1) ? mN[ThisDimension]-1 : index;
}

int calculateIndex(double4 ThisPoint, 
		    __global double * mN, 
		    __global double4 * mMinPoint, 
		    __global double * mInvCellSize
		  ) 
{
    int Index = 0;
    Index += calculatePosition(ThisPoint.z,2,mN,mMinPoint,mInvCellSize);
    Index *= mN[1];
    Index += calculatePosition(ThisPoint.y,1,mN,mMinPoint,mInvCellSize);
    Index *= mN[0];
    Index += calculatePosition(ThisPoint.x,0,mN,mMinPoint,mInvCellSize);
    return Index;
}

int calculateIndexForCell( double4 ThisIndex, 
			   __global double * mN 
			  )
{
    int Index = 0;
    Index += ThisIndex.z;
    Index *= mN[1];
    Index += ThisIndex.y;
    Index *= mN[0];
    Index += ThisIndex.x;
    return Index;
}

double4 calculateCell(double4 ThisPoint, 
		      double Radius, 
		      __global double * mN, 
		      __global double4 * mMinPoint, 
		      __global double * mInvCellSize )
{
    double4 Cell;
    Cell.x = calculatePosition(ThisPoint.x+Radius,0,mN,mMinPoint,mInvCellSize);
    Cell.y = calculatePosition(ThisPoint.y+Radius,1,mN,mMinPoint,mInvCellSize);
    Cell.z = calculatePosition(ThisPoint.z+Radius,2,mN,mMinPoint,mInvCellSize);
    return Cell;
}

__kernel void GenerateBins(__global double4 * points, 
			   __global double4 * cell,
			   __global int * IndexCell,
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
	IndexCell[i] = 0;
	IndexCellReference[i] = 0;
    }
  
    //Update storage counter, storing ahead
    for(int i = 0; i < POINT_SIZE; i++) 
    {
	size_t index = calculateIndex(points[i],N,MinPoint,InvCellSize);
	IndexCell[index+1]++; 
    }

    IndexCellReference[0] = 0;
    for(int i = 1; i < CELL_SIZE; i++) 
    {
	IndexCell[i] += IndexCell[i-1];
	IndexCellReference[i] = IndexCell[i];
    }

    //Storing in bins
    size_t k = 0;
    size_t pos,cellIndex;

    while(k < POINT_SIZE) 
    {
	cellIndex = calculateIndex(points[k],N,MinPoint,InvCellSize);
	pos = IndexCell[cellIndex]++;
	double4 point = points[k];
	BinsContainer[pos] = point;
	k++;
    }

    for(int i = 0; i < CELL_SIZE; i++)
    {
	cell[i] = BinsContainer[IndexCellReference[i]];
    }
}

}

__kernel void SearchInRadiusMultiple(__global int * IndexCellReference,
			     __global double4 * BinsContainer,
			     __global double * mInvCellSize,
			     __global double * mN,
			     __global double * radius,
			     __global double * radius2,
			     __global double4 * point,
			     __global double4 * mMinPoint,
			     __global double4 * outdata,
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
	cellBegin = calculateCell(point[ip],-radius[0],mN,mMinPoint,mInvCellSize);
	cellEnd   = calculateCell(point[ip], radius[0],mN,mMinPoint,mInvCellSize);
	results[ip] = 0;

	for(size_t i = cellBegin.x; i <= cellEnd.x; i++) {
	    for(size_t j = cellBegin.y; j <= cellEnd.y; j++) {
		for(size_t k = cellBegin.z; k <= cellEnd.z; k++) {
		  
		    index = calculateIndexForCell((double4)(i,j,k,-1),mN);
		    
		    for(int l = 0; (l + IndexCellReference[index]) < IndexCellReference[index+1]; l++) 
		    {
			if(functionDistance(point[ip],BinsContainer[(int)IndexCellReference[index]+l]) < radius2[0]) 
			{
			    if(results[ip] < maxResults[0])
				outdata[ip*maxResults[0]+results[ip]] = BinsContainer[(int)IndexCellReference[index]+l];
			    results[ip]++;
			}
		    }
		}
	    }
	}
    }
}