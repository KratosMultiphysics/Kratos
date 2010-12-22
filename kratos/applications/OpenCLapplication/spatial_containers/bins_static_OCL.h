//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-27 17:02:19 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_BINS_OCL_CONTAINER_H_INCLUDE)
#define KRATOS_BINS_OCL_CONTAINER_H_INCLUDE

#define WORKGROUP_SIZE 2

#include "../../../kratos/spatial_containers/tree.h"
#include "opencl_interface.h"
#include <malloc.h>
#include <stdio.h>
#include <string.h>

namespace Kratos {



template<  std::size_t TDimension,
           class TPointType,
		   class TContainerType,
		   class TPointerType = typename TContainerType::value_type,
		   class TIteratorType = typename TContainerType::iterator,
		   class TDistanceIteratorType = typename std::vector<double>::iterator,
		   class TDistanceFunction = Kratos::SearchUtils::SquaredDistanceFunction<TDimension,TPointType> >
class BinsOCL : public TreeNode<TDimension,TPointType, TPointerType, TIteratorType, TDistanceIteratorType> 
{


   public:

      /// Pointer definition of Bins
      KRATOS_CLASS_POINTER_DEFINITION(BinsOCL);

      typedef TPointType                         PointType;
      typedef TContainerType                     ContainerType;
      typedef TIteratorType                      IteratorType;
      typedef TDistanceIteratorType              DistanceIteratorType;
      typedef TPointerType                       PointerType;
      typedef TDistanceFunction                  DistanceFunction;
      enum { Dimension = TDimension };

      typedef TreeNode<Dimension,PointType,PointerType,IteratorType,DistanceIteratorType> TreeNodeType;

      typedef typename TreeNodeType::SizeType        SizeType;
      typedef typename TreeNodeType::IndexType       IndexType;
      typedef typename TreeNodeType::CoordinateType  CoordinateType;

      typedef TreeNodeType LeafType;

      // not always PointVector == ContainerType ( if ContainerType = C array )
      typedef std::vector<PointerType>        PointVector;
      typedef typename PointVector::iterator  PointIterator;

      typedef std::vector<IteratorType>          IteratorVector;
      typedef typename IteratorVector::iterator  IteratorIterator;
      typedef typename IteratorVector::const_iterator IteratorConstIterator;

      typedef Tvector<IndexType,TDimension>   CellType;
    
      typedef typename TreeNodeType::IteratorIteratorType IteratorIteratorType;
      typedef typename TreeNodeType::SearchStructureType SearchStructureType;

      typedef Kratos::SearchUtils::SearchNearestInRange<PointType,PointerType,IteratorType,DistanceFunction,CoordinateType> SearchNearestInRange;
      typedef Kratos::SearchUtils::SearchRadiusInRange<PointType,IteratorType,DistanceIteratorType,DistanceFunction,SizeType,CoordinateType> SearchRadiusInRange;
      typedef Kratos::SearchUtils::SearchBoxInRange<PointType,IteratorType,SizeType,TDimension> SearchBoxInRange;
      
    public:

	//************************************************************************

	//Constructor
	BinsOCL() : mPointBegin(this->NullIterator()), mPointEnd(this->NullIterator()) {};

	//************************************************************************

	BinsOCL( IteratorType const& PointBegin, IteratorType const& PointEnd, SizeType BucketSize = 1 )
	    : mPointBegin(PointBegin), mPointEnd(PointEnd)
	{
	    if(mPointBegin==mPointEnd)
		return;
	    
	    mPointSize = SearchUtils::PointerDistance(mPointBegin, mPointEnd);
	    mTriangleSize = 0;
	    mProblemSize = mPointSize;
	    pointsToSearch = new cl_double4[mProblemSize];
	    
	    CalculateBoundingBox();
	    CalculateCellSize();
	    AllocateCellsContainer();
	    InitOCL();
	}

	//************************************************************************
	 
	BinsOCL( IteratorType const& PointBegin, IteratorType const& PointEnd, cl_int4 * const& indexArray, int TSize, SizeType BucketSize = 1 )
	    : mPointBegin(PointBegin), mPointEnd(PointEnd)
	{
	    if(mPointBegin==mPointEnd)
		return;
	    
	    mPointSize = SearchUtils::PointerDistance(mPointBegin, mPointEnd);
	    mTriangleSize = TSize;
	    mProblemSize = mTriangleSize;
	    pointsToSearch = new cl_double4[mProblemSize];
	    
	    CalculateBoundingBox(indexArray);
	    CalculateCellSize();
	    AllocateCellsContainer();
	    InitOCL();
	    
	}

	//************************************************************************

	// destructor
	virtual ~BinsOCL(){ }

	//************************************************************************

	IteratorType Begin() { return mPointBegin; }

	//************************************************************************

	IteratorType End() { return mPointBegin; }

	//************************************************************************

	CoordinateType CellSize( SizeType const& iDim ) { return mCellSize[iDim]; }

	//************************************************************************

	SizeType NumCell( SizeType const& iDim ) { return mN[iDim]; }

	//************************************************************************
	
	void GenerateSampleInput(int newSize)
	{
	    free(pointsToSearch);
	    mProblemSize = newSize * newSize;
	    pointsToSearch = new cl_double4[mProblemSize];
	    double InvNewSize = 1/(double)newSize;
	    
	    int index = 0;
	    for(int i = 0; i < newSize; i++)
	    {
		for(int j = 0; j < newSize; j++)
		{
		    //for(int k = 0; k < mProblemSize; k++)
		    //{
			pointsToSearch[index].x = (i+1) * InvNewSize;
			pointsToSearch[index].y = (j+1) * InvNewSize;
			pointsToSearch[index].z = 0;//(k+1)*InvProblemSize;
			pointsToSearch[index].w = index++;//(k+1)*InvProblemSize;
		    //}
		}
	    }
	    
	}
	
	//************************************************************************
	 
    private:

	//************************************************************************
     
	void InitOCL() 
	{
	   
	   //"Advanced Micro Devices, Inc."
	   //"NVIDIA Corporation"
	   
	    OCLDeviceGroup = new Kratos::OpenCL::DeviceGroup(CL_DEVICE_TYPE_ALL,true,"Advanced Micro Devices, Inc.");
	    std::cout << "Found " << OCLDeviceGroup->DeviceNo << " device(s)." << std::endl;
	    for (cl_uint i = 0; i < OCLDeviceGroup->DeviceNo; i++)
	    {
		std::cout << "  Device " << i << ": " << Kratos::OpenCL::DeviceTypeString(OCLDeviceGroup->DeviceTypes[i]) << std::endl;
	    }
	    
	    char params[512] = "";
	    sprintf(params,"-D POINT_SIZE=%d -D TDIMENSION=%Zu -D CELL_SIZE=%d -D WORKGROUP_SIZE=%d",mProblemSize,TDimension,mCellSizeInt+1,WORKGROUP_SIZE);
	      
	    OCL_program = OCLDeviceGroup->BuildProgramFromFile("binshashOPT.cl",params);  
	}
	 
	//************************************************************************
	 
	//Calculate the barycenter of the triangle/tetrahedre formed by index in index4
	PointType BaryCenter(cl_int4 const& index4)
	{
	    IteratorType Point = mPointBegin;
	    PointType a,b,c,d,r;
	    a = (**(Point = mPointBegin + index4.x-1));
	    b = (**(Point = mPointBegin + index4.y-1));
	    c = (**(Point = mPointBegin + index4.z-1));
	    
	    //2D
	    if(index4.w-1 >= 0)
	    {
		d = (**(Point = mPointBegin + index4.w-1));
		
		for(int i = 0; i < TDimension; i++)
		{
		    r[i] = (a[i] + b[i] + c[i] + d[i]) / 4; // * 0.25f ? 
		}
	    }
	    //3D
	    else
	    {
		//std::cout << "Iteration: " << (**(Point = mPointBegin + index4.x-1)) << "" << a.id << " " << b.id << " " << c.id << std::endl;
		for(int i = 0; i < TDimension; i++)
		{
		    r[i] = (a[i] + b[i] + c[i]) / 3; // * 0.33333333333333334f ? 
		}
	    }
	    return r;
	}

	//Calculates the bounding box for standar input of points.
	void CalculateBoundingBox() 
	{
	    for(SizeType i = 0 ; i < TDimension ; i++)
	    {
		mMinPoint[i] = (**mPointBegin)[i];
		mMaxPoint[i] = (**mPointBegin)[i];
	    }
	   
	    for(IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
	    {
		for(SizeType i = 0 ; i < TDimension ; i++)
		{
		    if( (**Point)[i] < mMinPoint[i] ) mMinPoint[i] = (**Point)[i];
		    if( (**Point)[i] > mMaxPoint[i] ) mMaxPoint[i] = (**Point)[i];
		}
	    }

	    // Check que se ha echo todo bien
	    for(SizeType i = 0 ; i < TDimension ; i++)
		std::cout << mMinPoint[i] << "  ";
	    std::cout << std::endl;
	}
	 	
	//Calculates the bounding box for barycenters of tringle/tetrahedre input
	void CalculateBoundingBox(cl_int4 * const& indexArray) 
	{
	    mBaryCenters = (PointType *)malloc(sizeof(PointType) * mTriangleSize); //new PointType[mTriangleSize];
	    mTriangles = (cl_int4 *)malloc(sizeof(cl_int4) * mTriangleSize); 
	    
	    PointType barycenter = BaryCenter(indexArray[0]);
	    
	    for(SizeType i = 0 ; i < TDimension ; i++)
	    {
		mMinPoint[i] = barycenter[i];
		mMaxPoint[i] = barycenter[i];
	    }
	    
	    mBaryCenters[0] = barycenter;
	    mBaryCenters[0].id = 0;
	    
// 	    pointsToSearch[0].x = barycenter[0];
// 	    pointsToSearch[0].y = barycenter[1];
// 	    pointsToSearch[0].z = barycenter[2];
	    
	    mTriangles[0] = indexArray[0];
	   
	    for(int j = 1; j < mProblemSize; j++)
	    {
		barycenter = BaryCenter(indexArray[j]);
		barycenter.id = j;
		
		for(SizeType i = 0 ; i < TDimension ; i++){
		    if( barycenter[i] < mMinPoint[i] ) mMinPoint[i] = barycenter[i];
		    if( barycenter[i] > mMaxPoint[i] ) mMaxPoint[i] = barycenter[i];
		}
		//std::cout << "BarycenterPoint: " << barycenter << std::endl;
		mBaryCenters[j] = barycenter;
		mBaryCenters[j].id = j;
		
// 		pointsToSearch[j].x = barycenter[0];
// 		pointsToSearch[j].y = barycenter[1];
// 		pointsToSearch[j].z = barycenter[2];
		
		mTriangles[j] = indexArray[j];
	    }
	   
	    // Check que se ha echo todo bien
	    for(SizeType i = 0 ; i < TDimension ; i++)
		std::cout << mMinPoint[i] << "  ";
	    std::cout << std::endl;
	}

	//************************************************************************

	void CalculateCellSize() 
	{ 
	    CoordinateType delta[TDimension];
	    CoordinateType alpha[TDimension];
	    CoordinateType mult_delta = 1.00;
	    SizeType index = 0;
	  
	    for(SizeType i = 0 ; i < TDimension ; i++) 
	    {
		delta[i] = mMaxPoint[i] - mMinPoint[i];
		if ( delta[i] > delta[index] )
		    index = i;
		delta[i] = (delta[i] == 0.00) ? 1.00 : delta[i];
	    }

	    for(SizeType i = 0 ; i < TDimension ; i++){
		alpha[i] = delta[i] / delta[index];
		mult_delta *= alpha[i];
	    }

	    mN[index] = static_cast<SizeType>( pow(static_cast<CoordinateType>(mProblemSize/mult_delta), 1.00/TDimension)+1 );
	   
	    for(SizeType i = 0 ; i < TDimension ; i++)
	    {
		if(i!=index) 
		{
		    mN[i] = static_cast<SizeType>(alpha[i] * mN[index]);
		    mN[i] = ( mN[i] == 0 ) ? 1 : mN[i];
		}
	    }

	    for(SizeType i = 0 ; i < TDimension ; i++)
	    {
		mCellSize[i] = delta[i] / mN[i];
		mInvCellSize[i] = 1.00 / mCellSize[i];
	    }
	 
	}

	//************************************************************************

	void CalculateCellSize( CoordinateType BoxSize ) 
	{
	    for(SizeType i = 0 ; i < TDimension ; i++)
	    {
		mCellSize[i] = BoxSize;
		mInvCellSize[i] = 1.00 / mCellSize[i];
		mN[i] = static_cast<SizeType>( (mMaxPoint[i]-mMinPoint[i]) / mCellSize[i]) + 1;
	    }
	}

	//************************************************************************

	void AllocateCellsContainer() 
	{
	    mCellSizeInt = 1;
	    for(SizeType i = 0 ; i < TDimension ; i++)
		  mCellSizeInt *= mN[i];
	    std::cout << "Size: \t" << mCellSizeInt << std::endl;
	    mIndexCell.resize(mCellSizeInt+1);
	    mIndexCellBegin = mIndexCell.begin();
	    mIndexCellEnd   = mIndexCell.end();
	}

	//************************************************************************

	void GenerateBins()
	{
	    struct timespec begin;
	    struct timespec end;
	    
	    clock_gettime( CLOCK_REALTIME, &begin );
	    clock_gettime( CLOCK_REALTIME, &end );
	   
    	    clock_gettime( CLOCK_REALTIME, &begin );
	    OCLGenerateBinsA      = OCLDeviceGroup->RegisterKernel(OCL_program,"GenerateBinsA");
	    OCLGenerateBinsB1     = OCLDeviceGroup->RegisterKernel(OCL_program,"scanExclusiveLocal1");
	    OCLGenerateBinsB2     = OCLDeviceGroup->RegisterKernel(OCL_program,"scanExclusiveLocal2");
	    OCLGenerateBinsB3  	  = OCLDeviceGroup->RegisterKernel(OCL_program,"uniformUpdate");
	    OCLGenerateBinsC      = OCLDeviceGroup->RegisterKernel(OCL_program,"GenerateBinsC");
	    OCLSearchInRadius     = OCLDeviceGroup->RegisterKernel(OCL_program,"SearchInRadiusMultiple");
	    OCLSearchTriangles    = OCLDeviceGroup->RegisterKernel(OCL_program,"SearchTriangle");
	    OCLSearchNearest      = OCLDeviceGroup->RegisterKernel(OCL_program,"SearchNearestMultiple");
	    OCLSearchNearestCubic = OCLDeviceGroup->RegisterKernel(OCL_program,"SearchNearestMultipleCubic");
	    clock_gettime( CLOCK_REALTIME, &end );
	    
	    std::cout << "Kernel Reg:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
	    
	    int partialReduction = 0;
	    int Offset = 0;
	    int PowSize = 1;
	    int pr = 0;
	    int wiz = 512; // TODO: change this, now is the max number of work items in my card!
	    
	    while(PowSize < mCellSizeInt+1)
		PowSize <<= 1;
	    
	    std::cout << "POW:\t\t" << PowSize << std::endl;
	    
	    int * IndexCellReference   = (int *)malloc( (PowSize) * sizeof(int) );
	    int * IndexCellReferenceO  = (int *)malloc( (PowSize) * sizeof(int) );
	    cl_double4 * BinsContainer = (cl_double4 *)malloc(mProblemSize * sizeof(cl_double4));
	    cl_double4 * PointsBins    = (cl_double4 *)malloc(mProblemSize * sizeof(cl_double4));
	    cl_double4 * MinPoint      = new cl_double4();
	    
	    double InvCellSize[TDimension];
	    double N[TDimension];
	        	    
	    clock_gettime( CLOCK_REALTIME, &begin );
	    
	    if(mTriangleSize)
	    {
		for(int i = 0; i < mTriangleSize; i++)
		{
		  //std::cout << "MODE T  " <<  mProblemPoints[i].id << std::endl;
		    PointsBins[i].x = mBaryCenters[i][0];
		    PointsBins[i].y = mBaryCenters[i][1];
		    PointsBins[i].z = mBaryCenters[i][2];
		    PointsBins[i].w = mBaryCenters[i].id;
		}
	    }
	    else
	    {
		for( IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
		{
		    static int k = 0;
		    PointsBins[k].x = (**Point)[0];
		    PointsBins[k].y = (**Point)[1];
		    PointsBins[k].z = (**Point)[2];
		    PointsBins[k].w = (**Point).id;
		    k++;
		}
	    }

	    MinPoint->x = mMinPoint[0];
	    MinPoint->y = mMinPoint[1];
	    MinPoint->z = mMinPoint[2];
	    
	    for(int i = 0; i < PowSize; i++)
		IndexCellReference[i] = 0;
	    
	    for(int i = 0; i < TDimension; i++) 
	    {
		InvCellSize[i] = mInvCellSize[i];
		N[i] = mN[i];
	    }
	    clock_gettime( CLOCK_REALTIME, &end );
	    
	    std::cout << "Var filling:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
	    
	    OCL_PointsBins         = OCLDeviceGroup->CreateBuffer(sizeof(cl_double4) * mProblemSize, CL_MEM_READ_ONLY);   
	    OCL_IndexCellReference = OCLDeviceGroup->CreateBuffer(sizeof(int) * PowSize, CL_MEM_READ_WRITE);
	    OCL_IndexCellReferenceO= OCLDeviceGroup->CreateBuffer(sizeof(int) * PowSize, CL_MEM_READ_WRITE);
	    OCL_IndexCellReferenceU= OCLDeviceGroup->CreateBuffer(sizeof(int) * PowSize, CL_MEM_READ_WRITE);
	    OCL_BinsContainer      = OCLDeviceGroup->CreateBuffer(sizeof(cl_double4) * mProblemSize, CL_MEM_READ_WRITE);
	    OCL_InvCellSize        = OCLDeviceGroup->CreateBuffer(sizeof(double) * TDimension, CL_MEM_READ_ONLY);
	    OCL_N                  = OCLDeviceGroup->CreateBuffer(sizeof(double) * TDimension, CL_MEM_READ_ONLY);
	    OCL_MinPoint           = OCLDeviceGroup->CreateBuffer(sizeof(cl_double4), CL_MEM_READ_ONLY);
	    OCL_PowSize            = OCLDeviceGroup->CreateBuffer(sizeof(int), CL_MEM_READ_ONLY);
	    OCL_Offset             = OCLDeviceGroup->CreateBuffer(sizeof(int), CL_MEM_READ_ONLY);
	    OCL_PartialReduction   = OCLDeviceGroup->CreateBuffer(sizeof(int), CL_MEM_READ_WRITE);
	    
	    OCLDeviceGroup->CopyBuffer(OCL_PointsBins        , OpenCL::HostToDevice, OpenCL::VoidPList(1,PointsBins));
	    OCLDeviceGroup->CopyBuffer(OCL_IndexCellReference, OpenCL::HostToDevice, OpenCL::VoidPList(1,IndexCellReference));
	    OCLDeviceGroup->CopyBuffer(OCL_InvCellSize	     , OpenCL::HostToDevice, OpenCL::VoidPList(1,InvCellSize));
	    OCLDeviceGroup->CopyBuffer(OCL_N                 , OpenCL::HostToDevice, OpenCL::VoidPList(1,N));
	    OCLDeviceGroup->CopyBuffer(OCL_MinPoint          , OpenCL::HostToDevice, OpenCL::VoidPList(1,MinPoint));
	    OCLDeviceGroup->CopyBuffer(OCL_PowSize           , OpenCL::HostToDevice, OpenCL::VoidPList(1,&wiz));
	    OCLDeviceGroup->CopyBuffer(OCL_Offset            , OpenCL::HostToDevice, OpenCL::VoidPList(1,&Offset));
	    OCLDeviceGroup->CopyBuffer(OCL_PartialReduction  , OpenCL::HostToDevice, OpenCL::VoidPList(1,&pr));
	
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsA, 0,  OCL_PointsBins);
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsA, 1,  OCL_IndexCellReference);
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsA, 2,  OCL_InvCellSize);
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsA, 3,  OCL_N);
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsA, 4,  OCL_MinPoint);
	    
	    clock_gettime( CLOCK_REALTIME, &begin );
	    OCLDeviceGroup->ExecuteKernel(OCLGenerateBinsA, mTriangleSize);
	    clock_gettime( CLOCK_REALTIME, &end );
	    
	    std::cout << "KERNEL A EXECUTED:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
	    
	    OCLDeviceGroup->CopyBuffer(OCL_IndexCellReference, OpenCL::DeviceToHost, OpenCL::VoidPList(1,IndexCellReference));
	    /*
	    for(int i = 0; i < PowSize; i++) 
	    {
	      std::cout << "IndexCellRe " << i << "  " << IndexCellReference[i] << std::endl;
	    }
	    std::cout << "------------------" << std::endl;
	    */
	    
	     //SERIAL VERSION
	    OCLDeviceGroup->CopyBuffer(OCL_IndexCellReference, OpenCL::DeviceToHost, OpenCL::VoidPList(1,IndexCellReferenceO));
	    for(int i = 1; i < mCellSizeInt+1; i++)
	    {
	        IndexCellReferenceO[i]+=IndexCellReferenceO[i-1];
	    }
	    OCLDeviceGroup->CopyBuffer(OCL_IndexCellReferenceO, OpenCL::HostToDevice, OpenCL::VoidPList(1,IndexCellReferenceO));
	    
	    //PARALLEL VERSION
	    
// 	    int aux = 4 * WORKGROUP_SIZE;
// 	    OCLDeviceGroup->CopyBuffer(OCL_PowSize, OpenCL::HostToDevice, OpenCL::VoidPList(1,&aux));
// 	        	    
// 	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsB1, 0, OCL_IndexCellReferenceO);
// 	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsB1, 1, OCL_IndexCellReference);
// 	    OCLDeviceGroup->SetLocalMemAsKernelArg(OCLGenerateBinsB1, 2, 2 * WORKGROUP_SIZE * sizeof(int));
// 	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsB1, 3, OCL_PowSize);
// 	    
// 	    clock_gettime( CLOCK_REALTIME, &begin );
// 	    
// 	    OCLDeviceGroup->ExecuteKernel(OCLGenerateBinsB1, PowSize / 4);
// 	    std::cout << "KERNEL B1 EXECUTED:\t\t" << std::endl;
// 	     
// 	    OCLDeviceGroup->CopyBuffer(OCL_IndexCellReferenceO, OpenCL::DeviceToHost, OpenCL::VoidPList(1,IndexCellReferenceO));
// // 	    for(int i = 0; i < PowSize; i++) {
// // 	      std::cout << "IndexCellReB1 " << i << "  " << IndexCellReferenceO[i] << std::endl;
// // 	    } 
// // 	    std::cout << "------------------" << std::endl;
// 	    
// 	    aux = (8 * WORKGROUP_SIZE) / (4 * WORKGROUP_SIZE);
// 	    OCLDeviceGroup->CopyBuffer(OCL_PowSize, OpenCL::HostToDevice, OpenCL::VoidPList(1,&aux));
// 	    
// 	    aux = (PowSize / (8 * WORKGROUP_SIZE)) * ((8 * WORKGROUP_SIZE) / (4 * WORKGROUP_SIZE));
// 	    OCLDeviceGroup->CopyBuffer(OCL_Offset, OpenCL::HostToDevice, OpenCL::VoidPList(1,&aux));
// 	    
//             OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsB2, 0, OCL_IndexCellReferenceU);
//     	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsB2, 1, OCL_IndexCellReferenceO);
// 	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsB2, 2, OCL_IndexCellReference);
// 	    OCLDeviceGroup->SetLocalMemAsKernelArg(OCLGenerateBinsB2, 3, 2 * WORKGROUP_SIZE * sizeof(int));
// 	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsB2, 4, OCL_Offset);
// 	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsB2, 5, OCL_PowSize);
// 	    
// 	    aux = ((aux % WORKGROUP_SIZE) == 0) ? aux : (aux - aux % WORKGROUP_SIZE + WORKGROUP_SIZE);
// 	    
// 	    OCLDeviceGroup->ExecuteKernel(OCLGenerateBinsB2, aux);
// 	    std::cout << "KERNEL B2 EXECUTED:\t\t" << std::endl;
// 	    
// 	    OCLDeviceGroup->CopyBuffer(OCL_IndexCellReferenceO, OpenCL::DeviceToHost, OpenCL::VoidPList(1,IndexCellReferenceO));
// // 	    for(int i = 0; i < PowSize; i++) {
// // 	      std::cout << "IndexCellReB2 " << i << "  " << IndexCellReferenceO[i] << std::endl;
// // 	    } 
// // 	    std::cout << "------------------" << std::endl;
// 	    
// 	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsB3, 0, OCL_IndexCellReferenceO);
// 	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsB3, 1, OCL_IndexCellReferenceU);
// 	    
// 	    OCLDeviceGroup->ExecuteKernel(OCLGenerateBinsB3, (PowSize / (4 * WORKGROUP_SIZE)) * WORKGROUP_SIZE);
// 	    std::cout << "KERNEL B3 EXECUTED:\t\t" << std::endl;
// 	    
// 	    OCLDeviceGroup->CopyBuffer(OCL_IndexCellReferenceO, OpenCL::DeviceToHost, OpenCL::VoidPList(1,IndexCellReferenceO));
// // 	    for(int i = 0; i < PowSize; i++) {
// // 	      std::cout << "IndexCellReB3 " << i << "  " << IndexCellReferenceO[i] << std::endl;
// // 	    }
// // 	    std::cout << "------------------" << std::endl;
// 	    
// 	    clock_gettime( CLOCK_REALTIME, &end );
	    
	    
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsC, 0,  OCL_PointsBins);
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsC, 1,  OCL_IndexCellReferenceO);
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsC, 2,  OCL_InvCellSize);
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsC, 3,  OCL_N);
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsC, 4,  OCL_MinPoint);
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsC, 5,  OCL_BinsContainer);   
	    
	    clock_gettime( CLOCK_REALTIME, &begin );
	    OCLDeviceGroup->ExecuteKernel(OCLGenerateBinsC, mProblemSize);
	    clock_gettime( CLOCK_REALTIME, &end );
	    
	    std::cout << "KernelC Ex:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
	    
	    clock_gettime( CLOCK_REALTIME, &begin );
	    OCLDeviceGroup->CopyBuffer(OCL_BinsContainer, OpenCL::DeviceToHost, OpenCL::VoidPList(1,BinsContainer));
	    clock_gettime( CLOCK_REALTIME, &end );
	    
	    std::cout << "Kernel Copy:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
	    
	    OCLDeviceGroup->CopyBuffer(OCL_IndexCellReferenceO, OpenCL::DeviceToHost, OpenCL::VoidPList(1,IndexCellReferenceO));
	    /*
	    for(int i = 0; i < PowSize; i++) {
	      std::cout << "IndexCellReKE " << i << "  " << IndexCellReferenceO[i] << std::endl;
	    }
	    std::cout << "------------------" << std::endl;
	    
	    for(int i = 0; i < SearchUtils::PointerDistance(mPointBegin, mPointEnd); i++ ) 
	    {
	      std::cout << BinsContainer[i].w << std::endl;
	    }
	    */
	}

	//************************************************************************

	IndexType CalculatePosition( CoordinateType const& ThisCoord, SizeType ThisDimension )
	{
	    CoordinateType d_index = (ThisCoord - mMinPoint[ThisDimension]) * mInvCellSize[ThisDimension];
	    IndexType index = static_cast<SizeType>( (d_index < 0.00) ? 0.00 : d_index );
	    
	    return  (index > mN[ThisDimension]-1) ? mN[ThisDimension]-1 : index;
	}

	//************************************************************************

	IndexType CalculateIndex( PointType const& ThisPoint )
	{
	    IndexType Index = 0;
	    for(SizeType iDim = TDimension-1 ; iDim > 0 ; iDim--)
	    {
		Index += CalculatePosition(ThisPoint[iDim],iDim);
		Index *= mN[iDim-1];
	    }
	    Index += CalculatePosition(ThisPoint[0],0);
	    
	    return Index;
	}

	//************************************************************************

	IndexType CalculateIndex( CellType const& ThisIndex )
	{
	    IndexType Index = 0;
	    for(SizeType iDim = TDimension-1 ; iDim > 0 ; iDim--)
	    {
		Index += ThisIndex[iDim];
		Index *= mN[iDim-1];
	    }
	    Index += ThisIndex[0];
	    
	    return Index;
	}

	//************************************************************************

	CellType CalculateCell( PointType const& ThisPoint )
	{
	    CellType Cell;
	    for(SizeType i = 0 ; i < TDimension ; i++)
		Cell[i] = CalculatePosition(ThisPoint[i],i);
	    
	    return Cell;
	}

	CellType CalculateCell( PointType const& ThisPoint, CoordinateType Radius )
	{
	    CellType Cell;
	    for(SizeType i = 0 ; i < TDimension ; i++)
		Cell[i] = CalculatePosition(ThisPoint[i]+Radius,i);
	    
	    return Cell;
	}

	//************************************************************************

    public:

	//************************************************************************
     
	void generateBins() 
	{
	    GenerateBins();
	}

	//************************************************************************
     
	void prepareData(PointType const& ThisPoint) 
	{
	    static int index = 0;
	  
	    pointsToSearch[index].x = ThisPoint[0];
	    pointsToSearch[index].y = ThisPoint[1];
	    pointsToSearch[index].z = ThisPoint[2];
	  
	    index++;
	    if (index == mProblemSize) 
		index = 0;
	}
  
	//************************************************************************
	  
	void allocateOCLBuffers(int ConcurrentPoints, int maxResults)
	{
	    int ConcurrentPointsReal = ConcurrentPoints;
	    if(ConcurrentPoints > mProblemSize) 
	      ConcurrentPointsReal = mProblemSize;
	    OCL_Radius         = OCLDeviceGroup->CreateBuffer(sizeof(double), CL_MEM_READ_ONLY);
	    OCL_Radius2        = OCLDeviceGroup->CreateBuffer(sizeof(double), CL_MEM_READ_ONLY);
	    OCL_maxResults     = OCLDeviceGroup->CreateBuffer(sizeof(int), CL_MEM_READ_ONLY);
	    OCL_PointsToSearch = OCLDeviceGroup->CreateBuffer(sizeof(cl_double4) * ConcurrentPoints, CL_MEM_READ_ONLY);
	    OCL_results        = OCLDeviceGroup->CreateBuffer(sizeof(int) * ConcurrentPointsReal, CL_MEM_WRITE_ONLY);
	    OCL_outData        = OCLDeviceGroup->CreateBuffer(sizeof(int) * ConcurrentPointsReal * maxResults, CL_MEM_WRITE_ONLY);
	    OCL_distance       = OCLDeviceGroup->CreateBuffer(sizeof(double) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_NFunction      = OCLDeviceGroup->CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_PointsTriangle = OCLDeviceGroup->CreateBuffer(sizeof(cl_double4) * mPointSize, CL_MEM_READ_ONLY);
	    OCL_TriangleList   = OCLDeviceGroup->CreateBuffer(sizeof(cl_int4) * mTriangleSize, CL_MEM_READ_ONLY);
	    OCL_uVector        = OCLDeviceGroup->CreateBuffer(sizeof(double) * mProblemSize, CL_MEM_READ_ONLY);
	    std::cout << "DEBUG: allocateOCLBuffers - OK" << std::endl;
	}
         
	void searchInRadiusOCL(double Radius, int ConcurrentPoints, int maxResults) 
	{
	    int processed = 0;
	    
	    int result = 0;
	    
	    double HOST_memRadius  = Radius;
	    double HOST_memRadius2 = Radius * Radius;
	    
	    int amount;
	
	    while (processed < mProblemSize)
	    {	  
		amount = (ConcurrentPoints > mProblemSize) ? mProblemSize : (processed + ConcurrentPoints) < mProblemSize ? ConcurrentPoints : mProblemSize - processed;
		
		int HOST_w_size = amount;
		int * results;
		int * resultPoints;
		
		results      = (int *)malloc(sizeof(int) * amount);
		resultPoints = (int *)malloc(sizeof(int) * amount * maxResults);

		OCLDeviceGroup->CopyBuffer(OCL_Radius        , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius));
		OCLDeviceGroup->CopyBuffer(OCL_Radius2       , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius2));
		OCLDeviceGroup->CopyBuffer(OCL_PointsToSearch, OpenCL::HostToDevice, OpenCL::VoidPList(1,&pointsToSearch[processed]));
		OCLDeviceGroup->CopyBuffer(OCL_maxResults    , OpenCL::HostToDevice, OpenCL::VoidPList(1,&maxResults));
				
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 0,  OCL_IndexCellReferenceO);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 1,  OCL_BinsContainer);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 2 , OCL_InvCellSize);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 3 , OCL_N);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 4,  OCL_Radius);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 5,  OCL_Radius2);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 6,  OCL_PointsToSearch);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 7,  OCL_MinPoint);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 8,  OCL_outData);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 9,  OCL_results);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 10, OCL_maxResults);

		OCLDeviceGroup->ExecuteKernel(OCLSearchInRadius, amount);
		
		OCLDeviceGroup->CopyBuffer(OCL_results, OpenCL::DeviceToHost, OpenCL::VoidPList(1,results));
		OCLDeviceGroup->CopyBuffer(OCL_outData, OpenCL::DeviceToHost, OpenCL::VoidPList(1,resultPoints));
		
		for(int j = 0; j < amount; j++)
		    result += results[j];
		
		free(results);
		free(resultPoints);
		
		processed += amount;
		//std::cout << "Total Results while " << processed << " points done: " << result << std::endl;
	    }
	    
	    std::cout << "Total Results: " << result << " (" << maxResults * SearchUtils::PointerDistance(mPointBegin, mPointEnd) << " stored)"<< std::endl;
	}
	
	void searchTriangles(double Radius, int ConcurrentPoints, int maxResults) 
	{
	    int processed = 0;
	    int result = 0;
	    
	    double HOST_memRadius  = Radius;
	    double HOST_memRadius2 = Radius * Radius;
	    
	    int amount;
	    int searchSpaceSize = mProblemSize;
	    
	    cl_double4 * PointsTriangles = (cl_double4 *)malloc(mPointSize * sizeof(cl_double4));
	    double * uVector = (double *)malloc(mProblemSize * sizeof(double));
	    
	    for( IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
	    {
		static int k = 0;
		PointsTriangles[k].x = (**Point)[0];
		PointsTriangles[k].y = (**Point)[1];
		PointsTriangles[k].z = (**Point)[2];
		PointsTriangles[k].w = (**Point).id;
		k++;
	    }
	    
	    for(int i = 0; i < mProblemSize; i++)
	    {
		uVector[i] = pointsToSearch[i].y;
	    }
	    
	    OCLDeviceGroup->CopyBuffer(OCL_Radius        , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius));
	    OCLDeviceGroup->CopyBuffer(OCL_Radius2       , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius2));
	    OCLDeviceGroup->CopyBuffer(OCL_PointsTriangle, OpenCL::HostToDevice, OpenCL::VoidPList(1,PointsTriangles));
	    OCLDeviceGroup->CopyBuffer(OCL_TriangleList  , OpenCL::HostToDevice, OpenCL::VoidPList(1,mTriangles));
	    OCLDeviceGroup->CopyBuffer(OCL_uVector       , OpenCL::HostToDevice, OpenCL::VoidPList(1,uVector));
	
	    while (processed < mProblemSize)
	    {	 
		amount = (ConcurrentPoints > mProblemSize) ? mProblemSize : (processed + ConcurrentPoints) < mProblemSize ? ConcurrentPoints : mProblemSize - processed;

		cl_double4 * testResults = (cl_double4 *)malloc(sizeof(cl_double4) * amount);
		
		OCLDeviceGroup->CopyBuffer(OCL_PointsToSearch, OpenCL::HostToDevice, OpenCL::VoidPList(1,&pointsToSearch[processed]));

		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTriangles, 0,  OCL_IndexCellReferenceO);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTriangles, 1,  OCL_BinsContainer);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTriangles, 2,  OCL_PointsTriangle);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTriangles, 3,  OCL_TriangleList);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTriangles, 4,  OCL_InvCellSize);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTriangles, 5 , OCL_N);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTriangles, 6,  OCL_Radius);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTriangles, 7,  OCL_Radius2);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTriangles, 8,  OCL_PointsToSearch);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTriangles, 9,  OCL_MinPoint);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTriangles, 10, OCL_uVector);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTriangles, 11, OCL_NFunction); //Interpolation of N; for each point

		OCLDeviceGroup->ExecuteKernel(OCLSearchTriangles, amount);

		OCLDeviceGroup->CopyBuffer(OCL_NFunction, OpenCL::DeviceToHost, OpenCL::VoidPList(1,testResults));

		for(int j = 0; j < 100; j++)
		{
		    //if(uVector[j] != (testResults[j].x + testResults[j].y + testResults[j].z))
		    std::cout << "Expected:\t" << uVector[j] << "\t" << (testResults[j].x + testResults[j].y + testResults[j].z + testResults[j].w) << "\t\tResutl" << std::endl;
		}
			
		//std::cout << "Processed: " << processed << " Amount: " << amount << std::endl;

		processed += amount;
		
		std::cout << processed << " points processed" << std::endl;
		
		free(testResults);

		//std::cout << "Total Results while " << processed << " points done: " << result << std::endl;
	    }
	    
	    free(PointsTriangles);
	    free(uVector);
	    
	    //std::cout << "Total Results: " << result << " (" << maxResults * SearchUtils::PointerDistance(mPointBegin, mPointEnd) << " stored)"<< std::endl;
	}
         
	void searchNearestOCL(double Radius, int ConcurrentPoints) 
	{
	    int processed = 0;
	    
	    int result = 0;
	    
	    double HOST_memRadius  = Radius;
	    double HOST_memRadius2 = Radius * Radius;
	    
	    int amount;
	
	    while (processed < mProblemSize)
	    {	  
		amount = (ConcurrentPoints > mProblemSize) ? mProblemSize : ((processed + ConcurrentPoints) < mProblemSize) ? ConcurrentPoints : mProblemSize - processed;
		
		int HOST_w_size = amount;
		int * resultPoints;
		
		resultPoints = (int *)malloc(sizeof(int) * amount);

		OCLDeviceGroup->CopyBuffer(OCL_Radius        , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius));
		OCLDeviceGroup->CopyBuffer(OCL_Radius2       , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius2));
		OCLDeviceGroup->CopyBuffer(OCL_PointsToSearch, OpenCL::HostToDevice, OpenCL::VoidPList(1,&pointsToSearch[processed]));

		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearest, 0,  OCL_IndexCellReferenceO);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearest, 1,  OCL_BinsContainer);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearest, 2 , OCL_InvCellSize);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearest, 3 , OCL_N);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearest, 4,  OCL_Radius);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearest, 5,  OCL_Radius2);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearest, 6,  OCL_PointsToSearch);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearest, 7,  OCL_MinPoint);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearest, 8,  OCL_distance);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearest, 9,  OCL_results);

		OCLDeviceGroup->ExecuteKernel(OCLSearchNearest, amount);
	
		OCLDeviceGroup->CopyBuffer(OCL_results, OpenCL::DeviceToHost, OpenCL::VoidPList(1,resultPoints));

		result = resultPoints[amount-1];

		free(resultPoints);

		processed += amount;
	    }
	    
	    std::cout << "Nearest Point ID: " << result << std::endl;
	}
         
         
	void searchNearestOCLCubic(double Radius, int ConcurrentPoints) 
	{
	    int processed = 0;
	    
	    int result = 0;
	    
	    double HOST_memRadius  = Radius;
	    double HOST_memRadius2 = Radius * Radius;
	    
	    int amount;

	    while (processed < mProblemSize)
	    {	  
		amount = (ConcurrentPoints > mProblemSize) ? mProblemSize : (processed + ConcurrentPoints) < mProblemSize ? ConcurrentPoints : mProblemSize - processed;
		
		int HOST_w_size = amount;
		int * resultPoints;
		
		resultPoints = (int *)malloc(sizeof(int) * amount);

		OCLDeviceGroup->CopyBuffer(OCL_Radius        , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius));
		OCLDeviceGroup->CopyBuffer(OCL_Radius2       , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius2));
		OCLDeviceGroup->CopyBuffer(OCL_PointsToSearch, OpenCL::HostToDevice, OpenCL::VoidPList(1,&pointsToSearch[processed]));

		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearestCubic, 0,  OCL_IndexCellReferenceO);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearestCubic, 1,  OCL_BinsContainer);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearestCubic, 2 , OCL_InvCellSize);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearestCubic, 3 , OCL_N);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearestCubic, 4,  OCL_Radius);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearestCubic, 5,  OCL_Radius2);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearestCubic, 6,  OCL_PointsToSearch);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearestCubic, 7,  OCL_MinPoint);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearestCubic, 8,  OCL_distance);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchNearestCubic, 9,  OCL_results);

		OCLDeviceGroup->ExecuteKernel(OCLSearchNearestCubic, amount);

		OCLDeviceGroup->CopyBuffer(OCL_results, OpenCL::DeviceToHost, OpenCL::VoidPList(1,resultPoints));

		result = resultPoints[amount-1];

		free(resultPoints);

		processed += amount;
	    }
	    
	    std::cout << "Nearest Point ID: " << result << std::endl;
	}

	 //************************************************************************
	 //************************************************************************

	 /// Turn back information as a string.
	 virtual std::string Info() const
	 {
	   return "BinsContainer";
	 }

	 /// Print information about this object.
	 virtual void PrintInfo(std::ostream& rOStream) const
	 {
	   rOStream << "BinsContainer";
	 }

	 /// Print object's data.
	 virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const
	 {
	   rOStream << Perfix << "Bin[" << SearchUtils::PointerDistance(mPointBegin, mPointEnd) << "] : " << std::endl;
	   for(IteratorConstIterator i_cell = mIndexCell.begin() ; i_cell != mIndexCell.end()-1 ; i_cell++)
	   {
		 rOStream << Perfix << "[ " ;
		 for(IteratorType i_point = *i_cell ; i_point != *(i_cell+1) ; i_point++)
		   rOStream << **i_point << " ";
		 rOStream << "]" << std::endl;
	   }
	   rOStream << std::endl;
	 }

	 /// Print Size of Container
	 void PrintSize( std::ostream& rout ){
	   rout << " BinsSize: ";
	   for(SizeType i = 0 ; i < TDimension ; i++)
		 rout << "[" << mN[i] << "]";
	   rout << std::endl;
	 }

	 /// Print Limits Points of the Container
	 void PrintBox( std::ostream& rout ){
	   rout << " BinsBox: Min [";  mMinPoint.Print(rout);
	   rout <<       "];  Max [";  mMaxPoint.Print(rout);
	   rout <<       "];  Size ["; mCellSize.Print(rout);
	   rout << "]" << std::endl;
	 }

	 /// Assignment operator.
	 BinsOCL& operator=(BinsOCL const& rOther);

	 /// Copy constructor.
	 BinsOCL(BinsOCL const& rOther);

	 TPointType GetMinPoint()
	 {
	   TPointType point;
	   for(SizeType i = 0 ; i < TDimension ; i++)
		 point[i] = mMinPoint[i];
	   return point;
	 }

	 TPointType GetMaxPoint()
	 {
	   TPointType point;
	   for(SizeType i = 0 ; i < TDimension ; i++)
		 point[i] = mMaxPoint[i];
	   return point;
	 }

    private:
      
    // OCL mem and auxiliars
	Kratos::OpenCL::DeviceGroup *OCLDeviceGroup;
	
	cl_int Err;
	
	cl_uint OCL_program;
	
	//Kernels
	cl_uint OCLGenerateBinsA;
	cl_uint OCLGenerateBinsB; 
	cl_uint OCLGenerateBinsB1; 
	cl_uint OCLGenerateBinsB2; 
	cl_uint OCLGenerateBinsB3;
	cl_uint OCLGenerateBinsC; 
	cl_uint OCLSearchInRadius;
	cl_uint OCLSearchTriangles;
	cl_uint OCLSearchNearest;
	cl_uint OCLSearchNearestCubic;
	      
	cl_uint OCL_PointsBins;
	cl_uint OCL_IndexCell;
	cl_uint OCL_IndexCellReference; 
	cl_uint OCL_IndexCellReferenceO;
	cl_uint OCL_InvCellSize;
	cl_uint OCL_N;
	cl_uint OCL_MinPoint; 
	cl_uint OCL_BinsContainer; 
	cl_uint OCL_IndexCellReferenceU;
	cl_uint OCL_PowSize;
	cl_uint OCL_Offset;
	cl_uint OCL_PartialReduction;
      
	cl_uint OCL_Radius; 
	cl_uint OCL_Radius2; 
	cl_uint OCL_PointsToSearch;
	cl_uint OCL_PointsTriangle;
	cl_uint OCL_TriangleList;
	cl_uint OCL_outData;
	cl_uint OCL_results; 
	cl_uint OCL_distance; 
	cl_uint OCL_resultsNum; 
	cl_uint OCL_w_size;
	cl_uint OCL_maxResults;
	cl_uint OCL_NFunction;
	cl_uint OCL_uVector;
	
	cl_double4 * pointsToSearch; 
	
	int mCellSizeInt;
	int mPointSize;
	int mTriangleSize;
	int mProblemSize;
    
    // Point and Barycenter Access Iterators (vector reordered!!)
	IteratorType    mPointBegin;
	IteratorType	mPointEnd;

    // Bin Parameters (Sizes,BoundingBox,...)
	PointType  mMinPoint;
	PointType  mMaxPoint;
	PointType  mCellSize;
	PointType  mInvCellSize;
	
	cl_int4 * mTriangles;
	PointType * mBaryCenters;
	
	Tvector<SizeType,TDimension>  mN;

    // Bins Access Vector ( vector<Iterator> )
	IteratorVector           mIndexCell;
	IteratorVector           mIndexCellOCL;
	IteratorIterator         mIndexCellBegin;
	IteratorIterator         mIndexCellEnd;

    // Work Variables ( For non-copy of Search Variables )
    //SearchStructureType mBox;

   public:

	 static TreeNodeType* Construct(IteratorType PointsBegin, IteratorType PointsEnd, PointType MaxPoint, PointType MinPoint, SizeType BucketSize)
	 {
	   SizeType number_of_points = SearchUtils::PointerDistance(PointsBegin,PointsEnd);
	   if (number_of_points == 0)
		 return NULL;
	   else 
	   {
		 return new BinsOCL( PointsBegin, PointsEnd, MinPoint, MaxPoint, BucketSize );
	   }
	 }


};

template< std::size_t TDimension, class TPointType, class TContainerType, class TPointerType,
          class TIteratorType, class TDistanceIteratorType, class TDistanceFunction >
std::ostream & operator<<( std::ostream& rOStream, BinsOCL<TDimension,TPointType,TContainerType,TPointerType,TIteratorType,TDistanceIteratorType,TDistanceFunction>& rThis)
{
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintSize(rOStream);
	rThis.PrintData(rOStream);
	return rOStream;
}



}

#endif // KRATOS_BINS_CONTAINER_H_INCLUDE
