//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-27 17:02:19 $
//   Revision:            $Revision: 1.1.1.1 $
//
//

#if !defined(KRATOS_BINS_OCL_CONTAINER_H_INCLUDE)
#define KRATOS_BINS_OCL_CONTAINER_H_INCLUDE

#define WORKGROUP_SIZE 128

#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
#define KRATOS_CL_4_ARRAY_X(Arr,p) (Arr[p].x)
#define KRATOS_CL_4_ARRAY_Y(Arr,p) (Arr[p].y)
#define KRATOS_CL_4_ARRAY_Z(Arr,p) (Arr[p].z)
#define KRATOS_CL_4_ARRAY_W(Arr,p) (Arr[p].w)

#define KRATOS_CL_4_X(Arr) (Arr.x)
#define KRATOS_CL_4_Y(Arr) (Arr.y)
#define KRATOS_CL_4_Z(Arr) (Arr.z)
#define KRATOS_CL_4_W(Arr) (Arr.w)
#else
#define KRATOS_CL_4_ARRAY_X(Arr,p) (Arr[p].s[0])
#define KRATOS_CL_4_ARRAY_Y(Arr,p) (Arr[p].s[1])
#define KRATOS_CL_4_ARRAY_Z(Arr,p) (Arr[p].s[2])
#define KRATOS_CL_4_ARRAY_W(Arr,p) (Arr[p].s[3])

#define KRATOS_CL_4_X(Arr) (Arr.s[0])
#define KRATOS_CL_4_Y(Arr) (Arr.s[1])
#define KRATOS_CL_4_Z(Arr) (Arr.s[2])
#define KRATOS_CL_4_W(Arr) (Arr.s[3])
#endif

#include "../../../kratos/spatial_containers/tree.h"
#include "../custom_utilities/opencl_interface.h"
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
	 
	BinsOCL( IteratorType const& PointBegin, IteratorType const& PointEnd, cl_int4 * const& indexArray, int TSize, int sampleSize, SizeType BucketSize = 1 )
	    : mPointBegin(PointBegin), mPointEnd(PointEnd)
	{
	    if(mPointBegin==mPointEnd)
		return;
	    
	    mPointSize = SearchUtils::PointerDistance(mPointBegin, mPointEnd);
	    mTriangleSize = TSize;

	    GenerateSampleInput(sampleSize);   
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
	    int numReps = 1;
	    
	    mProblemSize = newSize * newSize * numReps;
	    pointsToSearch = new cl_double4[mProblemSize];
	    
	    
	    int index = 0;
	    for(int rep = 0; rep < numReps; rep++)
	    {
	      	double baseX = -0.98600;
		double baseY = -0.98500;
		double incX = (1.9940 + 0.99600) / newSize;
		double incY = (0.99518 + 0.99500) / newSize;

// 	      	double baseX = 0.01f;
// 		double baseY = 0.01f;
// 		double incX = (0.995f - 0.01f) / newSize;
// 		double incY = (0.995f - 0.01f) / newSize;

		double X = baseX;
		double Y = baseY;
		for(int i = 0; X < 1.9940f; i++)
// 		for(int i = 0; X < 0.995f; i++)
		{
		    Y = baseY;
		    for(int j = 0; Y < 0.99518f; j++)
// 		    for(int j = 0; Y < 0.995f; j++)
		    {
			//for(int k = 0; k < mProblemSize; k++)
			//{
			if(index == 0){

			    KRATOS_CL_4_ARRAY_X(pointsToSearch,index) = baseX;
			    KRATOS_CL_4_ARRAY_Y(pointsToSearch,index) = baseY;
			    KRATOS_CL_4_ARRAY_Z(pointsToSearch,index) = 0.0003f;
// 			    KRATOS_CL_4_ARRAY_Z(pointsToSearch,index) = 0;
			    KRATOS_CL_4_ARRAY_W(pointsToSearch,index) = index++;
			} else 
			{

			    KRATOS_CL_4_ARRAY_X(pointsToSearch,index) = X;
			    KRATOS_CL_4_ARRAY_Y(pointsToSearch,index) = Y;
			    KRATOS_CL_4_ARRAY_Z(pointsToSearch,index) = 0.0003f;
// 			    KRATOS_CL_4_ARRAY_Z(pointsToSearch,index) = 0;
			    KRATOS_CL_4_ARRAY_W(pointsToSearch,index) = index++;
			}
			Y += incY;
			//}
		    }
		    X += incX;
		}
	    }
	    
	    std::cout << "Sample of:" << newSize << "x" << newSize << " (" << mProblemSize << ") elements" << std::endl;
	    
	}
	
	//************************************************************************
	 
    private:

	//************************************************************************
     
	void InitOCL() 
	{
	   
	   //"Advanced Micro Devices, Inc."
	   //"NVIDIA Corporation"
	   
	    OCLDeviceGroup = new Kratos::OpenCL::DeviceGroup(CL_DEVICE_TYPE_GPU,true,"");
	    std::cout << "Found " << OCLDeviceGroup->DeviceNo << " device(s)." << std::endl;
	    for (cl_uint i = 0; i < OCLDeviceGroup->DeviceNo; i++)
	    {
		std::cout << "  Device " << i << ": " << Kratos::OpenCL::DeviceTypeString(OCLDeviceGroup->DeviceTypes[i]) << std::endl;
	    }
	    
	    char params[512] = "";
	    sprintf(params,"-D POINT_SIZE=%d -D TDIMENSION=%Zu -D CELL_SIZE=%d -D WORKGROUP_SIZE=%d",mTriangleSize ? mTriangleSize : mPointSize,TDimension,mCellSizeInt+1,WORKGROUP_SIZE);
	      
	    OCL_program = OCLDeviceGroup->BuildProgramFromFile("binshashOPT.cl",params);  
	}
	 
	//************************************************************************
	 
	//Calculate the barycenter of the triangle/tetraedre formed by index in index4
	PointType BaryCenter(cl_int4 const& index4)
	{
	    IteratorType Point = mPointBegin;
	    PointType a,b,c,d,r;

 	    a = (**(Point = mPointBegin + KRATOS_CL_4_X(index4)-1));
 	    b = (**(Point = mPointBegin + KRATOS_CL_4_Y(index4)-1));
 	    c = (**(Point = mPointBegin + KRATOS_CL_4_Z(index4)-1));
	    
	    //3D
// 	    if(index4.w > 0)

	    if(KRATOS_CL_4_W(index4) > 0)
	    {	
// 		d = (**(Point = mPointBegin + index4.w-1));

		d = (**(Point = mPointBegin + KRATOS_CL_4_W(index4)-1));
		
		for(int i = 0; i < TDimension; i++)
		{
		  r[i] = (a[i] + b[i] + c[i] + d[i]) / 4; // * 0.25f ? 
		}
	    }
	    //2D
	    else
	    {
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
	 	
	//Calculates the bounding box for barycenters of tringle/tetraedre input
	void CalculateBoundingBox(cl_int4 * const& indexArray) 
	{
	    mBaryCenters = (PointType *)malloc(sizeof(PointType) * mTriangleSize);
	    mTriangles = (cl_int4 *)malloc(sizeof(cl_int4) * mTriangleSize); 
	    
	    PointType barycenter = BaryCenter(indexArray[0]);
	    
	    for(SizeType i = 0 ; i < TDimension ; i++)
	    {
		mMinPoint[i] = barycenter[i];
		mMaxPoint[i] = barycenter[i];
	    }
	    
	    mBaryCenters[0] = barycenter;
	    mBaryCenters[0].id = 0;
	    
	    mTriangles[0] = indexArray[0];
	   
	    for(int j = 1; j < mTriangleSize; j++)
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
		
		mTriangles[j] = indexArray[j];
	    }
	   
	    // Check que se ha echo todo bien
	    std::cout << "Min: ";
	    for(SizeType i = 0 ; i < TDimension ; i++)
		std::cout << mMinPoint[i] << "  ";
	    std::cout << std::endl;
	    
	    std::cout << "Max: ";
	    for(SizeType i = 0 ; i < TDimension ; i++)
		std::cout << mMaxPoint[i] << "  ";
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
		delta[i] = (mMaxPoint[i] - mMinPoint[i]);
		if ( delta[i] > delta[index] )
		    index = i;
		delta[i] = (delta[i] == 0.00) ? 0.01 : delta[i]; //*** 1.00 -> 0.01
	    }
	    
	    std::cout << "Delta: ";
	    for(SizeType i = 0 ; i < TDimension ; i++)
		std::cout << delta[i] << "  ";
	    std::cout << std::endl;

	    for(SizeType i = 0 ; i < TDimension ; i++){
		alpha[i] = delta[i] / delta[index];
		mult_delta *= alpha[i];
	    }
	    
	    std::cout << "MultDelta: " << mult_delta << std::endl;
	    
	    std::cout << "Alpha: ";
	    for(SizeType i = 0 ; i < TDimension ; i++)
		std::cout << alpha[i] << "  ";
	    std::cout << std::endl;

	    mN[index] = static_cast<SizeType>( pow(static_cast<CoordinateType>((mTriangleSize ? mTriangleSize : mPointSize)/mult_delta), 1.00/TDimension)+1 );
	   
	    std::cout << "Index: " << index << std::endl;
	    std::cout << "mN[Index]: " << mN[index] << std::endl;
 	    
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
	    
	    std::cout << "N: " << mN[0] << " " << mN[1] << " " << mN[2] << std::endl;
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
	    OCLSearchTriangles   = OCLDeviceGroup->RegisterKernel(OCL_program,"SearchTriangle3D");
// 	    OCLSearchTrianglesB   = OCLDeviceGroup->RegisterKernel(OCL_program,"SearchTriangle3DB");
	    OCLSearchNearest      = OCLDeviceGroup->RegisterKernel(OCL_program,"SearchNearestMultiple");
	    OCLSearchNearestCubic = OCLDeviceGroup->RegisterKernel(OCL_program,"SearchNearestMultipleCubic");
	    clock_gettime( CLOCK_REALTIME, &end );
	    
	    std::cout << "Kernel Reg:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
	    
	    int partialReduction = 0;
	    int Offset = 0;
	    int PowSize = 1;
	    int pr = 0;
	    
	    while(PowSize < mCellSizeInt+1)
		PowSize <<= 1;
	    
	    std::cout << "POW:\t\t" << PowSize << std::endl;
	    
	    cl_double4 * BinsContainer;
	    cl_double4 * PointsBins;
	    cl_double4 * MinPoint      = new cl_double4();
	    
	    int * IndexCellReference   = (int *)malloc( (PowSize) * sizeof(int) );
	    int * IndexCellReferenceO  = (int *)malloc( (PowSize) * sizeof(int) );
	    
	    double InvCellSize[TDimension];
	    double N[TDimension];
	        	    
	    clock_gettime( CLOCK_REALTIME, &begin );
	    
	    if(mTriangleSize)
	    {
	        BinsContainer = (cl_double4 *)malloc(mTriangleSize * sizeof(cl_double4));
	        PointsBins    = (cl_double4 *)malloc(mTriangleSize * sizeof(cl_double4));
		for(int i = 0; i < mTriangleSize; i++)
		{ 
		    KRATOS_CL_4_ARRAY_X(PointsBins,i) = mBaryCenters[i][0];
		    KRATOS_CL_4_ARRAY_Y(PointsBins,i) = mBaryCenters[i][1];
		    KRATOS_CL_4_ARRAY_Z(PointsBins,i) = mBaryCenters[i][2];
		    KRATOS_CL_4_ARRAY_W(PointsBins,i) = mBaryCenters[i].id;
		}
		OCL_PointsBins         = OCLDeviceGroup->CreateBuffer(sizeof(cl_double4) * mTriangleSize, CL_MEM_READ_ONLY);
		OCL_BinsContainer      = OCLDeviceGroup->CreateBuffer(sizeof(cl_double4) * mTriangleSize, CL_MEM_READ_WRITE);
	    }
	    else
	    {
	      	BinsContainer = (cl_double4 *)malloc(mPointSize * sizeof(cl_double4));
	        PointsBins    = (cl_double4 *)malloc(mPointSize * sizeof(cl_double4));
		for( IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
		{
		    static int k = 0;	    
		    KRATOS_CL_4_ARRAY_X(PointsBins,k) = (**Point)[0];
		    KRATOS_CL_4_ARRAY_Y(PointsBins,k) = (**Point)[1];
		    KRATOS_CL_4_ARRAY_Z(PointsBins,k) = (**Point)[2];
		    KRATOS_CL_4_ARRAY_W(PointsBins,k) = (**Point).id;
		    k++;
		}
		OCL_PointsBins         = OCLDeviceGroup->CreateBuffer(sizeof(cl_double4) * mPointSize, CL_MEM_READ_ONLY);
		OCL_BinsContainer      = OCLDeviceGroup->CreateBuffer(sizeof(cl_double4) * mPointSize, CL_MEM_READ_WRITE);
	    }

	    KRATOS_CL_4_X(MinPoint[0]) = mMinPoint[0];
	    KRATOS_CL_4_Y(MinPoint[0]) = mMinPoint[1];
	    KRATOS_CL_4_Z(MinPoint[0]) = mMinPoint[2];
	    
	    for(int i = 0; i < PowSize; i++)
		IndexCellReference[i] = 0;
	    
	    for(int i = 0; i < TDimension; i++) 
	    {
		InvCellSize[i] = mInvCellSize[i];
		N[i] = mN[i];
	    }
	    clock_gettime( CLOCK_REALTIME, &end );
	    
	    std::cout << "Var filling:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
	    
	    std::cout << "mProblemSize: " << mProblemSize << std::endl;
	    OCL_IndexCellReference = OCLDeviceGroup->CreateBuffer(sizeof(int) * PowSize, CL_MEM_READ_WRITE);
	    OCL_IndexCellReferenceO= OCLDeviceGroup->CreateBuffer(sizeof(int) * PowSize, CL_MEM_READ_WRITE);
	    OCL_IndexCellReferenceU= OCLDeviceGroup->CreateBuffer(sizeof(int) * PowSize, CL_MEM_READ_WRITE);
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
// 	    OCLDeviceGroup->CopyBuffer(OCL_PowSize           , OpenCL::HostToDevice, OpenCL::VoidPList(1,&wiz));
	    OCLDeviceGroup->CopyBuffer(OCL_Offset            , OpenCL::HostToDevice, OpenCL::VoidPList(1,&Offset));
	    OCLDeviceGroup->CopyBuffer(OCL_PartialReduction  , OpenCL::HostToDevice, OpenCL::VoidPList(1,&pr));
	
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsA, 0,  OCL_PointsBins);
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsA, 1,  OCL_IndexCellReference);
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsA, 2,  OCL_InvCellSize);
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsA, 3,  OCL_N);
	    OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBinsA, 4,  OCL_MinPoint);
	    
	    clock_gettime( CLOCK_REALTIME, &begin );
	    OCLDeviceGroup->ExecuteKernel(OCLGenerateBinsA, mTriangleSize ? mTriangleSize : mPointSize);
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
	    OCLDeviceGroup->ExecuteKernel(OCLGenerateBinsC, mTriangleSize ? mTriangleSize : mPointSize);
	    clock_gettime( CLOCK_REALTIME, &end );
	    
	    std::cout << "KernelC Ex:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
	    
	    clock_gettime( CLOCK_REALTIME, &begin );
	    OCLDeviceGroup->CopyBuffer(OCL_BinsContainer, OpenCL::DeviceToHost, OpenCL::VoidPList(1,BinsContainer));
	    clock_gettime( CLOCK_REALTIME, &end );
	    
	    std::cout << "Kernel Copy:\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
	    
	    OCLDeviceGroup->CopyBuffer(OCL_IndexCellReferenceO, OpenCL::DeviceToHost, OpenCL::VoidPList(1,IndexCellReferenceO));
	    
	    
// 	    for(int i = 0; i < mCellSizeInt+1; i++) {
// 	      std::cout << "IndexCellReKE " << i << "  " << IndexCellReferenceO[i] << std::endl;
// 	    }
// 	    std::cout << "------------------" << std::endl;
	    
	    
// 	    for(int i = 0; i < mTriangleSize; i++ ) 
// 	    {
// 	      std::cout << i << " " << BinsContainer[i].w << std::endl;
// 	    }
	    
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
	    
	    KRATOS_CL_4_ARRAY_X(pointsToSearch,index) = ThisPoint[0];
	    KRATOS_CL_4_ARRAY_Y(pointsToSearch,index) = ThisPoint[1];
	    KRATOS_CL_4_ARRAY_Z(pointsToSearch,index) = ThisPoint[2];
	  
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
	    OCL_amount 	       = OCLDeviceGroup->CreateBuffer(sizeof(int), CL_MEM_READ_ONLY);
	    OCL_PointsToSearch = OCLDeviceGroup->CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_uVector	       = OCLDeviceGroup->CreateBuffer(sizeof(double) * ConcurrentPointsReal, CL_MEM_READ_ONLY);
	    OCL_results        = OCLDeviceGroup->CreateBuffer(sizeof(int) * ConcurrentPointsReal, CL_MEM_WRITE_ONLY);
	    OCL_outData        = OCLDeviceGroup->CreateBuffer(sizeof(int) * ConcurrentPointsReal * maxResults, CL_MEM_WRITE_ONLY);
	    OCL_distance       = OCLDeviceGroup->CreateBuffer(sizeof(double) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_NFunction      = OCLDeviceGroup->CreateBuffer(sizeof(cl_double4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
	    OCL_PointsTriangle = OCLDeviceGroup->CreateBuffer(sizeof(cl_double4) * mPointSize, CL_MEM_READ_ONLY);
	    OCL_TriangleList   = OCLDeviceGroup->CreateBuffer(sizeof(cl_int4) * mTriangleSize, CL_MEM_READ_ONLY);
	    OCL_searchIndex    = OCLDeviceGroup->CreateBuffer(sizeof(cl_int4) * ConcurrentPointsReal, CL_MEM_READ_WRITE);
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
// 	    std::cout << "SearchTriangle start" << std::endl;
	    
	    int processed = 0;
	    int result = 0;
	    
	    double HOST_memRadius  = Radius;
	    double HOST_memRadius2 = Radius * Radius;
	    
	    int amount;
	    
// 	    std::cout << "Normal var declared" << std::endl;
	    
	    cl_double4 * PointsTriangles = (cl_double4 *)malloc(sizeof(cl_double4) * mPointSize);
	    double * uVector = (double *)malloc(sizeof(double) * mProblemSize);
	    
// 	    std::cout << "Dynamic Var declared and allocated" << std::endl;
	    
	    int k = 0;
	    for( IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
	    {
		KRATOS_CL_4_ARRAY_X(PointsTriangles,k) = (**Point)[0];
		KRATOS_CL_4_ARRAY_Y(PointsTriangles,k) = (**Point)[1];
		KRATOS_CL_4_ARRAY_Z(PointsTriangles,k) = (**Point)[2];
		KRATOS_CL_4_ARRAY_W(PointsTriangles,k) = (**Point).id;
		k++;
	    }
	    
// 	    std::cout << "PointsTriangles filled" << std::endl;
// 	    std::cout << "Problem size " <<  mProblemSize << std::endl;
	    for(int i = 0; i < mProblemSize; i++)
	    {
// 	      	uVector[i] = pointsToSearch[i].x;

		uVector[i] = KRATOS_CL_4_ARRAY_Y(pointsToSearch,i);
	    }
	    
// 	    std::cout << "uVector filled" << std::endl;
	    
	    OCLDeviceGroup->CopyBuffer(OCL_Radius        , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius));
	    OCLDeviceGroup->CopyBuffer(OCL_Radius2       , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius2));
	    OCLDeviceGroup->CopyBuffer(OCL_PointsTriangle, OpenCL::HostToDevice, OpenCL::VoidPList(1,PointsTriangles));
	    OCLDeviceGroup->CopyBuffer(OCL_TriangleList  , OpenCL::HostToDevice, OpenCL::VoidPList(1,mTriangles));
	
// 	    std::cout << "OCL buffers filled" << std::endl;
	    
	    while (processed < mProblemSize)
	    {	 
		amount = (ConcurrentPoints > mProblemSize) ? mProblemSize : (processed + ConcurrentPoints) < mProblemSize ? ConcurrentPoints : mProblemSize - processed;

		cl_double4 * testResults = (cl_double4 *)malloc(sizeof(cl_double4) * amount);
	
		
		//Prevent from copying elements out of bounds
		cl_double4 * pointBuffer = (cl_double4 *)malloc(sizeof(cl_double4) * ConcurrentPoints);
		
		int pointCounter = 0;
		while(pointCounter < amount && pointCounter < ConcurrentPoints ) 
		{
		  pointBuffer[pointCounter] = pointsToSearch[processed + pointCounter];
		  pointCounter++;
		}
	
		OCLDeviceGroup->CopyBuffer(OCL_PointsToSearch, OpenCL::HostToDevice, OpenCL::VoidPList(1,&pointBuffer[0]));
		OCLDeviceGroup->CopyBuffer(OCL_uVector	     , OpenCL::HostToDevice, OpenCL::VoidPList(1,&uVector[processed]));
		OCLDeviceGroup->CopyBuffer(OCL_amount        , OpenCL::HostToDevice, OpenCL::VoidPList(1,&amount));

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
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTriangles, 11, OCL_NFunction);
		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTriangles, 12, OCL_amount);
		
   struct timespec begin;
   struct timespec end;
   
   clock_gettime( CLOCK_REALTIME, &begin );
		
		OCLDeviceGroup->ExecuteKernel(OCLSearchTriangles, amount);
		
   clock_gettime( CLOCK_REALTIME, &end );

// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesA, 0,  OCL_IndexCellReferenceO);
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesA, 1,  OCL_InvCellSize);
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesA, 2 , OCL_N);
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesA, 3,  OCL_Radius);
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesA, 4,  OCL_PointsToSearch);
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesA, 5,  OCL_MinPoint);
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesA, 6,  OCL_searchIndex);
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesA, 7,  OCL_amount);
// 
// 		OCLDeviceGroup->ExecuteKernel(OCLSearchTrianglesA, amount);
// 		
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesB, 0,  OCL_BinsContainer);
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesB, 1,  OCL_PointsTriangle);
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesB, 2,  OCL_TriangleList);
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesB, 3,  OCL_Radius2);
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesB, 4,  OCL_PointsToSearch);
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesB, 5,  OCL_uVector);		
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesB, 6,  OCL_NFunction);
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesB, 7,  OCL_searchIndex);
// 		OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchTrianglesB, 8,  OCL_amount);
// 		
// 		OCLDeviceGroup->ExecuteKernel(OCLSearchTrianglesB, amount);
		
		OCLDeviceGroup->CopyBuffer(OCL_NFunction, OpenCL::DeviceToHost, OpenCL::VoidPList(1,testResults));

// 		std::cout << "Check results..." << std::endl;
 		int errors = 0;
 		for(int j = 0; j < amount; j++)
		{
		    double diff = uVector[j+processed] - (testResults[j].x + testResults[j].y + testResults[j].z + testResults[j].w);
		    if((diff > 0.001 || diff < -0.001))
		    {
// 			std::cout << "(" <<pointsToSearch[j+processed].w << " " << pointsToSearch[j+processed].x << " " << pointsToSearch[j+processed].y << " " << pointsToSearch[j+processed].z << ") " << 
// 				     "Expected:\t" << " " << uVector[j+processed] << "\t" << (testResults[j].x + testResults[j].y + testResults[j].z + testResults[j].w ) << "\t\tResutl  - (" << 
// 				     testResults[j].x << " " << testResults[j].y << " " << testResults[j].z << " " << testResults[j].w << ") " << std::endl;
  		    	errors++;
		    }
		}
		std::cout << "Processed: " << amount << " Errors: " << errors << "  ";// <<  std::endl;
		 std::cout << "\t\t\t" << ((float)(end.tv_sec - begin.tv_sec) + (float)(end.tv_nsec-begin.tv_nsec)/1000000000) << std::endl;
		
		//std::cout << "Processed: " << processed << " Amount: " << amount << std::endl;

		processed += amount;
		
		//std::cout << processed << " points processed" << std::endl;
		
		free(testResults);

		//std::cout << "Total Results while " << processed << " points done: " << result << std::endl;
	    }
	    
	    k = 0;
	    for( IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
	    {
// 		PointsTriangles[k].x = NULL;
// 		PointsTriangles[k].y = NULL;
// 		PointsTriangles[k].z = NULL;
// 		PointsTriangles[k].w = NULL;
		
		KRATOS_CL_4_ARRAY_X(PointsTriangles,k) = NULL;
		KRATOS_CL_4_ARRAY_Y(PointsTriangles,k) = NULL;
		KRATOS_CL_4_ARRAY_Z(PointsTriangles,k) = NULL;
		KRATOS_CL_4_ARRAY_W(PointsTriangles,k) = NULL;
		
		k++;
	    }
	    
	    free(PointsTriangles);
	    //free(uVector);
	    
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
	cl_uint OCLSearchTrianglesB;
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
	cl_uint OCL_amount;
	cl_uint OCL_NFunction;
	cl_uint OCL_uVector;
	cl_uint OCL_searchIndex;
	
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
