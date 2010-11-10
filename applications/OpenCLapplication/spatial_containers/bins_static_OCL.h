//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-27 17:02:19 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_BINS_OCL_CONTAINER_H_INCLUDE)
#define KRATOS_BINS_OCL_CONTAINER_H_INCLUDE

#include "tree.h"
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
      
      Kratos::OpenCL::DeviceGroup *OCLDeviceGroup;
      
      cl_int Err;
      int Size;
      
      cl_uint OCLGenerateBins, OCLSearchInRadius;
      cl_mem OCL_Points, OCL_Cell, OCL_IndexCell, OCL_IndexCellReference, OCL_InvCellSize, OCL_N, OCL_MinPoint, OCL_BinsContainer;
      
      cl_double4 * pointsToSearch; 

   public:


	 //************************************************************************

	 // constructor 1
	 BinsOCL() : mPointBegin(this->NullIterator()), mPointEnd(this->NullIterator()) {};

	 //************************************************************************

	 BinsOCL( IteratorType const& PointBegin, IteratorType const& PointEnd, SizeType BucketSize = 1 )
	   : mPointBegin(PointBegin), mPointEnd(PointEnd)
	 {
	   if(mPointBegin==mPointEnd)
		 return;
	   InitOCL();
	   CalculateBoundingBox();
	   CalculateCellSize();
	   AllocateCellsContainer();
	   GenerateBins();
	 }

	 //************************************************************************

	 BinsOCL( IteratorType const& PointBegin, IteratorType const& PointEnd, PointType const& MinPoint, PointType const& MaxPoint, SizeType BucketSize = 1 )
	   : mPointBegin(PointBegin), mPointEnd(PointEnd)
	 {
	   if(mPointBegin==mPointEnd)
		 return;

	   for(SizeType i = 0 ; i < TDimension ; i++)
	   {
		 mMinPoint[i] = MinPoint[i];
		 mMaxPoint[i] = MaxPoint[i];
	   }
	   InitOCL();
	   CalculateCellSize();
	   AllocateCellsContainer();
	   GenerateBins();
	 }

	 //************************************************************************

	 BinsOCL( IteratorType const& PointBegin, IteratorType const& PointEnd, CoordinateType BoxSize, SizeType BucketSize = 1 )
	   : mPointBegin(PointBegin), mPointEnd(PointEnd)
	 {
	   if(mPointBegin==mPointEnd)
		 return;
	   InitOCL();
	   CalculateBoundingBox();
	   CalculateCellSize(BoxSize);
	   AllocateCellsContainer();
	   GenerateBins();
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
	 
   private:

	 //************************************************************************
     
	 void InitOCL() {
	   //Kratos::OpenCL::DeviceGroup OCLDeviceGroup(CL_DEVICE_TYPE_ALL);
	   //"Advanced Micro Devices, Inc."
	   //"NVIDIA Corporation"
	   OCLDeviceGroup = new Kratos::OpenCL::DeviceGroup(CL_DEVICE_TYPE_ALL,true,"NVIDIA Corporation");
	   std::cout << "Found " << OCLDeviceGroup->DeviceNo << " device(s)." << std::endl;
	   for (cl_uint i = 0; i < OCLDeviceGroup->DeviceNo; i++)
	   {
               std::cout << "  Device " << i << ": " << Kratos::OpenCL::DeviceTypeString(OCLDeviceGroup->DeviceTypes[i]) << std::endl;
	   }
	 }
	 
	 //************************************************************************

	 void CalculateBoundingBox() {

	   for(SizeType i = 0 ; i < TDimension ; i++){
		 mMinPoint[i] = (**mPointBegin)[i];
		 mMaxPoint[i] = (**mPointBegin)[i];
	   }
	   
	   for(IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
		 for(SizeType i = 0 ; i < TDimension ; i++){
		   if( (**Point)[i] < mMinPoint[i] ) mMinPoint[i] = (**Point)[i];
		   if( (**Point)[i] > mMaxPoint[i] ) mMaxPoint[i] = (**Point)[i];
		 }

	   
	   // Check que se ha echo todo bien
	   for(SizeType i = 0 ; i < TDimension ; i++)
	      std::cout << mMinPoint[i] << "  ";
           std::cout << std::endl;
	  
	 }

	 //************************************************************************

	 void CalculateCellSize() {
	   
	  CoordinateType delta[TDimension];
	  CoordinateType alpha[TDimension];
	  CoordinateType mult_delta = 1.00;
	  SizeType index = 0;
	  
	  for(SizeType i = 0 ; i < TDimension ; i++) {
	    delta[i] = mMaxPoint[i] - mMinPoint[i];
	    if ( delta[i] > delta[index] )
	      index = i;
	    delta[i] = (delta[i] == 0.00) ? 1.00 : delta[i];
	  }

	  for(SizeType i = 0 ; i < TDimension ; i++){
	    alpha[i] = delta[i] / delta[index];
	    mult_delta *= alpha[i];
	  }

	  mN[index] = static_cast<SizeType>( pow(static_cast<CoordinateType>(SearchUtils::PointerDistance(mPointBegin,mPointEnd)/mult_delta), 1.00/TDimension)+1 );
	   
	  for(SizeType i = 0 ; i < TDimension ; i++){
	    if(i!=index) {
	      mN[i] = static_cast<SizeType>(alpha[i] * mN[index]);
	      mN[i] = ( mN[i] == 0 ) ? 1 : mN[i];
	    }
	  }

	  for(SizeType i = 0 ; i < TDimension ; i++){
	    mCellSize[i] = delta[i] / mN[i];
	    mInvCellSize[i] = 1.00 / mCellSize[i];
	  }
	 
	 }

	 //************************************************************************

	 void CalculateCellSize( CoordinateType BoxSize ) {
	   for(SizeType i = 0 ; i < TDimension ; i++){
		 mCellSize[i] = BoxSize;
		 mInvCellSize[i] = 1.00 / mCellSize[i];
		 mN[i] = static_cast<SizeType>( (mMaxPoint[i]-mMinPoint[i]) / mCellSize[i]) + 1;
	   }
	 }

	 //************************************************************************

	 void AllocateCellsContainer() {
	   Size = 1;
	   for(SizeType i = 0 ; i < TDimension ; i++)
		 Size *= mN[i];
	   std::cout << "Size: \t" << Size << std::endl;
	   mIndexCell.resize(Size+1);
	   mIndexCellBegin = mIndexCell.begin();
	   mIndexCellEnd   = mIndexCell.end();
	   
	   std::cout << "Building Ocl program.." << std::endl;
	   
	   pointsToSearch = new cl_double4[mPointEnd - mPointBegin];
	   
	   std::cout << "Generate Params..." << std::endl;
	   
	   char params[512] = "";
	   sprintf(params,"-g -cl-mad-enable -D POINT_SIZE=%ld -D T_DIMENSION=%Zu -D CELL_SIZE=%d",(mPointEnd-mPointBegin),TDimension,Size+1);
	   
	   std::cout << "Building programs" << std::endl;	   
	   OCLDeviceGroup->BuildProgramFromFile("binshashOPT.cl",params);
	   
	   std::cout << "Register Kernels..." << std::endl;	   
	   OCLGenerateBins = OCLDeviceGroup->RegisterKernel(0,"GenerateBins");
	   OCLSearchInRadius = OCLDeviceGroup->RegisterKernel(0,"SearchInRadiusMultiple");
	   
	   int j = 0;
	   
	   cl_double4 * points = new cl_double4[(mPointEnd-mPointBegin)];
	   cl_double4 * MinPoint = new cl_double4();
	   
	   double InvCellSize[TDimension];
	   double N[TDimension];

	   for(IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++) {
		  points[j].x = (**Point)[0];
		  points[j].y = (**Point)[1];
		  points[j].z = (**Point)[2];
		  j++;
	   }
	   
	   MinPoint->x = mMinPoint[0];
	   MinPoint->y = mMinPoint[1];
	   MinPoint->z = mMinPoint[2];
	   
	   for(int i = 0; i < TDimension; i++) 
	   {
	      InvCellSize[i] = mInvCellSize[i];
	      N[i] = mN[i];
	   }
	   
	   int host_IndexCell[Size+1];
	   for(int i = 0; i < Size+1; i++)
	       host_IndexCell[i] = 0;
	   
	   // Prepare buffers
	   OCL_Points             = clCreateBuffer(OCLDeviceGroup->Contexts[0], CL_MEM_READ_ONLY, (mPointEnd-mPointBegin) * sizeof(cl_double4), NULL, &Err);
	   KRATOS_OCL_CHECK(Err);
	   
	   OCL_InvCellSize        = clCreateBuffer(OCLDeviceGroup->Contexts[0], CL_MEM_READ_ONLY, TDimension * sizeof(double), NULL, &Err);
	   KRATOS_OCL_CHECK(Err);
	   
	   OCL_N                  = clCreateBuffer(OCLDeviceGroup->Contexts[0], CL_MEM_READ_ONLY, TDimension * sizeof(double), NULL, &Err);
	   KRATOS_OCL_CHECK(Err);
	   
	   OCL_MinPoint           = clCreateBuffer(OCLDeviceGroup->Contexts[0], CL_MEM_READ_ONLY, sizeof(cl_double4), NULL, &Err);
	   KRATOS_OCL_CHECK(Err);
	   
	   OCL_Cell               = clCreateBuffer(OCLDeviceGroup->Contexts[0], CL_MEM_READ_WRITE, (Size+1) * sizeof(cl_double4), NULL, &Err);
	   KRATOS_OCL_CHECK(Err);
	   
	   OCL_IndexCell          = clCreateBuffer(OCLDeviceGroup->Contexts[0], CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, (Size+1) * sizeof(int), &host_IndexCell, &Err);
	   KRATOS_OCL_CHECK(Err);
	   
	   OCL_IndexCellReference = clCreateBuffer(OCLDeviceGroup->Contexts[0], CL_MEM_READ_WRITE, (Size+1) * sizeof(int), NULL, &Err);
	   KRATOS_OCL_CHECK(Err);
	   
	   OCL_BinsContainer      = clCreateBuffer(OCLDeviceGroup->Contexts[0], CL_MEM_READ_WRITE, (mPointEnd-mPointBegin) * sizeof(cl_double4), NULL, &Err);
	   KRATOS_OCL_CHECK(Err);
	     
	   // Load data into the input buffer
	   std::cout << "Loading data into inputbuffer..." << std::endl;
	   Err = clEnqueueWriteBuffer(OCLDeviceGroup->CommandQueues[0], OCL_Points,      CL_TRUE, 0, (mPointEnd-mPointBegin) * sizeof(cl_double4), points, 0, NULL, NULL);
	   KRATOS_OCL_CHECK(Err);
	   
	   Err = clEnqueueWriteBuffer(OCLDeviceGroup->CommandQueues[0], OCL_InvCellSize, CL_TRUE, 0, TDimension * sizeof(double), InvCellSize, 0, NULL, NULL);
	   KRATOS_OCL_CHECK(Err);
	   
	   Err = clEnqueueWriteBuffer(OCLDeviceGroup->CommandQueues[0], OCL_N,           CL_TRUE, 0, TDimension * sizeof(double), N, 0, NULL, NULL);
	   KRATOS_OCL_CHECK(Err);
	   
	   Err = clEnqueueWriteBuffer(OCLDeviceGroup->CommandQueues[0], OCL_MinPoint,    CL_TRUE, 0, sizeof(cl_double4), MinPoint, 0, NULL, NULL);
	   KRATOS_OCL_CHECK(Err);
	   
	   // Set arguments
	   std::cout << "Setting arguments..." << std::endl;
	   OCLDeviceGroup->SetKernelArg(OCLGenerateBins, 0,  OCL_Points);
	   OCLDeviceGroup->SetKernelArg(OCLGenerateBins, 1,  OCL_Cell);
	   OCLDeviceGroup->SetKernelArg(OCLGenerateBins, 2,  OCL_IndexCell);
	   OCLDeviceGroup->SetKernelArg(OCLGenerateBins, 3,  OCL_IndexCellReference);
	   OCLDeviceGroup->SetKernelArg(OCLGenerateBins, 4,  OCL_InvCellSize);
	   
	   OCLDeviceGroup->SetKernelArg(OCLGenerateBins, 5,  OCL_N);
	   OCLDeviceGroup->SetKernelArg(OCLGenerateBins, 6,  OCL_MinPoint);
	   OCLDeviceGroup->SetKernelArg(OCLGenerateBins, 7,  OCL_BinsContainer);
	   
	   std::cout << "Executing kernel..." << std::endl;
	   OCLDeviceGroup->ExecuteKernel(OCLGenerateBins, 1);
	   
	   // Release mem
	   clReleaseMemObject(OCL_Points);
	   clReleaseMemObject(OCL_IndexCell);
	 }

	 //************************************************************************

	 void GenerateBins( ){
	   PointVector TempPoint(mPointBegin,mPointEnd);

	   // Reset index vector
	   for( IteratorIterator Iter = mIndexCell.begin(); Iter != mIndexCell.end(); Iter++)
		 *Iter = mPointBegin;

	   // Update storage counter, storing ahead
	   for( IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
		 mIndexCell[ CalculateIndex(**Point) + 1 ]++;

	   // Storage/reshufing pass 1

	   // Update storage counter and store
	   for( IteratorIterator Iter = mIndexCell.begin()+1 ; Iter != mIndexCell.end() ; Iter++)
		 *Iter = *(Iter-1) + SearchUtils::PointerDistance(mPointBegin,*Iter);
	   
	   // Point pass 2
	   // Store the points in lbin1

	   // Update storage counter, storing in lbin1
	   for( PointIterator Point = TempPoint.begin() ; Point != TempPoint.end() ; Point++) {
		 //std::cout << "POINT:\t" << (**Point) << " \tHAVE INDEX:\t" << CalculateIndex(**Point) << std::endl;  
		 *(mIndexCell[CalculateIndex(**Point)]++) = *Point;
	   }
/* 
	   // TEST  !!! OJO -> No aumenta el contador del IteratorIterator !!!!
	   for( PointIterator Point = TempPoint.begin() ; Point != TempPoint.end() ; Point++)
	   {
	     Iter = mIndexCell[CalculateIndex(**Point)];
	     while( Iter != Point )
		 {
	       Iter2 = mIndexCell[CalculateIndex(**Iter)];
	       std::swap(*Iter,*Iter2);
	     }
	   }
*/  

	   // Storage/reshuffing pass 2

	   // Loop over bins, in reverse order
	   for(IteratorIterator Iter = mIndexCell.end()-1; Iter != mIndexCell.begin(); Iter--)
		 *Iter = *(Iter-1);
	   mIndexCell[0] = mPointBegin;
	   
	   //Esto checkea que se genere el bins bien
	   /*
	   cl_double4 OCL_CellHost[Size];
	   Err = clEnqueueReadBuffer(OCLDeviceGroup->CommandQueues[0], OCL_Cell, CL_TRUE, 0, sizeof(cl_double4) * (Size + 1), &OCL_CellHost, 0, NULL, NULL);

	   //Pointchecker
	   
	   int errors = 0;
	   std::cout << "Point checker" << std::endl;
	   for(int i = 0; i < Size; i++) 
	   {
	     if ((OCL_CellHost[i].x != (**mIndexCell[i])[0] || OCL_CellHost[i].y != (**mIndexCell[i])[1] || OCL_CellHost[i].z != (**mIndexCell[i])[2])) 
	     {
	       errors++;
	       std::cout << "!!!----------------------!!!" << std::endl;
	       std::cout << "Point D:  " << CalculateIndex(**mIndexCell[i]) << "\t" << (**mIndexCell[i])[0] << " " << (**mIndexCell[i])[1] << " " << (**mIndexCell[i])[2] << " " << std::endl;
	       std::cout << "Point CL: " << i << "\t" << OCL_CellHost[i].x << " " << OCL_CellHost[i].y << " " << OCL_CellHost[i].z << "  " << std::endl;
	     }
	   }
	   if(errors)
	     std::cout << "DETECTED: " << errors << " Error/s" << std::endl;
	   */
	   
	   // Release Memobjwcts that are not shared betwen kernels
	   clReleaseMemObject(OCL_Cell);
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
	   for(SizeType iDim = TDimension-1 ; iDim > 0 ; iDim--){
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
	   for(SizeType iDim = TDimension-1 ; iDim > 0 ; iDim--){
		 Index += ThisIndex[iDim];
		 Index *= mN[iDim-1];
	   }
	   Index += ThisIndex[0];
	   return Index;
	 }

	 //************************************************************************

	 CellType CalculateCell( PointType const& ThisPoint ){
	   CellType Cell;
	   for(SizeType i = 0 ; i < TDimension ; i++)
		 Cell[i] = CalculatePosition(ThisPoint[i],i);
	   return Cell;
	 }

	 CellType CalculateCell( PointType const& ThisPoint, CoordinateType Radius ){
	   CellType Cell;
	   for(SizeType i = 0 ; i < TDimension ; i++)
		 Cell[i] = CalculatePosition(ThisPoint[i]+Radius,i);
	   return Cell;
	 }

	 //************************************************************************

   public:

	 //************************************************************************
     //************************************************************************

     PointerType ExistPoint( PointerType const& ThisPoint, CoordinateType const Tolerance = static_cast<CoordinateType>(10.0*DBL_EPSILON) )
     {
       PointerType Nearest;
       CoordinateType Distance = static_cast<CoordinateType>(DBL_MAX);
       bool Found;
       SearchStructureType Box( CalculateCell(*ThisPoint,-Tolerance), CalculateCell(*ThisPoint,Tolerance), mN );
       SearchNearestInBox( *ThisPoint, Nearest, Distance, Box, Found );
       if(Found)
         return Nearest;
       return this->NullPointer();
     }

     //************************************************************************
     //************************************************************************

	 PointerType SearchNearestPoint( PointType const& ThisPoint )
	 {
	   PointerType Result            = *mPointBegin;                           //static_cast<PointerType>(NULL);
	   CoordinateType ResultDistance = static_cast<CoordinateType>(DBL_MAX);
       SearchStructureType Box( CalculateCell(ThisPoint), mN, mIndexCellBegin );
	   SearchNearestPointLocal( ThisPoint, Result, ResultDistance, Box );
	   return Result;
	 }

	 //************************************************************************

	 PointerType SearchNearestPoint( PointType const& ThisPoint, CoordinateType& rResultDistance )
	 {
	   PointerType Result = *mPointBegin;                           //static_cast<PointerType>(NULL);
	   rResultDistance    = static_cast<CoordinateType>(DBL_MAX);
       SearchStructureType Box( CalculateCell(ThisPoint), mN, mIndexCellBegin );
	   SearchNearestPointLocal( ThisPoint, Result, rResultDistance, Box);
	   return Result;
	 }

	 //************************************************************************
	 
     // New Thread Safe!!!
	 PointerType SearchNearestPoint( PointType const& ThisPoint, CoordinateType& rResultDistance, SearchStructureType& Box )
	 {
	   PointerType Result            = *mPointBegin;                           //static_cast<PointerType>(NULL);
	   CoordinateType ResultDistance = static_cast<CoordinateType>(DBL_MAX);
       Box.Set( CalculateCell(ThisPoint), mN, mIndexCellBegin );
	   SearchNearestPointLocal( ThisPoint, Result, rResultDistance, Box);
	   return Result;
	 }
     
	 //************************************************************************
	 //************************************************************************

     void SearchNearestPoint( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance )
     {
       SearchStructureType Box( CalculateCell(ThisPoint), mN, mIndexCellBegin );
       SearchNearestPointLocal(ThisPoint,rResult,rResultDistance,Box);
     }

     //************************************************************************

	 void SearchNearestPoint( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance, SearchStructureType& Box )
	 {
	   // This case is when BinStatic is a LeafType in Other Spacial Structure
	   // Then, it is possible a better Result before this search
       Box.Set( CalculateCell(ThisPoint), mN, mIndexCellBegin );
	   SearchNearestPointLocal( ThisPoint, rResult, rResultDistance, Box );
	 }

	 //************************************************************************
	 //************************************************************************

     // **** THREAD SAFE  -> The user pass the SearchStructure (BinBox)
	 void SearchNearestPointLocal( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance, SearchStructureType& Box ) 
	 {
	   if( mPointBegin == mPointEnd )
		 return;

       bool Found;

	   // initial search
	   ++Box;
	   SearchNearestInBox( ThisPoint, rResult, rResultDistance, Box, Found );
	   // increase mBox and try again
	   while(!Found){
		 ++Box;
		 SearchNearestInBox( ThisPoint, rResult, rResultDistance, Box, Found );
	   }
	 
     }

	 //************************************************************************
	 //************************************************************************

	 SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, IteratorType Results,
         DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults )
	 {	   
	   CoordinateType Radius2 = Radius * Radius;
	   SizeType NumberOfResults = 0;
	   SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   // CALCULATE BOUNING BOX WITH RADIUS
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
	   return NumberOfResults;
	 }

	 //************************************************************************
	 
     SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, IteratorType Results,
         DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
	 {
	   CoordinateType Radius2 = Radius * Radius;
       SizeType NumberOfResults = 0;
	   Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
       return NumberOfResults;
	 }

     //************************************************************************

	 void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults )
	 {
       SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
	 }

     //************************************************************************
	 
     void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
	 {
       Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
	 }
     
     
     //************************************************************************

     // **** THREAD SAFE
	
     // Dimension = 1
     void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
	 {
       SearchRadiusInRange()(*(Box.RowBegin),*(Box.RowEnd),ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
	 }
	
     // Dimension = 2
     void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
	 {
       for(IndexType I = Box.Axis[1].Begin() ; I <= Box.Axis[1].End() ; I += Box.Axis[1].Block )
         SearchRadiusInRange()(Box.RowBegin[I],Box.RowEnd[I],ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
	 }
	
     // Dimension = 3
     void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
	 {
	      static int index = 0;
	      
	      pointsToSearch[index].x = ThisPoint[0];
	      pointsToSearch[index].y = ThisPoint[1];
	      pointsToSearch[index].z = ThisPoint[2];
	      
	      index++;
	      
	      for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
		  for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
		      SearchRadiusInRange()(Box.RowBegin[I],Box.RowEnd[I],ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
	    
		  /*
	      static double total = 0;
	      double results = -1;
	      
	      //cl_double4 resultPoints[2000];
	      
	      
	      //Read number of results
	      Err = clEnqueueReadBuffer(OCLDeviceGroup->CommandQueues[0], OCL_results, CL_TRUE, 0, sizeof(double), &results, 0, NULL, NULL);
	      
	      //Read result points
	      //Err = clEnqueueReadBuffer(OCLDeviceGroup->CommandQueues[0], OCL_outData, CL_TRUE, 0, sizeof(cl_double4) * 20, &resultPoints, 0, NULL, NULL);
	     
	      //NumberOfResults = results;
	      total += results;
	      */
	      //Check if results are correct
	      /*if (NumberOfResults - results != 0 && false) {
		std::cout << "-----------------------------------------------" << std::endl;
		std::cout << "POINT: \t" << OCL_memPoint.x << " " << OCL_memPoint.y << " "  <<  OCL_memPoint.z << std::endl;
		
		total += results;
		if (NumberOfResults - results != 0)
		  std::cout << "*********************************************" << std::endl;
		std::cout << "DEF results: \t" << NumberOfResults << std::endl;
		//std::cout << Results << std::endl;
		std::cout << "OCL results: \t" << results << "\t total: \t" << total << std::endl;
		for( int i = 0; i < results; i++)
		{
		  std::cout << resultPoints[i].x << " " <<
				resultPoints[i].y << " " <<
				resultPoints[i].z << " " << " \t" << resultPoints[i].w <<
				
		  std::endl;
		}
		if (NumberOfResults - results != 0)
		  std::cout << "*********************************************" << std::endl;
	      }*/
	      
	      //Release execution-specific kernel buffers
	      /*clReleaseMemObject(OCL_Radius);
	      clReleaseMemObject(OCL_Point);
	      clReleaseMemObject(OCL_outData);
	      clReleaseMemObject(OCL_results);*/
	      
	      //NumberOfResults = results;
	      //std::cout << "DEF results \t" << NumberOfResults << "\t" << "OCL results: \t" << results << "\t total: \t" << total << std::endl;
	  }
         
         void computeresultsN(double Radius, int ConcurrentPoints) 
         {
	      cl_mem OCL_Radius;
	      cl_mem OCL_Radius2;
	      cl_mem OCL_Points;
	      cl_mem OCL_outData;
	      cl_mem OCL_results;
	      cl_mem OCL_resultsNum;
	      cl_mem OCL_w_size;
	      cl_mem OCL_maxResults;
	      
	      double total = 0;
	      
	      int pointSize = (mPointEnd - mPointBegin);
	      int processed = 0;
	      
	      int maxResults = 2;
	      int result = 0;
	      
	      double HOST_memRadius  = Radius;
	      double HOST_memRadius2 = Radius * Radius;
	      
	      while (processed < pointSize)
	      {
		  int amount;
		 
		  if (ConcurrentPoints > pointSize)
		    amount = pointSize;
		  else
		    amount = (processed + ConcurrentPoints) < pointSize ? ConcurrentPoints : pointSize - processed;
		  
		  std::cout << "amount " << amount << std::endl;
		  
		  int HOST_w_size = amount;
		  int results[amount];
		  
		  cl_double4 resultPoints[maxResults*amount]; // 20 resultados por punto
		
		  OCL_Radius  = clCreateBuffer(OCLDeviceGroup->Contexts[0], CL_MEM_READ_ONLY, sizeof(double), NULL, &Err);
		  OCL_Radius2 = clCreateBuffer(OCLDeviceGroup->Contexts[0], CL_MEM_READ_ONLY, sizeof(double), NULL, &Err);
		  
		  OCL_w_size  = clCreateBuffer(OCLDeviceGroup->Contexts[0]    , CL_MEM_READ_ONLY, sizeof(double), NULL, &Err);
		  OCL_maxResults  = clCreateBuffer(OCLDeviceGroup->Contexts[0], CL_MEM_READ_ONLY, sizeof(double), NULL, &Err);
		  
		  OCL_Points  = clCreateBuffer(OCLDeviceGroup->Contexts[0], CL_MEM_READ_ONLY, sizeof(cl_double4) * amount, NULL, &Err);

		  OCL_results = clCreateBuffer(OCLDeviceGroup->Contexts[0], CL_MEM_WRITE_ONLY, sizeof(int) * amount , NULL, &Err);
		  OCL_outData = clCreateBuffer(OCLDeviceGroup->Contexts[0], CL_MEM_WRITE_ONLY, sizeof(cl_double4) * amount * maxResults, NULL, &Err);

		  
		  Err = clEnqueueWriteBuffer(OCLDeviceGroup->CommandQueues[0], OCL_Radius , CL_TRUE, 0, sizeof(double),  &HOST_memRadius, 0, NULL, NULL);
		  KRATOS_OCL_CHECK(Err);

		  Err = clEnqueueWriteBuffer(OCLDeviceGroup->CommandQueues[0], OCL_Radius2, CL_TRUE, 0, sizeof(double),  &HOST_memRadius2, 0, NULL, NULL);
		  KRATOS_OCL_CHECK(Err);
		  
		  Err = clEnqueueWriteBuffer(OCLDeviceGroup->CommandQueues[0], OCL_w_size , CL_TRUE, 0, sizeof(int),     &HOST_w_size, 0, NULL, NULL);
		  KRATOS_OCL_CHECK(Err);
	
		  Err = clEnqueueWriteBuffer(OCLDeviceGroup->CommandQueues[0], OCL_Points , CL_TRUE, 0, sizeof(cl_double4) * amount, &pointsToSearch[processed],  0, NULL, NULL);
		  KRATOS_OCL_CHECK(Err);
		  
		  Err = clEnqueueWriteBuffer(OCLDeviceGroup->CommandQueues[0], OCL_maxResults , CL_TRUE, 0, sizeof(int),  &maxResults, 0, NULL, NULL);
		  KRATOS_OCL_CHECK(Err);
		  
		  OCLDeviceGroup->SetKernelArg(OCLSearchInRadius, 0,  OCL_IndexCellReference);
		  OCLDeviceGroup->SetKernelArg(OCLSearchInRadius, 1,  OCL_BinsContainer);
		
		  OCLDeviceGroup->SetKernelArg(OCLSearchInRadius, 2 , OCL_InvCellSize);
		  OCLDeviceGroup->SetKernelArg(OCLSearchInRadius, 3 , OCL_N);
		  OCLDeviceGroup->SetKernelArg(OCLSearchInRadius, 4,  OCL_Radius);
		  OCLDeviceGroup->SetKernelArg(OCLSearchInRadius, 5,  OCL_Radius2);
		  OCLDeviceGroup->SetKernelArg(OCLSearchInRadius, 6,  OCL_Points);
		  OCLDeviceGroup->SetKernelArg(OCLSearchInRadius, 7,  OCL_MinPoint);
		  OCLDeviceGroup->SetKernelArg(OCLSearchInRadius, 8,  OCL_outData);
		  OCLDeviceGroup->SetKernelArg(OCLSearchInRadius, 9, OCL_results);
		  OCLDeviceGroup->SetKernelArg(OCLSearchInRadius, 10, OCL_w_size);
		  OCLDeviceGroup->SetKernelArg(OCLSearchInRadius, 11, OCL_maxResults);
		  
		  //Execute kernel
		  OCLDeviceGroup->ExecuteKernel(OCLSearchInRadius, amount);
		  
		  //Read number of results and results
		  Err = clEnqueueReadBuffer(OCLDeviceGroup->CommandQueues[0], OCL_results, CL_TRUE, 0, sizeof(int) * amount, &results, 0, NULL, NULL);
		  Err = clEnqueueReadBuffer(OCLDeviceGroup->CommandQueues[0], OCL_outData, CL_TRUE, 0, sizeof(cl_double4) * amount * maxResults, &resultPoints, 0, NULL, NULL);
		  
		  for(int j = 0; j < amount; j++)
		      result += results[j];
		  
		  clReleaseMemObject(OCL_Radius);
		  clReleaseMemObject(OCL_Radius2);
		  clReleaseMemObject(OCL_Points);
		  clReleaseMemObject(OCL_outData);
		  clReleaseMemObject(OCL_results);
		  clReleaseMemObject(OCL_resultsNum);
		  clReleaseMemObject(OCL_w_size);
		  clReleaseMemObject(OCL_maxResults);
		  
		  processed += amount;
		  std::cout << "Total Results while " << processed << " points done: " << result << std::endl;
	      }
	      
	      std::cout << "Total Results: " << result << " (" << maxResults * (mPointEnd - mPointBegin) << " stored)"<< std::endl;
         }

	 //************************************************************************
	 //************************************************************************

	 SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType Radius, IteratorType Results, SizeType MaxNumberOfResults )
	 {
	   CoordinateType Radius2 = Radius * Radius;
       SizeType NumberOfResults = 0;
       SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
       return NumberOfResults;
	 }

	 //************************************************************************

	 SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType Radius, IteratorType Results,
         SizeType MaxNumberOfResults, SearchStructureType& Box )
	 {
	   CoordinateType Radius2 = Radius * Radius;
       SizeType NumberOfResults = 0;
	   Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
       return NumberOfResults;
	 }

	 //************************************************************************

	 void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults )
	 {
       SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
	 }

	 //************************************************************************

	 void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
	 {
       Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN, mIndexCellBegin );
	   SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
	 }

	 //************************************************************************

     // Dimension = 1
	 void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, 
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
	 {
       SearchRadiusInRange()(*(Box.RowBegin),*(Box.RowEnd),ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
	 }
	 
     // Dimension = 2
     void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, 
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
	 {
       for(IndexType I = Box.Axis[1].Begin() ; I <= Box.Axis[1].End() ; I += Box.Axis[1].Block )
         SearchRadiusInRange()(Box.RowBegin[I],Box.RowEnd[I],ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
	 }
     
     // Dimension = 3
     void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, 
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
	 {
       for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
         for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
           SearchRadiusInRange()(Box.RowBegin[I],Box.RowEnd[I],ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
	 }

	 //************************************************************************
	 //************************************************************************

     // Dimension = 1
	 void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance, 
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, bool& Found )
     {
       Found = false;
       SearchNearestInRange()( *(Box.RowBegin), *(Box.RowEnd), ThisPoint, ResultPoint, ResultDistance, Found );
     }

     // Dimension = 2
	 void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance, 
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, bool& Found )
     {
       Found = false;
       for(IndexType I = Box.Axis[1].Begin() ; I <= Box.Axis[1].End() ; I += Box.Axis[1].Block )
         SearchNearestInRange()( Box.RowBegin[I], Box.RowEnd[I], ThisPoint, ResultPoint, ResultDistance, Found );
     }
	 
     // Dimension = 3
     void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance, 
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, bool& Found )
     {
       Found = false;
       for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
         for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
           SearchNearestInRange()( Box.RowBegin[I], Box.RowEnd[I], ThisPoint, ResultPoint, ResultDistance, Found );
     }

	 //************************************************************************
	 //************************************************************************
	 
	 SizeType SearchInBox( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType Results,
         SizeType MaxNumberOfResults )
     {
       SizeType NumberOfResults = 0;
       SearchStructureType Box( CalculateCell(SearchMinPoint), CalculateCell(SearchMaxPoint), mN, mIndexCellBegin );
       SearchInBoxLocal( SearchMinPoint, SearchMaxPoint, Results, NumberOfResults, MaxNumberOfResults, Box );
       return NumberOfResults;
     }
     
     //************************************************************************
	 
     void SearchInBox(PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& Results, SizeType& NumberOfResults,
            SizeType const& MaxNumberOfResults )
     {
       NumberOfResults = 0;
       SearchStructureType Box( CalculateCell(SearchMinPoint), CalculateCell(SearchMaxPoint), mN, mIndexCellBegin );
       SearchInBoxLocal( SearchMinPoint, SearchMaxPoint, Results, NumberOfResults, MaxNumberOfResults, Box );
     }
     
     //************************************************************************

     // Dimension = 1
	 void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
     {
       SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,*(Box.RowBegin),*(Box.RowEnd),ResultsPoint,NumberOfResults,MaxNumberOfResults);
     }

     // Dimension = 2
	 void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
     {
       for(IndexType I = Box.Axis[1].Begin() ; I <= Box.Axis[1].End() ; I += Box.Axis[1].Block )
         SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,Box.RowBegin[I],Box.RowEnd[I],ResultsPoint,NumberOfResults,MaxNumberOfResults);
     }

     // Dimension = 3
	 void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
         SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
         SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
     {
       for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
         for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
           SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,Box.RowBegin[I],Box.RowEnd[I],ResultsPoint,NumberOfResults,MaxNumberOfResults);
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

     // Point Access Iterators (vector reordered!!)
	 IteratorType     mPointBegin;
	 IteratorType     mPointEnd;

     // Bin Parameters (Sizes,BoundingBox,...)
	 PointType  mMinPoint;
	 PointType  mMaxPoint;
	 PointType  mCellSize;
	 PointType  mInvCellSize;
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
