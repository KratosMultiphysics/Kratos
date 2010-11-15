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
      cl_uint OCL_Points, OCL_Cell, OCL_IndexCell, OCL_IndexCellReference, OCL_InvCellSize, OCL_N, OCL_MinPoint, OCL_BinsContainer;
      
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
	   
	   pointsToSearch = new cl_double4[SearchUtils::PointerDistance(mPointBegin, mPointEnd)];
	   
	   std::cout << "Generate Params..." << std::endl;
	   
	   char params[512] = "";
	   sprintf(params,"-D POINT_SIZE=%ld -D T_DIMENSION=%Zu -D CELL_SIZE=%d",SearchUtils::PointerDistance(mPointBegin, mPointEnd),TDimension,Size+1);
	   
	   std::cout << "Building programs" << std::endl;	   
	   cl_uint OCL_program = OCLDeviceGroup->BuildProgramFromFile("binshashOPT.cl",params);
	   
	   std::cout << "Register Kernels..." << std::endl;	   
	   OCLGenerateBins   = OCLDeviceGroup->RegisterKernel(OCL_program,"GenerateBins");
	   OCLSearchInRadius = OCLDeviceGroup->RegisterKernel(OCL_program,"SearchInRadiusMultiple");
	   
	   int j = 0;
	   
	   cl_double4 * points = new cl_double4[SearchUtils::PointerDistance(mPointBegin, mPointEnd)];
	   cl_double4 * MinPoint = new cl_double4();
	   
	   double InvCellSize[TDimension];
	   double N[TDimension];

	   for(IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++) {
		  points[j].x = (**Point)[0];
		  points[j].y = (**Point)[1];
		  points[j].z = (**Point)[2];
		  points[j].w = j;
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
	   
	   // Prepare buffers
	   OCL_Points             = OCLDeviceGroup->CreateBuffer(SearchUtils::PointerDistance(mPointBegin, mPointEnd) * sizeof(cl_double4), CL_MEM_READ_ONLY);
	   OCL_InvCellSize        = OCLDeviceGroup->CreateBuffer(TDimension * sizeof(double), CL_MEM_READ_ONLY);
	   OCL_N                  = OCLDeviceGroup->CreateBuffer(TDimension * sizeof(double), CL_MEM_READ_ONLY);
	   OCL_MinPoint           = OCLDeviceGroup->CreateBuffer(sizeof(cl_double4), CL_MEM_READ_ONLY);   
	   OCL_Cell               = OCLDeviceGroup->CreateBuffer((Size+1) * sizeof(cl_double4), CL_MEM_READ_WRITE);
	   OCL_IndexCellReference = OCLDeviceGroup->CreateBuffer((Size+1) * sizeof(int), CL_MEM_READ_WRITE);   
	   OCL_BinsContainer      = OCLDeviceGroup->CreateBuffer(SearchUtils::PointerDistance(mPointBegin, mPointEnd) * sizeof(cl_double4), CL_MEM_READ_WRITE);

	   // Load data into the input buffer
	   std::cout << "Loading data into inputbuffer..." << std::endl;
	   OCLDeviceGroup->CopyBuffer(OCL_Points     , OpenCL::HostToDevice, OpenCL::VoidPList(1,points));
	   OCLDeviceGroup->CopyBuffer(OCL_InvCellSize, OpenCL::HostToDevice, OpenCL::VoidPList(1,InvCellSize));
	   OCLDeviceGroup->CopyBuffer(OCL_N          , OpenCL::HostToDevice, OpenCL::VoidPList(1,N));
	   OCLDeviceGroup->CopyBuffer(OCL_MinPoint   , OpenCL::HostToDevice, OpenCL::VoidPList(1,MinPoint));
	   
	   // Set arguments
	   std::cout << "Setting arguments..." << std::endl;
	   OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBins, 0,  OCL_Points);
	   OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBins, 1,  OCL_IndexCellReference);
	   OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBins, 2,  OCL_InvCellSize);
	   
	   OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBins, 3,  OCL_N);
	   OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBins, 4,  OCL_MinPoint);
	   OCLDeviceGroup->SetBufferAsKernelArg(OCLGenerateBins, 5,  OCL_BinsContainer);
	   
	   std::cout << "Executing kernel..." << std::endl;
	   OCLDeviceGroup->ExecuteKernel(OCLGenerateBins, Size+1);
	 }

	 //************************************************************************

	 void GenerateBins( ){
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
     
     void prepareData(PointType const& ThisPoint) {
	  static int index = 0;
	  
	  pointsToSearch[index].x = ThisPoint[0];
	  pointsToSearch[index].y = ThisPoint[1];
	  pointsToSearch[index].z = ThisPoint[2];
	  
	  index++;
	  if (index == SearchUtils::PointerDistance(mPointBegin, mPointEnd)) index = 0;
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
	      for(IndexType II = Box.Axis[2].Begin() ; II <= Box.Axis[2].End() ; II += Box.Axis[2].Block )
		  for(IndexType I = II + Box.Axis[1].Begin() ; I <= II + Box.Axis[1].End() ; I += Box.Axis[1].Block )
		      SearchRadiusInRange()(Box.RowBegin[I],Box.RowEnd[I],ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
	  }
         
	  void computeresultsN(double Radius, int ConcurrentPoints, int maxResults) 
	  {
	      cl_uint OCL_Radius, OCL_Radius2, OCL_Points, OCL_outData, OCL_results, OCL_resultsNum, OCL_w_size, OCL_maxResults;
	      
	      int pointSize = SearchUtils::PointerDistance(mPointBegin, mPointEnd);
	      int processed = 0;
	      
	      int result = 0;
	      
	      double HOST_memRadius  = Radius;
	      double HOST_memRadius2 = Radius * Radius;
	      
	      int amount;
	      
	      OCL_Radius     = OCLDeviceGroup->CreateBuffer(sizeof(double), CL_MEM_READ_ONLY);
	      OCL_Radius2    = OCLDeviceGroup->CreateBuffer(sizeof(double), CL_MEM_READ_ONLY);
	      OCL_w_size     = OCLDeviceGroup->CreateBuffer(sizeof(double), CL_MEM_READ_ONLY);
	      OCL_maxResults = OCLDeviceGroup->CreateBuffer(sizeof(double), CL_MEM_READ_ONLY);
	      OCL_Points     = OCLDeviceGroup->CreateBuffer(sizeof(cl_double4) * ConcurrentPoints, CL_MEM_READ_ONLY);
	      OCL_results    = OCLDeviceGroup->CreateBuffer(sizeof(int) * ConcurrentPoints, CL_MEM_WRITE_ONLY);
	      OCL_outData    = OCLDeviceGroup->CreateBuffer(sizeof(int) * ConcurrentPoints * maxResults, CL_MEM_WRITE_ONLY);
	  
	      while (processed < pointSize)
	      {	  
		  amount = (ConcurrentPoints > pointSize) ? pointSize : (processed + ConcurrentPoints) < pointSize ? ConcurrentPoints : pointSize - processed;
		  
		  int HOST_w_size = amount;
		  int * results;
		  int * resultPoints;
		  
		  results = (int *)malloc(sizeof(int) * amount);
		  resultPoints = (int *)malloc(sizeof(int) * amount * maxResults);

		  OCLDeviceGroup->CopyBuffer(OCL_Radius     , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius));
		  OCLDeviceGroup->CopyBuffer(OCL_Radius2    , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_memRadius2));
		  OCLDeviceGroup->CopyBuffer(OCL_w_size     , OpenCL::HostToDevice, OpenCL::VoidPList(1,&HOST_w_size));
		  OCLDeviceGroup->CopyBuffer(OCL_Points     , OpenCL::HostToDevice, OpenCL::VoidPList(1,&pointsToSearch[processed]));
		  OCLDeviceGroup->CopyBuffer(OCL_maxResults , OpenCL::HostToDevice, OpenCL::VoidPList(1,&maxResults));
		  		  
		  OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 0,  OCL_IndexCellReference);
		  OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 1,  OCL_BinsContainer);
		  OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 2 , OCL_InvCellSize);
		  OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 3 , OCL_N);
		  OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 4,  OCL_Radius);
		  OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 5,  OCL_Radius2);
		  OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 6,  OCL_Points);
		  OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 7,  OCL_MinPoint);
		  OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 8,  OCL_outData);
		  OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 9,  OCL_results);
		  OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 10, OCL_w_size);
		  OCLDeviceGroup->SetBufferAsKernelArg(OCLSearchInRadius, 11, OCL_maxResults);
		  
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
