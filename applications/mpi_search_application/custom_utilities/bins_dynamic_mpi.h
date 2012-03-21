//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2007-03-27 17:02:19 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_BINS_DYNAMIC_MPI_CONTAINER_H_INCLUDE)
#define KRATOS_BINS_DYNAMIC_MPI_CONTAINER_H_INCLUDE

#include "mpi.h"
#include "spatial_containers/tree.h"
#include "containers/buffer.h"
#include "includes/serializer.h"

#define KRATOS_FILL_SERIALIZATOR_BUFFER(size,n,source,dest) \
for(int i = 0; i < size; i++) (*dest) << source[size*n+i]

namespace Kratos {

  template<
      std::size_t TDimension,
      class TPointType,
      class TContainerType,
      class TPointerType = typename TContainerType::value_type,
      class TIteratorType = typename TContainerType::iterator,
      class TDistanceIteratorType = typename std::vector<double>::iterator,
      class TDistanceFunction = Kratos::SearchUtils::SquaredDistanceFunction<TDimension,TPointType>
      >
   class BinsDynamicMpi : public TreeNode<TDimension,TPointType, TPointerType, TIteratorType, TDistanceIteratorType, typename std::vector<TPointerType>::iterator >
	  {
		
	public:

        /// Pointer definition of BinsDynamicMpi
        KRATOS_CLASS_POINTER_DEFINITION(BinsDynamicMpi);

	typedef TreeNode<TDimension,TPointType,TPointerType,TIteratorType,TDistanceIteratorType> TreeNodeType;
        typedef TPointType                         PointType;
	typedef TContainerType                     ContainerType;
	typedef TIteratorType                      IteratorType;
	typedef TDistanceIteratorType              DistanceIteratorType;
	typedef TPointerType                       PointerType;
	typedef TDistanceFunction                  DistanceFunction;
	
        enum { Dimension = TDimension };

	typedef typename TreeNodeType::CoordinateType  CoordinateType;  // double
	typedef typename TreeNodeType::SizeType        SizeType;        // std::size_t
	typedef typename TreeNodeType::IndexType       IndexType;       // std::size_t

	typedef TreeNodeType LeafType;
    
        typedef typename TreeNodeType::IteratorIteratorType IteratorIteratorType;
        typedef typename TreeNodeType::SearchStructureType SearchStructureType;


	// Local Container ( PointPointer Container per Cell )
	// can be different to ContainerType
        // not always PointVector == ContainerType ( if ContainerType = C array )
	typedef std::vector<PointerType>       PointVector;
	typedef typename PointVector::iterator PointIterator;
		
        // Global Container
	typedef std::vector<PointVector>        CellsContainerType;
        //typedef typename CellsContainerType::iterator IteratorIteratorType;

	typedef Tvector<IndexType,TDimension>   CellType;
        
        typedef Kratos::SearchUtils::SearchNearestInRange<PointType,PointerType,PointIterator,DistanceFunction,CoordinateType> SearchNearestInRange;
        typedef Kratos::SearchUtils::SearchRadiusInRange<PointType,PointIterator,DistanceIteratorType,DistanceFunction,SizeType,CoordinateType,IteratorType> SearchRadiusInRange;
	typedef Kratos::SearchUtils::SearchBoxInRange<PointType,PointIterator,SizeType,TDimension,IteratorType> SearchBoxInRange;


	public:

	//************************************************************************

	// constructor 1
	BinsDynamicMpi() : mPointBegin(this->NullIterator()), mPointEnd(this->NullIterator()), mNumPoints(0)
	{};

	///UNSUPPORTED CLASS CONSTRUCTOR
        //************************************************************************
         
//         BinsDynamicMpi( IteratorType const& PointBegin, IteratorType const& PointEnd, SizeType BucketSize = 1 )
//         : mPointBegin(PointBegin), mPointEnd(PointEnd)
//         {
//            if(mPointBegin==mPointEnd)
//               return;
//            mNumPoints = std::distance(mPointBegin,mPointEnd);
//            CalculateBoundingBox();
//            CalculateCellSize();
//            AllocateCellsContainer();
//            GenerateBins();
//         }
//         
//         //************************************************************************
//          
//         BinsDynamicMpi( IteratorType const& PointBegin, IteratorType const& PointEnd, PointType const& MinPoint, PointType const& MaxPoint, SizeType BucketSize = 1 )
//         : mPointBegin(PointBegin), mPointEnd(PointEnd)
//         {
//            if(mPointBegin==mPointEnd)
//               return;
// 
//            mNumPoints = std::distance(mPointBegin,mPointEnd);
// 		   for(SizeType i = 0 ; i < TDimension ; i++)
// 		   {
// 			 mMinPoint[i] = MinPoint[i];
// 			 mMaxPoint[i] = MaxPoint[i];
// 		   }
//            CalculateCellSize();
//            AllocateCellsContainer();
//            GenerateBins();
//         }
//         
//         //************************************************************************
//          
//         BinsDynamicMpi( PointType const& MinPoint, PointType const& MaxPoint, SizeType BucketSize )
//           : mNumPoints(0)
//         {
// 	for(SizeType i = 0 ; i < TDimension ; i++)
// 		   {
// 			 mMinPoint[i] = MinPoint[i];
// 			 mMaxPoint[i] = MaxPoint[i];
// 		   }
//            CalculateCellSize(BucketSize);
//            AllocateCellsContainer();
//         }
//         
//         //************************************************************************
//         
//         BinsDynamicMpi( IteratorType const& PointBegin, IteratorType const& PointEnd, CoordinateType BoxSize, SizeType BucketSize = 1 )
//         : mPointBegin(PointBegin), mPointEnd(PointEnd)
//         {
//            if(mPointBegin==mPointEnd)
//               return;
//            mNumPoints = std::distance(mPointBegin,mPointEnd);
//            CalculateBoundingBox();
//            CalculateCellSize(BoxSize);
//            AllocateCellsContainer();
//            GenerateBins();
//         }
        
        ///Constructor - MPI
        //************************************************************************
        
        BinsDynamicMpi( ModelPart * StaticMesh, ModelPart * ParticMesh, CoordinateType BoxSize, SizeType BucketSize = 1 )
//         : mPointBegin(StaticMesh.GetCommunicator().LocalMesh().NodesBegin()), mPointEnd(StaticMesh.GetCommunicator().LocalMesh().NodesEnd())
        {
	    //Initialize standard dynamic bins structure
	    IteratorType mPointIterator;
	   
	    mPointBegin = new PointType* [StaticMesh->GetCommunicator().LocalMesh().NumberOfNodes()];
	    mPointEnd = mPointBegin + StaticMesh->GetCommunicator().LocalMesh().NumberOfNodes();
	    
	    mPointIterator = mPointBegin;
	   
	    std::cout << "Parsing local mesh elements" << std::endl;
	    for( ModelPart::NodesContainerType::iterator inode = StaticMesh->GetCommunicator().LocalMesh().NodesBegin(); inode != StaticMesh->GetCommunicator().LocalMesh().NodesEnd(); inode++, mPointIterator++)
	    { 
		PointType auxPoint;

		auxPoint[0] = inode->X();
		auxPoint[1] = inode->Y();
		auxPoint[2] = inode->Z();

		(*mPointIterator) = new PointType(auxPoint);
	    }
	    
	    if(mPointBegin==mPointEnd)
		return;
	   
	    mNumPoints = std::distance(mPointBegin,mPointEnd);
	    CalculateBoundingBox();
	    CalculateCellSize(BoxSize);
	    AllocateCellsContainer();
	    GenerateBins();
	   
	    //Set up MPI interface for dynamic bins
	    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	    
	    std::cout << "Initialize done" << std::endl;
        }
        
        //************************************************************************

        // destructor
        virtual ~BinsDynamicMpi(){ }
        
        //************************************************************************
        
        IteratorType Begin() { return mPointBegin; }

        //************************************************************************
        
        IteratorType End() { return mPointBegin; }

        //************************************************************************
        
        CoordinateType CellSize( SizeType const& iDim ) { return mCellSize[iDim]; }
        
        //************************************************************************
        
        SizeType NumCell( SizeType const& iDim ) { return mN[iDim]; }
        
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
        }
        
        //************************************************************************
        
        void CalculateCellSize() 
        {

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
            SizeType Size = 1;
            for(SizeType i = 0 ; i < TDimension ; i++)
                Size *= mN[i];
            // Resize Global Container
            mPoints.resize(Size);
        }
        
        //************************************************************************
         
        void GenerateBins(){
           
           for(IteratorType i_point = mPointBegin ; i_point != mPointEnd ; i_point++)
              mPoints[CalculateIndex(**i_point)].push_back(*i_point); 
            
        }
        
        //************************************************************************

        IndexType CalculatePosition( CoordinateType const& ThisCoord, SizeType ThisDimension )
        {
            CoordinateType d_index = (ThisCoord - mMinPoint[ThisDimension]) * mInvCellSize[ThisDimension];
            IndexType index = static_cast<IndexType>( (d_index < 0.00) ? 0.00 : d_index );
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
        
         void AddPoint( PointerType const& ThisPoint ){
            mPoints[CalculateIndex(*ThisPoint)].push_back(ThisPoint);
            mNumPoints++;
         }
         
        //************************************************************************
         
		void MPI_ExistPoint( PointerType const& ThisPoint, PointerType ResultNearest, CoordinateType const Tolerance = static_cast<CoordinateType>(10.0*DBL_EPSILON) )
		{
			PointerType Nearest, remoteNearest[mpi_size], resultNearest[mpi_size], remoteThisPoint[mpi_size];
			CoordinateType Distance, remoteDistance[mpi_size], resultDistance[mpi_size];
			bool Found, remoteFound[mpi_size], resultFound[mpi_size];
			
			int msgSendSize = 0;
			int msgRecvSize = 0;
			
			int msgResSendSize = 0;
			int msgResRecvSize = 0;
			
			std::cout << "(" << mpi_rank << ") --- " << (*ThisPoint) << " --- " << std::endl;
			
			Serializer particleSerializer;
			particleSerializer.save("nodes",ThisPoint);  

			std::stringstream* serializer_buffer;
			
			serializer_buffer = (std::stringstream *)particleSerializer.pGetBuffer();
			msgSendSize = serializer_buffer->str().size();
			
			MPI_Allreduce(&msgSendSize,&msgRecvSize,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
			
			char mpi_send_buffer[(msgRecvSize+1)];
			char mpi_recv_buffer[(msgRecvSize+1) * mpi_size];
			
			strcpy (mpi_send_buffer, serializer_buffer->str().c_str());
			mpi_send_buffer[msgSendSize] = '\0';

			MPI_Allgather(mpi_send_buffer,(msgRecvSize+1),MPI_CHAR,mpi_recv_buffer,(msgRecvSize+1),MPI_CHAR,MPI_COMM_WORLD);
			
			for(int i = 0; i < mpi_size; i++) 
			{
				Serializer recvParticleSerializer;
			 	serializer_buffer = (std::stringstream *)recvParticleSerializer.pGetBuffer();
		
				for(int j = 0; mpi_recv_buffer[(msgRecvSize+1)*i+j] != '\0'; j++) 
				{
					(*serializer_buffer) << mpi_recv_buffer[(msgRecvSize+1)*i+j];
			  	}
				  
			  	remoteThisPoint[i]	= new PointType();
			  	remoteNearest[i]   	= new PointType();
			  	remoteDistance[i] 	= static_cast<CoordinateType>(DBL_MAX);
				  
				recvParticleSerializer.load("nodes",remoteThisPoint[i]);
				  
				std::cout << "(" << mpi_rank << ")" << " Restored Par: " << "(" << remoteThisPoint[i]->X() << " " << remoteThisPoint[i]->Y() << " " << remoteThisPoint[i]->Z() << ")" << std::endl;
				  
				SearchStructureType remote_Box( CalculateCell(*remoteThisPoint[i],-Tolerance), CalculateCell(*remoteThisPoint[i],Tolerance), mN );
				SearchNearestInBox( *remoteThisPoint[i], remoteNearest[i], remoteDistance[i], remote_Box, remoteFound[i] );
				  
				std::cout << "(" << mpi_rank << ") Found point: (" << remoteThisPoint[i]->X() << " " << remoteThisPoint[i]->Y() << " " << remoteThisPoint[i]->Z() << ") from process(" << i << "): " << (*(remoteNearest[i])) << " with dist: " << remoteDistance[i] << std::endl;
			}

			std::stringstream * res_serializer_buffer[mpi_size];
			
			for(int i = 0; i < mpi_size; i++) 
			{
				Serializer resSerializer;
				resSerializer.save("nodes",remoteNearest[i]);
			
				res_serializer_buffer[i] = (std::stringstream *)resSerializer.pGetBuffer();
				msgResSendSize = res_serializer_buffer[i]->str().size();
		
				msgResSendSize = msgResSendSize > msgResRecvSize ? msgResSendSize : msgResRecvSize;
		
				MPI_Allreduce(&msgResSendSize,&msgResRecvSize,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
			}
			
			char mpi_res_send_buffer[((msgResRecvSize + 1) * mpi_size)];
			char mpi_res_recv_buffer[((msgResRecvSize + 1) * mpi_size)];
			
			for(int i = 0; i < mpi_size; i++) 
			{
				strcpy(&mpi_res_send_buffer[(msgResRecvSize + 1) * i], res_serializer_buffer[i]->str().c_str());
				mpi_res_send_buffer[(msgResRecvSize + 1) * i + res_serializer_buffer[i]->str().size()] = '\0';
			}

			MPI_Alltoall(mpi_res_send_buffer,(msgResRecvSize+1),MPI_CHAR,mpi_res_recv_buffer,(msgResRecvSize+1),MPI_CHAR,MPI_COMM_WORLD);
			MPI_Alltoall(remoteDistance,1,MPI_DOUBLE,resultDistance,1,MPI_DOUBLE,MPI_COMM_WORLD);
			MPI_Alltoall(remoteFound,1,MPI_BYTE,resultFound,1,MPI_BYTE,MPI_COMM_WORLD);
			
			for (int i = 0; i < mpi_size; i++)
			{
				Serializer recvResParticleSerializer;
				serializer_buffer = (std::stringstream *)recvResParticleSerializer.pGetBuffer();
		
				for(int j = 0; mpi_res_recv_buffer[(msgResRecvSize+1)*i+j] != '\0'; j++) 
				{
			  		(*serializer_buffer) << mpi_res_recv_buffer[(msgResRecvSize+1)*i+j];
				}
		
				resultNearest[i] = new PointType();
				recvResParticleSerializer.load("nodes",resultNearest[i]);

				std::cout << "(" << mpi_rank << ") Result point from process (" << i << "): (" << resultNearest[i]->X() << " " << resultNearest[i]->Y() << " " << resultNearest[i]->Z() << ") with dist: " << resultDistance[i] << std::endl; 
			}

			Nearest 	= resultNearest[0];
			Distance 	= resultDistance[0];
			Found 		= resultFound[0];

			for(int i = 1; i < mpi_size; i++) 
			{
				if(resultFound[i] && resultDistance[i] < Distance) 
				{
					Nearest 	= resultNearest[0];
					Distance 	= resultDistance[0];
					Found 		= resultFound[0];
				}
			}
			
			ResultNearest = this->NullPointer();

	 	    if(Found)
	 	    	ResultNearest = Nearest;
		}
		
		void MPI_SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, IteratorType Results, 
             DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults, SizeType const& ResultsNumberOfResults )
        {
			CoordinateType Radius2 = Radius * Radius;
			SizeType NumberOfResults = 0;
			SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
			SearchInRadiusMpiWrapper( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
		   
//            return NumberOfResults;
        }

         //************************************************************************

//          void MPI_SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, IteratorType Results,
//              DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
//          {
//            CoordinateType Radius2 = Radius * Radius;
//            SizeType NumberOfResults = 0;
//            Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
//            SearchInRadiusMpiWrapper( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
// //            return NumberOfResults;
//          }
// 
//          //************************************************************************
// 
//          void MPI_SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
//              DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults )
//          {
//            SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
//            SearchInRadiusMpiWrapper( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
//          }
// 
//          //************************************************************************
// 
//          void MPI_SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
//              DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
//          {
//            Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
//            SearchInRadiusMpiWrapper( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
//          }

		void SearchInRadiusMpiWrapper( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
             DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,TDimension>& Box )
		{  
			PointType remoteThisPoint[mpi_size];
			PointerType remoteResults[mpi_size][MaxNumberOfResults], recvResults[mpi_size][MaxNumberOfResults];
			double remoteResultsDistances[mpi_size][MaxNumberOfResults], recvResultsDistances[mpi_size][MaxNumberOfResults];
			int messageSendNumberOfResults[mpi_size], messageRecvNumberOfResults[mpi_size];
			SizeType remoteNumberOfResults[mpi_size];

			int msgSendSize = 0;
			int msgRecvSize = 0;

			int msgResSendSize = 0;
			int msgResRecvSize = 0;

			std::cout << "(" << mpi_rank << ") --- " << ThisPoint << " --- " << std::endl;

			Serializer particleSerializer;
			particleSerializer.save("nodes",ThisPoint);  

			std::stringstream* serializer_buffer;

			serializer_buffer = (std::stringstream *)particleSerializer.pGetBuffer();
			msgSendSize = serializer_buffer->str().size();

			MPI_Allreduce(&msgSendSize,&msgRecvSize,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
			
			char mpi_send_buffer[(msgRecvSize+1)];
			char mpi_recv_buffer[(msgRecvSize+1) * mpi_size];

			strcpy (mpi_send_buffer, serializer_buffer->str().c_str());
			mpi_send_buffer[msgSendSize] = '\0';

			MPI_Allgather(mpi_send_buffer,(msgRecvSize+1),MPI_CHAR,mpi_recv_buffer,(msgRecvSize+1),MPI_CHAR,MPI_COMM_WORLD);

			for(int i = 0; i < mpi_size; i++) 
			{
				IteratorType remoteResultsPointer = &remoteResults[i][0];
				double * remoteResultsDistancesPointer = remoteResultsDistances[i];
				
				Serializer recvParticleSerializer;
				serializer_buffer = (std::stringstream *)recvParticleSerializer.pGetBuffer();

				for(int j = 0; mpi_recv_buffer[(msgRecvSize+1)*i+j] != '\0'; j++) 
				{
					(*serializer_buffer) << mpi_recv_buffer[(msgRecvSize+1)*i+j];
				}

				remoteNumberOfResults[i] = 0;
				
				std::cout << "Restoring Point" << std::endl;

				recvParticleSerializer.load("nodes",remoteThisPoint[i]);
					
				std::cout << "(" << mpi_rank << ")" << " Restored Par: " << "(" << remoteThisPoint[i].X() << " " << remoteThisPoint[i].Y() << " " << remoteThisPoint[i].Z() << ")" << std::endl;

				SearchInRadiusLocal(remoteThisPoint[i],Radius,Radius2,remoteResultsPointer,remoteResultsDistancesPointer,remoteNumberOfResults[i],MaxNumberOfResults,Box);

				std::cout << "(" << mpi_rank << ") Found points for: (" << remoteThisPoint[i].X() << " " << remoteThisPoint[i].Y() << " " << remoteThisPoint[i].Z() << ") from process(" << i << "): FOUND: " << remoteNumberOfResults[i] << std::endl;

				for(int j = 0; j < remoteNumberOfResults[i]; j++)
				{
					std::cout << "(" << mpi_rank << ")\t" << *(remoteResults[i][j]) << " " << remoteResultsDistances[i][j] << std::endl;	
				}
			}
			
			std::cout << "(" << mpi_rank << ") Send Back Results " << std::endl;
			
			//////////////////////////////////////////////////////////////////////////////////////////////
			// 	SEND BACK RESULTS 																		//
			//////////////////////////////////////////////////////////////////////////////////////////////

			std::stringstream * res_serializer_buffer[mpi_size][MaxNumberOfResults];
			std::string message[mpi_size][MaxNumberOfResults];
			int bufferSize[mpi_size][MaxNumberOfResults];
	
			for(int i = 0; i < mpi_size; i++) 
			{
				for(int j = 0; j < remoteNumberOfResults[i]; j++) 
				{
					Serializer resSerializer;
					resSerializer.save("nodes",remoteResults[i][j]);
					  
					res_serializer_buffer[i][j] = (std::stringstream *)resSerializer.pGetBuffer();
					message[i][j] = std::string(res_serializer_buffer[i][j]->str().c_str());
					bufferSize[i][j] = res_serializer_buffer[i][j]->str().size();
					
					msgResSendSize = msgResSendSize > bufferSize[i][j] ? msgResSendSize : bufferSize[i][j];
				}
			}
			
			MPI_Allreduce(&msgResSendSize,&msgResRecvSize,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
			
			std::cout << "(" << mpi_rank << ") Results Prepared " << std::endl;
			
			char mpi_res_send_buffer[((msgResRecvSize + 1) * mpi_size * MaxNumberOfResults)];
			char mpi_res_recv_buffer[((msgResRecvSize + 1) * mpi_size * MaxNumberOfResults)];
			
			std::cout << "(" << mpi_rank << ") Buffers allocated " << std::endl;
			
			for(int i = 0; i < mpi_size; i++) 
			{
				messageSendNumberOfResults[i] = remoteNumberOfResults[i];
				for(int j = 0; j < remoteNumberOfResults[i]; j++) 
				{
					strcpy(&mpi_res_send_buffer[(msgResRecvSize + 1)*MaxNumberOfResults*i+(msgResRecvSize + 1)*j], message[i][j].c_str());
					mpi_res_send_buffer[(msgResRecvSize + 1)*MaxNumberOfResults*i+(msgResRecvSize + 1)*j+bufferSize[i][j]] = '\0';
				}
			}
			
			std::cout << "(" << mpi_rank << ") Buffers Filled: " <<  ((msgResRecvSize + 1) * mpi_size * MaxNumberOfResults) * sizeof(char)<< std::endl;

			MPI_Alltoall(mpi_res_send_buffer,((msgResRecvSize+1) * MaxNumberOfResults),MPI_CHAR,mpi_res_recv_buffer,((msgResRecvSize+1) * MaxNumberOfResults),MPI_CHAR,MPI_COMM_WORLD);
			MPI_Alltoall(messageSendNumberOfResults,1,MPI_INT,messageRecvNumberOfResults,1,MPI_INT,MPI_COMM_WORLD);
			MPI_Alltoall(remoteResultsDistances[0],MaxNumberOfResults,MPI_DOUBLE,recvResultsDistances[0],MaxNumberOfResults,MPI_DOUBLE,MPI_COMM_WORLD);
			
			std::cout << "(" << mpi_rank << ") Number Of Results: " << messageRecvNumberOfResults[0] << " " << messageRecvNumberOfResults[1] << std::endl;
			
			for (int i = 0; i < mpi_size; i++)
			{
				for(int j = 0; j < messageRecvNumberOfResults[i]; j++)
				{
					Serializer recvResParticleSerializer;
					serializer_buffer = (std::stringstream *)recvResParticleSerializer.pGetBuffer();
			
					for(int k = 0; mpi_res_recv_buffer[(msgResRecvSize + 1)*MaxNumberOfResults*i+(msgResRecvSize + 1)*j+k] != '\0'; k++) 
					{
						(*serializer_buffer) << mpi_res_recv_buffer[(msgResRecvSize + 1)*MaxNumberOfResults*i+(msgResRecvSize + 1)*j+k];
					}
			
					recvResults[i][j] = new PointType();
					recvResParticleSerializer.load("nodes",recvResults[i][j]);

					std::cout << "(" << mpi_rank << ") Result point from process (" << i << "): (" << recvResults[i][j]->X() << " " << recvResults[i][j]->Y() << " " << recvResults[i][j]->Z() << ") with dist: " << recvResultsDistances[i][j] << std::endl;
				}
			}
			
			//////////////////////////////////////////////////////////////////////////////////////////////
			// 	RESULT REDUCTION 																		//
			//////////////////////////////////////////////////////////////////////////////////////////////
		   
		}
         
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////
        //// LOCAL FUCNTIONS	(own domain)									////
        ////////////////////////////////////////////////////////////////////////////

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
         
         PointerType SearchNearestPoint( PointType const& ThisPoint )
         {
           if( mPointBegin == mPointEnd )
             return this->NullPointer();

           PointerType Result            = *mPointBegin;
           CoordinateType ResultDistance = static_cast<CoordinateType>(DBL_MAX);
           SearchStructureType Box( CalculateCell(ThisPoint), mN );
           SearchNearestPointLocal( ThisPoint, Result, ResultDistance, Box );
           return Result;
         }
         
         //************************************************************************
         
         PointerType SearchNearestPoint( PointType const& ThisPoint, CoordinateType ResultDistance )
         {
           if( mPointBegin == mPointEnd )
             return this->NullPointer();

           PointerType Result = *mPointBegin;
           ResultDistance     = static_cast<CoordinateType>(DBL_MAX);
           SearchStructureType Box( CalculateCell(ThisPoint), mN );
           SearchNearestPointLocal( ThisPoint, Result, ResultDistance, Box);
           return Result;
         }

         //************************************************************************

         // New Thread Safe!!!
         PointerType SearchNearestPoint( PointType const& ThisPoint, CoordinateType& rResultDistance, SearchStructureType& Box )
         {
           PointerType Result = *mPointBegin; //static_cast<PointerType>(NULL);
           rResultDistance    = static_cast<CoordinateType>(DBL_MAX);
           Box.Set( CalculateCell(ThisPoint), mN );
           SearchNearestPointLocal( ThisPoint, Result, rResultDistance, Box);
           return Result;
         }
         
         //************************************************************************

         void SearchNearestPoint( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance )
         {
           SearchStructureType Box;
           Box.Set( CalculateCell(ThisPoint), mN );
           SearchNearestPointLocal(ThisPoint,rResult,rResultDistance,Box);
         }

         //************************************************************************

         void SearchNearestPoint( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance, SearchStructureType& Box )
         {
           // This case is when BinStatic is a LeafType in Other Spacial Structure
           // Then, it is possible a better Result before this search
           Box.Set( CalculateCell(ThisPoint), mN );
           SearchNearestPointLocal( ThisPoint, rResult, rResultDistance, Box );
         }

         //************************************************************************

         void SearchNearestPointLocal( PointType const& ThisPoint, PointerType& rResult, CoordinateType& rResultDistance, SearchStructureType& Box )
         {
           if( mPointBegin == mPointEnd )
             return;

           bool Found = false;

           // set mBox
           Box.Set( CalculateCell(ThisPoint), mN );

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
         
         SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, IteratorType Results, 
             DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults )
         {
           CoordinateType Radius2 = Radius * Radius;
           SizeType NumberOfResults = 0;
           SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
           SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
           return NumberOfResults;
         }

         //************************************************************************

         SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, IteratorType Results,
             DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
         {
           CoordinateType Radius2 = Radius * Radius;
           SizeType NumberOfResults = 0;
           Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
           SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
           return NumberOfResults;
         }

         //************************************************************************

         void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
             DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults )
         {
           SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
           SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
         }

         //************************************************************************

         void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
             DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
         {
           Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
           SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
         }

         //************************************************************************

         // **** THREAD SAFE

         // Dimension = 1
         void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
             DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
         {
           for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
             SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
         }

         // Dimension = 2
         void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
             DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
         {
           for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
             for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
               SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
         }

         // Dimension = 3
         void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
             DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
         {
           for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
             for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
               for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                 SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
		 }

         //************************************************************************
         
         SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType Radius, IteratorType Results, SizeType MaxNumberOfResults )
         {
           CoordinateType Radius2 = Radius * Radius;
           SizeType NumberOfResults = 0;
           SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
           SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
           return NumberOfResults;
         }

         //************************************************************************

         SizeType SearchInRadius( PointType const& ThisPoint, CoordinateType Radius, IteratorType Results,
             SizeType MaxNumberOfResults, SearchStructureType& Box )
         {
           CoordinateType Radius2 = Radius * Radius;
           SizeType NumberOfResults = 0;
           Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
           SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
           return NumberOfResults;
         }

         //************************************************************************

         void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
             SizeType& NumberOfResults, SizeType const& MaxNumberOfResults )
         {
           SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
           SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
         }

         //************************************************************************

         void SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
             SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
         {
           Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
           SearchInRadiusLocal( ThisPoint, Radius, Radius2, Results, NumberOfResults, MaxNumberOfResults, Box );
         }
         
         //************************************************************************
         
         // **** THREAD SAFE

         // Dimension = 1
         void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
             SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
         {
           for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I++ )
             SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
         }

         // Dimension = 2
         void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
             SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
         {
           for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
             for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I++ )
               SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
         }

         // Dimension = 3
         void SearchInRadiusLocal( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
             SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
         {
           for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
             for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
               for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I++ )
                 SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Results,NumberOfResults,MaxNumberOfResults);
         }
         
         //************************************************************************
         //************************************************************************

         // Dimension = 1
         void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance, 
             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, bool& Found )
         {
           Found = false;
           for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
             SearchNearestInRange()( mPoints[I].begin(), mPoints[I].end(), ThisPoint, ResultPoint, ResultDistance, Found );
         }

         // Dimension = 2
         void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance, 
             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, bool& Found )
         {
           Found = false;
           for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
             for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
               SearchNearestInRange()( mPoints[I].begin(), mPoints[I].end(), ThisPoint, ResultPoint, ResultDistance, Found );
         }

         // Dimension = 3
         void SearchNearestInBox( PointType const& ThisPoint, PointerType& ResultPoint, CoordinateType& ResultDistance, 
             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, bool& Found )
         {
           Found = false;
           for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
             for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
               for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                 SearchNearestInRange()( mPoints[I].begin(), mPoints[I].end(), ThisPoint, ResultPoint, ResultDistance, Found );
         }

         //************************************************************************
         //************************************************************************

         SizeType SearchInBox( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType Results,
             SizeType MaxNumberOfResults )
         {
           SizeType NumberOfResults = 0;
           SearchStructureType Box( CalculateCell(SearchMinPoint), CalculateCell(SearchMaxPoint), mN );
           SearchInBoxLocal( SearchMinPoint, SearchMaxPoint, Results, NumberOfResults, MaxNumberOfResults, Box );
           return NumberOfResults;
         }

         //************************************************************************

         void SearchInBox(PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& Results, SizeType& NumberOfResults,
             SizeType const& MaxNumberOfResults )
         {
           NumberOfResults = 0;
           SearchStructureType Box( CalculateCell(SearchMinPoint), CalculateCell(SearchMaxPoint), mN );
           SearchInBoxLocal( SearchMinPoint, SearchMaxPoint, Results, NumberOfResults, MaxNumberOfResults, Box );
         }

         //************************************************************************
         
         // Dimension = 1
         void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
             SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
         {
           for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
             SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,mPoints[I].begin(),mPoints[I].end(),ResultsPoint,NumberOfResults,MaxNumberOfResults);
         }

         // Dimension = 2
         void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
             SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
         {
           for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
             for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
               SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,mPoints[I].begin(),mPoints[I].end(),ResultsPoint,NumberOfResults,MaxNumberOfResults);
         }

         // Dimension = 3
         void SearchInBoxLocal( PointType const& SearchMinPoint, PointType const& SearchMaxPoint, IteratorType& ResultsPoint,
             SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
             SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
         {
           for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
             for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
               for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                 SearchBoxInRange()(SearchMinPoint,SearchMaxPoint,mPoints[I].begin(),mPoints[I].end(),ResultsPoint,NumberOfResults,MaxNumberOfResults);
         }
         
         //************************************************************************

         /// Turn back information as a string.
         virtual std::string Info() const
         {
            return "BinsDynamicMpi";
         }
      
         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const
         {
            rOStream << "BinsDynamicMpi";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const
         {
			rOStream << Perfix << "Bin[" << SearchUtils::PointerDistance(mPointBegin, mPointEnd) << "] : " << std::endl;
			for(typename CellsContainerType::const_iterator i_cell = mPoints.begin() ; i_cell != mPoints.end() ; i_cell++)
			{
			    rOStream << Perfix << "[ " ;
			    for(typename PointVector::const_iterator i_point = i_cell->begin() ; i_point != i_cell->end() ; i_point++)
				rOStream << **i_point << "    ";
			    rOStream << " ]" << std::endl;
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
         BinsDynamicMpi& operator=(BinsDynamicMpi const& rOther);

         /// Copy constructor.
         BinsDynamicMpi(BinsDynamicMpi const& rOther);

    private:

         IteratorType     mPointBegin;
         IteratorType     mPointEnd;
         
         Tvector<CoordinateType,TDimension>  mMinPoint;
         Tvector<CoordinateType,TDimension>  mMaxPoint;
         Tvector<CoordinateType,TDimension>  mCellSize;
         Tvector<CoordinateType,TDimension>  mInvCellSize;
         Tvector<SizeType,TDimension>        mN;
         SizeType                            mNumPoints;

         // Bins Access Vector ( vector<Iterator> )
         CellsContainerType mPoints;

         // Work Variables ( For non-copy of Search Variables )
         //BinBox SearchBox;
	 
	 //MPI_interface variables
	 int mpi_rank;
	 int mpi_size;

	public:
	  static TreeNodeType* Construct(IteratorType PointsBegin, IteratorType PointsEnd, PointType MaxPoint, PointType MinPoint, SizeType BucketSize)
	  {
		 
		SizeType number_of_points = SearchUtils::PointerDistance(PointsBegin,PointsEnd);
		if (number_of_points == 0)
		  return NULL;
		else 
		{
		  return new BinsDynamicMpi( PointsBegin, PointsEnd, MinPoint, MaxPoint, BucketSize );
		}

	  }

};

   template<
      std::size_t TDimension,
      class TPointType,
      class TContainerType,
      class TPointerType,
      class TIteratorType,
      class TDistanceIteratorType,
      class TDistanceFunction >
std::ostream & operator<<( std::ostream& rOStream,
      BinsDynamicMpi<TDimension,TPointType,TContainerType,TPointerType,TIteratorType,TDistanceIteratorType,TDistanceFunction>& rThis)
{
   rThis.PrintInfo(rOStream);
   rOStream << std::endl;
   rThis.PrintSize(rOStream);
   rThis.PrintData(rOStream);
   return rOStream;
};



};

#endif // KRATOS_BINS_DYNAMIC_MPI_CONTAINER_H_INCLUD
