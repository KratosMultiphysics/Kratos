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

namespace Kratos {

  /// This class its an implementation of BinsDynamic using MPI
  /**
   * Use the seam way you use the generic BinsDynamic
   */
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
	typedef Kratos::SearchUtils::SquaredDistanceFunction<TDimension,PointType> SquaredDistanceFunction;


	public:

		/// Default constructor.
		/**
		 * Empy constructor, you shouldn't use this unless you know what you are doing.
		 */
		BinsDynamicMpi() : mPointBegin(this->NullIterator()), mPointEnd(this->NullIterator()), mNumPoints(0)
		{};
        
		/// ModelPart Constructor.
		/**
		 * Creates and initializes BinsDynamic using the LocalMesh of the ModelPart provided as argument.
		 * @param StaticMesh The geometry used to generate the Bins
		 * @param ParticMesh Not used atm
		 * @param BoxSize Size of the box
		 * @param BucketSize default = 1
		 */
        BinsDynamicMpi( ModelPart * StaticMesh, ModelPart * ParticMesh, CoordinateType BoxSize, SizeType BucketSize = 1 )
        {
			IteratorType mPointIterator;
			
			mPointBegin = new PointType* [StaticMesh->GetCommunicator().LocalMesh().NumberOfNodes()];
			mPointEnd = mPointBegin + StaticMesh->GetCommunicator().LocalMesh().NumberOfNodes();
			
			mPointIterator = mPointBegin;
		  
			std::cout << "Parsing local mesh elements: " << StaticMesh->GetCommunicator().LocalMesh().NumberOfNodes() << std::endl;
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
		  
			MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
			MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
			
			std::cout << "Initialize done" << std::endl;
        }
        
        //************************************************************************

		/// Destructor.
        virtual ~BinsDynamicMpi(){ }
        
        //************************************************************************
        
        /// Pointer to first element.
        IteratorType Begin() { return mPointBegin; }

        //************************************************************************
        
        /// Pointer to last element.
        IteratorType End() { return mPointBegin; }

        //************************************************************************
        
        /// Size of specific dimension.
        CoordinateType CellSize( SizeType const& iDim ) { return mCellSize[iDim]; }
        
        //************************************************************************
        
        /// Number of cells of specific dimension.
        SizeType NumCell( SizeType const& iDim ) { return mN[iDim]; }
        
        //************************************************************************
        
        /// Calcutes bounding box of particles in bins
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
        
        /// Calcutes cell Size
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
        
        /// Calcutes cell Size gived the container box
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
	
		//************************************************************************
		
		///////////////////////////////////////////////////////////////////////////
		// MPI Single Input Search 
		///////////////////////////////////////////////////////////////////////////
		
		void MPISingleSearchInRadiusTest() 
		{
			//Parameters
			int NumberOfPoints 	= 5;
			int MaxNumberOfResults = 5;
			double Radius		= 0.1f;
		  
			//MultiSearch Test
			PointerType  PointInput = new PointType[NumberOfPoints];
			IteratorType Results    = new PointerType[MaxNumberOfResults];
			
			DistanceIteratorType Distances = new double[MaxNumberOfResults];
			
			for(int i = 0; i < NumberOfPoints; i++) 
			{
				PointType temp;

				temp[0] = (i+1)/NumberOfPoints;
				temp[1] = (i+1)/NumberOfPoints;
				temp[2] = 0;

				PointInput[i] = PointType(temp);
			}
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			MPI_SearchInRadius(PointInput[NumberOfPoints/2], Radius, Results, Distances, MaxNumberOfResults);
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			//Check Results
// 			for(int i = 0; i < res; i++) 
// 				std::cout << "(" << mpi_rank << ")" << (*Results[i]) << "\tDIST\t" << Distances[i] << std::endl;
		}
	
		SizeType MPI_SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, IteratorType Results, 
			DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults)
		{
			CoordinateType Radius2 = Radius * Radius;
			SizeType NumberOfResults = 0;
			SearchInRadiusMpiWrapper( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults );
		   
			return NumberOfResults;
		}

		//************************************************************************

// 		SizeType MPI_SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, IteratorType Results,
// 			DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
// 		{
// 			CoordinateType Radius2 = Radius * Radius;
// 			SizeType NumberOfResults = 0;
// 			Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
// 			SearchInRadiusMpiWrapper( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
// 			return NumberOfResults;
// 		}

		//************************************************************************

// 		SizeType MPI_SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
// 			DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults )
// 		{
// 			SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
// 			SearchInRadiusMpiWrapper( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
// 		}

		//************************************************************************

// 		SizeType MPI_SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
// 			DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
// 		{
// 			Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
// 			SearchInRadiusMpiWrapper( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
// 		}

		void SearchInRadiusMpiWrapper( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
             DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults )
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

				for(size_t j = 0; mpi_recv_buffer[(msgRecvSize+1)*i+j] != '\0'; j++) 
				{
					(*serializer_buffer) << mpi_recv_buffer[(msgRecvSize+1)*i+j];
				}

				remoteNumberOfResults[i] = 0;

				recvParticleSerializer.load("nodes",remoteThisPoint[i]);

				SearchStructureType Box( CalculateCell(remoteThisPoint[i],-Radius), CalculateCell(remoteThisPoint[i],Radius), mN );
				SearchInRadiusLocal(remoteThisPoint[i],Radius,Radius2,remoteResultsPointer,remoteResultsDistancesPointer,remoteNumberOfResults[i],MaxNumberOfResults,Box);

				for(size_t j = 0; j < remoteNumberOfResults[i]; j++)
				{
					std::cout << "(" << mpi_rank << ")\t" << *(remoteResults[i][j]) << " " << remoteResultsDistances[i][j] << std::endl;	
				}
			}

			std::stringstream * res_serializer_buffer[mpi_size][MaxNumberOfResults];
			std::string message[mpi_size][MaxNumberOfResults];
			int bufferSize[mpi_size][MaxNumberOfResults];
	
			for(int i = 0; i < mpi_size; i++) 
			{
				for(size_t j = 0; j < remoteNumberOfResults[i]; j++) 
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

			char mpi_res_send_buffer[((msgResRecvSize + 1) * mpi_size * MaxNumberOfResults)];
			char mpi_res_recv_buffer[((msgResRecvSize + 1) * mpi_size * MaxNumberOfResults)];

			for(int i = 0; i < mpi_size; i++) 
			{
				messageSendNumberOfResults[i] = remoteNumberOfResults[i];
				for(size_t j = 0; j < remoteNumberOfResults[i]; j++) 
				{
					strcpy(&mpi_res_send_buffer[(msgResRecvSize + 1)*MaxNumberOfResults*i+(msgResRecvSize + 1)*j], message[i][j].c_str());
					mpi_res_send_buffer[(msgResRecvSize + 1)*MaxNumberOfResults*i+(msgResRecvSize + 1)*j+bufferSize[i][j]] = '\0';
				}
			}

			MPI_Alltoall(mpi_res_send_buffer,((msgResRecvSize+1) * MaxNumberOfResults),MPI_CHAR,mpi_res_recv_buffer,((msgResRecvSize+1) * MaxNumberOfResults),MPI_CHAR,MPI_COMM_WORLD);
			MPI_Alltoall(messageSendNumberOfResults,1,MPI_INT,messageRecvNumberOfResults,1,MPI_INT,MPI_COMM_WORLD);
			MPI_Alltoall(remoteResultsDistances[0],MaxNumberOfResults,MPI_DOUBLE,recvResultsDistances[0],MaxNumberOfResults,MPI_DOUBLE,MPI_COMM_WORLD);

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

// 					std::cout << "(" << mpi_rank << ") Result point from process (" << i << "): (" << recvResults[i][j]->X() << " " << recvResults[i][j]->Y() << " " << recvResults[i][j]->Z() << ") with dist: " << recvResultsDistances[i][j] << std::endl;
				}
			}
		}
		
		///////////////////////////////////////////////////////////////////////////
		// MPI Single Input END
		///////////////////////////////////////////////////////////////////////////
		
		///////////////////////////////////////////////////////////////////////////
		// MPI Multiple Input Search 
		///////////////////////////////////////////////////////////////////////////
		
		void MPIMultiSearchInRadiusTest(const SizeType& NumberOfPoints, const SizeType& MaxNumberOfResults, const double& Radius, const SizeType& times)
		{
			//MultiSearch Test
			PointerType  PointInput = new PointType[NumberOfPoints];
			for(int i = 0; i < NumberOfPoints; i++) 
			{
				PointType temp;
				
				temp[0] = (double)(rand()%1000)/(double)1000;//PointBegin[i]->X();
				temp[1] = (double)(rand()%1000)/(double)1000;//PointBegin[i]->Y();
				temp[2] = (double)(rand()%1000)/(double)1000;//mPointBegin[i]->Z();

				PointInput[i] = PointType(temp);
			}
			
			for(int i = 0; i < times; i++)
			{
				vector<SizeType> NumberOfResults(NumberOfPoints);
				vector<vector<PointerType> > Results(NumberOfPoints, vector<PointerType>(MaxNumberOfResults));
				vector<vector<double> > ResultsDistances(NumberOfPoints, vector<double>(MaxNumberOfResults,0));
			
				MPI_SearchInRadius(PointInput, NumberOfPoints, Radius, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults);
				MPI_Barrier(MPI_COMM_WORLD);
			}
			
			//Check Results
// 			for(int i = 0; i < res; i++) 
// 				std::cout << "(" << mpi_rank << ")" << (*Results[i]) << "\tDIST\t" << Distances[i] << std::endl;
		}
	
		void MPI_SearchInRadius( PointerType const& ThisPoints, SizeType const& NumberOfPoints, CoordinateType const& Radius, vector<vector<PointerType> > Results, 
			vector<vector<double> > ResultsDistances, vector<SizeType> NumberOfResults, SizeType const& MaxNumberOfResults)
		{
			CoordinateType Radius2 = Radius * Radius;
			
			SearchInRadiusMpiWrapper( ThisPoints, NumberOfPoints, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults );
		}

		/// Act as wrapper between external function and its implementation
		/**
		 * This function provides all mpi functionality requiered to execute the parallel multi input searchInRaidus.
		 * the method implemented by this function is ALL-vs-ALL. It means all particles not found in the local
		 * processes are send to all other processes
		 * @param ThisPoints List of points to be search
		 * @param NumberOfPoints Number of points to be search 
		 * @param Radius Radius of search
		 * @param Radius2 Radius of search ^2
		 * @param Results List of results
		 * @param ResultsDistances Distance of the results
		 * @param NumberOfResults Number of results
		 * @param MaxNumberOfResults Maximum number of results returned for each point
		 */
		void SearchInRadiusMpiWrapper( PointerType const& ThisPoints, SizeType const& NumberOfPoints, CoordinateType const& Radius, CoordinateType const& Radius2, vector<vector<PointerType> > Results,
             vector<vector<double> > ResultsDistances, vector<SizeType> NumberOfResults, SizeType const& MaxNumberOfResults )
		{  
			PointType remoteThisPoints[mpi_size][NumberOfPoints];
			
			vector<vector<vector<PointerType> > > remoteResults(mpi_size, vector<vector<PointerType> >(NumberOfPoints, vector<PointerType>(MaxNumberOfResults)));
			vector<vector<vector<double> > > remoteResultsDistances(mpi_size, vector<vector<double> >(NumberOfPoints, vector<double>(MaxNumberOfResults)));
			
			SizeType remoteNumberOfResults[mpi_size][NumberOfPoints];

			int NumberOfSendPoints = 0;
			int MaxNumberOfSendPoints = 0;
			int RemoteNumberOfSendPoints[mpi_size];
			
			int msgSendSize = 0;
			int msgRecvSize = 0;
			
			//Local search
			for(int i = 0; i < NumberOfPoints; i++)
			{
			  	IteratorType ResultsPointer      = &Results[i][0];
				double * ResultsDistancesPointer = &ResultsDistances[i][0];
			  
				NumberOfResults[i] = 0;
				
				SearchStructureType Box( CalculateCell(ThisPoints[i],-Radius), CalculateCell(ThisPoints[i],Radius), mN );
				SearchInRadiusLocal(ThisPoints[i],Radius,Radius2,ResultsPointer,ResultsDistancesPointer,NumberOfResults[i],MaxNumberOfResults,Box);

				if(NumberOfResults[i] < MaxNumberOfResults) NumberOfSendPoints++;
			}
			
			std::cout << "(" << mpi_rank << ") N: " << NumberOfSendPoints << std::endl;
			
			//Only search points not found previously in local mesh
			std::stringstream * serializer_buffer[NumberOfSendPoints];
			std::string message[NumberOfSendPoints];
			
			for(int i = 0, j = 0; i < NumberOfPoints; i++) 
			{
				if(NumberOfResults[i] < MaxNumberOfResults)
				{
					Serializer particleSerializer;
					const PointType& ThisPoint = ThisPoints[i];
					
					particleSerializer.save("nodes",&ThisPoint);
					
					serializer_buffer[j] = (std::stringstream *)particleSerializer.pGetBuffer();
					message[j] = std::string(serializer_buffer[j]->str().c_str());
					msgSendSize = msgSendSize > serializer_buffer[j]->str().size() ? msgSendSize : serializer_buffer[j]->str().size();
					j++;
				}
			}

			//Message Size commuincation
			MPI_Allreduce(&msgSendSize,&msgRecvSize,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
			//Max number of particles to be search by any process
			MPI_Allreduce(&NumberOfSendPoints,&MaxNumberOfSendPoints,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
			//Number of particles requested to be search by one process
			MPI_Allgather(&NumberOfSendPoints,1,MPI_INT,RemoteNumberOfSendPoints,1,MPI_INT,MPI_COMM_WORLD);
			
			char mpi_send_buffer[(msgRecvSize+1) * MaxNumberOfSendPoints];
			char mpi_recv_buffer[(msgRecvSize+1) * MaxNumberOfSendPoints * mpi_size];
				
			int messageNumberOfResults[mpi_size][MaxNumberOfSendPoints];
			
			std::cout << "(" << mpi_rank << ") Max number Of Send Points: " << MaxNumberOfSendPoints << " : " << NumberOfSendPoints << std::endl;
			
			for(size_t i = 0; i < NumberOfSendPoints; i++) 
			{
				strcpy(&mpi_send_buffer[(msgRecvSize + 1)*i], message[i].c_str());
				
				//Some particles have serialization buffer with size < msgRecvSize, we need to fill the rest of string to parse it later
				for(int j = message[i].size(); j < (msgRecvSize + 1); j++)
					mpi_send_buffer[(msgRecvSize + 1)*i+j] = '\0';
			}
			
			std::cout << "(" << mpi_rank << ") Buffer formated ok" << std::endl;
			//Particle data
			MPI_Allgather(mpi_send_buffer,(msgRecvSize+1)*MaxNumberOfSendPoints,MPI_CHAR,mpi_recv_buffer,(msgRecvSize+1)*MaxNumberOfSendPoints,MPI_CHAR,MPI_COMM_WORLD);
			
			std::cout << "(" << mpi_rank << ") "<< "(" << getpid() << ") Begin Search" << std::endl;
			
			for(int i = 0; i < mpi_size; i++) 
			{
				size_t k = 0;
				if(i != mpi_rank)
				{
					//MUST resize vectors or declare them here in order to avoid NULL problems later in result serialization
					remoteResults[i].resize(RemoteNumberOfSendPoints[i]);
					remoteResultsDistances[i].resize(MaxNumberOfSendPoints);

					for(size_t j = 0; j < RemoteNumberOfSendPoints[i]; j++) 
					{
						remoteResults[i][j].resize(MaxNumberOfResults);
						remoteResultsDistances[i][j].resize(MaxNumberOfResults);
						
						Serializer particleSerializer;
						serializer_buffer[0] = (std::stringstream *)particleSerializer.pGetBuffer();

						for(; mpi_recv_buffer[(msgRecvSize+1)*MaxNumberOfSendPoints*i+k] != '\0'; k++)
							(*serializer_buffer[0]) << mpi_recv_buffer[(msgRecvSize+1)*MaxNumberOfSendPoints*i+k];
						while (mpi_recv_buffer[(msgRecvSize+1)*MaxNumberOfSendPoints*i+k] == '\0') k++;

						PointerType remoteThisPointsPointer = &remoteThisPoints[i][j];
						particleSerializer.load("nodes",remoteThisPointsPointer);
						
						IteratorType remoteResultsPointer      = &remoteResults[i][j][0];
						double * remoteResultsDistancesPointer = &remoteResultsDistances[i][j][0];

						remoteNumberOfResults[i][j] = 0;

						SearchStructureType Box( CalculateCell(remoteThisPoints[i][j],-Radius), CalculateCell(remoteThisPoints[i][j],Radius), mN );
						SearchInRadiusLocal(remoteThisPoints[i][j],Radius,Radius2,remoteResultsPointer,remoteResultsDistancesPointer,remoteNumberOfResults[i][j],MaxNumberOfResults,Box);
						
// 						remoteResults[i][j].resize(remoteNumberOfResults[i][j]);
						for(int l = remoteNumberOfResults[i][j]; l < MaxNumberOfResults; l++) 
							remoteResults[i][j][l] = new PointType();
						
/*						std::cout << "(" << mpi_rank << ") Found points for: (" << remoteThisPoints[i][j].X() << " " << remoteThisPoints[i][j].Y() << " " << remoteThisPoints[i][j].Z() << ") from process(" << i << "): FOUND: " << remoteNumberOfResults[i][j] << std::endl; 
						for(size_t k = 0; k < remoteNumberOfResults[i][j]; k++) 
							std::cout << "(" << mpi_rank << ")\t" << *(remoteResults[i][j][k]) << " " << remoteResultsDistances[i][j][k] << std::endl;    */   
					}
				}
			}
			
			std::cout << "(" << mpi_rank << ") Found elements from foregin proceeses" << std::endl;
			std::stringstream * res_serializer_buffer[mpi_size];
			std::string res_message[mpi_size];
			int res_bufferSize[mpi_size];
			
			msgSendSize = 0;
			msgRecvSize = 0;
	
			for(int i = 0; i < mpi_size; i++)
			{
				if(RemoteNumberOfSendPoints[i] && mpi_rank != i)
				{
					Serializer resSerializer;
					resSerializer.save("nodes",remoteResults[i]);
					  
// 					std::cout << "(" << mpi_rank << ") Serializing result buffer of size: " << remoteResults[i].size() << " for process: " << i << std::endl;
					res_serializer_buffer[i] = (std::stringstream *)resSerializer.pGetBuffer();
					res_message[i] = std::string(res_serializer_buffer[i]->str().c_str());
					res_bufferSize[i] = res_serializer_buffer[i]->str().size();
					
					msgSendSize = msgSendSize > res_bufferSize[i] ? msgSendSize : res_bufferSize[i];
				}
			}
			
			//Result buffer size
			MPI_Allreduce(&msgSendSize,&msgRecvSize,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

			char mpi_res_send_buffer[(msgRecvSize+1) * mpi_size];
			char mpi_res_recv_buffer[(msgRecvSize+1) * mpi_size];

			for(int i = 0; i < mpi_size; i++)
			{
				if(RemoteNumberOfSendPoints[i] && mpi_rank != i)
				{
					strcpy(&mpi_res_send_buffer[(msgRecvSize+1)*i],res_message[i].c_str());
					mpi_res_send_buffer[(msgRecvSize+1)*i+res_bufferSize[i]] = '\0';
				}
				
				for(size_t j = 0; j < RemoteNumberOfSendPoints[i]; j++) messageNumberOfResults[i][j] = remoteNumberOfResults[i][j];
				for(size_t j = RemoteNumberOfSendPoints[i]; j < MaxNumberOfSendPoints; j++) messageNumberOfResults[i][j] = 0;
			}

			//Number of results of each point of each process
			MPI_Alltoall(messageNumberOfResults,MaxNumberOfSendPoints,MPI_INT,messageNumberOfResults,MaxNumberOfSendPoints,MPI_INT,MPI_COMM_WORLD);
			//Results data
			MPI_Alltoall(mpi_res_send_buffer,(msgRecvSize+1),MPI_CHAR,mpi_res_recv_buffer,(msgRecvSize+1),MPI_CHAR,MPI_COMM_WORLD);

			for (int i = 0; i < mpi_size; i++)
			{
				if (i != mpi_rank)
				{
					Serializer particleSerializer;
					serializer_buffer[0] = (std::stringstream *)particleSerializer.pGetBuffer();
					
					for(int l = (msgRecvSize+1)*i; mpi_res_recv_buffer[l] != '\0'; l++) 
						(*serializer_buffer[0]) << mpi_res_recv_buffer[l];

					vector<vector<PointerType> > tResults(NumberOfSendPoints, vector<PointerType>(MaxNumberOfResults));
					for(int j = 0; j < NumberOfSendPoints; j++)
						for(int k = 0; k < MaxNumberOfResults; k++)
							tResults[j][k] = new PointType();
					
					particleSerializer.load("nodes",tResults);
					
					for(size_t j = 0; j < NumberOfPoints; j++)
					{
						if(NumberOfResults[i] < MaxNumberOfResults)
						{
							for(int k = 0; k < messageNumberOfResults[i][j]; k++)
							{
								Results[j][NumberOfResults[j]] = tResults[j][k];
								
								PointType& a = (ThisPoints[j]);
								PointType& b = (*Results[j][NumberOfResults[j]]);
								
								ResultsDistances[j][NumberOfResults[j]] = DistanceFunction()(a,b);
								NumberOfResults[j]++;
							}
						}
// 						std::cout << "(" << mpi_rank << ") Found " <<  NumberOfResults[j] << " results for point: " << ThisPoints[j] << std::endl;
// 						for(int k = 0; k < NumberOfResults[j]; k++)
// 						{
// 							std::cout << "(" << mpi_rank << ")\t(" << Results[j][k]->X() << " " << Results[j][k]->Y() << " " << Results[j][k]->Z() << ") with dist: " << ResultsDistances[j][k] << std::endl;
// 						}
					}
				}
			}
		}

		///////////////////////////////////////////////////////////////////////////
		// MPI END Multi Input End
		///////////////////////////////////////////////////////////////////////////

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
			while(!Found)
			{
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
         
		///////////////////////////////////////////////////////////////////////////
		// Multiple Input Search 
		///////////////////////////////////////////////////////////////////////////
		
		void MultiSearchInRadiusTest() 
		{
			//Parameters
			int NumberOfPoints 	= 5;
			int NumberOfResults = 5;
			double Radius		= 0.1f;
		  
			//MultiSearch Test
			PointerType  PointInput = new PointType[NumberOfPoints];
			IteratorType Results    = new PointerType[NumberOfPoints * NumberOfResults];
			
			DistanceIteratorType Distances = new double[NumberOfPoints  * NumberOfResults];
			
			for(int i = 0; i < NumberOfPoints; i++) 
			{
				PointType temp;

				temp[0] = (i+1)/NumberOfPoints;
				temp[1] = (i+1)/NumberOfPoints;
				temp[2] = 0;

				PointInput[i] = PointType(temp);
			}
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			int res = SearchInRadius(PointInput, NumberOfPoints, Radius, Results, Distances, NumberOfResults);
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			//Check Results
			for(int i = 0; i < res; i++) 
				std::cout << "(" << mpi_rank << ")" << (*Results[i]) << "\tDIST\t" << Distances[i] << std::endl;
		}
 
		SizeType SearchInRadius( PointerType const& ThisPoints, SizeType const& NumberOfPoints, CoordinateType const& Radius, IteratorType Results, 
			DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults )
		{
			CoordinateType Radius2 = Radius * Radius;
			SizeType NumberOfResults = 0;
			SearchStructureType Box[NumberOfPoints];
			
			for(size_t i = 0; i < NumberOfPoints; i++) 
				Box[i] = SearchStructureType( CalculateCell(ThisPoints[i],-Radius), CalculateCell(ThisPoints[i],Radius), mN );
			
			SearchInRadiusLocal( ThisPoints, NumberOfPoints, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
			
			return NumberOfResults;
		}

		//************************************************************************

		SizeType SearchInRadius( PointerType const& ThisPoints, SizeType const& NumberOfPoints, CoordinateType const& Radius, IteratorType Results,
			DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults, SearchStructureType Box[] )
		{
			CoordinateType Radius2 = Radius * Radius;
			SizeType NumberOfResults = 0;
			
			for(size_t i = 0; i < NumberOfPoints; i++) 
				Box[i].Set( CalculateCell(ThisPoints[i],-Radius), CalculateCell(ThisPoints[i],Radius), mN );
			
			SearchInRadiusLocal( ThisPoints, NumberOfPoints, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
			return NumberOfResults;
		}

		//************************************************************************

		void SearchInRadius( PointerType const& ThisPoints, SizeType const& NumberOfPoints, CoordinateType const& Radius, CoordinateType const& Radius2, 
			IteratorType& Results, DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults )
		{
			SearchStructureType Box[NumberOfPoints];
			
			for(size_t i = 0; i < NumberOfPoints; i++) 
				Box[i] = SearchStructureType( CalculateCell(ThisPoints[i],-Radius), CalculateCell(ThisPoints[i],Radius), mN );
			
			SearchInRadiusLocal( ThisPoints, NumberOfPoints, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
		}

		//************************************************************************

		void SearchInRadius( PointerType const& ThisPoints, SizeType const& NumberOfPoints, CoordinateType const& Radius, CoordinateType const& Radius2, 
			IteratorType& Results, DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType Box[] )
		{
			for(size_t i = 0; i < NumberOfPoints; i++) 
				Box[i].Set( CalculateCell(ThisPoints[i],-Radius), CalculateCell(ThisPoints[i],Radius), mN );
			SearchInRadiusLocal( ThisPoints, NumberOfPoints, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
		}

		//************************************************************************

		// Dimension = 1
		void SearchInRadiusLocal( PointerType const& ThisPoint, SizeType const& NumberOfPoints, CoordinateType const& Radius, CoordinateType const& Radius2, 
			IteratorType& Results, DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
			SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
		{
			for(size_t i = 0; i < NumberOfPoints; i++)
			{
			  	SizeType thisNumberOfResults = 0;
				IteratorType ResulPointer = &Results[NumberOfResults];
				DistanceIteratorType ResultsDistancesPointer = &ResultsDistances[NumberOfResults];
				for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
					SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
				NumberOfResults += thisNumberOfResults;
			}
		}

		// Dimension = 2
		void SearchInRadiusLocal( PointerType const& ThisPoint, SizeType const& NumberOfPoints, CoordinateType const& Radius, CoordinateType const& Radius2, 
			IteratorType& Results, DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
			SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
		{
		  	for(size_t i = 0; i < NumberOfPoints; i++)
			{
				SizeType thisNumberOfResults = 0;
				IteratorType ResulPointer = &Results[NumberOfResults];
				DistanceIteratorType ResultsDistancesPointer = &ResultsDistances[NumberOfResults];
				for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
					for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
						SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint,Radius2,Results,ResultsDistances,NumberOfResults,MaxNumberOfResults);
				NumberOfResults += thisNumberOfResults;
			}
		}

		// Dimension = 3
		void SearchInRadiusLocal( PointerType const& ThisPoint, SizeType const& NumberOfPoints, CoordinateType const& Radius, CoordinateType const& Radius2, 
			IteratorType& Results, DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults,
			SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3> * const& Box )
		{
			for(size_t i = 0; i < NumberOfPoints; i++)
			{
				SizeType thisNumberOfResults = 0;
				IteratorType ResulPointer = &Results[NumberOfResults];
				DistanceIteratorType ResultsDistancesPointer = &ResultsDistances[NumberOfResults];
				for(IndexType III = Box[i].Axis[2].Begin() ; III <= Box[i].Axis[2].End() ; III += Box[i].Axis[2].Block )
					for(IndexType II = III + Box[i].Axis[1].Begin() ; II <= III + Box[i].Axis[1].End() ; II += Box[i].Axis[1].Block )
						for(IndexType I = II + Box[i].Axis[0].Begin() ; I <= II + Box[i].Axis[0].End() ; I += Box[i].Axis[0].Block )
							SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint[i],Radius2,ResulPointer,ResultsDistancesPointer,thisNumberOfResults,MaxNumberOfResults);
				NumberOfResults += thisNumberOfResults;
			}
		}         
		 	 
		///////////////////////////////////////////////////////////////////////////
		// Thread Safe ?
		///////////////////////////////////////////////////////////////////////////
		        
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

         IteratorType mPointBegin;
         IteratorType mPointEnd;
         
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
