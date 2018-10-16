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

#include "includes/serializer.h"
#include "utilities/timer.h"

namespace Kratos
{
  /// This class its an implementation of BinsDynamic using MPI
  /**
   * Use the seam way you use the generic BinsDynamic
   */
template<class TConfigure>
class BinsDynamicMpi
{
public:

    enum { Dimension = TConfigure::Dimension };

    typedef TConfigure                                      Configure;
    typedef typename TConfigure::PointType                  PointType;
    typedef typename TConfigure::PointVector                ContainerType;
    typedef typename TConfigure::PointIterator              IteratorType;
    typedef typename TConfigure::DistanceIterator           DistanceIteratorType;
    typedef typename TConfigure::PtrPointType               PointerType;
    typedef typename TConfigure::DistanceFunction           DistanceFunction;

    typedef std::vector<PointerType>                        PointVector;
    typedef std::vector<PointVector>                        CellsContainerType;
    typedef typename PointVector::iterator                  PointIterator;

    typedef TreeNode<Dimension,PointType,PointerType,IteratorType,DistanceIteratorType> TreeNodeType;

    typedef typename TreeNodeType::CoordinateType  CoordinateType;  // double
    typedef typename TreeNodeType::SizeType        SizeType;        // std::size_t
    typedef typename TreeNodeType::IndexType       IndexType;       // std::size_t

    typedef TreeNodeType LeafType;

    typedef typename TreeNodeType::IteratorIteratorType IteratorIteratorType;
    typedef typename TreeNodeType::SearchStructureType SearchStructureType;

    typedef Tvector<IndexType,Dimension>   CellType;
    typedef Kratos::SearchUtils::SearchNearestInRange<PointType,PointerType,PointIterator,DistanceFunction,CoordinateType> SearchNearestInRange;
    typedef Kratos::SearchUtils::SearchRadiusInRange<PointType,PointIterator,DistanceIteratorType,DistanceFunction,SizeType,CoordinateType,IteratorType> SearchRadiusInRange;
    typedef Kratos::SearchUtils::SearchBoxInRange<PointType,PointIterator,SizeType,Dimension,IteratorType> SearchBoxInRange;
    typedef Kratos::SearchUtils::SquaredDistanceFunction<Dimension,PointType> SquaredDistanceFunction;

    /// Pointer definition of BinsDynamicMpi
    KRATOS_CLASS_POINTER_DEFINITION(BinsDynamicMpi);

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
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        IteratorType mPointIterator;

        this->StaticMesh = StaticMesh;

        if (mpi_size != 1)
        {
            mPointBegin = new PointType* [StaticMesh->GetCommunicator().LocalMesh().NumberOfNodes()];
            mPointEnd = mPointBegin + StaticMesh->GetCommunicator().LocalMesh().NumberOfNodes();

            std::cout << "Local Mesh has: " << StaticMesh->GetCommunicator().LocalMesh().NumberOfNodes() << " Nodes" << std::endl;
            std::cout << "Ghost Mesh has: " << StaticMesh->GetCommunicator().GhostMesh().NumberOfNodes() << " Nodes" << std::endl;
            std::cout << "Full  Mesh has: " << StaticMesh->NumberOfNodes() << " Nodes" << std::endl;

            mPointIterator = mPointBegin;
            for( ModelPart::NodesContainerType::iterator inode = StaticMesh->GetCommunicator().LocalMesh().NodesBegin(); inode != StaticMesh->GetCommunicator().LocalMesh().NodesEnd(); inode++, mPointIterator++)
            {
                PointType auxPoint;

                auxPoint[0] = inode->X();
                auxPoint[1] = inode->Y();
                auxPoint[2] = inode->Z();

                (*mPointIterator) = new PointType(auxPoint);

//                     std::cout << "(" << mpi_rank << ") "<< inode->Id() << " - " << inode->X() << " " << inode->Y() << " " << inode->Z() << std::endl;
            }
        }
        else
        {
            mPointBegin = new PointType* [StaticMesh->NumberOfNodes()];
            mPointEnd = mPointBegin + StaticMesh->NumberOfNodes();

            mPointIterator = mPointBegin;
            for( ModelPart::NodesContainerType::iterator inode = StaticMesh->NodesBegin(); inode != StaticMesh->NodesEnd(); inode++, mPointIterator++)
            {
                PointType auxPoint;

                auxPoint[0] = inode->X();
                auxPoint[1] = inode->Y();
                auxPoint[2] = inode->Z();

                (*mPointIterator) = new PointType(auxPoint);
            }
        }


        if(mPointBegin==mPointEnd)
          return;

        mNumPoints = std::distance(mPointBegin,mPointEnd);
        CalculateBoundingBox();
        CalculateCellSize(BoxSize);
        AllocateCellsContainer();
        GenerateBins();
        GenerateCommunicationGraph();
    }

    //************************************************************************

    /// Destructor.
    virtual ~BinsDynamicMpi()
    {
        char msg[12] = {'b','i','n','s','_','X','.','t','i','m','e','\0'};
        msg[5] = '0' + mpi_rank;
        Timer::SetOuputFile(msg);
        Timer::PrintTimingInformation();
    }

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
    void CalculateBoundingBox()
    {
        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            mMinPoint[i] = (**mPointBegin)[i];
            mMaxPoint[i] = (**mPointBegin)[i];
        }

        for(IteratorType Point = mPointBegin ; Point != mPointEnd ; Point++)
        {
            for(SizeType i = 0 ; i < Dimension ; i++)
            {
                if( (**Point)[i] < mMinPoint[i] ) mMinPoint[i] = (**Point)[i];
                if( (**Point)[i] > mMaxPoint[i] ) mMaxPoint[i] = (**Point)[i];
            }
        }
    }

    //************************************************************************

    /// Calcutes cell Size
    void CalculateCellSize()
    {

      CoordinateType delta[Dimension];
      CoordinateType alpha[Dimension];
      CoordinateType mult_delta = 1.00;
      SizeType index = 0;
      for(SizeType i = 0 ; i < Dimension ; i++) {
        delta[i] = mMaxPoint[i] - mMinPoint[i];
        if ( delta[i] > delta[index] )
          index = i;
        delta[i] = (delta[i] == 0.00) ? 1.00 : delta[i];
      }

      for(SizeType i = 0 ; i < Dimension ; i++){
        alpha[i] = delta[i] / delta[index];
        mult_delta *= alpha[i];
      }

      mN[index] = static_cast<SizeType>( pow(static_cast<CoordinateType>(SearchUtils::PointerDistance(mPointBegin,mPointEnd)/mult_delta), 1.00/Dimension)+1 );

      for(SizeType i = 0 ; i < Dimension ; i++){
        if(i!=index) {
          mN[i] = static_cast<SizeType>(alpha[i] * mN[index]);
          mN[i] = ( mN[i] == 0 ) ? 1 : mN[i];
        }
      }

//       for(SizeType i = 0 ; i < Dimension ; i++){
//         mN[i] *= mpi_size;
//       }

      for(SizeType i = 0 ; i < Dimension ; i++){
        mCellSize[i] = delta[i] / mN[i];
        mInvCellSize[i] = 1.00 / mCellSize[i];
      }

    }

    //************************************************************************

    /// Calcutes cell Size gived the container box
    void CalculateCellSize( CoordinateType BoxSize ) {
        for(SizeType i = 0 ; i < Dimension ; i++){
            mCellSize[i] = BoxSize;
            mInvCellSize[i] = 1.00 / mCellSize[i];
            mN[i] = static_cast<SizeType>( (mMaxPoint[i]-mMinPoint[i]) / mCellSize[i]) + 1;
        }
      }

    //************************************************************************

    void AllocateCellsContainer() {
        SizeType Size = 1;
        for(SizeType i = 0 ; i < Dimension ; i++)
            Size *= mN[i];
        // Resize Global Container
        mPoints.resize(Size);
    }

    //************************************************************************

    void GenerateBins(){

        for(IteratorType i_point = mPointBegin ; i_point != mPointEnd ; i_point++)
          mPoints[CalculateIndex(**i_point)].push_back(*i_point);

    }

    void GenerateCommunicationGraph()
    {
        double * MpiMinPoints = new double[mpi_size * Dimension];
        double * MpiMaxPoints = new double[mpi_size * Dimension];

        double * MyMinPoint = new double[Dimension];
        double * MyMaxPoint = new double[Dimension];

        for(size_t i = 0; i < Dimension; i++)
        {
            MyMinPoint[i] = mMinPoint[i];
            MyMaxPoint[i] = mMaxPoint[i];
        }

        mpi_connectivity = vector<int>(mpi_size);
        mpi_MinPoints = vector<vector<double> >(mpi_size, vector<double>(Dimension));
        mpi_MaxPoints = vector<vector<double> >(mpi_size, vector<double>(Dimension));

        MPI_Allgather(MyMinPoint,Dimension,MPI_DOUBLE,MpiMinPoints,Dimension,MPI_DOUBLE,MPI_COMM_WORLD);
        MPI_Allgather(MyMaxPoint,Dimension,MPI_DOUBLE,MpiMaxPoints,Dimension,MPI_DOUBLE,MPI_COMM_WORLD);

        for(int i = 0; i < mpi_size; i++)
        {
            mpi_connectivity[i] = 0;

            for(size_t j = 0; j < Dimension; j++)
            {
                mpi_MinPoints[i][j] = MpiMinPoints[i * Dimension + j];
                mpi_MaxPoints[i][j] = MpiMaxPoints[i * Dimension + j];
            }
        }

        delete [] MpiMinPoints;
        delete [] MpiMaxPoints;
        delete [] MyMinPoint;
        delete [] MyMaxPoint;
    }

    void PrepareCommunications(int * NumberOfSendElements,int * NumberOfRecvElements,int * msgSendSize,int * msgRecvSize)
    {
        MPI_Alltoall(msgSendSize,1,MPI_INT,msgRecvSize,1,MPI_INT,MPI_COMM_WORLD);
        MPI_Alltoall(NumberOfSendElements,1,MPI_INT,NumberOfRecvElements,1,MPI_INT,MPI_COMM_WORLD);
    }

    void AsyncSendAndRecive(std::string messages[],int * msgSendSize,int * msgRecvSize)
    {
        int NumberOfCommunicationEvents = 0;
        int NumberOfCommunicationEventsIndex = 0;

        char * recvBuffers[mpi_size];

        for(int j = 0; j < mpi_size; j++)
        {
            if(j != mpi_rank && msgRecvSize[j]) NumberOfCommunicationEvents++;
            if(j != mpi_rank && msgSendSize[j]) NumberOfCommunicationEvents++;
        }

        MPI_Request * reqs = new MPI_Request[NumberOfCommunicationEvents];
        MPI_Status * stats = new MPI_Status[NumberOfCommunicationEvents];

        //Set up all receive and send events
        for(int i = 0; i < mpi_size; i++)
        {
            if(i != mpi_rank && msgRecvSize[i])
            {
                recvBuffers[i] = (char *)malloc(sizeof(char) * msgRecvSize[i]);

                MPI_Irecv(recvBuffers[i],msgRecvSize[i],MPI_CHAR,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
            }

            if(i != mpi_rank && msgSendSize[i])
            {
                char * mpi_send_buffer = (char *)malloc(sizeof(char) * msgSendSize[i]);

                memcpy(mpi_send_buffer,messages[i].c_str(),msgSendSize[i]);
                MPI_Isend(mpi_send_buffer,msgSendSize[i],MPI_CHAR,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
            }
        }

        //wait untill all communications finish
        MPI_Waitall(NumberOfCommunicationEvents, reqs, stats);

        for(int i = 0; i < mpi_size; i++)
        {
            if(i != mpi_rank && msgRecvSize[i])
              messages[i] = std::string(recvBuffers[i],msgRecvSize[i]);
        }

        delete [] reqs;
        delete [] stats;
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
        for(SizeType iDim = Dimension-1 ; iDim > 0 ; iDim--){
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
        for(SizeType iDim = Dimension-1 ; iDim > 0 ; iDim--){
            Index += ThisIndex[iDim];
            Index *= mN[iDim-1];
        }
        Index += ThisIndex[0];
        return Index;
      }

      //************************************************************************

      CellType CalculateCell( PointType const& ThisPoint ){
        CellType Cell;
        for(SizeType i = 0 ; i < Dimension ; i++)
            Cell[i] = CalculatePosition(ThisPoint[i],i);
        return Cell;
      }

      CellType CalculateCell( PointType const& ThisPoint, CoordinateType Radius ){
        CellType Cell;
        for(SizeType i = 0 ; i < Dimension ; i++)
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

        char * mpi_send_buffer = new char[(msgRecvSize+1)];
        char * mpi_recv_buffer = new char[(msgRecvSize+1) * mpi_size];

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

              remoteThisPoint[i]    = new PointType();
              remoteNearest[i]       = new PointType();
              remoteDistance[i]     = static_cast<CoordinateType>(DBL_MAX);

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

        Nearest     = resultNearest[0];
        Distance     = resultDistance[0];
        Found         = resultFound[0];

        for(int i = 1; i < mpi_size; i++)
        {
            if(resultFound[i] && resultDistance[i] < Distance)
            {
                Nearest     = resultNearest[0];
                Distance     = resultDistance[0];
                Found         = resultFound[0];
            }
        }

        ResultNearest = this->NullPointer();

        if(Found)
            ResultNearest = Nearest;

        delete [] mpi_send_buffer;
        delete [] mpi_recv_buffer;
    }

    //************************************************************************

    ///////////////////////////////////////////////////////////////////////////
    // MPI Single Input Search
    ///////////////////////////////////////////////////////////////////////////

    void MPISingleSearchInRadiusTest(const SizeType& NumberOfPoints, const SizeType& MaxNumberOfResults, const double& Radius, const SizeType& times)
    {
        PointerType PointInput = new PointType[NumberOfPoints];

        for(int i = 0; i < NumberOfPoints; i++)
        {
            PointType temp;

            temp[0] = (i+1)/NumberOfPoints;
            temp[1] = (i+1)/NumberOfPoints;
            temp[2] = 0;

            PointInput[i] = PointType(temp);
        }

        std::vector<SizeType>                  NumberOfResults(NumberOfPoints);
        std::vector<std::vector<PointerType> > Results(NumberOfPoints, std::vector<PointerType>(MaxNumberOfResults));
        std::vector<std::vector<double> >      ResultsDistances(NumberOfPoints, std::vector<double>(MaxNumberOfResults,0));

        MPI_SearchInRadius(&PointInput[NumberOfPoints/2], Radius, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    void MPI_SearchInRadius( PointerType const& ThisPoints, CoordinateType const& Radius, std::vector<std::vector<PointerType> >& Results,
        std::vector<std::vector<double> >& ResultsDistances, std::vector<SizeType>& NumberOfResults, SizeType const& MaxNumberOfResults)
    {
        CoordinateType Radius2 = Radius * Radius;
        SearchInRadiusMpiWrapperSingle( ThisPoints, 1, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults );
    }

    //************************************************************************

//         SizeType MPI_SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, IteratorType Results,
//             DistanceIteratorType ResultsDistances, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
//         {
//             CoordinateType Radius2 = Radius * Radius;
//             SizeType NumberOfResults = 0;
//             Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
//             SearchInRadiusMpiWrapper( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box );
//             return NumberOfResults;
//         }

    //************************************************************************

//         SizeType MPI_SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
//             DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults )
//         {
//             SearchStructureType Box( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
//             SearchInRadiusMpiWrapper( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
//         }

    //************************************************************************

//         SizeType MPI_SearchInRadius( PointType const& ThisPoint, CoordinateType const& Radius, CoordinateType const& Radius2, IteratorType& Results,
//             DistanceIteratorType& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults, SearchStructureType& Box )
//         {
//             Box.Set( CalculateCell(ThisPoint,-Radius), CalculateCell(ThisPoint,Radius), mN );
//             SearchInRadiusMpiWrapper( ThisPoint, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults, Box);
//         }

    ///////////////////////////////////////////////////////////////////////////
    // MPI Single Input END
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // MPI Multiple Input Search
    ///////////////////////////////////////////////////////////////////////////

    void MPIMultiSearchInRadiusTest(const SizeType& NumberOfPoints, const SizeType& MaxNumberOfResults, const double& Radius, const SizeType& times)
    {
        //MultiSearch Test
        Timer::Start("ALL");
        PointerType PointInput = new PointType[NumberOfPoints];

        for (ModelPart::ElementsContainerType::iterator el_it = StaticMesh->ElementsBegin();
            el_it != StaticMesh->ElementsEnd() && el_it-StaticMesh->ElementsBegin() < NumberOfPoints; el_it++)
        {
            PointType temp;
            Geometry<Node < 3 > >& geom = el_it->GetGeometry();

            temp[0] = (geom[0].X() + geom[1].X() + geom[2].X() + geom[3].X())/4;
            temp[1] = (geom[0].Y() + geom[1].Y() + geom[2].Y() + geom[3].Y())/4;
            temp[2] = (geom[0].Z() + geom[1].Z() + geom[2].Z() + geom[3].Z())/4;

            PointInput[el_it-StaticMesh->ElementsBegin()] = PointType(temp);
        }

        int rest = 0;
        for(size_t i = 0; i < times; i++)
        {
            std::vector<SizeType>                  NumberOfResults(NumberOfPoints);
            std::vector<std::vector<PointerType> > Results(NumberOfPoints, std::vector<PointerType>(MaxNumberOfResults));
            std::vector<std::vector<double> >      ResultsDistances(NumberOfPoints, std::vector<double>(MaxNumberOfResults,0));

            MPI_SearchInRadius(PointInput, NumberOfPoints, Radius, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults);
            MPI_Barrier(MPI_COMM_WORLD);

            if(i == times-1)
            {
                size_t max = 0;
                size_t min = MaxNumberOfResults;
                //Check Results
                for(size_t i = 0; i < NumberOfPoints; i++)
                {
                    rest += NumberOfResults[i];
                    max = NumberOfResults[i] > max ? NumberOfResults[i] : max;
                    min = NumberOfResults[i] < min ? NumberOfResults[i] : min;
                }
                std::cout << "(" << mpi_rank << ") Found: " << rest << " results,  aprox " << rest/NumberOfPoints << " results per point. Max: " << max << " Min: " << min << std::endl;
            }
        }

        int total_size = 0;
        int total_resu = 0;
        int local_size = NumberOfPoints;
        int local_resu = rest;

        MPI_Allreduce(&local_size,&total_size,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&local_resu,&total_resu,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

        if(mpi_rank == 0)
          std::cout << "Total Point search: " << total_size << " " << total_resu << std::endl;

        Timer::Stop("ALL");
    }

    void MPI_SearchInRadius( PointerType const& ThisPoints, SizeType const& NumberOfPoints, CoordinateType const& Radius, std::vector<std::vector<PointerType> >& Results,
        std::vector<std::vector<double> >& ResultsDistances, std::vector<SizeType>& NumberOfResults, SizeType const& MaxNumberOfResults)
    {
        CoordinateType Radius2 = Radius * Radius;

        SearchInRadiusMpiWrapperSingle( ThisPoints, NumberOfPoints, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults );
    }

    /// Act as wrapper between external function and its implementation
    /**
      * This function provides all mpi functionality requiered to execute the parallel multi input searchInRaidus.
      * the method implemented by this function is one-to-many. It means all particles not found in the local
      * processes are send to all process intersecting the search radius of the particle
      * @param ThisPoints List of points to be search
      * @param NumberOfPoints Number of points to be search
      * @param Radius Radius of search
      * @param Radius2 Radius of search^2
      * @param Results List of results
      * @param ResultsDistances Distance of the results
      * @param NumberOfResults Number of results
      * @param MaxNumberOfResults Maximum number of results returned for each point
      */
    void SearchInRadiusMpiWrapperSingle( PointerType const& ThisPoints, SizeType const& NumberOfPoints, CoordinateType const& Radius, CoordinateType const& Radius2, std::vector<std::vector<PointerType> >& Results,
          std::vector<std::vector<double> >& ResultsDistances, std::vector<SizeType>& NumberOfResults, SizeType const& MaxNumberOfResults )
    {
        std::vector<std::vector<PointerType> >  remoteResults(mpi_size, std::vector<PointerType>(0));
        std::vector<std::vector<PointerType> >  SearchPetitions(mpi_size, std::vector<PointerType>(0));
        std::vector<std::vector<PointerType> >  SearchResults(mpi_size, std::vector<PointerType>(0));
        std::vector<std::vector<PointerType> >  SendPointToProcess(mpi_size, std::vector<PointerType>(0));

        std::string messages[mpi_size];

        int NumberOfSendPoints[mpi_size];
        int NumberOfRecvPoints[mpi_size];

        std::vector<bool> SendPoint(NumberOfPoints*mpi_size);

        int msgSendSize[mpi_size];
        int msgRecvSize[mpi_size];

        PointerType CommunicationToken = new PointType();
        CommunicationToken->X() = std::numeric_limits<double>::max();

        for(int i = 0; i < mpi_size; i++)
        {
            NumberOfSendPoints[i] = 0;
            msgSendSize[i] = 0;
        }

        //Local search
        Timer::Start("Calculate Local");
        for(size_t i = 0; i < NumberOfPoints; i++)
        {
            IteratorType ResultsPointer      = &Results[i][0];
            double * ResultsDistancesPointer = &ResultsDistances[i][0];

            NumberOfResults[i] = 0;

            SearchStructureType Box( CalculateCell(ThisPoints[i],-Radius), CalculateCell(ThisPoints[i],Radius), mN );
            SearchInRadiusLocal(ThisPoints[i],Radius,Radius2,ResultsPointer,ResultsDistancesPointer,NumberOfResults[i],MaxNumberOfResults,Box);

            //For each point with results < MaxResults and each process excluding ourself
            if(NumberOfResults[i] < MaxNumberOfResults)
            {
                for(int j = 0; j < mpi_size; j++)
                {
                    if(j != mpi_rank)
                    {
                        int intersect = 0;
                        for(size_t k = 0; k < Dimension; k++)
                            if((ThisPoints[i][k]+Radius >= mpi_MaxPoints[j][k] && ThisPoints[i][k]-Radius <= mpi_MaxPoints[j][k]) ||
                               (ThisPoints[i][k]+Radius >= mpi_MinPoints[j][k] && ThisPoints[i][k]-Radius <= mpi_MinPoints[j][k]) ||
                               (ThisPoints[i][k]-Radius >= mpi_MinPoints[j][k] && ThisPoints[i][k]+Radius <= mpi_MaxPoints[j][k])
                            ) intersect++;

                        SendPoint[j*NumberOfPoints+i] = 0;
                        if(intersect == Dimension)
                        {
                            SendPoint[j*NumberOfPoints+i]=1;
                            NumberOfSendPoints[j]++;
                        }
                    }
                }
            }
        }

        for(int i = 0; i < mpi_size; i++)
        {
            if(i != mpi_rank && NumberOfSendPoints[i])
            {
                int k = 0;
                SendPointToProcess[i].resize(NumberOfSendPoints[i]);
                for(size_t j = 0; j < NumberOfPoints; j++)
                    if(SendPoint[i*NumberOfPoints+j])
                        SendPointToProcess[i][k++] = &ThisPoints[j];
            }
        }
        Timer::Stop("Calculate Local");

        Timer::Start("Transfer Particles");
//         for(int i = 0; i < mpi_size; i++)
//         {
//             if(mpi_rank != i)
//                 TConfigure::Save(SendPointToProcess[i],messages[i]);
//             msgSendSize[i] = messages[i].size();
//         }
//
//         PrepareCommunications(msgSendSize,msgRecvSize,NumberOfSendPoints,NumberOfRecvPoints);
//         AsyncSendAndRecive(messages,msgSendSize,msgRecvSize);
//
//         for(int i = 0; i < mpi_size; i++)
//         {
//             if(mpi_rank != i && messages[i].size())
//                 TConfigure::Load(SearchPetitions[i],messages[i]);
//         }

        TConfigure::AsyncSendAndRecive(SendPointToProcess,SearchPetitions,msgSendSize,msgRecvSize);
        Timer::Stop("Transfer Particles");

        Timer::Start("Calculate Remote");
        //Calculate remote points
        for(int i = 0; i < mpi_size; i++)
        {
            if(i != mpi_rank && msgRecvSize[i])
            {
                int accum_results = 0;
                NumberOfRecvPoints[i] = SearchPetitions[i].size();
                std::vector<PointerType>& remoteSearchPetitions = SearchPetitions[i];
                remoteResults[i].resize((MaxNumberOfResults+1)*NumberOfRecvPoints[i]);
                for(int j = 0; j < NumberOfRecvPoints[i]; j++)
                {
                    IteratorType remoteResultsPointer      = &remoteResults[i][accum_results];
                    PointType remotePointPointer           = *remoteSearchPetitions[j];
                    SizeType thisNumberOfResults = 0;

                    SearchStructureType Box( CalculateCell(remotePointPointer,-Radius), CalculateCell(remotePointPointer,Radius), mN );
                    SearchInRadiusLocal(remotePointPointer,Radius,Radius2,remoteResultsPointer,thisNumberOfResults,MaxNumberOfResults,Box);
                    accum_results += thisNumberOfResults;
                    remoteResults[i][accum_results++] = CommunicationToken;
                }
                remoteResults[i].resize(accum_results);
                NumberOfSendPoints[i] = accum_results;
            }
        }
        Timer::Stop("Calculate Remote");


        Timer::Start("Transfer Results");
        for(int i = 0; i < mpi_size; i++)
        {
            if(mpi_rank != i)
              TConfigure::Save(remoteResults[i],messages[i]);
            msgSendSize[i] = messages[i].size();
        }

        PrepareCommunications(msgSendSize,msgRecvSize,NumberOfSendPoints,NumberOfRecvPoints);
        AsyncSendAndRecive(messages,msgSendSize,msgRecvSize);

        for(int i = 0; i < mpi_size; i++)
        {
            if(mpi_rank != i && messages[i].size())
                TConfigure::Load(SearchResults[i],messages[i]);
        }
        Timer::Stop("Transfer Results");

        Timer::Start("Prepare-C");
        for (int i = 0; i < mpi_size; i++)
        {
            if (i != mpi_rank)
            {
                std::vector<PointerType>& remoteSearchResults = SearchResults[i];

                int result_iterator = 0;
                for(size_t j = 0; j < NumberOfPoints; j++)
                {
                    if(SendPoint[i*NumberOfPoints+j])
                    {
                        int token = 0;

                        for(; !token && result_iterator < NumberOfRecvPoints[i]; result_iterator++)
                        {
                            PointType& a = ThisPoints[j];
                            PointType& b = *remoteSearchResults[result_iterator];

                            if(b.X() == std::numeric_limits<double>::max())
                              token = 1;

                            if (!token)
                            {
                                double dist = DistanceFunction()(a,b);

                                if (dist <= Radius2)
                                {
                                    if (NumberOfResults[j] < MaxNumberOfResults)
                                    {
                                        Results[j][NumberOfResults[j]] = remoteSearchResults[result_iterator];

                                        ResultsDistances[j][NumberOfResults[j]] = dist;
                                        NumberOfResults[j]++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        Timer::Stop("Prepare-C");
    }

    ///////////////////////////////////////////////////////////////////////////
    // MPI Search In Radius
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

    void MultiSearchInRadiusTest(const SizeType& NumberOfPoints, const SizeType& MaxNumberOfResults, const double& Radius, const SizeType& times)
    {
        //MultiSearch Test
        Timer::Start("ALL");
        PointerType  PointInput = new PointType[NumberOfPoints];
//             std::vector<IteratorType Results>(NumberOfPoints * MaxNumberOfResults);

//             std::vector<std::vector<PointerType> > Results(NumberOfPoints, std::vector<PointerType>(MaxNumberOfResults));
//             std::vector<std::vector<double> >      ResultsDistances(NumberOfPoints, std::vector<double>(MaxNumberOfResults,0));

        for (ModelPart::ElementsContainerType::iterator el_it = StaticMesh->ElementsBegin();
            el_it != StaticMesh->ElementsEnd() && el_it-StaticMesh->ElementsBegin() < NumberOfPoints; el_it++)
        {
            PointType temp;
            Geometry<Node < 3 > >& geom = el_it->GetGeometry();

            temp[0] = (geom[0].X() + geom[1].X() + geom[2].X() + geom[3].X())/4;
            temp[1] = (geom[0].Y() + geom[1].Y() + geom[2].Y() + geom[3].Y())/4;
            temp[2] = (geom[0].Z() + geom[1].Z() + geom[2].Z() + geom[3].Z())/4;

            PointInput[el_it-StaticMesh->ElementsBegin()] = PointType(temp);
        }

        int rest = 0;
        #pragma omp parallel for reduction(+:rest)
        for(size_t i = 0; i < NumberOfPoints; i++)
        {
            IteratorType Results    = new PointerType[MaxNumberOfResults];
            DistanceIteratorType Distances = new double[MaxNumberOfResults];
            PointType mypoint = PointInput[i];
//                 IteratorType
            rest += SearchInRadius(mypoint, Radius, Results, Distances, MaxNumberOfResults);
        }
        std::cout << "Rest: " << rest << std::endl;
        Timer::Stop("ALL");
    }

    SizeType SearchInRadius( PointerType const& ThisPoints, SizeType const& NumberOfPoints, CoordinateType const& Radius, std::vector<std::vector<PointerType> > Results,
        std::vector<std::vector<double> >  ResultsDistances, SizeType const& MaxNumberOfResults )
    {
        CoordinateType Radius2 = Radius * Radius;
        SizeType NumberOfResults = 0;
        /*SearchStructureType Box[NumberOfPoints]*/;

//             for(size_t i = 0; i < NumberOfPoints; i++)
//                 Box[i] = SearchStructureType( CalculateCell(ThisPoints[i],-Radius), CalculateCell(ThisPoints[i],Radius), mN );
//
        SearchInRadiusLocal( ThisPoints, NumberOfPoints, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults/*, Box*/ );

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
        std::vector<std::vector<PointerType> >& Results, std::vector<std::vector<double> >& ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults /*, SearchStructureType Box[]*/ )
    {
//             for(size_t i = 0; i < NumberOfPoints; i++)
//                 Box[i].Set( CalculateCell(ThisPoints[i],-Radius), CalculateCell(ThisPoints[i],Radius), mN );
        SearchInRadiusLocal( ThisPoints, NumberOfPoints, Radius, Radius2, Results, ResultsDistances, NumberOfResults, MaxNumberOfResults);
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
        std::vector<std::vector<PointerType> > Results, std::vector<std::vector<double> > ResultsDistances, SizeType& NumberOfResults, SizeType const& MaxNumberOfResults
        /*, SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3> * const& Box*/ )
    {
        int ompSafe_NumberOfResults = 0;
//             #pragma omp parallel for reduction(+:ompSafe_NumberOfResults)
        for(size_t i = 0; i < NumberOfPoints; i++)
        {
            SizeType thisNumberOfResults = 0;
            IteratorType ResulPointer = &Results[i][0];
            DistanceIteratorType ResultsDistancesPointer = &ResultsDistances[i][0];
            SearchStructureType Box;
            Box.Set( CalculateCell(ThisPoint[i],-Radius), CalculateCell(ThisPoint[i],Radius), mN );
            for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
                for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
                    for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                    {
                        SearchRadiusInRange()(mPoints[I].begin(),mPoints[I].end(),ThisPoint[i],Radius2,ResulPointer,ResultsDistancesPointer,thisNumberOfResults,MaxNumberOfResults);
                    }
            ompSafe_NumberOfResults += thisNumberOfResults;
        }
        NumberOfResults = ompSafe_NumberOfResults;
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
      for(SizeType i = 0 ; i < Dimension ; i++)
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

    Tvector<CoordinateType,Dimension>  mMinPoint;
    Tvector<CoordinateType,Dimension>  mMaxPoint;
    Tvector<CoordinateType,Dimension>  mCellSize;
    Tvector<CoordinateType,Dimension>  mInvCellSize;
    Tvector<SizeType,Dimension>        mN;
    SizeType                            mNumPoints;

    ModelPart * StaticMesh;

    // Bins Access Vector ( vector<Iterator> )
    CellsContainerType mPoints;

    // Work Variables ( For non-copy of Search Variables )
    //BinBox SearchBox;

    //MPI interface
    int mpi_rank;
    int mpi_size;

    //MPI Communication
    vector<int> mpi_connectivity;
    vector<vector<double> > mpi_MinPoints;
    vector<vector<double> > mpi_MaxPoints;

public:
//     static TreeNodeType* Construct(IteratorType PointsBegin, IteratorType PointsEnd, PointType MaxPoint, PointType MinPoint, SizeType BucketSize)
//     {
//
//       SizeType number_of_points = SearchUtils::PointerDistance(PointsBegin,PointsEnd);
//       if (number_of_points == 0)
//         return NULL;
//       else
//       {
//         return new BinsDynamicMpi( PointsBegin, PointsEnd, MinPoint, MaxPoint, BucketSize );
//       }
//
//     }
};

template<class TConfigure>
std::ostream & operator<<( std::ostream& rOStream,
      BinsDynamicMpi<TConfigure>& rThis)
{
   rThis.PrintInfo(rOStream);
   rOStream << std::endl;
   rThis.PrintSize(rOStream);
   rThis.PrintData(rOStream);
   return rOStream;
}



}

#endif // KRATOS_BINS_DYNAMIC_MPI_CONTAINER_H_INCLUD
