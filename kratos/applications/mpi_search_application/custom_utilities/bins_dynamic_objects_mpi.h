/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/


//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Nelson Lafontaine $
//   Date:                $Date: 29-09-2010 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_BINS_DYNAMIC_OBJECTS_MPI_CONTAINER_H_INCLUDED)
#define  KRATOS_BINS_DYNAMIC_OBJECTS_MPI_CONTAINER_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <time.h>

// Project includes
#include "mpi.h"
#include "spatial_containers/tree.h"
#include "spatial_containers/cell.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define CUSTOMTIMER 1  // ACTIVATES AND DISABLES ::TIMER:::::

/* Timer defines */
#ifdef CUSTOMTIMER
#define KRATOS_TIMER_START(t) Timer::Start(t);
#define KRATOS_TIMER_STOP(t) Timer::Stop(t);
#else
#define KRATOS_TIMER_START(t)
#define KRATOS_TIMER_STOP(t)
#endif

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{


/// Short class definition.
/** Detail class definition.
*/
template<class TConfigure>
class BinsObjectDynamicMpi
{
public:
    ///@name Type Definitions
    ///@{

    enum { Dimension = TConfigure::Dimension };

    //Point
    typedef TConfigure                                      Configure;
    typedef typename TConfigure::PointType                  PointType;

    //Container
    typedef typename TConfigure::PointerType                PointerType;
    typedef typename TConfigure::ContainerType              ContainerType;
    typedef typename TConfigure::IteratorType               IteratorType;
    typedef typename TConfigure::DistanceIteratorType       DistanceIteratorType;
    typedef typename TConfigure::ResultContainerType        ResultContainerType;
    typedef typename TConfigure::ElementsContainerType      ElementsContainerType;
    typedef typename TConfigure::ResultIteratorType         ResultIteratorType;
    typedef typename TConfigure::PointerContactType         PointerContactType;

    //Search Structures
    typedef Cell<Configure> CellType;
    typedef std::vector<CellType> CellContainerType;
    typedef std::vector<int> DomainContainerType;
    typedef typename CellContainerType::iterator CellContainerIterator;

    typedef TreeNode<Dimension, PointType, PointerType, IteratorType,  typename TConfigure::DistanceIteratorType> TreeNodeType;
    typedef typename TreeNodeType::CoordinateType       CoordinateType;  // double
    typedef typename TreeNodeType::SizeType             SizeType;        // std::size_t
    typedef typename TreeNodeType::IndexType            IndexType;       // std::size_t


    typedef Tvector<IndexType,Dimension>                IndexArray;
    typedef Tvector<SizeType,Dimension>                 SizeArray;
    typedef Tvector<CoordinateType,Dimension>           CoordinateArray;

    ///Contact Pair
    typedef typename TConfigure::ContainerContactType   ContainerContactType;
    typedef typename TConfigure::IteratorContactType    IteratorContactType;

    ///typedef TreeNodeType LeafType;
    typedef typename TreeNodeType::IteratorIteratorType IteratorIteratorType;
    typedef typename TreeNodeType::SearchStructureType  SearchStructureType;

    /// Pointer definition of BinsObjectDynamic
    KRATOS_CLASS_POINTER_DEFINITION(BinsObjectDynamicMpi);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BinsObjectDynamicMpi() {}
    /// Constructor de bins a bounding box


    BinsObjectDynamicMpi (IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd)
        : mObjectsBegin(ObjectsBegin), mObjectsEnd(ObjectsEnd)
    {
        mObjectsSize = SearchUtils::PointerDistance(mObjectsBegin,mObjectsEnd);

        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        BinsExtension = 2;

        CalculateBoundingBox();             // Calculate mMinPoint, mMaxPoint
        GenerateCommunicationGraph();       // Share The min and max
        CalculateCellSize(mObjectsSize);    // Calculate number of Cells
        AllocateContainer();                // Allocate cell list
        GenerateBins();                     // Fill Cells with objects

    }

    BinsObjectDynamicMpi (IteratorType const& ObjectsBegin, IteratorType const& ObjectsEnd, CoordinateType CellSize)
        : mObjectsBegin(ObjectsBegin), mObjectsEnd(ObjectsEnd)
    {
        mObjectsSize = SearchUtils::PointerDistance(mObjectsBegin,mObjectsEnd);
/*
        CommunicationToken = PointerType(ObjectsBegin[0]);
        CommunicationToken->GetValue(NUMBER_OF_NEIGHBOURS) = -2;*/

        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        BinsExtension = 2;

        CalculateBoundingBox();             // Calculate mMinPoint, mMaxPoint
        GenerateCommunicationGraph();
        CalculateCellSize(CellSize);        // Calculate number of Cells
        AllocateContainer();                // Allocate cell list
        GenerateBins();                     // Fill Cells with objects

    }

    BinsObjectDynamicMpi (const PointType& MinPoint, const PointType& MaxPoint, CoordinateType CellSize)
        : mObjectsSize(0), mObjectsBegin(0), mObjectsEnd(0)
    {

        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

//         CommunicationToken = PointerType(ObjectsBegin[0]);
//         CommunicationToken->GetValue(NUMBER_OF_NEIGHBOURS) = -1;

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            mMinPoint[i] = MinPoint[i];
            mMaxPoint[i] = MaxPoint[i];
        }
        CalculateCellSize(CellSize);
        AllocateContainer();
    }

    BinsObjectDynamicMpi (const PointType& MinPoint, const PointType& MaxPoint, SizeType NumPoints)
        : mObjectsSize(0), mObjectsBegin(0), mObjectsEnd(0)
    {

        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

//         CommunicationToken = PointerType(ObjectsBegin[0]);
//         CommunicationToken->GetValue(NUMBER_OF_NEIGHBOURS) = -1;

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            mMinPoint[i] = MinPoint[i];
            mMaxPoint[i] = MaxPoint[i];
        }
        CalculateCellSize(NumPoints);
        AllocateContainer();
    }


    /// Destructor.
    virtual ~BinsObjectDynamicMpi() {
    }


//************************************************************************
//************************************************************************
    SizeType SearchObjectsInCell(const PointType& ThisPoint, ResultIteratorType Result, const SizeType& MaxNumberOfResults)
    {
        IndexType icell = CalculateIndex(ThisPoint);

        /*for(IndexType I = 0; I< Dimension; I++)*/
        if(mCells[icell].Size() < MaxNumberOfResults)
        {
            for(IteratorType i_object = mCells[icell].Begin() ; i_object != mCells[icell].End(); i_object++, Result++)
                *Result = *i_object;
            return mCells[icell].Size();
        }
        else
            return -1;

    }
//************************************************************************
//************************************************************************


    SizeType SearchObjects(PointerType& ThisObject, ResultIteratorType& Result,  const SizeType& MaxNumberOfResults )
    {
        PointType Low, High;
        SearchStructureType Box;
        SizeType NumberOfResults = 0;
        TConfigure::CalculateBoundingBox(ThisObject, Low, High);
        Box.Set( CalculateCell(Low), CalculateCell(High), mN );
        SearchInBoxLocal(ThisObject, Result, NumberOfResults, MaxNumberOfResults, Box );
        return NumberOfResults;
    }


//************************************************************************
//************************************************************************

    SizeType SearchObjects(PointerType& ThisObject, ResultContainerType& Result)
    {
        PointType Low, High;
        SearchStructureType Box;
        TConfigure::CalculateBoundingBox(ThisObject, Low, High);
        Box.Set( CalculateCell(Low), CalculateCell(High), mN );
        SearchInBoxLocal(ThisObject, Result, Box );
        return Result.size();
    }


//************************************************************************
//************************************************************************

   /// Act as a wrapper between external function and its implementation
    /**
      * This function provides all mpi functionality requiered to execute the parallel multi input searchInRaidus.
      * The method implemented by this function is one-to-one what means that all particles not found in the local
      * processes are send only to the processes intersecting the search radius of the particle
      * @param ThisObjects List of objects to be search
      * @param NumberOfObjects Number of points to be search
      * @param Radius Radius of search
      * @param Radius2 Radius of search^2
      * @param Results List of results
      * @param ResultsDistances Distance of the results
      * @param NumberOfResults Number of results
      * @param MaxNumberOfResults Maximum number of results returned for each point
      **/
    void SearchObjectsMpi(const ElementsContainerType & ThisObjects, SizeType const& NumberOfObjects, std::vector<double> const& Radius, std::vector<std::vector<PointerType> >& Results,
          std::vector<std::vector<double> >& ResultsDistances, std::vector<SizeType>& NumberOfResults, SizeType const& MaxNumberOfResults, Communicator& Communicator, bool useRealData = true)
    {
        std::vector<ElementsContainerType> remoteResults(mpi_size);
        std::vector<ElementsContainerType> SearchPetitions(mpi_size);
        std::vector<ElementsContainerType> SearchResults(mpi_size);
        std::vector<ElementsContainerType> SendObjectToProcess(mpi_size);

        std::vector<std::vector<double> >  SendRadiusToProcess(mpi_size, std::vector<double>(0));
        std::vector<std::vector<double> >  SearchPetitionsRadius(mpi_size, std::vector<double>(0));
        std::vector<std::vector<double> >  SendResultsPerPoint(mpi_size, std::vector<double>(0));
        std::vector<std::vector<double> >  RecvResultsPerPoint(mpi_size, std::vector<double>(0));

        int * NumberOfSendPoints = new int[mpi_size];
        int * NumberOfRecvPoints = new int[mpi_size];

        std::vector<bool> SendPoint(NumberOfObjects*mpi_size);

        PointType Low, High;
        SearchStructureType Box;

        for(int i = 0; i < mpi_size; i++)
        {
            NumberOfSendPoints[i] = 0;
        }

        IteratorType it_begin = const_cast<typename ElementsContainerType::ContainerType& >(ThisObjects.GetContainer()).begin();
        IteratorType it_end   = const_cast<typename ElementsContainerType::ContainerType& >(ThisObjects.GetContainer()).end();

        int objectCounter = 0;

        for (IteratorType object_pointer_it = it_begin; object_pointer_it != it_end; ++object_pointer_it)
        {
            ResultIteratorType   ResultsPointer          = Results[objectCounter].begin();
            DistanceIteratorType ResultsDistancesPointer = ResultsDistances[objectCounter].begin();

            NumberOfResults[objectCounter] = 0;

            TConfigure::CalculateBoundingBox((*object_pointer_it), Low, High, Radius[objectCounter]);
            Box.Set( CalculateCell(Low), CalculateCell(High), mN );

            SearchInRadiusExclusive((*object_pointer_it), Radius[objectCounter], ResultsPointer, ResultsDistancesPointer, NumberOfResults[objectCounter], MaxNumberOfResults, Box );

            //For each point with results < MaxResults and each process excluding ourself
            if(NumberOfResults[objectCounter] < MaxNumberOfResults)
            {
                for(int j = 0; j < mpi_size; j++)
                {
                    if((((*object_pointer_it)->GetGeometry()[0].FastGetSolutionStepValue(PARTITION_MASK)) & (1<<j)))
                    //if(true)
                    {
                        SendPoint[j*NumberOfObjects+objectCounter]=1;
                        NumberOfSendPoints[j]++;
                    }
                }
            }

            objectCounter++;
        }

        for(int i = 0; i < mpi_size; i++)
        {
            if(i != mpi_rank && NumberOfSendPoints[i])
            {
                int k = 0;

                SendObjectToProcess[i].reserve(NumberOfSendPoints[i]);
                SendRadiusToProcess[i].resize(NumberOfSendPoints[i]);

                for(size_t j = 0; j < NumberOfObjects; j++)
                {
                    if( SendPoint[i*NumberOfObjects+j])
                    {
                        IteratorType itrObject = const_cast<typename ElementsContainerType::ContainerType& >(ThisObjects.GetContainer()).begin() + j;

                        SendObjectToProcess[i].push_back(*itrObject);
                        SendRadiusToProcess[i][k] = Radius[j];

                        k++;
                    }
                }
            }
        }

        TConfigure::TransferObjects(Communicator,SendObjectToProcess,SearchPetitions);
        TConfigure::TransferObjects(SendRadiusToProcess,SearchPetitionsRadius);

        Communicator::NeighbourIndicesContainerType communicator_ranks = Communicator.NeighbourIndices();

        //Calculate remote points
        for(int i = 0; i < mpi_size; i++)
        {
            int NumberOfRanks = Communicator.GetNumberOfColors();
            if(i != mpi_rank)
            {
                int destination = -1;

                for(int j = 0; j < NumberOfRanks; j++)
                    if(i == communicator_ranks[j])
                        destination = j;

                int accum_results = 0;

                NumberOfRecvPoints[i] = SearchPetitions[i].size();
                remoteResults[i].reserve((MaxNumberOfResults+1)*NumberOfRecvPoints[i]);
                SendResultsPerPoint[i].resize(SearchPetitionsRadius[i].size());

                std::vector<double> TempResultsDistances(MaxNumberOfResults);

                for(int j = 0; j < NumberOfRecvPoints[i]; j++)
                {
                    DistanceIteratorType ResultsDistancesPointer = TempResultsDistances.begin(); //Useless
                    ResultContainerType  TempResults(MaxNumberOfResults);
                    ResultIteratorType   remoteResultsPointer = TempResults.begin();

                    SizeType thisNumberOfResults = 0;
                    *ResultsDistancesPointer = 0;

                    TConfigure::CalculateBoundingBox((SearchPetitions[i].GetContainer())[j], Low, High, SearchPetitionsRadius[i][j]);

                    Box.Set( CalculateCell(Low), CalculateCell(High), mN );

                    SearchInRadiusExclusive((SearchPetitions[i].GetContainer())[j], SearchPetitionsRadius[i][j], remoteResultsPointer, ResultsDistancesPointer, thisNumberOfResults, MaxNumberOfResults, Box );

                    for(ResultIteratorType result_it = TempResults.begin(); result_it != remoteResultsPointer; ++result_it)
                    {
                        Communicator.LocalMesh(destination).Elements().push_back((*result_it));
                        Communicator.LocalMesh(destination).Nodes().push_back((*result_it)->GetGeometry()(0));

                        (remoteResults[i].GetContainer()).push_back(*result_it);

                        accum_results++;
                    }

                    SendResultsPerPoint[i][j] = thisNumberOfResults;
                }

                NumberOfSendPoints[i] = accum_results;
            }
        }

        TConfigure::TransferObjects(Communicator,remoteResults,SearchResults);
        TConfigure::TransferObjects(SendResultsPerPoint,RecvResultsPerPoint);

        for(int i = 0; i < mpi_size; i++) //for all ranks
        {
            int NumberOfRanks = Communicator.GetNumberOfColors();
            if(i != mpi_rank) //not being myself
            {
                int ParticleCounter = 0;
                int ResultCounter = 0;

                int origin = -1;

                for(int j = 0; j < NumberOfRanks; j++)
                    if(i == communicator_ranks[j])
                        origin = j;

                for(size_t j = 0; j < NumberOfObjects; j++)
                {
                    if(SendPoint[i*NumberOfObjects+j])
                    {
                        for(size_t k = 0; k < RecvResultsPerPoint[i][ParticleCounter]; k++)
                        {
                            double dist;
                            IteratorType itrObject = const_cast<typename ElementsContainerType::ContainerType& >(ThisObjects.GetContainer()).begin() + j;

                            Results[j][NumberOfResults[j]] = (SearchResults[i].GetContainer())[ResultCounter];
                            TConfigure::Distance((*itrObject),(SearchResults[i].GetContainer())[ResultCounter],dist);
                            ResultsDistances[j][NumberOfResults[j]] = dist;
                            NumberOfResults[j]++;

                            Communicator.GhostMesh().Elements().push_back((SearchResults[i].GetContainer())[ResultCounter]);
                            Communicator.GhostMesh(origin).Elements().push_back((SearchResults[i].GetContainer())[ResultCounter]);

                            bool repeat = 0;
                            for(ModelPart::NodesContainerType::iterator it = Communicator.GhostMesh(origin).Nodes().begin(); !repeat && it != Communicator.GhostMesh(origin).Nodes().end(); ++it)
                              if((SearchResults[i].GetContainer())[ResultCounter]->GetGeometry()[0].Id() == it->Id())
                                repeat = 0;
                            if(!repeat)
                            {
                                Communicator.GhostMesh().Nodes().push_back((SearchResults[i].GetContainer())[ResultCounter]->GetGeometry()(0));
                                Communicator.GhostMesh(origin).Nodes().push_back((SearchResults[i].GetContainer())[ResultCounter]->GetGeometry()(0));
                            }

                            ResultCounter++;
                        }
                        ParticleCounter++;
                    }
                }
            }
        }

        delete [] NumberOfSendPoints;
        delete [] NumberOfRecvPoints;
    }

//************************************************************************
//************************************************************************

    SizeType SearchObjectsInner(PointerType& ThisObject, ResultIteratorType& Result,  const SizeType& MaxNumberOfResults )
    {
        PointType Low, High;
        SearchStructureType Box;
        SizeType NumberOfResults = 0;
        TConfigure::CalculateBoundingBox(ThisObject, Low, High);
        Box.Set( CalculateCell(Low), CalculateCell(High), mN );
        SearchObjectLocalInner(ThisObject, Result, NumberOfResults, MaxNumberOfResults, Box );
        return NumberOfResults;
    }


//************************************************************************
//************************************************************************

    SizeType SearchObjectsInner(PointerType& ThisObject, ResultContainerType& Result)
    {
        PointType Low, High;
        SearchStructureType Box;
        TConfigure::CalculateBoundingBox(ThisObject, Low, High);
        Box.Set( CalculateCell(Low), CalculateCell(High), mN );
        SearchObjectLocalInner(ThisObject, Result, Box );
        return Result.size();
    }

//************************************************************************
//************************************************************************

    SizeType SearchObjectsInRadius(PointerType& ThisObject, const double& Radius, ResultIteratorType& Results,
                                   DistanceIteratorType ResultDistances, const SizeType& MaxNumberOfResults)
    {
        PointType Low, High;
        SearchStructureType Box;
        SizeType NumberOfResults = 0;
        TConfigure::CalculateBoundingBox(ThisObject, Low, High, Radius);
        Box.Set( CalculateCell(Low), CalculateCell(High), mN );
        SearchInRadius(ThisObject, Radius, Results, ResultDistances, NumberOfResults, MaxNumberOfResults, Box );
        return NumberOfResults;
    }

//************************************************************************
//************************************************************************

    SizeType SearchObjectsInRadiusInner(PointerType& ThisObject, const double& Radius, ResultIteratorType& Results,
                                   DistanceIteratorType ResultDistances, const SizeType& MaxNumberOfResults)
    {
        PointType Low, High;
        SearchStructureType Box;
        SizeType NumberOfResults = 0;
        TConfigure::CalculateBoundingBox(ThisObject, Low, High, Radius);
        Box.Set( CalculateCell(Low), CalculateCell(High), mN );
        SearchInRadiusExclusive(ThisObject, Radius, Results, ResultDistances, NumberOfResults, MaxNumberOfResults, Box );
        return NumberOfResults;
    }

//************************************************************************
//************************************************************************

    void SearchContact(ContainerContactType& Result)
    {
        for (CellContainerIterator icell = mCells.begin() ; icell!= mCells.end(); icell++)
            icell->SearchContact(Result);
    }

//************************************************************************
//************************************************************************

    SizeType SearchContact(IteratorContactType& Result, const SizeType& MaxNumberOfResults )
    {
        SizeType NumberOfResults = 0;
        for (CellContainerIterator icell = mCells.begin() ; icell!= mCells.end(); icell++)
        {
            icell->SearchContact(Result, NumberOfResults, MaxNumberOfResults);
        }
        return NumberOfResults;
    }

//************************************************************************
//************************************************************************

    void AddObject(const PointerType& ThisObject)
    {
        PointType Low, High;
        SearchStructureType Box;
        TConfigure::CalculateBoundingBox(ThisObject, Low, High);
        Box.Set( CalculateCell(Low), CalculateCell(High), mN );
        FillObject(Box,ThisObject);
        mObjectsSize++;
    }

    void RemoveObject(const PointerType& ThisObject)
    {
        PointType Low, High;
        SearchStructureType Box;
        TConfigure::CalculateBoundingBox(ThisObject, Low, High);
        Box.Set( CalculateCell(Low), CalculateCell(High), mN );
        RemoveObjectLocal(Box,ThisObject);
        mObjectsSize--;
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "BinsObjectDynamicMpi" ;
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream, std::string const& Perfix = std::string()) const
    {
        rOStream << " BinsSize: ";
        for(SizeType i = 0 ; i < Dimension ; i++)
            rOStream << "[" << mN[i] << "]";
        rOStream << std::endl;
        rOStream << "  CellSize: ";
        for(SizeType i = 0 ; i < Dimension ; i++)
            rOStream << "[" << mCellSize[i] << "]";
        rOStream << std::endl;
        SizeType nn = 0;
        for(SizeType i = 0 ; i < mCells.size(); i++)
            nn += mCells[i].Size();
        rOStream << "NumPointers: " << nn << std::endl;
    }

    /// Print Size of Container
    void PrintSize( std::ostream& rout )
    {
        rout << " BinsSize: ";
        for(SizeType i = 0 ; i < Dimension ; i++)
            rout << "[" << mN[i] << "]";
        rout << std::endl;
    }

    /// Print Limits Points of the Container
    void PrintBox( std::ostream& rout )
    {
        rout << " BinsBox: Min [";
        mMinPoint.Print(rout);
        rout <<       "];  Max [";
        mMaxPoint.Print(rout);
        rout <<       "];  Size [";
        mCellSize.Print(rout);
        rout << "]" << std::endl;
    }


    CellContainerType& GetCellContainer()
    {
        return mCells;
    }

    SizeArray& GetDivisions()
    {
        return mN;
    }

    CoordinateArray& GetCellSize()
    {
        return mCellSize;
    }



//************************************************************************
//************************************************************************


private:


    ///@}
    ///@name Friends
    ///@{


    ///@}

    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{




    /// Computa los boxes de cada uno de los elementos del model part
    void CalculateBoundingBox()
    {

        PointType Low, High;
        TConfigure::CalculateBoundingBox(*mObjectsBegin,mMinPoint,mMaxPoint);

#ifdef _OPENMP
        SizeType number_of_threads = omp_get_max_threads();
#else
        SizeType number_of_threads = 1;
#endif

        std::vector<SizeType> node_partition;
        CreatePartition(number_of_threads, mObjectsSize, node_partition);

        std::vector<PointType> Max(number_of_threads);
        std::vector<PointType> Min(number_of_threads);

        for(SizeType k=0; k<number_of_threads; k++ )
        {
            Max[k] = mMaxPoint;
            Min[k] = mMinPoint;
        }

        #pragma omp parallel for  private(High, Low)
        for(int k=0; k<static_cast<int>(number_of_threads); k++)
        {
            IteratorType i_begin = mObjectsBegin + node_partition[k];
            IteratorType i_end   = mObjectsBegin + node_partition[k+1];

            for (IteratorType i_object  = i_begin ; i_object != i_end ; i_object++ )
            {
                TConfigure::CalculateBoundingBox(*i_object, Low, High);
                for(SizeType i = 0 ; i < Dimension ; i++)
                {
                    Max[k][i] = (Max[k][i]  < High[i]) ? High[i] : Max[k][i];
                    Min[k][i] = (Min[k][i]  > Low[i])  ? Low[i]  : Min[k][i];
                }
            }
        }

        for(SizeType k=0; k<number_of_threads; k++)
        {
            for(SizeType i = 0 ; i < Dimension ; i++)
            {
                mMaxPoint[i]  = (mMaxPoint[i]  < Max[k][i]) ? Max[k][i] : mMaxPoint[i];
                mMinPoint[i]  = (mMinPoint[i]  > Min[k][i]) ? Min[k][i] : mMinPoint[i];
            }
        }

        PointType Epsilon = mMaxPoint - mMinPoint;

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            mMaxPoint[i] += Epsilon[i] * 0.1;
            mMinPoint[i] -= Epsilon[i] * 0.1;
        }
    }


//************************************************************************
//************************************************************************

    void CalculateCellSize()
    {

        CoordinateType delta[Dimension];
        CoordinateType alpha[Dimension];
        CoordinateType mult_delta = 1.00;
        SizeType index = 0;
        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            delta[i] = mMaxPoint[i] - mMinPoint[i];
            if ( delta[i] > delta[index] )
                index = i;
            delta[i] = (delta[i] == 0.00) ? 1.00 : delta[i];
        }

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            alpha[i] = delta[i] / delta[index];
            mult_delta *= alpha[i];
        }

        mN[index] = static_cast<SizeType>( pow(static_cast<CoordinateType>(mObjectsSize/mult_delta), 1.00/Dimension) +1 );

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            if(i!=index)
            {
                mN[i] = static_cast<SizeType>(alpha[i] * mN[index]);
                mN[i] = ( mN[i] == 0 ) ? 1 : mN[i];
            }
        }

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            mCellSize[i] = delta[i] / mN[i];
            mInvCellSize[i] = 1.00 / mCellSize[i];
        }
    }

    void CalculateCellSize(const CoordinateType& CellSize)
    {
        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            mCellSize[i] = CellSize;
            mInvCellSize[i] = 1.00 / mCellSize[i];
            mN[i] = static_cast<SizeType>( (mMaxPoint[i]-mMinPoint[i]) / mCellSize[i]) + 1;
        }
    }

    void CalculateCellSize( const SizeType& NumPoints )
    {

        CoordinateType delta[Dimension];
        CoordinateType alpha[Dimension];
        CoordinateType mult_delta = 1.00;
        SizeType index = 0;

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            delta[i] = mMaxPoint[i] - mMinPoint[i];
            if ( delta[i] > delta[index] )
                index = i;
            delta[i] = (delta[i] == 0.00) ? 1.00 : delta[i];
        }

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            alpha[i] = delta[i] / delta[index];
            mult_delta *= alpha[i];
        }

        mN[index] = static_cast<SizeType>( pow(static_cast<CoordinateType>(NumPoints/mult_delta), 1.00/Dimension) +1 );

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            if(i!=index)
            {
                mN[i] = static_cast<SizeType>(alpha[i] * mN[index]);
                mN[i] = ( mN[i] == 0 ) ? 1 : mN[i];
            }
        }

        for(SizeType i = 0 ; i < Dimension ; i++)
        {
            mCellSize[i] = delta[i] / mN[i];
            mInvCellSize[i] = 1.00 / mCellSize[i];
        }
    }

//************************************************************************
//************************************************************************

    void GenerateBins()
    {
        PointType Low, High, Center;
        SearchStructureType Box;
        /// Fill container with objects
        for(IteratorType i_object = mObjectsBegin ; i_object != mObjectsEnd ; i_object++)
        {
            TConfigure::CalculateBoundingBox(*i_object, Low, High);
            Box.Set( CalculateCell(Low), CalculateCell(High), mN );
            FillObject(Box, *i_object);
        }
    }


//************************************************************************
//************************************************************************

    void GenerateCommunicationGraph()
    {
        double * MpiMinPoints = new double[mpi_size * Dimension];
        double * MpiMaxPoints = new double[mpi_size * Dimension];

        double * MyMinPoint = new double[Dimension];
        double * MyMaxPoint = new double[Dimension];

        mMinBoundingBox = vector<PointType>(mpi_size);
        mMaxBoundingBox = vector<PointType>(mpi_size);

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
            if(mpi_rank != i)
            {
                mpi_connectivity[i] = 0;

                for(size_t j = 0; j < Dimension; j++)
                {
                    mpi_MinPoints[i][j] = MpiMinPoints[i * Dimension + j];
                    mpi_MaxPoints[i][j] = MpiMaxPoints[i * Dimension + j];
                }

                mMinBoundingBox[i].X() = mpi_MinPoints[i][0];
                mMaxBoundingBox[i].X() = mpi_MaxPoints[i][0];

                mMinBoundingBox[i].Y() = mpi_MinPoints[i][1];
                mMaxBoundingBox[i].Y() = mpi_MaxPoints[i][1];

                mMinBoundingBox[i].Z() = mpi_MinPoints[i][2];
                mMaxBoundingBox[i].Z() = mpi_MaxPoints[i][2];
            }
        }

        delete [] MpiMinPoints;
        delete [] MpiMaxPoints;

        delete [] MyMinPoint;
        delete [] MyMaxPoint;
    }


//************************************************************************
//************************************************************************

// **** THREAD SAFE
// Dimension = 1
    void SearchInBoxLocal(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;
        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0])
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                mCells[I].SearchObjects(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchInBoxLocal(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
            {
                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    mCells[I].SearchObjects(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
            }
        }
    }

    // Dimension = 3
    void SearchInBoxLocal(PointerType& ThisObject, ResultIteratorType& Result,
                          SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {

        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2] += mCellSize[2], MaxCell[2] += mCellSize[2] )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
                {
                    if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    {
                        mCells[I].SearchObjects(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
                    }
                }
            }
        }
    }

//************************************************************************
//************************************************************************

// **** THREAD SAFE
// Dimension = 1
    void SearchObjectLocalInner(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                                SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;
        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0])
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                mCells[I].SearchObjectsInner(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchObjectLocalInner(PointerType& ThisObject, ResultIteratorType& Result, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                                SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
            {
                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    mCells[I].SearchObjectsInner(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
            }
        }
    }

    // Dimension = 3
    void SearchObjectLocalInner(PointerType& ThisObject, ResultIteratorType& Result,
                                SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                                SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {

        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2] += mCellSize[2], MaxCell[2] += mCellSize[2] )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
                {
                    if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    {
                        mCells[I].SearchObjectsInner(ThisObject, Result, NumberOfResults, MaxNumberOfResults);
                    }
                }
            }
        }
    }

//************************************************************************
//************************************************************************


    // **** THREAD SAFE

    // Dimension = 1
    void SearchInBoxLocal(PointerType& ThisObject, ResultContainerType& Result,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {

        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;
        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
        {
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                mCells[I].SearchObjects(ThisObject, Result);
        }
    }

    // Dimension = 2
    void SearchInBoxLocal(PointerType& ThisObject, ResultContainerType& Result,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];

        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
            {
                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    mCells[I].SearchObjects(ThisObject, Result);
            }
        }
    }

    // Dimension = 3
    void SearchInBoxLocal(PointerType& ThisObject, ResultContainerType& Result,
                          SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];

        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                {
                    if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    {
                        mCells[I].SearchObjects(ThisObject, Result);
                    }
                }
            }
        }
    }



//************************************************************************
//************************************************************************


    // **** THREAD SAFE

    // Dimension = 1
    void SearchObjectLocalInner(PointerType& ThisObject, ResultContainerType& Result,
                                SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;
        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block )
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                mCells[I].SearchObjectsInner(ThisObject, Result);

    }

    // Dimension = 2
    void SearchObjectLocalInner(PointerType& ThisObject, ResultContainerType& Result,
                                SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];

        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
            {
                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    mCells[I].SearchObjectsInner(ThisObject, Result);
            }
        }
    }

    // Dimension = 3
    void SearchObjectLocalInner(PointerType& ThisObject, ResultContainerType& Result,
                                SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];

        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block )
        {
            MinCell[2] = MinBox[2];
            MaxCell[2] = MaxBox[2];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block )
                {
                    if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell))
                    {
                        mCells[I].SearchObjectsInner(ThisObject, Result);
                    }
                }
            }
        }
    }


//************************************************************************
//************************************************************************


    // **** THREAD SAFE

    // Dimension = 1
    void SearchInRadius(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, DistanceIteratorType ResultDistances, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];

        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0])
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                mCells[I].SearchObjectsInRaius(ThisObject, Radius, Result, NumberOfResults, MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchInRadius(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, DistanceIteratorType ResultDistances, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
            {
                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                    mCells[I].SearchObjectsInRaius(ThisObject, Radius, Result, NumberOfResults, MaxNumberOfResults);
            }
        }
    }

    // Dimension = 3
    void SearchInRadius(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, DistanceIteratorType ResultDistances, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {

        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2] += mCellSize[2], MaxCell[2] += mCellSize[2] )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
                {
                    if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                    {
                        mCells[I].SearchObjectsInRadius(ThisObject, Radius, Result, ResultDistances, NumberOfResults, MaxNumberOfResults);
                    }
                }
            }
        }
    }

//************************************************************************
//************************************************************************


    // **** THREAD SAFE

    // Dimension = 1
    void SearchInRadiusExclusive(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, DistanceIteratorType ResultDistances, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];

        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0])
            if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                mCells[I].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, NumberOfResults, MaxNumberOfResults);
    }

    // Dimension = 2
    void SearchInRadiusExclusive(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, DistanceIteratorType ResultDistances, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box )
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
            {
                if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                    mCells[I].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, NumberOfResults, MaxNumberOfResults);
            }
        }
    }

    // Dimension = 3
    void SearchInRadiusExclusive(PointerType& ThisObject, CoordinateType const& Radius, ResultIteratorType& Result, DistanceIteratorType ResultDistances, SizeType& NumberOfResults, const SizeType& MaxNumberOfResults,
                        SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {

        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2] += mCellSize[2], MaxCell[2] += mCellSize[2] )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
                {
                    if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                    {
                        mCells[I].SearchObjectsInRadiusExclusive(ThisObject, Radius, Result, ResultDistances, NumberOfResults, MaxNumberOfResults);
                    }
                }
            }
        }
    }

    //TEMP!!!!
    int EmptyCellsArround(PointerType& ThisObject, CoordinateType const& Radius,SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box )
    {
        int found = 0;

        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2] += mCellSize[2], MaxCell[2] += mCellSize[2] )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1] += mCellSize[1], MaxCell[1] += mCellSize[1] )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0] += mCellSize[0], MaxCell[0] += mCellSize[0] )
                {
                    if(TConfigure::IntersectionBox(ThisObject, MinCell, MaxCell, Radius))
                    {
                        found += (mCells[I].Size() == 0) ? 1 : 0;
                    }
                }
            }
        }

        return found;
    }

//************************************************************************
//************************************************************************

    IndexArray  CalculateCell( const PointType& ThisPoint )
    {
        IndexArray IndexCell;
        for(SizeType i = 0 ; i < Dimension ; i++)
            IndexCell[i] = CalculatePosition(ThisPoint[i],i);
        return IndexCell;
    }

    IndexType CalculateIndex( const PointType& ThisPoint )
    {
        IndexType Index = 0;
        for(SizeType iDim = Dimension-1 ; iDim > 0 ; iDim--)
        {
            Index += CalculatePosition(ThisPoint[iDim],iDim);
            Index *= mN[iDim-1];
        }
        Index += CalculatePosition(ThisPoint[0],0);
        return Index;
    }

//************************************************************************
//************************************************************************

    // Dimension = 1
    void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, const PointerType& i_object)
    {
        PointType  MinCell, MaxCell;

        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
        {
            if(TConfigure::IntersectionBox(i_object, MinCell, MaxCell))
                mCells[I].Add(i_object);
        }
    }


//************************************************************************
//************************************************************************

    // Dimension = 2
    void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, const PointerType& i_object)
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
            {
                if(TConfigure::IntersectionBox(i_object,MinCell,MaxCell))
                    mCells[I].Add(i_object);
            }
        }
    }


//************************************************************************
//************************************************************************

    // Dimension = 3
    void FillObject( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, const PointerType& i_object)
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2]+=mCellSize[2], MaxCell[2]+=mCellSize[2] )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
                {
                    if(TConfigure::IntersectionBox(i_object,MinCell,MaxCell))
                        mCells[I].Add(i_object);
                }
            }
        }
    }

//     void FillDomain( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, const PointerType& i_object)
//     {
//         PointType  MinCell, MaxCell;
//         PointType  MinBox, MaxBox;
//
//         for(SizeType i = 0; i < 3; i++)
//         {
//             MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
//             MaxBox[i] = MinBox[i] + mCellSize[i];
//         }
//
//         MinCell[2] = MinBox[2];
//         MaxCell[2] = MaxBox[2];
//         for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2]+=mCellSize[2], MaxCell[2]+=mCellSize[2] )
//         {
//             MinCell[1] = MinBox[1];
//             MaxCell[1] = MaxBox[1];
//             for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
//             {
//                 MinCell[0] = MinBox[0];
//                 MaxCell[0] = MaxBox[0];
//                 for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
//                 {
//                     if(TConfigure::IntersectionBox(i_object,MinCell,MaxCell))
//                         mCellsDom[I]=1;
//                 }
//             }
//         }
//     }

//************************************************************************
//************************************************************************

    // Dimension = 1
    void RemoveObjectLocal( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,1>& Box, const PointerType& i_object)
    {
        PointType  MinCell, MaxCell;

        MinCell[0] = static_cast<CoordinateType>(Box.Axis[0].Min) * mCellSize[0] + mMinPoint[0];  //
        MaxCell[0] = MinCell[0] + mCellSize[0];
        for(IndexType I = Box.Axis[0].Begin() ; I <= Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
        {
            if(TConfigure::IntersectionBox(i_object, MinCell, MaxCell))
                mCells[I].Remove(i_object);
        }
    }


//************************************************************************
//************************************************************************

    // Dimension = 2
    void RemoveObjectLocal( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,2>& Box, const PointerType& i_object)
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 2; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[1] = MinBox[1];
        MaxCell[1] = MaxBox[1];
        for(IndexType II = Box.Axis[1].Begin() ; II <= Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
        {
            MinCell[0] = MinBox[0];
            MaxCell[0] = MaxBox[0];
            for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
            {
                if(TConfigure::IntersectionBox(i_object,MinCell,MaxCell))
                    mCells[I].Remove(i_object);
            }
        }
    }


//************************************************************************
//************************************************************************

    // Dimension = 3
    void RemoveObjectLocal( SearchStructure<IndexType,SizeType,CoordinateType,IteratorType,IteratorIteratorType,3>& Box, const PointerType& i_object)
    {
        PointType  MinCell, MaxCell;
        PointType  MinBox, MaxBox;

        for(SizeType i = 0; i < 3; i++)
        {
            MinBox[i] = static_cast<CoordinateType>(Box.Axis[i].Min) * mCellSize[i] + mMinPoint[i];  //
            MaxBox[i] = MinBox[i] + mCellSize[i];
        }

        MinCell[2] = MinBox[2];
        MaxCell[2] = MaxBox[2];
        for(IndexType III = Box.Axis[2].Begin() ; III <= Box.Axis[2].End() ; III += Box.Axis[2].Block, MinCell[2]+=mCellSize[2], MaxCell[2]+=mCellSize[2] )
        {
            MinCell[1] = MinBox[1];
            MaxCell[1] = MaxBox[1];
            for(IndexType II = III + Box.Axis[1].Begin() ; II <= III + Box.Axis[1].End() ; II += Box.Axis[1].Block, MinCell[1]+=mCellSize[1], MaxCell[1]+=mCellSize[1] )
            {
                MinCell[0] = MinBox[0];
                MaxCell[0] = MaxBox[0];
                for(IndexType I = II + Box.Axis[0].Begin() ; I <= II + Box.Axis[0].End() ; I += Box.Axis[0].Block, MinCell[0]+=mCellSize[0], MaxCell[0]+=mCellSize[0] )
                {
                    if(TConfigure::IntersectionBox(i_object,MinCell,MaxCell))
                        mCells[I].Remove(i_object);
                }
            }
        }
    }

//************************************************************************
//************************************************************************

    IndexType CalculatePosition( CoordinateType const& ThisCoord, const SizeType& ThisDimension )
    {
        CoordinateType d_index = (ThisCoord - mMinPoint[ThisDimension]) * mInvCellSize[ThisDimension];
        IndexType index = static_cast<IndexType>( (d_index < 0.00) ? 0.00 : d_index );
        return  (index > mN[ThisDimension]-1) ? mN[ThisDimension]-1 : index;

    }


//************************************************************************
//************************************************************************

    void AllocateContainer()
    {
        SizeType Size = mN[0];
        for(SizeType i = 1 ; i < Dimension ; i++)
            Size *= mN[i];
        mCells.resize(Size);
//         mCellsDom.resize(Size);
    }


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    PointType    mMinPoint;
    PointType    mMaxPoint;

    IteratorType mObjectsBegin;
    IteratorType mObjectsEnd;
    SizeType     mObjectsSize;

    CoordinateArray  mCellSize;
    CoordinateArray  mInvCellSize;
    SizeArray        mN;

    CellContainerType   mCells;  ///The bin
//     DomainContainerType mCellsDom;  ///The domain bin

    ///@}
    ///@name MPI Variables
    ///@{

    int mpi_rank;
    int mpi_size;

    int BinsExtension;

    vector<int> mpi_connectivity;
    vector<vector<double> > mpi_MinPoints;
    vector<vector<double> > mpi_MaxPoints;

    vector<PointType> mMinBoundingBox;
    vector<PointType> mMaxBoundingBox;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    inline void CreatePartition(SizeType number_of_threads, const SizeType number_of_rows, std::vector<SizeType>& partitions)
    {
        partitions.resize(number_of_threads+1);
        SizeType partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for(SizeType i = 1; i<number_of_threads; i++)
            partitions[i] = partitions[i-1] + partition_size ;
    }


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}


public:
    /// Assignment operator.
    BinsObjectDynamicMpi<TConfigure> & operator=(const BinsObjectDynamicMpi<TConfigure> & rOther)
    {
        mMinPoint            = rOther.mMinPoint;
        mMaxPoint            = rOther.mMaxPoint;
        mObjectsBegin        = rOther.mObjectsBegin;
        mObjectsEnd          = rOther.mObjectsEnd;
        mObjectsSize         = rOther.mObjectsSize;
        mCellSize            = rOther.mCellSize;
        mInvCellSize         = rOther.mInvCellSize;
        mN                   = rOther.mN;
        mCells               = rOther.mCells;
//         mCellsDom            = rOther.mCellsDom;
        return *this;
    }

    /// Copy constructor.
    BinsObjectDynamicMpi(const BinsObjectDynamicMpi& rOther)
    {
        *this =  rOther;
    }





}; // Class BinsObjectDynamic

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TConfigure>
inline std::istream& operator >> (std::istream& rIStream,
                                  BinsObjectDynamicMpi<TConfigure>& rThis)
{
    return rIStream;
}


/// output stream function
template<class TConfigure>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BinsObjectDynamicMpi<TConfigure> & rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined
