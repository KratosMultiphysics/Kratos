//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//



#if !defined(KRATOS_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_COMMUNICATOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


// External includes

// Project includes
#include "includes/data_communicator.h"
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/element.h"
#include "includes/mesh.h"
#include "includes/parallel_environment.h"

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
class Communicator
{
public:
    ///@name  Enum's
    ///@{


    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Communicator
    KRATOS_CLASS_POINTER_DEFINITION(Communicator);

    typedef unsigned int IndexType;

    typedef unsigned int SizeType;

    typedef Node < 3 > NodeType;

    typedef Properties PropertiesType;

    typedef Element ElementType;

    typedef Condition ConditionType;

    typedef DenseVector<int> NeighbourIndicesContainerType;

    typedef Mesh<NodeType, PropertiesType, ElementType, ConditionType> MeshType;

    typedef PointerVector<MeshType> MeshesContainerType;

    /// Nodes container. Which is a vector set of nodes with their Id's as key.
    typedef MeshType::NodesContainerType NodesContainerType;

    /** Iterator over the nodes. This iterator is an indirect
        iterator over Node::Pointer which turn back a reference to
        node by * operator and not a pointer for more convenient
        usage. */
    typedef MeshType::NodeIterator NodeIterator;

    /** Const iterator over the nodes. This iterator is an indirect
        iterator over Node::Pointer which turn back a reference to
        node by * operator and not a pointer for more convenient
        usage. */
    typedef MeshType::NodeConstantIterator NodeConstantIterator;

    /** Iterator over the properties. This iterator is an indirect
        iterator over Properties::Pointer which turn back a reference to
        properties by * operator and not a pointer for more convenient
        usage. */

    /// Properties container. Which is a vector set of Properties with their Id's as key.
    typedef MeshType::PropertiesContainerType PropertiesContainerType;

    /** Iterator over the Properties. This iterator is an indirect
        iterator over Node::Pointer which turn back a reference to
        node by * operator and not a pointer for more convenient
        usage. */
    typedef MeshType::PropertiesIterator PropertiesIterator;

    /** Const iterator over the Properties. This iterator is an indirect
        iterator over Properties::Pointer which turn back a reference to
        Properties by * operator and not a pointer for more convenient
        usage. */
    typedef MeshType::PropertiesConstantIterator PropertiesConstantIterator;

    /** Iterator over the properties. This iterator is an indirect
        iterator over Properties::Pointer which turn back a reference to
        properties by * operator and not a pointer for more convenient
        usage. */

    /// Element container. A vector set of Elements with their Id's as key.
    typedef MeshType::ElementsContainerType ElementsContainerType;

    /** Iterator over the Elements. This iterator is an indirect
        iterator over Elements::Pointer which turn back a reference to
        Element by * operator and not a pointer for more convenient
        usage. */
    typedef MeshType::ElementIterator ElementIterator;

    /** Const iterator over the Elements. This iterator is an indirect
        iterator over Elements::Pointer which turn back a reference to
        Element by * operator and not a pointer for more convenient
        usage. */
    typedef MeshType::ElementConstantIterator ElementConstantIterator;

    /// Condintions container. A vector set of Conditions with their Id's as key.
    typedef MeshType::ConditionsContainerType ConditionsContainerType;

    /** Iterator over the Conditions. This iterator is an indirect
       iterator over Conditions::Pointer which turn back a reference to
       Condition by * operator and not a pointer for more convenient
       usage. */
    typedef MeshType::ConditionIterator ConditionIterator;

    /** Const iterator over the Conditions. This iterator is an indirect
        iterator over Conditions::Pointer which turn back a reference to
        Condition by * operator and not a pointer for more convenient
        usage. */
    typedef MeshType::ConditionConstantIterator ConditionConstantIterator;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Communicator() : mNumberOfColors(1)
        , mpLocalMesh(MeshType::Pointer(new MeshType))
        , mpGhostMesh(MeshType::Pointer(new MeshType))
        , mpInterfaceMesh(MeshType::Pointer(new MeshType))
        , mrDataCommunicator(ParallelEnvironment::GetDataCommunicator("Serial"))

    {
        MeshType mesh;
        mLocalMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
        mGhostMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
        mInterfaceMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
    }

    /// Constructor using a custom DataCommunicator.
    /** This constructor is intended for use from derived classes,
     *  since the base Communicator class will often not use the communicator at all.
     *  @param rDataCommunicator Reference to a DataCommunicator.
     */
    Communicator(const DataCommunicator& rDataCommunicator)
        : mNumberOfColors(1)
        , mpLocalMesh(MeshType::Pointer(new MeshType))
        , mpGhostMesh(MeshType::Pointer(new MeshType))
        , mpInterfaceMesh(MeshType::Pointer(new MeshType))
        , mrDataCommunicator(rDataCommunicator)
    {
        MeshType mesh;
        mLocalMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
        mGhostMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
        mInterfaceMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
    }

    /// Copy constructor.

    Communicator(Communicator const& rOther)
        : mNumberOfColors(rOther.mNumberOfColors)
        , mNeighbourIndices(rOther.mNeighbourIndices)
        , mpLocalMesh(MeshType::Pointer(rOther.mpLocalMesh))
        , mpGhostMesh(MeshType::Pointer(rOther.mpGhostMesh))
        , mpInterfaceMesh(MeshType::Pointer(rOther.mpInterfaceMesh))
        , mLocalMeshes(rOther.mLocalMeshes)
        , mGhostMeshes(rOther.mGhostMeshes)
        , mInterfaceMeshes(rOther.mInterfaceMeshes)
        , mrDataCommunicator(rOther.mrDataCommunicator)
    {
    }

    virtual Communicator::Pointer Create(const DataCommunicator& rDataCommunicator) const
    {
        KRATOS_TRY

        return Kratos::make_shared<Communicator>(rDataCommunicator);

        KRATOS_CATCH("");
    }

    virtual Communicator::Pointer Create() const
    {
        return Create(ParallelEnvironment::GetDataCommunicator("Serial"));
    }

    /// Destructor.
    virtual ~Communicator() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Communicator & operator=(Communicator const& rOther) = delete;

    ///@}
    ///@name Access
    ///@{

    virtual bool IsDistributed() const
    {
        return false;
    }

    virtual int MyPID() const
    {
        return mrDataCommunicator.Rank();
    }

    virtual int TotalProcesses() const
    {
        return mrDataCommunicator.Size();
    }

    SizeType GetNumberOfColors() const
    {
        return mNumberOfColors;
    }

    void SetNumberOfColors(SizeType NewNumberOfColors)
    {
        if (mNumberOfColors == NewNumberOfColors)
            return;

        mNumberOfColors = NewNumberOfColors;
        MeshType mesh;

        mLocalMeshes.clear();
        mGhostMeshes.clear();
        mInterfaceMeshes.clear();

        for (IndexType i = 0; i < mNumberOfColors; i++)
        {
            mLocalMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
            mGhostMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
            mInterfaceMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
        }
    }

    NeighbourIndicesContainerType& NeighbourIndices()
    {
        return mNeighbourIndices;
    }

    NeighbourIndicesContainerType const& NeighbourIndices() const
    {
        return mNeighbourIndices;
    }

    // Set the local mesh pointer to the given mesh
    void SetLocalMesh(MeshType::Pointer pGivenMesh)
    {
        mpLocalMesh = pGivenMesh;
    }

    // Returns pointer to the mesh storing all local entites

    MeshType::Pointer pLocalMesh()
    {
        return mpLocalMesh;
    }

    // Returns pointer to the mesh storing all ghost entites

    MeshType::Pointer pGhostMesh()
    {
        return mpGhostMesh;
    }

    // Returns pointer to the mesh storing all interface entites

    MeshType::Pointer pInterfaceMesh()
    {
        return mpInterfaceMesh;
    }

    // Returns a constant pointer to the mesh storing all local entites

    const MeshType::Pointer pLocalMesh() const
    {
        return mpLocalMesh;
    }

    // Returns a constant pointer to the mesh storing all ghost entites

    const MeshType::Pointer pGhostMesh() const
    {
        return mpGhostMesh;
    }

    // Returns a constant pointer to the mesh storing all interface entites

    const MeshType::Pointer pInterfaceMesh() const
    {
        return mpInterfaceMesh;
    }

    MeshType::Pointer pLocalMesh(IndexType ThisIndex)
    {
        return mLocalMeshes(ThisIndex);
    }

    MeshType::Pointer pGhostMesh(IndexType ThisIndex)
    {
        return mGhostMeshes(ThisIndex);
    }

    MeshType::Pointer pInterfaceMesh(IndexType ThisIndex)
    {
        return mInterfaceMeshes(ThisIndex);
    }

    const MeshType::Pointer pLocalMesh(IndexType ThisIndex) const
    {
        return mLocalMeshes(ThisIndex);
    }

    const MeshType::Pointer pGhostMesh(IndexType ThisIndex) const
    {
        return mGhostMeshes(ThisIndex);
    }

    const MeshType::Pointer pInterfaceMesh(IndexType ThisIndex) const
    {
        return mInterfaceMeshes(ThisIndex);
    }

    // Returns the reference to the mesh storing all local entites

    MeshType& LocalMesh()
    {
        return *mpLocalMesh;
    }

    // Returns the reference to the mesh storing all ghost entites

    MeshType& GhostMesh()
    {
        return *mpGhostMesh;
    }

    // Returns the reference to the mesh storing all interface entites

    MeshType& InterfaceMesh()
    {
        return *mpInterfaceMesh;
    }

    // Returns a constant reference to the mesh storing all local entites

    MeshType const& LocalMesh() const
    {
        return *mpLocalMesh;
    }

    // Returns a constant reference to the mesh storing all ghost entites

    MeshType const& GhostMesh() const
    {
        return *mpGhostMesh;
    }

    // Returns a constant reference to the mesh storing all interface entites

    MeshType const& InterfaceMesh() const
    {
        return *mpInterfaceMesh;
    }

    MeshType& LocalMesh(IndexType ThisIndex)
    {
        return mLocalMeshes[ThisIndex];
    }

    MeshType& GhostMesh(IndexType ThisIndex)
    {
        return mGhostMeshes[ThisIndex];
    }

    MeshType& InterfaceMesh(IndexType ThisIndex)
    {
        return mInterfaceMeshes[ThisIndex];
    }

    MeshType const& LocalMesh(IndexType ThisIndex) const
    {
        return mLocalMeshes[ThisIndex];
    }

    MeshType const& GhostMesh(IndexType ThisIndex) const
    {
        return mGhostMeshes[ThisIndex];
    }

    MeshType const& InterfaceMesh(IndexType ThisIndex) const
    {
        return mInterfaceMeshes[ThisIndex];
    }

    MeshesContainerType& LocalMeshes()
    {
        return mLocalMeshes;
    }

    MeshesContainerType& GhostMeshes()
    {
        return mGhostMeshes;
    }

    MeshesContainerType& InterfaceMeshes()
    {
        return mInterfaceMeshes;
    }

    MeshesContainerType const& LocalMeshes() const
    {
        return mLocalMeshes;
    }

    MeshesContainerType const& GhostMeshes() const
    {
        return mGhostMeshes;
    }

    MeshesContainerType const& InterfaceMeshes() const
    {
        return mInterfaceMeshes;
    }

    ///@}
    ///@name Operations
    ///@{

    KRATOS_DEPRECATED_MESSAGE("This function is deprecated, please retrieve the DataCommunicator with GetDataCommunicator and use it directly.")
    void Barrier() const
    {
        mrDataCommunicator.Barrier();
    }

    KRATOS_DEPRECATED_MESSAGE("This function is deprecated, please retrieve the DataCommunicator with GetDataCommunicator and use it directly.")
    bool SumAll(int& rValue) const
    {
        rValue = mrDataCommunicator.SumAll(rValue);
        return true;
    }

    KRATOS_DEPRECATED_MESSAGE("This function is deprecated, please retrieve the DataCommunicator with GetDataCommunicator and use it directly.")
    bool SumAll(double& rValue) const
    {
        rValue = mrDataCommunicator.SumAll(rValue);
        return true;
    }

    KRATOS_DEPRECATED_MESSAGE("This function is deprecated, please retrieve the DataCommunicator with GetDataCommunicator and use it directly.")
    bool SumAll(array_1d<double, 3>& rValue) const
    {
        rValue = mrDataCommunicator.SumAll(rValue);
        return true;
    }

    KRATOS_DEPRECATED_MESSAGE("This function is deprecated, please retrieve the DataCommunicator with GetDataCommunicator and use it directly.")
    bool MinAll(int& rValue) const
    {
        rValue = mrDataCommunicator.MinAll(rValue);
        return true;
    }

    KRATOS_DEPRECATED_MESSAGE("This function is deprecated, please retrieve the DataCommunicator with GetDataCommunicator and use it directly.")
    bool MinAll(double& rValue) const
    {
        rValue = mrDataCommunicator.MinAll(rValue);
        return true;
    }

    KRATOS_DEPRECATED_MESSAGE("This function is deprecated, please retrieve the DataCommunicator with GetDataCommunicator and use it directly.")
    bool MaxAll(int& rValue) const
    {
        rValue = mrDataCommunicator.MaxAll(rValue);
        return true;
    }

    KRATOS_DEPRECATED_MESSAGE("This function is deprecated, please retrieve the DataCommunicator with GetDataCommunicator and use it directly.")
    bool MaxAll(double& rValue) const
    {
        rValue = mrDataCommunicator.MaxAll(rValue);
        return true;
    }

    KRATOS_DEPRECATED_MESSAGE("This function is deprecated, please retrieve the DataCommunicator with GetDataCommunicator and use it directly.")
    bool ScanSum(const double& send_partial, double& receive_accumulated) const
    {
        receive_accumulated = mrDataCommunicator.ScanSum(send_partial);
        return true;
    }

    KRATOS_DEPRECATED_MESSAGE("This function is deprecated, please retrieve the DataCommunicator with GetDataCommunicator and use it directly.")
    bool ScanSum(const int& send_partial, int& receive_accumulated) const
    {
        receive_accumulated = mrDataCommunicator.ScanSum(send_partial);
        return true;
    }

    virtual bool SynchronizeNodalSolutionStepsData()
    {
        // #if defined(KRATOS_USING_MPI )
        // 	std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        // #endif
        return true;

    }

    virtual bool SynchronizeDofs()
    {
        // #if defined(KRATOS_USING_MPI )
        // 	std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        // #endif
        return true;

    }

    virtual bool SynchronizeVariable(Variable<int> const& rThisVariable)
    {
        // #if defined(KRATOS_USING_MPI )
        //  std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        // #endif
        return true;
    }

    virtual bool SynchronizeVariable(Variable<double> const& rThisVariable)
    {
        // #if defined(KRATOS_USING_MPI )
        //  std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        // #endif
        return true;
    }

    virtual bool SynchronizeVariable(Variable<bool> const& rThisVariable)
    {
        return true;
    }

    virtual bool SynchronizeVariable(Variable<array_1d<double, 3 > > const& rThisVariable)
    {
        // #if defined(KRATOS_USING_MPI )
        //  std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        // #endif
        return true;
    }

    virtual bool SynchronizeVariable(Variable<array_1d<double, 4 > > const& rThisVariable)
    {
        return true;
    }

    virtual bool SynchronizeVariable(Variable<array_1d<double, 6 > > const& rThisVariable)
    {
        return true;
    }

    virtual bool SynchronizeVariable(Variable<array_1d<double, 9 > > const& rThisVariable)
    {
        return true;
    }

    virtual bool SynchronizeVariable(Variable<Vector> const& rThisVariable)
    {
        // #if defined(KRATOS_USING_MPI )
        //  std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        // #endif
        return true;
    }

    virtual bool SynchronizeVariable(Variable<Matrix> const& rThisVariable)
    {
        // #if defined(KRATOS_USING_MPI )
        //  std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        // #endif
        return true;
    }

    virtual bool SynchronizeNonHistoricalVariable(Variable<int> const& rThisVariable)
    {
        return true;
    }

    virtual bool SynchronizeNonHistoricalVariable(Variable<double> const& rThisVariable)
    {
        return true;
    }

    virtual bool SynchronizeNonHistoricalVariable(Variable<bool> const& rThisVariable)
    {
        return true;
    }

    virtual bool SynchronizeNonHistoricalVariable(Variable<array_1d<double, 3 > > const& rThisVariable)
    {
        return true;
    }

    virtual bool SynchronizeNonHistoricalVariable(Variable<array_1d<double, 4 > > const& rThisVariable)
    {
        return true;
    }

    virtual bool SynchronizeNonHistoricalVariable(Variable<array_1d<double, 6 > > const& rThisVariable)
    {
        return true;
    }

    virtual bool SynchronizeNonHistoricalVariable(Variable<array_1d<double, 9 > > const& rThisVariable)
    {
        return true;
    }

    virtual bool SynchronizeNonHistoricalVariable(Variable<Vector> const& rThisVariable)
    {
        return true;
    }

    virtual bool SynchronizeNonHistoricalVariable(Variable<Matrix> const& rThisVariable)
    {
        return true;
    }

    /// Synchronize variable in nodal solution step data to the minimum value across all processes.
    /** @param ThisVariable The variable to be synchronized.
     */
    virtual bool SynchronizeCurrentDataToMin(Variable<double> const& ThisVariable)
    {
        return true;
    }

    /// Synchronize variable in nodal data to the minimum value across all processes.
    /** @param ThisVariable The variable to be synchronized.
     */
    virtual bool SynchronizeNonHistoricalDataToMin(Variable<double> const& ThisVariable)
    {
        return true;
    }

    virtual bool SynchronizeElementalFlags()
    {
        return true;
    }

    virtual bool AssembleCurrentData(Variable<int> const& ThisVariable)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;

    }

    virtual bool AssembleCurrentData(Variable<double> const& ThisVariable)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;

    }

    virtual bool AssembleCurrentData(Variable<array_1d<double, 3 > > const& ThisVariable)
    {
        // #if defined(KRATOS_USING_MPI )
        // 	std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        // #endif
        return true;

    }

    virtual bool AssembleCurrentData(Variable<Vector> const& ThisVariable)
    {
        // #if defined(KRATOS_USING_MPI )
        // 	std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        // #endif
        return true;

    }

    virtual bool AssembleCurrentData(Variable<Matrix> const& ThisVariable)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;

    }

    virtual bool AssembleNonHistoricalData(Variable<int> const& ThisVariable)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;

    }

    virtual bool AssembleNonHistoricalData(Variable<double> const& ThisVariable)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;

    }

    virtual bool AssembleNonHistoricalData(Variable<array_1d<double, 3 > > const& ThisVariable)
    {
        // #if defined(KRATOS_USING_MPI )
        // 	std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        // #endif
        return true;

    }


    virtual bool AssembleNonHistoricalData(Variable<DenseVector<array_1d<double,3> > > const& ThisVariable)
    {
        // #if defined(KRATOS_USING_MPI )
        //  std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        // #endif
        return true;

    }

    virtual bool AssembleNonHistoricalData(Variable<Vector> const& ThisVariable)
    {
        // #if defined(KRATOS_USING_MPI )
        // 	std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        // #endif
        return true;

    }

    virtual bool AssembleNonHistoricalData(Variable<Matrix> const& ThisVariable)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;

    }

    virtual bool SynchronizeElementalNonHistoricalVariable(Variable<int> const& ThisVariable)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;
    }

    virtual bool SynchronizeElementalNonHistoricalVariable(Variable<double> const& ThisVariable)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;
    }

    virtual bool SynchronizeElementalNonHistoricalVariable(Variable<array_1d<double, 3 > > const& ThisVariable)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;
    }

    virtual bool SynchronizeElementalNonHistoricalVariable(Variable<DenseVector<array_1d<double,3> > > const& ThisVariable)
    {
    /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;
    }

    virtual bool SynchronizeElementalNonHistoricalVariable(Variable<DenseVector<int> > const& ThisVariable)
    {
    /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;
    }

    virtual bool SynchronizeElementalNonHistoricalVariable(Variable<Vector> const& ThisVariable)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;
    }

    virtual bool SynchronizeElementalNonHistoricalVariable(Variable<Matrix> const& ThisVariable)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;
    }

    virtual bool TransferObjects(std::vector<NodesContainerType>& SendObjects, std::vector<NodesContainerType>& RecvObjects)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;
    }

    virtual bool TransferObjects(std::vector<ElementsContainerType>& SendObjects, std::vector<ElementsContainerType>& RecvObjects)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;
    }

    virtual bool TransferObjects(std::vector<ConditionsContainerType>& SendObjects, std::vector<ConditionsContainerType>& RecvObjects)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;
    }

    virtual bool TransferObjects(std::vector<NodesContainerType>& SendObjects, std::vector<NodesContainerType>& RecvObjects,Kratos::Serializer& particleSerializer)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;
    }

    virtual bool TransferObjects(std::vector<ElementsContainerType>& SendObjects, std::vector<ElementsContainerType>& RecvObjects,Kratos::Serializer& particleSerializer)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;
    }

    virtual bool TransferObjects(std::vector<ConditionsContainerType>& SendObjects, std::vector<ConditionsContainerType>& RecvObjects,Kratos::Serializer& particleSerializer)
    {
        /*#if defined(KRATOS_USING_MPI )
                std::cout << "WARNING: Using serial communicator with MPI defined. Use ModelPart::SetCommunicator to set its communicator to MPICommunicator" << std::endl;
        #endif*/
        return true;
    }

    virtual bool SynchronizeOrNodalFlags(const Flags& TheFlags)
    {
        return true;
    }

    virtual bool SynchronizeAndNodalFlags(const Flags& TheFlags)
    {
        return true;
    }

    virtual bool SynchronizeNodalFlags()
    {
        return true;
    }

    void Clear()
    {
        mNumberOfColors = 0;
        mNeighbourIndices.clear();
        mpLocalMesh->MeshType::Clear();
        mpGhostMesh->MeshType::Clear();
        mpInterfaceMesh->MeshType::Clear();
        mLocalMeshes.clear();
        mGhostMeshes.clear();
        mInterfaceMeshes.clear();
    }

    ///@}
    ///@name Access
    ///@{

    virtual const DataCommunicator& GetDataCommunicator() const
    {
        return mrDataCommunicator;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.

    virtual std::string Info() const
    {
        return "Communicator";
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const
    {
        for (IndexType i = 0; i < mLocalMeshes.size(); i++)
        {
            rOStream << "    Local Mesh " << i << " : " << std::endl;
            LocalMesh(i).PrintData(rOStream);
            rOStream << "    Ghost Mesh " << i << " : " << std::endl;
            GhostMesh(i).PrintData(rOStream);
            rOStream << "    Interface Mesh " << i << " : " << std::endl;
            InterfaceMesh(i).PrintData(rOStream);
        }
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{
    SizeType mNumberOfColors;

    NeighbourIndicesContainerType mNeighbourIndices;

    // To store all local entities
    MeshType::Pointer mpLocalMesh;

    // To store all ghost entities
    MeshType::Pointer mpGhostMesh;

    // To store all interface entities
    MeshType::Pointer mpInterfaceMesh;

    // To store interfaces local entities
    MeshesContainerType mLocalMeshes;

    // To store interfaces ghost entities
    MeshesContainerType mGhostMeshes;

    // To store interfaces ghost+local entities
    MeshesContainerType mInterfaceMeshes;

    // Interface to MPI communication
    const DataCommunicator& mrDataCommunicator;

    ///@}
    ///@name Protected Operators
    ///@{


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




    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


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

}; // Class Communicator

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  Communicator& rThis);

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const Communicator& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


} // namespace Kratos.

#endif // KRATOS_COMMUNICATOR_H_INCLUDED  defined


