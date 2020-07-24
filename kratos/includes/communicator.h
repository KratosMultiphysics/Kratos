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
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/element.h"
#include "includes/mesh.h"

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

// Forward declaration of DataCommunicator
class DataCommunicator;

/// The Commmunicator class manages communication for distributed ModelPart instances.
/** The base Communicator class only holds the required data (local and remote mesh interfaces)
 *  for communication. The actual communication is implemented in the derived MPICommunicator.
 */
class KRATOS_API(KRATOS_CORE) Communicator
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
    Communicator();

    /// Constructor using a custom DataCommunicator.
    /** This constructor is intended for use from derived classes,
     *  since the base Communicator class will often not use the communicator at all.
     *  @param rDataCommunicator Reference to a DataCommunicator.
     */
    Communicator(const DataCommunicator& rDataCommunicator);

    /// Copy constructor.
    Communicator(Communicator const& rOther);

    /// Destructor.
    virtual ~Communicator() = default;

    virtual Communicator::Pointer Create(const DataCommunicator& rDataCommunicator) const;

    virtual Communicator::Pointer Create() const;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Communicator & operator=(Communicator const& rOther) = delete;

    ///@}
    ///@name Access
    ///@{

    virtual bool IsDistributed() const;

    virtual int MyPID() const;

    virtual int TotalProcesses() const;

    SizeType GlobalNumberOfNodes() const;

    SizeType GlobalNumberOfElements() const;

    SizeType GlobalNumberOfConditions() const;

    SizeType GetNumberOfColors() const;

    void SetNumberOfColors(SizeType NewNumberOfColors);

    void AddColors(SizeType NumberOfAddedColors);

    NeighbourIndicesContainerType& NeighbourIndices();

    NeighbourIndicesContainerType const& NeighbourIndices() const;

    /// Set the local mesh pointer to the given mesh
    void SetLocalMesh(MeshType::Pointer pGivenMesh);

    /// Returns pointer to the mesh storing all local entites
    MeshType::Pointer pLocalMesh();

    /// Returns pointer to the mesh storing all ghost entites
    MeshType::Pointer pGhostMesh();

    /// Returns pointer to the mesh storing all interface entites
    MeshType::Pointer pInterfaceMesh();

    /// Returns a constant pointer to the mesh storing all local entites
    const MeshType::Pointer pLocalMesh() const;

    /// Returns a constant pointer to the mesh storing all ghost entites
    const MeshType::Pointer pGhostMesh() const;

    /// Returns a constant pointer to the mesh storing all interface entites
    const MeshType::Pointer pInterfaceMesh() const;

    MeshType::Pointer pLocalMesh(IndexType ThisIndex);

    MeshType::Pointer pGhostMesh(IndexType ThisIndex);

    MeshType::Pointer pInterfaceMesh(IndexType ThisIndex);

    const MeshType::Pointer pLocalMesh(IndexType ThisIndex) const;

    const MeshType::Pointer pGhostMesh(IndexType ThisIndex) const;

    const MeshType::Pointer pInterfaceMesh(IndexType ThisIndex) const;

    /// Returns the reference to the mesh storing all local entites
    MeshType& LocalMesh();

    /// Returns the reference to the mesh storing all ghost entites
    MeshType& GhostMesh();

    /// Returns the reference to the mesh storing all interface entites
    MeshType& InterfaceMesh();

    /// Returns a constant reference to the mesh storing all local entites
    MeshType const& LocalMesh() const;

    /// Returns a constant reference to the mesh storing all ghost entites
    MeshType const& GhostMesh() const;

    /// Returns a constant reference to the mesh storing all interface entites
    MeshType const& InterfaceMesh() const;

    MeshType& LocalMesh(IndexType ThisIndex);

    MeshType& GhostMesh(IndexType ThisIndex);

    MeshType& InterfaceMesh(IndexType ThisIndex);

    MeshType const& LocalMesh(IndexType ThisIndex) const;

    MeshType const& GhostMesh(IndexType ThisIndex) const;

    MeshType const& InterfaceMesh(IndexType ThisIndex) const;

    MeshesContainerType& LocalMeshes();

    MeshesContainerType& GhostMeshes();

    MeshesContainerType& InterfaceMeshes();

    MeshesContainerType const& LocalMeshes() const;

    MeshesContainerType const& GhostMeshes() const;

    MeshesContainerType const& InterfaceMeshes() const;

    virtual const DataCommunicator& GetDataCommunicator() const;

    ///@}
    ///@name Operations
    ///@{

    virtual bool SynchronizeNodalSolutionStepsData();

    virtual bool SynchronizeDofs();

    virtual bool SynchronizeVariable(Variable<int> const& rThisVariable);

    virtual bool SynchronizeVariable(Variable<double> const& rThisVariable);

    virtual bool SynchronizeVariable(Variable<bool> const& rThisVariable);

    virtual bool SynchronizeVariable(Variable<array_1d<double, 3 > > const& rThisVariable);

    virtual bool SynchronizeVariable(Variable<array_1d<double, 4 > > const& rThisVariable);

    virtual bool SynchronizeVariable(Variable<array_1d<double, 6 > > const& rThisVariable);

    virtual bool SynchronizeVariable(Variable<array_1d<double, 9 > > const& rThisVariable);

    virtual bool SynchronizeVariable(Variable<Vector> const& rThisVariable);

    virtual bool SynchronizeVariable(Variable<Matrix> const& rThisVariable);

    virtual bool SynchronizeNonHistoricalVariable(Variable<int> const& rThisVariable);

    virtual bool SynchronizeNonHistoricalVariable(Variable<double> const& rThisVariable);

    virtual bool SynchronizeNonHistoricalVariable(Variable<bool> const& rThisVariable);

    virtual bool SynchronizeNonHistoricalVariable(Variable<array_1d<double, 3 > > const& rThisVariable);

    virtual bool SynchronizeNonHistoricalVariable(Variable<array_1d<double, 4 > > const& rThisVariable);

    virtual bool SynchronizeNonHistoricalVariable(Variable<array_1d<double, 6 > > const& rThisVariable);

    virtual bool SynchronizeNonHistoricalVariable(Variable<array_1d<double, 9 > > const& rThisVariable);

    virtual bool SynchronizeNonHistoricalVariable(Variable<Vector> const& rThisVariable);

    virtual bool SynchronizeNonHistoricalVariable(Variable<Matrix> const& rThisVariable);

    /// Synchronize variable in nodal solution step data to the minimum value across all processes.
    /** @param ThisVariable The variable to be synchronized.
     */
    virtual bool SynchronizeCurrentDataToMin(Variable<double> const& ThisVariable);

    /// Synchronize variable in nodal data to the minimum value across all processes.
    /** @param ThisVariable The variable to be synchronized.
     */
    virtual bool SynchronizeNonHistoricalDataToMin(Variable<double> const& ThisVariable);

    virtual bool SynchronizeElementalFlags();

    virtual bool AssembleCurrentData(Variable<int> const& ThisVariable);

    virtual bool AssembleCurrentData(Variable<double> const& ThisVariable);

    virtual bool AssembleCurrentData(Variable<array_1d<double, 3 > > const& ThisVariable);

    virtual bool AssembleCurrentData(Variable<Vector> const& ThisVariable);

    virtual bool AssembleCurrentData(Variable<Matrix> const& ThisVariable);

    virtual bool AssembleNonHistoricalData(Variable<int> const& ThisVariable);

    virtual bool AssembleNonHistoricalData(Variable<double> const& ThisVariable);

    virtual bool AssembleNonHistoricalData(Variable<array_1d<double, 3 > > const& ThisVariable);

    virtual bool AssembleNonHistoricalData(Variable<DenseVector<array_1d<double,3> > > const& ThisVariable);

    virtual bool AssembleNonHistoricalData(Variable<Vector> const& ThisVariable);

    virtual bool AssembleNonHistoricalData(Variable<Matrix> const& ThisVariable);

    virtual bool SynchronizeElementalNonHistoricalVariable(Variable<int> const& ThisVariable);

    virtual bool SynchronizeElementalNonHistoricalVariable(Variable<double> const& ThisVariable);

    virtual bool SynchronizeElementalNonHistoricalVariable(Variable<array_1d<double, 3 > > const& ThisVariable);

    virtual bool SynchronizeElementalNonHistoricalVariable(Variable<DenseVector<array_1d<double,3> > > const& ThisVariable);

    virtual bool SynchronizeElementalNonHistoricalVariable(Variable<DenseVector<int> > const& ThisVariable);

    virtual bool SynchronizeElementalNonHistoricalVariable(Variable<Vector> const& ThisVariable);

    virtual bool SynchronizeElementalNonHistoricalVariable(Variable<Matrix> const& ThisVariable);

    virtual bool TransferObjects(std::vector<NodesContainerType>& SendObjects, std::vector<NodesContainerType>& RecvObjects);

    virtual bool TransferObjects(std::vector<ElementsContainerType>& SendObjects, std::vector<ElementsContainerType>& RecvObjects);

    virtual bool TransferObjects(std::vector<ConditionsContainerType>& SendObjects, std::vector<ConditionsContainerType>& RecvObjects);

    virtual bool TransferObjects(std::vector<NodesContainerType>& SendObjects, std::vector<NodesContainerType>& RecvObjects,Kratos::Serializer& particleSerializer);

    virtual bool TransferObjects(std::vector<ElementsContainerType>& SendObjects, std::vector<ElementsContainerType>& RecvObjects,Kratos::Serializer& particleSerializer);

    virtual bool TransferObjects(std::vector<ConditionsContainerType>& SendObjects, std::vector<ConditionsContainerType>& RecvObjects,Kratos::Serializer& particleSerializer);

    virtual bool SynchronizeOrNodalFlags(const Flags& TheFlags);

    virtual bool SynchronizeAndNodalFlags(const Flags& TheFlags);

    virtual bool SynchronizeNodalFlags();

    void Clear();

    ///@}
    ///@name Inquiry
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream, std::string const& rPrefixString="") const;

    ///@}
    ///@name Input and output
    ///@{

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


