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

#include "includes/communicator.h"
#include "includes/data_communicator.h"
#include "includes/parallel_environment.h"

namespace Kratos
{

// Life cycle and creation ////////////////////////////////////////////////////

// Default constructor
Communicator::Communicator()
    : mNumberOfColors(1)
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

// Constructor with custom DataCommunicator.
Communicator::Communicator(const DataCommunicator& rDataCommunicator)
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

// Copy constructor.
Communicator::Communicator(Communicator const& rOther)
    : mNumberOfColors(rOther.mNumberOfColors)
    , mNeighbourIndices(rOther.mNeighbourIndices)
    , mpLocalMesh(MeshType::Pointer(rOther.mpLocalMesh))
    , mpGhostMesh(MeshType::Pointer(rOther.mpGhostMesh))
    , mpInterfaceMesh(MeshType::Pointer(rOther.mpInterfaceMesh))
    , mLocalMeshes(rOther.mLocalMeshes)
    , mGhostMeshes(rOther.mGhostMeshes)
    , mInterfaceMeshes(rOther.mInterfaceMeshes)
    , mrDataCommunicator(rOther.mrDataCommunicator)
{}

Communicator::Pointer Communicator::Create(const DataCommunicator& rDataCommunicator) const
{
    KRATOS_TRY

    return Kratos::make_shared<Communicator>(rDataCommunicator);

    KRATOS_CATCH("");
}

Communicator::Pointer Communicator::Create() const
{
    return Kratos::make_shared<Communicator>();
}

// Public Access //////////////////////////////////////////////////////////////

bool Communicator::IsDistributed() const
{
    return false;
}

int Communicator::MyPID() const
{
    return mrDataCommunicator.Rank();
}

int Communicator::TotalProcesses() const
{
    return mrDataCommunicator.Size();
}

Communicator::SizeType Communicator::GlobalNumberOfNodes() const
{
    return mrDataCommunicator.SumAll(static_cast<unsigned int>(mpLocalMesh->NumberOfNodes()));
}

Communicator::SizeType Communicator::GlobalNumberOfElements() const
{
    return mrDataCommunicator.SumAll(static_cast<unsigned int>(mpLocalMesh->NumberOfElements()));
}

Communicator::SizeType Communicator::GlobalNumberOfConditions() const
{
    return mrDataCommunicator.SumAll(static_cast<unsigned int>(mpLocalMesh->NumberOfConditions()));
}

Communicator::SizeType Communicator::GetNumberOfColors() const
{
    return mNumberOfColors;
}

void Communicator::SetNumberOfColors(SizeType NewNumberOfColors)
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

void Communicator::AddColors(SizeType NumberOfAddedColors)
{
    if (NumberOfAddedColors < 1)
        return;

    mNumberOfColors += NumberOfAddedColors;
    MeshType mesh;

    for (IndexType i = 0; i < NumberOfAddedColors; i++)
    {
        mLocalMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
        mGhostMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
        mInterfaceMeshes.push_back(Kratos::make_shared<MeshType>(mesh.Clone()));
    }
}

Communicator::NeighbourIndicesContainerType& Communicator::NeighbourIndices()
{
    return mNeighbourIndices;
}

Communicator::NeighbourIndicesContainerType const& Communicator::NeighbourIndices() const
{
    return mNeighbourIndices;
}

// Set the local mesh pointer to the given mesh
void Communicator::SetLocalMesh(MeshType::Pointer pGivenMesh)
{
    mpLocalMesh = pGivenMesh;
}

// Returns pointer to the mesh storing all local entites
Communicator::MeshType::Pointer Communicator::pLocalMesh()
{
    return mpLocalMesh;
}

// Returns pointer to the mesh storing all ghost entites
Communicator::MeshType::Pointer Communicator::pGhostMesh()
{
    return mpGhostMesh;
}

// Returns pointer to the mesh storing all interface entites
Communicator::MeshType::Pointer Communicator::pInterfaceMesh()
{
    return mpInterfaceMesh;
}

// Returns a constant pointer to the mesh storing all local entites
const Communicator::MeshType::Pointer Communicator::pLocalMesh() const
{
    return mpLocalMesh;
}

// Returns a constant pointer to the mesh storing all ghost entites
const Communicator::MeshType::Pointer Communicator::pGhostMesh() const
{
    return mpGhostMesh;
}

// Returns a constant pointer to the mesh storing all interface entites
const Communicator::MeshType::Pointer Communicator::pInterfaceMesh() const
{
    return mpInterfaceMesh;
}

Communicator::MeshType::Pointer Communicator::pLocalMesh(IndexType ThisIndex)
{
    return mLocalMeshes(ThisIndex);
}

Communicator::MeshType::Pointer Communicator::pGhostMesh(IndexType ThisIndex)
{
    return mGhostMeshes(ThisIndex);
}

Communicator::MeshType::Pointer Communicator::pInterfaceMesh(IndexType ThisIndex)
{
    return mInterfaceMeshes(ThisIndex);
}

const Communicator::MeshType::Pointer Communicator::pLocalMesh(IndexType ThisIndex) const
{
    return mLocalMeshes(ThisIndex);
}

const Communicator::MeshType::Pointer Communicator::pGhostMesh(IndexType ThisIndex) const
{
    return mGhostMeshes(ThisIndex);
}

const Communicator::MeshType::Pointer Communicator::pInterfaceMesh(IndexType ThisIndex) const
{
    return mInterfaceMeshes(ThisIndex);
}

// Returns the reference to the mesh storing all local entites
Communicator::MeshType& Communicator::LocalMesh()
{
    return *mpLocalMesh;
}

// Returns the reference to the mesh storing all ghost entites
Communicator::MeshType& Communicator::GhostMesh()
{
    return *mpGhostMesh;
}

// Returns the reference to the mesh storing all interface entites
Communicator::MeshType& Communicator::InterfaceMesh()
{
    return *mpInterfaceMesh;
}

// Returns a constant reference to the mesh storing all local entites
Communicator::MeshType const& Communicator::LocalMesh() const
{
    return *mpLocalMesh;
}

// Returns a constant reference to the mesh storing all ghost entites
Communicator::MeshType const& Communicator::GhostMesh() const
{
    return *mpGhostMesh;
}

// Returns a constant reference to the mesh storing all interface entites
Communicator::MeshType const& Communicator::InterfaceMesh() const
{
    return *mpInterfaceMesh;
}

Communicator::MeshType& Communicator::LocalMesh(IndexType ThisIndex)
{
    return mLocalMeshes[ThisIndex];
}

Communicator::MeshType& Communicator::GhostMesh(IndexType ThisIndex)
{
    return mGhostMeshes[ThisIndex];
}

Communicator::MeshType& Communicator::InterfaceMesh(IndexType ThisIndex)
{
    return mInterfaceMeshes[ThisIndex];
}

Communicator::MeshType const& Communicator::LocalMesh(IndexType ThisIndex) const
{
    return mLocalMeshes[ThisIndex];
}

Communicator::MeshType const& Communicator::GhostMesh(IndexType ThisIndex) const
{
    return mGhostMeshes[ThisIndex];
}

Communicator::MeshType const& Communicator::InterfaceMesh(IndexType ThisIndex) const
{
    return mInterfaceMeshes[ThisIndex];
}

Communicator::MeshesContainerType& Communicator::LocalMeshes()
{
    return mLocalMeshes;
}

Communicator::MeshesContainerType& Communicator::GhostMeshes()
{
    return mGhostMeshes;
}

Communicator::MeshesContainerType& Communicator::InterfaceMeshes()
{
    return mInterfaceMeshes;
}

Communicator::MeshesContainerType const& Communicator::LocalMeshes() const
{
    return mLocalMeshes;
}

Communicator::MeshesContainerType const& Communicator::GhostMeshes() const
{
    return mGhostMeshes;
}

Communicator::MeshesContainerType const& Communicator::InterfaceMeshes() const
{
    return mInterfaceMeshes;
}

const DataCommunicator& Communicator::GetDataCommunicator() const
{
    return mrDataCommunicator;
}

// Public Operatrions /////////////////////////////////////////////////////////

bool Communicator::SynchronizeNodalSolutionStepsData()
{
    return true;
}

bool Communicator::SynchronizeDofs()
{
    return true;
}

bool Communicator::SynchronizeVariable(Variable<int> const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeVariable(Variable<double> const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeVariable(Variable<bool> const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeVariable(Variable<array_1d<double, 3 > > const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeVariable(Variable<array_1d<double, 4 > > const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeVariable(Variable<array_1d<double, 6 > > const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeVariable(Variable<array_1d<double, 9 > > const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeVariable(Variable<Vector> const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeVariable(Variable<Matrix> const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeNonHistoricalVariable(Variable<int> const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeNonHistoricalVariable(Variable<double> const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeNonHistoricalVariable(Variable<bool> const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeNonHistoricalVariable(Variable<array_1d<double, 3 > > const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeNonHistoricalVariable(Variable<array_1d<double, 4 > > const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeNonHistoricalVariable(Variable<array_1d<double, 6 > > const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeNonHistoricalVariable(Variable<array_1d<double, 9 > > const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeNonHistoricalVariable(Variable<Vector> const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeNonHistoricalVariable(Variable<Matrix> const& rThisVariable)
{
    return true;
}

bool Communicator::SynchronizeCurrentDataToMin(Variable<double> const& ThisVariable)
{
    return true;
}

bool Communicator::SynchronizeNonHistoricalDataToMin(Variable<double> const& ThisVariable)
{
    return true;
}

bool Communicator::SynchronizeElementalFlags()
{
    return true;
}

bool Communicator::AssembleCurrentData(Variable<int> const& ThisVariable)
{
    return true;
}

bool Communicator::AssembleCurrentData(Variable<double> const& ThisVariable)
{
    return true;
}

bool Communicator::AssembleCurrentData(Variable<array_1d<double, 3 > > const& ThisVariable)
{
    return true;
}

bool Communicator::AssembleCurrentData(Variable<Vector> const& ThisVariable)
{
    return true;
}

bool Communicator::AssembleCurrentData(Variable<Matrix> const& ThisVariable)
{
    return true;
}

bool Communicator::AssembleNonHistoricalData(Variable<int> const& ThisVariable)
{
    return true;
}

bool Communicator::AssembleNonHistoricalData(Variable<double> const& ThisVariable)
{
    return true;
}

bool Communicator::AssembleNonHistoricalData(Variable<array_1d<double, 3 > > const& ThisVariable)
{
    return true;
}


bool Communicator::AssembleNonHistoricalData(Variable<DenseVector<array_1d<double,3> > > const& ThisVariable)
{
    return true;
}

bool Communicator::AssembleNonHistoricalData(Variable<Vector> const& ThisVariable)
{
    return true;
}

bool Communicator::AssembleNonHistoricalData(Variable<Matrix> const& ThisVariable)
{
    return true;
}

bool Communicator::SynchronizeElementalNonHistoricalVariable(Variable<int> const& ThisVariable)
{
    return true;
}

bool Communicator::SynchronizeElementalNonHistoricalVariable(Variable<double> const& ThisVariable)
{
    return true;
}

bool Communicator::SynchronizeElementalNonHistoricalVariable(Variable<array_1d<double, 3 > > const& ThisVariable)
{
    return true;
}

bool Communicator::SynchronizeElementalNonHistoricalVariable(Variable<DenseVector<array_1d<double,3> > > const& ThisVariable)
{
    return true;
}

bool Communicator::SynchronizeElementalNonHistoricalVariable(Variable<DenseVector<int> > const& ThisVariable)
{
    return true;
}

bool Communicator::SynchronizeElementalNonHistoricalVariable(Variable<Vector> const& ThisVariable)
{
    return true;
}

bool Communicator::SynchronizeElementalNonHistoricalVariable(Variable<Matrix> const& ThisVariable)
{
    return true;
}

bool Communicator::TransferObjects(std::vector<NodesContainerType>& SendObjects, std::vector<NodesContainerType>& RecvObjects)
{
    return true;
}

bool Communicator::TransferObjects(std::vector<ElementsContainerType>& SendObjects, std::vector<ElementsContainerType>& RecvObjects)
{
    return true;
}

bool Communicator::TransferObjects(std::vector<ConditionsContainerType>& SendObjects, std::vector<ConditionsContainerType>& RecvObjects)
{
    return true;
}

bool Communicator::TransferObjects(std::vector<NodesContainerType>& SendObjects, std::vector<NodesContainerType>& RecvObjects,Kratos::Serializer& particleSerializer)
{
    return true;
}

bool Communicator::TransferObjects(std::vector<ElementsContainerType>& SendObjects, std::vector<ElementsContainerType>& RecvObjects,Kratos::Serializer& particleSerializer)
{
    return true;
}

bool Communicator::TransferObjects(std::vector<ConditionsContainerType>& SendObjects, std::vector<ConditionsContainerType>& RecvObjects,Kratos::Serializer& particleSerializer)
{
    return true;
}

bool Communicator::SynchronizeOrNodalFlags(const Flags& TheFlags)
{
    return true;
}

bool Communicator::SynchronizeAndNodalFlags(const Flags& TheFlags)
{
    return true;
}

bool Communicator::SynchronizeNodalFlags()
{
    return true;
}

void Communicator::Clear()
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

// Inquiry ////////////////////////////////////////////////////////////////////

std::string Communicator::Info() const
{
    return "Communicator";
}

void Communicator::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

void Communicator::PrintData(std::ostream& rOStream, std::string const& rPrefixString) const
{
    rOStream << rPrefixString << "    Local Mesh " << " : " << std::endl;
    mpLocalMesh->PrintData(rOStream, rPrefixString + "    ");
    rOStream << rPrefixString << "    Ghost Mesh " << " : " << std::endl;
    mpGhostMesh->PrintData(rOStream, rPrefixString + "    ");
    rOStream << rPrefixString << "    Interface Mesh " << " : " << std::endl;
    mpInterfaceMesh->PrintData(rOStream, rPrefixString + "    ");
}

} // namespace Kratos
