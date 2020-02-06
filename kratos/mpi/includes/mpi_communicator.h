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

#if !defined(KRATOS_MPI_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_MPI_COMMUNICATOR_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>

// External includes
#include "mpi.h"

// Project includes
#include "includes/data_communicator.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/mpi_serializer.h"
#include "mpi/includes/mpi_data_communicator.h"

#define CUSTOMTIMER 1

/* Timer defines */
#include "utilities/timer.h"
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

namespace MPIInternals
{

template<class ValueType> struct SendTraits;

template<> struct SendTraits<int>
{
    using SendType = int;
    using BufferType = std::vector<SendType>;
    constexpr static bool IsFixedSize = true;
    constexpr static std::size_t BlockSize = 1;
};

template<> struct SendTraits<bool>
{
    using SendType = int; // Note: sending bool as int to avoid the problems of std::vector<bool>
    using BufferType = std::vector<SendType>;
    constexpr static bool IsFixedSize = true;
    constexpr static std::size_t BlockSize = 1;
};

template<> struct SendTraits<double>
{
    using SendType = double;
    using BufferType = std::vector<SendType>;
    constexpr static bool IsFixedSize = true;
    constexpr static std::size_t BlockSize = 1;
};

template<std::size_t TDim> struct SendTraits< array_1d<double,TDim> >
{
    using SendType = double;
    using BufferType = std::vector<SendType>;
    constexpr static bool IsFixedSize = true;
    constexpr static std::size_t BlockSize = TDim;
};

template<> struct SendTraits< Kratos::Flags >
{
    using SendType = double;
    using BufferType = std::vector<SendType>;
    constexpr static bool IsFixedSize = true;
    constexpr static std::size_t BlockSize = sizeof(Kratos::Flags)/sizeof(double);
};

template<> struct SendTraits< Vector >
{
    using SendType = double;
    using BufferType = std::vector<SendType>;
    constexpr static bool IsFixedSize = false;

    static inline std::size_t GetMessageSize(const Vector& rValue)
    {
        return rValue.size();
    }
};

template<> struct SendTraits< Matrix >
{
    using SendType = double;
    using BufferType = std::vector<SendType>;
    constexpr static bool IsFixedSize = false;

    static inline std::size_t GetMessageSize(const Matrix& rValue)
    {
        return rValue.data().size();
    }
};


template<typename TVectorValue> struct SendTraits< DenseVector<TVectorValue> >
{
    using SendType = double;
    using BufferType = std::vector<SendType>;
    constexpr static bool IsFixedSize = false;

    static inline std::size_t GetMessageSize(const DenseVector<TVectorValue>& rValue)
    {
        return rValue.size()*sizeof(TVectorValue)/sizeof(double);
    }
};


template<> struct SendTraits< Kratos::VariablesListDataValueContainer >
{
    using SendType = double;
    using BufferType = std::string;
    constexpr static bool IsFixedSize = false;

    static inline std::size_t GetMessageSize(const Kratos::VariablesListDataValueContainer& rValue)
    {
        return (rValue.TotalSize()+1)*sizeof(char);
    }
};

template<> struct SendTraits< Node<3>::DofsContainerType >
{
    using SendType = int;
    using BufferType = std::vector<SendType>;
    constexpr static bool IsFixedSize = false;

    static inline std::size_t GetMessageSize(const Node<3>::DofsContainerType& rValue)
    {
        return rValue.size();
    }
};

template<typename ValueType> struct DirectCopyTransfer
{
    using SendType = typename SendTraits<ValueType>::SendType;

    static inline void WriteBuffer(const ValueType& rValue, SendType* pBuffer)
    {
        *(ValueType*)(pBuffer) = rValue;
    }

    static inline void ReadBuffer(const SendType* pBuffer, ValueType& rValue)
    {
        rValue = *(reinterpret_cast<const ValueType*>(pBuffer));
    }
};

template<typename ValueType> struct DynamicArrayTypeTransfer
{
    using SendType = typename SendTraits<ValueType>::SendType;

    static inline void WriteBuffer(const ValueType& rValue, SendType* pBuffer)
    {
        std::memcpy(pBuffer, &(rValue.data()[0]), rValue.data().size()*sizeof(double));
    }

    static inline void ReadBuffer(const SendType* pBuffer, ValueType& rValue)
    {
        std::memcpy(&(rValue.data()[0]), pBuffer, rValue.data().size()*sizeof(double));
    }
};

template<class TValue> struct SendTools;

template<>
struct SendTools<int> : public DirectCopyTransfer<int> {};

template<>
struct SendTools<double>: public DirectCopyTransfer<double> {};

template<>
struct SendTools<bool>: public DirectCopyTransfer<bool> {};

template<std::size_t TDim>
struct SendTools< array_1d<double,TDim> >: public DirectCopyTransfer<array_1d<double,TDim>> {};

template<>
struct SendTools< Kratos::Flags >: public DirectCopyTransfer<Kratos::Flags> {};

template<>
struct SendTools< Vector >: public DynamicArrayTypeTransfer<Vector> {};

template<>
struct SendTools< Matrix >: public DynamicArrayTypeTransfer<Matrix> {};

template<typename TVectorValue>
struct SendTools< DenseVector<TVectorValue> >
{
    using SendType = typename SendTraits< DenseVector<TVectorValue> >::SendType;

    static inline void WriteBuffer(const DenseVector<TVectorValue>& rValue, SendType* pBuffer)
    {
        std::memcpy(pBuffer, &(rValue.data()[0]), rValue.size()*sizeof(TVectorValue));
    }

    static inline void ReadBuffer(const SendType* pBuffer, DenseVector<TVectorValue>& rValue)
    {
        std::size_t position = 0;
        for (unsigned int i = 0; i < rValue.size(); i++)
        {
            rValue[i] = *(reinterpret_cast<const TVectorValue*>(pBuffer + position));
            position += sizeof(TVectorValue) / sizeof(double);
        }
    }
};

template<> struct SendTools< Kratos::VariablesListDataValueContainer >
{
    using SendType = typename SendTraits< Kratos::VariablesListDataValueContainer >::SendType;

    static inline void WriteBuffer(const Kratos::VariablesListDataValueContainer& rValue, SendType* pBuffer)
    {
        std::memcpy(pBuffer, rValue.Data(), rValue.TotalSize() * sizeof(double));
    }

    static inline void ReadBuffer(const SendType* pBuffer, Kratos::VariablesListDataValueContainer& rValue)
    {
        std::memcpy(rValue.Data(), pBuffer, rValue.TotalSize() * sizeof(double));
    }
};

template<> struct SendTools< Node<3>::DofsContainerType >
{
    using SendType = SendTraits< Node<3>::DofsContainerType >::SendType;

    static inline void WriteBuffer(const Node<3>::DofsContainerType& rValue, SendType* pBuffer)
    {
        unsigned int i = 0;
        for (auto i_dof = rValue.begin(); i_dof != rValue.end(); ++i_dof)
        {
            *(pBuffer + i) = (*i_dof)->EquationId();
            ++i;
        }
    }

    static inline void ReadBuffer(const SendType* pBuffer, Node<3>::DofsContainerType& rValue)
    {
        unsigned int i = 0;
        for (auto i_dof = rValue.begin(); i_dof != rValue.end(); ++i_dof)
        {
            (*i_dof)->SetEquationId(*(pBuffer + i));
            ++i;
        }
    }
};

template<
    class TDatabaseAccess,
    bool IsFixedSize = SendTraits<typename TDatabaseAccess::ValueType>::IsFixedSize >
struct BufferAllocation {
    using ValueType = typename TDatabaseAccess::ValueType;
    static inline std::size_t GetSendSize(TDatabaseAccess& rAccess, const Communicator::MeshType& rSourceMesh);
    static inline std::size_t GetSendSize(const ValueType& rValue);
};

template<class TDatabaseAccess>
struct BufferAllocation<TDatabaseAccess, true> {
    using ValueType = typename TDatabaseAccess::ValueType;
    static inline std::size_t GetSendSize(TDatabaseAccess& rAccess, const Communicator::MeshType& rSourceMesh)
    {
        const auto& r_container = rAccess.GetContainer(rSourceMesh);
        std::size_t num_objects = r_container.size();
        return num_objects * SendTraits<ValueType>::BlockSize;
    }

    static inline std::size_t GetSendSize(const ValueType&)
    {
        return SendTraits<ValueType>::BlockSize;
    }
};

template<class TDatabaseAccess>
struct BufferAllocation<TDatabaseAccess, false> {
    using ValueType = typename TDatabaseAccess::ValueType;
    static inline std::size_t GetSendSize(TDatabaseAccess& rAccess, const Communicator::MeshType& rSourceMesh)
    {
        const auto& r_container = rAccess.GetContainer(rSourceMesh);
        std::size_t buffer_size = 0;
        for (auto iter = r_container.begin(); iter != r_container.end(); ++iter)
        {
            buffer_size += MPIInternals::SendTraits<ValueType>::GetMessageSize(rAccess.GetValue(iter));
        }

        return buffer_size;
    }

    static inline std::size_t GetSendSize(const ValueType& rValue)
    {
        return MPIInternals::SendTraits<ValueType>::GetMessageSize(rValue);
    }
};

class NodalContainerAccess {
public:

    using ContainerType = Communicator::MeshType::NodesContainerType;
    using IteratorType = Communicator::MeshType::NodesContainerType::iterator;
    using ConstIteratorType = Communicator::MeshType::NodesContainerType::const_iterator;

    ContainerType& GetContainer(Communicator::MeshType& rMesh)
    {
        return rMesh.Nodes();
    }

    const ContainerType& GetContainer(const Communicator::MeshType& rMesh)
    {
        return rMesh.Nodes();
    }
};

class ElementalContainerAccess {
public:

    using ContainerType = Communicator::MeshType::ElementsContainerType;
    using IteratorType = Communicator::MeshType::ElementsContainerType::iterator;
    using ConstIteratorType = Communicator::MeshType::ElementsContainerType::const_iterator;

    ContainerType& GetContainer(Communicator::MeshType& rMesh)
    {
        return rMesh.Elements();
    }

    const ContainerType& GetContainer(const Communicator::MeshType& rMesh)
    {
        return rMesh.Elements();
    }
};

template<class TValue> class NodalSolutionStepValueAccess: public NodalContainerAccess
{
    const Variable<TValue>& mrVariable;

public:

    using ValueType = TValue;
    using SendType = typename SendTraits<TValue>::SendType;

    NodalSolutionStepValueAccess(const Variable<TValue>& mrVariable):
        NodalContainerAccess(),
        mrVariable(mrVariable)
    {}

    TValue& GetValue(IteratorType& iter)
    {
        return iter->FastGetSolutionStepValue(mrVariable);
    }

    const TValue& GetValue(const ConstIteratorType& iter)
    {
        return iter->FastGetSolutionStepValue(mrVariable);
    }
};

template<class TValue> class NodalDataAccess: public NodalContainerAccess
{
    const Variable<TValue>& mrVariable;

public:

    using ValueType = TValue;
    using SendType = typename SendTraits<TValue>::SendType;

    NodalDataAccess(const Variable<TValue>& mrVariable):
        NodalContainerAccess(),
        mrVariable(mrVariable)
    {}

    TValue& GetValue(IteratorType& iter)
    {
        return iter->GetValue(mrVariable);
    }

    const TValue& GetValue(const ConstIteratorType& iter)
    {
        return iter->GetValue(mrVariable);
    }
};

class NodalFlagsAccess: public NodalContainerAccess
{
public:

    const Kratos::Flags& mrMask;

    using ValueType = Kratos::Flags;
    using SendType = typename SendTraits<ValueType>::SendType;

    NodalFlagsAccess(const Kratos::Flags& rMask):
        NodalContainerAccess(),
        mrMask(rMask)
    {}

    Kratos::Flags& GetValue(IteratorType& iter)
    {
        return *iter;
    }

    const Kratos::Flags& GetValue(const ConstIteratorType& iter)
    {
        return *iter;
    }
};

class NodalSolutionStepDataAccess: public NodalContainerAccess
{
public:

    using ValueType = Kratos::VariablesListDataValueContainer;
    using SendType = typename SendTraits<ValueType>::SendType;

    ValueType& GetValue(IteratorType& iter)
    {
        return iter->SolutionStepData();
    }

    const ValueType& GetValue(const ConstIteratorType& iter)
    {
        return iter->SolutionStepData();
    }
};

class DofIdAccess: public NodalContainerAccess
{
public:

    using ValueType = Node<3>::DofsContainerType;
    using SendType = typename SendTraits<ValueType>::SendType;

    ValueType& GetValue(IteratorType& iter)
    {
        return iter->GetDofs();
    }

    const ValueType& GetValue(const ConstIteratorType& iter)
    {
        return iter->GetDofs();
    }
};

template<class TValue> class ElementalDataAccess: public ElementalContainerAccess
{
    const Variable<TValue>& mrVariable;

public:

    using ValueType = TValue;
    using SendType = typename SendTraits<TValue>::SendType;

    ElementalDataAccess(const Variable<TValue>& mrVariable):
        ElementalContainerAccess(),
        mrVariable(mrVariable)
    {}

    TValue& GetValue(IteratorType& iter)
    {
        return iter->GetValue(mrVariable);
    }

    const TValue& GetValue(const ConstIteratorType& iter)
    {
        return iter->GetValue(mrVariable);
    }
};

class ElementalFlagsAccess: public ElementalContainerAccess
{
public:

    const Kratos::Flags& mrMask;

    using ValueType = Kratos::Flags;
    using SendType = typename SendTraits<ValueType>::SendType;

    ElementalFlagsAccess(const Kratos::Flags& rMask):
        ElementalContainerAccess(),
        mrMask(rMask)
    {}

    Kratos::Flags& GetValue(IteratorType& iter)
    {
        return *iter;
    }

    const Kratos::Flags& GetValue(const ConstIteratorType& iter)
    {
        return *iter;
    }
};

}

/// MPICommunicator manages the transfer of ModelPart data in MPI distributed memory environment.
class MPICommunicator : public Communicator
{

enum class DistributedType {
    Local,
    Ghost
};

// Auxiliary type for compile-time dispatch of local/ghost mesh access
template<DistributedType TDistributed> struct MeshAccess
{
    constexpr static DistributedType Value = TDistributed;
};

enum class OperationType {
    Replace,
    SumValues,
    MinValues,
    OrAccessedFlags,
    AndAccessedFlags,
    ReplaceAccessedFlags
};

// Auxiliary type for compile-time dispatch of the reduction operation in data transfer methods
template<OperationType TOperation> struct Operation
{
    constexpr static OperationType Value = TOperation;
};

public:
    ///@name  Enum's
    ///@{


    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MPICommunicator
    KRATOS_CLASS_POINTER_DEFINITION(MPICommunicator);

    typedef Communicator BaseType;

    typedef BaseType::IndexType IndexType;

    typedef BaseType::SizeType SizeType;

    typedef BaseType::NodeType NodeType;

    typedef BaseType::PropertiesType PropertiesType;

    typedef BaseType::ElementType ElementType;

    typedef BaseType::ConditionType ConditionType;

    typedef BaseType::NeighbourIndicesContainerType NeighbourIndicesContainerType;

    typedef BaseType::MeshType MeshType;

    typedef BaseType::MeshesContainerType MeshesContainerType;

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

    /// Constructor using the VariablesList of the ModelPart that will use this communicator.
    MPICommunicator(VariablesList* Variables_list) : BaseType(DataCommunicator::GetDefault()), mpVariables_list(Variables_list)
    {}

    /// Constructor using the VariablesList and a custom DataCommunicator.
    /** @param pVariablesList Pointer to the VariablesList of the ModelPart that will use this communicator.
     *  @param rDataCommunicator Reference to a DataCommunicator, which will be used for MPI calls.
     */
    MPICommunicator(VariablesList* pVariablesList, const DataCommunicator& rDataCommunicator)
    : BaseType(rDataCommunicator)
    , mpVariables_list(pVariablesList)
    {}

    /// Copy constructor.

    MPICommunicator(MPICommunicator const& rOther)
    : BaseType(rOther)
    , mpVariables_list(rOther.mpVariables_list)
    {}



    Communicator::Pointer Create() const override
    {
        return Create(DataCommunicator::GetDefault());
    }

    Communicator::Pointer Create(const DataCommunicator& rDataCommunicator) const override
    {
        KRATOS_TRY

        return Kratos::make_shared<MPICommunicator>(mpVariables_list, rDataCommunicator);

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~MPICommunicator() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MPICommunicator & operator=(MPICommunicator const& rOther) = delete;

    ///@}
    ///@name Access
    ///@{

    bool IsDistributed() const override
    {
        return true;
    }

    ///@}
    ///@name Operations
    ///@{

    bool SynchronizeNodalSolutionStepsData() override
    {
        constexpr MeshAccess<DistributedType::Local> local_meshes;
        constexpr MeshAccess<DistributedType::Ghost> ghost_meshes;
        constexpr Operation<OperationType::Replace> replace;
        MPIInternals::NodalSolutionStepDataAccess nodal_solution_step_access;

        TransferDistributedValuesUnknownSize(local_meshes, ghost_meshes, nodal_solution_step_access, replace);

        return true;
    }

    bool SynchronizeDofs() override
    {
        constexpr MeshAccess<DistributedType::Local> local_meshes;
        constexpr MeshAccess<DistributedType::Ghost> ghost_meshes;
        constexpr Operation<OperationType::Replace> replace;
        MPIInternals::DofIdAccess dof_id_access;

        TransferDistributedValues(local_meshes, ghost_meshes, dof_id_access, replace);

        return true;
    }

    bool SynchronizeVariable(Variable<int> const& rThisVariable) override
    {
        MPIInternals::NodalSolutionStepValueAccess<int> solution_step_value_access(rThisVariable);
        SynchronizeFixedSizeValues(solution_step_value_access);
        return true;
    }

    bool SynchronizeVariable(Variable<double> const& rThisVariable) override
    {
        MPIInternals::NodalSolutionStepValueAccess<double> solution_step_value_access(rThisVariable);
        SynchronizeFixedSizeValues(solution_step_value_access);
        return true;
    }

    bool SynchronizeVariable(Variable<bool> const& rThisVariable) override
    {
        MPIInternals::NodalSolutionStepValueAccess<bool> solution_step_value_access(rThisVariable);
        SynchronizeFixedSizeValues(solution_step_value_access);
        return true;
    }

    bool SynchronizeVariable(Variable<array_1d<double, 3 > > const& rThisVariable) override
    {
        MPIInternals::NodalSolutionStepValueAccess<array_1d<double,3>> solution_step_value_access(rThisVariable);
        SynchronizeFixedSizeValues(solution_step_value_access);
        return true;
    }

    bool SynchronizeVariable(Variable<array_1d<double, 4 > > const& rThisVariable) override
    {
        MPIInternals::NodalSolutionStepValueAccess<array_1d<double,4>> solution_step_value_access(rThisVariable);
        SynchronizeFixedSizeValues(solution_step_value_access);
        return true;
    }

    bool SynchronizeVariable(Variable<array_1d<double, 6 > > const& rThisVariable) override
    {
        MPIInternals::NodalSolutionStepValueAccess<array_1d<double,6>> solution_step_value_access(rThisVariable);
        SynchronizeFixedSizeValues(solution_step_value_access);
        return true;
    }

    bool SynchronizeVariable(Variable<array_1d<double, 9 > > const& rThisVariable) override
    {
        MPIInternals::NodalSolutionStepValueAccess<array_1d<double,9>> solution_step_value_access(rThisVariable);
        SynchronizeFixedSizeValues(solution_step_value_access);
        return true;
    }

    bool SynchronizeVariable(Variable<Vector> const& rThisVariable) override
    {
        MPIInternals::NodalSolutionStepValueAccess<Vector> solution_step_value_access(rThisVariable);
        SynchronizeDynamicVectorValues(solution_step_value_access);
        return true;
    }

    bool SynchronizeVariable(Variable<Matrix> const& rThisVariable) override
    {
        MPIInternals::NodalSolutionStepValueAccess<Matrix> solution_step_value_access(rThisVariable);
        SynchronizeDynamicMatrixValues(solution_step_value_access);
        return true;
    }

    bool SynchronizeNonHistoricalVariable(Variable<int> const& rThisVariable) override
    {
        MPIInternals::NodalDataAccess<int> nodal_data_access(rThisVariable);
        SynchronizeFixedSizeValues(nodal_data_access);
        return true;
    }

    bool SynchronizeNonHistoricalVariable(Variable<double> const& rThisVariable) override
    {
        MPIInternals::NodalDataAccess<double> nodal_data_access(rThisVariable);
        SynchronizeFixedSizeValues(nodal_data_access);
        return true;
    }

    bool SynchronizeNonHistoricalVariable(Variable<bool> const& rThisVariable) override
    {
        MPIInternals::NodalDataAccess<bool> nodal_data_access(rThisVariable);
        SynchronizeFixedSizeValues(nodal_data_access);
        return true;
    }

    bool SynchronizeNonHistoricalVariable(Variable<array_1d<double, 3 > > const& rThisVariable) override
    {
        MPIInternals::NodalDataAccess<array_1d<double,3>> nodal_data_access(rThisVariable);
        SynchronizeFixedSizeValues(nodal_data_access);
        return true;
    }

    bool SynchronizeNonHistoricalVariable(Variable<array_1d<double, 4 > > const& rThisVariable) override
    {
        MPIInternals::NodalDataAccess<array_1d<double,4>> nodal_data_access(rThisVariable);
        SynchronizeFixedSizeValues(nodal_data_access);
        return true;
    }

    bool SynchronizeNonHistoricalVariable(Variable<array_1d<double, 6 > > const& rThisVariable) override
    {
        MPIInternals::NodalDataAccess<array_1d<double,6>> nodal_data_access(rThisVariable);
        SynchronizeFixedSizeValues(nodal_data_access);
        return true;
    }

    bool SynchronizeNonHistoricalVariable(Variable<array_1d<double, 9 > > const& rThisVariable) override
    {
        MPIInternals::NodalDataAccess<array_1d<double,9>> nodal_data_access(rThisVariable);
        SynchronizeFixedSizeValues(nodal_data_access);
        return true;
    }

    bool SynchronizeNonHistoricalVariable(Variable<Vector> const& rThisVariable) override
    {
        MPIInternals::NodalDataAccess<Vector> nodal_data_access(rThisVariable);
        SynchronizeDynamicVectorValues(nodal_data_access);
        return true;
    }

    bool SynchronizeNonHistoricalVariable(Variable<Matrix> const& rThisVariable) override
    {
        MPIInternals::NodalDataAccess<Matrix> nodal_data_access(rThisVariable);
        SynchronizeDynamicMatrixValues(nodal_data_access);
        return true;
    }

    bool SynchronizeCurrentDataToMin(Variable<double> const& ThisVariable) override
    {
        constexpr MeshAccess<DistributedType::Local> local_meshes;
        constexpr MeshAccess<DistributedType::Ghost> ghost_meshes;
        constexpr Operation<OperationType::Replace> replace;
        constexpr Operation<OperationType::MinValues> min;
        MPIInternals::NodalSolutionStepValueAccess<double> nodal_solution_step_access(ThisVariable);

        // Calculate min on owner rank
        TransferDistributedValues(ghost_meshes, local_meshes, nodal_solution_step_access, min);

        // Synchronize result on ghost copies
        TransferDistributedValues(local_meshes, ghost_meshes, nodal_solution_step_access, replace);

        return true;
    }

    bool SynchronizeNonHistoricalDataToMin(Variable<double> const& ThisVariable) override
    {
        constexpr MeshAccess<DistributedType::Local> local_meshes;
        constexpr MeshAccess<DistributedType::Ghost> ghost_meshes;
        constexpr Operation<OperationType::Replace> replace;
        constexpr Operation<OperationType::MinValues> min;
        MPIInternals::NodalDataAccess<double> nodal_data_access(ThisVariable);

        // Calculate min on owner rank
        TransferDistributedValues(ghost_meshes, local_meshes, nodal_data_access, min);

        // Synchronize result on ghost copies
        TransferDistributedValues(local_meshes, ghost_meshes, nodal_data_access, replace);

        return true;
    }

    bool AssembleCurrentData(Variable<int> const& ThisVariable) override
    {
        MPIInternals::NodalSolutionStepValueAccess<int> solution_step_access(ThisVariable);
        AssembleFixedSizeValues(solution_step_access);
        return true;
    }

    bool AssembleCurrentData(Variable<double> const& ThisVariable) override
    {
        MPIInternals::NodalSolutionStepValueAccess<double> solution_step_access(ThisVariable);
        AssembleFixedSizeValues(solution_step_access);
        return true;
    }

    bool AssembleCurrentData(Variable<array_1d<double, 3 > > const& ThisVariable) override
    {
        MPIInternals::NodalSolutionStepValueAccess<array_1d<double,3>> solution_step_access(ThisVariable);
        AssembleFixedSizeValues(solution_step_access);
        return true;
    }

    bool AssembleCurrentData(Variable<Vector> const& ThisVariable) override
    {
        MPIInternals::NodalSolutionStepValueAccess<Vector> solution_step_access(ThisVariable);
        AssembleDynamicVectorValues(solution_step_access);
        return true;
    }

    bool AssembleCurrentData(Variable<Matrix> const& ThisVariable) override
    {
        MPIInternals::NodalSolutionStepValueAccess<Matrix> solution_step_access(ThisVariable);
        AssembleDynamicMatrixValues(solution_step_access);
        return true;
    }

    bool AssembleNonHistoricalData(Variable<int> const& ThisVariable) override
    {
        MPIInternals::NodalDataAccess<int> nodal_data_access(ThisVariable);
        AssembleFixedSizeValues(nodal_data_access);
        return true;
    }

    bool AssembleNonHistoricalData(Variable<double> const& ThisVariable) override
    {
        MPIInternals::NodalDataAccess<double> nodal_data_access(ThisVariable);
        AssembleFixedSizeValues(nodal_data_access);
        return true;
    }

    bool AssembleNonHistoricalData(Variable<array_1d<double, 3 > > const& ThisVariable) override
    {
        MPIInternals::NodalDataAccess<array_1d<double,3>> nodal_data_access(ThisVariable);
        AssembleFixedSizeValues(nodal_data_access);
        return true;
    }

    bool AssembleNonHistoricalData(Variable<DenseVector<array_1d<double,3> > > const& ThisVariable) override
    {
        MPIInternals::NodalDataAccess<DenseVector<array_1d<double,3>>> nodal_data_access(ThisVariable);
        AssembleDynamicVectorValues(nodal_data_access);
        return true;
    }

    bool AssembleNonHistoricalData(Variable<Vector> const& ThisVariable) override
    {
        MPIInternals::NodalDataAccess<Vector> nodal_data_access(ThisVariable);
        AssembleDynamicVectorValues(nodal_data_access);
        return true;
    }

    bool AssembleNonHistoricalData(Variable<Matrix> const& ThisVariable) override
    {
        MPIInternals::NodalDataAccess<Matrix> nodal_data_access(ThisVariable);
        AssembleDynamicMatrixValues(nodal_data_access);
        return true;
    }

    /////////////////////////////////////////////////////////////////////////////

    bool SynchronizeElementalNonHistoricalVariable(Variable<int> const& ThisVariable) override
    {
        MPIInternals::ElementalDataAccess<int> elemental_data_access(ThisVariable);
        SynchronizeFixedSizeValues(elemental_data_access);
        return true;
    }

    bool SynchronizeElementalNonHistoricalVariable(Variable<double> const& ThisVariable) override
    {
        MPIInternals::ElementalDataAccess<double> elemental_data_access(ThisVariable);
        SynchronizeFixedSizeValues(elemental_data_access);
        return true;
    }

    bool SynchronizeElementalNonHistoricalVariable(Variable<array_1d<double, 3 > > const& ThisVariable) override
    {
        MPIInternals::ElementalDataAccess<array_1d<double,3>> elemental_data_access(ThisVariable);
        SynchronizeFixedSizeValues(elemental_data_access);
        return true;
    }

    bool SynchronizeElementalNonHistoricalVariable(Variable<DenseVector<array_1d<double,3> > > const& ThisVariable) override
    {
        MPIInternals::ElementalDataAccess<DenseVector<array_1d<double,3>>> elemental_data_access(ThisVariable);
        AssembleDynamicVectorValues(elemental_data_access);
        return true;
    }

    bool SynchronizeElementalNonHistoricalVariable(Variable<DenseVector<int> > const& ThisVariable) override
    {
        MPIInternals::ElementalDataAccess<DenseVector<int>> elemental_data_access(ThisVariable);
        AssembleDynamicVectorValues(elemental_data_access);
        return true;
    }

    bool SynchronizeElementalNonHistoricalVariable(Variable<Vector> const& ThisVariable) override
    {
        MPIInternals::ElementalDataAccess<Vector> elemental_data_access(ThisVariable);
        SynchronizeDynamicVectorValues(elemental_data_access);
        return true;
    }

    bool SynchronizeElementalNonHistoricalVariable(Variable<Matrix> const& ThisVariable) override
    {
        MPIInternals::ElementalDataAccess<Matrix> elemental_data_access(ThisVariable);
        SynchronizeDynamicMatrixValues(elemental_data_access);
        return true;
    }

    /////////////////////////////////////////////////////////////////////////////

    /**
     * Transfer objects from a given process to a destination process
     * @param SendObjects list of objects to be send.      SendObjects[i] -> Objects to   process i
     * @param RecvObjects list of objects to be received.  RecvObjects[i] -> objects from process i
     **/
    bool TransferObjects(std::vector<NodesContainerType>& SendObjects, std::vector<NodesContainerType>& RecvObjects) override
    {
        AsyncSendAndReceiveObjects<NodesContainerType>(SendObjects,RecvObjects);
        return true;
    }

    /**
    * Transfer objects from a given process to a destination process
    * @param SendObjects list of objects to be send.      SendObjects[i] -> Objects to   process i
    * @param RecvObjects list of objects to be received.  RecvObjects[i] -> objects from process i
    **/
    bool TransferObjects(std::vector<ElementsContainerType>& SendObjects, std::vector<ElementsContainerType>& RecvObjects) override
    {
        AsyncSendAndReceiveObjects<ElementsContainerType>(SendObjects,RecvObjects);
        return true;
    }

    /**
    * Transfer objects from a given process to a destination process
    * @param SendObjects list of objects to be send.      SendObjects[i] -> Objects to   process i
    * @param RecvObjects list of objects to be received.  RecvObjects[i] -> objects from process i
    **/
    bool TransferObjects(std::vector<ConditionsContainerType>& SendObjects, std::vector<ConditionsContainerType>& RecvObjects) override
    {
        AsyncSendAndReceiveObjects<ConditionsContainerType>(SendObjects,RecvObjects);
        return true;
    }

    /**
     * Transfer objects from a given process to a destination process
     * @param SendObjects list of objects to be send.      SendObjects[i] -> Objects to   process i
     * @param RecvObjects list of objects to be received.  RecvObjects[i] -> objects from process i
     **/
    bool TransferObjects(std::vector<NodesContainerType>& SendObjects, std::vector<NodesContainerType>& RecvObjects,Kratos::Serializer& particleSerializer) override
    {
        AsyncSendAndReceiveObjects<NodesContainerType>(SendObjects,RecvObjects);
        return true;
    }

    /**
    * Transfer objects from a given process to a destination process
    * @param SendObjects list of objects to be send.      SendObjects[i] -> Objects to   process i
    * @param RecvObjects list of objects to be received.  RecvObjects[i] -> objects from process i
    **/
    bool TransferObjects(std::vector<ElementsContainerType>& SendObjects, std::vector<ElementsContainerType>& RecvObjects,Kratos::Serializer& particleSerializer) override
    {
        AsyncSendAndReceiveObjects<ElementsContainerType>(SendObjects,RecvObjects);
        return true;
    }

    /**
    * Transfer objects from a given process to a destination process
    * @param SendObjects list of objects to be send.      SendObjects[i] -> Objects to   process i
    * @param RecvObjects list of objects to be received.  RecvObjects[i] -> objects from process i
    **/
    bool TransferObjects(std::vector<ConditionsContainerType>& SendObjects, std::vector<ConditionsContainerType>& RecvObjects,Kratos::Serializer& particleSerializer) override
    {
        AsyncSendAndReceiveObjects<ConditionsContainerType>(SendObjects,RecvObjects);
        return true;
    }

    /**
     * Assemble the values of the chosen flags on each node with an OR operation.
     * @param[in] TheFlags Names of the flags to be synchronized.
     */
    bool SynchronizeOrNodalFlags(const Flags& TheFlags) override
    {
        constexpr MeshAccess<DistributedType::Local> local_meshes;
        constexpr MeshAccess<DistributedType::Ghost> ghost_meshes;
        MPIInternals::NodalFlagsAccess nodal_flag_access(TheFlags);
        constexpr Operation<OperationType::OrAccessedFlags> or_accessed;
        constexpr Operation<OperationType::ReplaceAccessedFlags> replace_accessed;

        TransferDistributedValues(ghost_meshes, local_meshes, nodal_flag_access, or_accessed);

        TransferDistributedValues(local_meshes, ghost_meshes, nodal_flag_access, replace_accessed);
        return true;
    }

    /**
     * Assemble the values of the chosen flags on each node with an AND operation.
     * @param[in] TheFlags Names of the flags to be synchronized.
     */
    bool SynchronizeAndNodalFlags(const Flags& TheFlags) override
    {
        constexpr MeshAccess<DistributedType::Local> local_meshes;
        constexpr MeshAccess<DistributedType::Ghost> ghost_meshes;
        MPIInternals::NodalFlagsAccess nodal_flag_access(TheFlags);
        constexpr Operation<OperationType::AndAccessedFlags> and_accessed;
        constexpr Operation<OperationType::ReplaceAccessedFlags> replace_accessed;

        TransferDistributedValues(ghost_meshes, local_meshes, nodal_flag_access, and_accessed);

        TransferDistributedValues(local_meshes, ghost_meshes, nodal_flag_access, replace_accessed);
        return true;
    }

    bool SynchronizeNodalFlags() override
    {
        constexpr MeshAccess<DistributedType::Local> local_meshes;
        constexpr MeshAccess<DistributedType::Ghost> ghost_meshes;
        MPIInternals::NodalFlagsAccess nodal_flags_access(Flags::AllDefined());
        constexpr Operation<OperationType::Replace> replace;

        TransferDistributedValues(local_meshes, ghost_meshes, nodal_flags_access, replace);
        return true;
    }

    bool SynchronizeElementalFlags() override
    {
        constexpr MeshAccess<DistributedType::Local> local_meshes;
        constexpr MeshAccess<DistributedType::Ghost> ghost_meshes;
        MPIInternals::ElementalFlagsAccess elemental_flags_access(Flags::AllDefined());
        constexpr Operation<OperationType::Replace> replace;

        TransferDistributedValues(local_meshes, ghost_meshes, elemental_flags_access, replace);
        return true;
    }

    /////////////////////////////////////////////////////////////////////////////

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

    std::string Info() const override
    {
        return "MPICommunicator";
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.

    void PrintData(std::ostream& rOStream) const override
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

    VariablesList* mpVariables_list;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void PrintNodesId()
    {
        NeighbourIndicesContainerType& neighbours_indices = NeighbourIndices();

        int nproc = mrDataCommunicator.Size();
        int rank = mrDataCommunicator.Rank();

        mrDataCommunicator.Barrier();
        for (int proc_id = 0; proc_id < nproc; proc_id++)
        {
            if (proc_id == rank)
            {

                for (int i_color = 0; i_color < static_cast<int>(neighbours_indices.size()); i_color++)
                {
                    if ((neighbours_indices[i_color]) >= 0)
                    {
                        NodesContainerType& r_local_nodes = LocalMesh(i_color).Nodes();
                        NodesContainerType& r_ghost_nodes = GhostMesh(i_color).Nodes();
                        std::string tag = "Local nodes in rank ";
                        PrintNodesId(r_local_nodes, tag, i_color);
                        tag = "Ghost nodes in rank ";
                        PrintNodesId(r_ghost_nodes, tag, i_color);
                        tag = "Interface nodes in rank ";
                        PrintNodesId(InterfaceMesh(i_color).Nodes(), tag, i_color);

                    }
                }
            }
            mrDataCommunicator.Barrier();
        }
    }

    template<class TNodesArrayType>
    void PrintNodesId(TNodesArrayType& rNodes, std::string Tag, int color)
    {
        int rank = mrDataCommunicator.Rank();
        std::cout << Tag << rank << " with color " << color << ":";
        for (typename TNodesArrayType::iterator i_node = rNodes.begin(); i_node != rNodes.end(); i_node++)
            std::cout << i_node->Id() << ", ";

        std::cout << std::endl;
    }

    template<class TObjectType>
    bool AsyncSendAndReceiveObjects(std::vector<TObjectType>& SendObjects, std::vector<TObjectType>& RecvObjects)
    {
        int mpi_rank = mrDataCommunicator.Rank();
        int mpi_size = mrDataCommunicator.Size();

        MPI_Comm comm = MPIDataCommunicator::GetMPICommunicator(mrDataCommunicator);

        int * msgSendSize = new int[mpi_size];
        int * msgRecvSize = new int[mpi_size];

        char ** message = new char * [mpi_size];
        char ** mpi_send_buffer = new char * [mpi_size];

        for(int i = 0; i < mpi_size; i++)
        {
            msgSendSize[i] = 0;
            msgRecvSize[i] = 0;
        }

        for(int i = 0; i < mpi_size; i++)
        {
            if(mpi_rank != i)
            {
                Kratos::MpiSerializer serializer;

                serializer.save("VariableList",mpVariables_list);
                serializer.save("ObjectList",SendObjects[i].GetContainer());

                std::stringstream * stream = (std::stringstream *)serializer.pGetBuffer();
                const std::string & stream_str = stream->str();
                const char * cstr = stream_str.c_str();

                msgSendSize[i] = sizeof(char) * (stream_str.size()+1);
                mpi_send_buffer[i] = (char *)malloc(msgSendSize[i]);
                memcpy(mpi_send_buffer[i],cstr,msgSendSize[i]);
            }
        }

        MPI_Alltoall(msgSendSize,1,MPI_INT,msgRecvSize,1,MPI_INT,comm);

        int NumberOfCommunicationEvents      = 0;
        int NumberOfCommunicationEventsIndex = 0;

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
                message[i] = (char *)malloc(sizeof(char) * msgRecvSize[i]);

                MPI_Irecv(message[i],msgRecvSize[i],MPI_CHAR,i,0,comm,&reqs[NumberOfCommunicationEventsIndex++]);
            }

            if(i != mpi_rank && msgSendSize[i])
            {
                MPI_Isend(mpi_send_buffer[i],msgSendSize[i],MPI_CHAR,i,0,comm,&reqs[NumberOfCommunicationEventsIndex++]);
            }
        }

        //wait untill all communications finish
        int err = MPI_Waitall(NumberOfCommunicationEvents, reqs, stats);

        KRATOS_ERROR_IF(err != MPI_SUCCESS) << "Error in MPICommunicator asynchronous data transfer" << std::endl;

        mrDataCommunicator.Barrier();

        for(int i = 0; i < mpi_size; i++)
        {
            if (i != mpi_rank && msgRecvSize[i])
            {
                Kratos::MpiSerializer serializer;
                std::stringstream * serializer_buffer;

                serializer_buffer = (std::stringstream *)serializer.pGetBuffer();
                serializer_buffer->write(message[i], msgRecvSize[i]);

                VariablesList* tmp_mpVariables_list = NULL;

                serializer.load("VariableList",tmp_mpVariables_list);

                if(tmp_mpVariables_list != NULL)
                  delete tmp_mpVariables_list;
                tmp_mpVariables_list = mpVariables_list;

                serializer.load("ObjectList",RecvObjects[i].GetContainer());
            }

            mrDataCommunicator.Barrier();
        }

        // Free buffers
        for(int i = 0; i < mpi_size; i++)
        {
            if(msgRecvSize[i])
                free(message[i]);

            if(msgSendSize[i])
                free(mpi_send_buffer[i]);
        }

        delete [] reqs;
        delete [] stats;

        delete [] message;
        delete [] mpi_send_buffer;

        delete [] msgSendSize;
        delete [] msgRecvSize;

        return true;
    }

    template<class TDatabaseAccess>
    void SynchronizeFixedSizeValues(TDatabaseAccess& rVariableAccess)
    {
        constexpr MeshAccess<DistributedType::Local> local_meshes;
        constexpr MeshAccess<DistributedType::Ghost> ghost_meshes;
        constexpr Operation<OperationType::Replace> replace;

        TransferDistributedValues(local_meshes, ghost_meshes, rVariableAccess, replace);
    }

    template<class TDatabaseAccess>
    void SynchronizeDynamicVectorValues(TDatabaseAccess& rVariableAccess)
    {
        constexpr MeshAccess<DistributedType::Local> local_meshes;
        constexpr MeshAccess<DistributedType::Ghost> ghost_meshes;
        constexpr Operation<OperationType::Replace> replace;

        // Communicate vector sizes to ghost copies
        MatchDynamicVectorSizes(local_meshes, ghost_meshes, rVariableAccess);

        // Synchronize values on ghost copies
        TransferDistributedValues(local_meshes, ghost_meshes, rVariableAccess, replace);
    }

    template<class TDatabaseAccess>
    void SynchronizeDynamicMatrixValues(TDatabaseAccess& rVariableAccess)
    {
        constexpr MeshAccess<DistributedType::Local> local_meshes;
        constexpr MeshAccess<DistributedType::Ghost> ghost_meshes;
        constexpr Operation<OperationType::Replace> replace;

        // Communicate matrix sizes to ghost copies
        MatchDynamicMatrixSizes(local_meshes, ghost_meshes, rVariableAccess);

        // Synchronize values on ghost copies
        TransferDistributedValues(local_meshes, ghost_meshes, rVariableAccess, replace);
    }

    template<class TDatabaseAccess>
    void AssembleFixedSizeValues(TDatabaseAccess& rVariableAccess)
    {
        constexpr MeshAccess<DistributedType::Local> local_meshes;
        constexpr MeshAccess<DistributedType::Ghost> ghost_meshes;
        constexpr Operation<OperationType::Replace> replace;
        constexpr Operation<OperationType::SumValues> sum;

        // Assemble results on owner rank
        TransferDistributedValues(ghost_meshes, local_meshes, rVariableAccess, sum);

        // Synchronize result on ghost copies
        TransferDistributedValues(local_meshes, ghost_meshes, rVariableAccess, replace);
    }

    template<class TDatabaseAccess>
    void AssembleDynamicVectorValues(TDatabaseAccess& rVariableAccess)
    {
        constexpr MeshAccess<DistributedType::Local> local_meshes;
        constexpr MeshAccess<DistributedType::Ghost> ghost_meshes;
        constexpr Operation<OperationType::Replace> replace;
        constexpr Operation<OperationType::SumValues> sum;

        // Combine vector sizes on owner rank
        MatchDynamicVectorSizes(ghost_meshes, local_meshes, rVariableAccess);

        // Communicate vector sizes to ghost copies
        MatchDynamicVectorSizes(local_meshes, ghost_meshes, rVariableAccess);

        // From this point on, we can assume buffer sizes will always match for all ranks

        // Assemble results on owner rank
        TransferDistributedValues(ghost_meshes, local_meshes, rVariableAccess, sum);

        // Synchronize result on ghost copies
        TransferDistributedValues(local_meshes, ghost_meshes, rVariableAccess, replace);
    }

    template<class TDatabaseAccess>
    void AssembleDynamicMatrixValues(TDatabaseAccess& rVariableAccess)
    {
        constexpr MeshAccess<DistributedType::Local> local_meshes;
        constexpr MeshAccess<DistributedType::Ghost> ghost_meshes;
        constexpr Operation<OperationType::Replace> replace;
        constexpr Operation<OperationType::SumValues> sum;

        // Combine matrix sizes on owner rank
        MatchDynamicMatrixSizes(ghost_meshes, local_meshes, rVariableAccess);

        // Communicate matrix sizes to ghost copies
        MatchDynamicMatrixSizes(local_meshes, ghost_meshes, rVariableAccess);

        // From this point on, we can assume buffer sizes will always match for all ranks

        // Assemble results on owner rank
        TransferDistributedValues(ghost_meshes, local_meshes, rVariableAccess, sum);

        // Synchronize result on ghost copies
        TransferDistributedValues(local_meshes, ghost_meshes, rVariableAccess, replace);
    }

    MeshType& GetMesh(IndexType Color, const MeshAccess<DistributedType::Local>)
    {
        return LocalMesh(Color);
    }

    MeshType& GetMesh(IndexType Color, const MeshAccess<DistributedType::Ghost>)
    {
        return GhostMesh(Color);
    }

    template<class TDatabaseAccess>
    std::size_t ReduceValues(
        const typename TDatabaseAccess::SendType* pBuffer,
        TDatabaseAccess& rAccess,
        typename TDatabaseAccess::IteratorType ContainerIterator,
        Operation<OperationType::Replace>)
    {
        // Replace the value by reading the sent buffer in place.
        // Note that dynamic types are assumed to already have the
        // correct size (communicated in advance, if necessary).
        using ValueType = typename TDatabaseAccess::ValueType;
        auto& r_destination = rAccess.GetValue(ContainerIterator);
        MPIInternals::SendTools<ValueType>::ReadBuffer(pBuffer, r_destination);

        return MPIInternals::BufferAllocation<TDatabaseAccess>::GetSendSize(r_destination);
    }

    template<class TDatabaseAccess>
    std::size_t ReduceValues(
        const typename TDatabaseAccess::SendType* pBuffer,
        TDatabaseAccess& rAccess,
        typename TDatabaseAccess::IteratorType ContainerIterator,
        Operation<OperationType::SumValues>)
    {
        using ValueType = typename TDatabaseAccess::ValueType;
        ValueType& r_current = rAccess.GetValue(ContainerIterator);
        ValueType recv_value(r_current); // creating by copy to have the correct size in dynamic types
        MPIInternals::SendTools<ValueType>::ReadBuffer(pBuffer, recv_value);
        r_current += recv_value;

        return MPIInternals::BufferAllocation<TDatabaseAccess>::GetSendSize(recv_value);
    }

    template<class TDatabaseAccess>
    std::size_t ReduceValues(
        const typename TDatabaseAccess::SendType* pBuffer,
        TDatabaseAccess& rAccess,
        typename TDatabaseAccess::IteratorType ContainerIterator,
        Operation<OperationType::MinValues>)
    {
        using ValueType = typename TDatabaseAccess::ValueType;
        ValueType& r_current = rAccess.GetValue(ContainerIterator);
        ValueType recv_value(r_current); // creating by copy to have the correct size in dynamic types
        MPIInternals::SendTools<ValueType>::ReadBuffer(pBuffer, recv_value);
        if (recv_value < r_current) r_current = recv_value;

        return MPIInternals::BufferAllocation<TDatabaseAccess>::GetSendSize(recv_value);
    }

    template<class TDatabaseAccess>
    std::size_t ReduceValues(
        const typename TDatabaseAccess::SendType* pBuffer,
        TDatabaseAccess& rAccess,
        typename TDatabaseAccess::IteratorType ContainerIterator,
        Operation<OperationType::AndAccessedFlags>)
    {
        using ValueType = typename TDatabaseAccess::ValueType;
        ValueType recv_value;
        MPIInternals::SendTools<ValueType>::ReadBuffer(pBuffer, recv_value);
        rAccess.GetValue(ContainerIterator) &= recv_value | ~rAccess.mrMask;

        return MPIInternals::BufferAllocation<TDatabaseAccess>::GetSendSize(recv_value);
    }

    template<class TDatabaseAccess>
    std::size_t ReduceValues(
        const typename TDatabaseAccess::SendType* pBuffer,
        TDatabaseAccess& rAccess,
        typename TDatabaseAccess::IteratorType ContainerIterator,
        Operation<OperationType::OrAccessedFlags>)
    {
        using ValueType = typename TDatabaseAccess::ValueType;
        ValueType recv_value;
        MPIInternals::SendTools<ValueType>::ReadBuffer(pBuffer, recv_value);
        rAccess.GetValue(ContainerIterator) |= recv_value & rAccess.mrMask;

        return MPIInternals::BufferAllocation<TDatabaseAccess>::GetSendSize(recv_value);
    }

    template<class TDatabaseAccess>
    std::size_t ReduceValues(
        const typename TDatabaseAccess::SendType* pBuffer,
        TDatabaseAccess& rAccess,
        typename TDatabaseAccess::IteratorType ContainerIterator,
        Operation<OperationType::ReplaceAccessedFlags>)
    {
        using ValueType = typename TDatabaseAccess::ValueType;
        ValueType recv_value;
        MPIInternals::SendTools<ValueType>::ReadBuffer(pBuffer, recv_value);

        ValueType& r_current = rAccess.GetValue(ContainerIterator);
        r_current.AssignFlags( (recv_value & rAccess.mrMask) | (r_current & ~rAccess.mrMask) );

        return MPIInternals::BufferAllocation<TDatabaseAccess>::GetSendSize(recv_value);
    }

    template<
        typename TSourceAccess,
        typename TDestinationAccess,
        class TDatabaseAccess,
        typename TReductionOperation>
    void TransferDistributedValues(
        TSourceAccess SourceType,
        TDestinationAccess DestinationType,
        TDatabaseAccess& rAccess,
        TReductionOperation Reduction)
    {
        using DataType = typename TDatabaseAccess::ValueType;
        using BufferType = typename MPIInternals::SendTraits<DataType>::BufferType;
        int destination = 0;

        NeighbourIndicesContainerType& neighbour_indices = NeighbourIndices();

        BufferType send_values;
        BufferType recv_values;

        for (unsigned int i_color = 0; i_color < neighbour_indices.size(); i_color++)
        {
            if ( (destination = neighbour_indices[i_color]) >= 0)
            {
                MeshType& r_source_mesh = GetMesh(i_color, SourceType);
                AllocateBuffer(send_values, r_source_mesh, rAccess);

                MeshType& r_destination_mesh = GetMesh(i_color, DestinationType);
                AllocateBuffer(recv_values, r_destination_mesh, rAccess);

                if ( (send_values.size() == 0) && (recv_values.size() == 0) )
                {
                    continue; // nothing to transfer, skip communication step
                }

                FillBuffer(send_values, r_source_mesh, rAccess);

                mrDataCommunicator.SendRecv(
                    send_values, destination, i_color,
                    recv_values, destination, i_color);

                UpdateValues(recv_values, r_destination_mesh, rAccess, Reduction);
            }
        }
    }

    template<
        typename TSourceAccess,
        typename TDestinationAccess,
        class TDatabaseAccess,
        typename TReductionOperation>
    void TransferDistributedValuesUnknownSize(
        TSourceAccess SourceType,
        TDestinationAccess DestinationType,
        TDatabaseAccess& rAccess,
        TReductionOperation Reduction)
    {
        using DataType = typename TDatabaseAccess::ValueType;
        using BufferType = typename MPIInternals::SendTraits<DataType>::BufferType;
        int destination = 0;

        NeighbourIndicesContainerType& neighbour_indices = NeighbourIndices();

        BufferType send_values;
        BufferType recv_values;

        for (unsigned int i_color = 0; i_color < neighbour_indices.size(); i_color++)
        {
            if ( (destination = neighbour_indices[i_color]) >= 0)
            {
                MeshType& r_source_mesh = GetMesh(i_color, SourceType);
                MeshType& r_destination_mesh = GetMesh(i_color, DestinationType);

                FillBuffer(send_values, r_source_mesh, rAccess);

                std::vector<int> send_size{(int)send_values.size()};
                std::vector<int> recv_size{0};

                mrDataCommunicator.SendRecv(
                    send_size, destination, i_color,
                    recv_size, destination, i_color);

                recv_values.resize(recv_size[0]);

                if ( (send_values.size() == 0) && (recv_values.size() == 0) )
                {
                    continue; // nothing to transfer, skip communication step
                }

                mrDataCommunicator.SendRecv(
                    send_values, destination, i_color,
                    recv_values, destination, i_color);

                UpdateValues(recv_values, r_destination_mesh, rAccess, Reduction);
            }
        }
    }

    template<
        class TDatabaseAccess,
        typename TValue = typename TDatabaseAccess::ValueType,
        typename TSendType = typename MPIInternals::SendTraits<TValue>::SendType>
    void AllocateBuffer(std::vector<TSendType>& rBuffer, const MeshType& rSourceMesh, TDatabaseAccess& rAccess)
    {
        const std::size_t buffer_size = MPIInternals::BufferAllocation<TDatabaseAccess>::GetSendSize(rAccess, rSourceMesh);

        if (rBuffer.size() != buffer_size)
        {
            rBuffer.resize(buffer_size);
        }
    }

    template<
        class TDatabaseAccess,
        typename TValue = typename TDatabaseAccess::ValueType,
        typename TSendType = typename MPIInternals::SendTraits<TValue>::SendType>
    void FillBuffer(std::vector<TSendType>& rBuffer, MeshType& rSourceMesh, TDatabaseAccess& rAccess)
    {
        auto& r_container = rAccess.GetContainer(rSourceMesh);
        TSendType* p_buffer = rBuffer.data();
        std::size_t position = 0;
        for (auto iter = r_container.begin(); iter != r_container.end(); ++iter)
        {
            TValue& r_value = rAccess.GetValue(iter);
            MPIInternals::SendTools<TValue>::WriteBuffer(r_value, p_buffer + position);
            position += MPIInternals::BufferAllocation<TDatabaseAccess>::GetSendSize(r_value);
        }
    }

    template<
        class TDatabaseAccess,
        typename TValue = typename TDatabaseAccess::ValueType>
    void FillBuffer(std::string& rBuffer, MeshType& rSourceMesh, TDatabaseAccess& rAccess)
    {
        StreamSerializer serializer;
        auto& r_container = rAccess.GetContainer(rSourceMesh);

        for (auto iter = r_container.begin(); iter != r_container.end(); ++iter)
        {
            TValue& r_value = rAccess.GetValue(iter);
            serializer.save("Value", r_value);
        }

        rBuffer = serializer.GetStringRepresentation();
    }

    template<
        class TDatabaseAccess,
        typename TReductionOperation,
        typename TValue = typename TDatabaseAccess::ValueType,
        typename TSendType = typename MPIInternals::SendTraits<TValue>::SendType>
    void UpdateValues(
        const std::vector<TSendType>& rBuffer,
        MeshType& rSourceMesh,
        TDatabaseAccess& rAccess,
        TReductionOperation Operation)
    {
        auto& r_container = rAccess.GetContainer(rSourceMesh);
        const TSendType* p_buffer = rBuffer.data();
        std::size_t position = 0;

        for (auto iter = r_container.begin(); iter != r_container.end(); ++iter)
        {
            position += ReduceValues(p_buffer + position, rAccess, iter, Operation);
        }

        KRATOS_WARNING_IF_ALL_RANKS("MPICommunicator", position > rBuffer.size())
        << GetDataCommunicator()
        << "Error in estimating receive buffer size." << std::endl;
    }

    template<
        class TDatabaseAccess,
        typename TValue = typename TDatabaseAccess::ValueType>
    void UpdateValues(
        const std::string& rBuffer,
        MeshType& rSourceMesh,
        TDatabaseAccess& rAccess,
        Operation<OperationType::Replace>)
    {
        StreamSerializer serializer;
        std::stringstream* serializer_buffer = (std::stringstream *)serializer.pGetBuffer();
        serializer_buffer->write(rBuffer.data(), rBuffer.size());

        auto& r_container = rAccess.GetContainer(rSourceMesh);
        for (auto iter = r_container.begin(); iter != r_container.end(); ++iter)
        {
            serializer.load("Value", rAccess.GetValue(iter));
        }
    }

    template<
        typename TSourceAccess,
        typename TDestinationAccess,
        class TDatabaseAccess>
    void MatchDynamicVectorSizes(
        TSourceAccess SourceType,
        TDestinationAccess DestinationType,
        TDatabaseAccess& rAccess)
    {
        using TVectorType = typename TDatabaseAccess::ValueType;
        int destination = 0;

        NeighbourIndicesContainerType& neighbour_indices = NeighbourIndices();

        std::vector<int> send_sizes;
        std::vector<int> recv_sizes;

        bool resize_error = false;
        std::stringstream error_detail;

        for (unsigned int i_color = 0; i_color < neighbour_indices.size(); i_color++)
        {
            if ( (destination = neighbour_indices[i_color]) >= 0)
            {
                MeshType& r_source_mesh = GetMesh(i_color, SourceType);
                const auto& r_source_container = rAccess.GetContainer(r_source_mesh);
                const std::size_t num_values_to_send = r_source_container.size();

                MeshType& r_destination_mesh = GetMesh(i_color, DestinationType);
                auto& r_destination_container = rAccess.GetContainer(r_destination_mesh);
                const std::size_t num_values_to_recv = r_destination_container.size();

                if ( (num_values_to_send == 0) && (num_values_to_recv == 0) )
                {
                    continue; // nothing to transfer, skip communication step
                }

                if (send_sizes.size() != num_values_to_send)
                {
                    send_sizes.resize(num_values_to_send);
                }

                if (recv_sizes.size() != num_values_to_recv)
                {
                    recv_sizes.resize(num_values_to_recv);
                }

                int position = 0;
                for (auto iter = r_source_container.begin(); iter != r_source_container.end(); ++iter)
                {
                    const TVectorType& r_value = rAccess.GetValue(iter);
                    send_sizes[position++] = r_value.size();
                }

                mrDataCommunicator.SendRecv(
                    send_sizes, destination, i_color,
                    recv_sizes, destination, i_color);

                position = 0;
                for (auto iter = r_destination_container.begin(); iter != r_destination_container.end(); ++iter)
                {
                    std::size_t source_size = recv_sizes[position++];
                    if (source_size != 0)
                    {
                        TVectorType& r_value = rAccess.GetValue(iter);
                        if (r_value.size() == source_size)
                        {
                            continue; // everything ok!
                        }
                        else if (r_value.size() == 0)
                        {
                            r_value.resize(source_size, false);
                        }
                        else
                        {
                            resize_error = true;
                            error_detail
                            << "On rank " << mrDataCommunicator.Rank() << ": "
                            << "local size: " << r_value.size() << " "
                            << "source size: " << source_size << "." << std::endl;
                        }
                    }
                }
            }
        }

        KRATOS_ERROR_IF(mrDataCommunicator.ErrorIfTrueOnAnyRank(resize_error))
        << "Size mismatch in Vector size synchronization." << std::endl
        << error_detail.str();
    }

    template<
        typename TSourceAccess,
        typename TDestinationAccess,
        class TDatabaseAccess>
    void MatchDynamicMatrixSizes(
        TSourceAccess SourceType,
        TDestinationAccess DestinationType,
        TDatabaseAccess& rAccess)
    {
        using TMatrixType = typename TDatabaseAccess::ValueType;
        int destination = 0;

        NeighbourIndicesContainerType& neighbour_indices = NeighbourIndices();

        std::vector<int> send_sizes;
        std::vector<int> recv_sizes;

        bool resize_error = false;
        std::stringstream error_detail;

        for (unsigned int i_color = 0; i_color < neighbour_indices.size(); i_color++)
        {
            if ( (destination = neighbour_indices[i_color]) >= 0)
            {
                MeshType& r_source_mesh = GetMesh(i_color, SourceType);
                const auto& r_source_container = rAccess.GetContainer(r_source_mesh);
                const std::size_t num_values_to_send = 2*r_source_container.size();

                MeshType& r_destination_mesh = GetMesh(i_color, DestinationType);
                auto& r_destination_container = rAccess.GetContainer(r_destination_mesh);
                const std::size_t num_values_to_recv = 2*r_destination_container.size();

                if ( (num_values_to_send == 0) && (num_values_to_recv == 0) )
                {
                    continue; // nothing to transfer, skip communication step
                }

                if (send_sizes.size() != num_values_to_send)
                {
                    send_sizes.resize(num_values_to_send);
                }

                if (recv_sizes.size() != num_values_to_recv)
                {
                    recv_sizes.resize(num_values_to_recv);
                }

                int position = 0;
                for (auto iter = r_source_container.begin(); iter != r_source_container.end(); ++iter)
                {
                    const TMatrixType& r_value = rAccess.GetValue(iter);
                    send_sizes[position++] = r_value.size1();
                    send_sizes[position++] = r_value.size2();
                }

                mrDataCommunicator.SendRecv(
                    send_sizes, destination, i_color,
                    recv_sizes, destination, i_color);

                position = 0;
                for (auto iter = r_destination_container.begin(); iter != r_destination_container.end(); ++iter)
                {
                    std::size_t source_size_1 = recv_sizes[position++];
                    std::size_t source_size_2 = recv_sizes[position++];
                    if (source_size_1 != 0  && source_size_2 != 0)
                    {
                        TMatrixType& r_value = rAccess.GetValue(iter);
                        if (r_value.size1() == source_size_1 && r_value.size2() == source_size_2)
                        {
                            continue; // everything ok!
                        }
                        else if (r_value.size1() == 0 && r_value.size2() == 0)
                        {
                            r_value.resize(source_size_1, source_size_2, false);
                        }
                        else
                        {
                            resize_error = true;
                            error_detail
                            << "On rank " << mrDataCommunicator.Rank() << ": "
                            << "local size: (" << r_value.size1() << "," << r_value.size2() << ") "
                            << "source size: (" << source_size_1 << "," << source_size_2 << ")." << std::endl;
                        }
                    }
                }
            }
        }

        KRATOS_ERROR_IF(mrDataCommunicator.ErrorIfTrueOnAnyRank(resize_error))
        << "Size mismatch in Matrix size synchronization." << std::endl
        << error_detail.str();
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

}; // Class MPICommunicator

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  MPICommunicator& rThis);

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const MPICommunicator& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_MPI_COMMUNICATOR_H_INCLUDED  defined
