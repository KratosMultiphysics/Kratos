//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/lock_object.h"
#include "geometries/point.h"
#include "includes/dof.h"
#include "containers/pointer_vector_set.h"
#include "containers/variables_list_data_value_container.h"
#include "containers/flags.h"
#include "intrusive_ptr/intrusive_ptr.hpp"
#include "containers/global_pointers_vector.h"
#include "containers/data_value_container.h"
#include "containers/nodal_data.h"
#include "includes/kratos_flags.h"

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

/**
 * @ingroup KratosCore
 * @class Node
 * @brief This class defines the node
 * @details The node class from Kratos is defined in this class. Represents a node in a finite element model.
 * The Node class provides various functionalities for handling nodal data,
 * degrees of freedom (DOFs), etc...
 * @author Pooyan Dadvand
 * @author Riccardo Rossi
 */
class Node final
    : public Point, public Flags
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Node
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Node);

    /// Base type
    using BaseType = Point;

    /// Point type definition
    using PointType = Point;

    /// Dof type
    using DofType = Dof<double>;

    /// The index type
    using IndexType = std::size_t;

    /// The size type
    using SizeType = std::size_t;

    /// The DoF container type definition
    using DofsContainerType = std::vector<std::unique_ptr<Dof<double>>>;

    /// The solution step data container type
    using SolutionStepsNodalDataContainerType = VariablesListDataValueContainer;

    /// The block type
    using BlockType = VariablesListDataValueContainer::BlockType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief @brief Default constructor.
     */
    Node()
        : BaseType()
        , Flags()
        , mNodalData(0)
        , mDofs()
        , mData()
        , mInitialPosition()
        , mNodeLock()
    {
        CreateSolutionStepData();
    }

    /**
     * @brief Constructor with a given node ID.
     * @param NewId The unique index identifier of the node.
     */
    explicit Node(IndexType NewId)
        : BaseType()
        , Flags()
        , mNodalData(NewId)
        , mDofs()
        , mData()
        , mInitialPosition()
        , mNodeLock()
    {
        KRATOS_ERROR <<  "Calling the default constructor for the node ... illegal operation!!" << std::endl;
        CreateSolutionStepData();
    }

    /**
     * @brief Constructor for a 1D node.
     * @param NewId The unique index identifier of the node.
     * @param NewX The X-coordinate of the node.
     */
    Node(IndexType NewId, const double NewX)
        : BaseType(NewX)
        , Flags()
        , mNodalData(NewId)
        , mDofs()
        , mData()
        , mInitialPosition(NewX)
        , mNodeLock()
    {
        CreateSolutionStepData();
    }

    /**
     * @brief Constructor for a 2D node.
     * @param NewId The unique index identifier of the node.
     * @param NewX The X-coordinate of the node.
     * @param NewY The Y-coordinate of the node.
     */
    Node(IndexType NewId, const double NewX, const double NewY)
        : BaseType(NewX, NewY)
        , Flags()
        , mNodalData(NewId)
        , mDofs()
        , mData()
        , mInitialPosition(NewX, NewY)
        , mNodeLock()
    {
        CreateSolutionStepData();
    }

    /**
     * @brief Constructor for a 3D node.
     * @param NewId The unique index identifier of the node.
     * @param NewX The X-coordinate of the node.
     * @param NewY The Y-coordinate of the node.
     * @param NewZ The Z-coordinate of the node.
     */
    Node(IndexType NewId, const double NewX, const double NewY, const double NewZ)
        : BaseType(NewX, NewY, NewZ)
        , Flags()
        , mNodalData(NewId)
        , mDofs()
        , mData()
        , mInitialPosition(NewX, NewY, NewZ)
        , mNodeLock()
    {
        CreateSolutionStepData();
    }

    /**
     * @brief Constructor from an existing point.
     * @param NewId The unique index identifier of the node.
     * @param rThisPoint The point from which to initialize the node.
     */
    Node(IndexType NewId, Point const& rThisPoint)
        : BaseType(rThisPoint)
        , Flags()
        , mNodalData(NewId)
        , mDofs()
        , mData()
        , mInitialPosition(rThisPoint)
        , mNodeLock()
    {
        CreateSolutionStepData();
    }

    /** Copy constructor. Initialize this node with given node.*/
    explicit Node(const Node& rOtherNode)
        : BaseType(rOtherNode)
        , Flags(rOtherNode)
        , mNodalData(rOtherNode.Id())
        , mData(rOtherNode.mData)
        , mInitialPosition(rOtherNode.mInitialPosition)
        , mNodeLock()
    {
        // Deep copying the nodal data
        this->mNodalData = rOtherNode.mNodalData;

        // Deep copying the dofs
        for (const auto& it_dof : rOtherNode.mDofs) {
            mDofs.push_back(Kratos::make_unique<DofType>(*it_dof));
            mDofs.back()->SetNodalData(&mNodalData);
        }
    }

    /**
     * @brief Constructor using coordinates stored in a vector expression.
     * @details Initializes the node with the coordinates in the given vector.
     * @tparam TVectorType The type of the vector expression.
     * @param NewId The unique index identifier of the node.
     * @param rOtherCoordinates The vector expression containing the node coordinates.
     */
    template<class TVectorType>
    Node(IndexType NewId, vector_expression<TVectorType> const& rOtherCoordinates)
        : BaseType(rOtherCoordinates)
        , Flags()
        , mNodalData(NewId)
        , mDofs()
        , mData()
        , mInitialPosition(rOtherCoordinates)
        , mNodeLock()
    {
        CreateSolutionStepData();
    }

    /**
     * @brief Constructor using coordinates stored in a standard vector.
     * @param NewId The unique index identifier of the node.
     * @param rOtherCoordinates The std::vector containing the node coordinates.
     */
    Node(IndexType NewId, std::vector<double> const& rOtherCoordinates)
        : BaseType(rOtherCoordinates)
        , Flags()
        , mNodalData(NewId)
        , mDofs()
        , mData()
        , mInitialPosition()
        , mNodeLock()
    {
        CreateSolutionStepData();
    }

    /**
     * @brief Constructor for a 3D node with variables list and data.
     * @param NewId The unique index identifier of the node.
     * @param NewX The X-coordinate of the node.
     * @param NewY The Y-coordinate of the node.
     * @param NewZ The Z-coordinate of the node.
     * @param pVariablesList Pointer to the variables list associated with the node.
     * @param ThisData Pointer to the data block associated with the node.
     * @param NewQueueSize The queue size for storing historical data (default = 1).
     */
    Node(IndexType NewId, const double NewX, const double NewY, const double NewZ, 
         VariablesList::Pointer pVariablesList, BlockType const * ThisData, SizeType NewQueueSize = 1)
         : BaseType(NewX, NewY, NewZ)
         , Flags()
         , mNodalData(NewId, pVariablesList,ThisData,NewQueueSize)
         , mDofs()
         , mData()
         , mInitialPosition(NewX, NewY, NewZ)
         , mNodeLock()
     {
     }

    /**
     * @brief Creates a clone of the current node with a new ID.
     * @details This method creates a new node that is a copy of the current node,
     * but with a different ID specified by the user.
     * @param NewId The ID to be assigned to the new node.
     * @return Node::Pointer A pointer to the newly created node.
     */
    typename Node::Pointer Clone(IndexType NewId)
    {
        Node::Pointer p_new_node = Kratos::make_intrusive<Node>(*this);
        p_new_node->SetId(NewId);
        return p_new_node;
    }

    /**
     * @brief Destructor for the Node class.
     * @details Clears any stored solution step data before destruction.
     */
    ~Node() override
    {
        ClearSolutionStepsData();
    }

    /**
     * @brief Gets the reference count of the node.
     * @details Public API of intrusive_ptr
     * @return The current reference count.
     */
    unsigned int use_count() const noexcept
    {
        return mReferenceCounter;
    }

    /**
     * @brief Retrieves the node ID.
     * @return The ID of the node.
     */
    IndexType Id() const
    {
        return mNodalData.Id();
    }

    /**
     * @brief Retrieves the node ID (alternative method).
     * @return The ID of the node.
     */
    IndexType GetId() const
    {
        return mNodalData.Id();
    }

    /**
     * @brief Sets a new ID for the node.
     * @param NewId The new ID to assign.
     */
    void SetId(IndexType NewId)
    {
        mNodalData.SetId(NewId);
    }

    /**
     * @brief Retrieves the lock object for thread safety.
     * @return A reference to the lock object.
     */
    LockObject& GetLock()
    {
        return mNodeLock;
    }

    /**
     * @brief Locks the node for thread-safe operations.
     */
    inline void SetLock()
    {
        mNodeLock.lock();
    }

    /**
     * @brief Unlocks the node.
     */
    inline void UnSetLock()
    {
        mNodeLock.unlock();
    }

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Assignment operator for deep copying a node.
     * @details Copies the nodal data, DOFs, and initial position from another node.
     * @param rOther The node to copy from.
     * @return A reference to the copied node.
     */
    Node& operator=(const Node& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator =(rOther);
    
        mNodalData = rOther.mNodalData;
    
        // Deep copying the dofs
        for(typename DofsContainerType::const_iterator it_dof = rOther.mDofs.begin() ; it_dof != rOther.mDofs.end() ; it_dof++) {
            pAddDof(**it_dof);
        }
    
        mData = rOther.mData;
        mInitialPosition = rOther.mInitialPosition;
    
        return *this;
    }

    /**
     * @brief Equality operator.
     * @details Compares two nodes based on their position.
     * @param rOther The node to compare with.
     * @return True if nodes are equal, false otherwise.
     */
    bool operator==(const Node& rOther)
    {
        return Point::operator ==(rOther);
    }    

    /**
     * @brief Accesses a solution step value for a given variable and step index.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable to access.
     * @param SolutionStepIndex The solution step index.
     * @return Reference to the stored value.
     */
    template<class TVariableType>
    typename TVariableType::Type& operator()(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
    {
        return GetSolutionStepValue(rThisVariable, SolutionStepIndex);
    }

    /**
     * @brief Accesses the current solution step value for a given variable.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable to access.
     * @return Reference to the stored value.
     */
    template<class TVariableType>
    typename TVariableType::Type& operator()(const TVariableType& rThisVariable)
    {
        return GetSolutionStepValue(rThisVariable);
    }

    /**
     * @brief Accessor for node data using a variable.
     * @tparam TDataType The type of data stored.
     * @param rThisVariable The variable representing the data.
     * @return Reference to the stored data.
     */
    template<class TDataType>
    TDataType& operator[](const Variable<TDataType>& rThisVariable)
    {
        return GetValue(rThisVariable);
    }

    /**
     * @brief Const accessor for node data using a variable.
     * @tparam TDataType The type of data stored.
     * @param rThisVariable The variable representing the data.
     * @return Const reference to the stored data.
     */
    template<class TDataType>
    const TDataType& operator[](const Variable<TDataType>& rThisVariable) const
    {
        return GetValue(rThisVariable);
    }

    /**
     * @brief Accessor for node coordinates.
     * @param ThisIndex The index of the coordinate (0 = X, 1 = Y, 2 = Z).
     * @return Reference to the coordinate value.
     */
    double& operator[](IndexType ThisIndex)
    {
        return BaseType::operator[](ThisIndex);
    }    

    /**
     * @brief Const accessor for node coordinates.
     * @param ThisIndex The index of the coordinate (0 = X, 1 = Y, 2 = Z).
     * @return The coordinate value.
     */
    double operator[](IndexType ThisIndex) const
    {
        return BaseType::operator[](ThisIndex);
    }

    ///@}
    ///@name Nodal Data
    ///@{

    /**
     * @brief Creates a new solution step data entry.
     * @details Adds a new solution step at the front of the solution step data container.
     */
    void CreateSolutionStepData()
    {
        SolutionStepData().PushFront();
    }

    /**
     * @brief Clones the current solution step data.
     * @details Copies the front solution step data to create a new step.
     */
    void CloneSolutionStepData()
    {
        SolutionStepData().CloneFront();
    }

    /**
     * @brief Overwrites solution step data at a given index.
     * @details Copies the data from a source solution step index to a destination index.
     * @param SourceSolutionStepIndex The index of the source solution step.
     * @param DestinationSourceSolutionStepIndex The index of the destination step where the data will be assigned.
     */
    void OverwriteSolutionStepData(IndexType SourceSolutionStepIndex, IndexType DestinationSourceSolutionStepIndex)
    {
        SolutionStepData().AssignData(SolutionStepData().Data(SourceSolutionStepIndex), DestinationSourceSolutionStepIndex);
    }

    /**
     * @brief Clears all stored solution step data.
     * @details Removes all solution step data entries.
     */
    void ClearSolutionStepsData()
    {
        SolutionStepData().Clear();
    }

    /**
     * @brief Sets the list of solution step variables.
     * @details Assigns a new variables list to the solution step data.
     * @param pVariablesList Pointer to the new list of variables.
     */
    void SetSolutionStepVariablesList(VariablesList::Pointer pVariablesList)
    {
        SolutionStepData().SetVariablesList(pVariablesList);
    }

    /**
     * @brief Retrieves the solution step data container.
     * @return Reference to the solution step data container.
     */
    VariablesListDataValueContainer& SolutionStepData()
    {
        return mNodalData.GetSolutionStepData();
    }

    /**
     * @brief Retrieves the solution step data container (const version).
     * @return Const reference to the solution step data container.
     */
    const VariablesListDataValueContainer& SolutionStepData() const
    {
        return mNodalData.GetSolutionStepData();
    }

    /**
     * @deprecated This method is deprecated. Use GetData() instead.
     * @brief Retrieves the data container for the node.
     * @return Reference to the data container.
     */
    KRATOS_DEPRECATED_MESSAGE("This method is deprecated. Use 'GetData()' instead.")
    DataValueContainer& Data()
    {
        return mData;
    }

    /**
     * @brief Retrieves the data container for the node.
     * @return Reference to the data container.
     */
    DataValueContainer& GetData()
    {
        return mData;
    }

    /**
     * @brief Retrieves the data container for the node (const version).
     * @return Const reference to the data container.
     */
    const DataValueContainer& GetData() const
    {
        return mData;
    }

    /**
     * @brief Retrieves the solution step value for a given variable.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable whose value is to be retrieved.
     * @return Reference to the stored value.
     */
    template<class TVariableType>
    typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable)
    {
        return SolutionStepData().GetValue(rThisVariable);
    }

    /**
     * @brief Retrieves the solution step value for a given variable (const version).
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable whose value is to be retrieved.
     * @return Const reference to the stored value.
     */
    template<class TVariableType>
    typename TVariableType::Type const& GetSolutionStepValue(const TVariableType& rThisVariable) const
    {
        return SolutionStepData().GetValue(rThisVariable);
    }

    /**
     * @brief Retrieves the solution step value for a given variable at a specific step.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable whose value is to be retrieved.
     * @param SolutionStepIndex The index of the solution step.
     * @return Reference to the stored value.
     */
    template<class TVariableType>
    typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
    {
        return SolutionStepData().GetValue(rThisVariable, SolutionStepIndex);
    }

    /**
     * @brief Retrieves the solution step value for a given variable at a specific step (const version).
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable whose value is to be retrieved.
     * @param SolutionStepIndex The index of the solution step.
     * @return Const reference to the stored value.
     */
    template<class TVariableType>
    typename TVariableType::Type const& GetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex) const
    {
        return SolutionStepData().GetValue(rThisVariable, SolutionStepIndex);
    }

    /**
     * @brief Checks if the solution step data contains a given variable.
     * @param rThisVariable The variable to check.
     * @return True if the variable exists in the solution step data, false otherwise.
     */
    bool SolutionStepsDataHas(const VariableData& rThisVariable) const
    {
        return SolutionStepData().Has(rThisVariable);
    }

    /**
     * @brief Retrieves the solution step value for a given variable quickly.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable whose value is to be retrieved.
     * @return Reference to the stored value.
     */
    template<class TVariableType>
    typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable)
    {
        return SolutionStepData().FastGetValue(rThisVariable);
    }

    /**
     * @brief Retrieves the solution step value for a given variable quickly (const version).
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable whose value is to be retrieved.
     * @return Const reference to the stored value.
     */
    template<class TVariableType>
    const typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable) const
    {
        return SolutionStepData().FastGetValue(rThisVariable);
    }

    /**
     * @brief Retrieves the solution step value for a given variable at a specific step quickly.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable whose value is to be retrieved.
     * @param SolutionStepIndex The index of the solution step.
     * @return Reference to the stored value.
     */
    template<class TVariableType>
    typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
    {
        return SolutionStepData().FastGetValue(rThisVariable, SolutionStepIndex);
    }

    /**
     * @brief Retrieves the solution step value for a given variable at a specific step quickly (const version).
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable whose value is to be retrieved.
     * @param SolutionStepIndex The index of the solution step.
     * @return Const reference to the stored value.
     */
    template<class TVariableType>
    const typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex) const
    {
        return SolutionStepData().FastGetValue(rThisVariable, SolutionStepIndex);
    }

    /**
     * @brief Retrieves the solution step value for a given variable at a specific step and position quickly.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable whose value is to be retrieved.
     * @param SolutionStepIndex The index of the solution step.
     * @param ThisPosition The position in the solution step data.
     * @return Reference to the stored value.
     */
    template<class TVariableType>
    typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex, IndexType ThisPosition)
    {
        return SolutionStepData().FastGetValue(rThisVariable, SolutionStepIndex, ThisPosition);
    }

    /**
     * @brief Retrieves the current solution step value for a given variable at a specific position quickly.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable whose value is to be retrieved.
     * @param ThisPosition The position in the solution step data.
     * @return Reference to the stored value.
     */
    template<class TVariableType>
    typename TVariableType::Type& FastGetCurrentSolutionStepValue(const TVariableType& rThisVariable, IndexType ThisPosition)
    {
        return SolutionStepData().FastGetCurrentValue(rThisVariable, ThisPosition);
    }

    /**
     * @brief Retrieves the stored value of a given variable.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable whose value is to be retrieved.
     * @return Reference to the stored value.
     */
    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rThisVariable)
    {
        return mData.GetValue(rThisVariable);
    }

    /**
     * @brief Retrieves the stored value of a given variable (const version).
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable whose value is to be retrieved.
     * @return Const reference to the stored value.
     */
    template<class TVariableType>
    typename TVariableType::Type const& GetValue(const TVariableType& rThisVariable) const
    {
        return mData.GetValue(rThisVariable);
    }

    /**
     * @brief Retrieves the value of a variable, checking both node data and solution step data.
     * @details If the variable is not present in the node's data, it fetches from the solution step data.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable to retrieve.
     * @param SolutionStepIndex The solution step index.
     * @return Reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
    {
        if(!mData.Has(rThisVariable))
            return SolutionStepData().GetValue(rThisVariable, SolutionStepIndex);

        return mData.GetValue(rThisVariable);
    }

    /**
     * @brief Retrieves the value of a variable, checking both node data and solution step data (const version).
     * @details If the variable is not present in the node's data, it fetches from the solution step data.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable to retrieve.
     * @param SolutionStepIndex The solution step index.
     * @return Const reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type const& GetValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex) const
    {
        if(!mData.Has(rThisVariable))
            return SolutionStepData().GetValue(rThisVariable, SolutionStepIndex);

        return mData.GetValue(rThisVariable);
    }

    /**
     * @brief Sets the value of a given variable.
     * @details Stores a new value for the given variable in the node's data.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable to store.
     * @param rValue The value to assign.
     */
    template<class TVariableType>
    void SetValue(const TVariableType& rThisVariable, typename TVariableType::Type const& rValue)
    {
        mData.SetValue(rThisVariable, rValue);
    }

    /**
     * @brief Checks if a variable exists in the node's data.
     * @tparam TDataType The type of the variable.
     * @param rThisVariable The variable to check.
     * @return True if the variable exists, false otherwise.
     */
    template<class TDataType>
    bool Has(const Variable<TDataType>& rThisVariable) const
    {
        return mData.Has(rThisVariable);
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Fix a Degree of Freedom (DoF) associated with a given variable.
     * @details This function searches for the DoF corresponding to @p rDofVariable
     * in the internal list of DoFs (@p mDofs) and fixes it if found. If the
     * variable is not found, it adds a new DoF and fixes it.
     * @tparam TVariableType Type of the variable associated with the DoF.
     * @param rDofVariable Reference to the variable to fix.
     */
    template<class TVariableType>
    inline void Fix(const TVariableType& rDofVariable)
    {
        for(auto it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == rDofVariable){
                (*it_dof)->FixDof();
                return;
            }
        }

#ifdef KRATOS_DEBUG
        if(OpenMPUtils::IsInParallel() != 0) {
            KRATOS_ERROR << "Attempting to Fix the variable: " << rDofVariable << " within a parallel region. This is not permitted. Create the Dof first by pAddDof" << std::endl;
        }
#endif
        pAddDof(rDofVariable)->FixDof();
    }

    /**
     * @brief Free a Degree of Freedom (DoF) associated with a given variable.
     * @details This function searches for the DoF corresponding to @p rDofVariable
     * in the internal list of DoFs (@p mDofs) and frees it if found. If the
     * variable is not found, it adds a new DoF and frees it.
     * @tparam TVariableType Type of the variable associated with the DoF.
     * @param rDofVariable Reference to the variable to free.
     */
    template<class TVariableType>
    inline void Free(const TVariableType& rDofVariable)
    {
         for(auto it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == rDofVariable){
                (*it_dof)->FreeDof();
                return;
            }
        }

#ifdef KRATOS_DEBUG
        if(OpenMPUtils::IsInParallel() != 0) {
            KRATOS_ERROR << "Attempting to Free the variable: " << rDofVariable << " within a parallel region. This is not permitted. Create the Dof first by pAddDof" << std::endl;
        }
#endif
        pAddDof(rDofVariable)->FreeDof();
    }

    /**
     * @brief Get the buffer size.
     * @details This function returns the size of the buffer used for storing
     * solution step data.
     * @return IndexType The size of the buffer.
     */
    IndexType GetBufferSize() const
    {
        return SolutionStepData().QueueSize();
    }

    /**
     * @brief Set the buffer size.
     * @details This function resizes the buffer used for storing solution step data
     * to the given size @p NewBufferSize.
     * @param NewBufferSize New size for the buffer.
     */
    void SetBufferSize(IndexType NewBufferSize)
    {
        SolutionStepData().Resize(NewBufferSize);
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Returns the initial position of the node (const version).
     * @details This function provides read-only access to the initial position of the node.
     * @return A const reference to the initial position of the node.
     */
    const Point& GetInitialPosition() const
    {
        return mInitialPosition;
    }

    /**
     * @brief Returns the initial position of the node (non-const version).
     * @details This function provides access to the initial position of the node and allows modifications.
     * @return A reference to the initial position of the node.
     */
    Point& GetInitialPosition()
    {
        return mInitialPosition;
    }

    /**
     * @brief Returns the X coordinate of the initial position.
     * @details This function provides access to the X coordinate of the initial position and allows modifications.
     * @return A reference to the X coordinate of the initial position.
     */
    double& X0()
    {
        return mInitialPosition.X();
    }

    /**
     * @brief Returns the Y coordinate of the initial position.
     * @details This function provides access to the Y coordinate of the initial position and allows modifications.
     * @return A reference to the Y coordinate of the initial position.
     */
    double& Y0()
    {
        return mInitialPosition.Y();
    }

    /**
     * @brief Returns the Z coordinate of the initial position.
     * @details This function provides access to the Z coordinate of the initial position and allows modifications.
     * @return A reference to the Z coordinate of the initial position.
     */
    double& Z0()
    {
        return mInitialPosition.Z();
    }

    /**
     * @brief Returns the X coordinate of the initial position (const version).
     * @details This function provides read-only access to the X coordinate of the initial position.
     * @return The X coordinate of the initial position.
     */
    double X0() const
    {
        return mInitialPosition.X();
    }

    /**
     * @brief Returns the Y coordinate of the initial position (const version).
     * @details This function provides read-only access to the Y coordinate of the initial position.
     * @return The Y coordinate of the initial position.
     */
    double Y0() const
    {
        return mInitialPosition.Y();
    }

    /**
     * @brief Returns the Z coordinate of the initial position (const version).
     * @details This function provides read-only access to the Z coordinate of the initial position.
     * @return The Z coordinate of the initial position.
     */
    double Z0() const
    {
        return mInitialPosition.Z();
    }

    /**
     * @brief Sets the initial position of the node using a PointType.
     * @details This function sets the initial position of the node by copying the values from the given `NewInitialPosition`.
     * @param NewInitialPosition The new initial position to be set.
     */
    void SetInitialPosition(const Point& NewInitialPosition)
    {
        mInitialPosition.X() = NewInitialPosition.X();
        mInitialPosition.Y() = NewInitialPosition.Y();
        mInitialPosition.Z() = NewInitialPosition.Z();
    }

    /**
     * @brief Sets the initial position of the node using X, Y, and Z coordinates.
     * @details This function sets the initial position of the node by directly assigning values to the X, Y, and Z coordinates.
     * @param X The X coordinate of the new initial position.
     * @param Y The Y coordinate of the new initial position.
     * @param Z The Z coordinate of the new initial position.
     */
    void SetInitialPosition(double X,double Y, double Z)    
    {
        mInitialPosition.X() = X;
        mInitialPosition.Y() = Y;
        mInitialPosition.Z() = Z;
    }

    /**
     * @brief Returns the pointer to the list of variables associated with the node.
     * @details This function provides access to the list of variables for the current solution step.
     * @return A pointer to the list of variables.
     */
    VariablesList::Pointer pGetVariablesList()
    {
        return SolutionStepData().pGetVariablesList();
    }    

    /**
     * @brief Returns the pointer to the list of variables associated with the node (const version).
     * @details This function provides read-only access to the list of variables for the current solution step.
     * @return A const pointer to the list of variables.
     */
    const VariablesList::Pointer pGetVariablesList() const
    {
        return SolutionStepData().pGetVariablesList();
    }

    ///@}
    ///@name Dofs
    ///@{

    /**
     * @brief Get DoF  position with a given position
     * @param rDofVariable Name of the variable to search the position
     * @tparam TVariableType The variable type template argument
     * @return The position of the given DoF variable
     */
    template<class TVariableType>
    inline unsigned int GetDofPosition(TVariableType const& rDofVariable) const
    {
        typename DofsContainerType::const_iterator it_dof = mDofs.end();
        for(it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == rDofVariable){
                break;
            }
        }

        return it_dof - mDofs.begin();
    }

    /**
     * @brief Get dof with a given position. If not found it is search
     * @param rDofVariable Name of the variable
     * @param pos Position of the DoF
     * @tparam TVariableType The variable type template argument
     * @return The DoF associated to the given variable
     */
    template<class TVariableType>
    inline const DofType& GetDof(TVariableType const& rDofVariable, int pos) const
    {
        auto it_begin = mDofs.begin();
        auto it_end = mDofs.end();
        typename DofsContainerType::const_iterator it;
        // If the guess is exact return the guess
        if(pos < it_end-it_begin) {
            it = it_begin + pos;
            if( (*it)->GetVariable() == rDofVariable) {
                return **it;
            }
        }

        // Otherwise do a find
        for(auto it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == rDofVariable){
                return **it_dof;
            }
        }

        KRATOS_ERROR <<  "Non-existent DOF in node #" << Id() << " for variable : " << rDofVariable.Name() << std::endl;
    }

    /**
     * @brief Get DoF for a given variable
     * @param rDofVariable Name of the variable
     * @tparam TVariableType The variable type template argument
     * @return The DoF associated to the given variable
     */
    template<class TVariableType>
    inline const DofType& GetDof(TVariableType const& rDofVariable) const
    {
        for(auto it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == rDofVariable){
                return **it_dof;
            }
        }

        KRATOS_ERROR <<  "Non-existent DOF in node #" << Id() << " for variable : " << rDofVariable.Name() << std::endl;

    }

    /**
     * @brief Returns all the degrees of freedom (DOFs) associated with the node.
     * @details This function provides access to the container holding all the degrees of freedom
     * (DOFs) for the node. The returned container is mutable, allowing modifications to the DOFs.
     * @return A reference to the container of DOFs.
     */
    DofsContainerType& GetDofs()
    {
        return mDofs;
    }

    /**
     * @brief Returns all the degrees of freedom (DOFs) associated with the node (const version).
     * @details This function provides a read-only access to the container holding all the degrees of freedom
     * (DOFs) for the node. The returned container cannot be modified.
     * @return A const reference to the container of DOFs.
     */
    const DofsContainerType& GetDofs() const
    {
        return mDofs;
    }

    /**
     * @brief Get DoF counted pointer for a given variable
     * @param rDofVariable Name of the variable
     * @tparam TVariableType The variable type template argument
     * @return The DoF associated to the given variable
     */
    template<class TVariableType>
    inline typename DofType::Pointer pGetDof(TVariableType const& rDofVariable) const
    {
        for(auto it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == rDofVariable){
                return (*it_dof).get();
            }
        }

        KRATOS_ERROR <<  "Non-existent DOF in node #" << Id() << " for variable : " << rDofVariable.Name() << std::endl;
    }

    /**
     * @brief Get DoF counted pointer with a given position. If not found it is search
     * @param rDofVariable Name of the variable
     * @param Position Position of the DoF
     * @tparam TVariableType The variable type template argument
     * @return The DoF associated to the given variable
     */
    template<class TVariableType>
    inline typename DofType::Pointer pGetDof(
        TVariableType const& rDofVariable,
        int Position
        ) const
    {
        const auto it_begin = mDofs.begin();
        const auto it_end = mDofs.end();
        // If the guess is exact return the guess
        if(Position < it_end-it_begin) {
            auto it_dof = it_begin + Position;
            if( (*it_dof)->GetVariable() == rDofVariable) {
                return (*it_dof).get();
            }
        }

        // Otherwise do a find
        for(auto it_dof = it_begin; it_dof != it_end; ++it_dof){
            if((*it_dof)->GetVariable() == rDofVariable){
                return (*it_dof).get();
            }
        }

        KRATOS_ERROR <<  "Non-existent DOF in node #" << Id() << " for variable : " << rDofVariable.Name() << std::endl;
    }

    /**
     * @brief Adds a degree of freedom (DOF) to the node, returning the newly added DOF or the existing one if it already exists.
     * @details This function checks if the DOF for the given variable already exists in the node's DOF list. If it exists, 
     * it returns the existing DOF. Otherwise, it creates a new DOF for the given variable and adds it to the list.
     * The DOF list is then sorted to maintain order.
     * @tparam TVariableType The type of the variable associated with the degree of freedom.
     * @param rDofVariable The variable associated with the degree of freedom to be added.
     * @return A pointer to the added or existing DOF.
     */
    template<class TVariableType>
    inline typename DofType::Pointer pAddDof(TVariableType const& rDofVariable)
    {
        KRATOS_TRY

        for(auto it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == rDofVariable){
                return (*it_dof).get();
            }
        }

        mDofs.push_back(Kratos::make_unique<DofType>(&mNodalData, rDofVariable));

        DofType* p_new_dof = mDofs.back().get();

        SortDofs();

        return p_new_dof;

        KRATOS_CATCH(*this);
    }

    /**
     * @brief Adds a degree of freedom (DOF) to the node, returning the newly added DOF or the existing one if it already exists.
     * @details This function checks if the DOF for the given variable already exists in the node's DOF list. If it exists, 
     * it updates the reaction and returns the existing DOF. If not, it creates a new DOF based on the provided `SourceDof`.
     * The DOF list is then sorted to maintain order.
     * @param SourceDof The DOF to be added or whose reaction should be updated.
     * @return A pointer to the added or existing DOF.
     */
    inline typename DofType::Pointer pAddDof(DofType const& SourceDof)
    {
        KRATOS_TRY

        for(auto it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == SourceDof.GetVariable()) {
                if((*it_dof)->GetReaction() != SourceDof.GetReaction()) {
                    **it_dof = SourceDof;
                    (*it_dof)->SetNodalData(&mNodalData);
                }
                return (*it_dof).get();
            }
        }

        mDofs.push_back(Kratos::make_unique<DofType>(SourceDof));
        mDofs.back()->SetNodalData(&mNodalData);

        DofType* p_new_dof = mDofs.back().get();

        SortDofs();

        return p_new_dof;

        KRATOS_CATCH(*this);
    }

    /**
     * @brief Adds a degree of freedom (DOF) to the node, returning the newly added DOF or the existing one if it already exists, with a reaction.
     * @details This function checks if the DOF for the given variable already exists in the node's DOF list. If it exists, 
     * it updates the reaction and returns the existing DOF. Otherwise, it creates a new DOF with both the variable 
     * and the reaction, and adds it to the list. The DOF list is then sorted to maintain order.
     * @tparam TVariableType The type of the variable associated with the degree of freedom.
     * @tparam TReactionType The type of the reaction associated with the degree of freedom.
     * @param rDofVariable The variable associated with the degree of freedom to be added.
     * @param rDofReaction The reaction associated with the degree of freedom to be added.
     * @return A pointer to the added or existing DOF.
     */
    template<class TVariableType, class TReactionType>
    inline typename DofType::Pointer pAddDof(TVariableType const& rDofVariable, TReactionType const& rDofReaction)
    {
        KRATOS_TRY

        for(auto it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == rDofVariable){
                (*it_dof)->SetReaction(rDofReaction);
                return (*it_dof).get();
            }
        }

        mDofs.push_back(Kratos::make_unique<DofType>(&mNodalData, rDofVariable, rDofReaction));

        DofType* p_new_dof = mDofs.back().get();

        SortDofs();

        return p_new_dof;

        KRATOS_CATCH(*this);
    }

    /**
     * @brief Adds a degree of freedom (DOF) to the node, returning the newly added DOF or the existing one if it already exists.
     * @details This function checks if the DOF for the given variable already exists in the node's DOF list. If it exists, 
     * it returns the existing DOF. Otherwise, it creates a new DOF for the given variable and adds it to the list.
     * The DOF list is then sorted to maintain order.
     * @tparam TVariableType The type of the variable associated with the degree of freedom.
     * @param rDofVariable The variable associated with the degree of freedom to be added.
     * @return A reference to the added or existing DOF.
     */
    template<class TVariableType>
    inline DofType& AddDof(TVariableType const& rDofVariable)
    {
        KRATOS_TRY

        for(auto it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == rDofVariable){
                return **it_dof;
            }
        }

        mDofs.push_back(Kratos::make_unique<DofType>(&mNodalData, rDofVariable));

        DofType* p_new_dof = mDofs.back().get();

        SortDofs();

        return *p_new_dof;

        KRATOS_CATCH(*this);
    }

    /**
     * @brief Adds a degree of freedom (DOF) to the node, returning the newly added DOF or the existing one if it already exists, with a reaction.
     * @details This function checks if the DOF for the given variable already exists in the node's DOF list. If it exists, 
     * it updates the reaction and returns the existing DOF. Otherwise, it creates a new DOF with both the variable 
     * and the reaction, and adds it to the list. The DOF list is then sorted to maintain order.
     * @tparam TVariableType The type of the variable associated with the degree of freedom.
     * @tparam TReactionType The type of the reaction associated with the degree of freedom.
     * @param rDofVariable The variable associated with the degree of freedom to be added.
     * @param rDofReaction The reaction associated with the degree of freedom to be added.
     * @return A reference to the added or existing DOF.
     */
    template<class TVariableType, class TReactionType>
    inline DofType& AddDof(TVariableType const& rDofVariable, TReactionType const& rDofReaction)
    {
        KRATOS_TRY

        for(auto it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == rDofVariable){
                (*it_dof)->SetReaction(rDofReaction);
                return **it_dof;
            }
        }

        mDofs.push_back(Kratos::make_unique<DofType>(&mNodalData, rDofVariable, rDofReaction));

        DofType* p_new_dof = mDofs.back().get();

        SortDofs();

        return *p_new_dof;

        KRATOS_CATCH(*this);
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Checks if the degree of freedom (DOF) for a given variable is present on the node.
     * @details This function iterates over the DOFs associated with the node and checks whether
     * any of the DOFs correspond to the provided `rDofVariable`.
     * @param rDofVariable The variable associated with the degree of freedom to be checked.
     * @return `true` if the degree of freedom for the given variable is present, otherwise `false`.
     */
    inline bool HasDofFor(const VariableData& rDofVariable) const
    {
        for(auto it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == rDofVariable){
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Checks if the degree of freedom (DOF) for a given variable is fixed on the node.
     * @details This function iterates over the DOFs associated with the node and checks whether
     * any of the DOFs corresponding to the given `rDofVariable` are marked as fixed.
     * @param rDofVariable The variable associated with the degree of freedom to be checked.
     * @return `true` if the degree of freedom for the given variable is fixed, otherwise `false`.
     */
    inline bool IsFixed(const VariableData& rDofVariable) const
    {
        for(auto it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == rDofVariable){
                return (*it_dof)->IsFixed();
            }
        }
        return false;
    }

    /**
     * @brief Checks if the GeometricalObject is active
     * @return True by default, otherwise depending on the ACTIVE flag
     */
    inline bool IsActive() const 
    {
        return IsDefined(ACTIVE) ? Is(ACTIVE) : true;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Node #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
        if(!mDofs.empty())
            rOStream << std::endl << "    Dofs :" << std::endl;
        for(typename DofsContainerType::const_iterator i = mDofs.begin() ; i != mDofs.end() ; i++)
            rOStream << "        " << (*i)->Info() << std::endl;
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

    /// The nodal data (historical variables)
    NodalData mNodalData;

    /// Storage for the dof of the node
    DofsContainerType  mDofs;

    /// A container with data related to this node (non-historical variables) 
    DataValueContainer mData;

    /// Initial Position of the node
    Point mInitialPosition;

    /// The lock object of the node
    LockObject mNodeLock;

    ///@}
    ///@name Private Operators
    ///@{

    // This block is needed for refcounting
    mutable std::atomic<int> mReferenceCounter{0};

    /**
     * @brief Increments the reference counter of the given NodeType object.
     * @details This function increments the reference counter of the NodeType object `x`
     * by 1 using relaxed memory ordering.
     * @param x Pointer to the NodeType object whose reference counter is to be incremented.
     */
    friend void intrusive_ptr_add_ref(const Node* x)
    {
        x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
    }

    /**
     * @brief Decrements the reference counter of the given NodeType object and deletes it if zero.
     * @details This function decrements the reference counter of the NodeType object `x`
     * by 1 using release memory ordering. If the counter reaches zero after decrement,
     * the function ensures memory visibility using acquire memory ordering and deletes `x`.
     * @param x Pointer to the NodeType object whose reference counter is to be decremented.
     */
    friend void intrusive_ptr_release(const Node* x)
    {
        if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
        std::atomic_thread_fence(std::memory_order_acquire);
        delete x;
        }
    }

    /**
     * @brief Sorts the degrees of freedom (DOFs) based on their variable keys.
     * @details This function sorts the `mDofs` container (which holds unique pointers to `DofType` objects)
     * in ascending order of the keys of their associated variables. The sorting is done using
     * the `std::sort` algorithm and compares the `Key()` of the `GetVariable()` for each DOF.
     */
    void SortDofs()
    {
        std::sort(mDofs.begin(), mDofs.end(), [](Kratos::unique_ptr<DofType> const& First, Kratos::unique_ptr<DofType> const& Second)->bool{
            return First->GetVariable().Key() < Second->GetVariable().Key();
        });
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    /**
     * @brief The save operation which copies the database of the class
     * @param rSerializer The serializer used to preserve the information
     */
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Point );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
        rSerializer.save("NodalData", &mNodalData); // Storing it as pointer to be shared by Dof pointer
        rSerializer.save("Data", mData);
        rSerializer.save("Initial Position", mInitialPosition);
        rSerializer.save("Data", mDofs);
    }    

    /**
     * The load operation which restores the database of the class
     * @param rSerializer The serializer used to preserve the information
     */
    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Point );
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
        NodalData* p_nodal_data = &mNodalData;
        rSerializer.load("NodalData", p_nodal_data);
        rSerializer.load("Data", mData);
        rSerializer.load("Initial Position", mInitialPosition);
        rSerializer.load("Data", mDofs);
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

}; // Class Node

///@}

// template class KRATOS_API(KRATOS_CORE) KratosComponents<Node >;
// template class KRATOS_API(KRATOS_CORE) KratosComponents<Node::Pointer >;

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Node& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Node& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : ";
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.