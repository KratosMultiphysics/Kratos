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
class Node 
    : public Point, public Flags
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Node
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Node);

    /// Base type
    using BaseType = Point;

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

    /// Default constructor.
    Node();

    explicit Node(IndexType NewId )
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

    /// 1d constructor.
    Node(IndexType NewId, const double NewX);

    /// 2d constructor.
    Node(IndexType NewId, const double NewX, const double NewY);

    /// 3d constructor.
    Node(IndexType NewId, const double NewX, const double NewY, const double NewZ);

    /// Point constructor.
    Node(IndexType NewId, Point const& rThisPoint);

    /** Copy constructor. Initialize this node with given node.*/
    Node(Node const& rOtherNode) = delete;

    /**
     * Constructor using coordinates stored in given array. Initialize
    this point with the coordinates in the array. */
    template<class TVectorType>
    Node(IndexType NewId, vector_expression<TVectorType> const&  rOtherCoordinates)
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

    /** Constructor using coordinates stored in given std::vector. Initialize
    this point with the coordinates in the array. */
    Node(IndexType NewId, std::vector<double> const&  rOtherCoordinates);

    /// 3d with variables list and data constructor.
    Node(IndexType NewId, const double NewX, const double NewY, const double NewZ, VariablesList::Pointer  pVariablesList, BlockType const * ThisData, SizeType NewQueueSize = 1);

    /**
     * @brief Creates a clone of the current node.
     * @details This function creates a new node with the same ID and coordinates as
     * the current node, copies its nodal data, DOFs, and flags.
     * @return Pointer to the cloned node.
     */
    typename Node::Pointer Clone();

    /**
     * @brief Destructor for the Node class.
     * @details Clears any stored solution step data before destruction.
     */
    ~Node() override;

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
    IndexType Id() const;

    /**
     * @brief Retrieves the node ID (alternative method).
     * @return The ID of the node.
     */
    IndexType GetId() const;

    /**
     * @brief Sets a new ID for the node.
     * @param NewId The new ID to assign.
     */
    void SetId(IndexType NewId);

    /**
     * @brief Retrieves the lock object for thread safety.
     * @return A reference to the lock object.
     */
    LockObject& GetLock();

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
    Node& operator=(const Node& rOther);

    /**
     * @brief Equality operator.
     * @details Compares two nodes based on their position.
     * @param rOther The node to compare with.
     * @return True if nodes are equal, false otherwise.
     */
    bool operator==(const Node& rOther);

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
    double& operator[](IndexType ThisIndex);

    /**
     * @brief Const accessor for node coordinates.
     * @param ThisIndex The index of the coordinate (0 = X, 1 = Y, 2 = Z).
     * @return The coordinate value.
     */
    double operator[](IndexType ThisIndex) const;

    ///@}
    ///@name Nodal Data
    ///@{

    /**
     * @brief Creates a new solution step data entry.
     * @details Adds a new solution step at the front of the solution step data container.
     */
    void CreateSolutionStepData();

    /**
     * @brief Clones the current solution step data.
     * @details Copies the front solution step data to create a new step.
     */
    void CloneSolutionStepData();

    /**
     * @brief Overwrites solution step data at a given index.
     * @details Copies the data from a source solution step index to a destination index.
     * @param SourceSolutionStepIndex The index of the source solution step.
     * @param DestinationSourceSolutionStepIndex The index of the destination step where the data will be assigned.
     */
    void OverwriteSolutionStepData(IndexType SourceSolutionStepIndex, IndexType DestinationSourceSolutionStepIndex);

    /**
     * @brief Clears all stored solution step data.
     * @details Removes all solution step data entries.
     */
    void ClearSolutionStepsData();

    /**
     * @brief Sets the list of solution step variables.
     * @details Assigns a new variables list to the solution step data.
     * @param pVariablesList Pointer to the new list of variables.
     */
    void SetSolutionStepVariablesList(VariablesList::Pointer pVariablesList);

    /**
     * @brief Retrieves the solution step data container.
     * @return Reference to the solution step data container.
     */
    VariablesListDataValueContainer& SolutionStepData();

    /**
     * @brief Retrieves the solution step data container (const version).
     * @return Const reference to the solution step data container.
     */
    const VariablesListDataValueContainer& SolutionStepData() const;

    /**
     * @deprecated This method is deprecated. Use GetData() instead.
     * @brief Retrieves the data container for the node.
     * @return Reference to the data container.
     */
    KRATOS_DEPRECATED_MESSAGE("This method is deprecated. Use 'GetData()' instead.")
    DataValueContainer& Data();

    /**
     * @brief Retrieves the data container for the node.
     * @return Reference to the data container.
     */
    DataValueContainer& GetData();

    /**
     * @brief Retrieves the data container for the node (const version).
     * @return Const reference to the data container.
     */
    const DataValueContainer& GetData() const;

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
    bool SolutionStepsDataHas(const VariableData& rThisVariable) const;

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

    IndexType GetBufferSize() const;

    void SetBufferSize(IndexType NewBufferSize);

    ///@}
    ///@name Access
    ///@{

    const Point& GetInitialPosition() const;

    Point& GetInitialPosition();

    double& X0();

    double& Y0();
    
    double& Z0();

    double X0() const;
    
    double Y0() const;
    
    double Z0() const;

    void SetInitialPosition(const Point& NewInitialPosition);

    void SetInitialPosition(double X,double Y, double Z);

    VariablesList::Pointer pGetVariablesList();

    const VariablesList::Pointer pGetVariablesList() const;

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
        typename DofsContainerType::const_iterator it_begin = mDofs.begin();
        typename DofsContainerType::const_iterator it_end = mDofs.end();
        typename DofsContainerType::const_iterator it;
        //if the guess is exact return the guess
        if(pos < it_end-it_begin)
        {
            it = it_begin + pos;
            if( (*it)->GetVariable() == rDofVariable)
            {
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

    /** returns all of the Dofs  */
    DofsContainerType& GetDofs();

    const DofsContainerType& GetDofs() const;

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
    inline  bool IsFixed(const VariableData& rDofVariable) const
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
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

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
    void SortDofs();

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
    void save(Serializer& rSerializer) const override;

    /**
     * The load operation which restores the database of the class
     * @param rSerializer The serializer used to preserve the information
     */
    void load(Serializer& rSerializer) override;

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
