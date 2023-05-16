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
#include "containers/variables_list_data_value_container.h"
#include "containers/flags.h"
#include "intrusive_ptr/intrusive_ptr.hpp"
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
 * @class Node
 * @brief This class defines the node
 * @details The node class from Kratos is defined in this class
 * @ingroup KratosCore
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) Node

    : public Point, public Flags
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Node
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Node);

    /// Base type
    using BaseType = Point;

    /// Point type
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
     * @brief Default constructor. Constructs a Node object.
     */
    Node();

    /**
     * @brief Constructs a Node object with a given ID.
     * @param NewId The ID of the new node.
     */
    explicit Node(const IndexType NewId );

    /**
     * @brief Constructs a 1D Node object with a given ID and X coordinate.
     * @param NewId The ID of the new node.
     * @param NewX The X coordinate of the new node.
     */
    Node(const IndexType NewId, const double NewX);

    /**
     * @brief Constructs a 2D Node object with a given ID, X and Y coordinates.
     * @param NewId The ID of the new node.
     * @param NewX The X coordinate of the new node.
     * @param NewY The Y coordinate of the new node.
     */
    Node(const IndexType NewId, const double NewX, const double NewY);

    /**
     * @brief Constructs a 3D Node object with a given ID, X, Y, and Z coordinates.
     * @param NewId The ID of the new node.
     * @param NewX The X coordinate of the new node.
     * @param NewY The Y coordinate of the new node.
     * @param NewZ The Z coordinate of the new node.
     */
    Node(const IndexType NewId, const double NewX, const double NewY, const double NewZ);

    /**
     * @brief Constructs a Node object with a given ID and Point.
     * @param NewId The ID of the new node.
     * @param rThisPoint The point of the new node.
     */
    Node(const IndexType NewId, PointType const& rThisPoint);

    /**
     * @brief Copy constructor. This constructor is deleted to prevent copying.
     * @param rOtherNode The node to be copied.
     */
    Node(Node const& rOtherNode) = delete;

    /**
     * @brief Constructs a Node object with a given ID and coordinates stored in a given vector.
     * @tparam TVectorType The type of the vector.
     * @param NewId The ID of the new node.
     * @param rOtherCoordinates The coordinates of the new node.
     */
    template<class TVectorType>
    Node(const IndexType NewId, vector_expression<TVectorType> const&  rOtherCoordinates)
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
     * @brief Constructs a Node object with a given ID and coordinates stored in a given std::vector.
     * @param NewId The ID of the new node.
     * @param rOtherCoordinates The coordinates of the new node.
     */
    Node(const IndexType NewId, std::vector<double> const&  rOtherCoordinates);

    /**
     * @brief Constructs a 3D Node object with a given ID, X, Y, and Z coordinates, variables list, data, and queue size.
     * @param NewId The ID of the new node.
     * @param NewX The X coordinate of the new node.
     * @param NewY The Y coordinate of the new node.
     * @param NewZ The Z coordinate of the new node.
     * @param pVariablesList The variables list of the new node.
     * @param ThisData The data of the new node.
     * @param NewQueueSize The queue size of the new node.
     */
    Node(const IndexType NewId, const double NewX, const double NewY, const double NewZ, VariablesList::Pointer  pVariablesList, BlockType const * ThisData, SizeType NewQueueSize = 1);

    /**
     * @brief Clones the current Node and returns a pointer to the new Node.
     * @return typename Node::Pointer A pointer to the new Node.
     */
    typename Node::Pointer Clone();

    /// Destructor.
    ~Node() override;

    //*********************************************
    //public API of intrusive_ptr
    inline unsigned int use_count() const noexcept
    {
        return mReferenceCounter;
    }

    //*********************************************

    /**
     * @brief This method returns the current Id of the node.
     * @return The Id of the node.
     */
    inline IndexType Id() const
    {
        return mNodalData.Id();
    }

    /**
     * @brief This method returns the current Id of the node.
     * @return The Id of the node.
     */
    inline IndexType GetId() const
    {
        return mNodalData.Id();
    }

    /**
     * @brief This method sets the Id of the node.
     * @param NewId The new Id of the node.
     */
    inline void SetId(const IndexType NewId)
    {
        mNodalData.SetId(NewId);
    }

    /**
     * @brief This method returns the lock of the node.
     * @return The lock of the node.
    */
    inline LockObject& GetLock()
    {
        return mNodeLock;
    }

    /**
     * @brief Locks the node lock for thread-safe access. This function is inline to reduce function call overhead.
     */
    inline void SetLock()
    {
        mNodeLock.lock();
    }

    /**
     * @brief Unlocks the mutex lock for the calling thread that is currently locked by 
     * the calling thread. This function is inline to reduce function call overhead.
     */
    inline void UnSetLock()
    {
        mNodeLock.unlock();
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Node& operator=(const Node& rOther);

    /// Comparison operator. 
    bool operator==(const Node& rOther);

    /**
     * @brief Overloads the function call operator to return the value of the given variable at the specified solution step index.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @param SolutionStepIndex The index of the solution step.
     * @return typename TVariableType::Type& A reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type& operator()(const TVariableType& rThisVariable, const IndexType SolutionStepIndex)
    {
        return GetSolutionStepValue(rThisVariable, SolutionStepIndex);
    }

    /**
     * @brief Overloads the function call operator to return the value of the given variable from the current solution step.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @return typename TVariableType::Type& A reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type& operator()(const TVariableType& rThisVariable)
    {
        return GetSolutionStepValue(rThisVariable);
    }

    /**
     * @brief Overloads the subscript operator to return the value of the given variable.
     * @tparam TDataType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @return TDataType& A reference to the value of the variable.
     */
    template<class TDataType>
    TDataType& operator[](const Variable<TDataType>& rThisVariable)
    {
        return GetValue(rThisVariable);
    }

    /**
     * @brief Overloads the subscript operator to return a const reference to the value of the given variable.
     * @tparam TDataType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @return const TDataType& A const reference to the value of the variable.
     */
    template<class TDataType>
    const TDataType& operator[](const Variable<TDataType>& rThisVariable) const
    {
        return GetValue(rThisVariable);
    }

    /**
     * @brief Overloads the subscript operator to return a reference to the value at the specified index.
     * @param ThisIndex The index of the value.
     * @return double& A reference to the value at the index.
     */
    double& operator[](const IndexType ThisIndex);

    /**
     * @brief Overloads the subscript operator to return a const reference to the value at the specified index.
     * @param ThisIndex The index of the value.
     * @return double A const reference to the value at the index.
     */
    double operator[](const IndexType ThisIndex) const;

    ///@}
    ///@name Nodal Data
    ///@{

    /**
     * @brief This method creates the solution step data for the node.
    */
    void CreateSolutionStepData();

    /**
     * @brief This method clones the solution step data for the node.
    */
    void CloneSolutionStepData();

    /**
     * @brief This overrides the current solution step data with the data from the given node.
     * @param SourceSolutionStepIndex The index of the solution step to be copied.
     * @param DestinationSourceSolutionStepIndex The index of the solution step to be overwritten.
    */
    void OverwriteSolutionStepData(
        const IndexType SourceSolutionStepIndex,
        const IndexType DestinationSourceSolutionStepIndex
        )
    {
        SolutionStepData().AssignData(SolutionStepData().Data(SourceSolutionStepIndex), DestinationSourceSolutionStepIndex);
    }

    /**
     * @brief This method clears the solution step data for the node.
    */
    void ClearSolutionStepsData()
    {
        SolutionStepData().Clear();
    }

    /**
     * @brief This method sets the solution step data for the node.
    */
    void SetSolutionStepVariablesList(VariablesList::Pointer pVariablesList)
    {
        SolutionStepData().SetVariablesList(pVariablesList);
    }

    /**
     * @brief This method returns the solution step data for the node.
    */
    VariablesListDataValueContainer& SolutionStepData()
    {
        return mNodalData.GetSolutionStepData();
    }

    /**
     * @brief This method returns the solution step data for the node.
     * @details This method is const and should be used only when the data is not going to be modified.
    */
    const VariablesListDataValueContainer& SolutionStepData() const
    {
        return mNodalData.GetSolutionStepData();
    }

    /**
     * @brief Returns a reference to the value of the given variable from the current solution step.
     * @details This method accesses the value of the given variable from the SolutionStepData.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @return typename TVariableType::Type& A reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable)
    {
        return SolutionStepData().GetValue(rThisVariable);
    }

    /**
     * @brief Returns a const reference to the value of the given variable from the current solution step.
     * @details This method accesses the value of the given variable from the SolutionStepData.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @return typename TVariableType::Type const& A const reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type const& GetSolutionStepValue(const TVariableType& rThisVariable) const
    {
        return SolutionStepData().GetValue(rThisVariable);
    }

    /**
     * @brief Returns a reference to the value of the given variable at the specified solution step index.
     * @details This method accesses the value of the given variable from the SolutionStepData.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @param SolutionStepIndex The index of the solution step.
     * @return typename TVariableType::Type& A reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable, const IndexType SolutionStepIndex)
    {
        return SolutionStepData().GetValue(rThisVariable, SolutionStepIndex);
    }

    /**
     * @brief Returns a const reference to the value of the given variable at the specified solution step index.
     * @details This method accesses the value of the given variable from the SolutionStepData.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @param SolutionStepIndex The index of the solution step.
     * @return typename TVariableType::Type const& A const reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type const& GetSolutionStepValue(const TVariableType& rThisVariable, const IndexType SolutionStepIndex) const
    {
        return SolutionStepData().GetValue(rThisVariable, SolutionStepIndex);
    }

    /**
     * @brief Checks if the SolutionStepData has the given variable.
     * @details This method checks if the SolutionStepData contains the given variable.
     * @param rThisVariable The variable to check.
     * @return bool True if the SolutionStepData contains the variable, false otherwise.
     */
    bool SolutionStepsDataHas(const VariableData& rThisVariable) const
    {
        return SolutionStepData().Has(rThisVariable);
    }

    //*******************************************************************************************
    //By Riccardo
    //very similar to the one before BUT throws an error if the variable does not exist
    /**
     * @brief Returns a reference to the value of the given variable from the current solution step.
     * @details This method accesses the value of the given variable using a fast access method.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @return typename TVariableType::Type& A reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable)
    {
        return SolutionStepData().FastGetValue(rThisVariable);
    }

    /**
     * @brief Returns a const reference to the value of the given variable from the current solution step.
     * @details This method accesses the value of the given variable using a fast access method.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @return const typename TVariableType::Type& A const reference to the value of the variable.
     */
    template<class TVariableType>
    const typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable) const
    {
        return SolutionStepData().FastGetValue(rThisVariable);
    }

    /**
     * @brief Returns a reference to the value of the given variable at the specified solution step index.
     * @details This method accesses the value of the given variable using a fast access method.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @param SolutionStepIndex The index of the solution step.
     * @return typename TVariableType::Type& A reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, const IndexType SolutionStepIndex)
    {
        return SolutionStepData().FastGetValue(rThisVariable, SolutionStepIndex);
    }

    /**
     * @brief Returns a const reference to the value of the given variable at the specified solution step index.
     * @details This method accesses the value of the given variable using a fast access method.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @param SolutionStepIndex The index of the solution step.
     * @return const typename TVariableType::Type& A const reference to the value of the variable.
     */
    template<class TVariableType>
    const typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, const IndexType SolutionStepIndex) const
    {
        return SolutionStepData().FastGetValue(rThisVariable, SolutionStepIndex);
    }

    /**
     * @brief Returns a reference to the value of the given variable at the specified solution step index and position.
     * @details This method accesses the value of the given variable using a fast access method.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @param SolutionStepIndex The index of the solution step.
     * @param ThisPosition The position within the solution step.
     * @return typename TVariableType::Type& A reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, const IndexType SolutionStepIndex, IndexType ThisPosition)
    {
        return SolutionStepData().FastGetValue(rThisVariable, SolutionStepIndex, ThisPosition);
    }

    /**
     * @brief Returns a reference to the value of the given variable at the current solution step and specified position.
     * @details This method accesses the value of the given variable using a fast access method.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @param ThisPosition The position within the solution step.
     * @return typename TVariableType::Type& A reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type& FastGetCurrentSolutionStepValue(const TVariableType& rThisVariable, IndexType ThisPosition)
    {
        return SolutionStepData().FastGetCurrentValue(rThisVariable, ThisPosition);
    }

//*******************************************************************************************

    /**
     * @brief Returns a reference to the value of the given variable.
     * @details This method accesses the value of the given variable from the data container.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @return typename TVariableType::Type& A reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rThisVariable)
    {
        return mData.GetValue(rThisVariable);
    }

    /**
     * @brief Returns a const reference to the value of the given variable.
     * @details This method accesses the value of the given variable from the data container.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @return typename TVariableType::Type const& A const reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type const& GetValue(const TVariableType& rThisVariable) const
    {
        return mData.GetValue(rThisVariable);
    }

    /**
     * @brief Returns a reference to the value of the given variable at the specified solution step index.
     * @details If the data container does not have the variable, the method gets the value from SolutionStepData.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @param SolutionStepIndex The index of the solution step.
     * @return typename TVariableType::Type& A reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rThisVariable, const IndexType SolutionStepIndex)
    {
        if(!mData.Has(rThisVariable))
            return SolutionStepData().GetValue(rThisVariable, SolutionStepIndex);

        return mData.GetValue(rThisVariable);
    }

    /**
     * @brief Returns a const reference to the value of the given variable at the specified solution step index.
     * @details If the data container does not have the variable, the method gets the value from SolutionStepData.
     * @tparam TVariableType The type of the variable.
     * @param rThisVariable The variable of which the value is needed.
     * @param SolutionStepIndex The index of the solution step.
     * @return typename TVariableType::Type const& A const reference to the value of the variable.
     */
    template<class TVariableType>
    typename TVariableType::Type const& GetValue(const TVariableType& rThisVariable, const IndexType SolutionStepIndex) const
    {
        if(!mData.Has(rThisVariable))
            return SolutionStepData().GetValue(rThisVariable, SolutionStepIndex);

        return mData.GetValue(rThisVariable);
    }

    /**
     * @brief Sets the value of a given variable type.
     * @param rThisVariable the variable type to set the value of
     * @param rValue the value of the variable type to set
     */
    template<class TVariableType>
    void SetValue(const TVariableType& rThisVariable, typename TVariableType::Type const& rValue)
    {
        mData.SetValue(rThisVariable, rValue);
    }

    /**
     * @brief Checks if the specified variable is present in the data.
     * @param rThisVariable the variable to check for presence
     * @return true if the variable is present, false otherwise
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
     * @brief Fixes the degree of freedom associated with the given variable type.
     * @param rDofVariable the variable type for which the degree of freedom is to be fixed
     * @throws If KRATOS_DEBUG is defined, the function checks whether it is called from within a parallel region and
     * throws an error if it is.
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
        if(OpenMPUtils::IsInParallel() != 0)
        {
            KRATOS_ERROR << "Attempting to Fix the variable: " << rDofVariable << " within a parallel region. This is not permitted. Create the Dof first by pAddDof" << std::endl;
        }
#endif
        pAddDof(rDofVariable)->FixDof();
    }

    /**
     * @brief Free the degree of freedom (DOF) associated with the given variable type.
     * @param rDofVariable the variable type of the DOF to be freed.
     * @details This function iterates over the DOFs stored in mDofs and frees the DOF with the given variable type
     * (rDofVariable). If the DOF is found and freed, the function returns immediately. If the DOF is not
     * found, a new DOF is added to mDofs with the given variable type and then freed.
     * @throws If KRATOS_DEBUG is defined, the function checks whether it is called from within a parallel region and
     * throws an error if it is.
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
        if(OpenMPUtils::IsInParallel() != 0)
        {
            KRATOS_ERROR << "Attempting to Free the variable: " << rDofVariable << " within a parallel region. This is not permitted. Create the Dof first by pAddDof" << std::endl;
        }
#endif
        pAddDof(rDofVariable)->FreeDof();
    }

    /**
     * @brief Get the buffer size
     * @return The buffer size
     */
    IndexType GetBufferSize() const;

    /**
     * @brief Set the buffer size
     * @param NewBufferSize The new buffer size
     */
    void SetBufferSize(const IndexType NewBufferSize);

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Get the non-historical data container
     * @details Deprecated method. Use 'GetData()' instead.
     * @return The non-historical data container
     */
    KRATOS_DEPRECATED_MESSAGE("This method is deprecated. Use 'GetData()' instead.")
    inline DataValueContainer& Data()
    {
        return mData;
    }

    /**
     * @brief Get the non-historical data container
     * @return The non-historical data container
     */
    inline DataValueContainer& GetData()
    {
        return mData;
    }

    /**
     * @brief Get the non-historical data container
     * @details const version
     * @return The non-historical data container
     */
    inline const DataValueContainer& GetData() const
    {
        return mData;
    }

    /**
     * @brief Get the initial coordinates as a point
     * @return The initial position as a point
     * @details const version
     */
    inline const PointType& GetInitialPosition() const
    {
        return mInitialPosition;
    }

    /**
     * @brief Get the initial coordinates as a point
     * @return The initial position as a point
     */
    inline PointType& GetInitialPosition()
    {
        return mInitialPosition;
    }

    /**
     * @brief Get the initial X coordinate
     */
    inline double& X0()
    {
        return mInitialPosition.X();
    }

    /**
     * @brief Get the initial Y coordinate
     */
    inline double& Y0()
    {
        return mInitialPosition.Y();
    }

    /**
     * @brief Get the initial Z coordinate
     */
    inline double& Z0()
    {
        return mInitialPosition.Z();
    }

    /**
     * @brief Get the initial X coordinate
     * @details const version
     */
    inline double X0() const
    {
        return mInitialPosition.X();
    }

    /**
     * @brief Get the initial Y coordinate
     * @details const version
     */
    inline double Y0() const
    {
        return mInitialPosition.Y();
    }

    /**
     * @brief Get the initial Z coordinate
     * @details const version
     */
    inline double Z0() const
    {
        return mInitialPosition.Z();
    }

    /**
     * @brief Set the initial coordinates
     * @param rNewInitialPosition The new initial position
     */
    void SetInitialPosition(const PointType& rNewInitialPosition);

    /**
     * @brief Set the initial coordinates
     * @param X The new initial X coordinate
     * @param Y The new initial Y coordinate
     * @param Z The new initial Z coordinate
     */
    void SetInitialPosition(
        const double X,
        const double Y,
        const double Z
        );

    /**
     * @brief Get the list of variables
     * @return The list of variables
     */
    VariablesList::Pointer pGetVariablesList();

    /**
     * @brief Get the list of variables
     * @details const version
     * @return The list of variables
     */
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

    /**
     * @brief Returns all of the Dofs
     */
    DofsContainerType& GetDofs();

    /**
     * @brief Returns all of the Dofs
     * @details const version
     */
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
     * @brief Adds a Degree of Freedom (DoF) with the given variable to the current node if it does not already exist.
     * @param rDofVariable the variable of the DoF to add
     * @return a pointer to the newly added DoF
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
     * @brief Adds a degree of freedom (DOF) to the mesh node and returns the pointer to it.
     * @param SourceDof the DOF to be added
     * @return a pointer to the added DOF
     */
    inline typename DofType::Pointer pAddDof(DofType const& SourceDof)
    {
        KRATOS_TRY

        for(auto it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == SourceDof.GetVariable()){
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
     * @brief Adds a degree of freedom to the nodal data.
     * @param rDofVariable the degree of freedom variable
     * @param rDofReaction the reaction type of the degree of freedom
     * @return a pointer to the new degree of freedom
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
     * @brief Add a Dof to the node and return new added dof or existed one.
     * @param rDofVariable The variable
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
     * @brief Add a Dof to the node and return new added dof or existed one.
     * @param rDofVariable The variable
     * @param rDofReaction The reaction variable
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
     * @brief Determine if the provided `VariableData` has a degree of freedom (DOF).
     * @param rDofVariable the `VariableData` to check for a DOF.
     * @return `true` if a DOF exists for the provided `VariableData`, `false` otherwise.
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
     * @brief Determines whether a given VariableData object is fixed or not by searching
     * for it in the list of DOFs.
     * @param rDofVariable the VariableData object to search for
     * @return true if the VariableData object is fixed, false otherwise
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

    /// Non-historical storage for the nodal data
    NodalData mNodalData;

    /// storage for the dof of the node
    DofsContainerType  mDofs;

    /// A pointer to data related to this node.
    DataValueContainer mData;

    ///Initial Position of the node
    PointType mInitialPosition;

    /// The lock object for node
    LockObject mNodeLock;

    ///@}
    ///@name Private Operators
    ///@{
    //*********************************************
    //this block is needed for refcounting
    mutable std::atomic<int> mReferenceCounter{0};

    /**
     * @brief A friend function that increments the reference counter of a node type.
     * @param x a pointer to the node type for which the reference counter is incremented.
     */
    friend void intrusive_ptr_add_ref(const Node* x)
    {
        x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
    }

    /**
     * @brief Decrements the reference count of the given Node object, and deletes it if
     * the reference count becomes zero.
     * @param x Pointer to the Node object whose reference count needs to be decremented
     */
    friend void intrusive_ptr_release(const Node* x)
    {
        if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
        std::atomic_thread_fence(std::memory_order_acquire);
        delete x;
        }
    }
    //*********************************************

    /**
     * @brief This method sorts the dofs of the node
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
     * @brief This method is used to serialize the node
     */
    void save(Serializer& rSerializer) const override;

    /**
     * @brief This method is used to deserialize the node
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
