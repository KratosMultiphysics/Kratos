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
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>
#include <atomic>

// External includes

// Project includes
#include "includes/define.h"
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

/// This class defines the node
/** The node class from Kratos is defined in this class
*/
class Node : public Point, public Flags
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Node
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Node);

    /// Node definition
    using NodeType = Node;

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

    /// Default constructor.
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
    Node(IndexType NewId, double const& NewX)
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

    /// 2d constructor.
    Node(IndexType NewId, double const& NewX, double const& NewY)
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

    /// 3d constructor.
    Node(IndexType NewId, double const& NewX, double const& NewY, double const& NewZ)
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

    /// Point constructor.
    Node(IndexType NewId, PointType const& rThisPoint)
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
    Node(IndexType NewId, std::vector<double> const&  rOtherCoordinates)
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

    /// 3d with variables list and data constructor.
    Node(IndexType NewId, double const& NewX, double const& NewY, double const& NewZ, VariablesList::Pointer  pVariablesList, BlockType const * ThisData, SizeType NewQueueSize = 1)
        : BaseType(NewX, NewY, NewZ)
        , Flags()
        , mNodalData(NewId, pVariablesList,ThisData,NewQueueSize)
        , mDofs()
        , mData()
        , mInitialPosition(NewX, NewY, NewZ)
        , mNodeLock()
    {
    }

    typename Node::Pointer Clone()
    {
        Node::Pointer p_new_node = Kratos::make_intrusive<Node >( this->Id(), (*this)[0], (*this)[1], (*this)[2]);
        p_new_node->mNodalData = this->mNodalData;

        Node::DofsContainerType& my_dofs = (this)->GetDofs();
        for (typename DofsContainerType::const_iterator it_dof = my_dofs.begin(); it_dof != my_dofs.end(); it_dof++)
        {
            p_new_node->pAddDof(**it_dof);
        }

        p_new_node->mData = this->mData;
        p_new_node->mInitialPosition = this->mInitialPosition;

        p_new_node->Set(Flags(*this));

        return p_new_node;
    }

    /// Destructor.
    ~Node() override
    {
        ClearSolutionStepsData();
    }

    //*********************************************
    //public API of intrusive_ptr
    unsigned int use_count() const noexcept
    {
        return mReferenceCounter;
    }
    //*********************************************

    IndexType Id() const
    {
        return mNodalData.Id();
    }

    IndexType GetId() const
    {
        return mNodalData.Id();
    }

    void SetId(IndexType NewId)
    {
        mNodalData.SetId(NewId);
    }

    LockObject& GetLock()
    {
        return mNodeLock;
    }

    inline void SetLock()
    {
        mNodeLock.lock();
    }

    inline void UnSetLock()
    {
        mNodeLock.unlock();
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Node& operator=(const Node& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator =(rOther);

        mNodalData = rOther.mNodalData;

        // Deep copying the dofs
        for(typename DofsContainerType::const_iterator it_dof = rOther.mDofs.begin() ; it_dof != rOther.mDofs.end() ; it_dof++)
        {
            pAddDof(**it_dof);
        }

        mData = rOther.mData;
        mInitialPosition = rOther.mInitialPosition;

        return *this;
    }

    bool operator==(const Node& rOther)
    {
        return PointType::operator ==(rOther);
    }

    template<class TVariableType> typename TVariableType::Type& operator()(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
    {
        return GetSolutionStepValue(rThisVariable, SolutionStepIndex);
    }

    template<class TVariableType> typename TVariableType::Type& operator()(const TVariableType& rThisVariable)
    {
        return GetSolutionStepValue(rThisVariable);
    }

    template<class TDataType> TDataType& operator[](const Variable<TDataType>& rThisVariable)
    {
        return GetValue(rThisVariable);
    }

    template<class TDataType> const TDataType& operator[](const Variable<TDataType>& rThisVariable) const
    {
        return GetValue(rThisVariable);
    }

    double& operator[](IndexType ThisIndex)
    {
        return BaseType::operator[](ThisIndex);
    }

    double operator[](IndexType ThisIndex) const
    {
        return BaseType::operator[](ThisIndex);
    }

    ///@}
    ///@name Nodal Data
    ///@{

    void CreateSolutionStepData()
    {
        SolutionStepData().PushFront();
    }

    void CloneSolutionStepData()
    {
        SolutionStepData().CloneFront();
    }

    void OverwriteSolutionStepData(IndexType SourceSolutionStepIndex, IndexType DestinationSourceSolutionStepIndex)
    {
        SolutionStepData().AssignData(SolutionStepData().Data(SourceSolutionStepIndex), DestinationSourceSolutionStepIndex);
    }

    void ClearSolutionStepsData()
    {
        SolutionStepData().Clear();
    }

    void SetSolutionStepVariablesList(VariablesList::Pointer pVariablesList)
    {
        SolutionStepData().SetVariablesList(pVariablesList);
    }

    VariablesListDataValueContainer& SolutionStepData()
    {
        return mNodalData.GetSolutionStepData();
    }

    const VariablesListDataValueContainer& SolutionStepData() const
    {
        return mNodalData.GetSolutionStepData();
    }

    KRATOS_DEPRECATED_MESSAGE("This method is deprecated. Use 'GetData()' instead.")
    DataValueContainer& Data()
    {
        return mData;
    }

    DataValueContainer& GetData()
    {
        return mData;
    }

    const DataValueContainer& GetData() const
    {
        return mData;
    }

    template<class TVariableType> typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable)
    {
        return SolutionStepData().GetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type const& GetSolutionStepValue(const TVariableType& rThisVariable) const
    {
        return SolutionStepData().GetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
    {
        return SolutionStepData().GetValue(rThisVariable, SolutionStepIndex);
    }

    template<class TVariableType> typename TVariableType::Type const& GetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex) const
    {
        return SolutionStepData().GetValue(rThisVariable, SolutionStepIndex);
    }


    bool SolutionStepsDataHas(const VariableData& rThisVariable) const
    {
        return SolutionStepData().Has(rThisVariable);
    }

    //*******************************************************************************************
    //By Riccardo
    //very similar to the one before BUT throws an error if the variable does not exist
    template<class TVariableType> typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable)
    {
        return SolutionStepData().FastGetValue(rThisVariable);
    }

    template<class TVariableType> const typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable) const
    {
        return SolutionStepData().FastGetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
    {
        return SolutionStepData().FastGetValue(rThisVariable, SolutionStepIndex);
    }

    template<class TVariableType> const typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex) const
    {
        return SolutionStepData().FastGetValue(rThisVariable, SolutionStepIndex);
    }

    template<class TVariableType> typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex, IndexType ThisPosition)
    {
        return SolutionStepData().FastGetValue(rThisVariable, SolutionStepIndex, ThisPosition);
    }

    template<class TVariableType> typename TVariableType::Type& FastGetCurrentSolutionStepValue(const TVariableType& rThisVariable, IndexType ThisPosition)
    {
        return SolutionStepData().FastGetCurrentValue(rThisVariable, ThisPosition);
    }
//*******************************************************************************************

    template<class TVariableType> typename TVariableType::Type& GetValue(const TVariableType& rThisVariable)
    {
        return mData.GetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type const& GetValue(const TVariableType& rThisVariable) const
    {
        return mData.GetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type& GetValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
    {
        if(!mData.Has(rThisVariable))
            return SolutionStepData().GetValue(rThisVariable, SolutionStepIndex);

        return mData.GetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type const& GetValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex) const
    {
        if(!mData.Has(rThisVariable))
            return SolutionStepData().GetValue(rThisVariable, SolutionStepIndex);

        return mData.GetValue(rThisVariable);
    }

    template<class TVariableType>
    void SetValue(const TVariableType& rThisVariable, typename TVariableType::Type const& rValue)
    {
        mData.SetValue(rThisVariable, rValue);
    }

    template<class TDataType> bool Has(const Variable<TDataType>& rThisVariable) const
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
        if(OpenMPUtils::IsInParallel() != 0)
        {
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
        if(OpenMPUtils::IsInParallel() != 0)
        {
            KRATOS_ERROR << "Attempting to Free the variable: " << rDofVariable << " within a parallel region. This is not permitted. Create the Dof first by pAddDof" << std::endl;
        }
#endif
        pAddDof(rDofVariable)->FreeDof();
    }

    IndexType GetBufferSize() const
    {
        return SolutionStepData().QueueSize();
    }

    void SetBufferSize(IndexType NewBufferSize)
    {
        SolutionStepData().Resize(NewBufferSize);
    }

    ///@}
    ///@name Access
    ///@{

    const PointType& GetInitialPosition() const
    {
        return mInitialPosition;
    }
    PointType& GetInitialPosition()
    {
        return mInitialPosition;
    }

    double& X0()
    {
        return mInitialPosition.X();
    }
    double& Y0()
    {
        return mInitialPosition.Y();
    }
    double& Z0()
    {
        return mInitialPosition.Z();
    }

    double X0() const
    {
        return mInitialPosition.X();
    }
    double Y0() const
    {
        return mInitialPosition.Y();
    }
    double Z0() const
    {
        return mInitialPosition.Z();
    }

    void SetInitialPosition(const PointType& NewInitialPosition)
    {
        mInitialPosition.X() = NewInitialPosition.X();
        mInitialPosition.Y() = NewInitialPosition.Y();
        mInitialPosition.Z() = NewInitialPosition.Z();
    }

    void SetInitialPosition(double X,double Y, double Z)
    {
        mInitialPosition.X() = X;
        mInitialPosition.Y() = Y;
        mInitialPosition.Z() = Z;
    }

    VariablesList::Pointer pGetVariablesList()
    {
        return SolutionStepData().pGetVariablesList();
    }

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
    DofsContainerType& GetDofs()
    {
        return mDofs;
    }

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


    /** adds a Dof to the node and return new added dof or existed one. */
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

    /** adds a Dof to the node and return new added dof or existed one. */
    inline typename DofType::Pointer pAddDof(DofType const& SourceDof)
    {
        KRATOS_TRY

        for(auto it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == SourceDof.GetVariable()){
                if((*it_dof)->GetReaction() != SourceDof.GetReaction())
                {
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

    /** adds a Dof to the node and return new added dof or existed one. */
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

    /** adds a Dof to the node and return new added dof or existed one. */
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

    /** adds a Dof to the node and return new added dof or existed one. */
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

    /** Return true if the dof of freedom is present on the node */
    inline bool HasDofFor(const VariableData& rDofVariable) const
    {
        for(auto it_dof = mDofs.begin() ; it_dof != mDofs.end() ; it_dof++){
            if((*it_dof)->GetVariable() == rDofVariable){
                return true;
            }
        }
        return false;
    }

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

    NodalData mNodalData;

    /** storage for the dof of the node */
    DofsContainerType  mDofs;

    /** A pointer to data related to this node. */
    DataValueContainer mData;


    ///Initial Position of the node
    PointType mInitialPosition;

    LockObject mNodeLock;

    ///@}
    ///@name Private Operators
    ///@{
    //*********************************************
    //this block is needed for refcounting
    mutable std::atomic<int> mReferenceCounter{0};

    friend void intrusive_ptr_add_ref(const NodeType* x)
    {
        x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
    }

    friend void intrusive_ptr_release(const NodeType* x)
    {
        if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
        std::atomic_thread_fence(std::memory_order_acquire);
        delete x;
        }
    }
    //*********************************************


    void SortDofs(){
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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Point );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
        rSerializer.save("NodalData", &mNodalData); // Storing it as pointer to be shared by Dof pointer
        rSerializer.save("Data", mData);
        rSerializer.save("Initial Position", mInitialPosition);
        rSerializer.save("Data", mDofs);

    }

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

//     namespace Globals
//     {
// 	extern Node DefaultNode3;
//     }

}  // namespace Kratos.
