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


#if !defined(KRATOS_NODE_H_INCLUDED )
#define  KRATOS_NODE_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>


// External includes


// Project includes
#include "includes/define.h"
#include "geometries/point.h"
#include "includes/dof.h"
#include "containers/vector_map.h"
#include "containers/pointer_vector_set.h"
#include "containers/variables_list_data_value_container.h"
#include "utilities/indexed_object.h"
#include "containers/flags.h"

#include "containers/weak_pointer_vector.h"

#ifdef _OPENMP
#include "omp.h"
#endif


namespace Kratos
{

class Element;

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
template<std::size_t TDimension, class TDofType = Dof<double> >
class Node : public Point,  public IndexedObject, public Flags
{
    class GetDofKey : public std::unary_function<TDofType, VariableData::KeyType>
    {
    public:
        VariableData::KeyType operator()(TDofType  const & This)
        {
            return This.GetVariable().Key();
        }
    };

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Node
    KRATOS_CLASS_POINTER_DEFINITION(Node);

    typedef Node<TDimension, TDofType> NodeType;

    typedef Point BaseType;

    typedef Point PointType;

    typedef TDofType DofType;

    typedef std::size_t IndexType;

    typedef typename std::size_t SizeType;

    typedef PointerVectorSet<TDofType, GetDofKey> DofsContainerType;

    typedef VariablesListDataValueContainer SolutionStepsNodalDataContainerType;

    typedef VariablesListDataValueContainer::BlockType BlockType;

    typedef Variable<double> DoubleVariableType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Node()
        : BaseType()
        , IndexedObject(0)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition()
#ifdef _OPENMP
        , mNodeLock()
#endif
    {

        CreateSolutionStepData();
// 	mDofs.SetMaxBufferSize(0);

#ifdef _OPENMP
        omp_init_lock(&mNodeLock);
#endif



    }

    Node(IndexType NewId )
        : BaseType()
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition()
#ifdef _OPENMP
        , mNodeLock()
#endif
    {
        KRATOS_ERROR <<  "Calling the default constructor for the node ... illegal operation!!" << std::endl;
        CreateSolutionStepData();

    }

    /// 1d constructor.
    Node(IndexType NewId, double const& NewX)
        : BaseType(NewX)
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition(NewX)
#ifdef _OPENMP
        , mNodeLock()
#endif
    {
#ifdef _OPENMP
        omp_init_lock(&mNodeLock);
#endif
        CreateSolutionStepData();

    }

    /// 2d constructor.
    Node(IndexType NewId, double const& NewX, double const& NewY)
        : BaseType(NewX, NewY)
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition(NewX, NewY)
#ifdef _OPENMP
        , mNodeLock()
#endif
    {
#ifdef _OPENMP
        omp_init_lock(&mNodeLock);
#endif
        CreateSolutionStepData();

    }

    /// 3d constructor.
    Node(IndexType NewId, double const& NewX, double const& NewY, double const& NewZ)
        : BaseType(NewX, NewY, NewZ)
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition(NewX, NewY, NewZ)
#ifdef _OPENMP
        , mNodeLock()
#endif
    {
        CreateSolutionStepData();
// 	mDofs.SetMaxBufferSize(0);

#ifdef _OPENMP
        omp_init_lock(&mNodeLock);
#endif


    }

    /// Point constructor.
    Node(IndexType NewId, PointType const& rThisPoint)
        : BaseType(rThisPoint)
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition(rThisPoint)
#ifdef _OPENMP
        , mNodeLock()
#endif
    {

        CreateSolutionStepData();
// 	mDofs.SetMaxBufferSize(0);

#ifdef _OPENMP
        omp_init_lock(&mNodeLock);
#endif


    }

    /** Copy constructor. Initialize this node with given node.*/
    Node(Node const& rOtherNode) = delete;

    /** Copy constructor from a node with different dimension.*/
    template<SizeType TOtherDimension>
    Node(Node<TOtherDimension> const& rOtherNode) = delete;


    /**
     * Constructor using coordinates stored in given array. Initialize
    this point with the coordinates in the array. */
    template<class TVectorType>
    Node(IndexType NewId, vector_expression<TVectorType> const&  rOtherCoordinates)
        : BaseType(rOtherCoordinates)
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition(rOtherCoordinates)
#ifdef _OPENMP
        , mNodeLock()
#endif
    {

        CreateSolutionStepData();
// 	mDofs.SetMaxBufferSize(0);

#ifdef _OPENMP
        omp_init_lock(&mNodeLock);
#endif

    }



    /** Constructor using coordinates stored in given std::vector. Initialize
    this point with the coordinates in the array. */
    Node(IndexType NewId, std::vector<double> const&  rOtherCoordinates)
        : BaseType(rOtherCoordinates)
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData()
        , mInitialPosition()
#ifdef _OPENMP
        , mNodeLock()
#endif
    {
        CreateSolutionStepData();
// 	mDofs.SetMaxBufferSize(0);

#ifdef _OPENMP
        omp_init_lock(&mNodeLock);
#endif


    }

    /// 3d with variables list and data constructor.
    Node(IndexType NewId, double const& NewX, double const& NewY, double const& NewZ, VariablesList*  pVariablesList, BlockType const * ThisData, SizeType NewQueueSize = 1)
        : BaseType(NewX, NewY, NewZ)
        , IndexedObject(NewId)
        , Flags()
        , mDofs()
        , mData()
        , mSolutionStepsNodalData(pVariablesList,ThisData,NewQueueSize)
        , mInitialPosition(NewX, NewY, NewZ)
#ifdef _OPENMP
        , mNodeLock()
#endif
    {
// 	mDofs.SetMaxBufferSize(0);
#ifdef _OPENMP
        omp_init_lock(&mNodeLock);
#endif
    }

    typename Node<TDimension>::Pointer Clone()
    {
        Node<3>::Pointer p_new_node = Kratos::make_shared<Node<3> >( this->Id(), (*this)[0], (*this)[1], (*this)[2]);
        p_new_node->mSolutionStepsNodalData = this->mSolutionStepsNodalData;

        Node<3>::DofsContainerType& my_dofs = (this)->GetDofs();
        for (typename DofsContainerType::const_iterator it_dof = my_dofs.begin(); it_dof != my_dofs.end(); it_dof++)
        {
            p_new_node->pAddDof(*it_dof);
        }

        p_new_node->mData = this->mData;
        p_new_node->mInitialPosition = this->mInitialPosition;

        p_new_node->Set(Flags(*this));
        //KRATOS_ERROR << "Must implement correctly the copy of the flags" << std::endl;
        return p_new_node;
    }

    /// Destructor.
    ~Node() override
    {
#ifdef _OPENMP
        omp_destroy_lock(&mNodeLock);
#endif
    }

    void SetId(IndexType NewId) override
    {
        IndexedObject::SetId(NewId);
        Node<3>::DofsContainerType& my_dofs = (this)->GetDofs();
        for(Node<3>::DofsContainerType::iterator it_dof = my_dofs.begin(); it_dof != my_dofs.end(); it_dof++)
        {
            it_dof->SetId(NewId);
        }
    }

#ifdef _OPENMP
    omp_lock_t& GetLock()
    {
        return mNodeLock;
    }
#endif

    inline void SetLock()
    {
        //does nothing if openMP is not present
#ifdef _OPENMP
        omp_set_lock(&mNodeLock);
#endif
    }

    inline void UnSetLock()
    {
        //does nothing if openMP is not present
#ifdef _OPENMP
        omp_unset_lock(&mNodeLock);
#endif
    }

    ///@}
    ///@name Operators
    ///@{

    Node& operator=(const Node& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator =(rOther);
        IndexedObject::operator=(rOther);

        // Deep copying the dofs
        for(typename DofsContainerType::const_iterator it_dof = rOther.mDofs.begin() ; it_dof != rOther.mDofs.end() ; it_dof++)
        {
            pAddDof(*it_dof);
        }

        mData = rOther.mData;
        mSolutionStepsNodalData = rOther.mSolutionStepsNodalData;
        mInitialPosition = rOther.mInitialPosition;

        return *this;
    }

    /// Assignment operator.
    template<SizeType TOtherDimension>
    Node& operator=(const Node<TOtherDimension>& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator =(rOther);
        IndexedObject::operator=(rOther);
        for(typename DofsContainerType::const_iterator it_dof = rOther.mDofs.begin() ; it_dof != rOther.mDofs.end() ; it_dof++)
        {
            pAddDof(*it_dof);
        }

        mData = rOther.mData;
        mSolutionStepsNodalData = rOther.mSolutionStepsNodalData;
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

    template<class TAdaptorType> typename TAdaptorType::Type& operator[](const VariableComponent<TAdaptorType>& rThisVariable)
    {
        return GetValue(rThisVariable);
    }

    template<class TAdaptorType> const typename TAdaptorType::Type& operator[](const VariableComponent<TAdaptorType>& rThisVariable) const
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
        mSolutionStepsNodalData.PushFront();
    }

    void CloneSolutionStepData()
    {
        mSolutionStepsNodalData.CloneFront();
    }

    void OverwriteSolutionStepData(IndexType SourceSolutionStepIndex, IndexType DestinationSourceSolutionStepIndex)
    {
        mSolutionStepsNodalData.AssignData(mSolutionStepsNodalData.Data(SourceSolutionStepIndex), DestinationSourceSolutionStepIndex);
    }

    void ClearSolutionStepsData()
    {
        mSolutionStepsNodalData.Clear();
    }

    void SetSolutionStepVariablesList(VariablesList* pVariablesList)
    {
        mSolutionStepsNodalData.SetVariablesList(pVariablesList);
    }

    VariablesListDataValueContainer& SolutionStepData()
    {
        return mSolutionStepsNodalData;
    }

    const VariablesListDataValueContainer& SolutionStepData() const
    {
        return mSolutionStepsNodalData;
    }

    DataValueContainer& Data()
    {
        return mData;
    }

    const DataValueContainer& Data() const
    {
        return mData;
    }

    template<class TVariableType> typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable)
    {
        return mSolutionStepsNodalData.GetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type const& GetSolutionStepValue(const TVariableType& rThisVariable) const
    {
        return mSolutionStepsNodalData.GetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
    {
        return mSolutionStepsNodalData.GetValue(rThisVariable, SolutionStepIndex);
    }

    template<class TVariableType> typename TVariableType::Type const& GetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex) const
    {
        return mSolutionStepsNodalData.GetValue(rThisVariable, SolutionStepIndex);
    }


    template<class TDataType> bool SolutionStepsDataHas(const Variable<TDataType>& rThisVariable) const
    {
        return mSolutionStepsNodalData.Has(rThisVariable);
    }
    template<class TAdaptorType> bool SolutionStepsDataHas(const VariableComponent<TAdaptorType>& rThisVariable) const
    {
        return mSolutionStepsNodalData.Has(rThisVariable);
    }

    //*******************************************************************************************
    //By Riccardo
    //very similar to the one before BUT throws an error if the variable does not exist
    template<class TVariableType> typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable)
    {
        return mSolutionStepsNodalData.FastGetValue(rThisVariable);
    }

    template<class TVariableType> const typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable) const
    {
        return mSolutionStepsNodalData.FastGetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex)
    {
        return mSolutionStepsNodalData.FastGetValue(rThisVariable, SolutionStepIndex);
    }

    template<class TVariableType> const typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex) const
    {
        return mSolutionStepsNodalData.FastGetValue(rThisVariable, SolutionStepIndex);
    }

    template<class TVariableType> typename TVariableType::Type& FastGetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex, IndexType ThisPosition)
    {
        return mSolutionStepsNodalData.FastGetValue(rThisVariable, SolutionStepIndex, ThisPosition);
    }

    template<class TVariableType> typename TVariableType::Type& FastGetCurrentSolutionStepValue(const TVariableType& rThisVariable, IndexType ThisPosition)
    {
        return mSolutionStepsNodalData.FastGetCurrentValue(rThisVariable, ThisPosition);
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
            return mSolutionStepsNodalData.GetValue(rThisVariable, SolutionStepIndex);

        return mData.GetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type const& GetValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex) const
    {
        if(!mData.Has(rThisVariable))
            return mSolutionStepsNodalData.GetValue(rThisVariable, SolutionStepIndex);

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
    template<class TAdaptorType> bool Has(const VariableComponent<TAdaptorType>& rThisVariable) const
    {
        return mData.Has(rThisVariable);
    }

    ///@}
    ///@name Operations
    ///@{

    template<class TVariableType>
    inline void Fix(const TVariableType& rDofVariable)
    {
        typename DofsContainerType::iterator it_dof = mDofs.find(rDofVariable);
        if(it_dof != mDofs.end())
        {
            it_dof->FixDof();
        }
        else
        {
#ifdef KRATOS_DEBUG
            if(OpenMPUtils::IsInParallel() != 0)
            {
                KRATOS_ERROR << "Attempting to Fix the variable: " << rDofVariable << " within a parallel region. This is not permitted. Create the Dof first by pAddDof" << std::endl;
            }
#endif
            pAddDof(rDofVariable)->FixDof();
        }
    }

    template<class TVariableType>
    inline void Free(const TVariableType& rDofVariable)
    {
        typename DofsContainerType::iterator it_dof = mDofs.find(rDofVariable);
        if(it_dof != mDofs.end())
        {
            it_dof->FreeDof();
        }
        else
        {
#ifdef KRATOS_DEBUG
            if(OpenMPUtils::IsInParallel() != 0)
            {
                KRATOS_ERROR << "Attempting to Fix the variable: " << rDofVariable << " within a parallel region. This is not permitted. Create the Dof first by pAddDof" << std::endl;
            }
#endif
            pAddDof(rDofVariable)->FreeDof();
        }
    }

    IndexType GetBufferSize() const
    {
        return mSolutionStepsNodalData.QueueSize();
    }

    void SetBufferSize(IndexType NewBufferSize)
    {
        mSolutionStepsNodalData.Resize(NewBufferSize);
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

    VariablesList * pGetVariablesList()
    {
        return mSolutionStepsNodalData.pGetVariablesList();
    }

    const VariablesList * pGetVariablesList() const
    {
        return mSolutionStepsNodalData.pGetVariablesList();
    }

    ///@}
    ///@name Dofs
    ///@{

    //advanced functions by Riccardo
    template<class TVariableType>
    inline unsigned int GetDofPosition(TVariableType const& rDofVariable) const
    {
        typename DofsContainerType::const_iterator it=mDofs.find(rDofVariable.Key());
        return it - mDofs.begin();
    }

    template<class TVariableType>
    inline DofType& GetDof(TVariableType const& rDofVariable, int pos)
    {
        typename DofsContainerType::iterator it_begin = mDofs.begin();
        typename DofsContainerType::iterator it_end = mDofs.end();
        typename DofsContainerType::iterator it;
        // If the guess is exact return the guess
        if(pos < it_end-it_begin) {
            it = it_begin + pos;
            if( (it)->GetVariable() == rDofVariable) {
                return *it;
            }
        }

        // Otherwise do a find
        return GetDof(rDofVariable);
    }

    /** returns the Dof asociated with variable  */
    template<class TVariableType>
    inline DofType& GetDof(TVariableType const& rDofVariable)
    {
        typename DofsContainerType::iterator it=mDofs.find(rDofVariable.Key());
        if ( it!= mDofs.end() ) {
            return *it;
        }

        KRATOS_ERROR <<  "Not existant DOF in node #" << Id() << " for variable : "
            << rDofVariable.Name() << std::endl;
    }

    /** returns all of the Dofs  */
    DofsContainerType& GetDofs()
    {
        return mDofs;
    }

    /** returns a counted pointer to the Dof asociated with variable  */
    template<class TVariableType>
    inline typename DofType::Pointer pGetDof(TVariableType const& rDofVariable)
    {
        typename DofsContainerType::iterator it=mDofs.find(rDofVariable.Key());
        if ( it!= mDofs.end() ) {
            return *(it.base());
        }

        KRATOS_ERROR <<  "Not existant DOF in node #" << Id() << " for variable : "
            << rDofVariable.Name() << std::endl;
    }

    /** adds a Dof to the node and return new added dof or existed one. */
    template<class TVariableType>
    inline typename DofType::Pointer pAddDof(TVariableType const& rDofVariable)
    {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF(rDofVariable.Key() == 0) << "Variable  " << rDofVariable
            << " has key zero key when adding Dof for node " << this->Id() << std::endl;

        typename DofsContainerType::iterator it_dof = mDofs.find(rDofVariable);
        if(it_dof != mDofs.end())
        {
            return *(it_dof.base());
        }

        typename DofType::Pointer p_new_dof =  Kratos::make_shared<DofType>(Id(), &mSolutionStepsNodalData, rDofVariable);
        mDofs.insert(mDofs.begin(), p_new_dof);

//         if(!mDofs.IsSorted())
        mDofs.Sort();

        return p_new_dof;

        KRATOS_CATCH(*this);
    }

    /** adds a Dof to the node and return new added dof or existed one. */
    inline typename DofType::Pointer pAddDof(DofType const& SourceDof)
    {
        KRATOS_TRY

        typename DofsContainerType::iterator it_dof = mDofs.find(SourceDof.GetVariable());
        if(it_dof != mDofs.end())
        {
            if(it_dof->GetReaction() != SourceDof.GetReaction())
            {
                *it_dof = SourceDof;
                it_dof->SetId(Id());
                it_dof->SetSolutionStepsData(&mSolutionStepsNodalData);
            }
            return *(it_dof.base());
        }

        typename DofType::Pointer p_new_dof =  Kratos::make_shared<DofType>(SourceDof);
        mDofs.insert(mDofs.begin(), p_new_dof);

        p_new_dof->SetId(Id());

        p_new_dof->SetSolutionStepsData(&mSolutionStepsNodalData);

        mDofs.Sort();

        return p_new_dof;

        KRATOS_CATCH(*this);
    }

    /** adds a Dof to the node and return new added dof or existed one. */
    template<class TVariableType, class TReactionType>
    inline typename DofType::Pointer pAddDof(TVariableType const& rDofVariable, TReactionType const& rDofReaction)
    {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF_NOT(rDofVariable.Key() == 0) << "Variable  " << rDofVariable
            << " has key zero key when adding Dof for node " << this->Id() << std::endl;
        KRATOS_DEBUG_ERROR_IF_NOT(rDofReaction.Key() == 0) << "Reaction  " << rDofReaction
            << " has key zero when adding reactions for node " << this->Id() << std::endl;

        typename DofsContainerType::iterator it_dof = mDofs.find(rDofVariable);
        if(it_dof != mDofs.end())
        {
            it_dof->SetReaction(rDofReaction);
            return *(it_dof.base());
        }

        typename DofType::Pointer p_new_dof =  Kratos::make_shared<DofType>(Id(), &mSolutionStepsNodalData, rDofVariable, rDofReaction);
        mDofs.insert(mDofs.begin(), p_new_dof);

//         if(!mDofs.IsSorted())
        mDofs.Sort();

        return p_new_dof;

        KRATOS_CATCH(*this);

    }

    /** adds a Dof to the node and return new added dof or existed one. */
    template<class TVariableType>
    inline DofType& AddDof(TVariableType const& rDofVariable)
    {
        KRATOS_TRY

#ifdef KRATOS_DEBUG
        if(rDofVariable.Key() == 0)
        {
            KRATOS_ERROR << "Variable  " << rDofVariable << " has key zero key when adding Dof for node " << this->Id() << std::endl;
        }
#endif

        typename DofsContainerType::iterator it_dof = mDofs.find(rDofVariable);
        if(it_dof != mDofs.end())
        {
            return *it_dof;
        }

        typename DofType::Pointer p_new_dof =  Kratos::make_shared<DofType>(Id(), &mSolutionStepsNodalData, rDofVariable);
        mDofs.insert(mDofs.begin(), p_new_dof);

//         if(!mDofs.IsSorted())
        mDofs.Sort();

        return *p_new_dof;

        KRATOS_CATCH(*this);

    }

    /** adds a Dof to the node and return new added dof or existed one. */
    template<class TVariableType, class TReactionType>
    inline DofType& AddDof(TVariableType const& rDofVariable, TReactionType const& rDofReaction)
    {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF(rDofVariable.Key() == 0) << "Variable  " << rDofVariable
            << " has key zero key when adding Dof for node " << this->Id() << std::endl;
        KRATOS_DEBUG_ERROR_IF(rDofReaction.Key() == 0) << "Reaction  " << rDofReaction
            << " has key zero when adding reactions for node " << this->Id() << std::endl;

        typename DofsContainerType::iterator it_dof = mDofs.find(rDofVariable);
        if(it_dof != mDofs.end())
        {
            it_dof->SetReaction(rDofReaction);
            return *it_dof;
        }

        typename DofType::Pointer p_new_dof =  Kratos::make_shared<DofType>(Id(), &mSolutionStepsNodalData, rDofVariable, rDofReaction);
        mDofs.insert(mDofs.begin(), p_new_dof);

//         if(!mDofs.IsSorted())
        mDofs.Sort();

        return *p_new_dof;

        KRATOS_CATCH(*this);

    }

    ///@}
    ///@name Inquiry
    ///@{

    /** Return true if the dof of freedom is present on the node */
    inline bool HasDofFor(const VariableData& rDofVariable) const
    {
        return (mDofs.find(rDofVariable) != mDofs.end());
    }

    inline  bool IsFixed(const VariableData& rDofVariable) const
    {
        typename DofsContainerType::const_iterator i;
        return (((i= mDofs.find(rDofVariable)) == mDofs.end()) ? false : i->IsFixed());
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
            rOStream << "        " << i->Info() << std::endl;
// 	  rOStream << "        " << "solution steps  : " << *mSolutionStepsNodalData;
// 	  rOStream << "        " << "solution steps capacity : " << mSolutionStepsNodalData.GetContainer().capacity();
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

    /** storage for the dof of the node */
    DofsContainerType  mDofs;

    /** A pointer to data related to this node. */
    DataValueContainer mData;

    SolutionStepsNodalDataContainerType mSolutionStepsNodalData;

    ///Initial Position of the node
    PointType mInitialPosition;

#ifdef _OPENMP
    omp_lock_t mNodeLock;
#endif

    ///@}
    ///@name Private Operators
    ///@{



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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
        rSerializer.save("Data", mData);
        const SolutionStepsNodalDataContainerType* pSolutionStepsNodalData = &mSolutionStepsNodalData;
        // I'm saving it as pointer so the dofs pointers will point to it as stored pointer. Pooyan.
        rSerializer.save("Solution Steps Nodal Data", pSolutionStepsNodalData);
        rSerializer.save("Initial Position", mInitialPosition);
        rSerializer.save("Data", mDofs);

    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Point );
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject );
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
        rSerializer.load("Data", mData);
        SolutionStepsNodalDataContainerType* pSolutionStepsNodalData = &mSolutionStepsNodalData;
        rSerializer.load("Solution Steps Nodal Data", pSolutionStepsNodalData);
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

// template class KRATOS_API(KRATOS_CORE) KratosComponents<Node<3,double> >;
// template class KRATOS_API(KRATOS_CORE) KratosComponents<Node<3,double>::Pointer >;

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<std::size_t TDimension, class TDofType>
inline std::istream& operator >> (std::istream& rIStream,
                                  Node<TDimension, TDofType>& rThis);

/// output stream function
template<std::size_t TDimension, class TDofType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Node<TDimension, TDofType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : ";
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

//*********************************************************************************
//*********************************************************************************
//*********************************************************************************
//definition of the NEIGHBOUR_NODES variable
//*********************************************************************************
//*********************************************************************************
//*********************************************************************************

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_API

KRATOS_DEFINE_VARIABLE(WeakPointerVector<Node<3> >, NEIGHBOUR_NODES)
KRATOS_DEFINE_VARIABLE(WeakPointerVector<Node<3> >, FATHER_NODES)

#undef  KRATOS_EXPORT_MACRO
#define KRATOS_EXPORT_MACRO KRATOS_NO_EXPORT


//     namespace Globals
//     {
// 	extern Node<3> DefaultNode3;
//     }



}  // namespace Kratos.

#endif // KRATOS_NODE_H_INCLUDED  defined
