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
//
//

#if !defined(KRATOS_DOF_H_INCLUDED )
#define  KRATOS_DOF_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>

// External includes
#include <boost/variant.hpp>

// Project includes
#include "includes/define.h"
#include "containers/data_value_container.h"
#include "containers/nodal_data.h"
#include "containers/vector_component_adaptor.h"
#include "containers/array_1d.h"


namespace Kratos
{


#define KRATOS_DOF_TRAITS \
        KRATOS_MAKE_DOF_TRAIT(0) Variable<TDataType> KRATOS_END_DOF_TRAIT(0); \
        KRATOS_MAKE_DOF_TRAIT(1) VariableComponent<VectorComponentAdaptor<array_1d<TDataType, 3> > > KRATOS_END_DOF_TRAIT(1); \
        KRATOS_MAKE_DOF_TRAIT(2) VariableComponent<VectorComponentAdaptor<array_1d<TDataType, 4> > > KRATOS_END_DOF_TRAIT(2); \
        KRATOS_MAKE_DOF_TRAIT(3) VariableComponent<VectorComponentAdaptor<array_1d<TDataType, 6> > > KRATOS_END_DOF_TRAIT(3); \
        KRATOS_MAKE_DOF_TRAIT(4) VariableComponent<VectorComponentAdaptor<array_1d<TDataType, 9> > > KRATOS_END_DOF_TRAIT(4);



template<class TDataType, class TVariableType = Variable<TDataType> >
struct DofTrait
{
    static const int Id;
};


#define KRATOS_MAKE_DOF_TRAIT(id) \
        template<class TDataType> \
        struct DofTrait<TDataType,


#define KRATOS_END_DOF_TRAIT(id) \
        >{\
        static const int Id = id;\
        }


KRATOS_DOF_TRAITS


#undef KRATOS_MAKE_DOF_TRAIT
#undef KRATOS_END_DOF_TRAIT


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

/// Dof represents a degree of freedom (DoF).
/** It is a lightweight object which holds its variable, like TEMPERATURE, its
state of freedom, and a reference to its value in the data structure.
This class enables the system to work with different set of dofs and also
represents the Dirichlet condition assigned to each dof.
*/
template<class TDataType>
class Dof{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Dof
    //KRATOS_CLASS_POINTER_DEFINITION(Dof);
    using Pointer=Dof*;

    typedef std::size_t IndexType;

    typedef std::size_t EquationIdType;

    typedef VariablesListDataValueContainer SolutionStepsDataContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor. This constructor takes all necessary
    informations to construct a degree of freedom. Also default
    values are used to make it easier to define for simple cases.

    @param rThisVariable Variable which this degree of freedom
    holds. This variable considered as unknown of problem to solved
    and fixing by Fix() method also applied to it. It must be a
    TDataType variable or component not a vector. For example
    DISPLACEMENT_X in structural element.

    @see Node
    @see Variable
    @see VariableComponent
    */
    template<class TVariableType>
    Dof(NodalData* pThisNodalData,
        const TVariableType& rThisVariable)
        : mIsFixed(false),
          mVariableType(DofTrait<TDataType, TVariableType>::Id),
          mReactionType(DofTrait<TDataType, Variable<TDataType> >::Id),
          mEquationId(IndexType()),
          mpNodalData(pThisNodalData)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(pThisNodalData->GetSolutionStepData().Has(rThisVariable))
            << "The Dof-Variable " << rThisVariable.Name() << " is not "
            << "in the list of variables" << std::endl;

        mIndex = mpNodalData->GetSolutionStepData().pGetVariablesList()->AddDof(&rThisVariable);
    }

    /** Constructor. This constructor takes the same input
    as the previous one, but add the reaction on the DoF
    declaration


    @param rThisVariable Variable which this degree of freedom
    holds. This variable considered as unknown of problem to solved
    and fixing by Fix() method also applied to it. It must be a
    TDataType variable or component not a vector. For example
    DISPLACEMENT_X in structural element.


    @param rThisReaction This is the right hand side variable in
    the system of equation correspounding to variable this dof
    holding. For example THERMAL_FLOW in thermal element. It will
    be none as default.


    @see Node
    @see Variable
    @see VariableComponent
    */
    template<class TVariableType, class TReactionType>
    Dof(NodalData* pThisNodalData,
        const TVariableType& rThisVariable,
        const TReactionType& rThisReaction)
        : mIsFixed(false),
          mVariableType(DofTrait<TDataType, TVariableType>::Id),
          mReactionType(DofTrait<TDataType, TReactionType>::Id),
          mEquationId(IndexType()),
          mpNodalData(pThisNodalData)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(pThisNodalData->GetSolutionStepData().Has(rThisVariable))
            << "The Dof-Variable " << rThisVariable.Name() << " is not "
            << "in the list of variables" << std::endl;

        KRATOS_DEBUG_ERROR_IF_NOT(pThisNodalData->GetSolutionStepData().Has(rThisReaction))
            << "The Reaction-Variable " << rThisReaction.Name() << " is not "
            << "in the list of variables" << std::endl;
        
        mIndex = mpNodalData->GetSolutionStepData().pGetVariablesList()->AddDof(&rThisVariable, &rThisReaction);

    }

    //This default constructor is needed for serializer
    Dof()
        : mIsFixed(false),
          mVariableType(DofTrait<TDataType, Variable<TDataType> >::Id),
          mReactionType(DofTrait<TDataType, Variable<TDataType> >::Id),
          mIndex(),
          mEquationId(IndexType()),
          mpNodalData()
    {
    }

    /// Copy constructor.
    Dof(Dof const& rOther)
        : mIsFixed(rOther.mIsFixed),
          mVariableType(rOther.mVariableType),
          mReactionType(rOther.mReactionType),
          mIndex(rOther.mIndex),
          mEquationId(rOther.mEquationId),
          mpNodalData(rOther.mpNodalData)
    {
    }


    /// Destructor.
    ~Dof() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Dof& operator=(Dof const& rOther)
    {
        mIsFixed = rOther.mIsFixed;
        mEquationId = rOther.mEquationId;
        mpNodalData = rOther.mpNodalData;
        mIndex = rOther.mIndex;
        mVariableType = rOther.mVariableType;
        mReactionType = rOther.mReactionType;
        // mData = rOther.mData;

        return *this;
    }

    template<class TVariableType>
    typename TVariableType::Type& operator()(const TVariableType& rThisVariable, IndexType SolutionStepIndex = 0)
    {
        return GetSolutionStepValue(rThisVariable, SolutionStepIndex);
    }

    template<class TVariableType>
    typename TVariableType::Type const& operator()(const TVariableType& rThisVariable, IndexType SolutionStepIndex = 0) const
    {
        return GetSolutionStepValue(rThisVariable, SolutionStepIndex);
    }

    TDataType& operator()(IndexType SolutionStepIndex = 0)
    {
        return GetSolutionStepValue(SolutionStepIndex);
    }

    TDataType const& operator()(IndexType SolutionStepIndex = 0) const
    {
        return GetSolutionStepValue(SolutionStepIndex);
    }

    TDataType& operator[](IndexType SolutionStepIndex)
    {
        return GetSolutionStepValue(SolutionStepIndex);
    }

    TDataType const& operator[](IndexType SolutionStepIndex) const
    {
        return GetSolutionStepValue(SolutionStepIndex);
    }


    ///@}
    ///@name Operations
    ///@{



    TDataType& GetSolutionStepValue(IndexType SolutionStepIndex = 0)
    {
        return GetReference(GetVariable(), mpNodalData->GetSolutionStepData(), SolutionStepIndex, mVariableType);
    }


    TDataType const& GetSolutionStepValue(IndexType SolutionStepIndex = 0) const
    {
        return GetReference(GetVariable(), mpNodalData->GetSolutionStepData(), SolutionStepIndex, mVariableType);
    }


    template<class TVariableType>
    typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex = 0)
    {
        return mpNodalData->GetSolutionStepData().GetValue(rThisVariable, SolutionStepIndex);
    }


    template<class TVariableType>
    typename TVariableType::Type const& GetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex = 0) const
    {
        return mpNodalData->GetSolutionStepData().GetValue(rThisVariable, SolutionStepIndex);
    }


    TDataType& GetSolutionStepReactionValue(IndexType SolutionStepIndex = 0)
    {
        return GetReference(GetReaction(), mpNodalData->GetSolutionStepData(), SolutionStepIndex, mReactionType);
    }

    TDataType const& GetSolutionStepReactionValue(IndexType SolutionStepIndex = 0) const
    {
        return GetReference(GetReaction(), mpNodalData->GetSolutionStepData(), SolutionStepIndex, mReactionType);
    }


    ///@}
    ///@name Access
    ///@{

    IndexType Id() const
    {
        return mpNodalData->GetId();
    }

    IndexType GetId() const
    {
        return mpNodalData->GetId();
    }

    /** Returns variable assigned to this degree of freedom. */
    const VariableData& GetVariable() const
    {
        return mpNodalData->GetSolutionStepData().pGetVariablesList()->GetDofVariable(mIndex);
    }

    /** Returns reaction variable of this degree of freedom. */
    const VariableData& GetReaction() const
    {
        auto p_reaction = mpNodalData->GetSolutionStepData().pGetVariablesList()->pGetDofReaction(mIndex);
        return (p_reaction == nullptr) ? msNone : *p_reaction;
    }

    template<class TReactionType>
    void SetReaction(TReactionType const& rReaction)
    {
        mReactionType = DofTrait<TDataType, TReactionType>::Id;
        mpNodalData->GetSolutionStepData().pGetVariablesList()->SetDofReaction(&rReaction, mIndex);
    }

    /** Return the Equation Id related to this degree eof freedom.
     */
    EquationIdType EquationId() const
    {
        return mEquationId;
    }

    /** Sets the Equation Id to the desired value
     */
    void SetEquationId(EquationIdType NewEquationId)
    {
        mEquationId = NewEquationId;
    }

    /** Fixes the Dof
     */
    void FixDof()
    {
        mIsFixed=true;
    }


    /** Frees the degree of freedom
     */
    void FreeDof()
    {
        mIsFixed=false;
    }


    SolutionStepsDataContainerType* GetSolutionStepsData()
    {
        return &(mpNodalData->GetSolutionStepData());
    }

    void SetNodalData(NodalData* pNewNodalData)
    {
        auto p_variable = &GetVariable();
        auto p_reaction = mpNodalData->GetSolutionStepData().pGetVariablesList()->pGetDofReaction(mIndex);
        mpNodalData = pNewNodalData;
        if(p_reaction != nullptr){
            mIndex = mpNodalData->GetSolutionStepData().pGetVariablesList()->AddDof(p_variable, p_reaction);
        }
        else{
            mIndex = mpNodalData->GetSolutionStepData().pGetVariablesList()->AddDof(p_variable);
        }
    }

    bool HasReaction() const
    {
        return (mpNodalData->GetSolutionStepData().pGetVariablesList()->pGetDofReaction(mIndex) != nullptr);
    }

    ///@}
    ///@name Inquiry
    ///@{

    bool IsFixed() const
    {
        return mIsFixed;
    }


    bool IsFree() const
    {
        return !IsFixed();
    }

    ///@}
    ///@name Input and output
    ///@{


    /// Turn back information as a string.
    std::string Info() const 
    {
        std::stringstream buffer;


        if(IsFixed())
            buffer << "Fix " << GetVariable().Name() << " degree of freedom";
        else
            buffer << "Free " << GetVariable().Name() << " degree of freedom";


        return buffer.str();
    }



    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const 
    {
        rOStream << Info();
    }


    /// Print object's data.
    void PrintData(std::ostream& rOStream) const 
    {
        rOStream << "    Variable               : " << GetVariable().Name() << std::endl;
        rOStream << "    Reaction               : " << GetReaction().Name() << std::endl;
        if(IsFixed())
            rOStream << "    IsFixed                : True" << std::endl;
        else
            rOStream << "    IsFixed                : False" << std::endl;
        rOStream << "    Equation Id            : " << mEquationId << std::endl;
    }




    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{

    static const Variable<TDataType> msNone;
    static constexpr int msIsFixedPosition = 63;


    ///@}
    ///@name Member Variables
    ///@{


    /** True is is fixed */
    int mIsFixed : 1;

    int mVariableType : 4;

    int mReactionType : 4;

    int mIndex : 6;

    /** Equation identificator of the degree of freedom */
#ifdef KRATOS_ENV32BIT // Required to avoid overflow on 32 bit systems
    EquationIdType mEquationId : 32;
#else
    EquationIdType mEquationId : 48;
#endif

    /** A pointer to nodal data stored in node which is corresponded to this dof */
    NodalData* mpNodalData;

    ///@}
    ///@name Private Operators
    ///@{
#define KRATOS_MAKE_DOF_TRAIT(id) \
        case id : \
                return rData.GetValue(static_cast<


#define KRATOS_END_DOF_TRAIT(id) \
 const&>(ThisVariable), SolutionStepIndex);

    TDataType& GetReference(VariableData const& ThisVariable, VariablesListDataValueContainer& rData, IndexType SolutionStepIndex, int ThisId)
    {
        switch(ThisId)
        {
            KRATOS_DOF_TRAITS
        }
        KRATOS_ERROR << "Not supported type for Dof" << std::endl;
    }

    TDataType const& GetReference(VariableData const& ThisVariable, VariablesListDataValueContainer const& rData, IndexType SolutionStepIndex, int ThisId) const
    {
        switch(ThisId)
        {
            KRATOS_DOF_TRAITS
        }
        KRATOS_ERROR << "Not supported type for Dof" << std::endl;
    }

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const 
    {
        rSerializer.save("IsFixed", static_cast<bool>(mIsFixed));
        rSerializer.save("EquationId", static_cast<EquationIdType>(mEquationId));
        rSerializer.save("NodalData", mpNodalData);
        rSerializer.save("VariableType", static_cast<int>(mVariableType));
        rSerializer.save("ReactionType", static_cast<int>(mReactionType));
        rSerializer.save("Index", static_cast<int>(mIndex));

    }

    void load(Serializer& rSerializer) 
    {
        std::string name;
        bool is_fixed;
        rSerializer.load("IsFixed", is_fixed);
        mIsFixed=is_fixed;
        EquationIdType equation_id;
        rSerializer.load("EquationId", equation_id);
        mEquationId = equation_id;
        rSerializer.load("NodalData", mpNodalData);

        int variable_type;
        int reaction_type;
        rSerializer.load("VariableType", variable_type);
        rSerializer.load("ReactionType", reaction_type);

        mVariableType = variable_type;
        mReactionType = reaction_type;

        int index;
        rSerializer.load("Index", index);
        mIndex = index;
    }
    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Default constructor.
    //      Dof();



    ///@}

}; // Class Dof
template<class TDataType> const Variable<TDataType> Dof<TDataType>::msNone("NONE");


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TDataType>
inline std::istream& operator >> (std::istream& rIStream,
                                  Dof<TDataType>& rThis);


/// output stream function
template<class TDataType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Dof<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);


    return rOStream;
}
///@}
///@name Operations
///@{


/// Greater than operator
template<class TDataType>
inline bool operator > ( Dof<TDataType> const& First,
                         Dof<TDataType> const& Second)
{
    if(First.Id() == Second.Id())
        return (First.GetVariable().Key() > Second.GetVariable().Key());

    return (First.Id() > Second.Id());
}

/// Less than operator
template<class TDataType>
inline bool operator < ( Dof<TDataType> const& First,
                         Dof<TDataType> const& Second)
{
    if(First.Id() == Second.Id())
        return (First.GetVariable().Key() < Second.GetVariable().Key());

    return (First.Id() < Second.Id());
}

/// Greater equal operator
template<class TDataType>
inline bool operator >= ( Dof<TDataType> const& First,
                          Dof<TDataType> const& Second)
{
    if(First.Id() == Second.Id())
        return (First.GetVariable().Key() >= Second.GetVariable().Key());

    return (First.Id() > Second.Id());
}

/// Less equal operator
template<class TDataType>
inline bool operator <= ( Dof<TDataType> const& First,
                          Dof<TDataType> const& Second)
{
    if(First.Id() == Second.Id())
        return (First.GetVariable().Key() <= Second.GetVariable().Key());

    return (First.Id() < Second.Id());
}

/// Equal operator
template<class TDataType>
inline bool operator == ( Dof<TDataType> const& First,
                          Dof<TDataType> const& Second)
{
    return ((First.Id() == Second.Id()) && (First.GetVariable().Key() == Second.GetVariable().Key()));
}

///@}


}  // namespace Kratos.


#undef KRATOS_DOF_TRAITS
#undef KRATOS_MAKE_DOF_TRAIT
#undef KRATOS_END_DOF_TRAIT


#endif // KRATOS_DOF_H_INCLUDED  defined













