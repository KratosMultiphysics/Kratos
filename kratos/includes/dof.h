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
//
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>
#include <bitset>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/data_value_container.h"
#include "containers/nodal_data.h"
#include "containers/flags.h"
#include "includes/kratos_flags.h"

namespace Kratos
{

#define KRATOS_DOF_TRAITS \
        KRATOS_MAKE_DOF_TRAIT(0) Variable<TDataType> KRATOS_END_DOF_TRAIT(0);

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
    information to construct a degree of freedom. Also default
    values are used to make it easier to define for simple cases.

    @param rThisVariable Variable which this degree of freedom
    holds. This variable considered as unknown of problem to solved
    and fixing by Fix() method also applied to it. It must be a
    TDataType variable or component not a vector. For example
    DISPLACEMENT_X in structural element.

    @see Node
    @see Variable
    */
    template<class TVariableType>
    Dof(NodalData* pThisNodalData,
        const TVariableType& rThisVariable)
        : mIsFixed(false),
          mIsActive(false),
          mVariableType(DofTrait<TDataType, TVariableType>::Id),
          mReactionType(DofTrait<TDataType, Variable<TDataType> >::Id),
          mEquationId(IndexType()),
          mpNodalData(pThisNodalData)
    {
        KRATOS_DEBUG_ERROR_IF_NOT(pThisNodalData->GetSolutionStepData().Has(rThisVariable))
            << "The Dof-Variable " << rThisVariable.Name() << " is not "
            << "in the list of variables" << std::endl;

        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        mIndex = p_nodal_data->GetSolutionStepData().pGetVariablesList()->AddDof(&rThisVariable);
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
    */
    template<class TVariableType, class TReactionType>
    Dof(NodalData* pThisNodalData,
        const TVariableType& rThisVariable,
        const TReactionType& rThisReaction)
        : mIsFixed(false),
          mIsActive(false),
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

        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        mIndex = p_nodal_data->GetSolutionStepData().pGetVariablesList()->AddDof(&rThisVariable, &rThisReaction);
    }

    //This default constructor is needed for serializer
    Dof() noexcept
        : mIsFixed(false),
          mIsActive(false),
          mVariableType(DofTrait<TDataType, Variable<TDataType> >::Id),
          mReactionType(DofTrait<TDataType, Variable<TDataType> >::Id),
          mIndex(),
          mEquationId(IndexType()),
          mpNodalData()
    {
    }

    /// Copy constructor.
    Dof(Dof const& rOther) noexcept = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Dof& operator=(Dof const& rOther) noexcept = default;

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
        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        return GetReference(GetVariable(), p_nodal_data->GetSolutionStepData(), SolutionStepIndex, mVariableType);
    }

    TDataType const& GetSolutionStepValue(IndexType SolutionStepIndex = 0) const
    {
        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        return GetReference(GetVariable(), p_nodal_data->GetSolutionStepData(), SolutionStepIndex, mVariableType);
    }

    template<class TVariableType>
    typename TVariableType::Type& GetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex = 0)
    {
        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        return p_nodal_data->GetSolutionStepData().GetValue(rThisVariable, SolutionStepIndex);
    }

    template<class TVariableType>
    typename TVariableType::Type const& GetSolutionStepValue(const TVariableType& rThisVariable, IndexType SolutionStepIndex = 0) const
    {
        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        return p_nodal_data->GetSolutionStepData().GetValue(rThisVariable, SolutionStepIndex);
    }

    TDataType& GetSolutionStepReactionValue(IndexType SolutionStepIndex = 0)
    {
        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        return GetReference(GetReaction(), p_nodal_data->GetSolutionStepData(), SolutionStepIndex, mReactionType);
    }

    TDataType const& GetSolutionStepReactionValue(IndexType SolutionStepIndex = 0) const
    {
        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        return GetReference(GetReaction(), p_nodal_data->GetSolutionStepData(), SolutionStepIndex, mReactionType);
    }

    ///@}
    ///@name Access
    ///@{

    IndexType Id() const
    {
        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        return p_nodal_data->GetId();
    }

    IndexType GetId() const
    {
        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        return p_nodal_data->GetId();
    }

    /** Returns variable assigned to this degree of freedom. */
    const VariableData& GetVariable() const
    {
        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        return p_nodal_data->GetSolutionStepData().GetVariablesList().GetDofVariable(mIndex);
    }

    /** Returns reaction variable of this degree of freedom. */
    const VariableData& GetReaction() const
    {
        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        auto p_reaction = p_nodal_data->GetSolutionStepData().GetVariablesList().pGetDofReaction(mIndex);
        return (p_reaction == nullptr) ? msNone : *p_reaction;
    }

    template<class TReactionType>
    void SetReaction(TReactionType const& rReaction)
    {
        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        mReactionType = DofTrait<TDataType, TReactionType>::Id;
        p_nodal_data->GetSolutionStepData().pGetVariablesList()->SetDofReaction(&rReaction, mIndex);
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

    /// @brief Mark this @p Dof active.
    /// @details A @p Dof is active if it is present in the linear system that
    ///          a @ref BuilderAndSolver built, and has a valid equation ID
    ///          accessible at @ref Dof::EquationId.
    void Activate() noexcept
    {
        mIsActive = true;
    }

    /// @brief Mark this @p Dof inactive.
    /// @details A @p Dof is inactive if the @ref BuilderAndSolver did not use it during the
    ///          assembly of the linear system, in which case @ref Dof::EquationId is invalid.
    void Deactivate() noexcept
    {
        mIsActive = false;
    }

    void Set(Flags flags)
    {
        KRATOS_TRY
        auto [dof_flags, clean_pointer] = DofFlags::Decompose(mpNodalData);
        dof_flags.Set(flags);
        mpNodalData = dof_flags.Compose(clean_pointer);
        KRATOS_CATCH("")
    }

    void Set(Flags flags, bool value)
    {
        KRATOS_TRY
        auto [dof_flags, clean_pointer] = DofFlags::Decompose(mpNodalData);
        value ? dof_flags.Set(flags) : dof_flags.Unset(flags);
        mpNodalData = dof_flags.Compose(clean_pointer);
        KRATOS_CATCH("")
    }

    SolutionStepsDataContainerType* GetSolutionStepsData()
    {
        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        return &p_nodal_data->GetSolutionStepData();
    }

    void SetNodalData(NodalData* pNewNodalData)
    {
        auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        auto p_variable = &GetVariable();
        auto p_reaction = p_nodal_data->GetSolutionStepData().pGetVariablesList()->pGetDofReaction(mIndex);
        p_nodal_data = pNewNodalData;

        if(p_reaction != nullptr){
            mIndex = p_nodal_data->GetSolutionStepData().pGetVariablesList()->AddDof(p_variable, p_reaction);
        } else{
            mIndex = p_nodal_data->GetSolutionStepData().pGetVariablesList()->AddDof(p_variable);
        }

        mpNodalData = dof_flags.Compose(p_nodal_data);
    }

    bool HasReaction() const
    {
        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        return p_nodal_data->GetSolutionStepData().pGetVariablesList()->pGetDofReaction(mIndex) != nullptr;
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

    /// @brief Indicates whether the @p Dof is present in the system that the @ref BuilderAndSolver assembled.
    /// @details A @p Dof is active if it is present in the linear system that
    ///          a @ref BuilderAndSolver built, and has a valid equation ID
    ///          accessible at @ref Dof::EquationId.
    ///          On the other hand, a @p Dof is inactive if the @ref BuilderAndSolver did not
    ///          use it during the assembly of the linear system, in which case @ref Dof::EquationId
    ///          is invalid.
    bool IsActive() const noexcept
    {
        return mIsActive;
    }

    bool Is(Flags flags) const
    {
        KRATOS_TRY
        [[maybe_unused]] const auto [dof_flags, p_nodal_data] = DofFlags::Decompose(mpNodalData);
        return dof_flags.Is(flags);
        KRATOS_CATCH("")
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
    /// @brief 3-bit-wide bit array that hides its data in the 3 least significant bits of a pointer.
    /// @details List of supported flags:
    ///          - @ref MASTER
    ///          - @ref SLAVE
    ///          - there is space for one extra flag to support.
    class DofFlags
    {
    private:
        static constexpr std::size_t msMask = 0b111ul;

    public:
        /// @brief A pointer with access restricted to @ref DofFlags.
        /// @details This class is meant to catch attempts to dereference pointers
        ///          that were manipulated by @ref DofFlags at compile time.
        template <class T>
        class DirtyPointer
        {
        public:
            DirtyPointer() noexcept = default;

            explicit DirtyPointer(T* Pointer) noexcept : mPointer(Pointer) {}

            void save(Serializer& rSerializer) const
            {
                rSerializer.save("Pointer", mPointer);
            }

            void load(Serializer& rSerializer)
            {
                rSerializer.load("Pointer", mPointer);
            }

        private:
            friend class DofFlags;

            T* mPointer;
        };

        /// @brief Counterpart of the @ref MASTER @ref Flags "flag".
        static constexpr DofFlags Master = 0b001;

        /// @brief Counterpart of the @ref SLAVE @ref Flags "flag".
        static constexpr DofFlags Slave = 0b010;

        //constexpr static DofFlags LAST_SUPPORTED_FLAG_GOES_HERE = 0b100;

        /// @brief Default constructor that unsets all flags.
        constexpr DofFlags() noexcept : mData(0) {}

        /// @brief Convert a @ref Flags instance to @ref DofFlags.
        /// @throws if any unsupported flags are set.
        /// @details @see Dof for more information on supported flags.
        explicit DofFlags(Flags flags)
            : mData(0b000)
        {
            mData |= flags.Is(MASTER);
            mData |= flags.Is(SLAVE) << 1;
            //mData |= flags.Is(LAST_SUPPORTED_FLAG_GOES_HERE) << 2;

            // Check whether any flag is set other than
            // the supported ones, and throw if there are.
            flags.Set(MASTER, false);
            flags.Set(SLAVE, false);
            KRATOS_ERROR_IF(flags) << "unsupported flags set in " << flags;
        }

        /// @brief Break up a dirty pointer into a clean one and the @p DofFlags that corrupted it.
        /// @details @p DirtyPointer is assumed to be a pointer to an instance of @p T (whose alignment
        ///          is at least 8), but with the 3 least significant bits corrupted. The corrupted bits
        ///          should contain the data for a @p DofFlags instance. This function extracts the last
        ///          3 bits and returns the resulting @p DofFlags, as well as the @b valid pointer to @p T
        ///          (last 3 bits are unset). The resulting valid pointer can be safely dereferenced.
        /// @note Use @ref DofFlags::Compose to recover the original dirty pointer.
        template <class T>
        static constexpr std::pair<DofFlags,T*> Decompose(DirtyPointer<T> Dirty) noexcept
        {
            static_assert(1 << 3 <= std::alignment_of_v<T>);
            DofFlags flags(Data(reinterpret_cast<std::size_t>(Dirty.mPointer) & msMask));
            T* clean_pointer = reinterpret_cast<T*>(reinterpret_cast<std::size_t>(Dirty.mPointer) & ~msMask);
            return std::make_pair(flags, clean_pointer);
        }

        /// @brief Inject the 3 bits of data in this @p DofFlags into the 3 least significant bits of the provided pointer.
        /// @note Use @ref DofFlags::Decompose to recover the original valid pointer.
        /// @warning The returned pointer should not be dereferenced directly. Use @ref DofFlags::Decompose to recover the
        ///          original pointer.
        template <class T>
        constexpr DirtyPointer<T> Compose(T* CleanPointer) const
        {
            static_assert(1 << 3 <= std::alignment_of_v<T>);
            KRATOS_ERROR_IF(reinterpret_cast<std::size_t>(CleanPointer) & msMask)
                << "the provided pointer is not valid " << CleanPointer;
            T* dirty_pointer = reinterpret_cast<T*>(reinterpret_cast<std::size_t>(CleanPointer) | mData.to_ulong());
            return DirtyPointer<T>(dirty_pointer);
        }

        constexpr bool Is(DofFlags Other) const noexcept
        {
            return (mData & Other.mData).any();
        }

        bool Is(Flags flags) const
        {
            KRATOS_TRY
            return this->Is(DofFlags(flags));
            KRATOS_CATCH("")
        }

        void Set(Flags flags)
        {
            KRATOS_TRY
            DofFlags supported_flags(flags);
            mData = mData | supported_flags.mData;
            KRATOS_CATCH("")
        }

        void Unset(Flags flags)
        {
            KRATOS_TRY
            DofFlags supported_flags(flags);
            mData = mData & ~supported_flags.mData;
            KRATOS_CATCH("")
        }

    private:
        using Data = std::bitset<3>;

        explicit constexpr DofFlags(Data data) noexcept : mData(data) {}

        Data mData;
    }; // class DofFlags

    ///@name Static Member Variables
    ///@{

    static const Variable<TDataType> msNone;

    ///@}
    ///@name Member Variables
    ///@{

    /** True is is fixed */
    int mIsFixed : 1;

    /// @copydoc Dof::IsActive
    int mIsActive : 1;

    int mVariableType : 4;

    int mReactionType : 4;

    int mIndex : 6;

    /** Equation identificator of the degree of freedom */
    EquationIdType mEquationId : 48;

    /** A pointer to nodal data stored in node which is corresponded to this dof */
    typename DofFlags::template DirtyPointer<NodalData> mpNodalData;

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
        rSerializer.save("IsActive", static_cast<bool>(mIsActive));
        rSerializer.save("EquationId", static_cast<EquationIdType>(mEquationId));
        rSerializer.save("NodalData", mpNodalData);
        rSerializer.save("VariableType", static_cast<int>(mVariableType));
        rSerializer.save("ReactionType", static_cast<int>(mReactionType));
        rSerializer.save("Index", static_cast<int>(mIndex));

    }

    void load(Serializer& rSerializer)
    {
        bool is_fixed;
        rSerializer.load("IsFixed", is_fixed);
        mIsFixed=is_fixed;
        bool is_active;
        rSerializer.load("IsActive", is_active);
        mIsActive = is_active;
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
