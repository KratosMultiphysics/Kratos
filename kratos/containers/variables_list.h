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

#if !defined(KRATOS_VARIABLES_LIST_H_INCLUDED )
#define  KRATOS_VARIABLES_LIST_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <vector>
#include <atomic>

// External includes
#include <boost/iterator/indirect_iterator.hpp>


// Project includes
#include "includes/define.h"
#include "includes/kratos_components.h"
#include "containers/variable.h"

#ifdef KRATOS_DEBUG
#include "utilities/openmp_utils.h"
#endif


namespace Kratos
{

	///@name Kratos Classes
	///@{

	/// Holds a list of variables and their position in VariablesListDataValueContainer
	/** This class works tightly with VariablesListDataValueContainer and provides the
	the positions of variables for that containers
	*/
	class VariablesList
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of VariablesList
		KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(VariablesList);



		typedef std::size_t SizeType;

		typedef std::size_t IndexType;

		typedef double BlockType;

		typedef std::vector<IndexType> KeysContainerType;

		typedef std::vector<IndexType> PositionsContainerType;

		typedef std::vector<const VariableData*> VariablesContainerType;

		typedef VariableData data_type;
		typedef const VariableData* value_type;
		typedef const VariableData* const_pointer;
		typedef VariableData const& const_reference;

		typedef boost::indirect_iterator<VariablesContainerType::const_iterator>          const_iterator;
		typedef boost::indirect_iterator<VariablesContainerType::const_reverse_iterator>  const_reverse_iterator;

		typedef VariablesContainerType::size_type size_type;
		typedef VariablesContainerType::const_iterator ptr_const_iterator;
		typedef VariablesContainerType::const_reverse_iterator ptr_const_reverse_iterator;
		typedef VariablesContainerType::difference_type difference_type;


		///@}
		///@name Life Cycle
		///@{

		/// Default constructor. mPosition should have at least on entry.
        VariablesList()
		: mReferenceCounter(0)
		{}

        template <class TInputIteratorType>
		VariablesList(TInputIteratorType First, TInputIteratorType Last)
		: mReferenceCounter(0)
		{
			for (; First != Last; First++)
				Add(*First);
		}

		/// Copy constructor.
		VariablesList(VariablesList const& rOther) : mDataSize(rOther.mDataSize)
			, mHashFunctionIndex(rOther.mHashFunctionIndex)
			, mKeys(rOther.mKeys)
			, mPositions(rOther.mPositions)
			, mVariables(rOther.mVariables)
			, mDofVariables(rOther.mDofVariables)
			, mDofReactions(rOther.mDofReactions)
			, mReferenceCounter(0) 
			{}

		/// Destructor.
		virtual ~VariablesList()
		{
		}


		///@}
		///@name Operators
		///@{

		/// Assignment operator.
		VariablesList& operator=(VariablesList const& rOther)
		{
			mDataSize = rOther.mDataSize;
			mHashFunctionIndex = rOther.mHashFunctionIndex;
			mKeys = rOther.mKeys;
			mPositions = rOther.mPositions;
			mVariables = rOther.mVariables;
			mDofVariables=rOther.mDofVariables;
			mDofReactions=rOther.mDofReactions;

			return *this;
		}

		IndexType operator()(IndexType VariableKey) const
		{
			return GetPosition(VariableKey);
		}

		template<class TDataType>
		IndexType operator()(Variable<TDataType> const& ThisVariable) const
		{
			return GetPosition(ThisVariable.Key());
		}

		const VariableData* operator[](IndexType Index) const
		{
			return mVariables[Index];
		}


		bool operator==(const VariablesList& r) const // nothrow
		{
			if (size() != r.size())
				return false;
			else
				return std::equal(mPositions.begin(), mPositions.end(), r.mPositions.begin()) &&
				std::equal(mVariables.begin(), mVariables.end(), r.mVariables.begin());
		}

		//public API of intrusive_ptr
		unsigned int use_count() const noexcept
		{
			return mReferenceCounter;
		}
		///@}
		///@name Operations
		///@{

		const_iterator             begin() const
		{
			return const_iterator(mVariables.begin());
		}
		const_iterator             end() const
		{
			return const_iterator(mVariables.end());
		}
		const_reverse_iterator     rbegin() const
		{
			return const_reverse_iterator(mVariables.rbegin());
		}
		const_reverse_iterator     rend() const
		{
			return const_reverse_iterator(mVariables.rend());
		}
		ptr_const_iterator         ptr_begin() const
		{
			return mVariables.begin();
		}
		ptr_const_iterator         ptr_end() const
		{
			return mVariables.end();
		}
		ptr_const_reverse_iterator ptr_rbegin() const
		{
			return mVariables.rbegin();
		}
		ptr_const_reverse_iterator ptr_rend() const
		{
			return mVariables.rend();
		}

		const_reference  front() const /* nothrow */
		{
			assert(!IsEmpty());
			return *mVariables.front();
		}
		const_reference  back() const  /* nothrow */
		{
			assert(!IsEmpty());
			return *mVariables.back();
		}

		size_type size() const
		{
			return mVariables.size();
		}

		size_type max_size() const
		{
			return mVariables.max_size();
		}

		void swap(VariablesList& rOther)
		{
			SizeType temp = mDataSize;
			mDataSize = rOther.mDataSize;
			rOther.mDataSize = temp;

			temp = mHashFunctionIndex;
			mHashFunctionIndex = rOther.mHashFunctionIndex;
			rOther.mHashFunctionIndex = temp;

			mVariables.swap(rOther.mVariables);

			mDofVariables.swap(rOther.mDofVariables);
			mDofReactions.swap(rOther.mDofReactions);
			
			mKeys.swap(rOther.mKeys);
			mPositions.swap(rOther.mPositions);
		}

		template<class TOtherDataType>
		void push_back(TOtherDataType const& x)
		{
			Add(x);
		}

		void clear()
		{
			mDataSize = 0;
			mHashFunctionIndex = 0;
			mVariables.clear();
			mDofVariables.clear();
			mDofReactions.clear();
			mKeys = {static_cast<IndexType>(-1)};
			mPositions = {static_cast<IndexType>(-1)};
		}


		void Add(VariableData const& ThisVariable)
		{
			if (ThisVariable.Key() == 0)
				KRATOS_THROW_ERROR(std::logic_error,
					"Adding uninitialize variable to this variable list. Check if all variables are registered before kernel initialization", "");

			if (Has(ThisVariable))
				return;

			mVariables.push_back(&ThisVariable);
			SetPosition(ThisVariable.Key(), mDataSize);
			const SizeType block_size = sizeof(BlockType);
			mDataSize += static_cast<SizeType>(((block_size - 1) + ThisVariable.Size()) / block_size);
		}

		int AddDof(VariableData const* pThisDofVariable){
			
			for(std::size_t dof_index = 0 ; dof_index < mDofVariables.size() ; dof_index++){
				if(*mDofVariables[dof_index] == *pThisDofVariable){
					return static_cast<int>(dof_index);
				}
			}

#ifdef KRATOS_DEBUG
        if(OpenMPUtils::IsInParallel() != 0)
            KRATOS_ERROR << "attempting to call AddDof for: " << pThisDofVariable << ". Unfortunately the Dof was not added before and the operations is not threadsafe (this function is being called within a parallel region)" << std::endl;
#endif 
			mDofVariables.push_back(pThisDofVariable);
			mDofReactions.push_back(nullptr);

			KRATOS_DEBUG_ERROR_IF(mDofVariables.size()>64) << "Adding too many dofs to the node. Each node only can store 64 Dofs." << std::endl;

			return mDofVariables.size() - 1;
		}

		int AddDof(VariableData const* pThisDofVariable, VariableData const* pThisDofReaction){
			
			for(std::size_t dof_index = 0 ; dof_index < mDofVariables.size() ; dof_index++){
				if(*mDofVariables[dof_index] == *pThisDofVariable){
					mDofReactions[dof_index] = pThisDofReaction;
					return static_cast<int>(dof_index);
				}
			}

#ifdef KRATOS_DEBUG
        if(OpenMPUtils::IsInParallel() != 0)
            KRATOS_ERROR << "attempting to call AddDof for: " << pThisDofVariable << ". Unfortunately the Dof was not added before and the operations is not threadsafe (this function is being called within a parallel region)" << std::endl;
#endif 
			mDofVariables.push_back(pThisDofVariable);
			mDofReactions.push_back(pThisDofReaction);

			KRATOS_DEBUG_ERROR_IF(mDofVariables.size()>64) << "Adding too many dofs to the node. Each node only can store 64 Dofs." << std::endl;

			return mDofVariables.size() - 1;
		}

		const VariableData& GetDofVariable(int DofIndex) const {
			return *mDofVariables[DofIndex];
		}

		const VariableData* pGetDofReaction(int DofIndex) const {
			return mDofReactions[DofIndex];
		}

		void SetDofReaction(VariableData const* pThisDofReaction, int DofIndex) {
			KRATOS_DEBUG_ERROR_IF(static_cast<std::size_t>(DofIndex) >= mDofReactions.size()) << "The given dof with index = " << DofIndex  << " is not stored in this variables list" << std::endl;
			mDofReactions[DofIndex] = pThisDofReaction;
		}

		IndexType Index(IndexType VariableKey) const
		{
			return GetPosition(VariableKey);
		}

		template<class TDataType>
		IndexType Index(Variable<TDataType> const& ThisVariable) const
		{
			return GetPosition(ThisVariable.Key());
		}

		IndexType Index(const VariableData* pThisVariable) const
		{
			return GetPosition(pThisVariable->Key());
		}


		///@}
		///@name Access
		///@{

		SizeType DataSize() const
		{
			return mDataSize;
		}


		VariablesContainerType const& Variables()
		{
			return mVariables;
		}

		///@}
		///@name Inquiry
		///@{

		bool Has(const VariableData& rThisVariable) const
		{
			if (mPositions.empty())
				return false;

			if (rThisVariable.Key() == 0)
				return false;

			return mKeys[GetHashIndex(rThisVariable.Key(), mKeys.size(), mHashFunctionIndex)] == rThisVariable.Key();
		}

		bool IsEmpty() const
		{
			return mVariables.empty();
		}

		///@}
		///@name Input and output
		///@{

		/// Turn back information as a string.
		virtual std::string Info() const
		{
			return "variables list";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
		}

		/// Print object's data.
		virtual void PrintData(std::ostream& rOStream) const
		{
			rOStream << " with " << size() << " variables";
			rOStream << " (size : " << mDataSize << " blocks of " << sizeof(BlockType) << " bytes) " << std::endl;
			for (IndexType i = 0; i < mVariables.size(); ++i)
				rOStream << "    " << mVariables[i]->Name() << " \t-> " << GetPosition(mVariables[i]->Key()) << std::endl;

			rOStream << " with " << mDofVariables.size() << " Dofs:";	
			for (IndexType i = 0; i < mDofVariables.size(); ++i)
				rOStream << "    [" << mDofVariables[i]->Name() << " ,  " << ((mDofReactions[i] == nullptr) ? "NONE" : mDofReactions[i]->Name()) << "]" << std::endl;
		}


		///@}
	private:
		///@name Member Variables
		///@{
		SizeType mDataSize = 0;

		SizeType mHashFunctionIndex = 0;

		KeysContainerType mKeys = {static_cast<IndexType>(-1)};

		PositionsContainerType mPositions = {static_cast<IndexType>(-1)};

		VariablesContainerType mVariables;

		VariablesContainerType mDofVariables;

		VariablesContainerType mDofReactions;

		///@}
		///@name Private Operators
		///@{
		//this block is needed for refcounting
		mutable std::atomic<int> mReferenceCounter;

		friend void intrusive_ptr_add_ref(const VariablesList* x)
		{
			x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
		}

		friend void intrusive_ptr_release(const VariablesList* x)
		{
			if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
			std::atomic_thread_fence(std::memory_order_acquire);
			delete x;
			}
		}

		///@}
		///@name Private Operations
		///@{


		void SetPosition(IndexType Key, SizeType ThePosition) {
			if (mPositions.empty())
				ResizePositions();

			if (mPositions[GetHashIndex(Key,mPositions.size(),mHashFunctionIndex)] < mDataSize) // The position is ocupied and a resize  (as re hash) is needed
				ResizePositions();

			mKeys[GetHashIndex(Key, mPositions.size(), mHashFunctionIndex)] = Key;
			mPositions[GetHashIndex(Key, mPositions.size(), mHashFunctionIndex)] = ThePosition;
		}

		SizeType GetHashIndex(std::size_t Key, std::size_t TableSize, std::size_t HashFunctionIndex) const{
			return (Key >> HashFunctionIndex) & (TableSize - 1);
		}

		SizeType GetPosition(IndexType Key) const {
			SizeType index = GetHashIndex(Key,mPositions.size(),mHashFunctionIndex);
			return mPositions[index];
		}

		void ResizePositions() {
			bool size_is_ok = false;
			std::size_t new_size = mPositions.size();
			SizeType new_hash_function_index = 0;
			while (size_is_ok != true) {
				new_hash_function_index++;
				if (new_hash_function_index > 31) {
					new_hash_function_index = 0;
					new_size *= 2;
				}
				KeysContainerType new_keys(new_size, static_cast<IndexType>(-1));
				PositionsContainerType new_positions(new_size, static_cast<IndexType>(-1));
				size_is_ok = true;

					for (auto i_variable = mVariables.begin(); i_variable != mVariables.end(); i_variable++)
						if (new_positions[GetHashIndex((*i_variable)->Key(), new_size, new_hash_function_index)] > mDataSize) {
							new_positions[GetHashIndex((*i_variable)->Key(), new_size, new_hash_function_index)] = mPositions[GetHashIndex((*i_variable)->Key(), mPositions.size(), mHashFunctionIndex)];
							new_keys[GetHashIndex((*i_variable)->Key(), new_size, new_hash_function_index)] = (*i_variable)->Key();
						}
						else {
							size_is_ok = false;
							break;
						}

				if (size_is_ok) {
					mPositions.swap(new_positions);
					mKeys.swap(new_keys);
					mHashFunctionIndex = new_hash_function_index;
				}
			}
		}

		///@}
		///@name Serialization
		///@{

		friend class Serializer;


		virtual void save(Serializer& rSerializer) const
		{
			std::size_t size = mVariables.size();
			rSerializer.save("Size", size);
			for (std::size_t i = 0; i < size; i++)
			{
				rSerializer.save("VariableName", mVariables[i]->Name());
			}

			std::size_t dof_size = mDofVariables.size();
			rSerializer.save("DofSize", dof_size);
			for (std::size_t i = 0; i < dof_size; i++)
			{
				rSerializer.save("DofVariableName", mDofVariables[i]->Name());
				if(mDofReactions[i] == nullptr){
					rSerializer.save("HasReaction", false);
				}
				else{
					rSerializer.save("HasReaction", true);
					rSerializer.save("DofReactionName", mDofReactions[i]->Name());
				}
			}
		}

		virtual void load(Serializer& rSerializer)
		{
			std::size_t size;
			rSerializer.load("Size", size);
			std::string name;
			for (std::size_t i = 0; i < size; i++)
			{
				rSerializer.load("VariableName", name);
				Add(*KratosComponents<VariableData>::pGet(name));
			}
			rSerializer.load("DofSize", size);
			for (std::size_t i = 0; i < size; i++)
			{
				rSerializer.load("DofVariableName", name);
				bool has_reaction;
				rSerializer.load("HasReaction", has_reaction);
				
				if(has_reaction){
					std::string reaction_name;
					rSerializer.load("DofReactionName", reaction_name);
					AddDof(KratosComponents<VariableData>::pGet(name), KratosComponents<VariableData>::pGet(reaction_name));

				}
				else{
            		AddDof(KratosComponents<VariableData>::pGet(name), nullptr);
				}

			}

		}

	}; // Class VariablesList
	   ///@}
	   ///@name Input and output
	   ///@{


	   /// input stream function
	inline std::istream& operator >> (std::istream& rIStream,
		VariablesList& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream,
		const VariablesList& rThis)
	{
		rThis.PrintInfo(rOStream);
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@}


}  // namespace Kratos.

#endif // KRATOS_VARIABLES_LIST_H_INCLUDED  defined
