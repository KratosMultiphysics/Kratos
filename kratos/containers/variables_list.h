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


// External includes
#include <boost/iterator/indirect_iterator.hpp>


// Project includes
#include "includes/define.h"
#include "includes/kratos_components.h"
#include "containers/variable.h"


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
		KRATOS_CLASS_POINTER_DEFINITION(VariablesList);



		typedef std::size_t SizeType;

		typedef std::size_t IndexType;

		typedef double BlockType;

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
		//typedef typename VariablesContainerType::iterator iterator;
		//typedef typename VariablesContainerType::const_iterator const_iterator;
		//typedef typename VariablesContainerType::reverse_iterator reverse_iterator;
		//typedef typename VariablesContainerType::const_reverse_iterator const_reverse_iterator;


		///@}
		///@name Life Cycle
		///@{

		/// Default constructor. mPosition should have at least on entry.
		VariablesList() : mDataSize(0), mPositions(1, static_cast<IndexType>(-1)), mVariables()
		{
		}

		template <class TInputIteratorType>
		VariablesList(TInputIteratorType First, TInputIteratorType Last)
		{
			for (; First != Last; First++)
				push_back(*First);
		}

		/// Copy constructor.
		VariablesList(VariablesList const& rOther) : mDataSize(rOther.mDataSize)
			, mPositions(rOther.mPositions)
			, mVariables(rOther.mVariables) {}

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
			mPositions = rOther.mPositions;
			mVariables = rOther.mVariables;

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

			mVariables.swap(rOther.mVariables);
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
			mVariables.clear();
			mPositions.clear();
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

			for (auto i_variable = mVariables.begin(); i_variable != mVariables.end(); i_variable++)
				if ((*i_variable)->Key() == rThisVariable.Key())
					return true;

			return false;
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
		}


		///@}
	private:
		///@name Static Member Variables
		///@{


		///@}
		///@name Member Variables
		///@{

		SizeType mDataSize;

		double mSizeInverse;

		PositionsContainerType mPositions;

		VariablesContainerType mVariables;

		///@}
		///@name Private Operators
		///@{


		///@}
		///@name Private Operations
		///@{


		void SetPosition(IndexType TheIndex, SizeType ThePosition) {
			//if (mPositions.size() <= TheIndex)
			//	mPositions.resize(TheIndex + 1, static_cast<IndexType>(-1));
			if (mPositions.empty())
				ResizePositions();

			if (mPositions[TheIndex % mPositions.size()] < mDataSize) // The position is ocupied and a resize  (as re hash) is needed
				ResizePositions();

			mPositions[TheIndex % mPositions.size()] = ThePosition;

			mSizeInverse = 1.00 / mPositions.size();
		}

		SizeType GetPosition(IndexType TheIndex) const {
			// This is the fastest way I found to do the TheIndex % mPositions.size(); Pooyan.
			SizeType index = TheIndex - mPositions.size() * static_cast<SizeType>(TheIndex * mSizeInverse);
			//SizeType index = GetModulus(TheIndex);
			//SizeType index = ( ( TheIndex * 193) + ( ( TheIndex & 0xffff0000) >> 16) * 5373) & (mPositions.size() - 1);
			return mPositions[index];
		}

		void ResizePositions() {
			bool size_is_ok = false;
			SizeType new_size = mPositions.size();
			while (size_is_ok != true) {
				new_size++;
				PositionsContainerType new_positions(new_size, static_cast<IndexType>(-1));
				size_is_ok = true;

				for (auto i_variable = mVariables.begin(); i_variable != mVariables.end(); i_variable++)
					if (new_positions[(*i_variable)->Key() % new_size] > mDataSize)
						new_positions[(*i_variable)->Key() % new_size] = mPositions[(*i_variable)->Key() % mPositions.size()];
					else {
						size_is_ok = false;
						break;
					}
					if (size_is_ok)
						mPositions.swap(new_positions);
					KRATOS_WATCH(new_size);
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
				rSerializer.save("Variable Name", mVariables[i]->Name());
			}
		}

		virtual void load(Serializer& rSerializer)
		{
			std::size_t size;
			rSerializer.load("Size", size);
			std::string name;
			for (std::size_t i = 0; i < size; i++)
			{
				rSerializer.load("Variable Name", name);
				Add(*KratosComponents<VariableData>::pGet(name));
			}
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

	}; // Class VariablesList

	   ///@}

	   ///@name Type Definitions
	   ///@{


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


