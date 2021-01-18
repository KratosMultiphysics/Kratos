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
//                   Nelson Lafontaine
//
//

#if !defined(KRATOS_ALL_VARIABLES_DATA_VALUE_CONTAINER_H_INCLUDED )
#define KRATOS_ALL_VARIABLES_DATA_VALUE_CONTAINER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cstddef>
#include <cstring>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "includes/kratos_components.h"


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

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_EXPORT_DLL AllVariablesDataValueContainer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AllVariablesDataValueContainer
    KRATOS_CLASS_POINTER_DEFINITION(AllVariablesDataValueContainer);

    typedef double BlockType;

    /// Type of the container used for variables
    typedef BlockType* ContainerType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AllVariablesDataValueContainer()
        : mOffsets(), mpData(0)
    {
        // Allcating memory
        Allocate();

		// Assigning the zero of each variable to its position
		AssignZero();
    }

    /// Copy constructor.
    AllVariablesDataValueContainer(AllVariablesDataValueContainer const& rOther)
        : mOffsets(), mpData(0)
    {
        // Allcating memory
        Allocate();

		Copy(rOther);

    }


    /// Destructor.
    virtual ~AllVariablesDataValueContainer()
    {
        //Clear();
    }


    ///@}
    ///@name Operators
    ///@{

//       template<class TDataType> const TDataType& operator()(const VariableData& rThisVariable, SizeType QueueIndex = 0) const
// 	{
// 	  return GetValue(rThisVariable, QueueIndex);
// 	}

    template<class TDataType>
    TDataType& operator()(const Variable<TDataType>& rThisVariable)
    {
        return GetValue(rThisVariable);
    }

    template<class TDataType>
    const TDataType& operator()(const Variable<TDataType>& rThisVariable) const
    {
        return GetValue(rThisVariable);
    }

    template<class TAdaptorType>
    typename TAdaptorType::Type& operator()(const VariableComponent<TAdaptorType>& rThisVariable)
    {
        return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
    }

    template<class TAdaptorType>
    const typename TAdaptorType::Type& operator()(const VariableComponent<TAdaptorType>& rThisVariable) const
    {
        return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
    }

//       template<class TDataType> TDataType& operator[](const VariableData& rThisVariable)
// 	{
// 	  return GetValue(rThisVariable, 0);
// 	}

//       template<class TDataType> const TDataType& operator[](const VariableData& rThisVariable) const
// 	{
// 	  return GetValue(rThisVariable, 0);
// 	}

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
        return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable(), 0));
    }

    template<class TAdaptorType> const typename TAdaptorType::Type& operator[](const VariableComponent<TAdaptorType>& rThisVariable) const
    {
        return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable(), 0));
    }

    /// Assignment operator.
    AllVariablesDataValueContainer& operator=(const AllVariablesDataValueContainer& rOther)
    {
		Copy(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    template<class TDataType>
    TDataType& GetValue(const Variable<TDataType>& rThisVariable)
    {
        return *(TDataType*)Position(rThisVariable);
    }


    template<class TDataType>
    const TDataType& GetValue(const Variable<TDataType>& rThisVariable) const
    {
        return *(const TDataType*)Position(rThisVariable);
    }

    template<class TAdaptorType>
    typename TAdaptorType::Type& GetValue(const VariableComponent<TAdaptorType>& rThisVariable)
    {
        return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
    }
    template<class TAdaptorType>
    const typename TAdaptorType::Type& GetValue(const VariableComponent<TAdaptorType>& rThisVariable) const
    {
        return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
    }


    SizeType Size() const
    {
        return mOffsets[KratosComponents<VariableData>::GetComponents().size()];
    }


    template<class TDataType> void SetValue(const Variable<TDataType>& rThisVariable, TDataType const& rValue)
    {
        GetValue(rThisVariable) = rValue;
    }

    template<class TAdaptorType> void SetValue(const VariableComponent<TAdaptorType>& rThisVariable, typename TAdaptorType::Type const& rValue)
    {
        rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable())) = rValue;
    }

    void Clear()
    {
        //DestructAllElements(); //commented for a test. Pooyan.
        if(mpData)
            free(mpData);

        mpData = 0;
    }


    ///@}
    ///@name Access
    ///@{

    BlockType* Data()
    {
        return mpData;
    }

    BlockType* Data(VariableData const & rThisVariable)
    {
        return Position(rThisVariable);
    }

    //SizeType DataSize()
    //{
    //    return mpVariablesList->DataSize() * sizeof(BlockType);
    //}


    ///@}
    ///@name Inquiry
    ///@{

    template<class TDataType> bool Has(const Variable<TDataType>& rThisVariable) const
    {
		// This container should have all the variables. If it is not in the list something has gone wrong with registartion of the variable
        return true;
    }

    template<class TAdaptorType> bool Has(const VariableComponent<TAdaptorType>& rThisVariable) const
    {
		// This container should have all the variables. If it is not in the list something has gone wrong with registartion of the variable
        return true;
    }

    bool IsEmpty()
    {
        return false;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return std::string("All variables data value container");
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "All variables data value container";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {

        for(KratosComponents<VariableData>::ComponentsContainerType::const_iterator i_variable = KratosComponents<VariableData>::GetComponents().begin() ;
			i_variable != KratosComponents<VariableData>::GetComponents().end() ; i_variable++)
			if(i_variable->second->IsNotComponent())
			{
				rOStream <<"    ";
				rOStream << mOffsets[i_variable->second->Key()] << " : ";
				i_variable->second->Print(Position(*(i_variable->second)), rOStream);
				rOStream << std::endl;
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

    std::vector<int> mOffsets;

    ContainerType mpData;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    inline void Allocate()
    {
		int size = KratosComponents<VariableData>::GetComponents().size() + 1; // The + 1 is important to store the total size
		mOffsets.resize(size, double());
		int i = 0;
		const double a = 1.00 / sizeof(BlockType);
        for(KratosComponents<VariableData>::ComponentsContainerType::const_iterator i_variable = KratosComponents<VariableData>::GetComponents().begin() ;
                i_variable != KratosComponents<VariableData>::GetComponents().end() ; i_variable++)
				if (i_variable->second->IsNotComponent())
					mOffsets[++i] = static_cast<int>(i_variable->second->Size() * a);
				else
					mOffsets[++i] = 0;

		for(int j = 1 ; j < size ; j++)
			mOffsets[j] += mOffsets[j-1];

		mpData = (BlockType*)malloc(mOffsets[size-1] * sizeof(BlockType));
    }



    inline void AssignZero()
    {
        for(KratosComponents<VariableData>::ComponentsContainerType::const_iterator i_variable = KratosComponents<VariableData>::GetComponents().begin() ;
                i_variable != KratosComponents<VariableData>::GetComponents().end() ; i_variable++)
				if(i_variable->second->IsNotComponent())
					i_variable->second->AssignZero(Position(i_variable->second));
    }

    void DestructAllElements()
    {
        if(mpData == 0)
            return;
		if(mOffsets.size() != KratosComponents<VariableData>::GetComponents().size() + 1)
			return; //This would cause a memory leak but should only happen in the initialization where I cannot do anything better! Pooyan.

        for(KratosComponents<VariableData>::ComponentsContainerType::const_iterator i_variable = KratosComponents<VariableData>::GetComponents().begin() ;
                i_variable != KratosComponents<VariableData>::GetComponents().end() ; i_variable++)
		{
				if(i_variable->second->IsNotComponent())
					i_variable->second->Destruct(Position(i_variable->second));
		}
    }



	void Copy(AllVariablesDataValueContainer const& rOther)
	{
        for(KratosComponents<VariableData>::ComponentsContainerType::const_iterator i_variable = KratosComponents<VariableData>::GetComponents().begin() ;
                i_variable != KratosComponents<VariableData>::GetComponents().end() ; i_variable++)
		{
			if(i_variable->second->IsNotComponent())
			{
				SizeType offset = mOffsets[i_variable->second->Key()];
				i_variable->second->Copy(rOther.mpData + offset, mpData + offset);
			}
		}
	}

    void AssignData(BlockType* Source, BlockType* Destination)
    {
        for(KratosComponents<VariableData>::ComponentsContainerType::const_iterator i_variable = KratosComponents<VariableData>::GetComponents().begin() ;
                i_variable != KratosComponents<VariableData>::GetComponents().end() ; i_variable++)
				if(i_variable->second->IsNotComponent())
				{
					SizeType offset = mOffsets[i_variable->second->Key()];
					i_variable->second->Assign(Source + offset, Destination + offset);
				}
	}

	inline BlockType* Position(const VariableData * const pThisVariable) const
    {
        return mpData + mOffsets[pThisVariable->Key()];
    }

    inline BlockType* Position(VariableData const & rThisVariable) const
    {
        return mpData + mOffsets[rThisVariable.Key()];
    }

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save("Offsets", mOffsets);
		for(KratosComponents<VariableData>::ComponentsContainerType::const_iterator i_variable = KratosComponents<VariableData>::GetComponents().begin() ;
			i_variable != KratosComponents<VariableData>::GetComponents().end() ; i_variable++)
			if(i_variable->second->IsNotComponent())
			{
				i_variable->second->Save(rSerializer, Position(i_variable->second));
			}
	}

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load("Offsets", mOffsets);

		int size = KratosComponents<VariableData>::GetComponents().size();
		mpData = (BlockType*)malloc(mOffsets[size] * sizeof(BlockType));

		for(KratosComponents<VariableData>::ComponentsContainerType::const_iterator i_variable = KratosComponents<VariableData>::GetComponents().begin() ;
			i_variable != KratosComponents<VariableData>::GetComponents().end() ; i_variable++)
			if(i_variable->second->IsNotComponent())
			{
				i_variable->second->Load(rSerializer, Position(i_variable->second));
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

}; // Class AllVariablesDataValueContainer

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AllVariablesDataValueContainer& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AllVariablesDataValueContainer& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ALL_VARIABLES_DATA_VALUE_CONTAINER_H_INCLUDED  defined
