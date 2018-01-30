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
	           


#if !defined(KRATOS_VARIABLE_DATA_H_INCLUDED )
#define  KRATOS_VARIABLE_DATA_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <cstddef>


// External includes


// Project includes
#include "includes/define.h"
#include "utilities/counter.h"
#include "includes/serializer.h"


namespace Kratos
{
///@addtogroup Kratos
///@{

///@name Kratos Classes
///@{

/// This class is the base of variables and variable's components which contains their common data.
/** This class hold variables name and key and also adaptor type for variables components.
    It also has static method for generating a key based on the name of the variable
*/
class KRATOS_API(KRATOS_CORE) VariableData
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VariableData
    KRATOS_CLASS_POINTER_DEFINITION(VariableData);

    typedef std::size_t KeyType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Copy constructor
    VariableData(const VariableData& rOtherVariable);

    /// Destructor.
    virtual ~VariableData() {}


    ///@}
    ///@name Operators
    ///@{

    /** This operator return the key. by this method user can use Variable
    as argument for the places which variable key is needed. */
    operator size_t() const
    {
        return mKey;
    }

    ///@}
    ///@name Operations
    ///@{

    virtual void* Clone(const void* pSource) const;

    virtual void* Copy(const void* pSource, void* pDestination) const;

    virtual void Assign(const void* pSource, void* pDestination) const;

    virtual void AssignZero(void* pDestination) const;

    virtual void Destruct(void* pSource) const;

    virtual void Delete(void* pSource) const;

    virtual void Print(const void* pSource, std::ostream& rOStream) const;

    virtual void Allocate(void** pData) const;

    virtual void Save(Serializer& rSerializer, void* pData) const;

    virtual void Load(Serializer& rSerializer, void* pData) const;


    ///@}
    ///@name Access
    ///@{

    KeyType Key() const
    {
        return mKey;
    }

    /// NOTE: This function is for internal use and not
    /// to change arbitrary any variable's key
   void SetKey(KeyType NewKey);

    const std::string& Name() const
    {
        return mName;
    }

    std::size_t Size() const
    {
        return mSize;
    }

	bool IsComponent()
	{
		return mIsComponent;
	}

	bool IsNotComponent()
	{
		return !mIsComponent;
	}


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;

    ///@}
    ///@name Statics
    ///@{

	/// This static method generates a uinque key for given name and flags.
	/** The generated key contains a 32-bit uique hash and following information:
		- Copyable : if the is a value type and can be copied by memcopy
		- Component: for component of another variables
		- Component index: The index if is component
		- Size: size of the variable in number of double. 
		The order is as follow:

		64           size         32-bit hash                comp. index 0
		 |-----------|----|---------------------------------------------|-|
		*/
	static KeyType GenerateKey(const std::string& Name, std::size_t Size, std::size_t ComponentIndex);

    ///@}
    ///@name Friends
    ///@{

    friend bool operator==(const VariableData& rFirstVariable, const VariableData& rSecondVariable)
    {
        return (rFirstVariable.mKey == rSecondVariable.mKey);
    }

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{

    VariableData& operator=(const VariableData& rOtherVariable)
    {
        mName = rOtherVariable.mName;
        mKey = rOtherVariable.mKey;
        mSize = rOtherVariable.mSize;
		mIsComponent = rOtherVariable.mIsComponent;

        return *this;
    }

    ///@}
    ///@name Protected LifeCycle
    ///@{

    /// Constructor.
	VariableData(const std::string& NewName, std::size_t NewSize, bool Iscomponent = false);


    /** default constructor is to be used only with serialization due to the fact that
    each variable must have a name defined.*/
    VariableData() {}

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::string mName;

    /** Key value of this variable. Each variable will be locate by this
    value in each data structure. Variable constructor will initialize it. */
    KeyType mKey;

    std::size_t mSize;

	bool mIsComponent;

    ///@}
    ///@name Private Operations
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

    ///@}

}; // Class VariableData


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  VariableData& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const VariableData& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_VARIABLE_DATA_H_INCLUDED  defined 


