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
    
    /**
     * Clone creates a copy of the object using a copy constructor of the class. 
     * It is useful to avoid shallow copying of complex objects and also without 
     * actually having information about the variable type.
     * @param pSource The pointer of the variable to be cloned
     * @return A raw pointer of the variable
     */
    virtual void* Clone(const void* pSource) const;

    /**
     * Copy is very similar to Clone except that it also the destination 
     * pointer also passed to it. It is a helpful method specially 
     * to create a copy of heterogeneous data arrays
     * @param pSource The pointer of the variable to be copied
     * @param pDestination The pointer of the destination variable
     * @return A raw pointer of the variable
     */
    virtual void* Copy(const void* pSource, void* pDestination) const;

    /**
     * Assign is very similar to Copy. It just differs in using an assignment 
     * operator besides the copy constructor. Copy creates a new object while 
     * Assign does the assignment for two existing objects. 
     * @param pSource The pointer of the value to be assigned
     * @param pDestination The pointer of the destination value
     */
    virtual void Assign(const void* pSource, void* pDestination) const;

    /**
     * AssignZero is a special case of Assign for which variable zero value used as source. 
     * This method is useful for initializing arrays or resetting values in memory.
     * @param pDestination The pointer of the destination variable
     */
    virtual void AssignZero(void* pDestination) const;

    /**
     *  Delete removes an object of variable type from memory. It calls a 
     * destructor of objects to prevent memory leak and frees the memory 
     * allocated for this object assuming that the object is allocated in heap.
     * @param pSource The pointer of the variable to be deleted
     */
    virtual void Delete(void* pSource) const;
    
    /**
     *  Destruct eliminates an object maintaining the memory it is using. 
     * However, the unlike Delete it does nothing with the memory allocated to it. 
     * So it is very useful in case of reallocating a part of the memory.
     * @param pSource The pointer of the variable to be destructed
     */
    virtual void Destruct(void* pSource) const;

    /**
     *  Print is an auxiliary method to produce output of given variable 
     * knowing its address. For example writing an heterogenous container 
     * in an output stream can be done using this method. Point assumes 
     * that the streaming operator is defined for the variable type.
     * @param pSource The pointer of the variable to be printed
     * @param rOStream The stream used to print the information
     */
    virtual void Print(const void* pSource, std::ostream& rOStream) const;

    /**
     * PrintData is an auxiliary method to produce output only the value of given variable 
     * knowing its address. For example writing an heterogenous container 
     * in an output stream can be done using this method. Point assumes 
     * that the streaming operator is defined for the variable type.
     * @param pSource The pointer of the variable to be printed
     * @param rOStream The stream used to print the information
     */
    virtual void PrintData(const void* pSource, std::ostream& rOStream) const;

    /**
     * This method allocates the data of the variable
     * @param pData A pointer to the data to be allocated
     */
    virtual void Allocate(void** pData) const;
    
    /**
     * The save operation which backups the data of the class
     * @param rSerializer The serializer used to preserve the information
     * @param pData A pointer to the data to be saved
     */
    virtual void Save(Serializer& rSerializer, void* pData) const;

    /**
     * The load operation which restores the data of the class
     * @param rSerializer The serializer used to preserve the information
     * @param pData A pointer to the data to be loaded
     */
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

    bool IsComponent() const
    {
        return mIsComponent;
    }

    bool IsNotComponent() const
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
    static KeyType GenerateKey(const std::string& Name, std::size_t Size, bool IsComponent, char ComponentIndex);

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
    VariableData(const std::string& NewName, std::size_t NewSize, bool Iscomponent = false, char ComponentIndex = 0);


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


