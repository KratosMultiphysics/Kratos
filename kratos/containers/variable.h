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
//

#if !defined(KRATOS_VARIABLE_H_INCLUDED )
#define  KRATOS_VARIABLE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "variable_data.h"
#include "utilities/stl_io.h"


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

/// Variable class contains all information needed to store and retrive data from a data container.
/** Variable class contains all information needed to store and
    retrive data from a data container.  It contains key value which
    is needed for searching in data container. Also a zero value to
    use as a default value in the container. Finally it has the type of
    the Variable as its template parameter so the container can find it
    in its relative part.
*/
template<class TDataType>
class Variable : public VariableData
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Variable
    KRATOS_CLASS_POINTER_DEFINITION(Variable);

    /// type of this variable
    typedef TDataType Type;

    // Type used for key values which defined in VariableData
    typedef VariableData::KeyType KeyType;

    // Type of this varible with given TDataType
    typedef Variable<TDataType> VariableType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor with specific name and zero value */
    Variable(const std::string& NewName, const TDataType Zero = TDataType())
        : VariableData(NewName, sizeof(TDataType)), mZero(Zero)
    {
    }

    /// Copy constructor.
    Variable(const VariableType& rOtherVariable) : VariableData(rOtherVariable), mZero(rOtherVariable.mZero) {}

    /// Destructor.
    ~Variable() override {}

    ///@}
    ///@name Operators
    ///@{
    
    /// Assignment operator.
    VariableType& operator=(const VariableType& rOtherVariable)
    {
        VariableData::operator=(rOtherVariable);
        mZero = rOtherVariable.mZero;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    void* Clone(const void* pSource) const override
    {
        return new TDataType(*static_cast<const TDataType* >(pSource) );
    }

    void* Copy(const void* pSource, void* pDestination) const override
    {
        return new(pDestination) TDataType(*static_cast<const TDataType* >(pSource) );
    }

    void Assign(const void* pSource, void* pDestination) const override
    {
        (*static_cast<TDataType* >(pDestination) ) = (*static_cast<const TDataType* >(pSource) );
    }

    void AssignZero(void* pDestination) const override
    {
        //(*static_cast<TDataType* >(pDestination) ) = mZero;
        new (pDestination) TDataType(mZero);
    }

    void Delete(void* pSource) const override
    {
        delete static_cast<TDataType* >(pSource);
    }

    void Destruct(void* pSource) const override
    {
        static_cast<TDataType* >(pSource)->~TDataType();
    }

    void Print(const void* pSource, std::ostream& rOStream) const override
    {
        rOStream << Name() << " : " << *static_cast<const TDataType* >(pSource) ;
    }

    void Save(Serializer& rSerializer, void* pData) const override
    {
        // I'm saving by the value, it can be done by the pointer to detect shared data. Pooyan.
        rSerializer.save("Data",*static_cast<TDataType* >(pData));
    }

    void Allocate(void** pData) const override
    {
        *pData = new TDataType;
    }

    void Load(Serializer& rSerializer, void* pData) const override
    {
        rSerializer.load("Data",*static_cast<TDataType* >(pData));
    }

    static const VariableType& StaticObject()
    {
        return msStaticObject;
    }

    ///@}
    ///@name Access
    ///@{

    const TDataType& Zero() const
    {
        return mZero;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << Name() << " variable";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Name() << " variable";
    }

    /// Print object's data.
//       virtual void PrintData(std::ostream& rOStream) const;


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

    static const VariableType msStaticObject;

    ///@}
    ///@name Member Variables
    ///@{

    TDataType mZero;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, VariableData );
        rSerializer.save("Zero",mZero);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, VariableData );
        rSerializer.load("Zero",mZero);
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

    /** Default constructor is un accessible due to the fact that
    each variable must have a name defined.*/
    Variable() {}

    ///@}

}; // Class Variable

///@}

template<class TDataType>
const Variable<TDataType> Variable<TDataType>::msStaticObject("NONE");

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TDataType>
inline std::istream& operator >> (std::istream& rIStream,
                                  Variable<TDataType>& rThis);

/// output stream function
template<class TDataType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Variable<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_VARIABLE_H_INCLUDED  defined 


