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


#if !defined(KRATOS_VECTOR_COMPONENT_ADAPTOR_H_INCLUDED )
#define  KRATOS_VECTOR_COMPONENT_ADAPTOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "variable.h"


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
template<class TVectorType>
class VectorComponentAdaptor
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of VectorComponentAdaptor
    KRATOS_CLASS_POINTER_DEFINITION(VectorComponentAdaptor);

    typedef typename TVectorType::value_type Type;

    typedef TVectorType SourceType;

    typedef Variable<TVectorType>  SourceVariableType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    VectorComponentAdaptor(const SourceVariableType& rSourceVariable, int ComponentIndex)
        : mpSourceVariable(&rSourceVariable), mComponentIndex(ComponentIndex)
    {}

    /// Copy constructor.
    VectorComponentAdaptor(const VectorComponentAdaptor& rOther)
        : mpSourceVariable(rOther.mpSourceVariable), mComponentIndex(rOther.mComponentIndex)
    {}

    /// Destructor.
    virtual ~VectorComponentAdaptor() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Type& GetValue(SourceType& rValue) const
    {
        return rValue[mComponentIndex];
    }

    const Type& GetValue(const SourceType& rValue) const
    {
        return rValue[mComponentIndex];
    }


    static VectorComponentAdaptor const& StaticObject()
    {
        return msStaticObject;
    }

    ///@}
    ///@name Access
    ///@{

    const SourceVariableType& GetSourceVariable() const
    {
        return *mpSourceVariable;
    }

    int GetComponentIndex() const {
        return mComponentIndex;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << mpSourceVariable->Name() << " vector component " << mComponentIndex << " adaptor";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << mpSourceVariable->Name() << " vector component " << mComponentIndex << " adaptor";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
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

    static VectorComponentAdaptor const msStaticObject;

    ///@}
    ///@name Member Variables
    ///@{

    const SourceVariableType* mpSourceVariable;

    int mComponentIndex;

    ///@}
    ///@name Serialization
    ///@{

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Default constructor.
    VectorComponentAdaptor();

    ///@}

}; // Class VectorComponentAdaptor

///@}

template<class TVectorType>
const VectorComponentAdaptor<TVectorType> VectorComponentAdaptor<TVectorType>::msStaticObject = VectorComponentAdaptor<TVectorType>(VectorComponentAdaptor<TVectorType>::SourceVariableType::StaticObject(), 0);

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TDataType>
inline std::istream& operator >> (std::istream& IStream,
                                  VectorComponentAdaptor<TDataType>& rThis);

/// output stream function
template<class TDataType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const VectorComponentAdaptor<TDataType>& rThis)
{
    rThis.PrintInfo(OStream);
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_VECTOR_COMPONENT_ADAPTOR_H_INCLUDED  defined


