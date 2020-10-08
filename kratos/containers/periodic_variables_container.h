//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela Dalmau
//
//

#ifndef KRATOS_PERIODIC_VARIABLES_CONTAINER_H
#define	KRATOS_PERIODIC_VARIABLES_CONTAINER_H

// System includes
#include <string>
#include <iostream>


// External includes
#include <boost/iterator/indirect_iterator.hpp>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "containers/variable.h"


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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

/// A container of Kratos variables used to define a periodic boundary condition.
/** It is used by PeriodicCondition to identify the Dofs where the periodic condition applies.
 * @see PeriodicCondition
 */
class PeriodicVariablesContainer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PeriodicVariablesContainer
    KRATOS_CLASS_POINTER_DEFINITION(PeriodicVariablesContainer);

    /// Kratos double Variable
    typedef Variable<double> DoubleVariableType;

    /// Container of pointers to Kratos double variables
    typedef std::vector<const DoubleVariableType*> DoubleVariablesContainerType;

    /// Double Variable iterator
    typedef boost::indirect_iterator<DoubleVariablesContainerType::const_iterator> DoubleVariablesConstIterator;

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PeriodicVariablesContainer():
        mPeriodicDoubleVars()
    {
    }

    /// Copy constructor.
    PeriodicVariablesContainer(PeriodicVariablesContainer const& rOther):
        mPeriodicDoubleVars(rOther.mPeriodicDoubleVars)
    {
    }

    /// Destructor.
    virtual ~PeriodicVariablesContainer()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    PeriodicVariablesContainer& operator=(PeriodicVariablesContainer const& rOther)
    {
        this->mPeriodicDoubleVars = rOther.mPeriodicDoubleVars;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Add a scalar variable to the list of variables where periodic conditions will be imposed.
    void Add(DoubleVariableType const& rThisVariable)
    {
        if(rThisVariable.Key()== 0)
            KRATOS_THROW_ERROR(std::logic_error,
                         "Adding uninitialized variable to a list of periodic variables: ",rThisVariable.Name());

        mPeriodicDoubleVars.push_back(&rThisVariable);
    }

    /// Erase all information contained in this object.
    void Clear()
    {
        mPeriodicDoubleVars.clear();
    }

    ///@}
    ///@name Access
    ///@{

    /// Iterator for the list of scalar variables
    DoubleVariablesConstIterator DoubleVariablesBegin() const
    {
        return DoubleVariablesConstIterator( mPeriodicDoubleVars.begin() );
    }

    /// Iterator for the list of scalar variables
    DoubleVariablesConstIterator DoubleVariablesEnd() const
    {
        return DoubleVariablesConstIterator( mPeriodicDoubleVars.end() );
    }

    /// Total number of periodic variables (including scalars and vector components)
    SizeType size() const
    {
        return mPeriodicDoubleVars.size();
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
        buffer << "PeriodicVariablesContainer";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "PeriodicVariablesContainer";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "Double Variables:" << std::endl;
        for (DoubleVariablesContainerType::const_iterator it = mPeriodicDoubleVars.begin(); it != mPeriodicDoubleVars.end(); ++it)
        {
            (*it)->PrintInfo(rOStream);
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

    /// Container for double variables
    DoubleVariablesContainerType mPeriodicDoubleVars;


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

    virtual void save(Serializer& rSerializer) const
    {
        std::size_t DoubleVarSize= mPeriodicDoubleVars.size();
        rSerializer.save("DoubleVarSize",DoubleVarSize);
        for(std::size_t i = 0; i < DoubleVarSize; i++)
            rSerializer.save("Variable Name", mPeriodicDoubleVars[i]->Name());
    }

    virtual void load(Serializer& rSerializer)
    {
        std::string Name;
        std::size_t DoubleVarSize;
        rSerializer.load("DoubleVarSize",DoubleVarSize);
        for(std::size_t i = 0; i < DoubleVarSize; i++)
        {
            rSerializer.load("Variable Name", Name);
            Add( *(static_cast<DoubleVariableType*>(KratosComponents<VariableData>::pGet(Name))) );
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

}; // Class PeriodicVariablesContainer

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  PeriodicVariablesContainer& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const PeriodicVariablesContainer& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif	/* KRATOS_PERIODIC_VARIABLES_CONTAINER_H */

