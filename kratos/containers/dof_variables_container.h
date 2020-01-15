//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:         BSD License 
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela Dalmau
//                   Tobias Teschemacher
//

#ifndef KRATOS_DOF_VARIABLES_CONTAINER_H
#define	KRATOS_DOF_VARIABLES_CONTAINER_H

// System includes
#include <string>
#include <iostream>

// External includes
#include <boost/iterator/indirect_iterator.hpp>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// A container of Kratos double and component variables.
/** It can be filled using either Kratos::Variable<double> or
 * array_1d<double,3> components (VariableComponent< VectorComponentAdaptor< array_1d<double, 3 > > >).
 * It is used by PeriodicCondition to identify the Dofs where the periodic condition applies.
 * @see PeriodicCondition
 * It is also used by PenaltyCouplingCondition to identify the Dofs which shall be enforced.
 * @see PenaltyCouplingCondition
 */
class DofVariablesContainer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DofVariablesContainer
    KRATOS_CLASS_POINTER_DEFINITION(DofVariablesContainer);

    /// Kratos double Variable
    typedef Variable<double> DoubleVariableType;

    /// Container of pointers to Kratos double variables
    typedef std::vector<const DoubleVariableType*> DoubleVariablesContainerType;

    /// Double Variable iterator
    typedef boost::indirect_iterator<DoubleVariablesContainerType::const_iterator> DoubleVariableConstIterator;

    /// Component of a Kratos array_1d<double,3> Variable
    typedef VariableComponent< VectorComponentAdaptor< array_1d<double, 3 > > > VariableComponentType;

    /// Vector of pointers to Kratos double variables
    typedef std::vector<const VariableComponentType*> ComponentsVariableContainerType;

    /// Component variable iterator
    typedef boost::indirect_iterator<ComponentsVariableContainerType::const_iterator> ComponentsVariableConstIterator;

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DofVariablesContainer():
        mDofDoubleVariables(),
        mDofComponentsVariables()
    {
    }

    /// Copy constructor.
    DofVariablesContainer(DofVariablesContainer const& rOther):
        mDofDoubleVariables(rOther.mDofDoubleVariables),
        mDofComponentsVariables(rOther.mDofComponentsVariables)
    {
    }

    /// Destructor.
    virtual ~DofVariablesContainer()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    DofVariablesContainer& operator=(DofVariablesContainer const& rOther)
    {
        this->mDofDoubleVariables = rOther.mDofDoubleVariables;
        this->mDofComponentsVariables = rOther.mDofComponentsVariables;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Add a scalar variable to the list of variables
    void Add(DoubleVariableType const& rThisVariable)
    {
        KRATOS_ERROR_IF(rThisVariable.Key() == 0)
            << "Adding uninitialized variable to a list of dof variables container: "
            << rThisVariable.Name() << std::endl;

        mDofDoubleVariables.push_back(&rThisVariable);
    }

    /// Add a component of a vector variable to the list of variables
    void Add(VariableComponentType const& rThisVariable)
    {
        KRATOS_ERROR_IF(rThisVariable.Key() == 0)
            << "Adding uninitialized variable to dof variables container: "
            << rThisVariable.Name() << std::endl;

        mDofComponentsVariables.push_back(&rThisVariable);
    }

    /// Erase all information contained in this object.
    void Clear()
    {
        mDofDoubleVariables.clear();
        mDofComponentsVariables.clear();
    }

    ///@}
    ///@name Access
    ///@{

    /// Iterator for the list of scalar variables
    DoubleVariableConstIterator DoubleVariablesBegin() const
    {
        return DoubleVariableConstIterator( mDofDoubleVariables.begin() );
    }

    /// Iterator for the list of scalar variables
    DoubleVariableConstIterator DoubleVariablesEnd() const
    {
        return DoubleVariableConstIterator( mDofDoubleVariables.end() );
    }

    /// Iterator for the list of vector components
    ComponentsVariableConstIterator VariableComponentsBegin() const
    {
        return ComponentsVariableConstIterator( mDofComponentsVariables.begin() );
    }

    /// Iterator for the list of vector components
    ComponentsVariableConstIterator VariableComponentsEnd() const
    {
        return ComponentsVariableConstIterator( mDofComponentsVariables.end() );
    }

    /// Total number of Dof variables (including scalars and vector components)
    SizeType size() const
    {
        return mDofDoubleVariables.size() + mDofComponentsVariables.size();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "DofVariablesContainer";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DofVariablesContainer";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "Double Variables:" << std::endl;
        for (DoubleVariablesContainerType::const_iterator it = mDofDoubleVariables.begin();
            it != mDofDoubleVariables.end(); ++it)
        {
            (*it)->PrintInfo(rOStream);
            rOStream << std::endl;
        }
        rOStream << "Components Variables:" << std::endl;
        for (ComponentsVariableContainerType::const_iterator it = mDofComponentsVariables.begin();
            it != mDofComponentsVariables.end(); ++it)
        {
            (*it)->PrintInfo(rOStream);
            rOStream << std::endl;
        }
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    /// Container for double variables
    DoubleVariablesContainerType mDofDoubleVariables;

    /// Container for variable components
    ComponentsVariableContainerType mDofComponentsVariables;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        std::size_t DoubleVariablesSize= mDofDoubleVariables.size();
        rSerializer.save("DoubleVariablesSize",DoubleVariablesSize);
        for(std::size_t i = 0; i < DoubleVariablesSize; i++)
            rSerializer.save("Variable Name", mDofDoubleVariables[i]->Name());

        std::size_t ComponentsVariablesSize= mDofComponentsVariables.size();
        rSerializer.save("ComponentsVariablesSize", ComponentsVariablesSize);
        for(std::size_t i = 0; i < ComponentsVariablesSize; i++)
            rSerializer.save("Variable Name", mDofComponentsVariables[i]->Name());
    }

    virtual void load(Serializer& rSerializer)
    {
        std::string Name;
        std::size_t DoubleVariablesSize;
        rSerializer.load("DoubleVariablesSize",DoubleVariablesSize);
        for(std::size_t i = 0; i < DoubleVariablesSize; i++)
        {
            rSerializer.load("Variable Name", Name);
            Add( *(static_cast<DoubleVariableType*>(KratosComponents<VariableData>::pGet(Name))) );
        }

        std::size_t ComponentsVariablesSize;
        rSerializer.load("ComponentsVariablesSize", ComponentsVariablesSize);
        for(std::size_t i = 0; i < ComponentsVariablesSize; i++)
        {
            rSerializer.load("Variable Name", Name);
            Add( *(static_cast<VariableComponentType*>(KratosComponents<VariableData>::pGet(Name))) );
        }
    }

    ///@}

}; // Class DofVariablesContainer

///@name Input and output
///@{

/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  DofVariablesContainer& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const DofVariablesContainer& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif  /* KRATOS_DOF_VARIABLES_CONTAINER_H */

