//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Jordi Cotela
//


#if !defined(KRATOS_GENERIC_VARIABLES_LIST_H_INCLUDED )
#define  KRATOS_GENERIC_VARIABLES_LIST_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include "boost/variant.hpp"
#include "boost/visit_each.hpp"

// External includes


// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"

namespace Kratos
{

class VariableNamePrinter : public boost::static_visitor<>
{
    public:
        VariableNamePrinter(std::stringstream& rString) : mrString(rString)
        {
        }

        template <class T>
        void operator()(const T& rVar)
        {
            mrString << rVar.Name() << std::endl;
        }

    private:
        std::stringstream& mrString;
};


///@addtogroup KratosCore
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

/// Generic Container of Variables
/** This class defines a generic container to hold variables of any type.
 * it is based on the variant/visitor paradigm.
 * 
 * the usage is something like
 * 
 *  GenericVariablesList vars;
 *  vars.AddVariable(TEMPERATURE);
 *  vars.AddVariable(VELOCITY);
 *  vars.AddVariable(DISPLACEMENT_X);
 *  vars.AddVariable(std::string("DISPLACEMENT_Z")); //note that we can add by the name
 *
 * operations are defined by creating a custom visitor with an operator specialized per every type
 * which is then used as
 *
 *  MyTestVisitor aux(mp);
 *   vars.ApplyVisitor( aux );
*/
class GenericVariablesList
{
public:
    ///@name Type Definitions
    ///@{
    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > variable_component_type;

    typedef boost::variant<
        Variable<bool>,
        Variable<int>,
        Variable<unsigned int>,        
        Variable<double>,
        Variable<array_1d<double,3>>,
        Variable<Vector>,
        Variable<Matrix>,
        variable_component_type        
    > variable_variant_type;

    using visitor_base_type = boost::static_visitor<> ;

    /// Pointer definition of GenericVariablesList
    KRATOS_CLASS_POINTER_DEFINITION(GenericVariablesList);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GenericVariablesList()
    {}

    /// Destructor.
    virtual ~GenericVariablesList()
    {}

    template<class TVarType>
    void AddVariable(const TVarType& rVar){
        mVariables.push_back(rVar);
    }

    void AddVariable(const std::string& rVarName);

    ///@}
    ///@name Operators
    ///@{
    template<class TVistorType >
    void ApplyVisitor(TVistorType& rVisitor) const
    {
        for(const auto& rVarVariant : mVariables)
        {
            boost::apply_visitor(rVisitor,rVarVariant);
        }
    }
    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "GenericVariableList\n";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        std::stringstream msg;
        msg << "GenericVariableList containing : " << std::endl;
        
        VariableNamePrinter aux(msg);
        this->ApplyVisitor(aux);
        rOStream << msg.str();
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    std::vector< variable_variant_type > mVariables;



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

    ///@}

}; // Class GenericVariablesList

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                GenericVariablesList& rThis)
                {return rIStream;}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const GenericVariablesList& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GENERIC_VARIABLES_LIST_H_INCLUDED  defined
