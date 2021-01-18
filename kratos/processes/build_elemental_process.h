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


#if !defined(KRATOS_BUILD_ELEMENTAL_PROCESS_H_INCLUDED )
#define  KRATOS_BUILD_ELEMENTAL_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "containers/pointer_vector_set.h"


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
template<class TEquationSystemType>
class BuildElementalProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BuildElementalProcess
    KRATOS_CLASS_POINTER_DEFINITION(BuildElementalProcess);

    typedef Element ElementType;

    typedef PointerVectorSet<ElementType, IndexedObject> ElementsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    BuildElementalProcess(ElementsContainerType& rElements, TEquationSystemType & rEquationSystem, ProcessInfo& rProcessInfo)
        : mElements(rElements), mEquationSystem(rEquationSystem), mProcessInfo(rProcessInfo) {}

    /// Destructor.
    virtual ~BuildElementalProcess() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void Execute()
    {
        ElementType::VectorType rhs = ZeroVector(3);

        for(IO::ElementsContainerType::iterator i_element = mElements.begin() ; i_element != mElements.end() ; i_element++)
        {
            ElementType::MatrixType stiffness_matrix = ZeroMatrix(3,3);

            Element::EquationIdVectorType ids(3,0);

            i_element->EquationIdVector(ids, mProcessInfo);

            i_element->CalculateLocalSystem(stiffness_matrix, rhs, mProcessInfo);

            mEquationSystem.Assemble(stiffness_matrix, rhs, ids);
        }
    }


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
        return "BuildElementalProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "BuildElementalProcess";
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

    ElementsContainerType& mElements;
    TEquationSystemType& mEquationSystem;
    ProcessInfo& mProcessInfo;

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

    /// Assignment operator.
    BuildElementalProcess& operator=(BuildElementalProcess const& rOther);

    /// Copy constructor.
    BuildElementalProcess(BuildElementalProcess const& rOther);


    ///@}

}; // Class BuildElementalProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function

template<class TEquationSystemType>
inline std::istream& operator >> (std::istream& rIStream,
                                  BuildElementalProcess<TEquationSystemType>& rThis);

/// output stream function
template<class TEquationSystemType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BuildElementalProcess<TEquationSystemType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_BUILD_ELEMENTAL_PROCESS_H_INCLUDED  defined 


