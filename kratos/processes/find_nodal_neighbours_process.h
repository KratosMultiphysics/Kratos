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
//


#if !defined(KRATOS_FIND_NODAL_NEIGHBOURS_PROCESS_H_INCLUDED )
#define  KRATOS_FIND_NODAL_NEIGHBOURS_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "processes/find_global_nodal_elemental_neighbours_process.h"

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
class FindNodalNeighboursProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindNodalNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindNodalNeighboursProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    FindNodalNeighboursProcess(ModelPart& model_part)
        :   mr_model_part(model_part)
    {
        auto& r_comm = model_part.GetCommunicator().GetDataCommunicator();
        mpNodeNeighboursCalculator = Kratos::make_unique<FindGlobalNodalNeighboursProcess>(r_comm, model_part);
        mpElemNeighboursCalculator = Kratos::make_unique<FindGlobalNodalElementalNeighboursProcess>(r_comm, model_part);

        KRATOS_WARNING("FindNodalNeighboursProcess") << 
            R"(please call separetely FindGlobalNodalNeighboursProcess 
            and FindGlobalNodalElementalNeighboursProcess. 
            The two calculations are currently independent,
             hence memory savings can be achieved)" << std::endl;
    }

    FindNodalNeighboursProcess(ModelPart& model_part, unsigned int avg_elems, unsigned int avg_nodes)
        : FindNodalNeighboursProcess(model_part)
    {
        KRATOS_WARNING("FindNodalNeighboursProcess") << "parameters avg_elems and avg_nodes are currently ignored. This constructor will be removed on the 2 of april 2020" << std::endl;
    }
    /// Destructor.
    ~FindNodalNeighboursProcess() override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        mpNodeNeighboursCalculator->Execute();
        mpElemNeighboursCalculator->Execute();
    }

    void ClearNeighbours()
    {
        mpNodeNeighboursCalculator->ClearNeighbours();
        mpElemNeighboursCalculator->ClearNeighbours();
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
    std::string Info() const override
    {
        return "FindNodalNeighboursProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindNodalNeighboursProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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


    ///@}
    ///@name Member Variables
    ///@{
    ModelPart& mr_model_part;

    std::unique_ptr<FindGlobalNodalElementalNeighboursProcess> mpElemNeighboursCalculator = nullptr;
    std::unique_ptr<FindGlobalNodalNeighboursProcess> mpNodeNeighboursCalculator = nullptr;
    ///@}
    ///@name Private Operators
    ///@{

    //******************************************************************************************
    //******************************************************************************************
    template< class TDataType > void  AddUniqueWeakPointer
    (GlobalPointersVector<TDataType>& v, const GlobalPointer<TDataType> candidate)
    {
        auto i = v.begin();
        auto endit = v.end();
        while ( i != endit && (i)->Id() != (candidate)->Id())
        {
            i++;
        }
        if( i == endit )
        {
            v.push_back(candidate);
        }

    }

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
    FindNodalNeighboursProcess& operator=(FindNodalNeighboursProcess const& rOther);

    /// Copy constructor.
    //FindNodalNeighboursProcess(FindNodalNeighboursProcess const& rOther);


    ///@}

}; // Class FindNodalNeighboursProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FindNodalNeighboursProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FindNodalNeighboursProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_FIND_NODAL_NEIGHBOURS_PROCESS_H_INCLUDED  defined 


