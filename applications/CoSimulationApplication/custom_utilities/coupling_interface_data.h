//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#if !defined(KRATOS_COUPLING_INTERFACE_DATA_H_INCLUDED)
#define KRATOS_COUPLING_INTERFACE_DATA_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "containers/model.h"
#include "utilities/auxiliar_model_part_utilities.h"


namespace Kratos
{
///@addtogroup CoSimulationApplication
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(CO_SIMULATION_APPLICATION) CouplingInterfaceData
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CouplingInterfaceData
    KRATOS_CLASS_POINTER_DEFINITION(CouplingInterfaceData);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CouplingInterfaceData(
        Parameters Settings,
        Model& rModel,
        const std::string& rName="default",
        const std::string& rSolverName="default_solver");

    ///@}

    /// Assignment operator.
    CouplingInterfaceData& operator=(CouplingInterfaceData const& rOther) = delete;

    /// Copy constructor.
    CouplingInterfaceData(CouplingInterfaceData const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    static Parameters GetDefaultParameters();

    ///@}
    ///@name Access
    ///@{

    void GetData() const;

    void SetData();

    ModelPart& GetModelPart();

    ModelPart& GetModelPart() const;

    ///@}
    ///@name Inquiry
    ///@{

    bool IsDistributed() const;

    std::size_t Size() const;

    std::size_t GetBufferSize() const;

    std::string GetName() const;

    std::string GetSolverName() const;

    std::string GetModelPartName() const;

    std::string GetVariableName() const;

    bool UsesScalarVariable() const;

    DataLocation GetDataLocation() const

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart* mpModelpart;

    std::string mVariableName;

    DataLocation mDataLocation;

    ///@}
    ///@name Private Operations
    ///@{


    ///@}

}; // Class CouplingInterfaceData

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const CouplingInterfaceData& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_COUPLING_INTERFACE_DATA_H_INCLUDED defined
