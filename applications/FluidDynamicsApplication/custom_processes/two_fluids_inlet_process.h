//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Simon Wenczowski
//
//

#ifndef KRATOS_TWO_FLUIDS_INLET_PROCESS_H
#define KRATOS_TWO_FLUIDS_INLET_PROCESS_H

// System includes

// External includes

// Project includes
#include "processes/process.h"

// Application includes


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

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) TwoFluidsInletProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DistanceModificationProcess
    KRATOS_CLASS_POINTER_DEFINITION(TwoFluidsInletProcess);

    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Kratos parameters.

    /**
     * @brief Construct a new Two Fluids Inlet Process object
     * The constructor calculates the AUX_DISTANCE field that is needed for the smoothing of the distance field.
     * Also, it subdivides the inlet into to sub model parts for water and air.
     *
     * @param rModelPart Model part the inlet process is applied to
     * @param rParameters Checked parameters (no more checking mechanism applied)
     */
    TwoFluidsInletProcess(
        ModelPart& rModelPart,
        Parameters& rParameters,
        Process::Pointer pDistanceProcess );

    /// Destructor.
    ~TwoFluidsInletProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Function to smooth the distance field.
     * A continuous transition from the inlet distance the the distance field inside the domain shall be achieved.
     */
    void SmoothDistanceField();

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
        std::stringstream buffer;
        buffer << "TwoFluidsInletProcess";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "TwoFluidsInletProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


    ///@}
    ///@name Friends
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrInletModelPart;

    array_1d<double,3> mInterfaceNormal;

    array_1d<double,3> mInterfacePoint;

    double mInletRadius;

    ///@}
    ///@name Protected Operators
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

}; // Class DistanceModificationProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_TWO_FLUIDS_INLET_PROCESS_H
