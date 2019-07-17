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



#if !defined(KRATOS_COMPUTE_BDF_COEFFICIENTS_PROCESS_INCLUDED )
#define  KRATOS_COMPUTE_BDF_COEFFICIENTS_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <limits>
// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"


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

/**
 * @brief Auxiliary class to compute the BDF coefficients
 * This class computes the BDF coefficients for the time step values stored
 * in the ProcessInfo. It is valid for 1st and 2nd order BDF schemes and for
 * non-constant delta time values.
 */
class ComputeBDFCoefficientsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeBDFCoefficientsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeBDFCoefficientsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ComputeBDFCoefficientsProcess(ModelPart& model_part, unsigned int time_order)
        : mr_model_part(model_part), mtime_order(time_order)
    {
        KRATOS_WARNING_ONCE("DEPRECATION") << "\"ComputeBDFCoefficientsProcess\" is deprecated, please use the functionalities provided by /kratos/utilities/time_discretization.h" << std::endl;
        KRATOS_ERROR_IF(mtime_order == 0 || mtime_order > 2) << "Time order must be either \'1\' or \'2\'. Got " << mtime_order << std::endl;
        mr_model_part.GetProcessInfo()[BDF_COEFFICIENTS] = ZeroVector(mtime_order + 1);
    }

    /// Destructor.
    ~ComputeBDFCoefficientsProcess() = default;

    /// Assignment operator.
    ComputeBDFCoefficientsProcess &operator=(ComputeBDFCoefficientsProcess const &rOther) = delete;

    /// Copy constructor.
    ComputeBDFCoefficientsProcess(ComputeBDFCoefficientsProcess const& rOther) = delete;

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
    ///this function fills the vector BDF_COEFFICIENTS with the correct values
    ///depending on the time_order chosen
    void Execute() override
    {
        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo = mr_model_part.GetProcessInfo();

        if (mtime_order == 2)
        {
            //calculate the BDF coefficients
            double Dt = rCurrentProcessInfo[DELTA_TIME];
            double OldDt = rCurrentProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME];

            if(OldDt > 1e-10*Dt) //this should always be the case!!
            {
                double Rho = OldDt / Dt;
                double TimeCoeff = 1.0 / (Dt * Rho * Rho + Dt * Rho);

                Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
                BDFcoeffs.resize(3, false);

                BDFcoeffs[0] = TimeCoeff * (Rho * Rho + 2.0 * Rho); //coefficient for step n+1 (3/2Dt if Dt is constant)
                BDFcoeffs[1] = -TimeCoeff * (Rho * Rho + 2.0 * Rho + 1.0); //coefficient for step n (-4/2Dt if Dt is constant)
                BDFcoeffs[2] = TimeCoeff; //coefficient for step n-1 (1/2Dt if Dt is constant)
            }
            else
            {
                KRATOS_WATCH("Dt was zero at the previous time step!!!")
                Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
                BDFcoeffs.resize(3, false);

                BDFcoeffs[0] = 1.0/Dt; //coefficient for step n+1 (3/2Dt if Dt is constant)
                BDFcoeffs[1] = -1.0/Dt; //coefficient for step n (-4/2Dt if Dt is constant)
                BDFcoeffs[2] = 0.0; //coefficient for step n-1 (1/2Dt if Dt is constant)
            }
        }
        else if (mtime_order == 1)
        {
            double Dt = rCurrentProcessInfo[DELTA_TIME];
            double TimeCoeff = 1.0 / Dt;

            Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
            BDFcoeffs.resize(2, false);

            BDFcoeffs[0] = TimeCoeff; //coefficient for step n+1 (1/Dt)
            BDFcoeffs[1] = -TimeCoeff; //coefficient for step n (-1/Dt)
        }


        KRATOS_CATCH("")
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
        return "ComputeBDFCoefficientsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeBDFCoefficientsProcess";
    }

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

    ModelPart& mr_model_part;
    const unsigned int mtime_order;

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
}; // Class ComputeBDFCoefficientsProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ComputeBDFCoefficientsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ComputeBDFCoefficientsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_COMPUTE_BDF_COEFFICIENTS_PROCESS_INCLUDED  defined
