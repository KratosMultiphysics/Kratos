/*
==============================================================================
KratosPFEMApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:31 $
//   Revision:            $Revision: 1.2 $
//
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

/// Short class definition.
/** Detail class definition.
	calculate the nodal H for all the nodes depending on the min distance
	of the neighbouring nodes.

	lonely nodes are given the average value of the H
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
    }

    /// Destructor.
    virtual ~ComputeBDFCoefficientsProcess()
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
    ///this function fills the vector BDF_COEFFICIENTS with the correct values
    ///depending on the time_order chosen
    virtual void Execute() override
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
    virtual std::string Info() const override
    {
        return "ComputeBDFCoefficientsProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeBDFCoefficientsProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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

    /// Assignment operator.
    ComputeBDFCoefficientsProcess& operator=(ComputeBDFCoefficientsProcess const& rOther);

    /// Copy constructor.
    //ComputeBDFCoefficientsProcess(ComputeBDFCoefficientsProcess const& rOther);


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


