/*
==============================================================================
KratosULFApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
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
//   Date:                $Date: 2007-03-06 10:30:32 $
//   Revision:            $Revision: 1.2 $
//
//  this process marks distorted fluid elements and their direct neighbours

#if !defined(KRATOS_MARK_BAD_ELEMENTS_PROCESS_INCLUDED )
#define  KRATOS_MARK_BAD_ELEMENTS_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "custom_utilities/geometry_utilities2D.h"
#include "custom_elements/updated_lagrangian_fluid.h"


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
	Update the PRESSURE_FORCE on the nodes


*/

class MarkBadElementsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(MarkBadElementsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MarkBadElementsProcess(ModelPart& model_part)
        : mr_model_part(model_part)
    {
    }

    /// Destructor.
    ~MarkBadElementsProcess() override
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
        KRATOS_TRY
        for(ModelPart::ElementsContainerType::const_iterator im = mr_model_part.ElementsBegin(); im!=mr_model_part.ElementsEnd(); im++)
        {
            double x0 = im.X();
            double x1 = im.X();
            double x2 = pgeom[2].X();

            double y0 = pgeom[0].Y();
            double y1 = pgeom[1].Y();
            double y2 = pgeom[2].Y();

            /*if( ((y0<-0.1) || (y1<-0.1 ) || (y2<-0.1 )) && (x0>1.1 || x1>1.1 || x2>1.1 ))
            {
            	return false;
            }*/

            msJ(0,0)=2.0*(x1-x0);
            msJ(0,1)=2.0*(y1-y0);
            msJ(1,0)=2.0*(x2-x0);
            msJ(1,1)=2.0*(y2-y0);


            double detJ = msJ(0,0)*msJ(1,1)-msJ(0,1)*msJ(1,0);

            msJinv(0,0) =  msJ(1,1);
            msJinv(0,1) = -msJ(0,1);
            msJinv(1,0) = -msJ(1,0);
            msJinv(1,1) =  msJ(0,0);

            BoundedMatrix<double,2,2> check;


            if(detJ < 1e-12)
            {
                //std::cout << "detJ = " << detJ << std::endl;
                ////mark as boundary
                pgeom[0].GetSolutionStepValue(IS_BOUNDARY) = 1;
                pgeom[1].GetSolutionStepValue(IS_BOUNDARY) = 1;
                pgeom[2].GetSolutionStepValue(IS_BOUNDARY) = 1;
                return false;
            }

            else
            {

                double x0_2 = x0*x0 + y0*y0;
                double x1_2 = x1*x1 + y1*y1;
                double x2_2 = x2*x2 + y2*y2;

                //finalizing the calculation of the inverted matrix
                //std::cout<<"MATR INV"<<MatrixInverse(msJ)<<std::endl;
                msJinv /= detJ;
                //calculating the RHS
                ms_rhs[0] = (x1_2 - x0_2);
                ms_rhs[1] = (x2_2 - x0_2);

                //calculate position of the center
                noalias(msc) = prod(msJinv,ms_rhs);

                double radius = sqrt(pow(msc[0]-x0,2)+pow(msc[1]-y0,2));

                //calculate average h
                double h;
                h =  pgeom[0].FastGetSolutionStepValue(NODAL_H);
                h += pgeom[1].FastGetSolutionStepValue(NODAL_H);
                h += pgeom[2].FastGetSolutionStepValue(NODAL_H);
                h *= 0.333333333;
                if (radius < h*alpha_param)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }

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
        return "MarkBadElementsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MarkBadElementsProcess";
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
//		MarkBadElementsProcess& operator=(MarkBadElementsProcess const& rOther);

    /// Copy constructor.
//		MarkBadElementsProcess(MarkBadElementsProcess const& rOther);


    ///@}

}; // Class MarkBadElementsProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MarkBadElementsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MarkBadElementsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MARK_BAD_ELEMENTS_PROCESS_INCLUDED  defined 


