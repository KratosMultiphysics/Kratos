/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
/* *********************************************************
*
*   Last Modified by:    $Author:  $
*   Date:                $Date: 01-27-2010$
*   Revision:            $Revision: 1.00   $
*
* ***********************************************************/

#if !defined(NODAL_VALUES_UTILITY)
#define NODAL_VALUES_UTILITY



/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "containers/array_1d.h"

#include "lagrangian_mpm_application.h"
#include "processes/find_nodal_neighbours_process.h"
#include "utilities/math_utils.h"
//#include "custom_utilities/sd_math_utils.h"
#include <cmath>
#include <algorithm>



namespace Kratos
{



class Nodal_Values_Utility
{

public:


    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef boost::numeric::ublas::vector<Matrix> Matrix_Order_Tensor;
    typedef boost::numeric::ublas::vector<Vector> Vector_Order_Tensor;






    Nodal_Values_Utility( ModelPart& model_part, int domain_size ) : mr_model_part( model_part )
    {
        mdomain_size = domain_size;

    }

    ~Nodal_Values_Utility() {}



    virtual void CalculateNodalVarialbes( const Variable<Matrix >& rVariable)

    {

        KRATOS_TRY;



         //const Variable<double >& rVariable;
        Matrix Output;
        const ProcessInfo& rCurrentProcessInfo = mr_model_part.GetProcessInfo();


        if(rVariable == NODAL_STRESS)
        {




            for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode!=mr_model_part.NodesEnd(); inode++)
            {

                inode->FastGetSolutionStepValue(NODAL_STRESS) = ZeroMatrix(1,3);
                //inode->FastGetSolutionStepValue(NODAL_STRAIN) = ZeroMatrix(1,3);


            }


            for(ModelPart::ElementsContainerType::iterator iel = mr_model_part.ElementsBegin();
                iel!=mr_model_part.ElementsEnd(); iel++)
            {


                 iel->Calculate(rVariable, Output, rCurrentProcessInfo);





            }


            for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode!=mr_model_part.NodesEnd(); inode++)
            {

                //std::cout << inode->FastGetSolutionStepValue(NODAL_STRAIN) << std::endl;

                const double area = inode->FastGetSolutionStepValue(NODAL_AREA);
                //inode->FastGetSolutionStepValue(NODAL_STRAIN) /= area;
                inode->FastGetSolutionStepValue(NODAL_STRESS) /= area;

                //KRATOS_WATCH(inode->FastGetSolutionStepValue(NODAL_AREA));
                //KRATOS_WATCH(inode->FastGetSolutionStepValue(NODAL_STRESS));


            }




        }

        else if(rVariable == NODAL_STRAIN)
        {




            for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode!=mr_model_part.NodesEnd(); inode++)
            {


                inode->FastGetSolutionStepValue(NODAL_STRAIN) = ZeroMatrix(1,3);


            }


            for(ModelPart::ElementsContainerType::iterator iel = mr_model_part.ElementsBegin();
                iel!=mr_model_part.ElementsEnd(); iel++)
            {


                 iel->Calculate(rVariable, Output, rCurrentProcessInfo);



            }


            for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode!=mr_model_part.NodesEnd(); inode++)
            {

                //std::cout << inode->FastGetSolutionStepValue(NODAL_STRAIN) << std::endl;

                const double area = inode->FastGetSolutionStepValue(NODAL_AREA);
                inode->FastGetSolutionStepValue(NODAL_STRAIN) /= area;




            }




        }




        KRATOS_CATCH("");

    }


    virtual void CalculateNodalArea(const Variable<double >& rVariable)

    {

        KRATOS_TRY;


        //Process& neighbour_finder = FindNodalNeighboursProcess(mr_model_part,10, 10);
        //neighbour_finder.Execute();

        FindNodalNeighboursProcess  NodosVecinos(mr_model_part, 10, 10);
        NodosVecinos.Execute();
         //const Variable<double >& rVariable;
        double Output;
        const ProcessInfo& rCurrentProcessInfo = mr_model_part.GetProcessInfo();

        if(rVariable == NODAL_AREA)
        {


            for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode!=mr_model_part.NodesEnd(); inode++)
            {


                inode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;

            }


            for(ModelPart::ElementsContainerType::iterator iel = mr_model_part.ElementsBegin();
                iel!=mr_model_part.ElementsEnd(); iel++)
            {


                 iel->Calculate(rVariable, Output, rCurrentProcessInfo);



            }

            for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode!=mr_model_part.NodesEnd(); inode++)
            {

                //std::cout << inode->FastGetSolutionStepValue(NODAL_STRAIN) << std::endl;

//                 const double area = inode->FastGetSolutionStepValue(NODAL_AREA);
               // inode->FastGetSolutionStepValue(NODAL_H) = sqrt(area);




            }



        }

/*
        if(rVariable == NODE_EQUIVALENT_PLASTIC_STRAIN)
        {


            for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode!=mr_model_part.NodesEnd(); inode++)
            {


                inode->FastGetSolutionStepValue(NODE_EQUIVALENT_PLASTIC_STRAIN) = 0.0;
                //inode->FastGetSolutionStepValue(PLASTIC_STRAIN) = 0.0;

            }


            for(ModelPart::ElementsContainerType::iterator iel = mr_model_part.ElementsBegin();
                iel!=mr_model_part.ElementsEnd(); iel++)
            {


                 iel->Calculate(rVariable, Output, rCurrentProcessInfo);



            }

            for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode!=mr_model_part.NodesEnd(); inode++)
            {

                //std::cout << inode->FastGetSolutionStepValue(NODAL_STRAIN) << std::endl;

                //const double area = inode->FastGetSolutionStepValue(NODAL_AREA);
               // inode->FastGetSolutionStepValue(NODAL_H) = sqrt(area);




            }



        }



        else if ( rVariable == DAMAGE )
        {

/*

            for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode!=mr_model_part.NodesEnd(); inode++)
            {


                inode->FastGetSolutionStepValue(DAMAGE) = 0;

            }


            for(ModelPart::ElementsContainerType::iterator iel = mr_model_part.ElementsBegin();
                iel!=mr_model_part.ElementsEnd(); iel++)
            {


                 iel->Calculate(rVariable, Output, rCurrentProcessInfo);



            }



            for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
                inode!=mr_model_part.NodesEnd(); inode++)
            {

                //std::cout << inode->FastGetSolutionStepValue(NODAL_STRAIN) << std::endl;

                WeakPointerVector< Element >& rneigh_el = inode->GetValue(NEIGHBOUR_ELEMENTS);
                //const double area = inode->FastGetSolutionStepValue(NODAL_AREA);
                inode->FastGetSolutionStepValue(DAMAGE) /= rneigh_el.size();

               //KRATOS_WATCH(inode->FastGetSolutionStepValue(DAMAGE));


            }


        }*/


         KRATOS_CATCH("");


      }


private:

    ModelPart& mr_model_part;
    unsigned int mdomain_size;

};
}

#endif

