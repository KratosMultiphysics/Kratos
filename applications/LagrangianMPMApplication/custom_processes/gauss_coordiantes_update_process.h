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

#if !defined(GAUSS_COORDINATES_UPDATE_UTILITY)
#define GAUSS_COORDINATES_UPDATE_UTILITY



/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/variables.h"

#include "custom_utilities/compute_mls_shape_functions_utility.h"
#include "lagrangian_mpm_application_variables.h"
#include "utilities/math_utils.h"
#include <cmath>
#include <algorithm>



namespace Kratos
{



class Gauss_Coordinates_Update_Process : public Process
{

public:


    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;



    Gauss_Coordinates_Update_Process( ModelPart& model_part ) : mr_model_part( model_part )
    {
       /////////////

    }

    ~Gauss_Coordinates_Update_Process() {}



    virtual void Execute()

    {

        KRATOS_TRY;


        for(ModelPart::ElementsContainerType::iterator iel = mr_model_part.ElementsBegin();
            iel!=mr_model_part.ElementsEnd(); iel++)
        {

            //const unsigned int number_of_nodes = iel->GetGeometry().size();
            //const double& h = iel->GetValue(EFFECTIVE_RADIUS);

/*
            const double h = iel->GetValue(EFFECTIVE_RADIUS);

            const array_1d<double,3>& xg = iel->GetValue(GAUSS_POINT_COORDINATES);



            Matrix TempCoordinates(number_of_nodes,2);
            for ( unsigned int i = 0; i < iel->GetGeometry().size(); i++ )
            {
                row(TempCoordinates,i) = iel->GetGeometry()[i].GetInitialPosition()+
                         iel->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

             }
             Vector Ng( number_of_nodes);
             Matrix DN_DX( number_of_nodes, 2 );

             LinearMLSKernel::ComputeMLSKernel(Ng, DN_DX, TempCoordinates, xg, h);
*/
             Vector& Ng = iel->GetValue(SHAPE_FUNCTIONS);

             for ( unsigned int i = 0; i < iel->GetGeometry().size(); i++ )
              {
                 Vector delta_D = ZeroVector(2);
                 delta_D = iel->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT)
                         - iel->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
                 iel->GetValue(GAUSS_POINT_COORDINATES) += Ng[i] * delta_D;
                 //iel->GetValue(TEMP_POS) += Ng[i] * row(TempCoordinates,i);
                         //TempCoordinates(i,0);


              }
             //KRATOS_WATCH(Ng);
             //KRATOS_WATCH(iel->GetValue(DomainSize));
             //KRATOS_WATCH(iel->GetValue(TEMP_POS));

        }



        KRATOS_CATCH("");

    }





private:

    ModelPart& mr_model_part;


};
}

#endif

