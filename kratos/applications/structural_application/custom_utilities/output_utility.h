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
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2008-10-23 13:03:25 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

#if !defined(KRATOS_OUTPUT_UTILITY_INCLUDED )
#define  KRATOS_OUTPUT_UTILITY_INCLUDED
//System includes
//External includes
#include "boost/smart_ptr.hpp"

//Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "includes/element.h"
#include "integration/integration_point.h"
#include "geometries/geometry.h"
#include "spaces/ublas_space.h"

namespace Kratos
{
    class OutputUtility
    {
        public:
            typedef Geometry<Node<3> >::IntegrationPointsArrayType IntegrationPointsArrayType;
            typedef Geometry<Node<3> >::GeometryType GeometryType;
            typedef Geometry<Node<3> >::CoordinatesArrayType CoordinatesArrayType;
            
			/** 
             * Constructor.
                         */
            OutputUtility()
            {
                std::cout << "OutputUtility created" << std::endl;
            }
			
			/** 
             * Destructor.
                         */
            virtual ~OutputUtility()
            {}
            
            Vector GetStrain( ModelPart& rModelPart, Element& element, int gp_index )
            {
                std::vector<Matrix> output;
                Vector result = ZeroVector(6);
                element.CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_TENSOR, output, 
                        rModelPart.GetProcessInfo() );
                for( int i=0; i<6; i++ )
                    result[i] = output[gp_index](0,i);
                return( result );
            }
            
            Vector GetStress( ModelPart& rModelPart, Element& element, int gp_index )
            {
                std::vector<Matrix> output;
                Vector result = ZeroVector(6);
                element.CalculateOnIntegrationPoints(PK2_STRESS_TENSOR, output, 
                        rModelPart.GetProcessInfo() );
                for( int i=0; i<6; i++ )
                    result[i] = output[gp_index](0,i);
                return( result );
            }
            
            Vector GetInternalVariables( ModelPart& rModelPart, Element& element, int gp_index )
            {
                std::vector<Vector> output;
                element.CalculateOnIntegrationPoints(INTERNAL_VARIABLES, output, 
                        rModelPart.GetProcessInfo() );
                return( output[gp_index] );
            }
            
            
    };//Class OutputUtility
}//namespace Kratos.

#endif /* KRATOS_OUTPUT_UTILITY_INCLUDED defined */
