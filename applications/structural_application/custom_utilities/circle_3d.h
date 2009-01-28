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
*   Date:                $Date: 2008-04-29 11:41:31 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/

#if !defined(KRATOS_CIRCLE_3D_INCLUDED )
#define  KRATOS_CIRCLE_3D_INCLUDED
//System includes
//External includes
#include "boost/smart_ptr.hpp"
#include <cmath>

//Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
    class Circle3D
    {
        public:
            typedef Dof<double> TDofType;
            typedef PointerVectorSet<TDofType, IndexedObject> DofsArrayType;
            typedef ModelPart::ElementsContainerType ElementsArrayType;
            typedef double* ContainerType;
            typedef Element::DofsVectorType DofsVectorType;
            typedef Geometry<Node<3> >::IntegrationPointsArrayType IntegrationPointsArrayType;
            typedef Geometry<Node<3> >::GeometryType GeometryType;
            typedef Geometry<Node<3> >::CoordinatesArrayType CoordinatesArrayType;
            typedef UblasSpace<double, CompressedMatrix, Vector> SpaceType;
            
            /** 
             * Constructor.
             */
            Circle3D( Point<3> M, Point<3> e1, Point<3> e2 )
            :mM(M),mE1(e1),mE2(e2)
            {
                r1 = mE1-M;
                r2 = mE2-M;
                radius = MathUtils<double>::Norm(r1);
                radius2 = radius*radius;
                std::cout << "Circle created" << std::endl;
                KRATOS_WATCH( mM );
                KRATOS_WATCH( mE1 );
                KRATOS_WATCH( mE2 );
                KRATOS_WATCH( r1 );
                KRATOS_WATCH( r2 );
                KRATOS_WATCH( radius );
            }
            
            /**
             * 
             */
            Point<3>& GetCenter()
            {
                return mM;
            }
            
            /**
             * 
             */
            Point<3>& GetPoint( Point<3>& result, double s )
            {
                result = cos( s/radius )*r1 + sin( s/radius )*r2 + mM;
                return result;
            }
            
            /**
             * 
             */
            array_1d<double,3>& GetPoint1()
            {
                return mE1;
            }
            
            /**
             * 
             */
            array_1d<double,3>& GetPoint2()
            {
                return mE2;
            }
            
            /**
             * 
             */
            array_1d<double, 3> GetR1()
            {
                return r1;
            }
            
            /**
             * 
             */
            array_1d<double, 3> GetR2()
            {
                return r2;
            }
            
            /**
             * 
             */
            array_1d<double,3>& GetDerivative( array_1d<double,3>& result, double s )
            {
//                 result = r2-r1;
                result = (1.0/radius)*(-sin(s/radius)*r1+cos(s/radius)*r2);
                return result;
            }
            
            /**
             * 
             */
            double GetRadius()
            {
                return radius;
            }
            
            /**
             * 
             */
            double GetSqRadius()
            {
                return radius2;
            }
            
        private:
            
            Point<3> mM;
            Point<3> mE1;
            Point<3> mE2;
            array_1d<double, 3> r1;
            array_1d<double, 3> r2;
            double radius;
            double radius2;
    };//Class Circle3D
}//namespace Kratos.

#endif /* KRATOS_CIRCLE_3D_INCLUDED  defined */
