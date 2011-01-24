/*
==============================================================================
KratosIncompressibleFluidApplication
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
/*
 * File:   estimate_dt_utilities.h
 * Author: jcotela
 *
 * Created on January 18, 2011, 6:45 PM
 */

#ifndef KRATOS_ESTIMATE_DT_UTILITIES_H
#define	KRATOS_ESTIMATE_DT_UTILITIES_H

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"
#include "utilities/geometry_utilities.h"


namespace Kratos {

    ///@name Kratos Classes
    ///@{

    /// Estimate the time step in a fluid problem to obtain a given Courant number.
    template< unsigned int TDim >
    class EstimateDtUtil
    {
    public:

        ///@name Life Cycle
        ///@{

        /// Constructor
        /**
         @param rModelPart The model part containing the problem mesh
         */
        EstimateDtUtil(ModelPart& rModelPart):
            mrModelPart(rModelPart)
        {}

        /// Destructor
        ~EstimateDtUtil()
        {}
        
        ///@}
        ///@name Operations
        ///@{

        /// Calculate the maximum time step that satisfies the Courant-Friedrichs-Lewy (CFL) condition.
        /**
         @param CFL The upper limit for the Courant number
         @param dt_max Maximum admissible time step (upper bound to be used for situations with very low velocity fields)
         */
        double EstimateDt(double CFL, double dt_max)
        {
            KRATOS_TRY;

            const unsigned int NumNodes = TDim +1;

            double Area;
            array_1d<double, NumNodes> N;
            boost::numeric::ublas::bounded_matrix<double, NumNodes, TDim> DN_DX;

            double MaxDen = 0.0;
            double Den;

            for( ModelPart::ElementIterator itElem = mrModelPart.ElementsBegin(); itElem != mrModelPart.ElementsEnd(); ++itElem)
            {
                // Get the element's geometric parameters
                Geometry< Node<3> >& rGeom = itElem->GetGeometry();
                GeometryUtils::CalculateGeometryData(rGeom, DN_DX, N, Area);
                Den = 0;

                for(unsigned int i = 0; i < NumNodes; ++i)
                {
                    array_1d< double, 3 >& rVel = rGeom[i].FastGetSolutionStepValue(VELOCITY);

                    for(unsigned int j = 0; j < NumNodes; ++j)
                    {
                        double temp = 0;
                        for(unsigned int d = 0; d < TDim; ++d)
                        {
                            temp = fabs(rVel[d] * DN_DX(j,d));
                            if(temp > Den) Den = temp;
                        }
                    }
                }
                if(Den > MaxDen)
                    MaxDen = Den;
            }

            double dt = CFL / MaxDen;
            if(dt > dt_max)
                dt = dt_max;

            return dt;

            KRATOS_CATCH("")
        }
        
        ///@} // Operators

    private:
        
        ///@name Member Variables
        ///@{

        /// A reference to the problem's model part
        ModelPart& mrModelPart;
        
        ///@} // Member variables

    };

    ///@} // Kratos classes

} // namespace Kratos.


#endif	/* KRATOS_ESTIMATE_DT_UTILITIES_H */

