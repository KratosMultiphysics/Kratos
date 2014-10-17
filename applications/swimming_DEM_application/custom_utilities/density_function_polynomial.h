/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last Modified by:    $Author: gcasas $
//   Date:                $Date: 2014-10-09 10:34:00 $
//   Revision:            $Revision: 0.1 $
//
//

#if !defined(KRATOS_DENSITY_FUNCTION_POLYNOMIAL)
#define KRATOS_DENSITY_FUNCTION_POLYNOMIAL
// System includes
#include <string>
#include <iostream>
#include <stdlib.h>
#include "density_function.h"

namespace Kratos
{

template< unsigned int TDim>
class DensityFunctionPolynomial: public DensityFunction<TDim>
{
public:
KRATOS_CLASS_POINTER_DEFINITION(DensityFunctionPolynomial);

DensityFunctionPolynomial(const double range, const double shape_factor)
    : mR(range),
      mHeightOverR(shape_factor)
{
    m6 = 315 / (32 * KRATOS_M_PI * pow(mR, 9)) - 3 * mHeightOverR / pow(mR, 5);
    m4 = 7 * mHeightOverR / pow(mR, 3) - 315 / (16 * KRATOS_M_PI * pow(mR, 7));
    m2 = 315 / (32 * KRATOS_M_PI * pow(mR, 5)) - 5 * mHeightOverR / mR;
    m0 = mR * mHeightOverR;
}

virtual ~DensityFunctionPolynomial(){}

void ComputeWeights(std::vector<double> & distances, std::vector<double> & weights)
{
    double sum_of_weights_inv = 0.0;

    for (unsigned int i = 0; i != distances.size(); ++i){
        double radius_2 = distances[i] * distances[i];
        double weight = m6 * radius_2 * radius_2 * radius_2 + m4 * radius_2 * radius_2 + m2 * radius_2 + m0;
        weights[i] = weight;
        sum_of_weights_inv += weight;
    }

    sum_of_weights_inv = 1.0 / sum_of_weights_inv;

    // normalizing weights
    for (unsigned int i = 0; i != distances.size(); ++i){
        weights[i] *= sum_of_weights_inv;
    }
}

private:

double mR; // the radius of the support of the PDF;
double mHeightOverR;
double m6;
double m4;
double m2;
double m0;

}; // class DensityFunctionPolynomial

} // namespace Kratos.

#endif // KRATOS_DENSITY_FUNCTION_POLYNOMIAL
