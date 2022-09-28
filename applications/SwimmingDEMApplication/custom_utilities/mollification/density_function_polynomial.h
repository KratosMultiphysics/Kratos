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
#include <cmath>

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
    m6 = fm6(mR);
    m4 = fm4(mR);
    m2 = fm2(mR);
    m0 = fm0(mR);
}

virtual ~DensityFunctionPolynomial(){}

void ComputeWeights(std::vector<double> & distances, std::vector<double> & nodal_areas, std::vector<double> & weights)
{
    double total_nodal_area_inv = 0.0;

    for (unsigned int i = 0; i != distances.size(); ++i){
        total_nodal_area_inv += nodal_areas[i];
    }

    total_nodal_area_inv = 1.0 / total_nodal_area_inv;

    double sum_of_weights_inv = 0.0;

    for (unsigned int i = 0; i != distances.size(); ++i){
        double radius_2 = distances[i] * distances[i];
        double weight = nodal_areas[i] * (m6 * std::pow(radius_2, 3) + m6 * m2 * radius_2 + m0);
        weights[i] = weight;
        sum_of_weights_inv += weight;
    }

    sum_of_weights_inv = 1.0 / sum_of_weights_inv;

    // normalizing weights
    for (unsigned int i = 0; i != distances.size(); ++i){
        weights[i] *= sum_of_weights_inv;
    }
}

void ComputeWeights(std::vector<double> & distances, std::vector<double> & nodal_areas, const double max_nodal_area_inv, std::vector<double> & weights)
{
    double total_nodal_area_inv = 0.0;

    for (unsigned int i = 0; i != distances.size(); ++i){
        total_nodal_area_inv += nodal_areas[i];
    }

    total_nodal_area_inv = 1.0 / total_nodal_area_inv;

    double sum_of_weights_inv = 0.0;

    for (unsigned int i = 0; i != distances.size(); ++i){
        double radius_2 = distances[i] * distances[i];
        double r = mR;// * pow(nodal_areas[i] * max_nodal_area_inv, 0.3333333333333333333333333333333);
        double weight;

        weight = nodal_areas[i] * ((radius_2 > r * r) ? 0.0 : m6 * std::pow(radius_2, 3) + m6 * m2 * radius_2 + m0);
        weights[i] = weight;
        sum_of_weights_inv += weight;
    }

    if(std::abs(sum_of_weights_inv) < std::numeric_limits<double>::epsilon()) sum_of_weights_inv = 0.0;
    else sum_of_weights_inv = 1.0 / sum_of_weights_inv;

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

double fm6(double r)
{
    return 7/(16*std::pow(mR, 7));
}

double fm4(double r)
{
    return -15*(std::pow(mR, 2) + 2)*std::exp(std::pow(mR, 2)/2)/(2*std::pow(mR, 4)*(2*std::pow(mR, 3)*std::exp(std::pow(mR, 2)/2) + 14*mR*std::exp(std::pow(mR, 2)/2) - 15*std::sqrt(2)*std::sqrt(Globals::Pi)*std::exp(std::pow(mR, 2))*erf(std::sqrt(2)*mR/2)));
}

double fm2(double r)
{
    return -3*std::pow(mR, 4);
}

double fm0(double r)
{
    return 7/(8*mR);
}

}; // class DensityFunctionPolynomial

} // namespace Kratos.

#endif // KRATOS_DENSITY_FUNCTION_POLYNOMIAL
