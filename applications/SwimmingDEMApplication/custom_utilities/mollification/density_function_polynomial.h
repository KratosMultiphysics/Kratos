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
#include <cstdlib>
#include "density_function.h"
#include <cmath>
#include "includes/model_part.h"
//Database includes
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/search/point_point_search.h"

#include "utilities/binbased_fast_point_locator.h"
#include "utilities/binbased_nodes_in_element_locator.h"

namespace Kratos
{

template< unsigned int TDim>
class DensityFunctionPolynomial: public DensityFunction<TDim>
{
public:
KRATOS_CLASS_POINTER_DEFINITION(DensityFunctionPolynomial);

DensityFunctionPolynomial(const double range)
    : mR(range)
{
    m6 = fm6(mR);
    m4 = fm4(mR);
    m2 = fm2(mR);
    m0 = fm0(mR);
}

virtual ~DensityFunctionPolynomial(){}

void ComputeWeights(std::vector<double> & distances, std::vector<double> & nodal_areas, std::vector<double> & weights, ModelPart::NodesContainerType::ContainerType neighbor, BinBasedFastPointLocator<TDim> bin_of_objects_fluid, array_1d<double, 3> p_coord)
{
    double total_nodal_area_inv = 0.0;
    double min_distance = mR;
    int n_neigh;
    bool on_boundary = false;
    double sum_of_weights_inv = 0.0;
    array_1d<double, 3> virtual_coordinates;

    for (unsigned int i = 0; i != distances.size(); ++i){
        double distance = 0.0;
        total_nodal_area_inv += nodal_areas[i];
        if (neighbor[i]->Is(BOUNDARY)){
            on_boundary = true;
            for (unsigned int d = 0; d < TDim; ++d)
                distance += std::pow(neighbor[i]->Coordinates()[d] - p_coord[d],2);
            distance = std::sqrt(distance);
            if (distance <= min_distance){
                min_distance = distance;
                n_neigh = i;
            }
        }
    }
    if (on_boundary == true){
        // Calculate intersection between the normal direction of plane from particle and the plane
        array_1d<double, 3> normal = -neighbor[n_neigh]->FastGetSolutionStepValue(NORMAL);

        array_1d<double, 3> intersection_coordinates;
        double t = (-normal[0] * (p_coord[0] - neighbor[n_neigh]->Coordinates()[0]) - normal[1] * (p_coord[1] - neighbor[n_neigh]->Coordinates()[1]) - normal[2] * (p_coord[2] - neighbor[n_neigh]->Coordinates()[2])) / (-std::pow(normal[0],2)-std::pow(normal[1],2)-std::pow(normal[2],2));
        for (unsigned int d = 0; d < TDim; ++d)
            intersection_coordinates[d] = p_coord[d] - normal[d]*t;

        // Calculate virtual coordinates
        for (unsigned int d = 0; d < TDim; ++d)
            virtual_coordinates[d] = 2.0 * intersection_coordinates[d] - p_coord[d];

    }

    for (unsigned int i = 0; i != distances.size(); ++i){
        double virtual_weight = 0.0;
        if (on_boundary == true){
            double virtual_distance = 0.0;
            for (unsigned int d = 0; d < TDim; ++d)
                virtual_distance += std::pow(virtual_coordinates[d] - neighbor[i]->Coordinates()[d],2);
            virtual_distance = std::sqrt(virtual_distance);
            if (virtual_distance <= mR){
                double virtual_radius_2 = virtual_distance * virtual_distance;
                virtual_weight = nodal_areas[i] * (m6 * std::pow(virtual_radius_2, 3) + m6 * m2 * virtual_radius_2 + m0);
            }
        }
        double radius_2 = distances[i] * distances[i];
        double weight = nodal_areas[i] * (m6 * std::pow(radius_2, 3) + m6 * m2 * radius_2 + m0);
        weights[i] = weight + virtual_weight;
        sum_of_weights_inv += weight + virtual_weight;
    }

    sum_of_weights_inv = 1.0 / sum_of_weights_inv;

    // normalizing weights
    for (unsigned int i = 0; i != distances.size(); ++i){
        weights[i] *= sum_of_weights_inv;
    }
}

private:

double mR; // the radius of the support of the PDF;
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
