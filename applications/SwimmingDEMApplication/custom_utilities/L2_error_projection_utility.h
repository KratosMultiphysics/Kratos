//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//

#if !defined(KRATOS_L2_ERROR_PROJECTION_UTILITY)
#define KRATOS_L2_ERROR_PROJECTION_UTILITY

// External includes
#include <omp.h>
#include <vector>

// Project includes
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "includes/kratos_parameters.h"
#include "includes/process_info.h"

namespace Kratos
{
class L2ErrorProjection
{

public:

    typedef Node < 3 > NodeType;
    typedef Properties PropertiesType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    typedef ModelPart::ElementsContainerType::iterator  ElementIterator;
    typedef ModelPart::NodesContainerType::iterator     NodeIterator;
    typedef Kratos::Vector ShapeFunctionsType;

    KRATOS_CLASS_POINTER_DEFINITION(L2ErrorProjection);

L2ErrorProjection()
{}

virtual ~L2ErrorProjection(){}

double GetL2VectorProjection(ModelPart& r_model_part)
{
    const unsigned int n_elements = r_model_part.Elements().size();
    double interpolator = 0.0, result = 0.0, error_x = 0.0, error_y = 0.0, error_z = 0.0;
    const unsigned int dim = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
    Matrix NContainer;
    array_1d<double, 3> scalar_product;
    std::vector<double> error;

    for (unsigned int i = 0; i < n_elements; ++i){

        ElementIterator ielem = r_model_part.ElementsBegin() + i;
        GeometryType& rGeom = ielem->GetGeometry();
        const unsigned int NumNodes = rGeom.PointsNumber();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = ielem->GetGeometry().IntegrationPoints(GeometryData::GI_GAUSS_2);
        NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
        const SizeType& NumGauss = NContainer.size1();

        for (SizeType gss = 0; gss < NumGauss; ++gss){
            array_1d<double, 4> N;
            N = row(NContainer, gss);
            Vector DetJ = ZeroVector(NumGauss);
            ielem->GetGeometry().DeterminantOfJacobian(DetJ, GeometryData::GI_GAUSS_2);
            for (unsigned int j = 0; j < NumNodes; ++j){
                error_x += rGeom[j].FastGetSolutionStepValue(ERROR_X) * N[j];
                error_y += rGeom[j].FastGetSolutionStepValue(ERROR_Y) * N[j];
                if (dim == 3) error_z += rGeom[j].FastGetSolutionStepValue(ERROR_Z) * N[j];
            }
            error.push_back(error_x);
            error.push_back(error_y);
            if (dim == 3) error.push_back(error_z);
            // for (unsigned int d = 0; d < dim; ++d){
            //     scalar_product[d] = error[d];
            // }
            result += std::pow(SWIMMING_MODULUS_3(error),2) * DetJ[gss] * IntegrationPoints[gss].Weight();
            error.clear();
            error_x = 0.0;
            error_y = 0.0;
            if (dim == 3) error_z = 0.0;
        }
        interpolator += result;
        result = 0.0;
    }

    return std::sqrt(interpolator);
}


double GetL2ScalarProjection(ModelPart& r_model_part)
{
    const unsigned int n_elements = r_model_part.Elements().size();
    double interpolator = 0.0, result = 0.0, error = 0.0;
    Matrix NContainer;
    array_1d<double, 3> scalar_product;

    for (unsigned int i = 0; i < n_elements; ++i){

        ElementIterator ielem = r_model_part.ElementsBegin() + i;
        GeometryType& rGeom = ielem->GetGeometry();
        const unsigned int NumNodes = rGeom.PointsNumber();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = ielem->GetGeometry().IntegrationPoints(GeometryData::GI_GAUSS_2);
        NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
        const SizeType& NumGauss = NContainer.size1();

        for (SizeType gss = 0; gss < NumGauss; ++gss){
            array_1d<double, 4> N;
            N = row(NContainer, gss);
            Vector DetJ = ZeroVector(NumGauss);
            ielem->GetGeometry().DeterminantOfJacobian(DetJ, GeometryData::GI_GAUSS_2);
            for (unsigned int j = 0; j < NumNodes; ++j){
                error += rGeom[j].FastGetSolutionStepValue(SCALAR_ERROR) * N[j];
            }
            result += std::pow(error,2) * DetJ[gss] * IntegrationPoints[gss].Weight();
            error = 0.0;
        }
        interpolator += result;
        result = 0.0;
    }

    return std::sqrt(interpolator);
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

}; // Class L2ErrorProjection

} // namespace Kratos.

#endif // KRATOS_L2_ERROR_PROJECTION_UTILITY  defined