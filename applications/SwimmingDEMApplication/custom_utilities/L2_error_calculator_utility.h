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

#if !defined(KRATOS_L2_ERROR_CALCULATOR_UTILITY)
#define KRATOS_L2_ERROR_CALCULATOR_UTILITY

// External includes
#include <omp.h>
#include <vector>

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "includes/kratos_parameters.h"
#include "includes/process_info.h"

namespace Kratos
{
class L2ErrorNormCalculator
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

    KRATOS_CLASS_POINTER_DEFINITION(L2ErrorNormCalculator);

L2ErrorNormCalculator()
{}

virtual ~L2ErrorNormCalculator(){}

void ComputeDofsErrors(ModelPart& r_model_part)
{
    block_for_each(r_model_part.Nodes(), [&](Node<3>& r_node)
    {
        auto& r_vectorial_error = r_node.FastGetSolutionStepValue(VECTORIAL_ERROR);
        r_vectorial_error = r_node.FastGetSolutionStepValue(VELOCITY) - r_node.FastGetSolutionStepValue(EXACT_VELOCITY);

        auto& r_scalar_error = r_node.FastGetSolutionStepValue(SCALAR_ERROR);
        r_scalar_error = r_node.FastGetSolutionStepValue(PRESSURE) - r_node.FastGetSolutionStepValue(EXACT_PRESSURE);
    });
}

double GetL2VectorErrorNorm(ModelPart& r_model_part)
{
    const unsigned int n_elements = r_model_part.Elements().size();
    double squared_modulus, total_area = 0.0, sum_error = 0.0, result = 0.0, error_x = 0.0, error_y = 0.0, error_z = 0.0;
    const unsigned int dim = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
    Matrix NContainer;

    for (unsigned int i = 0; i < n_elements; ++i){

        ElementIterator ielem = r_model_part.ElementsBegin() + i;
        GeometryType& rGeom = ielem->GetGeometry();
        const unsigned int NumNodes = rGeom.PointsNumber();
        const auto& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
        const SizeType& NumGauss = NContainer.size1();
        const double& area = rGeom.Area();
        total_area += area;

        for (SizeType gss = 0; gss < NumGauss; ++gss){
            array_1d<double, 4> N;
            N = row(NContainer, gss);
            Vector DetJ = ZeroVector(NumGauss);
            rGeom.DeterminantOfJacobian(DetJ, GeometryData::GI_GAUSS_2);
            for (unsigned int j = 0; j < NumNodes; ++j){
                error_x += rGeom[j].FastGetSolutionStepValue(ERROR_X) * N[j];
                error_y += rGeom[j].FastGetSolutionStepValue(ERROR_Y) * N[j];
                error_z += 0.0;
                if (dim == 3) error_z += rGeom[j].FastGetSolutionStepValue(ERROR_Z) * N[j];
            }
            squared_modulus = std::pow(error_x,2)+std::pow(error_y,2)+std::pow(error_z,2);
            result += squared_modulus * DetJ[gss] * IntegrationPoints[gss].Weight();
            error_x = 0.0;
            error_y = 0.0;
            error_z = 0.0;
        }
        sum_error += result;
        result = 0.0;
    }

    return std::sqrt(sum_error/total_area);
}


double GetL2ScalarErrorNorm(ModelPart& r_model_part)
{
    const unsigned int n_elements = r_model_part.Elements().size();
    double sum_error = 0.0, result = 0.0, error = 0.0, total_area = 0.0;
    Matrix NContainer;

    for (unsigned int i = 0; i < n_elements; ++i){

        ElementIterator ielem = r_model_part.ElementsBegin() + i;
        GeometryType& rGeom = ielem->GetGeometry();
        const unsigned int NumNodes = rGeom.PointsNumber();
        const auto& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
        const SizeType& NumGauss = NContainer.size1();
        const double& area = rGeom.Area();
        total_area += area;


        for (SizeType gss = 0; gss < NumGauss; ++gss){
            array_1d<double, 4> N;
            N = row(NContainer, gss);
            Vector DetJ = ZeroVector(NumGauss);
            rGeom.DeterminantOfJacobian(DetJ, GeometryData::GI_GAUSS_2);
            for (unsigned int j = 0; j < NumNodes; ++j){
                error += rGeom[j].FastGetSolutionStepValue(SCALAR_ERROR) * N[j];
            }
            result += std::pow(error,2) * DetJ[gss] * IntegrationPoints[gss].Weight();
            error = 0.0;
        }
        sum_error += result;
        result = 0.0;
    }

    return std::sqrt(sum_error/total_area);
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

}; // Class L2ErrorNormCalculator

} // namespace Kratos.

#endif // KRATOS_L2_ERROR_CALCULATOR_UTILITY  defined