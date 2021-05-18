//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    
//

#if !defined(KRATOS_MESH_CONVERGENCE_ERROR_CALCULATOR_UTILITY)
#define KRATOS_MESH_CONVERGENCE_ERROR_CALCULATOR_UTILITY

// External includes
#include <omp.h>
#include <vector>

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "includes/kratos_parameters.h"
#include "includes/process_info.h"
#include "convection_diffusion_application_variables.h"

namespace Kratos
{
class MeshConvergenceErrorCalculator
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

    KRATOS_CLASS_POINTER_DEFINITION(MeshConvergenceErrorCalculator);

MeshConvergenceErrorCalculator()
{}

virtual ~MeshConvergenceErrorCalculator(){}

double GetL2VectorErrorNorm(ModelPart& r_model_part)
{
    double total_area = 0.0, result = 0.0;
    const unsigned int dim = r_model_part.GetProcessInfo()[DOMAIN_SIZE];

    for (const auto& r_node : r_model_part.Nodes())
    {
        double local_area = r_node.FastGetSolutionStepValue(NODAL_AREA);
        total_area += local_area;
        double error_x = r_node.FastGetSolutionStepValue(VECTORIAL_MESH_ERROR_X);
        double error_y = r_node.FastGetSolutionStepValue(VECTORIAL_MESH_ERROR_Y);
        double error_z = 0.0;
        if (dim == 3) error_z += r_node.FastGetSolutionStepValue(VECTORIAL_MESH_ERROR_Z);
        double squared_modulus = std::pow(error_x,2)+std::pow(error_y,2)+std::pow(error_z,2);
        result += squared_modulus * local_area;
    }

    return std::sqrt(result/total_area);
}


double GetL2ScalarErrorNorm(ModelPart& r_model_part)
{
    double total_area = 0.0, result = 0.0;

    for (const auto& r_node : r_model_part.Nodes())
    {
        double local_area = r_node.FastGetSolutionStepValue(NODAL_AREA);
        total_area += local_area;
        double error_p = r_node.FastGetSolutionStepValue(SCALAR_MESH_ERROR);
        result += std::pow(error_p,2) * local_area;
    }

    return std::sqrt(result/total_area);
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

}; // Class MeshConvergenceErrorCalculator

} // namespace Kratos.

#endif // KRATOS_MESH_CONVERGENCE_ERROR_CALCULATOR_UTILITY  defined