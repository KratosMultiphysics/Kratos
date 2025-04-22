//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aniol Sala
//

#if !defined(KRATOS_ERROR_ETHIER_FIELD)
#define KRATOS_ERROR_ETHIER_FIELD

// External includes
#include <omp.h>
#include <vector>

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "includes/kratos_parameters.h"
#include "includes/process_info.h"
#include "custom_utilities/fields/ethier_flow_field.h"

// Application includes
#include "swimming_DEM_application.h"
#include "swimming_dem_application_variables.h"

// Other applications includes
#include "fluid_dynamics_application_variables.h"
#include "fluid_dynamics_application.h"

namespace Kratos
{

class ErrorNormEthierFieldCalculator
{

public:

    typedef Node NodeType;
    typedef Properties PropertiesType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    typedef ModelPart::ElementsContainerType::iterator  ElementIterator;
    typedef ModelPart::NodesContainerType::iterator     NodeIterator;
    typedef Kratos::Vector ShapeFunctionsType;
    typedef Element::Pointer  ElementPointerType;
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    KRATOS_CLASS_POINTER_DEFINITION(ErrorNormEthierFieldCalculator);

    ErrorNormEthierFieldCalculator(const double a, const double b, const bool normalize_result): mA(a), mD(b), mNormalizeResult(normalize_result), mFlowField(a, b)
    {}

    ErrorNormEthierFieldCalculator(const double a, const double b): mA(a), mD(b), mNormalizeResult(false), mFlowField(a, b)
    {}

    double getL2NormFluidAccelWithoutRecoveryUsingGaussInterpolatedValues(ModelPart& r_model_part);   // Get the L2 norm using the shape functions (no recovered field)
    double getL2NormFluidAccelWithRecoveryUsingGaussInterpolatedValues(ModelPart& r_model_part);  // Get the L2 norm using the recovered field at the nodes, and then interpolated at the gauss points
    double getL2NormFluidAccelWithoutRecoveryUsingGaussExactValues(ModelPart& r_model_part);   // Get the L2 norm using the shape functions (no recovered field)
    double getL2NormFluidAccelWithRecoveryUsingGaussExactValues(ModelPart& r_model_part);  // Get the L2 norm using the exact values 

    virtual ~ErrorNormEthierFieldCalculator(){}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

private:
    const double mA;
    const double mD;
    const bool mNormalizeResult;
    EthierFlowField mFlowField;

}; // Class ErrorNormEthierFieldCalculator

} // namespace Kratos.

#endif // KRATOS_ERROR_ETHIER_FIELD defined