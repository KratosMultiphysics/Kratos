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

#if !defined(KRATOS_ERROR_TORUS)
#define KRATOS_ERROR_TORUS

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

class ErrorNormTorusCalculator
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

    KRATOS_CLASS_POINTER_DEFINITION(ErrorNormTorusCalculator);

    ErrorNormTorusCalculator(const double major_radius, const double minor_radius, const double center_velocity): mMajorRadius(major_radius), mMinorRadius(minor_radius), mCenterVelocity(center_velocity)
    {}

    double getDistanceToCenter(const array_1d<double, 3>& coor);
    double getVelocityModule(const array_1d<double, 3>& coor);
    array_1d<double, 3> getVelocity(const array_1d<double, 3>& coor);
    void CalculateMaterialAcceleration(const array_1d<double, 3>& coor, array_1d<double, 3>& accel);

    double getL2NormFluidAccelWithoutRecoveryUsingGaussInterpolatedValues(ModelPart& r_model_part);   // Get the L2 norm using the shape functions (no recovered field)
    double getL2NormFluidAccelWithRecoveryUsingGaussInterpolatedValues(ModelPart& r_model_part);  // Get the L2 norm using the recovered field at the nodes, and then interpolated at the gauss points
    double getL2NormFluidAccelWithoutRecoveryUsingGaussExactValues(ModelPart& r_model_part);   // Get the L2 norm using the shape functions (no recovered field)
    double getL2NormFluidAccelWithRecoveryUsingGaussExactValues(ModelPart& r_model_part);  // Get the L2 norm using the exact values 

    virtual ~ErrorNormTorusCalculator(){}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

private:
    const double mMajorRadius;
    const double mMinorRadius;
    const double mCenterVelocity;

}; // Class ErrorNormTorusCalculator

} // namespace Kratos.

#endif // KRATOS_ERROR_TORUS defined