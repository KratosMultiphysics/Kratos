//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "utilities/geometry_utilities.h"
#include "shallow_water_application_variables.h"
#include "wave_element.h"

namespace Kratos
{

template<std::size_t TNumNodes>
int WaveElement<TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int err = Element::Check(rCurrentProcessInfo);
    if (err != 0) return err;

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (const auto& r_node : this->GetGeometry())
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MANNING, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TOPOGRAPHY, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION, r_node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VERTICAL_VELOCITY, r_node)

        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(HEIGHT, r_node)
    }

    return err;

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rResult.size() != mLocalSize)
        rResult.resize(mLocalSize, false); // False says not to preserve existing storage!!

    const GeometryType& r_geom = GetGeometry();
    int counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rResult[counter++] = r_geom[i].GetDof(VELOCITY_X).EquationId();
        rResult[counter++] = r_geom[i].GetDof(VELOCITY_Y).EquationId();
        rResult[counter++] = r_geom[i].GetDof(HEIGHT).EquationId();
    }

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if(rElementalDofList.size() != mLocalSize)
        rElementalDofList.resize(mLocalSize);

    const GeometryType& r_geom = GetGeometry();
    int counter=0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rElementalDofList[counter++] = r_geom[i].pGetDof(VELOCITY_X);
        rElementalDofList[counter++] = r_geom[i].pGetDof(VELOCITY_Y);
        rElementalDofList[counter++] = r_geom[i].pGetDof(HEIGHT);
    }

    KRATOS_CATCH("")
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::GetValuesVector(Vector& rValues, int Step) const
{
    if (rValues.size() != mLocalSize)
        rValues.resize(mLocalSize, false);

    const GeometryType& r_geom = this->GetGeometry();
    IndexType counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(MOMENTUM_X, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(MOMENTUM_Y, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(HEIGHT, Step);
    }
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    if (rValues.size() != mLocalSize)
        rValues.resize(mLocalSize, false);

    const GeometryType& r_geom = this->GetGeometry();
    IndexType counter = 0;
    for (IndexType i = 0; i < TNumNodes; i++)
    {
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(ACCELERATION_X, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(ACCELERATION_Y, Step);
        rValues[counter++] = r_geom[i].FastGetSolutionStepValue(VERTICAL_VELOCITY, Step);
    }
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    KRATOS_ERROR << "ShallowWater2D3: This method is not supported by the formulation" << std::endl;
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    for (IndexType PointNumber = 0; PointNumber < 1; PointNumber++)
        rValues[PointNumber] = double(this->GetValue(rVariable));
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::InitializeData(ElementData& rData, const ProcessInfo& rCurrentProcessInfo)
{
    rData.stab_factor = rCurrentProcessInfo[STABILIZATION_FACTOR];
    rData.shock_stab_factor = rCurrentProcessInfo[SHOCK_STABILIZATION_FACTOR];
    rData.rel_dry_height = rCurrentProcessInfo[RELATIVE_DRY_HEIGHT];
    rData.gravity = rCurrentProcessInfo[GRAVITY_Z];
}

template<std::size_t TNumNodes>
void WaveElement<TNumNodes>::GetNodalData(ElementData& rData, const GeometryType& rGeometry)
{
    rData.length = rGeometry.Length();
    rData.height = 0.0;
    rData.velocity = ZeroVector(3);

    for (IndexType i = 0; i < 3; i++)
    {
        const IndexType block = 3 * i;

        const auto h = rGeometry[i].FastGetSolutionStepValue(HEIGHT);
        const array_1d<double,3> v = rGeometry[i].FastGetSolutionStepValue(VELOCITY);
        rData.height += h;
        rData.velocity += v;

        rData.topography[i] = rGeometry[i].FastGetSolutionStepValue(TOPOGRAPHY);

        rData.unknown[block]     = v[0];
        rData.unknown[block + 1] = v[1];
        rData.unknown[block + 2] = h;
    }
    const double lumping_factor = 1.0 / TNumNodes;
    rData.height = std::max(0.0, rData.height);
    rData.height *= lumping_factor;
    rData.velocity *= lumping_factor;
}

template<std::size_t TNumNodes>
double WaveElement<TNumNodes>::StabilizationParameter(const ElementData& rData) const
{
    return 0.0;
}

template<std::size_t TNumNodes>
typename WaveElement<TNumNodes>::LocalVectorType WaveElement<TNumNodes>::ToSystemVector(const array_1d<double,3>& rVector) const
{
    return ZeroVector(3*TNumNodes);
}

template<std::size_t TNumNodes>
typename WaveElement<TNumNodes>::LocalVectorType WaveElement<TNumNodes>::ToSystemVector(const double& rScalar) const
{
    return ZeroVector(3*TNumNodes);
}

template class WaveElement<3>;
template class WaveElement<4>;

} // namespace Kratos
