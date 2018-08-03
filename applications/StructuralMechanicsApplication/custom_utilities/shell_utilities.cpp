//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes

// External includes

// Project includes
#include "shell_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
namespace ShellUtilities
{
using SizeType = std::size_t;
using IndexType = std::size_t;

double dN_seren_dxi(const int nNode, const double Xi, const double Eta)
{
    // Natural derivatives of 8-node serendipity shape functions
    switch (nNode)
    {
    case 1:
        return -(-Eta + 1.0)*(-0.25*Xi + 0.25) -
            0.25*(-Eta + 1.0)*(-Eta - Xi - 1.0);
    case 2:
        return (-Eta + 1.0)*(0.25*Xi + 0.25) +
            0.25*(-Eta + 1.0)*(-Eta + Xi - 1.0);
    case 3:
        return (Eta + 1.0)*(0.25*Xi + 0.25) +
            0.25*(Eta + 1.0)*(Eta + Xi - 1.0);
    case 4:
        return -(Eta + 1.0)*(-0.25*Xi + 0.25) -
            0.25*(Eta + 1.0)*(Eta - Xi - 1.0);
    case 5:
        return -1.0*Xi*(-Eta + 1.0);
    case 6:
        return -0.5*Eta*Eta + 0.5;
    case 7:
        return -1.0*Xi*(Eta + 1.0);
    case 8:
        return 0.5*Eta*Eta - 0.5;
    default:
        KRATOS_ERROR <<
            "Error: ELEMENT ShellThinElement3D4N, METHOD dN_seren_dxi"
            << std::endl;
    }
}

double dN_seren_deta(const int nNode, const double Xi, const double Eta)
{
    // Natural derivatives of 8-node serendipity shape functions
    switch (nNode)
    {
    case 1:
        return -(-Eta + 1.0)*(-0.25*Xi + 0.25) -
            (-0.25*Xi + 0.25)*(-Eta - Xi - 1.0);
    case 2:
        return -(-Eta + 1.0)*(0.25*Xi + 0.25) -
            (0.25*Xi + 0.25)*(-Eta + Xi - 1.0);
    case 3:
        return (Eta + 1.0)*(0.25*Xi + 0.25) +
            (0.25*Xi + 0.25)*(Eta + Xi - 1.0);
    case 4:
        return (Eta + 1.0)*(-0.25*Xi + 0.25) +
            (-0.25*Xi + 0.25)*(Eta - Xi - 1.0);
    case 5:
        return 0.5*Xi*Xi - 0.5;
    case 6:
        return -1.0*Eta*(Xi + 1.0);
    case 7:
        return -0.5*Xi*Xi + 0.5;
    case 8:
        return -1.0*Eta*(-Xi + 1.0);
    default:
        KRATOS_ERROR <<
            "Error: ELEMENT ShellThinElement3D4N, METHOD dN_seren_dxi"
            << std::endl;
    }
}

void InterpToStandardGaussPoints(double& rV1, double& rV2, double& rV3)
{
    const double vg1 = rV1;
    const double vg2 = rV2;
    const double vg3 = rV3;
#ifdef OPT_AVERAGE_RESULTS
    rV1 = (vg1 + vg2 + vg3) / 3.0;
    rV2 = (vg1 + vg2 + vg3) / 3.0;
    rV3 = (vg1 + vg2 + vg3) / 3.0;
#else
    rV1 = (2.0*vg1) / 3.0 - vg2 / 3.0 + (2.0*vg3) / 3.0;
    rV2 = (2.0*vg1) / 3.0 + (2.0*vg2) / 3.0 - vg3 / 3.0;
    rV3 = (2.0*vg2) / 3.0 - vg1 / 3.0 + (2.0*vg3) / 3.0;
#endif // OPT_AVERAGE_RESULTS
}

void InterpToStandardGaussPoints(std::vector< double >& rV)
{
    if (rV.size() != 3) return;
    InterpToStandardGaussPoints(rV[0], rV[1], rV[2]);
}

void InterpToStandardGaussPoints(std::vector< array_1d<double, 3> >& rV)
{
    if (rV.size() != 3) return;
    for (IndexType i = 0; i < 3; ++i)
        InterpToStandardGaussPoints(rV[0][i], rV[1][i], rV[2][i]);
}

void InterpToStandardGaussPoints(std::vector< array_1d<double, 6> >& rV)
{
    if (rV.size() != 3) return;
    for (IndexType i = 0; i < 6; ++i)
        InterpToStandardGaussPoints(rV[0][i], rV[1][i], rV[2][i]);
}

void InterpToStandardGaussPoints(std::vector< Vector >& rV)
{
    if (rV.size() != 3) return;
    SizeType ncomp = rV[0].size();
    for (int i = 1; i < 3; ++i)
        if (rV[i].size() != ncomp)
            return;
    for (IndexType i = 0; i < ncomp; ++i)
        InterpToStandardGaussPoints(rV[0][i], rV[1][i], rV[2][i]);
}

void InterpToStandardGaussPoints(std::vector< Matrix >& rV)
{
    if (rV.size() != 3) return;
    SizeType nrows = rV[0].size1();
    SizeType ncols = rV[0].size2();
    for (int i = 1; i < 3; ++i)
        if (rV[i].size1() != nrows || rV[i].size2() != ncols)
            return;
    for (IndexType i = 0; i < nrows; ++i)
        for (IndexType j = 0; j < ncols; ++j)
            InterpToStandardGaussPoints
            (rV[0](i, j), rV[1](i, j), rV[2](i, j));
}

bool IsOrthotropic(const Properties& rProps)
{
    return rProps.Has(SHELL_ORTHOTROPIC_LAYERS);
}

double GetThickness(const Properties& rProps)
{
    if (IsOrthotropic(rProps))
    {
        double thickness = 0.0;
        const auto& orthotropic_layers = rProps.GetValue(SHELL_ORTHOTROPIC_LAYERS);
        for (IndexType i=0; i<orthotropic_layers.size1(); ++i)
            thickness += orthotropic_layers(i,0);
        return thickness;
    }
    else
    {
        KRATOS_DEBUG_ERROR_IF_NOT(rProps.Has(THICKNESS));
        return rProps.GetValue(THICKNESS);
    }
}

double GetThickness(const Properties& rProps, const IndexType Index)
{
    if (IsOrthotropic(rProps))
    {
        const auto& orthotropic_layers = rProps.GetValue(SHELL_ORTHOTROPIC_LAYERS);
        return orthotropic_layers(Index,0);
    }
    else
    {
        KRATOS_DEBUG_ERROR_IF_NOT(rProps.Has(THICKNESS));
        return rProps.GetValue(THICKNESS);
    }
}

double GetDensity(const Properties& rProps, const IndexType Index)
{
    if (IsOrthotropic(rProps))
    {
        const auto& orthotropic_layers = rProps.GetValue(SHELL_ORTHOTROPIC_LAYERS);
        return orthotropic_layers(Index,2);
    }
    else
    {
        KRATOS_DEBUG_ERROR_IF_NOT(rProps.Has(DENSITY));
        return rProps.GetValue(DENSITY);
    }
}

double GetOrientationAngle(const Properties& rProps, const IndexType Index)
{
    if (IsOrthotropic(rProps))
    {
        const auto& orthotropic_layers = rProps.GetValue(SHELL_ORTHOTROPIC_LAYERS);
        double orientation_angle = orthotropic_layers(Index,1);

        orientation_angle = std::fmod(orientation_angle, 360.0);
        if(orientation_angle < 0.0)
            orientation_angle += 360.0;
        return orientation_angle;
    }
    else
        return 0.0;
}


double GetOffset(const Properties& rProps)
{
if (rProps.Has(SHELL_OFFSET))
    return rProps.GetValue(SHELL_OFFSET);
else
    return 0.0;
}

} // namespace ShellUtilities

} // namespace Kratos.


