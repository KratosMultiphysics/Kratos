//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez

#if !defined(KRATOS_POTENTIAL_FLOW_UTILITIES_H_INCLUDED )
#define  KRATOS_POTENTIAL_FLOW_UTILITIES_H_INCLUDED


// Project includes
#include "containers/array_1d.h"
#include "includes/ublas_interface.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{
    // forward-declaring
    class Element;
    class ProcessInfo;
    class ModelPart;
    template< class TDataType >
    class GlobalPointersVector;
    template< class TDataType >
    class Dof;
    class Node;
    template< class TPointType >
    class Geometry;

namespace PotentialFlowUtilities
{
template <unsigned int TNumNodes, unsigned int TDim>
struct ElementalData
{
    template<typename TGeometryType>
    ElementalData(const TGeometryType& rGeometry)
    {
        GeometryUtils::CalculateGeometryData(rGeometry, DN_DX, N, vol);
    }

    array_1d<double, TNumNodes> potentials, distances;
    double vol;

    BoundedMatrix<double, TNumNodes, TDim> DN_DX;
    array_1d<double, TNumNodes> N;
};

typedef Node NodeType;
typedef Geometry<NodeType> GeometryType;

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) array_1d<double, NumNodes> GetWakeDistances(const Element& rElement);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) BoundedVector<double, NumNodes> GetPotentialOnNormalElement(const Element& rElement);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) BoundedVector<double, 2 * NumNodes> GetPotentialOnWakeElement(
    const Element& rElement, const array_1d<double, NumNodes>& rDistances);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) BoundedVector<double, NumNodes> GetPotentialOnUpperWakeElement(
    const Element& rElement, const array_1d<double, NumNodes>& rDistances);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) BoundedVector<double, NumNodes> GetPotentialOnLowerWakeElement(
    const Element& rElement, const array_1d<double, NumNodes>& rDistances);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) array_1d<double, Dim> ComputeVelocityNormalElement(const Element& rElement);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) array_1d<double, Dim> ComputeVelocityUpperWakeElement(const Element& rElement);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) array_1d<double, Dim> ComputeVelocityLowerWakeElement(const Element& rElement);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) array_1d<double, Dim> ComputeVelocity(const Element& rElement);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) array_1d<double, Dim> ComputePerturbedVelocity(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) array_1d<double, Dim> ComputePerturbedVelocityLowerElement(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeMaximumVelocitySquared(const ProcessInfo& rCurrentProcessInfo);

KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeVacuumVelocitySquared(const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeClampedVelocitySquared(const array_1d<double, Dim>& rVelocity, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeVelocityMagnitude(const double localMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeIncompressiblePressureCoefficient(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputePerturbationIncompressiblePressureCoefficient(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeCompressiblePressureCoefficient(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputePerturbationCompressiblePressureCoefficient(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeLocalSpeedOfSound(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeLocalSpeedofSoundSquared(const array_1d<double, Dim>& rVelocity,const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeSquaredSpeedofSoundFactor(const double localVelocitySquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputePerturbationLocalSpeedOfSound(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeLocalMachNumber(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeLocalMachNumberSquared(const array_1d<double, Dim>& rVelocity, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeDerivativeLocalMachSquaredWRTVelocitySquared(const array_1d<double, Dim>& rVelocity, const double localMachNumberSquared,const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputePerturbationLocalMachNumber(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeUpwindFactor(double localMachNumberSquared,const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double SelectMaxUpwindFactor(const array_1d<double, Dim>& rCurrentVelocity, const array_1d<double, Dim>& rUpwindVelocity, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) size_t ComputeUpwindFactorCase(array_1d<double, 3>& rUpwindFactorOptions);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeUpwindFactorDerivativeWRTMachSquared(const double localMachNumberSquared,const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeUpwindFactorDerivativeWRTVelocitySquared(const array_1d<double, Dim>& rVelocity,const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeDensity(const double localMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeUpwindedDensity(const array_1d<double, Dim>& rCurrentVelocity, const array_1d<double, Dim>& rUpwindVelocity, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeDensityDerivativeWRTVelocitySquared(const double localMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicAccelerating(const array_1d<double, Dim>& rCurrentVelocity, const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicDeaccelerating(const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicAccelerating(const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicDeaccelerating(const array_1d<double, Dim>& rUpwindVelocity, const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) bool CheckIfElementIsCutByDistance(const BoundedVector<double, NumNodes>& rNodalDistances);

KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) bool CheckIfElementIsTrailingEdge(const Element& rElement);

template <int Dim>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) void CheckIfWakeConditionsAreFulfilled(const ModelPart& rWakeModelPart, const double& rTolerance, const int& rEchoLevel);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) bool CheckWakeCondition(const Element& rElement, const double& rTolerance, const int& rEchoLevel);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) void GetSortedIds(std::vector<size_t>& Ids, const GeometryType& rGeom);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) void GetNodeNeighborElementCandidates(GlobalPointersVector<Element>& ElementCandidates, const GeometryType& rGeom);

template<int Dim>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) Vector ComputeKuttaNormal(const double angle);

template <class TContainerType>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) double CalculateArea(TContainerType& rContainer);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) void ComputePotentialJump(ModelPart& rWakeModelPart);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) void AddKuttaConditionPenaltyTerm(const Element& rElement,
                              Matrix& rLeftHandSideMatrix,
                              Vector& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) void AddKuttaConditionPenaltyPerturbationRHS(const Element& rElement,
                              Vector& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) void AddKuttaConditionPenaltyPerturbationLHS(const Element& rElement,
                              Matrix& rLeftHandSideMatrix,
                              const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) void AddPotentialGradientStabilizationTerm(Element& rElement,
                              Matrix& rLeftHandSideMatrix,
                              Vector& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo);

} // namespace PotentialFlow
} // namespace Kratos

#endif // KRATOS_POTENTIAL_FLOW_UTILITIES_H_INCLUDED  defined