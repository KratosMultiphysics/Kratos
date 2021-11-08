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
    template< std::size_t TDimension, class TDofType >
    class Node;
    template< class TPointType >
    class Geometry;

namespace PotentialFlowUtilities
{
template <unsigned int TNumNodes, unsigned int TDim>
struct ElementalData{
    array_1d<double, TNumNodes> potentials, distances;
    double vol;

    BoundedMatrix<double, TNumNodes, TDim> DN_DX;
    array_1d<double, TNumNodes> N;
};

typedef Node < 3, Dof<double> > NodeType;
typedef Geometry<NodeType> GeometryType;

template <int Dim, int NumNodes>
array_1d<double, NumNodes> GetWakeDistances(const Element& rElement);

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes> GetPotentialOnNormalElement(const Element& rElement);

template <int Dim, int NumNodes>
BoundedVector<double, 2 * NumNodes> GetPotentialOnWakeElement(
    const Element& rElement, const array_1d<double, NumNodes>& rDistances);

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes> GetPotentialOnUpperWakeElement(
    const Element& rElement, const array_1d<double, NumNodes>& rDistances);

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes> GetPotentialOnLowerWakeElement(
    const Element& rElement, const array_1d<double, NumNodes>& rDistances);

template <int Dim, int NumNodes>
array_1d<double, Dim> ComputeVelocityNormalElement(const Element& rElement);

template <int Dim, int NumNodes>
array_1d<double, Dim> ComputeVelocityUpperWakeElement(const Element& rElement);

template <int Dim, int NumNodes>
array_1d<double, Dim> ComputeVelocityLowerWakeElement(const Element& rElement);

template <int Dim, int NumNodes>
array_1d<double, Dim> ComputeVelocity(const Element& rElement);

template <int Dim, int NumNodes>
array_1d<double, Dim> ComputePerturbedVelocity(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
array_1d<double, Dim> ComputePerturbedVelocityLowerElement(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeMaximumVelocitySquared(const ProcessInfo& rCurrentProcessInfo);

double ComputeVacuumVelocitySquared(const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeClampedVelocitySquared(const array_1d<double, Dim>& rVelocity, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeVelocityMagnitude(const double localMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeIncompressiblePressureCoefficient(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputePerturbationIncompressiblePressureCoefficient(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeCompressiblePressureCoefficient(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputePerturbationCompressiblePressureCoefficient(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeLocalSpeedOfSound(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeLocalSpeedofSoundSquared(const array_1d<double, Dim>& rVelocity,const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeSquaredSpeedofSoundFactor(const double localVelocitySquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputePerturbationLocalSpeedOfSound(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeLocalMachNumber(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeLocalMachNumberSquared(const array_1d<double, Dim>& rVelocity, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeDerivativeLocalMachSquaredWRTVelocitySquared(const array_1d<double, Dim>& rVelocity, const double localMachNumberSquared,const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputePerturbationLocalMachNumber(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeUpwindFactor(double localMachNumberSquared,const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double SelectMaxUpwindFactor(const array_1d<double, Dim>& rCurrentVelocity, const array_1d<double, Dim>& rUpwindVelocity, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
size_t ComputeUpwindFactorCase(array_1d<double, 3>& rUpwindFactorOptions);

template <int Dim, int NumNodes>
double ComputeUpwindFactorDerivativeWRTMachSquared(const double localMachNumberSquared,const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeUpwindFactorDerivativeWRTVelocitySquared(const array_1d<double, Dim>& rVelocity,const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeDensity(const double localMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeUpwindedDensity(const array_1d<double, Dim>& rCurrentVelocity, const array_1d<double, Dim>& rUpwindVelocity, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeDensityDerivativeWRTVelocitySquared(const double localMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicAccelerating(const array_1d<double, Dim>& rCurrentVelocity, const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeUpwindedDensityDerivativeWRTVelocitySquaredSupersonicDeaccelerating(const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicAccelerating(const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeUpwindedDensityDerivativeWRTUpwindVelocitySquaredSupersonicDeaccelerating(const array_1d<double, Dim>& rUpwindVelocity, const double currentMachNumberSquared, const double upwindMachNumberSquared, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
bool CheckIfElementIsCutByDistance(const BoundedVector<double, NumNodes>& rNodalDistances);

bool CheckIfElementIsTrailingEdge(const Element& rElement);

template <int Dim>
void CheckIfWakeConditionsAreFulfilled(const ModelPart& rWakeModelPart, const double& rTolerance, const int& rEchoLevel);

template <int Dim, int NumNodes>
bool CheckWakeCondition(const Element& rElement, const double& rTolerance, const int& rEchoLevel);

template <int Dim, int NumNodes>
void GetSortedIds(std::vector<size_t>& Ids, const GeometryType& rGeom);

template <int Dim, int NumNodes>
void GetNodeNeighborElementCandidates(GlobalPointersVector<Element>& ElementCandidates, const GeometryType& rGeom);

template<int Dim>
Vector ComputeKuttaNormal(const double angle);

template <class TContainerType>
double CalculateArea(TContainerType& rContainer);

template <int Dim, int NumNodes>
void ComputePotentialJump(ModelPart& rWakeModelPart);

template <int Dim, int NumNodes>
void AddKuttaConditionPenaltyTerm(const Element& rElement,
                              Matrix& rLeftHandSideMatrix,
                              Vector& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
void AddPotentialGradientStabilizationTerm(Element& rElement,
                              Matrix& rLeftHandSideMatrix,
                              Vector& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo);

} // namespace PotentialFlow
} // namespace Kratos

#endif // KRATOS_POTENTIAL_FLOW_UTILITIES_H_INCLUDED  defined