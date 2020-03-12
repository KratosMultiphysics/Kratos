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
#include "includes/dof.h"

namespace Kratos
{
    class Element; // forward-declaring to not having to include it here
    class ProcessInfo; // forward-declaring to not having to include it here
    class ModelPart; // forward-declaring to not having to include it here

namespace PotentialFlowUtilities
{
template <unsigned int TNumNodes, unsigned int TDim>
struct ElementalData{
    array_1d<double, TNumNodes> potentials, distances;
    double vol;

    BoundedMatrix<double, TNumNodes, TDim> DN_DX;
    array_1d<double, TNumNodes> N;
};

typedef std::vector<std::size_t> EquationIdVectorType;
typedef std::vector< Dof<double>::Pointer > DofsVectorType;

template <int Dim, int NumNodes>
array_1d<double, NumNodes> GetWakeDistances(const Element& rElement);

template <int Dim, int NumNodes>
void GetEquationIdVectorNormalElement(const Element& rElement, EquationIdVectorType& rElementalIdList);

template <int Dim, int NumNodes>
void GetEquationIdVectorKuttaElement(const Element& rElement, EquationIdVectorType& rElementalIdList);

template <int Dim, int NumNodes>
void GetDofListNormalElement(const Element& rElement, DofsVectorType& rElementalDofList);

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
double ComputePerturbationLocalSpeedOfSound(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeMaximumVelocitySquared(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeLocalMachNumber(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputePerturbationLocalMachNumber(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputePerturbationDensity(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeDensityDerivative(const double& rDensity, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeUpwindDensity(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeSwitchingOperator(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeUpwindFactor(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeDerivativeUpwindFactorWRTMachNumberSquared(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
double ComputeDerivativeMachNumberSquaredWRTVelocitySquared(const Element& rElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes> ComputeDrhoDphiSupersonicAccelerating(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes> ComputeDrhoDphiUpSupersonicAccelerating(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes> ComputeDrhoDphiSupersonicDecelerating(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes> ComputeDrhoDphiUpSupersonicDecelerating(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes + 1> ComputeAndAssembleDrhoDphi(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes + 1> ComputeAndAssembleDrhoDphiSupersonicDecelerating(const Element& rElement, const Element& rUpstreamElement, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
BoundedVector<double, NumNodes + 1> AssembleDrhoDphi(const Element& rElement, const BoundedVector<double, NumNodes>& rDrhoDPhiCurrent, const Element& rUpstreamElement,  const BoundedVector<double, NumNodes>& rDrhoDPhiUpstream, const ProcessInfo& rCurrentProcessInfo);

template <int Dim, int NumNodes>
bool CheckIfElementIsCutByDistance(const BoundedVector<double, NumNodes>& rNodalDistances);

bool CheckIfElementIsTrailingEdge(const Element& rElement);

template <int Dim>
void CheckIfWakeConditionsAreFulfilled(const ModelPart& rWakeModelPart, const double& rTolerance, const int& rEchoLevel);

template <int Dim, int NumNodes>
bool CheckWakeCondition(const Element& rElement, const double& rTolerance, const int& rEchoLevel);

} // namespace PotentialFlow
} // namespace Kratos

#endif // KRATOS_POTENTIAL_FLOW_UTILITIES_H_INCLUDED  defined