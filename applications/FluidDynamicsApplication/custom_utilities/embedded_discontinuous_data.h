//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_EMBEDDED_DISCONTINUOUS_DATA_H)
#define KRATOS_EMBEDDED_DISCONTINUOUS_DATA_H

#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fluid_element_data.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template< class TFluidData >
class EmbeddedDiscontinuousData : public TFluidData
{
public:

///@name Type Definitions
///@{

using NodalScalarData = typename TFluidData::NodalScalarData;
using NodalVectorData = typename TFluidData::NodalVectorData;

typedef GeometryData::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;
typedef std::vector<array_1d<double,3>> InterfaceNormalsType;

/// Number of edges of the element (simplex elements are assumed)
constexpr static std::size_t NumEdges = (TFluidData::NumNodes == 3) ? 3 : 6;

///@}
///@name Public Members
///@{

double SlipLength;
double PenaltyCoefficient;

NodalScalarData ElementalDistances;
Vector ElementalEdgeDistances;
Vector ElementalEdgeDistancesExtrapolated;

Matrix PositiveSideN;
Matrix NegativeSideN;
ShapeFunctionsGradientsType PositiveSideDNDX;
ShapeFunctionsGradientsType NegativeSideDNDX;
Vector PositiveSideWeights;
Vector NegativeSideWeights;

Matrix PositiveInterfaceN;
Matrix NegativeInterfaceN;
ShapeFunctionsGradientsType PositiveInterfaceDNDX;
ShapeFunctionsGradientsType NegativeInterfaceDNDX;
Vector PositiveInterfaceWeights;
Vector NegativeInterfaceWeights;
InterfaceNormalsType PositiveInterfaceUnitNormals;
InterfaceNormalsType NegativeInterfaceUnitNormals;

std::vector< size_t > PositiveIndices;
std::vector< size_t > NegativeIndices;

std::size_t NumPositiveNodes;
std::size_t NumNegativeNodes;
std::size_t NumIntersectedEdges;
std::size_t NumIntersectedEdgesExtrapolated;

///@}
///@name Public Operations
///@{

/**
 * @brief Discontinuous embedded formulation data container initialization
 * This method initializes the discontinuous embedded formulation data container. This implies to intialize
 * the base formulation data container as well as to get the elemental distances from the elemental variable
 * ELEMENTAL_DISTANCES (note that this requires the ELEMENTAL_DISTANCES to be set before this operation) and
 * the elemental edge distances from ELEMENTAL_EDGE_DISTANCES and ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED.
 * @param rElement reference to the element that owns the data container
 * @param rProcessInfo reference to the current ProcessInfo container
 */
void Initialize(
    const Element& rElement,
    const ProcessInfo& rProcessInfo) override
{
    TFluidData::Initialize(rElement, rProcessInfo);
    this->FillFromElementData(ElementalDistances, ELEMENTAL_DISTANCES, rElement);
    this->FillFromElementData(ElementalEdgeDistances, ELEMENTAL_EDGE_DISTANCES, rElement);
    this->FillFromElementData(ElementalEdgeDistancesExtrapolated, ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED, rElement);

    NumPositiveNodes = 0;
    NumNegativeNodes = 0;
    NumIntersectedEdges = 0;
    NumIntersectedEdgesExtrapolated = 0;
}

/**
 * @brief Fills the boundary condition data fields
 * This method needs to be called in cut elements. It fills the data structure fields related to the boundary
 * condition imposition (slip length and penalty coefficient) by retrieving their value from the ProcessInfo.
 * @param rProcessInfo
 */
void InitializeBoundaryConditionData(const ProcessInfo& rProcessInfo)
{
    this->FillFromProcessInfo(SlipLength, SLIP_LENGTH, rProcessInfo);
    this->FillFromProcessInfo(PenaltyCoefficient, PENALTY_COEFFICIENT, rProcessInfo);
}

/**
 * @brief Discontinous embedded formulation data container check
 * Simple discontinuous embedded formulation data container check. The base formulation data container is
 * checked as well. Returns 0 if the check process succeeds.
 * @param rElement reference to the element that owns the data container
 * @param rProcessInfo reference to the current ProcessInfo container
 * @return int returns 0 if the check process succeeds
 */
static int Check(
    const Element& rElement,
    const ProcessInfo& rProcessInfo)
{
    int out = TFluidData::Check(rElement,rProcessInfo);
    return out;
}

/**
 * @brief Checks if the current element is intersected
 * Checks if the current element is intersected by checking the number of positive and negative distance nodes.
 * @return true if the element is intersected
 * @return false if the element is not intersected
 */
bool IsCut()
{
    if (IsIncised()) {
        return false;
    } else {
        return (NumPositiveNodes > 0) && (NumNegativeNodes > 0);
    }
}

/**
 * @brief Checks if the current element is partially intersected (incised) and if extrapolated intersections where calculated.
 * Checks if the current element is partially intersected by checking the number of extrapolated intersected edges
 * This number will only be non-zero if user provided flag to calculate extrapolated edge distances.
 * The case in which three edges of a tetrahedra are cut and element is only incised is also considered.
 * @return true if the element is incised
 * @return false if the element is not incised
 */
inline bool IsIncised()
{
    return NumIntersectedEdgesExtrapolated > 0 ? true : false;
}

///@}

};

///@}

///@}

}

#endif