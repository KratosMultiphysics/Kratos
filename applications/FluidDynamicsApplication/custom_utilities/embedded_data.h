//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#include <algorithm>
#include <cstddef>
#if !defined(KRATOS_EMBEDDED_DATA_H)
#define KRATOS_EMBEDDED_DATA_H

#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fluid_element_data.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template< class TFluidData >
class EmbeddedData : public TFluidData
{
public:

///@name Type Definitions
///@{

using NodalScalarData = typename TFluidData::NodalScalarData;
using NodalVectorData = typename TFluidData::NodalVectorData;

typedef GeometryData::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;
typedef std::vector<array_1d<double,3>> InterfaceNormalsType;

typedef Node NodeType;
typedef std::vector< NodeType::Pointer > NodePointerVectorType;

///@}
///@name Public Members
///@{

bool IsSlip;
bool ApplyNitscheBoundaryImposition;
bool ApplyConstraints;

double SlipLength;
double PenaltyCoefficient;

NodalScalarData Distance;

Matrix PositiveSideN;
ShapeFunctionsGradientsType PositiveSideDNDX;
Vector PositiveSideWeights;

Matrix PositiveInterfaceN;
ShapeFunctionsGradientsType PositiveInterfaceDNDX;
Vector PositiveInterfaceWeights;
InterfaceNormalsType PositiveInterfaceUnitNormals;

std::vector< size_t > PositiveIndices;
std::vector< size_t > NegativeIndices;

std::size_t NumPositiveNodes;
std::size_t NumNegativeNodes;

std::size_t LocalConstraintSize;
NodePointerVectorType ElementsNodesAndMasters;
Matrix ConstraintsRelationMatrix;

///@}
///@name Public Operations
///@{

void Initialize(
    const Element& rElement,
    const ProcessInfo& rProcessInfo) override
{
    TFluidData::Initialize(rElement, rProcessInfo);
    const Geometry<Node >& r_geometry = rElement.GetGeometry();
    this->FillFromHistoricalNodalData(Distance, DISTANCE, r_geometry);

    NumPositiveNodes = 0;
    NumNegativeNodes = 0;

    IsSlip = rElement.Is(SLIP) ? true : false;
}

/**
 * @brief Fills the boundary condition data fields
 * This method needs to be called in cut elements. It fills the data structure fields related to the boundary
 * condition imposition (slip length for the slip formulation and penalty coefficient for both formulations).
 * @param rProcessInfo
 */
void InitializeBoundaryConditionData(const ProcessInfo &rProcessInfo)
{
    if (IsSlip){
        this->FillFromProcessInfo(SlipLength, SLIP_LENGTH, rProcessInfo);
    }
    this->FillFromProcessInfo(PenaltyCoefficient, PENALTY_COEFFICIENT, rProcessInfo);
    this->FillFromProcessInfo(ApplyNitscheBoundaryImposition, APPLY_NITSCHE_BOUNDARY_IMPOSITION, rProcessInfo);
}

/** @brief
 * This method needs to be called in cut elements. It fills the data structure fields related to applying constraints
 * on negative nodes of cut elements.
 * @param rElement
 * @param rProcessInfo
 */
void InitializeConstraintData(const Element& rElement, const std::size_t BlockSize, const bool CalculateRelations)
{
    // Clear vector of pointers to nodes of the element plus master nodes of constrained negative nodes
    ElementsNodesAndMasters.clear();

    // Add nodes of the element to the nodes of the local system and check whether constraints need to be considered,
    // which is the case if one of the element's nodes has a negative distance value and is a slave node
    ApplyConstraints = false;
    for (auto &r_node : rElement.GetGeometry()) {
        ElementsNodesAndMasters.push_back(&r_node);

        if (r_node.GetValue(APPLY_EMBEDDED_CONSTRAINTS)) {
            ApplyConstraints = true;
        }
    }

    // If constraints need to be considered, add the corresponding master nodes to the local system and build relation matrix
    if (ApplyConstraints) {
        // Get additional master nodes, which are not nodes of the element already
        for (auto &r_node : rElement.GetGeometry()) {
            if (r_node.GetValue(APPLY_EMBEDDED_CONSTRAINTS)) {

                // Get vector of master node pointers for respective slave node
                const NodePointerVectorType NegNodeConstraintMasters = r_node.GetValue(EMBEDDED_CONSTRAINT_MASTERS);

                // Add master node if it was not already added because it is one of the element's nodes (unique entries)
                for (auto p_master : NegNodeConstraintMasters) {
                    auto it = std::find(ElementsNodesAndMasters.begin(), ElementsNodesAndMasters.end(), p_master);
                    if (it == ElementsNodesAndMasters.end()) {
                        ElementsNodesAndMasters.push_back(p_master);
                    }
                }
            }
        }
        KRATOS_WATCH("Adapting size of local system: " << ElementsNodesAndMasters.size());
        // Build the multi-point constraint's relation matrix if needed
        if (CalculateRelations) {
            //KRATOS_WATCH("Building constraint relation matrix ...");
            BuildRelationMatrix(rElement, BlockSize);
        }
    }

    // Calculate size of local system for resizing LHS and RHS, only differs from LocalSize
    // when element has constraints on its negative nodes, which master nodes are not nodes of the element
    LocalConstraintSize = ElementsNodesAndMasters.size() * BlockSize;
}

/** @brief
 * This method builds the constraint relation matrix one of the cut element's negative nodes is constrained.
 * @param rElement
 * @param BlockSize
 */
void BuildRelationMatrix(const Element& rElement, const std::size_t BlockSize)
{
    // Initialize relation matrix as identity matrix (so unconstrained dofs are not changed)
    const std::size_t size = ElementsNodesAndMasters.size() * BlockSize;
    if (ConstraintsRelationMatrix.size1() != size)
        ConstraintsRelationMatrix.resize(size, size, false);
    noalias(ConstraintsRelationMatrix) = IdentityMatrix(size);

    //Initialize the element's node counter for the position of the slave node in the local system
    std::size_t i_node = 0;

    // Fill relation matrix using constraints' master weights
    for (auto &r_node : rElement.GetGeometry()) {
        if (r_node.GetValue(APPLY_EMBEDDED_CONSTRAINTS)) {
            const std::size_t dim = BlockSize - 1;  //TODO: safe?

            // Get vector of master node pointers and master weights matrix for respective node
            const NodePointerVectorType NegNodeConstraintMasters = r_node.GetValue(EMBEDDED_CONSTRAINT_MASTERS);
            const Matrix NegNodeConstraintWeights = r_node.GetValue(EMBEDDED_CONSTRAINT_MASTER_WEIGHTS);

            //KRATOS_WATCH(NegNodeConstraintWeights);

            KRATOS_ERROR_IF_NOT(NegNodeConstraintWeights.size1() == dim)
                << "Size1 of master weights of an embedded constrained does not match the node's velocity dofs.";
            KRATOS_ERROR_IF_NOT(NegNodeConstraintWeights.size2() == NegNodeConstraintMasters.size()*dim)
                << "Number of master nodes of an embedded constrained node and its master weights do not match.";

            // Get position of first dof of slave node (VELOCITY_X)
            const std::size_t slave_x_position = i_node * BlockSize;

            // Find position of each master node and add weights for velocities
            std::size_t i_master_node = 0;
            for (auto p_master : NegNodeConstraintMasters) {
                // Find entry iterator of master node in the local system
                auto it = std::find(ElementsNodesAndMasters.begin(), ElementsNodesAndMasters.end(), p_master);
                // Get position of first dof of master node (VELOCITY_X)
                const std::size_t master_x_position = (it - ElementsNodesAndMasters.begin()) * BlockSize;
                // Add master weights of respective master dof (column) to slave dof (row)
                for (std::size_t d1 = 0; d1 < dim; ++d1) {
                    for (std::size_t d2 = 0; d2 < dim; ++d2) {
                        ConstraintsRelationMatrix(slave_x_position+d1, master_x_position+d2) = NegNodeConstraintWeights(d1, i_master_node*dim+d2);
                    }
                }
                i_master_node++;
            }

            // Set slave dofs diagonal values to zero (VELOCITY_X, VELOCITY_Y and if 3D VELOCITY_Z)
            for (std::size_t d = 0; d < dim; ++d) {
                ConstraintsRelationMatrix(slave_x_position+d, slave_x_position+d) = 0.0;
            }
        }
        i_node++;
    }
    //KRATOS_WATCH(ConstraintsRelationMatrix);
}

static int Check(const Element& rElement, const ProcessInfo& rProcessInfo)
{
    const Geometry< Node >& r_geometry = rElement.GetGeometry();
    for (unsigned int i = 0; i < TFluidData::NumNodes; i++) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE,r_geometry[i]);
    }

    int out = TFluidData::Check(rElement,rProcessInfo);
    return out;
}

bool IsCut() {
    return (NumPositiveNodes > 0) && (NumNegativeNodes > 0);
}

///@}

};

///@}

///@}

}

#endif