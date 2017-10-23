//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "modified_shape_functions/modified_shape_functions.h"

namespace Kratos
{
    
    /// ModifiedShapeFunctions implementation
    /// Default constructor
    ModifiedShapeFunctions::ModifiedShapeFunctions(GeometryPointerType pInputGeometry, Vector& rNodalDistances) :
        mpInputGeometry(pInputGeometry),
        mrNodalDistances(rNodalDistances) {
    };

    /// Destructor
    ModifiedShapeFunctions::~ModifiedShapeFunctions() {};

    /// Turn back information as a string.
    std::string ModifiedShapeFunctions::Info() const {
        return "Modified shape functions computation base class.";
    };
    
    /// Print information about this object.
    void ModifiedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
        rOStream << "Modified shape functions computation base class.";
    };
    
    /// Print object's data.
    void ModifiedShapeFunctions::PrintData(std::ostream& rOStream) const {
        const GeometryPointerType p_geometry = this->GetInputGeometry();
        const Vector nodal_distances = this->GetNodalDistances();
        rOStream << "Modified shape functions computation base class:\n";
        rOStream << "\tGeometry type: " << (*p_geometry).Info() << "\n";
        std::stringstream distances_buffer;
        for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
            distances_buffer << std::to_string(nodal_distances(i)) << " ";
        }
        rOStream << "\tDistance values: " << distances_buffer.str();
    };
    
    // Returns the input original geometry.
    ModifiedShapeFunctions::GeometryPointerType ModifiedShapeFunctions::GetInputGeometry() const {
        return mpInputGeometry;
    };
    
    // Returns the nodal distances vector.
    Vector ModifiedShapeFunctions::GetNodalDistances() const {
        return mrNodalDistances;
    };

    // Internally computes the splitting pattern and returns all the shape function values for both sides.
    void ModifiedShapeFunctions::GetShapeFunctionsAndGradientsValues(Matrix& rPositiveSideShapeFunctionsValues,
                                                                     Matrix& rNegativeSideShapeFunctionsValues,
                                                                     std::vector<Matrix>& rPositiveSideShapeFunctionsGradientsValues,
                                                                     std::vector<Matrix>& rNegativeSideShapeFunctionsGradientsValues,
                                                                     Vector& rPositiveSideWeightsValues,
                                                                     Vector& rNegativeWeightsValues,
                                                                     const IntegrationMethodType IntegrationMethod) {
        KRATOS_ERROR << "Calling the base class GetShapeFunctionsAndGradientsValues method. Call the specific geometry one.";
    };

    // Internally computes the splitting pattern and returns all the shape function values for the positive side.
    void ModifiedShapeFunctions::GetPositiveSideShapeFunctionsAndGradientsValues(Matrix &rPositiveSideShapeFunctionsValues,
                                                                                 std::vector<Matrix> &rPositiveSideShapeFunctionsGradientsValues,
                                                                                 Vector &rPositiveSideWeightsValues,
                                                                                 const IntegrationMethodType IntegrationMethod)
    {
        KRATOS_ERROR << "Calling the base class GetShapeFunctionsAndGradientsValuesPositiveSide method. Call the specific geometry one.";
    };

    // Internally computes the splitting pattern and returns all the shape function values for the negative side.
    void ModifiedShapeFunctions::GetNegativeSideShapeFunctionsAndGradientsValues(Matrix &rNegativeSideShapeFunctionsValues,
                                                                                 std::vector<Matrix> &rNegativeSideShapeFunctionsGradientsValues,
                                                                                 Vector &rNegativeSideWeightsValues,
                                                                                 const IntegrationMethodType IntegrationMethod)
    {
        KRATOS_ERROR << "Calling the base class GetShapeFunctionsAndGradientsValuesPositiveSide method. Call the specific geometry one.";
    };
    
    // Sets the condensation matrix to transform the subdivsion values to entire element ones.
    void ModifiedShapeFunctions::SetIntersectionPointsCondensationMatrix(Matrix& rIntPointCondMatrix,
        const int EdgeNodeI[],
        const int EdgeNodeJ[],
        const int SplitEdges[],
        const unsigned int SplitEdgesSize) {

        // Initialize intersection points condensation matrix
        const unsigned int nedges = mpInputGeometry->EdgesNumber();
        const unsigned int nnodes = mpInputGeometry->PointsNumber();

        rIntPointCondMatrix = ZeroMatrix(SplitEdgesSize, nnodes);

        // Fill the original geometry points main diagonal
        for (unsigned int i = 0; i < nnodes; ++i) {
            rIntPointCondMatrix(i,i) = 1.0;
        }

        // Compute the intersection points contributions
        unsigned int row = nnodes;
        for (unsigned int idedge = 0; idedge < nedges; ++idedge) {
            // Check if the edge has an intersection point
            if (SplitEdges[nnodes+idedge] != -1) {
                // Get the nodes that compose the edge
                const unsigned int edge_node_i = EdgeNodeI[idedge];
                const unsigned int edge_node_j = EdgeNodeJ[idedge];

                // Compute the relative coordinate of the intersection point over the edge
                const double aux_node_rel_location = std::abs (mrNodalDistances(edge_node_i)/(mrNodalDistances(edge_node_j)-mrNodalDistances(edge_node_i)));

                // Store the relative coordinate values as the original geometry nodes sh. function value in the intersections
                rIntPointCondMatrix(row, edge_node_i) = aux_node_rel_location;
                rIntPointCondMatrix(row, edge_node_j) = 1.0 - aux_node_rel_location;
            }
            row++;
        }
    }

    // Given the subdivision pattern of either the positive or negative side, computes the shape function values. 
    void ModifiedShapeFunctions::ComputeValuesOnOneSide(Matrix &rShapeFunctionsValues,
                                                        std::vector<Matrix> &rShapeFunctionsGradientsValues,
                                                        Vector &rWeightsValues,
                                                        const std::vector<IndexedPointGeometryPointerType> &rSubdivisionsVector,
                                                        const Matrix &p_matrix,
                                                        const IntegrationMethodType IntegrationMethod)
    {
        KRATOS_ERROR << "Calling the base class ComputeValuesOnOneSide method. Call the specific geometry one.";
    };

}; // namespace Kratos
