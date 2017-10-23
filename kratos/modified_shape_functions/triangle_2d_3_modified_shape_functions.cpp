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
#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

namespace Kratos
{
    
    /// Triangle2D3ModifiedShapeFunctions implementation
    /// Default constructor
    Triangle2D3ModifiedShapeFunctions::Triangle2D3ModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) :
        ModifiedShapeFunctions(pInputGeometry, rNodalDistances), 
        mpTriangleSplitter(new DivideTriangle2D3(*pInputGeometry, rNodalDistances)) {

        // Perform the element splitting
        mpTriangleSplitter->GenerateDivision();
        mpTriangleSplitter->GenerateIntersectionsSkin();
    };

    /// Destructor
    Triangle2D3ModifiedShapeFunctions::~Triangle2D3ModifiedShapeFunctions() {};

    /// Turn back information as a string.
    std::string Triangle2D3ModifiedShapeFunctions::Info() const {
        return "Triangle2D3N modified shape functions computation class.";
    };

    /// Print information about this object.
    void Triangle2D3ModifiedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
        rOStream << "Triangle2D3N modified shape functions computation class.";
    };

    /// Print object's data.
    void Triangle2D3ModifiedShapeFunctions::PrintData(std::ostream& rOStream) const {
        const GeometryPointerType p_geometry = this->GetInputGeometry();
        const Vector nodal_distances = this->GetNodalDistances();
        rOStream << "Triangle2D3N modified shape functions computation class:\n";
        rOStream << "\tGeometry type: " << (*p_geometry).Info() << "\n";
        std::stringstream distances_buffer;
        for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
            distances_buffer << std::to_string(nodal_distances(i)) << " ";
        }
        rOStream << "\tDistance values: " << distances_buffer.str();
    };

    // Internally computes the splitting pattern and returns all the shape function values for the positive side.
    void Triangle2D3ModifiedShapeFunctions::GetPositiveSideShapeFunctionsAndGradientsValues(Matrix &rPositiveSideShapeFunctionsValues,
                                                                                            std::vector<Matrix> &rPositiveSideShapeFunctionsGradientsValues,
                                                                                            Vector &rPositiveSideWeightsValues,
                                                                                            const IntegrationMethodType IntegrationMethod) {

        if (this->mIsSplit) {
            // Get the intersection points condensation matrix
            Matrix p_matrix;
            SetIntersectionPointsCondensationMatrix(p_matrix, 
                                                    mpTriangleSplitter->mEdgeNodeI, 
                                                    mpTriangleSplitter->mEdgeNodeJ, 
                                                    mpTriangleSplitter->mSplitEdges, 
                                                    mpTriangleSplitter->mSplitEdgesSize);

            // Compute the positive side values
            this->ComputeValuesOnOneSide(rPositiveSideShapeFunctionsValues, 
                                         rPositiveSideShapeFunctionsGradientsValues, 
                                         rPositiveSideWeightsValues,
                                         mpTriangleSplitter->mPositiveSubdivisions,
                                         p_matrix,
                                         IntegrationMethod);
        } else {
            KRATOS_ERROR << "Using the GetPositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
        }
    };

    // Internally computes the splitting pattern and returns all the shape function values for the negative side.
    void Triangle2D3ModifiedShapeFunctions::GetNegativeSideShapeFunctionsAndGradientsValues(Matrix &rNegativeSideShapeFunctionsValues,
                                                                                            std::vector<Matrix> &rNegativeSideShapeFunctionsGradientsValues,
                                                                                            Vector &rNegativeSideWeightsValues,
                                                                                            const IntegrationMethodType IntegrationMethod) {       

        if (this->mIsSplit) {
            // Get the intersection points condensation matrix
            Matrix p_matrix;
            SetIntersectionPointsCondensationMatrix(p_matrix, 
                                                    mpTriangleSplitter->mEdgeNodeI, 
                                                    mpTriangleSplitter->mEdgeNodeJ, 
                                                    mpTriangleSplitter->mSplitEdges, 
                                                    mpTriangleSplitter->mSplitEdgesSize);

            // Compute the negative side values
            this->ComputeValuesOnOneSide(rNegativeSideShapeFunctionsValues, 
                                         rNegativeSideShapeFunctionsGradientsValues, 
                                         rNegativeSideWeightsValues,
                                         mpTriangleSplitter->mNegativeSubdivisions,
                                         p_matrix,
                                         IntegrationMethod);
        } else {
            KRATOS_ERROR << "Using the GetNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
        }
    };

    // Internally computes the splitting pattern and returns all the shape function values for the positive interface side.
    void Triangle2D3ModifiedShapeFunctions::GetInterfacePositiveSideShapeFunctionsAndGradientsValues(Matrix &rInterfacePositiveSideShapeFunctionsValues,
                                                                                                     std::vector<Matrix> &rInterfacePositiveSideShapeFunctionsGradientsValues,
                                                                                                     Vector &rInterfacePositiveSideWeightsValues,
                                                                                                     const IntegrationMethodType IntegrationMethod) {

        // if (is_divided) {
        //     // Get the intersection points condensation matrix
        //     Matrix p_matrix;
        //     SetIntersectionPointsCondensationMatrix(p_matrix, 
        //                                             triangle_splitter.mEdgeNodeI, 
        //                                             triangle_splitter.mEdgeNodeJ, 
        //                                             triangle_splitter.mSplitEdges, 
        //                                             triangle_splitter.mSplitEdgesSize);

        //     // Compute the positive side values
        //     this->ComputeValuesOnOneSide(rPositiveSideShapeFunctionsValues, 
        //                                 rPositiveSideShapeFunctionsGradientsValues, 
        //                                 rPositiveSideWeightsValues,
        //                                 positive_subdivisions,
        //                                 p_matrix,
        //                                 IntegrationMethod);
        // } else {
        //     KRATOS_ERROR << "Using the GetPositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
        // }
    };

    // Internally computes the splitting pattern and returns all the shape function values for the negative interface side.
    void Triangle2D3ModifiedShapeFunctions::GetInterfaceNegativeSideShapeFunctionsAndGradientsValues(Matrix &rInterfaceNegativeSideShapeFunctionsValues,
                                                                                                     std::vector<Matrix> &rInterfaceNegativeSideShapeFunctionsGradientsValues,
                                                                                                     Vector &rInterfaceNegativeSideWeightsValues,
                                                                                                     const IntegrationMethodType IntegrationMethod) {

        // if (is_divided) {
        //     // Get the intersection points condensation matrix
        //     Matrix p_matrix;
        //     SetIntersectionPointsCondensationMatrix(p_matrix, 
        //                                             triangle_splitter.mEdgeNodeI, 
        //                                             triangle_splitter.mEdgeNodeJ, 
        //                                             triangle_splitter.mSplitEdges, 
        //                                             triangle_splitter.mSplitEdgesSize);

        //     // Compute the negative side values
        //     this->ComputeValuesOnOneSide(rNegativeSideShapeFunctionsValues, 
        //                                 rNegativeSideShapeFunctionsGradientsValues, 
        //                                 rNegativeSideWeightsValues,
        //                                 negative_subdivisions,
        //                                 p_matrix,
        //                                 IntegrationMethod);
        // } else {
        //     KRATOS_ERROR << "Using the GetNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
        // }
    };

};
