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
#include "modified_shape_functions/triangle_2d_3_ausas_incised_shape_functions.h"

namespace Kratos
{

/// Triangle2D3AusasIncisedShapeFunctions implementation
/// Default constructor
Triangle2D3AusasIncisedShapeFunctions::Triangle2D3AusasIncisedShapeFunctions(const GeometryPointerType pInputGeometry,
        const Vector& rExtrapolatedNodalDistances, const Vector& rExtrapolatedEdgeRatios) :
    Triangle2D3AusasModifiedShapeFunctions(pInputGeometry, rExtrapolatedNodalDistances), mExtraEdgeRatios(rExtrapolatedEdgeRatios) {};

/// Destructor
Triangle2D3AusasIncisedShapeFunctions::~Triangle2D3AusasIncisedShapeFunctions() {};

/// Turn back information as a string.
std::string Triangle2D3AusasIncisedShapeFunctions::Info() const {
    return "Triangle2D3N Ausas incised shape functions computation class.";
};

/// Print information about this object.
void Triangle2D3AusasIncisedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Triangle2D3N Ausas incised shape functions computation class.";
};

/// Print object's data.
void Triangle2D3AusasIncisedShapeFunctions::PrintData(std::ostream& rOStream) const {
    const GeometryPointerType p_geometry = this->GetInputGeometry();
    const Vector nodal_distances = this->GetNodalDistances();
    rOStream << "Triangle2D3N Ausas incised shape functions computation class:\n";
    rOStream << "\tGeometry type: " << (*p_geometry).Info() << "\n";
    std::stringstream distances_buffer;
    std::ostringstream stm;
    for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
        stm << nodal_distances(i);
        distances_buffer << stm.str() << " ";
    }
    rOStream << "\tExtrapolated distance values: " << distances_buffer.str();
};

// Returns the nodal distances vector.
const Vector& Triangle2D3AusasIncisedShapeFunctions::GetExtrapolatedEdgeRatios() const {
    return mExtraEdgeRatios;
};

}; //namespace Kratos
