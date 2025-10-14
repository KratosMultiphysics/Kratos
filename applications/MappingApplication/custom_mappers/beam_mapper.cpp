//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan Ignacio Camarotti, Philipp Bucher
//
// See PhD Thesis Tianyang Wang Chapter 5

// System includes

// External includes
#include "utilities/math_utils.h"

// Project includes
#include "beam_mapper.h"
#include "custom_utilities/mapping_matrix_utilities.h"
#include "mapping_application_variables.h"
#include <math.h>
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "custom_searching/interface_communicator_mpi.h"
#endif

namespace Kratos {

void BeamMapperInterfaceInfo::ProcessSearchResult(const InterfaceObject& rInterfaceObject)
{
    SaveSearchResult(rInterfaceObject,false);
}

void BeamMapperInterfaceInfo::ProcessSearchResultForApproximation(const InterfaceObject& rInterfaceObject)
{
    const auto p_geom = rInterfaceObject.pGetBaseGeometry();
    double proj_dist;
    double proj_dist_nodes;

    const Point point_to_proj(this->Coordinates());
    Point projection_point;
    Point local_coords;

    mCoordinates = point_to_proj;

    Vector linear_shape_function_values;
    Vector hermitian_shape_function_values;
    Vector hermitian_shape_function_derivatives_values;
    std::vector<int> eq_ids;

    std::abs(GeometricalProjectionUtilities::FastProjectOnLine(*(p_geom), point_to_proj, projection_point));
    
    Point min_distance_to_nodes;
    double distance_to_node1, distance_to_node2;
    distance_to_node1 = norm_2(projection_point - (*p_geom)[0]);
    distance_to_node2 = norm_2(projection_point - (*p_geom)[1]);
    
    if(distance_to_node1 < distance_to_node2)
        proj_dist_nodes = distance_to_node1;
    else
        proj_dist_nodes = distance_to_node2;
    
    if (proj_dist_nodes < mClosestProjectionDistance) {
        SetIsApproximation();

        mPairingIndex = ProjectionUtilities::PairingIndex::Closest_Point;
        p_geom->IsInside(projection_point, local_coords, 1e-14);
        p_geom->PointLocalCoordinates(local_coords, projection_point);
        if(local_coords[0] < -1.0 )
            local_coords[0] = -1.0;
        if(local_coords[0] > 1.0 )
            local_coords[0] = 1.0;

        bool ComputeApproximation = 0;
        ProjectionUtilities::ProjectOnLine(*p_geom, point_to_proj, mLocalCoordTol, linear_shape_function_values, eq_ids, proj_dist, ComputeApproximation); // Kust to get eq_ids
        //ProjectionUtilities::ProjectOnLineHermitian(*p_geom, point_to_proj, mLocalCoordTol, hermitian_shape_function_values, hermitian_shape_function_derivatives_values, proj_dist, projection_point);

        p_geom->ShapeFunctionsValues(linear_shape_function_values, local_coords);
        ProjectionUtilities::HermitianShapeFunctionsValues(hermitian_shape_function_values, hermitian_shape_function_derivatives_values, local_coords);
        mClosestProjectionDistance = proj_dist_nodes;
        mNodeIds = eq_ids;

        const std::size_t num_values_linear = linear_shape_function_values.size();
        const std::size_t num_values_hermitian = hermitian_shape_function_values.size();
        const std::size_t num_values_hermitian_der = hermitian_shape_function_derivatives_values.size();
        
        if (mLinearShapeFunctionValues.size() != num_values_linear) mLinearShapeFunctionValues.resize(num_values_linear);
        for (std::size_t i=0; i<num_values_linear; ++i) {
            mLinearShapeFunctionValues[i] = linear_shape_function_values[i];
        }

        if (mHermitianShapeFunctionValues.size() != num_values_linear) mHermitianShapeFunctionValues.resize(num_values_hermitian);
        for (std::size_t i=0; i<num_values_hermitian; ++i) {
            mHermitianShapeFunctionValues[i] = hermitian_shape_function_values[i];
        }

        if (mHermitianShapeFunctionValuesDerivatives.size() != num_values_hermitian_der) mHermitianShapeFunctionValuesDerivatives.resize(num_values_hermitian_der);
        for (std::size_t i=0; i<num_values_hermitian_der; ++i) {
            mHermitianShapeFunctionValuesDerivatives[i] = hermitian_shape_function_derivatives_values[i];
        }
        
        mProjectionOfPoint = projection_point;
        
        mpInterfaceObject = make_shared<InterfaceGeometryObject>(rInterfaceObject.pGetBaseGeometry());
    }
}

void BeamMapperInterfaceInfo::SaveSearchResult(const InterfaceObject& rInterfaceObject,
                                               const bool ComputeApproximation)
{
    const auto p_geom = rInterfaceObject.pGetBaseGeometry();

    double proj_dist;

    const Point point_to_proj(this->Coordinates());
    Point projection_point;

    mCoordinates = point_to_proj;

    Vector linear_shape_function_values;
    Vector hermitian_shape_function_values;
    Vector hermitian_shape_function_derivatives_values;
    std::vector<int> eq_ids;

    ProjectionUtilities::PairingIndex pairing_index;

    const auto geom_family = p_geom->GetGeometryFamily();
    KRATOS_ERROR_IF(geom_family != GeometryData::KratosGeometryFamily::Kratos_Linear) << "Invalid geometry of the Origin! The geometry should be a beam!";

    // Calculating and storing the shape function values
    // Linear shape functions
    pairing_index = ProjectionUtilities::ProjectOnLine(*p_geom, point_to_proj, mLocalCoordTol, linear_shape_function_values, eq_ids, proj_dist, ComputeApproximation);
    // Hermitian shape functions
    ProjectionUtilities::ProjectOnLineHermitian(*p_geom, point_to_proj, mLocalCoordTol, hermitian_shape_function_values, hermitian_shape_function_derivatives_values, proj_dist, projection_point);
    const bool is_full_projection = (pairing_index == ProjectionUtilities::PairingIndex::Line_Inside);
    

    if (is_full_projection) {
        SetLocalSearchWasSuccessful();
    } else {
        if (!ComputeApproximation) {
            return;
        } else {
            SetIsApproximation();
        }
    }

    const std::size_t num_values_linear = linear_shape_function_values.size();
    const std::size_t num_values_hermitian = hermitian_shape_function_values.size();
    const std::size_t num_values_hermitian_der = hermitian_shape_function_derivatives_values.size();


    if (pairing_index > mPairingIndex || (pairing_index == mPairingIndex && proj_dist < mClosestProjectionDistance)) {
        mPairingIndex = pairing_index;
        mClosestProjectionDistance = proj_dist;
        mNodeIds = eq_ids;

        if (mLinearShapeFunctionValues.size() != num_values_linear) mLinearShapeFunctionValues.resize(num_values_linear);
        for (std::size_t i=0; i<num_values_linear; ++i) {
            mLinearShapeFunctionValues[i] = linear_shape_function_values[i];
        }

        if (mHermitianShapeFunctionValues.size() != num_values_linear) mHermitianShapeFunctionValues.resize(num_values_hermitian);
        for (std::size_t i=0; i<num_values_hermitian; ++i) {
            mHermitianShapeFunctionValues[i] = hermitian_shape_function_values[i];
        }

        if (mHermitianShapeFunctionValuesDerivatives.size() != num_values_hermitian_der) mHermitianShapeFunctionValuesDerivatives.resize(num_values_hermitian_der);
        for (std::size_t i=0; i<num_values_hermitian_der; ++i) {
            mHermitianShapeFunctionValuesDerivatives[i] = hermitian_shape_function_derivatives_values[i];
        }

        mProjectionOfPoint = projection_point;
        
        mpInterfaceObject = make_shared<InterfaceGeometryObject>(rInterfaceObject.pGetBaseGeometry());
    }
    
}

void BeamMapperInterfaceInfo::ComputeRotationMatrix()
{
    std::vector<double> axisX;
    std::vector<double> axisY;
    std::vector<double> axisZ;

    axisX.resize(3);
    axisY.resize(3);
    axisZ.resize(3);
    
    const auto p_geom = mpInterfaceObject->pGetBaseGeometry();

    auto temp_v = (*p_geom)[1].Coordinates() - (*p_geom)[0].Coordinates();
    double lengthX = sqrt(temp_v[0]*temp_v[0] + temp_v[1]*temp_v[1] + temp_v[2]*temp_v[2]);
    
    KRATOS_ERROR_IF(lengthX < 0.000001) << "Lenght of the beam is 0.0" << std::endl;
    
    axisX[0] = temp_v[0] / lengthX;
    axisX[1] = temp_v[1] / lengthX;
    axisX[2] = temp_v[2] / lengthX;   

    if (axisX[0] == 1.0 && axisX[1] == 0.0 && axisX[2] == 0.0 ){
        axisY[0] = 0.0;
        axisY[1] = 1.0;
        axisY[2] = 0.0;
        axisZ[0] = 0.0;
        axisZ[1] = 0.0;
        axisZ[2] = 1.0;
    }
    else if (axisX[0] == 0.0 && axisX[1] == 1.0 && axisX[2] == 0.0 ){
        axisY[0] = 0.0;
        axisY[1] = 0.0;
        axisY[2] = 1.0;
        axisZ[0] = 1.0;
        axisZ[1] = 0.0;
        axisZ[2] = 0.0;
    }
    else if (axisX[0] == 0.0 && axisX[1] == 0.0 && axisX[2] == 1.0 ){
        axisY[0] = 0.0;
        axisY[1] = 1.0;
        axisY[2] = 0.0;
        axisZ[0] = 1.0;
        axisZ[1] = 0.0;
        axisZ[2] = 0.0;
    }
    else if (axisX[0] != 0.0 && axisX[1] != 0.0 && axisX[2] == 0.0 ){
        axisY[0] = -axisX[1];
        axisY[1] =  axisX[0];
        axisY[2] =  0.0;
        axisZ[0] = axisX[1]*axisY[2] - axisX[2]*axisY[1]; 
        axisZ[1] = axisX[2]*axisY[0] - axisX[0]*axisY[2];
        axisZ[2] = axisX[0]*axisY[1] - axisX[1]*axisY[0];
    }
    else if (axisX[0] != 0.0 && axisX[1] == 0.0 && axisX[2] != 0.0 ){
        axisY[0] = -axisX[2];
        axisY[1] =  0;
        axisY[2] =  axisX[0];
        axisZ[0] = axisX[1]*axisY[2] - axisX[2]*axisY[1]; 
        axisZ[1] = axisX[2]*axisY[0] - axisX[0]*axisY[2];
        axisZ[2] = axisX[0]*axisY[1] - axisX[1]*axisY[0];
    }
    else if (axisX[0] == 0.0 && axisX[1] != 0.0 && axisX[2] != 0.0){
        axisY[0] =  0;
        axisY[1] = -axisX[2];
        axisY[2] =  axisX[1];
        axisZ[0] = axisX[1]*axisY[2] - axisX[2]*axisY[1]; 
        axisZ[1] = axisX[2]*axisY[0] - axisX[0]*axisY[2];
        axisZ[2] = axisX[0]*axisY[1] - axisX[1]*axisY[0];
    }
    else{
        axisY[0] = 1;
        axisY[1] = 1;
        axisY[2] = (-axisX[0] - axisX[1]) / axisX[2];
        double lenghtY = sqrt(axisY[0]*axisY[0] + axisY[1]*axisY[1] + axisY[2]*axisY[2]);
        axisY[0] = axisY[0]/lenghtY;
        axisY[1] = axisY[1]/lenghtY;
        axisY[2] = axisY[2]/lenghtY;
        
        axisZ[0] = axisX[1]*axisY[2] - axisX[2]*axisY[1]; 
        axisZ[1] = axisX[2]*axisY[0] - axisX[0]*axisY[2];
        axisZ[2] = axisX[0]*axisY[1] - axisX[1]*axisY[0];
    }

    MatrixType _RotationMatrix(3, 3, 0.0);

    for(std::size_t j = 0; j < 3; j++)
    {
        _RotationMatrix(j, 0) = axisX[j];
        _RotationMatrix(j, 1) = axisY[j];
        _RotationMatrix(j, 2) = axisZ[j];
    }

    mRotationMatrixOfBeam = _RotationMatrix;

}

// This function does not do nothing for the beam mapper
void BeamMapperLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    KRATOS_WARNING("BeamMapperLocalSystem")
        << "CalculateAll() was called, but is not implemented for BeamMapperLocalSystem." << std::endl;
}

void BeamMapperLocalSystem::PairingInfo(std::ostream& rOStream, const int EchoLevel) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not initialized!" << std::endl;

    rOStream << "BeamMapperLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 1) {
        rOStream << " at Coordinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
        if (mPairingStatus == MapperLocalSystem::PairingStatus::Approximation) {
            mpNode->SetValue(PAIRING_STATUS, 0);
        } 
    }
}

}  // namespace Kratos.
