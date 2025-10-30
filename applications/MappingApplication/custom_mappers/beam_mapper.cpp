//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
// Main authors:    E. G. Loera Villeda
// Contributor:     Juan I. Camarotti
// See PhD Thesis Tianyang Wang Chapter 5

// System includes

// External includes
#include "utilities/math_utils.h"
#include "geometries/line_3d_2.h"

// Project includes
#include "beam_mapper.h"
#include "mappers/mapper_define.h"
#include "custom_utilities/mapping_matrix_utilities.h"
#include "mapping_application_variables.h"

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

    for (auto& node : p_geom->Points()) {
        node.X() = node.X0();
        node.Y() = node.Y0();
        node.Z() = node.Z0();
    }

    std::ignore = GeometricalProjectionUtilities::FastProjectOnLine(*(p_geom), point_to_proj, projection_point);
    
    Point min_distance_to_nodes;
    double distance_to_node1 = norm_2(projection_point - (*p_geom)[0]);
    double distance_to_node2 = norm_2(projection_point - (*p_geom)[1]);
    
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
        std::ignore = ProjectionUtilities::ProjectOnLine(*p_geom, point_to_proj, mLocalCoordTol, linear_shape_function_values, eq_ids, proj_dist, ComputeApproximation); // Kust to get eq_ids

        p_geom->ShapeFunctionsValues(linear_shape_function_values, local_coords);
        BeamMapperUtilities::HermitianShapeFunctionsValues(hermitian_shape_function_values, hermitian_shape_function_derivatives_values, local_coords);
        mClosestProjectionDistance = proj_dist_nodes;
        mNodeIds = eq_ids;

        const IndexType num_values_linear = linear_shape_function_values.size();
        const IndexType num_values_hermitian = hermitian_shape_function_values.size();
        const IndexType num_values_hermitian_der = hermitian_shape_function_derivatives_values.size();
        
        if (mLinearShapeFunctionValues.size() != num_values_linear) mLinearShapeFunctionValues.resize(num_values_linear);
        for (IndexType i=0; i<num_values_linear; ++i) {
            mLinearShapeFunctionValues[i] = linear_shape_function_values[i];
        }

        if (mHermitianShapeFunctionValues.size() != num_values_linear) mHermitianShapeFunctionValues.resize(num_values_hermitian);
        for (IndexType i=0; i<num_values_hermitian; ++i) {
            mHermitianShapeFunctionValues[i] = hermitian_shape_function_values[i];
        }

        if (mHermitianShapeFunctionValuesDerivatives.size() != num_values_hermitian_der) mHermitianShapeFunctionValuesDerivatives.resize(num_values_hermitian_der);
        for (IndexType i=0; i<num_values_hermitian_der; ++i) {
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

    ProjectionUtilities::PairingIndex pairing_index_linear;

    for (auto& node : p_geom->Points()) {
        node.X() = node.X0();
        node.Y() = node.Y0();
        node.Z() = node.Z0();
    }

    const auto geom_family = p_geom->GetGeometryFamily();
    KRATOS_ERROR_IF(geom_family != GeometryData::KratosGeometryFamily::Kratos_Linear) << "Invalid geometry of the Origin! The geometry should be a beam!";

    // Calculating and storing the shape function values
    // Linear shape functions
    pairing_index_linear = ProjectionUtilities::ProjectOnLine(*p_geom, point_to_proj, mLocalCoordTol, linear_shape_function_values, eq_ids, proj_dist, ComputeApproximation);
    // Hermitian shape functions
    std::ignore = BeamMapperUtilities::ProjectOnLineHermitian(*p_geom, point_to_proj, mLocalCoordTol, hermitian_shape_function_values, hermitian_shape_function_derivatives_values, proj_dist, projection_point);
    const bool is_full_projection = (pairing_index_linear == ProjectionUtilities::PairingIndex::Line_Inside);
    
    if (is_full_projection) {
        SetLocalSearchWasSuccessful();
    } else {
        if (!ComputeApproximation) {
            return;
        } else {
            SetIsApproximation();
        }
    }

    const IndexType num_values_linear = linear_shape_function_values.size();
    const IndexType num_values_hermitian = hermitian_shape_function_values.size();
    const IndexType num_values_hermitian_der = hermitian_shape_function_derivatives_values.size();


    if (pairing_index_linear > mPairingIndex || (pairing_index_linear == mPairingIndex && proj_dist < mClosestProjectionDistance)) {
        mPairingIndex = pairing_index_linear;
        mClosestProjectionDistance = proj_dist;
        mNodeIds = eq_ids;

        if (mLinearShapeFunctionValues.size() != num_values_linear) mLinearShapeFunctionValues.resize(num_values_linear);
        for (IndexType i=0; i<num_values_linear; ++i) {
            mLinearShapeFunctionValues[i] = linear_shape_function_values[i];
        }

        if (mHermitianShapeFunctionValues.size() != num_values_linear) mHermitianShapeFunctionValues.resize(num_values_hermitian);
        for (IndexType i=0; i<num_values_hermitian; ++i) {
            mHermitianShapeFunctionValues[i] = hermitian_shape_function_values[i];
        }

        if (mHermitianShapeFunctionValuesDerivatives.size() != num_values_hermitian_der) mHermitianShapeFunctionValuesDerivatives.resize(num_values_hermitian_der);
        for (IndexType i=0; i<num_values_hermitian_der; ++i) {
            mHermitianShapeFunctionValuesDerivatives[i] = hermitian_shape_function_derivatives_values[i];
        }

        mProjectionOfPoint = projection_point;
        
        mpInterfaceObject = make_shared<InterfaceGeometryObject>(rInterfaceObject.pGetBaseGeometry());
    }
    
}

void BeamMapperInterfaceInfo::ComputeRotationMatrix()
{
    std::vector<double> axis_X;
    std::vector<double> axis_Y;
    std::vector<double> axis_Z;

    axis_X.resize(3);
    axis_Y.resize(3);
    axis_Z.resize(3);
    
    const auto p_geom = mpInterfaceObject->pGetBaseGeometry();

    auto temp_v = (*p_geom)[1].Coordinates() - (*p_geom)[0].Coordinates();
    double length_X = sqrt(temp_v[0]*temp_v[0] + temp_v[1]*temp_v[1] + temp_v[2]*temp_v[2]);
    
    KRATOS_ERROR_IF(length_X < 0.000001) << "Lenght of the beam is 0.0" << std::endl;
    
    axis_X[0] = temp_v[0] / length_X;
    axis_X[1] = temp_v[1] / length_X;
    axis_X[2] = temp_v[2] / length_X;   

    if (axis_X[0] == 1.0 && axis_X[1] == 0.0 && axis_X[2] == 0.0 ){
        axis_Y[0] = 0.0;
        axis_Y[1] = 1.0;
        axis_Y[2] = 0.0;
        axis_Z[0] = 0.0;
        axis_Z[1] = 0.0;
        axis_Z[2] = 1.0;
    }
    else if (axis_X[0] == 0.0 && axis_X[1] == 1.0 && axis_X[2] == 0.0 ){
        axis_Y[0] = 0.0;
        axis_Y[1] = 0.0;
        axis_Y[2] = 1.0;
        axis_Z[0] = 1.0;
        axis_Z[1] = 0.0;
        axis_Z[2] = 0.0;
    }
    else if (axis_X[0] == 0.0 && axis_X[1] == 0.0 && axis_X[2] == 1.0 ){
        axis_Y[0] = 0.0;
        axis_Y[1] = 1.0;
        axis_Y[2] = 0.0;
        axis_Z[0] = 1.0;
        axis_Z[1] = 0.0;
        axis_Z[2] = 0.0;
    }
    else if (axis_X[0] != 0.0 && axis_X[1] != 0.0 && axis_X[2] == 0.0 ){
        axis_Y[0] = -axis_X[1];
        axis_Y[1] =  axis_X[0];
        axis_Y[2] =  0.0;
        axis_Z[0] = axis_X[1]*axis_Y[2] - axis_X[2]*axis_Y[1]; 
        axis_Z[1] = axis_X[2]*axis_Y[0] - axis_X[0]*axis_Y[2];
        axis_Z[2] = axis_X[0]*axis_Y[1] - axis_X[1]*axis_Y[0];
    }
    else if (axis_X[0] != 0.0 && axis_X[1] == 0.0 && axis_X[2] != 0.0 ){
        axis_Y[0] = -axis_X[2];
        axis_Y[1] =  0;
        axis_Y[2] =  axis_X[0];
        axis_Z[0] = axis_X[1]*axis_Y[2] - axis_X[2]*axis_Y[1]; 
        axis_Z[1] = axis_X[2]*axis_Y[0] - axis_X[0]*axis_Y[2];
        axis_Z[2] = axis_X[0]*axis_Y[1] - axis_X[1]*axis_Y[0];
    }
    else if (axis_X[0] == 0.0 && axis_X[1] != 0.0 && axis_X[2] != 0.0){
        axis_Y[0] =  0;
        axis_Y[1] = -axis_X[2];
        axis_Y[2] =  axis_X[1];
        axis_Z[0] = axis_X[1]*axis_Y[2] - axis_X[2]*axis_Y[1]; 
        axis_Z[1] = axis_X[2]*axis_Y[0] - axis_X[0]*axis_Y[2];
        axis_Z[2] = axis_X[0]*axis_Y[1] - axis_X[1]*axis_Y[0];
    }
    else{
        axis_Y[0] = 1;
        axis_Y[1] = 1;
        axis_Y[2] = (-axis_X[0] - axis_X[1]) / axis_X[2];
        double length_Y = sqrt(axis_Y[0]*axis_Y[0] + axis_Y[1]*axis_Y[1] + axis_Y[2]*axis_Y[2]);
        axis_Y[0] = axis_Y[0]/length_Y;
        axis_Y[1] = axis_Y[1]/length_Y;
        axis_Y[2] = axis_Y[2]/length_Y;
        
        axis_Z[0] = axis_X[1]*axis_Y[2] - axis_X[2]*axis_Y[1]; 
        axis_Z[1] = axis_X[2]*axis_Y[0] - axis_X[0]*axis_Y[2];
        axis_Z[2] = axis_X[0]*axis_Y[1] - axis_X[1]*axis_Y[0];
    }

    MatrixType rotation_matrix(3, 3, 0.0);

    for(IndexType j = 0; j < 3; j++)
    {
        rotation_matrix(j, 0) = axis_X[j];
        rotation_matrix(j, 1) = axis_Y[j];
        rotation_matrix(j, 2) = axis_Z[j];
    }

    mRotationMatrixOfBeam = rotation_matrix;

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

// ************ BeamMapper function definitions

template<class TSparseSpace, class TDenseSpace>
void BeamMapper<TSparseSpace, TDenseSpace>::InitializeInformationBeams(const Variable<array_1d<double, 3>>& rOriginVariablesDisplacements,
                                                                       const Variable<array_1d<double, 3>>& rOriginVariablesRotations,
                                                                       const Variable<array_1d<double, 3>>& rDestinationVariableDisplacement)
{   
    for( auto& r_local_sys : mMapperLocalSystems )
    {   
        if( r_local_sys->HasInterfaceInfo())
        {
            auto beam_sys = dynamic_cast<BeamMapperLocalSystem*>(r_local_sys.get());
            KRATOS_ERROR_IF_NOT(beam_sys) << "Expected BeamMapperLocalSystem!" << std::endl;

            MatrixType rotation_matrix_G_B(3, 3);
            VectorType translation_vector_B_P(3);
            VectorType linear_shape_values(2);
            VectorType hermitian_shape_values(4);
            VectorType hermitian_der_shape_values(4);
            GeometryType r_geom;
            NodePointerType p_node;

            beam_sys->CalculateRotationMatrixInterfaceInfos(rotation_matrix_G_B,
                                                               translation_vector_B_P,
                                                               linear_shape_values,
                                                               hermitian_shape_values,
                                                               hermitian_der_shape_values,
                                                               r_geom,
                                                               p_node);


            KRATOS_ERROR_IF_NOT(p_node) << "Node is a nullptr"<< std::endl;

            const std::vector<std::string> var_comps{"_X", "_Y", "_Z"};
            VectorType displacement_node_1_G(3); //Expresses in global coordinates
            VectorType displacement_node_2_G(3); //Expresses in global coordinates
            VectorType rotation_node_1_G(3); //Expresses in global coordinates
            VectorType rotation_node_2_G(3); //Expresses in global coordinates

            VectorType displacement_node_1_B(3); //Expresses in beam coordinates
            VectorType displacement_node_2_B(3); //Expresses in beam coordinates
            VectorType rotation_node_1_B(3); //Expresses in beam coordinates
            VectorType rotation_node_2_B(3); //Expresses in beam coordinates

            IndexType k = 0;

            for (const auto& var_ext : var_comps)
            {
                const auto& var_origin_disp = KratosComponents<ComponentVariableType>::Get(rOriginVariablesDisplacements.Name() + var_ext);
                displacement_node_1_G(k) = r_geom[0].FastGetSolutionStepValue(var_origin_disp);
                displacement_node_2_G(k) = r_geom[1].FastGetSolutionStepValue(var_origin_disp);


                const auto& var_origin_rot = KratosComponents<ComponentVariableType>::Get(rOriginVariablesRotations.Name() + var_ext);
                rotation_node_1_G(k) = r_geom[0].FastGetSolutionStepValue(var_origin_rot);
                rotation_node_2_G(k) = r_geom[1].FastGetSolutionStepValue(var_origin_rot);
                k++;
            }

            MatrixType rotation_matrix_B_G( 3, 3 );
            double determinant;
            MathUtils<double>::InvertMatrix3(rotation_matrix_G_B, rotation_matrix_B_G, determinant );
            
            TDenseSpace::Mult( rotation_matrix_B_G, displacement_node_1_G, displacement_node_1_B );
            TDenseSpace::Mult( rotation_matrix_B_G, rotation_node_1_G, rotation_node_1_B );
            TDenseSpace::Mult( rotation_matrix_B_G, displacement_node_2_G, displacement_node_2_B );
            TDenseSpace::Mult( rotation_matrix_B_G, rotation_node_2_G, rotation_node_2_B );
            
            // Initializing matrix of shape functions
            MatrixType shape_functions_matrix(6, 12, 0.0);

            VectorType beamVector = r_geom[1].Coordinates() - r_geom[0].Coordinates();
            double length_beam_vector = norm_2(beamVector);
            shape_functions_matrix(0 , 0) = linear_shape_values(0);
            shape_functions_matrix(0 , 6) = linear_shape_values(1);
            shape_functions_matrix(1 , 1) = hermitian_shape_values(0);
            shape_functions_matrix(1 , 5) = hermitian_shape_values(1) * length_beam_vector;
            shape_functions_matrix(1 , 7) = hermitian_shape_values(2);
            shape_functions_matrix(1 , 11) = hermitian_shape_values(3) * length_beam_vector;
            shape_functions_matrix(2 , 2) = hermitian_shape_values(0);
            shape_functions_matrix(2 , 4) = -hermitian_shape_values(1) * length_beam_vector;
            shape_functions_matrix(2 , 8) = hermitian_shape_values(2);
            shape_functions_matrix(2 , 10) = -hermitian_shape_values(3) * length_beam_vector;
            
            shape_functions_matrix(3 , 3) = linear_shape_values(0);
            shape_functions_matrix(3 , 9) = linear_shape_values(1);
            shape_functions_matrix(4 , 2) = -hermitian_der_shape_values(0) / length_beam_vector;
            shape_functions_matrix(4 , 4) = hermitian_der_shape_values(1);
            shape_functions_matrix(4 , 8) = -hermitian_der_shape_values(2) / length_beam_vector;
            shape_functions_matrix(4 , 10) = -hermitian_der_shape_values(3);
            shape_functions_matrix(5 , 1) = hermitian_der_shape_values(0) / length_beam_vector;
            shape_functions_matrix(5 , 5) = hermitian_der_shape_values(1);
            shape_functions_matrix(5 , 7) = hermitian_der_shape_values(2) / length_beam_vector;
            shape_functions_matrix(5 , 11) = hermitian_der_shape_values(3);
            
            VectorType DOF_vector(12);
            for (IndexType i = 0; i < 3; i++){
                DOF_vector(i) = displacement_node_1_B(i);
                DOF_vector(i + 3) = rotation_node_1_B(i);
                DOF_vector(i + 6) = displacement_node_2_B(i);
                DOF_vector(i + 9) = rotation_node_2_B(i);
            }

            VectorType displacements_rotations_P(6);
            TDenseSpace::Mult(shape_functions_matrix, DOF_vector, displacements_rotations_P);

            VectorType axis_X(3, 0.0);
            VectorType axis_Y(3, 0.0);
            VectorType axis_Z(3, 0.0);
            axis_X(0) = 1.0;
            axis_Y(1) = 1.0;
            axis_Z(2) = 1.0; 

            MatrixType rotation_X(3, 3);
            MatrixType rotation_Y(3, 3);
            MatrixType rotation_Z(3, 3);

            CalculateRotationMatrixWithAngle( axis_X, displacements_rotations_P(3), rotation_X);
            CalculateRotationMatrixWithAngle( axis_Y, displacements_rotations_P(4), rotation_Y);
            CalculateRotationMatrixWithAngle( axis_Z, displacements_rotations_P(5), rotation_Z);

            MatrixType rotation_matrix_P(3, 3);
            MatrixType tmp_rotation(3, 3);
            tmp_rotation = prod(rotation_X ,rotation_Z);
            rotation_matrix_P = prod(tmp_rotation, rotation_Y);

            VectorType TranslationVector_P(3);
            TranslationVector_P(0) = displacements_rotations_P(0);
            TranslationVector_P(1) = displacements_rotations_P(1);
            TranslationVector_P(2) = displacements_rotations_P(2);

            // Rigid body operation 
            // phi_G( X ) = R_G_B * R_P * (R_G_B ^ T) * X - R_G_B * R_P * (R_G_B ^ T) * t_B_P + R_G_B * t_B + t_B_P

            MatrixType matrix_product(3, 3); 
            MathUtils<double>::BDBtProductOperation(matrix_product, rotation_matrix_P, rotation_matrix_G_B);
            
            VectorType first_product(3), second_product(3), third_product(3), phi_G(3), displacement(3), X(3);

            X(0) = (r_local_sys->Coordinates())[0];
            X(1) = (r_local_sys->Coordinates())[1];
            X(2) = (r_local_sys->Coordinates())[2];
            TDenseSpace::Mult(matrix_product, X, first_product );
            TDenseSpace::Mult(matrix_product, translation_vector_B_P, second_product );
            TDenseSpace::Mult(rotation_matrix_G_B, TranslationVector_P, third_product );

            phi_G(0) = first_product(0) - second_product(0) + third_product(0) + translation_vector_B_P(0);
            phi_G(1) = first_product(1) - second_product(1) + third_product(1) + translation_vector_B_P(1);
            phi_G(2) = first_product(2) - second_product(2) + third_product(2) + translation_vector_B_P(2); 

            displacement(0) = phi_G(0 ) - X(0);
            displacement(1) = phi_G(1) - X(1);
            displacement(2) = phi_G(2) - X(2);

            IndexType c = 0;
            for( const auto& var_ext : var_comps )
            {
                const auto& var_destination_disp = KratosComponents<ComponentVariableType>::Get(rDestinationVariableDisplacement.Name() + var_ext);
                p_node->FastGetSolutionStepValue(var_destination_disp) = displacement( c );
                c++;
            }
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamMapper<TSparseSpace, TDenseSpace>::InitializeInformationBeamsCorotation(const Variable< array_1d<double, 3> >& rOriginVariablesDisplacements,
                                                                                 const Variable< array_1d<double, 3> >& rOriginVariablesRotations,
                                                                                 const Variable< array_1d<double, 3> >& rDestinationVariableDisplacement)
{
    for( auto& r_local_sys : mMapperLocalSystems )
    {   
        if( r_local_sys->HasInterfaceInfo())
        {
            auto beam_sys = dynamic_cast<BeamMapperLocalSystem*>(r_local_sys.get());
            KRATOS_ERROR_IF_NOT(beam_sys) << "Expected BeamMapperLocalSystem!" << std::endl;

            MatrixType rotation_matrix_G_B(3, 3);
            VectorType translation_vector_B_P(3);
            VectorType linear_shape_values(2);
            VectorType hermitian_shape_values(4);
            VectorType hermitian_der_shape_values(4);
            GeometryType r_geom;
            NodePointerType p_node;

            beam_sys->CalculateRotationMatrixInterfaceInfos(rotation_matrix_G_B,
                                                               translation_vector_B_P,
                                                               linear_shape_values,
                                                               hermitian_shape_values,
                                                               hermitian_der_shape_values,
                                                               r_geom,
                                                               p_node);

            KRATOS_ERROR_IF_NOT(p_node) << "Node is a nullptr"<< std::endl;

            const std::vector<std::string> var_comps{"_X", "_Y", "_Z"};
            VectorType displacement_node_1_G(3); // Expresses in global coordinates
            VectorType displacement_node_2_G(3); // Expresses in global coordinates
            VectorType rotation_node_1_G(3); // Expresses in global coordinates
            VectorType rotation_node_2_G(3); // Expresses in global coordinates

            VectorType displacement_node_1_B(3); // Expresses in beam coordinates
            VectorType displacement_node_2_B(3); // Expresses in beam coordinates
            VectorType rotation_node_1_B(3); // Expresses in beam coordinates
            VectorType rotation_node_2_B(3); // Expresses in beam coordinates

            IndexType k = 0;

            for (const auto& var_ext : var_comps)
            {
                const auto& var_origin_disp = KratosComponents<ComponentVariableType>::Get(rOriginVariablesDisplacements.Name() + var_ext);
                displacement_node_1_G(k) = r_geom[0].FastGetSolutionStepValue(var_origin_disp);
                displacement_node_2_G(k) = r_geom[1].FastGetSolutionStepValue(var_origin_disp);


                const auto& var_origin_rot = KratosComponents<ComponentVariableType>::Get(rOriginVariablesRotations.Name() + var_ext);
                rotation_node_1_G(k) = r_geom[0].FastGetSolutionStepValue(var_origin_rot);
                rotation_node_2_G(k) = r_geom[1].FastGetSolutionStepValue(var_origin_rot);
                k++;
            }

            MatrixType rotation_matrix_B_G( 3, 3 );
            double determinant;
            MathUtils<double>::InvertMatrix3(rotation_matrix_G_B, rotation_matrix_B_G, determinant );

            // Transforming the nodal displacements to the BCS
            TDenseSpace::Mult( rotation_matrix_B_G, displacement_node_1_G, displacement_node_1_B );
            TDenseSpace::Mult( rotation_matrix_B_G, displacement_node_2_G, displacement_node_2_B );
            
            // Transforming the nodal rotations to the BCS
            MatrixType rotation_G_1(3, 3);
            double angle_1 = norm_2(rotation_node_1_G);
            if (angle_1 != 0.0)
                rotation_node_1_G /= angle_1;
            CalculateRotationMatrixWithAngle(rotation_node_1_G, angle_1, rotation_G_1);
            MatrixType rotation_B_1(3, 3);
            MathUtils<double>::BtDBProductOperation(rotation_B_1, rotation_G_1,rotation_matrix_G_B);
            
            MatrixType rotation_G_2(3, 3);
            double _angle2 = norm_2(rotation_node_2_G);
            if (_angle2 != 0.0)
                rotation_node_2_G /= _angle2;
            CalculateRotationMatrixWithAngle(rotation_node_2_G, _angle2, rotation_G_2);
            MatrixType Rotation_B_2(3, 3);
            MathUtils<double>::BtDBProductOperation(Rotation_B_2, rotation_G_2,rotation_matrix_G_B);

            // Calculating R_d and t_d for Phi_d 
            VectorType e_x_d_G(3), e_x_d_B(3);
            e_x_d_G = r_geom[1].Coordinates() + displacement_node_2_G - (r_geom[0].Coordinates() + displacement_node_1_G); // this vector is described in global system 
            e_x_d_G /= norm_2(e_x_d_G);
            TDenseSpace::Mult(rotation_matrix_B_G, e_x_d_G, e_x_d_B); // transforming e_x_d to the beam coordinate system
            e_x_d_B /= norm_2(e_x_d_B);
            
            VectorType e_x(3, 0.0), n_d(3, 0.0);
            MatrixType R_d(3, 3), I(3, 3, 0.0), R_(3, 3, 0.0), I_2nT(3, 3, 0.0);
            e_x(0) = 1.0;
            I(0,0) = 1.0;
            I(1,1) = 1.0;
            I(2,2) = 1.0;
            R_(0, 0) = -1.0;
            R_(1, 1) = 1.0;
            R_(2, 2) = 1.0;
            n_d =  e_x + e_x_d_B;
            n_d /= norm_2(n_d);
            I_2nT = I - 2 * MathUtils<double>::TensorProduct3(n_d, n_d);
            R_d = prod(I_2nT, R_); // Equation 5.22 from Tianyang Wang's dissertation
            
            // Calculating the translation vector t_d of the section in the BCS
            VectorType t_d(3);
            t_d(0) = linear_shape_values(0) * displacement_node_1_B(0) + linear_shape_values(1) * displacement_node_2_B(0);
            t_d(1) = linear_shape_values(0) * displacement_node_1_B(1) + linear_shape_values(1) * displacement_node_2_B(1);
            t_d(2) = linear_shape_values(0) * displacement_node_1_B(2) + linear_shape_values(1) * displacement_node_2_B(2);

            // Calculating R_s (average rotation of the section)
            MatrixType Rs_Rl1(3, 3), Rs_Rl2(3, 3), R_d_T(3, 3);
            MathUtils<double>::InvertMatrix3(R_d, R_d_T, determinant);
            Rs_Rl1 = prod(R_d_T, rotation_B_1); 
            Rs_Rl2 = prod(R_d_T, Rotation_B_2);
            VectorType v_Rs_Rl1(3), v_Rs_Rl2(3);
            GetRotationVector(Rs_Rl1, v_Rs_Rl1);
            GetRotationVector(Rs_Rl2, v_Rs_Rl2);

            double theta_s = (v_Rs_Rl1(0) + v_Rs_Rl2(0)) / 2;
            MatrixType R_s(3, 3);
            VectorType e_x_s(3, 0.0);
            e_x_s(0) = 1.0;
            CalculateRotationMatrixWithAngle(e_x_s, theta_s, R_s);

            // Calculating rotation vector l on node 1 and 2 of the beam
            MatrixType Rl1(3, 3), Rl2(3, 3), R_s_T(3, 3);
            VectorType v_Rl1(3), v_Rl2(3);
            MathUtils<double>::InvertMatrix3(R_s, R_s_T, determinant);
            Rl1 = prod(R_s_T, Rs_Rl1);
            Rl2 = prod(R_s_T, Rs_Rl2);
            GetRotationVector(Rl1, v_Rl1);
            GetRotationVector(Rl2, v_Rl2); 

            MatrixType shape_functions_matrix(6, 12, 0.0);
            VectorType corBeamVector(3, 0.0);
            corBeamVector = r_geom[1].Coordinates() + displacement_node_2_G - (r_geom[0].Coordinates() + displacement_node_1_G); // this vector is described in global system 
            double length_beam_vector = norm_2(corBeamVector);            

            shape_functions_matrix(0 , 0) = linear_shape_values(0);
            shape_functions_matrix(0 , 6) = linear_shape_values(1);
            shape_functions_matrix(1 , 1) = hermitian_shape_values(0);
            shape_functions_matrix(1 , 5) = hermitian_shape_values(1) * length_beam_vector;
            shape_functions_matrix(1 , 7) = hermitian_shape_values(2);
            shape_functions_matrix(1 , 11) = hermitian_shape_values(3) * length_beam_vector;
            shape_functions_matrix(2 , 2) = hermitian_shape_values(0);
            shape_functions_matrix(2 , 4) = -hermitian_shape_values(1) * length_beam_vector;
            shape_functions_matrix(2 , 8) = hermitian_shape_values(2);
            shape_functions_matrix(2 , 10) = -hermitian_shape_values(3) * length_beam_vector;
            
            shape_functions_matrix(3 , 3) = linear_shape_values(0);
            shape_functions_matrix(3 , 9) = linear_shape_values(1);
            shape_functions_matrix(4 , 2) = -hermitian_der_shape_values(0) / length_beam_vector;
            shape_functions_matrix(4 , 4) = hermitian_der_shape_values(1);
            shape_functions_matrix(4 , 8) = -hermitian_der_shape_values(2) / length_beam_vector;
            shape_functions_matrix(4 , 10) = hermitian_der_shape_values(3);
            shape_functions_matrix(5 , 1) = hermitian_der_shape_values(0) / length_beam_vector;
            shape_functions_matrix(5 , 5) = hermitian_der_shape_values(1);
            shape_functions_matrix(5 , 7) = hermitian_der_shape_values(2) / length_beam_vector;
            shape_functions_matrix(5 , 11) = hermitian_der_shape_values(3); 

            VectorType DOF_vector_l(12);
            for (IndexType i = 0; i < 3; i++){
                DOF_vector_l(i) = 0.0;
                DOF_vector_l(i + 3) = v_Rl1(i);
                DOF_vector_l(i + 6) = 0.0;
                DOF_vector_l(i + 9) = v_Rl2(i);
            }

            VectorType displacements_rotations_L(6);
            TDenseSpace::Mult(shape_functions_matrix, DOF_vector_l, displacements_rotations_L);
            
            // Constructing phi_l
            VectorType axis_l(3), t_l(3);
            MatrixType R_l(3, 3), R_l_temp(3, 3), R_x(3, 3), R_y(3, 3), R_z(3, 3);
            t_l(0) = displacements_rotations_L(0); // t_l
            t_l(1) = displacements_rotations_L(1); // t_l
            t_l(2) = displacements_rotations_L(2); // t_l

            MatrixType Rx(3, 3), Ry(3, 3), Rz(3, 3), R(3, 3), R_temp(3, 3);
            VectorType e_1(3, 0.0), e_2(3, 0.0), e_3(3, 0.0);
            e_1(0) = 1.0;
            e_2(1) = 1.0;
            e_3(2) = 1.0;

            CalculateRotationMatrixWithAngle(e_1, displacements_rotations_L(3), Rx);
            CalculateRotationMatrixWithAngle(e_2, displacements_rotations_L(4), Ry);
            CalculateRotationMatrixWithAngle(e_3, displacements_rotations_L(5), Rz);
            R_temp = prod(Rx, Ry); 
            R_l = prod(Rz, R_temp);

            // Calculating R_P and t_P
            MatrixType R_P(3, 3), R_P_temp(3, 3);


            R_P_temp = prod(R_d, R_s);
            R_P = prod(R_P_temp, R_l);

            VectorType t_P(3), t_P_temp(3, 0.0);
            TDenseSpace::Mult(R_P_temp, t_l, t_P_temp);
            t_P = t_P_temp + t_d;

            // Rigid body operation 
            // phi_G( X ) = R_G_B * R_P * (R_G_B ^ T) * X - R_G_B * R_P * (R_G_B ^ T) * t_B_P + R_G_B * t_B + t_B_P

            MatrixType matrix_product(3, 3); 
            MathUtils<double>::BDBtProductOperation(matrix_product, R_P, rotation_matrix_G_B);
            
            VectorType first_product(3), second_product(3), third_product(3), phi_G(3), displacement(3), X(3);

            X(0) = (r_local_sys->Coordinates())[0];
            X(1) = (r_local_sys->Coordinates())[1];
            X(2) = (r_local_sys->Coordinates())[2];
            // Getting the first product of the rigid body operation
            TDenseSpace::Mult(matrix_product, X, first_product );
            // Getting the second product of the rigid body operation
            TDenseSpace::Mult(matrix_product, translation_vector_B_P, second_product );
            // Getting the third product of the rigid body operation
            TDenseSpace::Mult(rotation_matrix_G_B, t_P, third_product );

            phi_G(0) = first_product(0) - second_product(0) + third_product(0) + translation_vector_B_P(0);
            phi_G(1) = first_product(1) - second_product(1) + third_product(1) + translation_vector_B_P(1);
            phi_G(2) = first_product(2) - second_product(2) + third_product(2) + translation_vector_B_P(2); 

            displacement(0) = phi_G(0 ) - X(0);
            displacement(1) = phi_G(1) - X(1);
            displacement(2) = phi_G(2) - X(2);

            // For conservative mapping
            VectorType rotation_of_section(3);
            GetRotationVector(matrix_product, rotation_of_section);
            beam_sys->SaveRotationVectorValue(rotation_of_section);

            IndexType c = 0;
            for( const auto& var_ext : var_comps )
            {
                const auto& var_destination_disp = KratosComponents<ComponentVariableType>::Get(rDestinationVariableDisplacement.Name() + var_ext);
                p_node->FastGetSolutionStepValue(var_destination_disp) = displacement( c );
                c++;
            }
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamMapper<TSparseSpace, TDenseSpace>::InitializeInformationBeamsInverse(const Variable<array_1d<double, 3>>& rOriginVariablesForces,
                                                                              const Variable<array_1d<double, 3>>& rOriginVariablesMoments,
                                                                              const Variable<array_1d<double, 3>>& rDestinationVariableForces,
                                                                              const Kratos::Flags& rMappingOptions)
{   
    const double factor = rMappingOptions.Is(MapperFlags::SWAP_SIGN) ? -1.0 : 1.0;
    for( auto& r_local_sys : mMapperLocalSystems )
    {   
        if( r_local_sys->HasInterfaceInfo())
        {           
            auto beam_sys = dynamic_cast<BeamMapperLocalSystem*>(r_local_sys.get());
            KRATOS_ERROR_IF_NOT(beam_sys) << "Expected BeamMapperLocalSystem!" << std::endl;

            MatrixType rotation_matrix_G_B(3, 3);
            VectorType translation_vector_B_P(3);
            VectorType linear_shape_values(2);
            VectorType hermitian_shape_values(4);
            VectorType hermitian_der_shape_values(4);
            GeometryType r_geom; // beam line
            NodePointerType p_node; // surface node

            beam_sys->CalculateRotationMatrixInterfaceInfos(rotation_matrix_G_B,
                                                               translation_vector_B_P,
                                                               linear_shape_values,
                                                               hermitian_shape_values,
                                                               hermitian_der_shape_values,
                                                               r_geom,
                                                               p_node);

            KRATOS_ERROR_IF_NOT(p_node) << "Node is a nullptr"<< std::endl;

            const std::vector<std::string> var_comps{"_X", "_Y", "_Z"};
            VectorType surface_force(3), surface_moment(3); 
            VectorType point_P(3), point_Q(3), distance_PQ(3), distance_PQ_moved(3);

            IndexType k = 0;

            for (const auto& var_ext : var_comps)
            {
                const auto& var_destination_force = KratosComponents<ComponentVariableType>::Get(rDestinationVariableForces.Name() + var_ext);
                surface_force(k) = p_node->FastGetSolutionStepValue(var_destination_force);
                k++;
            }
            VectorType rotation_vector_of_section(3);
            MatrixType rotation_matrix_of_section(3, 3);
            beam_sys->GetValue(rotation_vector_of_section);

            double angle = norm_2(rotation_vector_of_section);
            if (angle != 0.0){
                rotation_vector_of_section(0) = rotation_vector_of_section(0) / angle;
                rotation_vector_of_section(1) = rotation_vector_of_section(1) / angle;
                rotation_vector_of_section(2) = rotation_vector_of_section(2) / angle;
            }

            CalculateRotationMatrixWithAngle(rotation_vector_of_section, angle, rotation_matrix_of_section);
            point_P = translation_vector_B_P;
            point_Q = p_node->Coordinates();
            distance_PQ = point_Q - point_P;
            TDenseSpace::Mult(rotation_matrix_of_section, distance_PQ, distance_PQ_moved);

            surface_moment(0) = distance_PQ_moved(1) * surface_force(2) - distance_PQ_moved(2) * surface_force(1);
            surface_moment(1) = distance_PQ_moved(2) * surface_force(0) - distance_PQ_moved(0) * surface_force(2);
            surface_moment(2) = distance_PQ_moved(0) * surface_force(1) - distance_PQ_moved(1) * surface_force(0);

            IndexType c = 0;
            for (const auto& var_ext : var_comps)
            {
                const auto& var_origin_force = KratosComponents<ComponentVariableType>::Get(rOriginVariablesForces.Name() + var_ext);
                const auto& var_origin_moment = KratosComponents<ComponentVariableType>::Get(rOriginVariablesMoments.Name() + var_ext);

                r_geom[0].FastGetSolutionStepValue(var_origin_force)  +=  linear_shape_values(0) * surface_force(c) * factor;
                r_geom[0].FastGetSolutionStepValue(var_origin_moment) +=  linear_shape_values(0) * surface_moment(c) * factor;
                r_geom[1].FastGetSolutionStepValue(var_origin_force)  +=  linear_shape_values(1) * surface_force(c) * factor;
                r_geom[1].FastGetSolutionStepValue(var_origin_moment) +=  linear_shape_values(1) * surface_moment(c) * factor; 
                c++;
            }
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
void BeamMapper<TSparseSpace, TDenseSpace>::CalculateRotationMatrixWithAngle( VectorType& rAxis, double& rAngle , MatrixType& rRotationMatrix)
{
    rRotationMatrix(0, 0) = cos( rAngle ) + pow(rAxis(0), 2) * (1 - cos( rAngle ));
    rRotationMatrix(0, 1) = rAxis(0) * rAxis(1) * (1 - cos( rAngle )) - rAxis(2)*sin(rAngle);
    rRotationMatrix(0, 2) = rAxis(0)*rAxis(2)*(1-cos(rAngle)) + rAxis(1)*sin(rAngle);

    rRotationMatrix(1, 0) = rAxis(0)*rAxis(1)*(1-cos(rAngle)) + rAxis(2)*sin(rAngle);
    rRotationMatrix(1, 1) = cos(rAngle) + pow(rAxis(1), 2)*(1 - cos(rAngle));
    rRotationMatrix(1, 2) = rAxis(1)*rAxis(2)*(1-cos(rAngle)) - rAxis(0)*sin(rAngle);

    rRotationMatrix(2, 0) = rAxis(0)*rAxis(2)*(1-cos(rAngle)) - rAxis(1)*sin(rAngle);
    rRotationMatrix(2, 1) = rAxis(1)*rAxis(2)*(1-cos(rAngle)) + rAxis(0)*sin(rAngle);
    rRotationMatrix(2, 2) = cos(rAngle) + pow(rAxis(2),2)*(1-cos(rAngle));
}

template<class TSparseSpace, class TDenseSpace>
void BeamMapper<TSparseSpace, TDenseSpace>::GetRotationVector(const MatrixType& rRotationMatrix, VectorType& rRotationVector) 
{
    // see Non-linear Modeling and Analysis of Solids and Structures (Steen Krenk 2009) P52
    double angle = rRotationMatrix(0, 0) + rRotationMatrix(1, 1) + rRotationMatrix(2, 2) - 1.0;
    const double pi = 3.14159265358979323846;

    angle /= 2.0;
    if (angle > 1.0)
        angle = 1.0;
    else if (angle < -1.0)
        angle = -1.0;

    angle = acos(angle); // between 0 and pi

    const double EPS = 1E-6;
    if (angle < EPS) {
        rRotationVector(0) = 0.0;
        rRotationVector(1) = 0.0;
        rRotationVector(2) = 0.0;

        return;
    } else if ((pi - angle) < EPS) {
        const double product11 = (rRotationMatrix(0,0) + 1.0) / 2.0;
        const double product22 = (rRotationMatrix(1,1) + 1.0) / 2.0;
        const double product33 = (rRotationMatrix(2,2) + 1.0) / 2.0;
        const double product12 = (rRotationMatrix(0,1) + 1.0) / 2.0;
        const double product23 = (rRotationMatrix(1,2) + 1.0) / 2.0;
        const double product13 = (rRotationMatrix(0,2) + 1.0) / 2.0;
        const double tmp1 = sqrt(product11);
        const double tmp2 = sqrt(product22);
        const double tmp3 = sqrt(product33);

        { // case 1 +++:
            rRotationVector(0) = tmp1;
            rRotationVector(1) = tmp2;
            rRotationVector(2) = tmp3;
            const double tmp12 = rRotationVector(0) * rRotationVector(1);
            const double tmp13 = rRotationVector(0) * rRotationVector(2);
            const double tmp23 = rRotationVector(1) * rRotationVector(2);
            if (std::abs(tmp12) < EPS || std::abs(tmp12 - product12) < std::abs(tmp12 + product12))
                if (std::abs(tmp13) < EPS || std::abs(tmp13 - product13) < std::abs(tmp13 + product13))
                    if (std::abs(tmp23) < EPS || std::abs(tmp23 - product23) < std::abs(tmp23 + product23)) {
                        rRotationVector(0) *= pi;
                        rRotationVector(1) *= pi;
                        rRotationVector(2) *= pi;
                        return;
                    }
        }
        { // case 2 +--:
            rRotationVector(0) = tmp1;
            rRotationVector(1) = -tmp2;
            rRotationVector(2) = -tmp3;
            const double tmp12 = rRotationVector[0] * rRotationVector[1];
            const double tmp13 = rRotationVector[0] * rRotationVector[2];
            const double tmp23 = rRotationVector[1] * rRotationVector[2];
            if (std::abs(tmp12) < EPS || std::abs(tmp12 - product12) < std::abs(tmp12 + product12))
                if (std::abs(tmp13) < EPS || std::abs(tmp13 - product13) < std::abs(tmp13 + product13))
                    if (std::abs(tmp23) < EPS || std::abs(tmp23 - product23) < std::abs(tmp23 + product23)) {
                        rRotationVector(0) *= pi;
                        rRotationVector(1) *= pi;
                        rRotationVector(2) *= pi;
                        return;
                    }
        }
        { // case 3 -+-:
            rRotationVector(0) = -tmp1;
            rRotationVector(1) = tmp2;
            rRotationVector(2) = -tmp3;
            const double tmp12 = rRotationVector(0) * rRotationVector(1);
            const double tmp13 = rRotationVector(0) * rRotationVector(2);
            const double tmp23 = rRotationVector(1) * rRotationVector(2);
            if (std::abs(tmp12) < EPS || std::abs(tmp12 - product12) < std::abs(tmp12 + product12))
                if (std::abs(tmp13) < EPS || std::abs(tmp13 - product13) < std::abs(tmp13 + product13))
                    if (std::abs(tmp23) < EPS || std::abs(tmp23 - product23) < std::abs(tmp23 + product23)) {
                        rRotationVector(0) *= pi;
                        rRotationVector(1) *= pi;
                        rRotationVector(2) *= pi;
                        return;
                    }
        }
        { // case 4 --+:
            rRotationVector(0) = -tmp1;
            rRotationVector(1) = -tmp2;
            rRotationVector(2) = tmp3;
            const double tmp12 = rRotationVector(0) * rRotationVector(1);
            const double tmp13 = rRotationVector(0) * rRotationVector(2);
            const double tmp23 = rRotationVector(1) * rRotationVector(2);
            if (std::abs(tmp12) < EPS || std::abs(tmp12 - product12) < std::abs(tmp12 + product12))
                if (std::abs(tmp13) < EPS || std::abs(tmp13 - product13) < std::abs(tmp13 + product13))
                    if (std::abs(tmp23) < EPS || std::abs(tmp23 - product23) < std::abs(tmp23 + product23)) {
                        rRotationVector(0) *= pi;
                        rRotationVector(1) *= pi;
                        rRotationVector(2) *= pi;
                        return;
                    }
        }
        assert(false);
    }

    double tmp = angle / 2.0 / sin(angle);
    rRotationVector(0) = -(rRotationMatrix(1,2) - rRotationMatrix(2,1)) * tmp;
    rRotationVector(1) =  (rRotationMatrix(0,2) - rRotationMatrix(2,0)) * tmp;
    rRotationVector(2) = -(rRotationMatrix(0,1) - rRotationMatrix(1,0)) * tmp;
}

template<class TSparseSpace, class TDenseSpace>
void BeamMapper<TSparseSpace, TDenseSpace>::InitializeOriginForcesAndMoments(const Variable<array_1d<double, 3>>& rOriginVariablesForces,
                                    const Variable<array_1d<double, 3>>& rOriginVariablesMoments)                                    
{   
    for( auto& r_local_sys : mMapperLocalSystems )
    {   
        if( r_local_sys->HasInterfaceInfo())
        {
            auto beam_sys = dynamic_cast<BeamMapperLocalSystem*>(r_local_sys.get());
            KRATOS_ERROR_IF_NOT(beam_sys) << "Expected BeamMapperLocalSystem!" << std::endl;

            MatrixType rotation_matrix_G_B(3, 3);
            VectorType translation_vector_B_P(3);
            VectorType linear_shape_values(2);
            VectorType hermitian_shape_values(4);
            VectorType hermitian_der_shape_values(4);
            GeometryType r_geom; // beam line
            NodePointerType p_node; // surface node

            beam_sys->CalculateRotationMatrixInterfaceInfos(rotation_matrix_G_B,
                                                               translation_vector_B_P,
                                                               linear_shape_values,
                                                               hermitian_shape_values,
                                                               hermitian_der_shape_values,
                                                               r_geom,
                                                               p_node);

            KRATOS_ERROR_IF_NOT(p_node) << "Node is a nullptr"<< std::endl;

            const std::vector<std::string> var_comps{"_X", "_Y", "_Z"};

            IndexType c = 0;
            for (const auto& var_ext : var_comps)
            {
                const auto& var_origin_force = KratosComponents<ComponentVariableType>::Get(rOriginVariablesForces.Name() + var_ext);
                const auto& var_origin_moment = KratosComponents<ComponentVariableType>::Get(rOriginVariablesMoments.Name() + var_ext);

                r_geom[0].FastGetSolutionStepValue(var_origin_force) =  0.0;
                r_geom[0].FastGetSolutionStepValue(var_origin_moment) =  0.0;
                r_geom[1].FastGetSolutionStepValue(var_origin_force) =  0.0;
                r_geom[1].FastGetSolutionStepValue(var_origin_moment) =  0.0; 
                c++;
            }
        }
    }

}

// template instantiation
template class BeamMapper<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>;

}  // namespace Kratos.
