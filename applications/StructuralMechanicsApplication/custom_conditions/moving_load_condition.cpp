// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Aron Noordam
//

// System includes


// External includes


// Project includes
#include "custom_conditions/moving_load_condition.h"
#include "includes/variables.h"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************
template< std::size_t TDim, std::size_t TNumNodes >
MovingLoadCondition< TDim, TNumNodes> ::MovingLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseLoadCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
template< std::size_t TDim, std::size_t TNumNodes >
MovingLoadCondition< TDim, TNumNodes>::MovingLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : BaseLoadCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************
template< std::size_t TDim, std::size_t TNumNodes >
Condition::Pointer MovingLoadCondition< TDim, TNumNodes>::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<MovingLoadCondition>(NewId, pGeom, pProperties);
}

//************************************************************************************
//************************************************************************************
template< std::size_t TDim, std::size_t TNumNodes >
Condition::Pointer MovingLoadCondition< TDim, TNumNodes>::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<MovingLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/
template< std::size_t TDim, std::size_t TNumNodes >
Condition::Pointer MovingLoadCondition< TDim, TNumNodes>::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<MovingLoadCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************
template< std::size_t TDim, std::size_t TNumNodes >
MovingLoadCondition< TDim, TNumNodes>::~MovingLoadCondition()
{
}

//************************************************************************************
//************************************************************************************

template< std::size_t TDim, std::size_t TNumNodes >
bool MovingLoadCondition<TDim, TNumNodes>::HasRotDof() const
{
    KRATOS_TRY
        return GetGeometry()[0].HasDofFor(ROTATION_Z) && GetGeometry().size() == 2;

    KRATOS_CATCH("")
}

template< std::size_t TDim, std::size_t TNumNodes >
void MovingLoadCondition<TDim, TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    const double local_x_coord = this->GetValue(MOVING_LOAD_LOCAL_DISTANCE);

    // check if cond should be calculated
    mIsMovingLoad = false;
    for (IndexType i = 0; i < TDim; ++i) {
        if (std::abs(this->GetValue(POINT_LOAD)[i]) > std::numeric_limits<double>::epsilon() && local_x_coord <= this->GetGeometry().Length() && local_x_coord >= 0.0) {
            mIsMovingLoad = true;
        }
    }
}


template< std::size_t TDim, std::size_t TNumNodes >
void MovingLoadCondition<TDim, TNumNodes>::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    if (mIsMovingLoad) {
        this->CalculateLoadPointDisplacementVector();
        this->CalculateLoadPointRotationVector();
    }
    else {
        this->SetValue(DISPLACEMENT, ZeroVector(3));
        this->SetValue(ROTATION, ZeroVector(3));
    }
}


template< std::size_t TDim, std::size_t TNumNodes >
void MovingLoadCondition< TDim, TNumNodes>::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType block_size = this->GetBlockSize();
    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * block_size;

    if ( CalculateStiffnessMatrixFlag == true ){ //calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size){
            rLeftHandSideMatrix.resize(mat_size, mat_size, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix(mat_size, mat_size); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ){ //calculation of the matrix is required
        if ( rRightHandSideVector.size( ) != mat_size){
            rRightHandSideVector.resize(mat_size, false );
        }

        noalias( rRightHandSideVector ) = ZeroVector(mat_size); //resetting RHS
    }

    // Vector with a loading applied to the condition
    array_1d<double, TDim > moving_load = ZeroVector(TDim);
    if( this->Has(POINT_LOAD) ){
        array_1d<double, 3 > point_load = this->GetValue(POINT_LOAD);
        for (IndexType i = 0; i < TDim; ++i)
        {
            moving_load[i] = point_load[i];
        }
    }

    // apply moving load if moving load is present
    if (mIsMovingLoad){

        const double local_x_coord = this->GetValue(MOVING_LOAD_LOCAL_DISTANCE);
        const GeometryType& r_geom = this->GetGeometry();

        bounded_matrix<double, TDim, TDim> rotation_matrix = ZeroMatrix(TDim, TDim);
        this->CalculateRotationMatrix(rotation_matrix, r_geom);

        // rotate load to local system
        array_1d<double, TDim> local_moving_load;
        noalias(local_moving_load) = prod(rotation_matrix, moving_load);

        VectorType normal_shape_functions_vector;
        VectorType shear_shape_functions_vector;
        VectorType rotational_shape_functions_vector;

        // if element has rotational degrees of freedom, shape functions are exact, thus no interpolation is required
        if (block_size > TDim){
            CalculateExactNormalShapeFunctions(normal_shape_functions_vector, local_x_coord);
            CalculateExactShearShapeFunctions(shear_shape_functions_vector, local_x_coord);
            CalculateExactRotationalShapeFunctions(rotational_shape_functions_vector, local_x_coord);
        } else {
            GeometryType::CoordinatesArrayType local_coordinates_array = ZeroVector(3);
            local_coordinates_array[0] = local_x_coord / r_geom.Length() * 2 - 1;

            r_geom.ShapeFunctionsValues(normal_shape_functions_vector, local_coordinates_array);
            r_geom.ShapeFunctionsValues(shear_shape_functions_vector, local_coordinates_array);

        }

        array_1d<double, TDim> load_first_node = ZeroVector(TDim);
        array_1d<double, TDim> load_end_node = ZeroVector(TDim);

        BoundedMatrix<double, TDim, TNumNodes> local_load_matrix = ZeroMatrix(TDim, TNumNodes);
        BoundedMatrix<double, TDim, TNumNodes> global_load_matrix = ZeroMatrix(TDim, TNumNodes);

        const Matrix global_moment_matrix = CalculateGlobalMomentMatrix(rotational_shape_functions_vector, local_moving_load);

        for (IndexType i_nod = 0; i_nod < TNumNodes; ++i_nod){
            local_load_matrix(0, i_nod) = normal_shape_functions_vector[i_nod] * local_moving_load[0];

            for (IndexType ii = 1; ii < TDim; ++ii){
                local_load_matrix(ii, i_nod) = shear_shape_functions_vector[i_nod] * local_moving_load[ii];
            }
        }

        // rotate load back to global
        noalias(global_load_matrix) = prod(trans(rotation_matrix), local_load_matrix);

        for (IndexType ii = 0; ii < TNumNodes; ++ii){
            const IndexType base = ii * block_size;

            // add load to rhs
            for (IndexType k = 0; k < TDim; ++k){
                const IndexType index = base + k;
                rRightHandSideVector[index] += global_load_matrix(k, ii);
            }

            // shape index is 0 or 1, rotation is only added to first and final node
            for (IndexType k = 0; k < block_size - TDim; ++k){

                const IndexType index = base + TDim + k;
                rRightHandSideVector[index] = global_moment_matrix(k, ii);
            }
        }
    }
    KRATOS_CATCH( "" )
}

template< std::size_t TDim, std::size_t TNumNodes >
Matrix MovingLoadCondition<TDim, TNumNodes>::CalculateGlobalMomentMatrix(const VectorType& rRotationalShapeFunctionVector, const array_1d<double, TDim>& rLocalMovingLoad) const
{
    KRATOS_TRY

    // check if condition has rotation dof 
    const bool has_rot_dof = this->HasRotDof();

    Matrix global_moment_matrix;

    if constexpr (TDim == 2) {
        global_moment_matrix.resize(1, TNumNodes, false);
    } else if constexpr (TDim == 3) {
        global_moment_matrix.resize(3, TNumNodes, false);
    } 

    if (has_rot_dof) {
        if constexpr (TDim == 2) {
            // rotation around z axis (2D)
            global_moment_matrix(0, 0) = rRotationalShapeFunctionVector[0] * rLocalMovingLoad[1];
            global_moment_matrix(0, 1) = rRotationalShapeFunctionVector[1] * rLocalMovingLoad[1];
        } else if constexpr (TDim == 3){
            // rotation around y and z axis (3D)
            global_moment_matrix(0, 0) = 0;
            global_moment_matrix(1, 0) = rRotationalShapeFunctionVector[0] * rLocalMovingLoad[2];
            global_moment_matrix(2, 0) = rRotationalShapeFunctionVector[0] * rLocalMovingLoad[1];

            global_moment_matrix(0, 1) = 0;
            global_moment_matrix(1, 1) = rRotationalShapeFunctionVector[1] * rLocalMovingLoad[2];
            global_moment_matrix(2, 1) = rRotationalShapeFunctionVector[1] * rLocalMovingLoad[1];
        }
    }
    return global_moment_matrix;
    KRATOS_CATCH("")

}

template< std::size_t TDim, std::size_t TNumNodes >
void MovingLoadCondition< TDim, TNumNodes>::CalculateExactNormalShapeFunctions(VectorType& rShapeFunctionsVector, const double LocalXCoord) const
{
    KRATOS_TRY
    if (rShapeFunctionsVector.size() != TNumNodes) {
        rShapeFunctionsVector.resize(TNumNodes, false);
    }
    const auto& r_geometry = this->GetGeometry();
    const double length = r_geometry.Length();

    rShapeFunctionsVector[0] = 1 - LocalXCoord / length;
    rShapeFunctionsVector[1] = LocalXCoord / length;

    KRATOS_CATCH("")

}


template< std::size_t TDim, std::size_t TNumNodes >
void MovingLoadCondition< TDim, TNumNodes>::CalculateExactShearShapeFunctions(VectorType& rShapeFunctionsVector, const double LocalXCoord) const
{
    KRATOS_TRY
    if (rShapeFunctionsVector.size() != TNumNodes) {
        rShapeFunctionsVector.resize(TNumNodes, false);
    }
    const auto& r_geometry = this->GetGeometry();
    const double length = r_geometry.Length();

    rShapeFunctionsVector[0] = 1 + 2 * std::pow((LocalXCoord / length), 3) - 3 * std::pow((LocalXCoord / length), 2);
    rShapeFunctionsVector[1] = -2 * std::pow((LocalXCoord / length), 3) + 3 * std::pow((LocalXCoord / length), 2);
    KRATOS_CATCH("")
}


template< std::size_t TDim, std::size_t TNumNodes >
void MovingLoadCondition< TDim, TNumNodes>::CalculateExactRotationalShapeFunctions(VectorType& rShapeFunctionsVector, const double LocalXCoord) const
{
    KRATOS_TRY
    if (rShapeFunctionsVector.size() != TNumNodes) {
        rShapeFunctionsVector.resize(TNumNodes, false);
    }
    const auto& r_geometry = this->GetGeometry();
    const double length = r_geometry.Length();

    rShapeFunctionsVector[0] = LocalXCoord + (std::pow(LocalXCoord, 3) / std::pow(length, 2)) - 2 * (std::pow(LocalXCoord, 2) / length);
    rShapeFunctionsVector[1] = (std::pow(LocalXCoord, 3) / std::pow(length, 2)) - (std::pow(LocalXCoord, 2) / length);
    KRATOS_CATCH("")
}

template< std::size_t TDim, std::size_t TNumNodes >
void MovingLoadCondition< TDim, TNumNodes>::CalculateExactShearShapeFunctionsDerivatives(VectorType& rShapeFunctionsDerivativesVector, const double LocalXCoord) const
{
    KRATOS_TRY
    if (rShapeFunctionsDerivativesVector.size() != TNumNodes) {
        rShapeFunctionsDerivativesVector.resize(TNumNodes, false);
    }

    const auto& r_geometry = this->GetGeometry();
    const double length = r_geometry.Length();

    rShapeFunctionsDerivativesVector[0] = 6 * std::pow(LocalXCoord, 2) / std::pow(length, 3) - 6 * LocalXCoord / std::pow(length, 2);
    rShapeFunctionsDerivativesVector[1] = -6 * std::pow(LocalXCoord, 2)/ std::pow(length, 3) + 6 * LocalXCoord / std::pow(length, 2);

    KRATOS_CATCH("")
}


template< std::size_t TDim, std::size_t TNumNodes >
void MovingLoadCondition< TDim, TNumNodes>::CalculateExactRotationalShapeFunctionsDerivatives(VectorType& rShapeFunctionsDerivativesVector, const double LocalXCoord) const
{
    KRATOS_TRY

    if (rShapeFunctionsDerivativesVector.size() != TNumNodes) {
        rShapeFunctionsDerivativesVector.resize(TNumNodes, false);
    }

    const auto& r_geometry = this->GetGeometry();
    const double length = r_geometry.Length();

    rShapeFunctionsDerivativesVector[0] = 1 + 3 * std::pow(LocalXCoord/length, 2)  - 4 * LocalXCoord / length;
    rShapeFunctionsDerivativesVector[1] = 3 * std::pow(LocalXCoord / length, 2) - 2 * LocalXCoord / length;

    KRATOS_CATCH("")
}

template< std::size_t TDim, std::size_t TNumNodes >
void MovingLoadCondition< TDim, TNumNodes>::CalculateRotationMatrix(BoundedMatrix<double, TDim, TDim>& rRotationMatrix, const GeometryType& rGeom)
{
    KRATOS_TRY

    constexpr double tolerance = 1e-8;
    //Unitary vector in local x direction
    array_1d<double, 3> vx = ZeroVector(3);
    noalias(vx) = rGeom.GetPoint(1) - rGeom.GetPoint(0);
    const double inv_norm_x = 1.0 / norm_2(vx);

    for (IndexType i = 0; i < TDim; ++i) {
        vx[i] *= inv_norm_x;
    }

    array_1d<double, 3> vy = ZeroVector(3);

    array_1d<double, 3> vy_tmp = ZeroVector(3);
    array_1d<double, 3> vz_tmp = ZeroVector(3);

    vy_tmp[1] = 1;
    vz_tmp[2] = 1;

    // Unitary vector in local y direction
    if (std::fabs(vx[0]) < tolerance && std::fabs(vx[1]) < tolerance) {
        MathUtils<double>::CrossProduct(vy, vy_tmp, vx);
    } else {
        MathUtils<double>::CrossProduct(vy, vz_tmp, vx);
    }

    // set 3d part of Rotation Matrix
    if constexpr (TDim == 3) {
        array_1d<double, 3> vz = ZeroVector(3);

        // Unitary vector in local y direction
        const double inv_norm_y = 1.0 / norm_2(vy);
        vy[0] *= inv_norm_y;
        vy[1] *= inv_norm_y;
        vy[2] *= inv_norm_y;

        // Unitary vector in local z direction
        MathUtils<double>::CrossProduct(vz, vx, vy);
        const double inv_norm_z = 1.0 / norm_2(vz);
        if (inv_norm_z > tolerance) {
            vz[0] *= inv_norm_z;
            vz[1] *= inv_norm_z;
            vz[2] *= inv_norm_z;
        }
   
        rRotationMatrix(0, 2) = vx[2];

        rRotationMatrix(1, 2) = vy[2];

        rRotationMatrix(2, 0) = vz[0];
        rRotationMatrix(2, 1) = vz[1];
        rRotationMatrix(2, 2) = vz[2];
    }

    // add 2d part of Rotation Matrix
    rRotationMatrix(0, 0) = vx[0];
    rRotationMatrix(0, 1) = vx[1];

    rRotationMatrix(1, 0) = vy[0];
    rRotationMatrix(1, 1) = vy[1];

    KRATOS_CATCH("")
}
    

template< std::size_t TDim, std::size_t TNumNodes >
Vector MovingLoadCondition< TDim, TNumNodes>::CalculateLoadPointDisplacementVector()
{

    KRATOS_TRY

    // Get global displacement vector
    Vector displacement_vector;
    this->GetValuesVector(displacement_vector);

    // check if rotation degrees of freedom are active
    const bool has_rot_dof = this->HasRotDof();

    // Convert displacement vector to a ndim by nnodes matrix
    bounded_matrix<double, TDim, TNumNodes> displacement_matrix = ZeroMatrix(TDim, TNumNodes);

    IndexType vector_index = 0;
    for (IndexType ii = 0; ii < TNumNodes; ++ii) {
        for (IndexType jj = 0; jj < TDim; ++jj) {
            displacement_matrix(jj, ii) = displacement_vector[vector_index];
            vector_index++;
        }
    }


    // Get global nodal rotation matrix
    bounded_matrix<double, 3, TNumNodes> nodal_rotation_matrix = ZeroMatrix(3, TNumNodes);
    if (has_rot_dof) {
        for (IndexType ii = 0; ii < TNumNodes; ++ii) {

            nodal_rotation_matrix(0, ii) = GetGeometry()[ii].FastGetSolutionStepValue(ROTATION_X);
            nodal_rotation_matrix(1, ii) = GetGeometry()[ii].FastGetSolutionStepValue(ROTATION_Y);
            nodal_rotation_matrix(2, ii) = GetGeometry()[ii].FastGetSolutionStepValue(ROTATION_Z);
        }
    }

    // calculate elemental rotation matrix
    bounded_matrix<double, TDim, TDim> rotation_matrix = ZeroMatrix(TDim, TDim);
    this->CalculateRotationMatrix(rotation_matrix, this->GetGeometry());

    // make sure the elemental rotation matrix is 3x3
    bounded_matrix<double, 3, 3> full_rotation_matrix = ZeroMatrix(3, 3);
    for (IndexType i = 0; i < TDim; i++) {
        for (IndexType j = 0; j < TDim; j++) {
            full_rotation_matrix(i, j) = rotation_matrix(i, j);
        }
    }
    if constexpr (TDim == 2) {
        full_rotation_matrix(2, 2) = 1;
    }

    // Calculate local nodal rotation matrix
    bounded_matrix<double, 3, TNumNodes> local_nodal_rotation_matrix = prod(full_rotation_matrix, nodal_rotation_matrix);

    // calculate local displacement matrix
    bounded_matrix<double, TDim, TNumNodes> local_disp_matrix = prod(rotation_matrix, displacement_matrix);

    // Get local coordinate of the moving load position
    const double local_x_coord = this->GetValue(MOVING_LOAD_LOCAL_DISTANCE);

    // calculate shape functions
    VectorType normal_shape_functions_vector;
    VectorType shear_shape_functions_vector;
    VectorType rotational_shape_functions_vector;

    auto& r_geom = this->GetGeometry();

    // if rotation DOFS are active, shape functions are the same as in euler beams, else use the condition element shape functions
    if (has_rot_dof) {

        this->CalculateExactNormalShapeFunctions(normal_shape_functions_vector, local_x_coord);
        this->CalculateExactShearShapeFunctions(shear_shape_functions_vector, local_x_coord);
        this->CalculateExactRotationalShapeFunctions(rotational_shape_functions_vector, local_x_coord);
    } else {
        GeometryType::CoordinatesArrayType local_coordinates_array = ZeroVector(3);
        local_coordinates_array[0] = local_x_coord / r_geom.Length() * 2 - 1;

        r_geom.ShapeFunctionsValues(normal_shape_functions_vector, local_coordinates_array);
        r_geom.ShapeFunctionsValues(shear_shape_functions_vector, local_coordinates_array);
    }

    // calculate local displacement vector at the location of the moving load
    VectorType local_disp_vector = ZeroVector(TDim);
    if constexpr (TDim == 2) {

        // calculate inner product local displacement and shape functions
        double disp_shear_axis = 0;
        double disp_normal_axis = 0;
        for (IndexType ii = 0; ii < TNumNodes; ++ii) {

            disp_normal_axis += local_disp_matrix(0, ii) * normal_shape_functions_vector[ii];
            disp_shear_axis += local_disp_matrix(1, ii) * shear_shape_functions_vector[ii];

            if (has_rot_dof) {
                disp_shear_axis += local_nodal_rotation_matrix(2, ii) * rotational_shape_functions_vector[ii];
            }
        }

        local_disp_vector(0) = disp_normal_axis;
        local_disp_vector(1) = disp_shear_axis;

    } else if constexpr (TDim == 3) {
        double disp_normal_axis = 0;
        double disp_shear_axis_1 = 0;
        double disp_shear_axis_2 = 0;


        // calculate inner product local displacement and shape functions
        for (IndexType ii = 0; ii < TNumNodes; ++ii) {

            disp_normal_axis += local_disp_matrix(0, ii) * normal_shape_functions_vector[ii];
            disp_shear_axis_1 += local_disp_matrix(1, ii) * shear_shape_functions_vector[ii];
            disp_shear_axis_2 += local_disp_matrix(2, ii) * shear_shape_functions_vector[ii];

            if (has_rot_dof) {
                disp_shear_axis_1 += local_nodal_rotation_matrix(2, ii) * rotational_shape_functions_vector[ii];
                disp_shear_axis_2 += local_nodal_rotation_matrix(1, ii) * rotational_shape_functions_vector[ii];
            }
        }

        local_disp_vector(0) = disp_normal_axis;
        local_disp_vector(1) = disp_shear_axis_1;
        local_disp_vector(2) = disp_shear_axis_2;
    }
    
    // calculate global displacement vector at the location of the moving load
    VectorType global_point_disp_vector = prod(trans(rotation_matrix), local_disp_vector);

    // make sure displacement is a vector with size 3
    Vector displacements = ZeroVector(3);
    for (IndexType ii = 0; ii < TDim; ++ii) {
        displacements(ii) = global_point_disp_vector(ii);
    }

    // Set Displacement at the location of the point load to the element
    this->SetValue(DISPLACEMENT, displacements);

    return displacements;

    KRATOS_CATCH("")
}

template< std::size_t TDim, std::size_t TNumNodes >
Vector MovingLoadCondition< TDim, TNumNodes>::CalculateLoadPointRotationVector()
{


    KRATOS_TRY

    // Get global displacement vector
    Vector displacement_vector;
    this->GetValuesVector(displacement_vector);

    // check if rotation degrees of freedom are active
    const bool has_rot_dof = this->HasRotDof();

    // Convert displacement vector to a ndim by nnodes matrix
    bounded_matrix<double, TDim, TNumNodes> displacement_matrix = ZeroMatrix(TDim, TNumNodes);

    IndexType vector_index = 0;
    for (IndexType ii = 0; ii < TNumNodes; ++ii) {

        for (IndexType jj = 0; jj < TDim; ++jj) {
            displacement_matrix(jj, ii) = displacement_vector[vector_index];
            vector_index++;
        }
    }

    // Get global nodal rotation matrix
    bounded_matrix<double, 3, TNumNodes> nodal_rotation_matrix = ZeroMatrix(3, TNumNodes);

    if (has_rot_dof) {
        for (IndexType ii = 0; ii < TNumNodes; ++ii) {

            nodal_rotation_matrix(0, ii) = GetGeometry()[ii].FastGetSolutionStepValue(ROTATION_X);
            nodal_rotation_matrix(1, ii) = GetGeometry()[ii].FastGetSolutionStepValue(ROTATION_Y);
            nodal_rotation_matrix(2, ii) = GetGeometry()[ii].FastGetSolutionStepValue(ROTATION_Z);
        }
    }


    // calculate elemental rotation matrix
    bounded_matrix<double, TDim, TDim> rotation_matrix = ZeroMatrix(TDim, TDim);
    this->CalculateRotationMatrix(rotation_matrix, this->GetGeometry());


    // define full rotation matrix which is required to back calculate the global rotation at the location of the moving load
    bounded_matrix<double, 3, 3> full_rotation_matrix = ZeroMatrix(3, 3);

    for (IndexType i = 0; i<TDim; ++i) {
        for (IndexType j = 0; j < TDim; ++j) {
            full_rotation_matrix(i, j) = rotation_matrix(i, j);
        }
    }
    if constexpr (TDim == 2) {
        full_rotation_matrix(2, 2) = 1;
    }


    // calculate local nodal rotation matrix
    bounded_matrix<double, 3, TNumNodes> local_nodal_rotation_matrix = prod(full_rotation_matrix, nodal_rotation_matrix);

    // calculate local displacement matrix
    bounded_matrix<double, TDim, TNumNodes> local_disp_matrix = prod(rotation_matrix, displacement_matrix);


    // get local location of the moving load
    const double local_x_coord = this->GetValue(MOVING_LOAD_LOCAL_DISTANCE);

    // calculate shape function derivatives vector
    VectorType shear_shape_functions_derivatives_vector = ZeroVector(TNumNodes);
    VectorType rotational_shape_functions_derivatives_vector = ZeroVector(TNumNodes);

    auto& r_geom = this->GetGeometry();

    // if rotation DOFS are active, shape functions derivatives are the same as in euler beams, else use the condition element shape function derivatives
    if (has_rot_dof) {

        this->CalculateExactShearShapeFunctionsDerivatives(shear_shape_functions_derivatives_vector, local_x_coord);
        this->CalculateExactRotationalShapeFunctionsDerivatives(rotational_shape_functions_derivatives_vector, local_x_coord);
    } else {
        GeometryType::CoordinatesArrayType local_coordinates_array = ZeroVector(3);
        local_coordinates_array[0] = local_x_coord / r_geom.Length() * 2 - 1;

        Matrix local_gradient_matrix;

        r_geom.ShapeFunctionsLocalGradients(local_gradient_matrix, local_coordinates_array);
        for (IndexType i = 0; i < TNumNodes; ++i) {
            shear_shape_functions_derivatives_vector(i) = local_gradient_matrix(i, 0);
        }
    }

    // calculate local rotation at the location of the moving load
    VectorType local_rot_vector = ZeroVector(3);
    if constexpr (TDim == 2) {
        
        double local_rotation = 0;
        // calculate inner product local displacement and shape functions
        for (IndexType ii = 0; ii < TNumNodes; ++ii) {

            local_rotation += local_disp_matrix(1, ii) * shear_shape_functions_derivatives_vector[ii];

            if (has_rot_dof) {
                local_rotation += local_nodal_rotation_matrix(2, ii) * rotational_shape_functions_derivatives_vector[ii];
            }
        }

        local_rot_vector(2) = local_rotation;
    } else if constexpr (TDim == 3) {

        double local_rotation_axis_1 = 0;
        double local_rotation_axis_2 = 0;

        // calculate inner product local displacement and shape functions
        for (IndexType ii = 0; ii < TNumNodes; ++ii) {

            local_rotation_axis_1 += local_disp_matrix(2, ii) * shear_shape_functions_derivatives_vector[ii];
            local_rotation_axis_2 += local_disp_matrix(1, ii) * shear_shape_functions_derivatives_vector[ii];

            if (has_rot_dof) {
                local_rotation_axis_1 += local_nodal_rotation_matrix(1, ii) * rotational_shape_functions_derivatives_vector[ii];
                local_rotation_axis_2 += local_nodal_rotation_matrix(2, ii) * rotational_shape_functions_derivatives_vector[ii];
            }
        }

        local_rot_vector(0) = 0;
        local_rot_vector(1) = local_rotation_axis_1;
        local_rot_vector(2) = local_rotation_axis_2;
    }

    // calculate global rotation vector at the location of the moving load
    VectorType global_point_rotation_vector = ZeroVector(3);
    if constexpr (TDim == 2) {
        global_point_rotation_vector(2) = local_rot_vector(2);
    } else if constexpr (TDim == 3) {
        global_point_rotation_vector = prod(trans(full_rotation_matrix), local_rot_vector);
    }

    // Set Displacement at the location of the point load to the element
    this->SetValue(ROTATION, global_point_rotation_vector);

    return global_point_rotation_vector;
    KRATOS_CATCH("")
}


template <std::size_t TDim, std::size_t TNumNodes>
void MovingLoadCondition<TDim, TNumNodes>::GetRotationsVector(Vector& rRotationsVector, const int Step) const
{

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dim = GetGeometry().WorkingSpaceDimension();

    SizeType mat_size;
    if constexpr (TDim == 2) {
        mat_size = number_of_nodes;
    } else {
        mat_size = number_of_nodes * dim;
    }


    if (rRotationsVector.size() != mat_size) {
        rRotationsVector.resize(mat_size, false);
    }

    if constexpr (TDim == 2) {
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            rRotationsVector[i] = GetGeometry()[i].FastGetSolutionStepValue(ROTATION_Z, Step);
        }
    } else {
        for (SizeType i = 0; i < number_of_nodes; ++i) {
            const array_1d<double, 3 >& r_rotation = GetGeometry()[i].FastGetSolutionStepValue(ROTATION, Step);
            const SizeType index = i * dim;
            for (SizeType k = 0; k < dim; ++k) {
                rRotationsVector[index + k] = r_rotation[k];
            }
        }
    }
}


template class MovingLoadCondition<2, 2>;
template class MovingLoadCondition<2, 3>;
template class MovingLoadCondition<3, 2>;
template class MovingLoadCondition<3, 3>;

} // Namespace Kratos


