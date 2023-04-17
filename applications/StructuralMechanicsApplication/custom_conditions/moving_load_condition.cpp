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

    const double local_x_coord = this->GetValue(MOVING_LOAD_LOCAL_DISTANCE);

    // check if cond should be calculated
    bool is_moving_load = false;
    for (IndexType i = 0; i < TDim; ++i){
        if (std::abs(moving_load[i]) > std::numeric_limits<double>::epsilon() && local_x_coord <= this->GetGeometry().Length() && local_x_coord >= 0.0){
            is_moving_load = true;
        }
    }

    // apply moving load if moving load is present
    if (is_moving_load){
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
Matrix MovingLoadCondition<TDim, TNumNodes>::CalculateGlobalMomentMatrix(const VectorType& rRotationalShapeFunctionVector, const array_1d<double, TDim> LocalMovingLoad) const
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
            global_moment_matrix(0, 0) = rRotationalShapeFunctionVector[0] * LocalMovingLoad[1];
            global_moment_matrix(0, 1) = rRotationalShapeFunctionVector[1] * LocalMovingLoad[1];
        } else if constexpr (TDim == 3){
            // rotation around y and z axis (3D)
            global_moment_matrix(0, 0) = 0;
            global_moment_matrix(1, 0) = rRotationalShapeFunctionVector[0] * LocalMovingLoad[2];
            global_moment_matrix(2, 0) = rRotationalShapeFunctionVector[0] * LocalMovingLoad[1];

            global_moment_matrix(0, 1) = 0;
            global_moment_matrix(1, 1) = rRotationalShapeFunctionVector[1] * LocalMovingLoad[2];
            global_moment_matrix(2, 1) = rRotationalShapeFunctionVector[1] * LocalMovingLoad[1];
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
    }
    else {
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
bool MovingLoadCondition<TDim, TNumNodes>::HasRotDof() const
{
    KRATOS_TRY
    return GetGeometry()[0].HasDofFor(ROTATION_Z) && GetGeometry().size() == 2;

    KRATOS_CATCH("")
}

template class MovingLoadCondition<2, 2>;
template class MovingLoadCondition<2, 3>;
template class MovingLoadCondition<3, 2>;
template class MovingLoadCondition<3, 3>;

} // Namespace Kratos


