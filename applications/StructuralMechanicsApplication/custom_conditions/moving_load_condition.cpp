// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

// System includes


// External includes


// Project includes
#include "custom_conditions/moving_load_condition.h"

#include <includes/element.h>

#include "includes/checks.h"

#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"
#include "utilities/beam_math_utilities.hpp"

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

    const unsigned int NumberOfNodes = GetGeometry().size();
    const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();

    const unsigned int block_size = this->GetBlockSize();
    // Resizing as needed the LHS
    unsigned int MatSize = NumberOfNodes * block_size;
   /* for (const auto& r_node : GetGeometry().Points()) {
        if (r_node.HasDofFor(ROTATION_X)) { MatSize += 1; }
        if (r_node.HasDofFor(ROTATION_Y)) { MatSize += 1; }
        if (r_node.HasDofFor(ROTATION_Z)) { MatSize += 1; }
    }*/


    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
        {
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size( ) != MatSize )
        {
            rRightHandSideVector.resize( MatSize, false );
        }

        noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
    }

    // Vector with a loading applied to the condition
    array_1d<double, TDim > MovingLoad = ZeroVector(TDim);
    if( this->Has(POINT_LOAD) )
    {
        noalias(MovingLoad) = this->GetValue(POINT_LOAD);
    }

    double local_x_coord = this->GetValue(MOVING_LOAD_LOCAL_DISTANCE);

	// check if cond should be calculated
    bool is_moving_load = false;
    for (int i = 0; i < TDim; ++i)
    {
        if (std::abs(MovingLoad[i]) > DBL_EPSILON && local_x_coord <=this->GetGeometry().Length() && local_x_coord >= 0.0)
        {
            is_moving_load = true;
        }
    }

    // apply moving load if moving load is present
    if (is_moving_load)
    {
        GeometryType& rGeom = this->GetGeometry();

        bounded_matrix<double, TDim, TDim> rotation_matrix= ZeroMatrix(TDim, TDim);
        CalculateRotationMatrix(rotation_matrix, rGeom);

        // rotate load to local system
        array_1d<double, TDim> localMovingLoad;
        noalias(localMovingLoad) = prod(rotation_matrix, MovingLoad);
       

        VectorType normalShapeFunctionsVector;
        VectorType shearShapeFunctionsVector;
        VectorType rotationalShapeFunctionsVector;



        // if element has rotational degrees of freedom, shape functions are exact, thus no interpolation is required
        if (block_size > TDim)
        {
            CalculateExactNormalShapeFunctions(normalShapeFunctionsVector, local_x_coord);
            CalculateExactShearShapeFunctions(shearShapeFunctionsVector, local_x_coord);
            CalculateExactRotationalShapeFunctions(rotationalShapeFunctionsVector, local_x_coord);
        }
        else
        {
            GeometryType::CoordinatesArrayType local_coordinates_array = ZeroVector(3);
            local_coordinates_array[0] = local_x_coord / rGeom.Length() * 2 - 1;

            rGeom.ShapeFunctionsValues(normalShapeFunctionsVector, local_coordinates_array);
            rGeom.ShapeFunctionsValues(shearShapeFunctionsVector, local_coordinates_array);

        }

        array_1d<double, TDim> load_first_node = ZeroVector(TDim);
        array_1d<double, TDim> load_end_node = ZeroVector(TDim);

        BoundedMatrix<double, TDim, TNumNodes> local_load_matrix = ZeroMatrix(TDim, TNumNodes);
        BoundedMatrix<double, TDim, TNumNodes> global_load_matrix = ZeroMatrix(TDim, TNumNodes);

        std::vector<double>  moment_first_node;
        std::vector<double>  moment_end_node;

        moment_first_node.resize(block_size - TDim, false);
        moment_end_node.resize(block_size - TDim, false);


        for (unsigned int nod = 0; nod < TNumNodes; ++nod)
        {
            local_load_matrix(0, nod) = normalShapeFunctionsVector[nod] * localMovingLoad[0];

            for (unsigned int ii = 1; ii < TDim; ++ii)
            {
                local_load_matrix(ii,nod) = shearShapeFunctionsVector[nod] * localMovingLoad[ii];
            }
        }

        // rotate load back to global
        noalias(global_load_matrix) = prod(trans(rotation_matrix), local_load_matrix);


        // rotation around y and z axis (3D)
        if (block_size > TDim + 1)
        {
            moment_first_node[0] = 0;
            moment_first_node[1] = rotationalShapeFunctionsVector[0] * localMovingLoad[2];
            moment_first_node[2] = rotationalShapeFunctionsVector[0] * localMovingLoad[1];

            moment_end_node[0] = 0;
            moment_end_node[1] = rotationalShapeFunctionsVector[1] * localMovingLoad[2];
            moment_end_node[2] = rotationalShapeFunctionsVector[1] * localMovingLoad[1];

        }
        // rotation around z axis (2D)
        else if (block_size > TDim)
        {
            moment_first_node[0] = rotationalShapeFunctionsVector[0] * localMovingLoad[1];
            moment_end_node[0] = rotationalShapeFunctionsVector[0] * localMovingLoad[1];
        }


        array_1d<int, 2> shape_indices;
        shape_indices[0] = 0;
        shape_indices[1] = NumberOfNodes - 1;

        for (unsigned int ii = 0; ii < TNumNodes; ++ii)
        //for (unsigned int ii = 0; ii < shape_indices.size(); ++ii)
        {
            const unsigned int base = shape_indices[ii] * block_size;

            // only add load and rotation to RHS if current node is first or final node
            for (unsigned int k = 0; k < TDim; ++k)
            {
                //const double load = global_load_matrix(k, ii);
                //const double load = (shape_indices[ii] == 0) ? load_first_node[k] : load_end_node[k];
                rRightHandSideVector[base + k] += global_load_matrix(k, ii);
            }


            // shape index is 0 or 1, rotation is only added to first and final node
            for (unsigned int k = 0; k < block_size - TDim; ++k)
            {
                //const double load = (shape_indices[ii] == 0) ? load_first_node[k] : load_end_node[ii];
                //double rotation = rotationalShapeFunctionsVector[shape_indices[ii]] * MovingLoad[0];

                const double moment = (shape_indices[ii] == 0) ? moment_first_node[k] : moment_end_node[k];

                rRightHandSideVector[base + TDim +k] = moment;
            }
        }
    }
    KRATOS_CATCH( "" )
}

template< std::size_t TDim, std::size_t TNumNodes >
void MovingLoadCondition< TDim, TNumNodes>::CalculateExactNormalShapeFunctions(VectorType& rShapeFunctionsVector, double local_x_coord) const
{
    if (rShapeFunctionsVector.size() != TNumNodes) {
        rShapeFunctionsVector.resize(TNumNodes, false);
    }
    auto& rGeom = this->GetGeometry();
    const double length = rGeom.Length();

    rShapeFunctionsVector[0] = 1 - local_x_coord / length;
    rShapeFunctionsVector[1] = local_x_coord / length;

    if (TNumNodes==3)
    {
        rShapeFunctionsVector[2] = 0;
    }

    //rResult[0] = 0.5 * (1.0 - rCoordinates[0]);
    //rResult[1] = 0.5 * (1.0 + rCoordinates[0]);
    //rShapeFunctionsVector[0] = 1 + 2 * std::pow((local_x_coord / length), 3) - 3 * std::pow((local_x_coord / length), 2);
    //rShapeFunctionsVector[1] = -2 * std::pow((local_x_coord / length), 3) + 3 * std::pow((local_x_coord / length), 2);
}


template< std::size_t TDim, std::size_t TNumNodes >
void MovingLoadCondition< TDim, TNumNodes>::CalculateExactShearShapeFunctions(VectorType& rShapeFunctionsVector, double local_x_coord) const
{
    if (rShapeFunctionsVector.size() != TNumNodes) {
        rShapeFunctionsVector.resize(TNumNodes, false);
    }
    auto& rGeom = this->GetGeometry();
    const double length = rGeom.Length();

    rShapeFunctionsVector[0] = 1 + 2 * std::pow((local_x_coord / length), 3) - 3 * std::pow((local_x_coord / length), 2);
    rShapeFunctionsVector[1] = -2 * std::pow((local_x_coord / length), 3) + 3 * std::pow((local_x_coord / length), 2);

    if (TNumNodes == 3)
    {
        rShapeFunctionsVector[2] = 0;
    }


}

template< std::size_t TDim, std::size_t TNumNodes >
void MovingLoadCondition< TDim, TNumNodes>::CalculateExactRotationalShapeFunctions(VectorType& rShapeFunctionsVector, double local_x_coord) const
{
    if (rShapeFunctionsVector.size() != TNumNodes) {
        rShapeFunctionsVector.resize(TNumNodes, false);
    }
    auto& rGeom = this->GetGeometry();
    double length = rGeom.Length();

    rShapeFunctionsVector[0] = local_x_coord + (pow(local_x_coord, 3) / pow(length, 2)) - 2 * (pow(local_x_coord, 2) / length);
    rShapeFunctionsVector[1] = (pow(local_x_coord, 3) / pow(length, 2)) - (pow(local_x_coord, 2) / length);

    if (TNumNodes == 3)
    {
        rShapeFunctionsVector[2] = 0;
    }
}

//************************************************************************************
//************************************************************************************

template<>
void MovingLoadCondition<2, 3>::CalculateRotationMatrix(BoundedMatrix<double, 2, 2>& rRotationMatrix, const Element::GeometryType& rGeom)
{
    constexpr double tolerance = 1e-8;

    //Unitary vector in local x direction
    array_1d<double, 3> Vx = ZeroVector(3);
    noalias(Vx) = rGeom.GetPoint(1) - rGeom.GetPoint(0);
    const double inv_norm_x = 1.0 / norm_2(Vx);
    Vx[0] *= inv_norm_x;
    Vx[1] *= inv_norm_x;

    rRotationMatrix(0, 0) = Vx[0];
    rRotationMatrix(0, 1) = Vx[1];

    array_1d<double, 3> Vy = ZeroVector(3);
    array_1d<double, 3> VyTmp = ZeroVector(3);
    array_1d<double, 3> VzTmp = ZeroVector(3);

    VyTmp[1] = 1;
    VzTmp[2] = 1;

    // Unitary vector in local y direction
    if (fabs(Vx[0]) < tolerance && fabs(Vx[1]) < tolerance) {
        MathUtils<double>::CrossProduct(Vy, VyTmp, Vx);
    }
    else {
        MathUtils<double>::CrossProduct(Vy, VzTmp, Vx);
    }

    //Rotation Matrix
    rRotationMatrix(0, 0) = Vx[0];
    rRotationMatrix(0, 1) = Vx[1];

    rRotationMatrix(1, 0) = Vy[0];
    rRotationMatrix(1, 1) = Vy[1];
}

template<>
/// <summary>
/// Rotation matrix of a line3D3 element
/// </summary>
/// <param name="rRotationMatrix"></param>
/// <param name="rGeom"></param>
void MovingLoadCondition<3, 3>::CalculateRotationMatrix(BoundedMatrix<double, 3, 3>& rRotationMatrix, const Element::GeometryType& rGeom)
{
    constexpr double tolerance = 1e-8;

    //Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx) = rGeom.GetPoint(1) - rGeom.GetPoint(0);
    const double inv_norm_x = 1.0 / norm_2(Vx);
    Vx[0] *= inv_norm_x;
    Vx[1] *= inv_norm_x;
    Vx[2] *= inv_norm_x;

    //initialise vectors in local y and z direction
    array_1d<double, 3> Vy = ZeroVector(3);
    array_1d<double, 3> Vz = ZeroVector(3);

    array_1d<double, 3> VyTmp = ZeroVector(3);
    array_1d<double, 3> VzTmp = ZeroVector(3);

    VyTmp[1] = 1;
    VzTmp[2] = 1;

    // Unitary vector in local y direction
    if (fabs(Vx[0]) < tolerance && fabs(Vx[1]) < tolerance) {
        MathUtils<double>::CrossProduct(Vy, VyTmp, Vx);
    }
    else {
        MathUtils<double>::CrossProduct(Vy, VzTmp, Vx);
    }

    double inv_norm_y = 1.0 / norm_2(Vy);
    Vy[0] *= inv_norm_y;
    Vy[1] *= inv_norm_y;
    Vy[2] *= inv_norm_y;

    // Unitary vector in local z direction
    MathUtils<double>::CrossProduct(Vz, Vx, Vy);
    const double inv_norm_z = 1.0 / norm_2(Vz);
    if (inv_norm_z > tolerance)
    {
        Vz[0] *= inv_norm_z;
        Vz[1] *= inv_norm_z;
        Vz[2] *= inv_norm_z;
    }

    //Rotation Matrix
    rRotationMatrix(0, 0) = Vx[0];
    rRotationMatrix(0, 1) = Vx[1];
    rRotationMatrix(0, 2) = Vx[2];

    rRotationMatrix(1, 0) = Vy[0];
    rRotationMatrix(1, 1) = Vy[1];
    rRotationMatrix(1, 2) = Vy[2];

    rRotationMatrix(2, 0) = Vz[0];
    rRotationMatrix(2, 1) = Vz[1];
    rRotationMatrix(2, 2) = Vz[2];
}

template<>
/// <summary>
/// Rotation matrix of a line2D2 element
/// </summary>
/// <param name="rRotationMatrix"></param>
/// <param name="rGeom"></param>
void MovingLoadCondition<2, 2>::CalculateRotationMatrix(BoundedMatrix<double, 2, 2>& rRotationMatrix, const Element::GeometryType& rGeom)
{
    constexpr double tolerance = 1e-8;

    //Unitary vector in local x direction
    array_1d<double, 3> Vx = ZeroVector(3);
    noalias(Vx) = rGeom.GetPoint(1) - rGeom.GetPoint(0);
    const double inv_norm_x = 1.0 / norm_2(Vx);
    Vx[0] *= inv_norm_x;
    Vx[1] *= inv_norm_x;

    rRotationMatrix(0, 0) = Vx[0];
    rRotationMatrix(0, 1) = Vx[1];

    array_1d<double, 3> Vy = ZeroVector(3);
    array_1d<double, 3> VyTmp = ZeroVector(3);
    array_1d<double, 3> VzTmp = ZeroVector(3);

    VyTmp[1] = 1;
    VzTmp[2] = 1;

    // Unitary vector in local y direction
    if (fabs(Vx[0]) < tolerance && fabs(Vx[1]) < tolerance) {
        MathUtils<double>::CrossProduct(Vy, VyTmp, Vx);
    }
    else {
        MathUtils<double>::CrossProduct(Vy, VzTmp, Vx);
    }

    //Rotation Matrix
    rRotationMatrix(0, 0) = Vx[0];
    rRotationMatrix(0, 1) = Vx[1];

    rRotationMatrix(1, 0) = Vy[0];
    rRotationMatrix(1, 1) = Vy[1];
}

template< >
/// <summary>
/// Rotation matrix of a line3D2 element
/// </summary>
/// <param name="rRotationMatrix"></param>
/// <param name="rGeom"></param>
void MovingLoadCondition<3, 2>::CalculateRotationMatrix(BoundedMatrix<double, 3, 3>& rRotationMatrix, const Element::GeometryType& rGeom)
{
     constexpr double tolerance = 1e-8;
  
    //Unitary vector in local x direction
    array_1d<double, 3> Vx;
    noalias(Vx) = rGeom.GetPoint(1) - rGeom.GetPoint(0);
    const double inv_norm_x = 1.0 / norm_2(Vx);
    Vx[0] *= inv_norm_x;
    Vx[1] *= inv_norm_x;
    Vx[2] *= inv_norm_x;

    //initialise vectors in local y and z direction
    array_1d<double, 3> Vy = ZeroVector(3);
    array_1d<double, 3> Vz = ZeroVector(3);

    array_1d<double, 3> VyTmp = ZeroVector(3);
    array_1d<double, 3> VzTmp = ZeroVector(3);

    VyTmp[1] = 1;
    VzTmp[2] = 1;

    // Unitary vector in local y direction
    if (fabs(Vx[0]) < tolerance && fabs(Vx[1]) < tolerance) {
        MathUtils<double>::CrossProduct(Vy, VyTmp, Vx);
    }
    else {
        MathUtils<double>::CrossProduct(Vy, VzTmp, Vx);
    }

    double inv_norm_y = 1.0 / norm_2(Vy);
    Vy[0] *= inv_norm_y;
    Vy[1] *= inv_norm_y;
    Vy[2] *= inv_norm_y;

    // Unitary vector in local z direction
    MathUtils<double>::CrossProduct(Vz, Vx, Vy);
    const double inv_norm_z = 1.0 / norm_2(Vz);
    if (inv_norm_z > tolerance)
    {
        Vz[0] *= inv_norm_z;
        Vz[1] *= inv_norm_z;
        Vz[2] *= inv_norm_z;
    }

    //Rotation Matrix
    rRotationMatrix(0, 0) = Vx[0];
    rRotationMatrix(0, 1) = Vx[1];
    rRotationMatrix(0, 2) = Vx[2];

    rRotationMatrix(1, 0) = Vy[0];
    rRotationMatrix(1, 1) = Vy[1];
    rRotationMatrix(1, 2) = Vy[2];

    rRotationMatrix(2, 0) = Vz[0];
    rRotationMatrix(2, 1) = Vz[1];
    rRotationMatrix(2, 2) = Vz[2];
}

template< std::size_t TDim, std::size_t TNumNodes >
bool MovingLoadCondition<TDim, TNumNodes>::HasRotDof() const
{
    return GetGeometry()[0].HasDofFor(ROTATION_Z) ;
}

template class MovingLoadCondition<2, 2>;
template class MovingLoadCondition<3, 2>;
template class MovingLoadCondition<2, 3>;
template class MovingLoadCondition<3, 3>;

} // Namespace Kratos


