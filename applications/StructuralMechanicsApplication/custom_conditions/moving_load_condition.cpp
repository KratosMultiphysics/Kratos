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

#include <includes/checks.h>

#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

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
    array_1d<double, 3 > MovingLoad = ZeroVector(3);
    if( this->Has(MOVING_LOAD ) )
    {
        noalias(MovingLoad) = this->GetValue( MOVING_LOAD );
    }

    // check if cond should be calculated
    bool is_moving_load = false;
    for (int i = 0; i < TDim; ++i)
    {
        if (abs(MovingLoad[i]) > DBL_EPSILON)
        {
            is_moving_load = true;
        }
    }

    // apply moving load if moving load is present
    if (is_moving_load)
    {

        double local_x_coord = this->GetValue(MOVING_LOAD_LOCAL_DISTANCE);

        VectorType translationalShapeFunctionsVector;
        VectorType rotationalShapeFunctionsVector;

        CalculateSecondOrderTranslationalShapeFunctions(translationalShapeFunctionsVector, local_x_coord);
        CalculateSecondOrderRotationalShapeFunctions(rotationalShapeFunctionsVector, local_x_coord);

        GeometryType& rGeom = this->GetGeometry();

        array_1d<double, TDim> load_first_node = ZeroVector(TDim);
        array_1d<double, TDim> load_end_node = ZeroVector(TDim);

        std::vector<double>  rotation_first_node;
        std::vector<double>  rotation_end_node;

        rotation_first_node.resize(block_size - TDim, false);
        rotation_end_node.resize(block_size - TDim, false);

        for (unsigned int ii = 1; ii < TDim; ++ii)
        {
            load_first_node[ii] = translationalShapeFunctionsVector[0] * MovingLoad[ii];
            load_end_node[ii] = translationalShapeFunctionsVector[1] * MovingLoad[ii];
        }

        
        if (block_size > TDim + 1)
        {
            rotation_first_node[0] = 0;
            rotation_first_node[1] = rotationalShapeFunctionsVector[0] * MovingLoad[2];
            rotation_first_node[2] = rotationalShapeFunctionsVector[0] * MovingLoad[1];

            rotation_end_node[0] = 0;
            rotation_end_node[1] = rotationalShapeFunctionsVector[1] * MovingLoad[2];
            rotation_end_node[2] = rotationalShapeFunctionsVector[1] * MovingLoad[1];

        }
        else if (block_size > TDim)
        {
            rotation_first_node[0] = rotationalShapeFunctionsVector[0] * MovingLoad[1];
            rotation_end_node[0] = rotationalShapeFunctionsVector[0] * MovingLoad[1];
        }
        

        array_1d<int, 2> shape_indices;
        shape_indices[0] = 0;
        shape_indices[1] = NumberOfNodes - 1;

        for (unsigned int ii = 0; ii < shape_indices.size(); ++ii)
        {
            const unsigned int base = shape_indices[ii] * block_size;

            // only add load and rotation to RHS if current node is first or final node
            for (unsigned int k = 0; k < Dimension; ++k)
            {
                const double load = (shape_indices[ii] == 0) ? load_first_node[k] : load_end_node[k];
                rRightHandSideVector[base + k] = load;
            }


            // shape index is 0 or 1, rotation is only added to first and final node
            for (unsigned int k = 0; k < block_size - Dimension; ++k)
            {
                //const double load = (shape_indices[ii] == 0) ? load_first_node[k] : load_end_node[ii];
                //double rotation = rotationalShapeFunctionsVector[shape_indices[ii]] * MovingLoad[0];

                const double moment = (shape_indices[ii] == 0) ? rotation_first_node[k] : rotation_first_node[k];

                rRightHandSideVector[base + TDim +k] = moment;
            }
        }
    }
    KRATOS_CATCH( "" )
}

template< std::size_t TDim, std::size_t TNumNodes >
void MovingLoadCondition< TDim, TNumNodes>::CalculateSecondOrderTranslationalShapeFunctions(VectorType& rShapeFunctionsVector, double local_x_coord) const
{
    if (rShapeFunctionsVector.size() != 2) {
        rShapeFunctionsVector.resize(2, false);
    }
    auto& rGeom = this->GetGeometry();
    double length = rGeom.Length();

    rShapeFunctionsVector[0] = 1 + 2 * std::pow((local_x_coord / length), 3) - 3 * std::pow((local_x_coord / length), 2);
    rShapeFunctionsVector[1] = -2 * std::pow((local_x_coord / length), 3) + 3 * std::pow((local_x_coord / length), 2);
}

template< std::size_t TDim, std::size_t TNumNodes >
void MovingLoadCondition< TDim, TNumNodes>::CalculateSecondOrderRotationalShapeFunctions(VectorType& rShapeFunctionsVector, double local_x_coord) const
{
    if (rShapeFunctionsVector.size() != 2) {
        rShapeFunctionsVector.resize(2, false);
    }
    auto& rGeom = this->GetGeometry();
    double length = rGeom.Length();

    rShapeFunctionsVector[0] = local_x_coord + (pow(local_x_coord, 3) / pow(length, 2)) - 2 * (pow(local_x_coord, 2) / length);
    rShapeFunctionsVector[1] = (pow(local_x_coord, 3) / pow(length, 2)) - (pow(local_x_coord, 2) / length);
}

//************************************************************************************
//************************************************************************************
template< std::size_t TDim, std::size_t TNumNodes >
double MovingLoadCondition< TDim, TNumNodes>::GetMovingLoadIntegrationWeight( ) const
{
    return 1.0;
}


template class MovingLoadCondition<2, 2>;
template class MovingLoadCondition<3, 2>;
template class MovingLoadCondition<2, 3>;
template class MovingLoadCondition<3, 3>;

} // Namespace Kratos


