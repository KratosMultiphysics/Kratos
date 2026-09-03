// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// --- STL Includes ---
#include <ranges>


// Project includes
#include "custom_conditions/point_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

PointLoadCondition::PointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseLoadCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

PointLoadCondition::PointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : BaseLoadCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer PointLoadCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<PointLoadCondition>(NewId, pGeom, pProperties);
}

//************************************************************************************
//************************************************************************************

Condition::Pointer PointLoadCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<PointLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PointLoadCondition::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<PointLoadCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

PointLoadCondition::~PointLoadCondition()
{
}

//************************************************************************************
//************************************************************************************

void PointLoadCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    const unsigned int NumberOfNodes = GetGeometry().size();
    const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    const unsigned int MatSize = NumberOfNodes * Dimension;

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
    array_1d<double, 3 > PointLoad = ZeroVector(3);
    if( this->Has( POINT_LOAD ) )
    {
        noalias(PointLoad) = this->GetValue( POINT_LOAD );
    }

    for (unsigned int ii = 0; ii < NumberOfNodes; ++ii)
    {
        const unsigned int base = ii*Dimension;

        if( GetGeometry()[ii].SolutionStepsDataHas( POINT_LOAD ) )
        {
            noalias(PointLoad) += GetGeometry()[ii].FastGetSolutionStepValue( POINT_LOAD );
        }

        for(unsigned int k = 0; k < Dimension; ++k)
        {
            rRightHandSideVector[base + k] += GetPointLoadIntegrationWeight() * PointLoad[k];
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

double PointLoadCondition::GetPointLoadIntegrationWeight() const
{
    return 1.0;
}


void PointLoadCondition::GetLoadInfluencingVariables(
    std::vector<IAdjoint::DynamicVariable>& rOutput,
    const ProcessInfo&) const {
        KRATOS_DEBUG_ERROR_IF_NOT(this->GetGeometry().size() == 1);
        KRATOS_TRY
            const std::array<const Variable<double>*,3> point_load_components {
                &POINT_LOAD_X,
                &POINT_LOAD_Y,
                &POINT_LOAD_Z};
            const auto dimension = GetGeometry().WorkingSpaceDimension();

            rOutput.clear();
            std::ranges::transform(
                point_load_components | std::views::take(dimension),
                std::back_inserter(rOutput),
                [] (const Variable<double>* p_variable) -> IAdjoint::DynamicVariable {
                    return {*p_variable, 0};
                });
        KRATOS_CATCH("")
}


void PointLoadCondition::ComputeLoadDerivative(
    Matrix& rOutput,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const ProcessInfo&,
    int) const {
        KRATOS_DEBUG_ERROR_IF_NOT(this->GetGeometry().size() == 1);
        KRATOS_TRY
            const std::array<IAdjoint::DynamicVariable,3> all_point_load_components {
                IAdjoint::DynamicVariable {POINT_LOAD_X, 0},
                IAdjoint::DynamicVariable {POINT_LOAD_Y, 0},
                IAdjoint::DynamicVariable {POINT_LOAD_Z, 0}};
            const auto dimension = GetGeometry().WorkingSpaceDimension();
            const auto point_load_components = all_point_load_components | std::views::take(dimension);

            rOutput = ZeroMatrix(dimension, Variables.size());
            for (std::size_t i_variable=0ul; i_variable<Variables.size(); ++i_variable) {
                // Find which point load component the current variable refers to.
                const auto it_point_load_component = std::ranges::find(
                    point_load_components,
                    Variables[i_variable]);
                KRATOS_ERROR_IF(it_point_load_component == point_load_components.end())
                    << "unexpected variable " << Variables[i_variable].Name() << " "
                    << Variables[i_variable].GetDynamicIndex();
                const auto i_dimension = std::distance(
                    point_load_components.begin(),
                    it_point_load_component);

                // Update the output matrix.
                rOutput(i_dimension, i_variable) += 1.0;
            } // for i_variable in range(Variables.size())
        KRATOS_CATCH("")
}


} // namespace Kratos
