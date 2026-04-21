// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

// System includes

// External includes

// Project includes
#include "finite_difference_utility.h"
#include "utilities/openmp_utils.h"

namespace Kratos {


template <IAdjoint::ResidualTerm Term, class TEntity>
requires (std::is_same_v<TEntity,Element> || std::is_same_v<TEntity,Condition>)
void FiniteDifferenceUtility::ComputeDerivative(
    const TEntity& rEntity,
    const std::conditional_t<
        Term == IAdjoint::ResidualTerm::Load,
        Vector,
        Matrix
    >& rValue,
    std::span<const std::pair<
        IAdjoint::DynamicVariable,
        Globals::DataLocation>
    > variables,
    Matrix& rOutput,
    double Perturbation,
    const ProcessInfo& rProcessInfo) {
        KRATOS_TRY

        // This wonderful piece of software does not care about petty concepts
        // such as const-correctness, so I have to do a disgusting const cast here.
        TEntity& r_mutable_entity = *const_cast<TEntity*>(&rEntity);

        // Make a copy of the entity, its geometry, nodes, properties and constitutive law.
        // Since the naming and design in Kratos is deliberately misleading, TEntity::Clone
        // does not actually perform a deep copy of the entity, so the above mentioned
        // components have to be copied separately.
        typename TEntity::Pointer p_entity_clone;

        {
            // Copy nodes.
            typename TEntity::NodesArrayType nodes;
            for (Node& r_node : r_mutable_entity.GetGeometry())
                nodes.push_back(r_node.Clone(r_node.Id()));

            // Geometry does not have clone function even though it's a base class.
            // The best I can do here is to construct a new one and hope the instance
            // does not have any exotic internal state.
            Geometry<Node>::Pointer p_geometry_clone = rEntity.GetGeometry().Create(
                rEntity.GetGeometry().Id(),
                rEntity.GetGeometry());

            // Copy the DataValueContainer attached to the geometry.
            p_geometry_clone->GetData() = rEntity.GetGeometry().GetData();

            // Copy the entity's properties.
            Properties::Pointer p_properties_clone(new Properties(rEntity.GetProperties()));

            // Clone the entity.
            //p_entity_clone = rEntity.Clone(
            //    rEntity.Id(),
            //    )
        }

        KRATOS_CATCH("")
}


#define KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(Term, TEntity)         \
    template void FiniteDifferenceUtility::ComputeDerivative<Term,TEntity>( \
        const TEntity&,                                                     \
        const std::conditional_t<                                           \
            Term == IAdjoint::ResidualTerm::Load,                           \
            Vector,                                                         \
            Matrix                                                          \
        >&,                                                                 \
        std::span<const std::pair<                                          \
            IAdjoint::DynamicVariable,                                      \
            Globals::DataLocation>>,                                        \
        Matrix&,                                                            \
        double,                                                             \
        const ProcessInfo&);


KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Load, Element)
KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Mass, Element)
KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Damping, Element)
KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Stiffness, Element)
KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Load, Condition)
KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Mass, Condition)
KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Damping, Condition)
KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY(IAdjoint::ResidualTerm::Stiffness, Condition)


#undef KRATOS_INSTANTIATE_FINITE_DIFFERENCE_UTILITY


    void FiniteDifferenceUtility::CalculateRightHandSideDerivative(Element& rElement,
                                                const Vector& rRHS,
                                                const Variable<double>& rDesignVariable,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if ( rElement.GetProperties().Has(rDesignVariable) )
        {

            // define working variables
            Vector RHS_perturbed;

            if ( (rOutput.size1() != 1) || (rOutput.size2() != rRHS.size() ) )
                rOutput.resize(1, rRHS.size(), false);

            // Save property pointer
            Properties::Pointer p_global_properties = rElement.pGetProperties();

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
            rElement.SetProperties(p_local_property);

            // perturb the design variable
            const double current_property_value = rElement.GetProperties()[rDesignVariable];
            p_local_property->SetValue(rDesignVariable, (current_property_value + rPertubationSize));

            // Compute RHS after perturbation
            rElement.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);

            // Compute derivative of RHS w.r.t. design variable with finite differences
            for(IndexType i = 0; i < RHS_perturbed.size(); ++i)
                rOutput(0, i) = (RHS_perturbed[i] - rRHS[i]) / rPertubationSize;

            // Give element original properties back
            rElement.SetProperties(p_global_properties);
        }
        else
            if ( (rOutput.size1() != 0) || (rOutput.size2() != 0) )
                rOutput.resize(0,0,false);

        KRATOS_CATCH("");
    }

    void FiniteDifferenceUtility::CalculateLeftHandSideDerivative(Element& rElement,
                                                const Matrix& rLHS,
                                                const array_1d_component_type& rDesignVariable,
                                                Node& rNode,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if( rDesignVariable == SHAPE_SENSITIVITY_X || rDesignVariable == SHAPE_SENSITIVITY_Y || rDesignVariable == SHAPE_SENSITIVITY_Z )
        {
            KRATOS_WARNING_IF("FiniteDifferenceUtility::CalculateLeftHandSideDerivative", OpenMPUtils::IsInParallel() != 0)
                << "The call of this non shared-memory-parallelized function within a parallel section should be avoided for efficiency reasons!" << std::endl;

            #pragma omp critical
            {
                const IndexType coord_dir = FiniteDifferenceUtility::GetCoordinateDirection(rDesignVariable);

                // define working variables
                Matrix LHS_perturbed;
                Vector dummy;

                if ( (rOutput.size1() != rLHS.size1()) || (rOutput.size2() != rLHS.size2() ) )
                    rOutput.resize(rLHS.size1(), rLHS.size2(), false);

                // perturb the design variable
                rNode.GetInitialPosition()[coord_dir] += rPertubationSize;
                rNode.Coordinates()[coord_dir] += rPertubationSize;

                // compute LHS after perturbation
                rElement.CalculateLocalSystem(LHS_perturbed, dummy, rCurrentProcessInfo);

                //compute derivative of RHS w.r.t. design variable with finite differences
                noalias(rOutput) = (LHS_perturbed - rLHS) / rPertubationSize;

                 // unperturb the design variable
                rNode.GetInitialPosition()[coord_dir] -= rPertubationSize;
                rNode.Coordinates()[coord_dir] -= rPertubationSize;
            }
        }
        else
        {
            KRATOS_WARNING("FiniteDifferenceUtility") << "Unsupported nodal design variable: " << rDesignVariable << std::endl;
            if ( (rOutput.size1() != 0) || (rOutput.size2() != 0) )
                rOutput.resize(0,0,false);
        }

        KRATOS_CATCH("");
    }

    void FiniteDifferenceUtility::CalculateMassMatrixDerivative(Element& rElement,
                                                const Matrix& rMassMatrix,
                                                const array_1d_component_type& rDesignVariable,
                                                Node& rNode,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if( rDesignVariable == SHAPE_SENSITIVITY_X || rDesignVariable == SHAPE_SENSITIVITY_Y || rDesignVariable == SHAPE_SENSITIVITY_Z )
        {
            KRATOS_WARNING_IF("FiniteDifferenceUtility::CalculateMassMatrixDerivative", OpenMPUtils::IsInParallel() != 0)
                << "The call of this non shared-memory-parallelized function within a parallel section should be avoided for efficiency reasons!" << std::endl;

            #pragma omp critical
            {
                const IndexType coord_dir = FiniteDifferenceUtility::GetCoordinateDirection(rDesignVariable);

                // define working variables
                Matrix perturbed_mass_matrix;

                if ( (rOutput.size1() != rMassMatrix.size1()) || (rOutput.size2() != rMassMatrix.size2() ) )
                    rOutput.resize(rMassMatrix.size1(), rMassMatrix.size2(), false);

                // perturb the design variable
                rNode.GetInitialPosition()[coord_dir] += rPertubationSize;
                rNode.Coordinates()[coord_dir] += rPertubationSize;

                // compute LHS after perturbation
                rElement.CalculateMassMatrix(perturbed_mass_matrix, rCurrentProcessInfo);

                //compute derivative of RHS w.r.t. design variable with finite differences
                noalias(rOutput) = (perturbed_mass_matrix - rMassMatrix) / rPertubationSize;

                 // unperturb the design variable
                rNode.GetInitialPosition()[coord_dir] -= rPertubationSize;
                rNode.Coordinates()[coord_dir] -= rPertubationSize;
            }
        }
        else
        {
            KRATOS_WARNING("FiniteDifferenceUtility") << "Unsupported nodal design variable: " << rDesignVariable << std::endl;
            if ( (rOutput.size1() != 0) || (rOutput.size2() != 0) )
                rOutput.resize(0,0,false);
        }

        KRATOS_CATCH("");
    }

    std::size_t FiniteDifferenceUtility::GetCoordinateDirection(const array_1d_component_type& rDesignVariable)
    {
        if( rDesignVariable == SHAPE_SENSITIVITY_X )
            return 0;
        else if( rDesignVariable == SHAPE_SENSITIVITY_Y )
            return 1;
        else if( rDesignVariable == SHAPE_SENSITIVITY_Z )
            return 2;
        else
            KRATOS_ERROR << "Invalid valiable component: " << rDesignVariable.Name() <<
                "Available is only 'SHAPE_SENSITIVITY_X','SHAPE_SENSITIVITY_Y' and 'SHAPE_SENSITIVITY_Z' " << std::endl;
    }

}  // namespace Kratos.

