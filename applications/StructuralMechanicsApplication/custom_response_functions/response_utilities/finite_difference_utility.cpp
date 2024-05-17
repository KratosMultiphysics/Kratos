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

namespace Kratos
{
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
        if( rDesignVariable == SHAPE_SENSITIVITY_X || rDesignVariable == TEMPERATURE)
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

