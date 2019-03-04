//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//


// System includes

// External includes

// Project includes
#include "element_finite_difference_utility.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{
    void ElementFiniteDifferenceUtility::CalculateRightHandSideDerivative(Element& rElement,
                                                const Vector& rRHS,
                                                const Variable<double>& rDesignVariable,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if ( rElement.GetProperties().Has(rDesignVariable) )
        {

            // define working variables
            Vector RHS_perturbed;

            if ( (rOutput.size1() != 1) || (rOutput.size2() != rRHS.size() ) )
                rOutput.resize(1, rRHS.size());

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

    void ElementFiniteDifferenceUtility::CalculateRightHandSideDerivative(Element& rElement,
                                                const Vector& rRHS,
                                                const array_1d_component_type& rDesignVariable,
                                                Node<3>& rNode,
                                                const double& rPertubationSize,
                                                Vector& rOutput,
                                                ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if( rDesignVariable == SHAPE_X || rDesignVariable == SHAPE_Y || rDesignVariable == SHAPE_Z )
        {
            KRATOS_WARNING_IF("ElementFiniteDifferenceUtility::CalculateRightHandSideDerivative", OpenMPUtils::IsInParallel() != 0)
                << "The call of this non omp-parallelized function within a parallel section should be avoided for efficiency reasons!" << std::endl;

            #pragma omp critical
            {
                const IndexType coord_dir = ElementFiniteDifferenceUtility::GetCoordinateDirection(rDesignVariable);

                // define working variables
                Vector RHS_perturbed;

                if ( rOutput.size() != rRHS.size() )
                    rOutput.resize(rRHS.size(), false);

                // perturb the design variable
                rNode.GetInitialPosition()[coord_dir] += rPertubationSize;
                rNode.Coordinates()[coord_dir] += rPertubationSize;

                // compute LHS after perturbation
                rElement.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);

                //compute derivative of RHS w.r.t. design variable with finite differences
                noalias(rOutput) = (RHS_perturbed - rRHS) / rPertubationSize;

                 // unperturb the design variable
                rNode.GetInitialPosition()[coord_dir] -= rPertubationSize;
                rNode.Coordinates()[coord_dir] -= rPertubationSize;
            }
        }
        else
        {
            KRATOS_WARNING("ElementFiniteDifferenceUtility") << "Unsupported nodal design variable: " << rDesignVariable << std::endl;
            if ( (rOutput.size() != 0) )
                rOutput.resize(0,false);
        }

        KRATOS_CATCH("");
    }

    void ElementFiniteDifferenceUtility::CalculateLeftHandSideDerivative(Element& rElement,
                                                const Matrix& rLHS,
                                                const array_1d_component_type& rDesignVariable,
                                                Node<3>& rNode,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if( rDesignVariable == SHAPE_X || rDesignVariable == SHAPE_Y || rDesignVariable == SHAPE_Z )
        {
            KRATOS_WARNING_IF("ElementFiniteDifferenceUtility::CalculateLeftHandSideDerivative", OpenMPUtils::IsInParallel() != 0)
                << "The call of this non omp-parallelized function within a parallel section should be avoided for efficiency reasons!" << std::endl;

            #pragma omp critical
            {
                const IndexType coord_dir = ElementFiniteDifferenceUtility::GetCoordinateDirection(rDesignVariable);

                // define working variables
                Matrix LHS_perturbed;
                Vector dummy;

                if ( (rOutput.size1() != rLHS.size1()) || (rOutput.size2() != rLHS.size2() ) )
                    rOutput.resize(rLHS.size1(), rLHS.size2());

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
            KRATOS_WARNING("ElementFiniteDifferenceUtility") << "Unsupported nodal design variable: " << rDesignVariable << std::endl;
            if ( (rOutput.size1() != 0) || (rOutput.size2() != 0) )
                rOutput.resize(0,0,false);
        }

        KRATOS_CATCH("");
    }

    void ElementFiniteDifferenceUtility::CalculateMassMatrixDerivative(Element& rElement,
                                                const Matrix& rMassMatrix,
                                                const array_1d_component_type& rDesignVariable,
                                                Node<3>& rNode,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if( rDesignVariable == SHAPE_X || rDesignVariable == SHAPE_Y || rDesignVariable == SHAPE_Z )
        {
            KRATOS_WARNING_IF("ElementFiniteDifferenceUtility::CalculateMassMatrixDerivative", OpenMPUtils::IsInParallel() != 0)
                << "The call of this non omp-parallelized function within a parallel section should be avoided for efficiency reasons!" << std::endl;

            #pragma omp critical
            {
                const IndexType coord_dir = ElementFiniteDifferenceUtility::GetCoordinateDirection(rDesignVariable);

                // define working variables
                Matrix perturbed_mass_matrix;

                if ( (rOutput.size1() != rMassMatrix.size1()) || (rOutput.size2() != rMassMatrix.size2() ) )
                    rOutput.resize(rMassMatrix.size1(), rMassMatrix.size2());

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
            KRATOS_WARNING("ElementFiniteDifferenceUtility") << "Unsupported nodal design variable: " << rDesignVariable << std::endl;
            if ( (rOutput.size1() != 0) || (rOutput.size2() != 0) )
                rOutput.resize(0,0,false);
        }

        KRATOS_CATCH("");
    }

    std::size_t ElementFiniteDifferenceUtility::GetCoordinateDirection(const array_1d_component_type& rDesignVariable)
    {
        if( rDesignVariable == SHAPE_X )
            return 0;
        else if( rDesignVariable == SHAPE_Y )
            return 1;
        else if( rDesignVariable == SHAPE_Z )
            return 2;
        else
            KRATOS_ERROR << "Invalid valiable component: " << rDesignVariable.Name() <<
                "Available is only 'SHAPE_X','SHAPE_Y' and 'SHAPE_Z' " << std::endl;
    }

}  // namespace Kratos.

