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

namespace Kratos
{
    void ElementFiniteDifferenceUtility::CalculateRightHandSideDerivative(Element& rElement,
                                                const Variable<double>& rDesignVariable,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if ( rElement.GetProperties().Has(rDesignVariable) )
        {

            // define working variables
            Vector RHS_unperturbed;
            Vector RHS_perturbed;

            // Compute RHS before perturbion
            rElement.CalculateRightHandSide(RHS_unperturbed, rCurrentProcessInfo);

            if ( (rOutput.size1() != 1) || (rOutput.size2() != RHS_unperturbed.size() ) )
                rOutput.resize(1, RHS_unperturbed.size());

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
                rOutput(0, i) = (RHS_perturbed[i] - RHS_unperturbed[i]) / rPertubationSize;

            // Give element original properties back
            rElement.SetProperties(p_global_properties);

            //call one last time to make sure everything is as it was before TODO improve this..
            rElement.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
        }
        else
            if ( (rOutput.size1() != 0) || (rOutput.size2() != 0) )
                rOutput.resize(0,0,false);

        KRATOS_CATCH("");
    }

    void ElementFiniteDifferenceUtility::CalculateRightHandSideDerivative(Element& rElement,
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
                << "This function is not thread safe for shape derivatives!" << std::endl;

            #pragma omp critical
			{
                const IndexType coord_dir = ElementFiniteDifferenceUtility::GetCoordinateDirection(rDesignVariable);

                // define working variables
                Vector RHS_unperturbed;
                Vector RHS_perturbed;

                // compute RHS before perturbion
                rElement.CalculateRightHandSide(RHS_unperturbed, rCurrentProcessInfo);

                if ( rOutput.size() != RHS_unperturbed.size() )
                    rOutput.resize(RHS_unperturbed.size(), false);

                // perturb the design variable
                rNode.GetInitialPosition()[coord_dir] += rPertubationSize;
                rNode.Coordinates()[coord_dir] += rPertubationSize;

                // compute LHS after perturbation
                rElement.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);

                //compute derivative of RHS w.r.t. design variable with finite differences
                noalias(rOutput) = (RHS_perturbed - RHS_unperturbed) / rPertubationSize;

                 // unperturb the design variable
                rNode.GetInitialPosition()[coord_dir] -= rPertubationSize;
                rNode.Coordinates()[coord_dir] -= rPertubationSize;

                //call one last time to make sure everything is as it was before TODO improve this..
                rElement.CalculateRightHandSide(RHS_perturbed, rCurrentProcessInfo);
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
                                                const array_1d_component_type& rDesignVariable,
                                                Node<3>& rNode,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if( rDesignVariable == SHAPE_X || rDesignVariable == SHAPE_Y || rDesignVariable == SHAPE_Z )
        {
            KRATOS_WARNING_IF("ElementFiniteDifferenceUtility::CalculateRightHandSideDerivative", OpenMPUtils::IsInParallel() != 0)
                << "This function is not thread safe for shape derivatives!" << std::endl;

            #pragma omp critical
			{
                const IndexType coord_dir = ElementFiniteDifferenceUtility::GetCoordinateDirection(rDesignVariable);

                // define working variables
                Matrix LHS_unperturbed;
                Matrix LHS_perturbed;
                Vector dummy;

                // compute LHS before perturbion
                rElement.CalculateLocalSystem(LHS_unperturbed, dummy, rCurrentProcessInfo);

                if ( (rOutput.size1() != LHS_unperturbed.size1()) || (rOutput.size2() != LHS_unperturbed.size2() ) )
                    rOutput.resize(LHS_unperturbed.size1(), LHS_unperturbed.size2());

                // perturb the design variable
                rNode.GetInitialPosition()[coord_dir] += rPertubationSize;
                rNode.Coordinates()[coord_dir] += rPertubationSize;

                // compute LHS after perturbation
                rElement.CalculateLocalSystem(LHS_perturbed, dummy, rCurrentProcessInfo);

                //compute derivative of RHS w.r.t. design variable with finite differences
                noalias(rOutput) = (LHS_perturbed - LHS_unperturbed) / rPertubationSize;

                 // unperturb the design variable
                rNode.GetInitialPosition()[coord_dir] -= rPertubationSize;
                rNode.Coordinates()[coord_dir] -= rPertubationSize;

                //call one last time to make sure everything is as it was before TODO improve this..
                rElement.CalculateLocalSystem(LHS_perturbed, dummy, rCurrentProcessInfo);
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
                                                const array_1d_component_type& rDesignVariable,
                                                Node<3>& rNode,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if( rDesignVariable == SHAPE_X || rDesignVariable == SHAPE_Y || rDesignVariable == SHAPE_Z )
        {
            KRATOS_WARNING_IF("ElementFiniteDifferenceUtility::CalculateRightHandSideDerivative", OpenMPUtils::IsInParallel() != 0)
                << "This function is not thread safe for shape derivatives!" << std::endl;

            #pragma omp critical
			{
                const IndexType coord_dir = ElementFiniteDifferenceUtility::GetCoordinateDirection(rDesignVariable);

                // define working variables
                Matrix unperturbed_mass_matrix;
                Matrix perturbed_mass_matrix;

                // compute mass matrix before perturbion
                rElement.CalculateMassMatrix(unperturbed_mass_matrix, rCurrentProcessInfo);

                if ( (rOutput.size1() != unperturbed_mass_matrix.size1()) || (rOutput.size2() != unperturbed_mass_matrix.size2() ) )
                    rOutput.resize(unperturbed_mass_matrix.size1(), unperturbed_mass_matrix.size2());

                // perturb the design variable
                rNode.GetInitialPosition()[coord_dir] += rPertubationSize;
                rNode.Coordinates()[coord_dir] += rPertubationSize;

                // compute LHS after perturbation
                rElement.CalculateMassMatrix(perturbed_mass_matrix, rCurrentProcessInfo);

                //compute derivative of RHS w.r.t. design variable with finite differences
                noalias(rOutput) = (perturbed_mass_matrix - unperturbed_mass_matrix) / rPertubationSize;

                 // unperturb the design variable
                rNode.GetInitialPosition()[coord_dir] -= rPertubationSize;
                rNode.Coordinates()[coord_dir] -= rPertubationSize;

                //call one last time to make sure everything is as it was before TODO improve this..
                rElement.CalculateMassMatrix(perturbed_mass_matrix, rCurrentProcessInfo);
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

