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
#include "differentiation_utility.h"

namespace Kratos
{
    void DifferentiationUtility::CalculateRigthHandSideDerivative(Element& rElement,
                                                const Variable<double>& rDesignVariable,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        if ( rElement.GetProperties().Has(rDesignVariable) )
        {

            // define working variables
            Vector RHS_unperturbed;
            Vector RHS_perturbed;

            ProcessInfo copy_process_info = rCurrentProcessInfo;

            // Compute RHS before perturbion
            rElement.CalculateRightHandSide(RHS_unperturbed, copy_process_info);

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
            rElement.CalculateRightHandSide(RHS_perturbed, copy_process_info);

            // Compute derivative of RHS w.r.t. design variable with finite differences
            for(IndexType i = 0; i < RHS_perturbed.size(); ++i)
                rOutput(0, i) = (RHS_perturbed[i] - RHS_unperturbed[i]) / rPertubationSize;

            // Give element original properties back
            rElement.SetProperties(p_global_properties);

            //call one last time to make sure everything is as it was before TODO improve this..
            rElement.CalculateRightHandSide(RHS_perturbed, copy_process_info);
        }
        else
            if ( (rOutput.size1() != 0) || (rOutput.size2() != 0) )
                rOutput.resize(0,0,false);

        KRATOS_CATCH("");
    }

    void DifferentiationUtility::CalculateRigthHandSideDerivative(Element& rElement,
                                                const Variable<array_1d<double,3>>& rDesignVariable,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        #ifdef _OPENMP
            KRATOS_ERROR_IF(omp_get_thread_num() > 0) <<
                "DifferentiationUtility::CalculateRigthHandSideDerivative " <<
                "is not thread safe for shape derivatives!" << omp_get_thread_num();
        #endif

        if(rDesignVariable == SHAPE)
        {
            // define working variables
            Vector RHS_unperturbed;
            Vector RHS_perturbed;
            ProcessInfo copy_process_info = rCurrentProcessInfo;

            const SizeType number_of_nodes = rElement.GetGeometry().PointsNumber();
            const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);

            // compute RHS before perturbion
            rElement.CalculateRightHandSide(RHS_unperturbed, copy_process_info);

            if ( (rOutput.size1() != dimension * number_of_nodes) || (rOutput.size2() != RHS_unperturbed.size() ) )
                rOutput.resize(dimension * number_of_nodes, RHS_unperturbed.size());

            IndexType index = 0;
            for(auto& node_i : rElement.GetGeometry())
            {
                for(IndexType coord_dir_i = 0; coord_dir_i < dimension; ++coord_dir_i)
                {
                    // perturb the design variable
                    node_i.GetInitialPosition()[coord_dir_i] += rPertubationSize;
                    node_i[coord_dir_i] += rPertubationSize;

                    // compute RHS after perturbation
                    rElement.CalculateRightHandSide(RHS_perturbed, copy_process_info);

                    //compute derivative of RHS w.r.t. design variable with finite differences
                    for(IndexType i = 0; i < RHS_perturbed.size(); ++i)
                        rOutput( (coord_dir_i + index*dimension), i) = (RHS_perturbed[i]-RHS_unperturbed[i])/rPertubationSize;

                    // Reset perturbed vector
                    noalias(RHS_perturbed) = ZeroVector(RHS_perturbed.size());

                    // unperturb the design variable
                    node_i.GetInitialPosition()[coord_dir_i] -= rPertubationSize;
                    node_i[coord_dir_i] -= rPertubationSize;

                }
                index++;

                //call one last time to make sure everything is as it was before TODO improve this..
                rElement.CalculateRightHandSide(RHS_perturbed, copy_process_info);

            }// end loop over element nodes
        }
        else
        {
            KRATOS_WARNING("DifferentiationUtility") << "Unsupported nodal design variable: " << rDesignVariable << std::endl;
            if ( (rOutput.size1() != 0) || (rOutput.size2() != 0) )
                rOutput.resize(0,0,false);
        }

        KRATOS_CATCH("");
    }

        void DifferentiationUtility::CalculateRigthHandSideDerivative(Element& rElement,
                                                const array_1d_component_type& rDesignVariable,
                                                Node<3>& rNode,
                                                const double& rPertubationSize,
                                                Vector& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        #ifdef _OPENMP
            KRATOS_ERROR_IF(omp_get_thread_num() > 0) <<
                "DifferentiationUtility::CalculateRigthHandSideDerivative " <<
                "is not thread safe for shape derivatives!" << omp_get_thread_num();
        #endif

        if( rDesignVariable == SHAPE_X || rDesignVariable == SHAPE_Y || rDesignVariable == SHAPE_Z )
        {
           const IndexType coord_dir = DifferentiationUtility::GetCoordinateDirection(rDesignVariable);

            // define working variables
            Vector RHS_unperturbed;
            Vector RHS_perturbed;
            ProcessInfo copy_process_info = rCurrentProcessInfo;

            // compute RHS before perturbion
            rElement.CalculateRightHandSide(RHS_unperturbed, copy_process_info);

            if ( rOutput.size() != RHS_unperturbed.size() )
                rOutput.resize(RHS_unperturbed.size(), false);

            // perturb the design variable
            rNode.GetInitialPosition()[coord_dir] += rPertubationSize;
            rNode.Coordinates()[coord_dir] += rPertubationSize;

            // compute LHS after perturbation
            rElement.CalculateRightHandSide(RHS_perturbed, copy_process_info);

            //compute derivative of RHS w.r.t. design variable with finite differences
            noalias(rOutput) = (RHS_perturbed - RHS_unperturbed) / rPertubationSize;

             // unperturb the design variable
            rNode.GetInitialPosition()[coord_dir] -= rPertubationSize;
            rNode.Coordinates()[coord_dir] -= rPertubationSize;

            //call one last time to make sure everything is as it was before TODO improve this..
            rElement.CalculateRightHandSide(RHS_perturbed, copy_process_info);
        }
        else
        {
            KRATOS_WARNING("DifferentiationUtility") << "Unsupported nodal design variable: " << rDesignVariable << std::endl;
            if ( (rOutput.size() != 0) )
                rOutput.resize(0,false);
        }

        KRATOS_CATCH("");
    }

    void DifferentiationUtility::CalculateLeftHandSideDerivative(Element& rElement,
                                                const array_1d_component_type& rDesignVariable,
                                                Node<3>& rNode,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        #ifdef _OPENMP
            KRATOS_ERROR_IF(omp_get_thread_num() > 0) <<
                "DifferentiationUtility::CalculateLeftHandSideDerivative " <<
                "is not thread safe for shape derivatives!" << omp_get_thread_num();
        #endif

        if( rDesignVariable == SHAPE_X || rDesignVariable == SHAPE_Y || rDesignVariable == SHAPE_Z )
        {
            const IndexType coord_dir = DifferentiationUtility::GetCoordinateDirection(rDesignVariable);

            // define working variables
            Matrix LHS_unperturbed;
            Matrix LHS_perturbed;
            Vector dummy;
            ProcessInfo copy_process_info = rCurrentProcessInfo;

            // compute LHS before perturbion
            rElement.CalculateLocalSystem(LHS_unperturbed, dummy ,copy_process_info);

            if ( (rOutput.size1() != LHS_unperturbed.size1()) || (rOutput.size2() != LHS_unperturbed.size2() ) )
                rOutput.resize(LHS_unperturbed.size1(), LHS_unperturbed.size2());

            // perturb the design variable
            rNode.GetInitialPosition()[coord_dir] += rPertubationSize;
            rNode.Coordinates()[coord_dir] += rPertubationSize;

            // compute LHS after perturbation
            rElement.CalculateLocalSystem(LHS_perturbed, dummy ,copy_process_info);

            //compute derivative of RHS w.r.t. design variable with finite differences
            noalias(rOutput) = (LHS_perturbed - LHS_unperturbed) / rPertubationSize;

             // unperturb the design variable
            rNode.GetInitialPosition()[coord_dir] -= rPertubationSize;
            rNode.Coordinates()[coord_dir] -= rPertubationSize;

            //call one last time to make sure everything is as it was before TODO improve this..
            rElement.CalculateLocalSystem(LHS_perturbed, dummy ,copy_process_info);
        }
        else
        {
            KRATOS_WARNING("DifferentiationUtility") << "Unsupported nodal design variable: " << rDesignVariable << std::endl;
            if ( (rOutput.size1() != 0) || (rOutput.size2() != 0) )
                rOutput.resize(0,0,false);
        }

        KRATOS_CATCH("");
    }

    void DifferentiationUtility::CalculateMassMatrixDerivative(Element& rElement,
                                                const array_1d_component_type& rDesignVariable,
                                                Node<3>& rNode,
                                                const double& rPertubationSize,
                                                Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        #ifdef _OPENMP
            KRATOS_ERROR_IF(omp_get_thread_num() > 0) <<
                "DifferentiationUtility::CalculateMassMatrixDerivative " <<
                "is not thread safe for shape derivatives!" << omp_get_thread_num();
        #endif

        if( rDesignVariable == SHAPE_X || rDesignVariable == SHAPE_Y || rDesignVariable == SHAPE_Z )
        {
            const IndexType coord_dir = DifferentiationUtility::GetCoordinateDirection(rDesignVariable);

            // define working variables
            Matrix unperturbed_mass_matrix;
            Matrix perturbed_mass_matrix;
            ProcessInfo copy_process_info = rCurrentProcessInfo;

            // compute mass matrix before perturbion
            rElement.CalculateMassMatrix(unperturbed_mass_matrix, copy_process_info);

            if ( (rOutput.size1() != unperturbed_mass_matrix.size1()) || (rOutput.size2() != unperturbed_mass_matrix.size2() ) )
                rOutput.resize(unperturbed_mass_matrix.size1(), unperturbed_mass_matrix.size2());

            // perturb the design variable
            rNode.GetInitialPosition()[coord_dir] += rPertubationSize;
            rNode.Coordinates()[coord_dir] += rPertubationSize;

            // compute LHS after perturbation
            rElement.CalculateMassMatrix(perturbed_mass_matrix, copy_process_info);

            //compute derivative of RHS w.r.t. design variable with finite differences
            noalias(rOutput) = (perturbed_mass_matrix - unperturbed_mass_matrix) / rPertubationSize;

             // unperturb the design variable
            rNode.GetInitialPosition()[coord_dir] -= rPertubationSize;
            rNode.Coordinates()[coord_dir] -= rPertubationSize;

            //call one last time to make sure everything is as it was before TODO improve this..
            rElement.CalculateMassMatrix(perturbed_mass_matrix, copy_process_info);
        }
        else
        {
            KRATOS_WARNING("DifferentiationUtility") << "Unsupported nodal design variable: " << rDesignVariable << std::endl;
            if ( (rOutput.size1() != 0) || (rOutput.size2() != 0) )
                rOutput.resize(0,0,false);
        }

        KRATOS_CATCH("");
    }

    std::size_t DifferentiationUtility::GetCoordinateDirection(const array_1d_component_type& rDesignVariable)
    {
        IndexType coord_dir = 0;
        if( rDesignVariable == SHAPE_X )
            coord_dir = 0;
        else if( rDesignVariable == SHAPE_Y )
            coord_dir = 1;
        else if( rDesignVariable == SHAPE_Z )
            coord_dir = 2;

        return coord_dir;
    }

}  // namespace Kratos.

