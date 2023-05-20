//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//  Kratos default license: kratos/license.txt
//
//  Main authors:    RAUL BRAVO
//

#if !defined( ROM_RESIDUALS_UTILITY_H_INCLUDED )
#define  ROM_RESIDUALS_UTILITY_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

/* Application includes */
#include "rom_application_variables.h"
#include "custom_utilities/rom_auxiliary_utilities.h"


namespace Kratos
{
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;

    // This utility returns the converged residuals projected onto the ROM basis Phi.
    class RomResidualsUtility
    {
        public:

        KRATOS_CLASS_POINTER_DEFINITION(RomResidualsUtility);

        RomResidualsUtility(
        ModelPart& rModelPart,
        Parameters ThisParameters,
        BaseSchemeType::Pointer pScheme
        ): mrModelPart(rModelPart), mpScheme(pScheme){
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "nodal_unknowns" : [],
            "number_of_rom_dofs" : 10,
            "petrov_galerkin_number_of_rom_dofs" : 10
        })" );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        mNodalVariablesNames = ThisParameters["nodal_unknowns"].GetStringArray();

        mNodalDofs = mNodalVariablesNames.size();
        mRomDofs = ThisParameters["number_of_rom_dofs"].GetInt();
        mPetrovGalerkinRomDofs = ThisParameters["petrov_galerkin_number_of_rom_dofs"].GetInt();

        // Setting up mapping: VARIABLE_KEY --> CORRECT_ROW_IN_BASIS
        for(int k=0; k<mNodalDofs; k++){
            if(KratosComponents<Variable<double>>::Has(mNodalVariablesNames[k]))
            {
                const auto& var = KratosComponents<Variable<double>>::Get(mNodalVariablesNames[k]);
                mMapPhi[var.Key()] = k;
            }
            else
                KRATOS_ERROR << "variable \""<< mNodalVariablesNames[k] << "\" not valid" << std::endl;

            }
        }

        ~RomResidualsUtility()= default;

        /**
        * Thread Local Storage containing dynamically allocated structures to avoid reallocating each iteration.
        */
        struct AssemblyTLS
        {
            Matrix phiE = {};                // Elemental Phi
            Matrix lhs = {};  // Elemental LHS
            Element::EquationIdVectorType eq_id = {}; // Elemental equation ID vector
            Element::DofsVectorType dofs = {};        // Elemental dof vector
            Matrix romA;        // reduced LHS
            Vector romB;        // reduced RHS
        };

        /**
         * Resizes a Matrix if it's not the right size
         */
        template<typename TMatrix>
        static void ResizeIfNeeded(TMatrix& rMat, const std::size_t Rows, const std::size_t Cols)

        {
            if(rMat.size1() != Rows || rMat.size2() != Cols) {
                rMat.resize(Rows, Cols, false);
            }
        };

        Matrix GetProjectedResidualsOntoPhi()
        {
            // Getting the number of elements and conditions from the model
            const int n_elements = static_cast<int>(mrModelPart.Elements().size());
            const int n_conditions = static_cast<int>(mrModelPart.Conditions().size());

            const auto& r_current_process_info = mrModelPart.GetProcessInfo();
            const auto el_begin = mrModelPart.ElementsBegin();
            const auto cond_begin = mrModelPart.ConditionsBegin();

            //contributions to the system
            Matrix lhs_contribution = ZeroMatrix(0, 0);
            Vector rhs_contribution = ZeroVector(0);

            //vector containing the localization in the system of the different terms
            Element::EquationIdVectorType equation_id;
            Matrix matrix_residuals( (n_elements + n_conditions), mRomDofs); // Matrix of reduced residuals.
            Matrix phi_elemental;

            //dofs container initialization
            Element::DofsVectorType elem_dofs;
            Condition::DofsVectorType cond_dofs;
            #pragma omp parallel firstprivate(n_elements, n_conditions, lhs_contribution, rhs_contribution, equation_id, phi_elemental, el_begin, cond_begin, elem_dofs, cond_dofs)
            {
                #pragma omp for nowait
                for (int k = 0; k < n_elements; k++){
                    const auto it_el = el_begin + k;
                    //detect if the element is active or not. If the user did not make any choice the element is active by default
                    bool element_is_active = true;
                    if ((it_el)->IsDefined(ACTIVE))
                        element_is_active = (it_el)->Is(ACTIVE);
                    if (element_is_active){
                        //calculate elemental contribution
                        mpScheme->CalculateSystemContributions(*it_el, lhs_contribution, rhs_contribution, equation_id, r_current_process_info);
                        it_el->GetDofList(elem_dofs, r_current_process_info);
                        //assemble the elemental contribution - here is where the ROM acts
                        //compute the elemental reduction matrix phi_elemental
                        const auto& r_geom = it_el->GetGeometry();
                        if(phi_elemental.size1() != elem_dofs.size() || phi_elemental.size2() != mRomDofs)
                            phi_elemental.resize(elem_dofs.size(), mRomDofs,false);
                        RomAuxiliaryUtilities::GetPhiElemental(phi_elemental, elem_dofs, r_geom, mMapPhi);
                        noalias(row(matrix_residuals, k)) = prod(trans(phi_elemental), rhs_contribution); // The size of the residual will vary only when using more ROM modes, one row per condition
                    }

                }

                #pragma omp for nowait
                for (int k = 0; k < n_conditions;  k++){
                    ModelPart::ConditionsContainerType::iterator it = cond_begin + k;
                    //detect if the condition is active or not. If the user did not make any choice the condition is active by default
                    bool condition_is_active = true;
                    if ((it)->IsDefined(ACTIVE))
                        condition_is_active = (it)->Is(ACTIVE);
                    if (condition_is_active){
                        it->GetDofList(cond_dofs, r_current_process_info);
                        //calculate elemental contribution
                        mpScheme->CalculateSystemContributions(*it, lhs_contribution, rhs_contribution, equation_id, r_current_process_info);
                        //assemble the elemental contribution - here is where the ROM acts
                        //compute the elemental reduction matrix phi_elemental
                        const auto& r_geom = it->GetGeometry();
                        if(phi_elemental.size1() != cond_dofs.size() || phi_elemental.size2() != mRomDofs)
                            phi_elemental.resize(cond_dofs.size(), mRomDofs,false);
                        RomAuxiliaryUtilities::GetPhiElemental(phi_elemental, cond_dofs, r_geom, mMapPhi);
                        noalias(row(matrix_residuals, k+n_elements)) = prod(trans(phi_elemental), rhs_contribution); // The size of the residual will vary only when using more ROM modes, one row per condition
                    }
                }
            }
        return matrix_residuals;
        }

        Matrix GetProjectedResidualsOntoPsi()
        {
            // Getting the number of elements and conditions from the model
            const int n_elements = static_cast<int>(mrModelPart.Elements().size());
            const int n_conditions = static_cast<int>(mrModelPart.Conditions().size());

            const auto& r_current_process_info = mrModelPart.GetProcessInfo();
            const auto el_begin = mrModelPart.ElementsBegin();
            const auto cond_begin = mrModelPart.ConditionsBegin();

            //contributions to the system
            Matrix lhs_contribution;
            Vector rhs_contribution;

            //vector containing the localization in the system of the different terms
            Element::EquationIdVectorType equation_id;
            Matrix matrix_residuals( (n_elements + n_conditions), mPetrovGalerkinRomDofs); // Matrix of reduced residuals.
            Matrix psi_elemental;
            
            //dofs container initialization
            Element::DofsVectorType elem_dofs;
            Condition::DofsVectorType cond_dofs;
            #pragma omp parallel firstprivate(n_elements, n_conditions, lhs_contribution, rhs_contribution, equation_id, psi_elemental, el_begin, cond_begin, elem_dofs, cond_dofs)
            {
                #pragma omp for nowait
                for (int k = 0; k < n_elements; k++){
                    const auto it_el = el_begin + k;
                    //detect if the element is active or not. If the user did not make any choice the element is active by default
                    const bool element_is_active = it_el->IsDefined(ACTIVE) ? it_el->Is(ACTIVE) : true;
                    if (element_is_active){
                        //calculate elemental contribution
                        mpScheme->CalculateSystemContributions(*it_el, lhs_contribution, rhs_contribution, equation_id, r_current_process_info);
                        it_el->GetDofList(elem_dofs, r_current_process_info);
                        //assemble the elemental contribution - here is where the ROM acts
                        //compute the elemental reduction matrix phi_elemental
                        const auto& r_geom = it_el->GetGeometry();
                        if(psi_elemental.size1() != elem_dofs.size() || psi_elemental.size2() != mPetrovGalerkinRomDofs)
                            psi_elemental.resize(elem_dofs.size(), mPetrovGalerkinRomDofs,false);
                        RomAuxiliaryUtilities::GetPsiElemental(psi_elemental, elem_dofs, r_geom, mMapPhi);
                        noalias(row(matrix_residuals, k)) = prod(trans(psi_elemental), rhs_contribution); // The size of the residual will vary only when using more ROM modes, one row per condition
                    }

                }

                #pragma omp for nowait
                for (int k = 0; k < n_conditions;  k++){
                    const auto it = cond_begin + k;
                    //detect if the condition is active or not. If the user did not make any choice the condition is active by default
                    const bool condition_is_active = it->IsDefined(ACTIVE) ? it->Is(ACTIVE) : true;
                    if (condition_is_active){
                        it->GetDofList(cond_dofs, r_current_process_info);
                        //calculate elemental contribution
                        mpScheme->CalculateSystemContributions(*it, lhs_contribution, rhs_contribution, equation_id, r_current_process_info);
                        //assemble the elemental contribution - here is where the ROM acts
                        //compute the elemental reduction matrix phi_elemental
                        const auto& r_geom = it->GetGeometry();
                        if(psi_elemental.size1() != cond_dofs.size() || psi_elemental.size2() != mPetrovGalerkinRomDofs)
                            psi_elemental.resize(cond_dofs.size(), mPetrovGalerkinRomDofs,false);
                        RomAuxiliaryUtilities::GetPsiElemental(psi_elemental, cond_dofs, r_geom, mMapPhi);
                        noalias(row(matrix_residuals, k+n_elements)) = prod(trans(psi_elemental), rhs_contribution); // The size of the residual will vary only when using more ROM modes, one row per condition
                    }
                }
            }
        return matrix_residuals;
        }

        Matrix GetProjectedGlobalLHSOntoPhi()
        {
            const auto& n_nodes = mrModelPart.NumberOfNodes();
            const int system_size = n_nodes*mNodalDofs;
            Matrix r_a = ZeroMatrix(system_size, mRomDofs);
            Vector r_b = ZeroVector(system_size);
            BuildGlobalLHS(mpScheme, mrModelPart, r_a, r_b);
            return r_a;
        }

        /**
         * Builds the reduced system of equations on rank 0 
         */
        void BuildGlobalLHS(
            BaseSchemeType::Pointer pScheme,
            ModelPart &rModelPart,
            Matrix &rA,
            Vector &rb) 
        {
            KRATOS_TRY

            // Build the system matrix by looping over elements and conditions and assembling to A
            KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

            // Get ProcessInfo from main model part
            const auto& r_current_process_info = rModelPart.GetProcessInfo();


            // Assemble all entities
            RomResidualsUtility::AssemblyTLS assembly_tls_container;

            const auto& r_elements = rModelPart.Elements();

            if(!r_elements.empty())
            {
                block_for_each(r_elements, assembly_tls_container, 
                    [&](Element& r_element, RomResidualsUtility::AssemblyTLS& r_thread_prealloc)
                {
                    CalculateLocalContributionLSPG(r_element, rA, rb, r_thread_prealloc, *pScheme, r_current_process_info);
                });
            }


            const auto& r_conditions = rModelPart.Conditions();

            if(!r_conditions.empty())
            {
                block_for_each(r_conditions, assembly_tls_container, 
                    [&](Condition& r_condition, RomResidualsUtility::AssemblyTLS& r_thread_prealloc)
                {
                    CalculateLocalContributionLSPG(r_condition, rA, rb, r_thread_prealloc, *pScheme, r_current_process_info);
                });
            }
            KRATOS_CATCH("")
        }

        template<typename TEntity>
        void CalculateLocalContributionLSPG(
            TEntity& rEntity,
            Matrix& rAglobal,
            Vector& rBglobal,
            RomResidualsUtility::AssemblyTLS& rPreAlloc,
            BaseSchemeType& rScheme,
            const ProcessInfo& rCurrentProcessInfo)
        {
            if (rEntity.IsDefined(ACTIVE) && rEntity.IsNot(ACTIVE)) return;

            // Calculate elemental contribution
            rScheme.CalculateSystemContributions(rEntity, rPreAlloc.lhs, rPreAlloc.romB, rPreAlloc.eq_id, rCurrentProcessInfo);
            rEntity.GetDofList(rPreAlloc.dofs, rCurrentProcessInfo);

            const std::size_t ndofs = rPreAlloc.dofs.size();
            ResizeIfNeeded(rPreAlloc.phiE, ndofs, mRomDofs);
            ResizeIfNeeded(rPreAlloc.romA, ndofs, mRomDofs);

            const auto &r_geom = rEntity.GetGeometry();
            RomAuxiliaryUtilities::GetPhiElemental(rPreAlloc.phiE, rPreAlloc.dofs, r_geom, mMapPhi);

            noalias(rPreAlloc.romA) = prod(rPreAlloc.lhs, rPreAlloc.phiE);


            // Assembly
            for(std::size_t row=0; row < ndofs; ++row)
            {
                const std::size_t global_row = rPreAlloc.eq_id[row];

                double& r_bi = rBglobal(global_row);
                AtomicAdd(r_bi, rPreAlloc.romB[row]);

                if(rPreAlloc.dofs[row]->IsFree())
                {
                    for(std::size_t col=0; col < mRomDofs; ++col)
                    {
                        // const std::size_t global_col = rPreAlloc.eq_id[col];
                        const std::size_t global_col = col;
                        double& r_aij = rAglobal(global_row, global_col);
                        AtomicAdd(r_aij, rPreAlloc.romA(row, col));
                    }
                }
            }
        }

    protected:
        std::vector< std::string > mNodalVariablesNames;
        int mNodalDofs;
        unsigned int mRomDofs;
        unsigned int mPetrovGalerkinRomDofs;
        ModelPart& mrModelPart;
        BaseSchemeType::Pointer mpScheme;
        std::unordered_map<Kratos::VariableData::KeyType, Matrix::size_type> mMapPhi;
    };


} // namespace Kratos



#endif // ROM_RESIDUALS_UTILITY_H_INCLUDED  defined