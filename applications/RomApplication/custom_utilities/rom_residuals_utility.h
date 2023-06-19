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
                MapPhi[var.Key()] = k;
            }
            else
                KRATOS_ERROR << "variable \""<< mNodalVariablesNames[k] << "\" not valid" << std::endl;

            }
        }

        ~RomResidualsUtility()= default;

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
                        RomAuxiliaryUtilities::GetPhiElemental(phi_elemental, elem_dofs, r_geom, MapPhi);
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
                        RomAuxiliaryUtilities::GetPhiElemental(phi_elemental, cond_dofs, r_geom, MapPhi);
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
                        RomAuxiliaryUtilities::GetPsiElemental(psi_elemental, elem_dofs, r_geom, MapPhi);
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
                        RomAuxiliaryUtilities::GetPsiElemental(psi_elemental, cond_dofs, r_geom, MapPhi);
                        noalias(row(matrix_residuals, k+n_elements)) = prod(trans(psi_elemental), rhs_contribution); // The size of the residual will vary only when using more ROM modes, one row per condition
                    }
                }
            }
        return matrix_residuals;
        }

        Matrix GetProjectedGlobalLHS()
        {
            const int n_elements = static_cast<int>(mrModelPart.Elements().size());
            const int n_conditions = static_cast<int>(mrModelPart.Conditions().size());
            const auto& n_nodes = mrModelPart.NumberOfNodes();

            const auto& r_current_process_info = mrModelPart.GetProcessInfo();
            
            const int system_size = n_nodes*mNodalDofs;

            const auto el_begin = mrModelPart.ElementsBegin();
            const auto cond_begin = mrModelPart.ConditionsBegin();

            //contributions to the system
            Matrix lhs_contribution = ZeroMatrix(0,0);

            //vector containing the localization in the system of the different terms
            Element::EquationIdVectorType equation_id;
            Matrix a_phi = ZeroMatrix(system_size, mRomDofs);

            //dofs container initialization
            Element::DofsVectorType elem_dofs;
            Condition::DofsVectorType cond_dofs;

            Matrix phi_elemental;
            Matrix temp_a_phi = ZeroMatrix(system_size,mRomDofs);
            Matrix aux;

            #pragma omp parallel firstprivate(n_elements, n_conditions, lhs_contribution, equation_id, el_begin, cond_begin, elem_dofs, cond_dofs)
            {

                #pragma omp for nowait
                for (int k = 0; k < static_cast<int>(n_elements); k++) {
                    const auto it_el = el_begin + k;

                    // Detect if the element is active or not. If the user did not make any choice the element is active by default
                    const bool element_is_active = it_el->IsDefined(ACTIVE) ? it_el->Is(ACTIVE) : true;

                    // Calculate elemental contribution
                    if (element_is_active){
                        mpScheme->CalculateLHSContribution(*it_el, lhs_contribution, equation_id, r_current_process_info);
                        it_el->GetDofList(elem_dofs, r_current_process_info);
                        const auto &r_geom = it_el->GetGeometry();
                        if(phi_elemental.size1() != elem_dofs.size() || phi_elemental.size2() != mRomDofs) {
                            phi_elemental.resize(elem_dofs.size(), mRomDofs,false);
                        }
                        if(aux.size1() != elem_dofs.size() || aux.size2() != mRomDofs) {
                            aux.resize(elem_dofs.size(), mRomDofs,false);
                        }
                        RomAuxiliaryUtilities::GetPhiElemental(phi_elemental, elem_dofs, r_geom, MapPhi);
                        noalias(aux) = prod(lhs_contribution, phi_elemental);
                        for(int d = 0; d < static_cast<int>(elem_dofs.size()); ++d){
                            if (elem_dofs[d]->IsFixed()==false){
                                row(temp_a_phi,elem_dofs[d]->EquationId()) += row(aux,d);// Add contributions to global system for free dofs.
                            }
                        }
                    }
                }

                #pragma omp for nowait
                for (int k = 0; k < static_cast<int>(n_conditions); k++){
                    const auto it = cond_begin + k;

                    // Detect if the element is active or not. If the user did not make any choice the condition is active by default
                    const bool condition_is_active = it->IsDefined(ACTIVE) ? it->Is(ACTIVE) : true;

                    // Calculate condition contribution
                    if (condition_is_active) {
                        it->GetDofList(cond_dofs, r_current_process_info);
                        mpScheme->CalculateLHSContribution(*it, lhs_contribution, equation_id, r_current_process_info);
                        const auto &r_geom = it->GetGeometry();
                        if(phi_elemental.size1() != cond_dofs.size() || phi_elemental.size2() != mRomDofs) {
                            phi_elemental.resize(cond_dofs.size(), mRomDofs,false);
                        }
                        if(aux.size1() != cond_dofs.size() || aux.size2() != mRomDofs) {
                            aux.resize(cond_dofs.size(), mRomDofs,false);
                        }
                        RomAuxiliaryUtilities::GetPhiElemental(phi_elemental, cond_dofs, r_geom, MapPhi);
                        noalias(aux) = prod(lhs_contribution, phi_elemental);
                        for(int d = 0; d < static_cast<int>(cond_dofs.size()); ++d){
                            if (cond_dofs[d]->IsFixed()==false){
                                row(temp_a_phi,cond_dofs[d]->EquationId()) += row(aux,d);
                            }
                        }
                    }
                }

                #pragma omp critical
                {
                    noalias(a_phi) += temp_a_phi;
                }

            }
            return a_phi;
        }

    protected:
        std::vector< std::string > mNodalVariablesNames;
        int mNodalDofs;
        unsigned int mRomDofs;
        unsigned int mPetrovGalerkinRomDofs;
        ModelPart& mrModelPart;
        BaseSchemeType::Pointer mpScheme;
        std::unordered_map<Kratos::VariableData::KeyType, Matrix::size_type> MapPhi;
    };


} // namespace Kratos



#endif // ROM_RESIDUALS_UTILITY_H_INCLUDED  defined