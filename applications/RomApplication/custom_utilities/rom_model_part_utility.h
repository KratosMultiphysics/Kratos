//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//  Kratos default license: kratos/license.txt
//
//  Main authors:    SEBASTIAN ARES DE PARGA REGALADO
//                   RAUL BRAVO
//

#if !defined( ROM_MODEL_PART_UTILITY_H_INCLUDED )
#define  ROM_MODEL_PART_UTILITY_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "spaces/ublas_space.h"
#include "utilities/math_utils.h"
#include "utilities/svd_utils.h"

/* Application includes */
#include "rom_application_variables.h"

namespace Kratos
{
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;

    // This utility returns the converged residuals projected onto the ROM basis Phi.
    class RomModelPartUtility
    {
        public:

        KRATOS_CLASS_POINTER_DEFINITION(RomModelPartUtility);

        RomModelPartUtility(
        ModelPart& rModelPart,
        Parameters ThisParameters,
        BaseSchemeType::Pointer pScheme
        ): mpModelPart(rModelPart), mpScheme(pScheme){
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "nodal_unknowns" : [],
            "number_of_dofs" : 10
        })" );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        mNodalVariablesNames = ThisParameters["nodal_unknowns"].GetStringArray();

        mNodalDofs = mNodalVariablesNames.size();
        mDofs = ThisParameters["number_of_dofs"].GetInt();

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

        ~RomModelPartUtility()= default;


        void GetPhiElemental(
            Matrix &PhiElemental,
            const Element::DofsVectorType &dofs,
            const Element::GeometryType &geom)
        {
            const auto *pcurrent_rom_nodal_basis = &(geom[0].GetValue(ROM_BASIS));
            int counter = 0;
            for(unsigned int k = 0; k < dofs.size(); ++k){
                auto variable_key = dofs[k]->GetVariable().Key();
                if(k==0)
                    pcurrent_rom_nodal_basis = &(geom[counter].GetValue(ROM_BASIS));
                else if(dofs[k]->Id() != dofs[k-1]->Id()){
                    counter++;
                    pcurrent_rom_nodal_basis = &(geom[counter].GetValue(ROM_BASIS));
                }
                if (dofs[k]->IsFixed())
                    noalias(row(PhiElemental, k)) = ZeroVector(PhiElemental.size2());
                else
                    noalias(row(PhiElemental, k)) = row(*pcurrent_rom_nodal_basis, MapPhi[variable_key]);
            }
        }


        Matrix Calculate()
        {
            // Getting the number of elements and conditions from the model
            const int nelements = static_cast<int>(mpModelPart.Elements().size());
            const int nconditions = static_cast<int>(mpModelPart.Conditions().size());

            const auto& CurrentProcessInfo = mpModelPart.GetProcessInfo();
            const auto el_begin = mpModelPart.ElementsBegin();
            const auto cond_begin = mpModelPart.ConditionsBegin();

            //contributions to the system
            Matrix LHS_Contribution = ZeroMatrix(0, 0);
            Vector RHS_Contribution = ZeroVector(0);

            //vector containing the localization in the system of the different terms
            Element::EquationIdVectorType EquationId;
            Matrix MatrixResiduals( (nelements + nconditions), mDofs); // Matrix of reduced residuals.
            Matrix PhiElemental;
            #pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, EquationId, PhiElemental, el_begin, cond_begin)
            {
                #pragma omp for nowait
                for (int k = 0; k < nelements; k++){
                    auto it_el = el_begin + k;
                    //detect if the element is active or not. If the user did not make any choice the element is active by default
                    bool element_is_active = true;
                    if ((it_el)->IsDefined(ACTIVE))
                        element_is_active = (it_el)->Is(ACTIVE);
                    if (element_is_active){
                        //calculate elemental contribution
                        mpScheme->CalculateSystemContributions(*it_el, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        Element::DofsVectorType dofs;
                        it_el->GetDofList(dofs, CurrentProcessInfo);
                        //assemble the elemental contribution - here is where the ROM acts
                        //compute the elemental reduction matrix PhiElemental
                        const auto& geom = it_el->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mDofs)
                            PhiElemental.resize(dofs.size(), mDofs,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        noalias(row(MatrixResiduals, k)) = prod(trans(PhiElemental), RHS_Contribution); // The size of the residual will vary only when using more ROM modes, one row per condition
                    }

                }

                #pragma omp for nowait
                for (int k = 0; k < nconditions;  k++){
                    ModelPart::ConditionsContainerType::iterator it = cond_begin + k;
                    //detect if the condition is active or not. If the user did not make any choice the condition is active by default
                    bool condition_is_active = true;
                    if ((it)->IsDefined(ACTIVE))
                        condition_is_active = (it)->Is(ACTIVE);
                    if (condition_is_active){
                        Condition::DofsVectorType dofs;
                        it->GetDofList(dofs, CurrentProcessInfo);
                        //calculate elemental contribution
                        mpScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        //assemble the elemental contribution - here is where the ROM acts
                        //compute the elemental reduction matrix PhiElemental
                        const auto& geom = it->GetGeometry();
                        if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mDofs)
                            PhiElemental.resize(dofs.size(), mDofs,false);
                        GetPhiElemental(PhiElemental, dofs, geom);
                        noalias(row(MatrixResiduals, k+nelements)) = prod(trans(PhiElemental), RHS_Contribution); // The size of the residual will vary only when using more ROM modes, one row per condition
                    }
                }
            }
        return MatrixResiduals;
        }

        Vector GetNonProjectedResiduals()
        {
            // Getting the number of elements and conditions from the model
            const int nelements = static_cast<int>(mpModelPart.Elements().size());
            const int nconditions = static_cast<int>(mpModelPart.Conditions().size());

            const auto& CurrentProcessInfo = mpModelPart.GetProcessInfo();
            const auto el_begin = mpModelPart.ElementsBegin();
            const auto cond_begin = mpModelPart.ConditionsBegin();

            int cond_node = cond_begin->GetGeometry().size();
            //vector containing the localization in the system of the different terms
            int TotalDofs = mNodalDofs*cond_node;

            //contributions to the system
            Matrix LHS_Contribution = ZeroMatrix(0, 0);
            Vector RHS_Contribution = ZeroVector(0);

            //vector containing the localization in the system of the different terms
            Element::EquationIdVectorType EquationId;
            Vector VectorResiduals((nelements + nconditions)*TotalDofs); // Matrix of reduced residuals.
            #pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, EquationId, el_begin, cond_begin)
            {
                #pragma omp for nowait
                for (int k = 0; k < nelements; k++){
                    auto it_el = el_begin + k;
                    //detect if the element is active or not. If the user did not make any choice the element is active by default
                    bool element_is_active = true;
                    if ((it_el)->IsDefined(ACTIVE))
                        element_is_active = (it_el)->Is(ACTIVE);
                    if (element_is_active){
                        //calculate elemental contribution
                        mpScheme->CalculateSystemContributions(*it_el, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        //assemble the elemental contribution - here is where the ROM acts
                        //compute the elemental reduction matrix PhiElemental
                        unsigned int local_size = RHS_Contribution.size();
                        unsigned int i_global = k*TotalDofs;
                        for (unsigned int i_local=0; i_local<local_size; i_local++){
                            VectorResiduals(i_global+i_local)=RHS_Contribution(i_local);
                        }
                    }

                }

                #pragma omp for nowait
                for (int k = 0; k < nconditions;  k++){
                    ModelPart::ConditionsContainerType::iterator it = cond_begin + k;
                    //detect if the condition is active or not. If the user did not make any choice the condition is active by default
                    bool condition_is_active = true;
                    if ((it)->IsDefined(ACTIVE))
                        condition_is_active = (it)->Is(ACTIVE);
                    if (condition_is_active){
                        //calculate elemental contribution
                        mpScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                        //assemble the elemental contribution - here is where the ROM acts
                        //compute the elemental reduction matrix PhiElemental
                        unsigned int local_size = RHS_Contribution.size();
                        unsigned int i_global = k*TotalDofs;
                        for (unsigned int i_local=0; i_local<local_size; i_local++){
                            VectorResiduals(i_global+i_local)=RHS_Contribution(i_local);
                        }
                    }
                }
            }
        return VectorResiduals;
        }

        Vector GetNonProjectedResidualsFromConditionList(Vector ConditionsList)
        {
            // Getting the number of elements and conditions from the model
            const int nconditions = static_cast<int>(ConditionsList.size());

            const auto& CurrentProcessInfo = mpModelPart.GetProcessInfo();
            const auto cond_begin = mpModelPart.ConditionsBegin();

            int cond_node = cond_begin->GetGeometry().size();
            //vector containing the localization in the system of the different terms
            int TotalDofs = mNodalDofs*cond_node;

            //contributions to the system
            Matrix LHS_Contribution = ZeroMatrix(0, 0);
            Vector RHS_Contribution = ZeroVector(0);

            //vector containing the localization in the system of the different terms
            Element::EquationIdVectorType EquationId;
            Vector VectorResiduals((nconditions)*TotalDofs); // Matrix of reduced residuals.
            #pragma omp parallel firstprivate( nconditions, LHS_Contribution, RHS_Contribution, EquationId)
            {
                #pragma omp for nowait
                for (int i=0; i<nconditions;i++) {
                    //detect if the element is active or not. If the user did not make any choice the element is active by default
                    //calculate elemental contribution
                    int cond_Id = ConditionsList(i);
                    mpScheme->CalculateSystemContributions(mpModelPart.GetCondition(cond_Id), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    //assemble the elemental contribution - here is where the ROM acts
                    //compute the elemental reduction matrix PhiElemental
                    unsigned int local_size = RHS_Contribution.size();
                    unsigned int i_global = i*TotalDofs;
                    for (unsigned int i_local=0; i_local<local_size; i_local++){
                        VectorResiduals(i_global+i_local)=RHS_Contribution(i_local);
                    }
                }
            }
        return VectorResiduals;
        }

        Vector GetNonProjectedResidualsFromElementList(Vector ElementList)
        {
            // Getting the number of elements and conditions from the model
            const int nelements = static_cast<int>(ElementList.size());

            const auto& CurrentProcessInfo = mpModelPart.GetProcessInfo();
            const auto elem_begin = mpModelPart.ElementsBegin();

            int elem_node = elem_begin->GetGeometry().size();
            //vector containing the localization in the system of the different terms
            int TotalDofs = mNodalDofs*elem_node;

            //contributions to the system
            Matrix LHS_Contribution = ZeroMatrix(0, 0);
            Vector RHS_Contribution = ZeroVector(0);

            //vector containing the localization in the system of the different terms
            Element::EquationIdVectorType EquationId;
            Vector VectorResiduals((nelements)*TotalDofs); // Matrix of reduced residuals.
            #pragma omp parallel firstprivate( nelements, LHS_Contribution, RHS_Contribution, EquationId)
            {
                #pragma omp for nowait
                for (int i=0; i<nelements;i++) {
                    //detect if the element is active or not. If the user did not make any choice the element is active by default
                    //calculate elemental contribution
                    int elem_Id = ElementList(i);
                    mpScheme->CalculateSystemContributions(mpModelPart.GetElement(elem_Id), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    //assemble the elemental contribution - here is where the ROM acts
                    //compute the elemental reduction matrix PhiElemental
                    unsigned int local_size = RHS_Contribution.size();
                    unsigned int i_global = i*TotalDofs;
                    for (unsigned int i_local=0; i_local<local_size; i_local++){
                        VectorResiduals(i_global+i_local)=RHS_Contribution(i_local);
                    }
                }
            }
        return VectorResiduals;
        }

        Vector GetConditionsList()
        {
            // Getting the number of elements and conditions from the model
            const int nconditions = static_cast<int>(mpModelPart.Conditions().size());

            const auto cond_begin = mpModelPart.ConditionsBegin();

            Vector VectorConditionList(nconditions); // Matrix of reduced residuals.
            #pragma omp parallel firstprivate(nconditions, cond_begin)
            {
                #pragma omp for nowait
                for (int k = 0; k < nconditions;  k++){
                    ModelPart::ConditionsContainerType::iterator it = cond_begin + k;
                    //detect if the condition is active or not. If the user did not make any choice the condition is active by default
                    bool condition_is_active = true;
                    if ((it)->IsDefined(ACTIVE))
                        condition_is_active = (it)->Is(ACTIVE);
                    if (condition_is_active){
                        VectorConditionList(k)=it->GetId();
                    }
                }
            }
        return VectorConditionList;
        }

        std::vector<IndexType> GetElementListFromNode(ModelPart& mrMainModelPart)
        {
            const int nNodes = static_cast<int>(mpModelPart.NumberOfNodes());

            std::vector<IndexType> ElementList;
            Vector NodeList(nNodes); 
            int i = 0;
            for (auto& node : mpModelPart.Nodes()){
                NodeList(i) = node.Id();
                i+=1;
            }
            for (auto& r_elem : mrMainModelPart.Elements()) {
  	            bool salir = 0;
                // auto& r_elem_Id = r_elem.Id();
                const auto& geom = r_elem.GetGeometry();
                for (auto& node : geom) {
                    // auto& node_Id = node.Id();
                    for (int k = 0; k < nNodes; k++){
                        if (NodeList(k)==node.Id()) {
                            ElementList.push_back(r_elem.Id());
                            salir = 1;
                            break;
                        }
                    }
                    if (salir){
                        break;
                    }
                }
            } 
        return ElementList;
        }



        // Matrix GetReactions_new(Matrix Snapshots_reactions, int scenarios)
        // {
        //     // Getting the number of elements and conditions from the model
        //     const int nconditions = static_cast<int>(mpModelPart.Conditions().size());

        //     const auto& CurrentProcessInfo = mpModelPart.GetProcessInfo();

        //     const auto cond_begin = mpModelPart.ConditionsBegin();
        //     int cond_node = cond_begin->GetGeometry().size();
        //     //vector containing the localization in the system of the different terms
        //     int TotalDofs = mNodalDofs*cond_node;
        //     Matrix MatrixReactions((nconditions*scenarios), TotalDofs); // Matrix of reduced residuals.
        //     Matrix PhiElemental;
        //     #pragma omp parallel firstprivate(nconditions, PhiElemental,scenarios)
        //     {
        //         #pragma omp for nowait
        //         for (int k = 0; k < nconditions;  k++){
        //             ModelPart::ConditionsContainerType::iterator it = cond_begin + k;
        //             //detect if the condition is active or not. If the user did not make any choice the condition is active by default
        //             bool condition_is_active = true;
        //             if ((it)->IsDefined(ACTIVE))
        //                 condition_is_active = (it)->Is(ACTIVE);
        //             if (condition_is_active){
        //                 Condition::DofsVectorType dofs;
        //                 it->GetDofList(dofs, CurrentProcessInfo);
        //                 //assemble the elemental contribution - here is where the ROM acts
        //                 //compute the elemental reduction matrix PhiElemental
        //                 const auto& geom = it->GetGeometry();
        //                 if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mDofs)
        //                     PhiElemental.resize(dofs.size(), mDofs,false);
        //                 GetPhiElemental(PhiElemental, dofs, geom);
        //                 Matrix PhiElemental_Inv (PhiElemental.size2(),PhiElemental.size1());
        //                 KRATOS_WATCH(PhiElemental)
        //                 PhiElemental_Inv = GetPseudoInverse(PhiElemental);
        //                 KRATOS_WATCH(PhiElemental_Inv)
        //                 Vector Snapshot_elemental_scenario( mDofs);
        //                 for (int i=0;i<scenarios;i++){
        //                     int begin = i*mDofs;
        //                     int RomDofs = mDofs;
        //                     for (int j=0;j<RomDofs;j++){
        //                         Snapshot_elemental_scenario(j) = Snapshots_reactions(k,begin+j);
        //                     }
        //                     noalias(row(MatrixReactions,i+(scenarios*k))) = trans(prod(PhiElemental_Inv,trans(Snapshot_elemental_scenario)));
        //                 }
        //             }
        //         }
        //     }
        // return MatrixReactions;
        // }

        // Matrix GetPseudoInverse(Matrix Matrix_to_invert)
        // {
        //     Matrix u_svd; // Orthogonal matrix (m x m)
        //     Matrix w_svd; // Rectangular diagonal matrix (m x n)
        //     Matrix v_svd; // Orthogonal matrix (n x n)
        //     std::string svd_type = "Jacobi"; // SVD decomposition type
        //     double svd_rel_tol = 1.0e-6; // Relative tolerance of the SVD decomposition (it will be multiplied by the input matrix norm)
        //     SVDUtils<double>::SingularValueDecomposition(Matrix_to_invert, u_svd, w_svd, v_svd, svd_type, svd_rel_tol);
        //     int data_rows = Matrix_to_invert.size1();
        //     int data_cols = Matrix_to_invert.size2();

        //     Matrix Matrix_pseudo_inverted (data_cols, data_rows);
        //     // Matrix aux (data_cols,data_cols);
        //     // Matrix aux_inv (data_cols,data_cols);
        //     // aux = prod(trans(Matrix_to_invert),Matrix_to_invert);
        //     // double det;
        //     // MathUtils<double>::InvertMatrix(aux,aux_inv,det);
        //     // Matrix_pseudo_inverted = prod(aux_inv,trans(Matrix_to_invert));
        //     Matrix_pseudo_inverted = prod(trans(v_svd),prod(w_svd,u_svd));
        // return Matrix_pseudo_inverted;  
        // }

        // Matrix GetDofArray()
        // {
        //     // Getting the number of elements and conditions from the model
        //     const int nconditions = static_cast<int>(mpModelPart.Conditions().size());

        //     const auto& CurrentProcessInfo = mpModelPart.GetProcessInfo();

        //     const auto cond_begin = mpModelPart.ConditionsBegin();
        //     int cond_node = cond_begin->GetGeometry().size();
        //     //vector containing the localization in the system of the different terms
        //     int TotalDofs = mNodalDofs*cond_node;
        //     Matrix dof_list(nconditions, TotalDofs); // Matrix of reduced residuals.
        //     #pragma omp parallel firstprivate(nconditions)
        //     {
        //         #pragma omp for nowait
        //         for (int k = 0; k < nconditions;  k++){
        //             ModelPart::ConditionsContainerType::iterator it = cond_begin + k;
        //             //detect if the condition is active or not. If the user did not make any choice the condition is active by default
        //             bool condition_is_active = true;
        //             if ((it)->IsDefined(ACTIVE))
        //                 condition_is_active = (it)->Is(ACTIVE);
        //             if (condition_is_active){
        //                 Condition::DofsVectorType dofs;
        //                 it->GetDofList(dofs, CurrentProcessInfo);
        //                 int counter = 0;
        //                 for (auto& p_dof : dofs){
        //                     const IndexType dof_id = p_dof->Id();
        //                     dof_list(k,counter) = dof_id;
        //                     counter += 1;
        //                 }
        //             }
        //         }
        //     }
        // return dof_list;
        // }
        

        protected:
            std::vector< std::string > mNodalVariablesNames;
            int mNodalDofs;
            unsigned int mDofs;
            BaseSchemeType::Pointer mpScheme;
            ModelPart& mpModelPart;
            std::unordered_map<Kratos::VariableData::KeyType,int> MapPhi;
        };



} // namespace Kratos



#endif // ROM_MODEL_PART_UTILITY_H_INCLUDED  defined