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
//            
//

#if !defined( ROM_MODEL_PART_UTILITY_H_INCLUDED )
#define  ROM_MODEL_PART_UTILITY_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "spaces/ublas_space.h"/*
#include "utilities/math_utils.h"
#include "utilities/sparse_matrix_multiplication_utility.h"
#include "utilities/svd_utils.h"*/

/* Application includes */
#include "rom_application_variables.h"


namespace Kratos
{
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;
    typedef ModelPart::DofType DofType;
    typedef ModelPart::DofsArrayType DofsArrayType;
    typedef typename std::unordered_set<DofType::Pointer, DofPointerHasher> DofSetType;
    typedef ModelPart::DofsVectorType DofsVectorType;

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

        // This function returns the residuals within the Initializing Model Part as a column vector.
        Vector GetNonProjectedResiduals()
        {
            // Getting the number of elements and conditions from the model
            const int nelements = static_cast<int>(mpModelPart.Elements().size());
            const int nconditions = static_cast<int>(mpModelPart.Conditions().size());

            const auto& CurrentProcessInfo = mpModelPart.GetProcessInfo();
            const auto el_begin = mpModelPart.ElementsBegin();
            const auto cond_begin = mpModelPart.ConditionsBegin();

            int cond_node = cond_begin->GetGeometry().size();
            int TotalDofs = mNodalDofs*cond_node;

            //contributions to the system
            Matrix LHS_Contribution = ZeroMatrix(0, 0);
            Vector RHS_Contribution = ZeroVector(0);

            //vector containing the localization in the system of the different terms
            Element::EquationIdVectorType EquationId;
            Vector VectorResiduals((nelements + nconditions)*TotalDofs); // Vector of residuals.
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

        // Returns a column vector containing the residuals of a given condition list (vector), the residuals vector contains all the DoFs residual values within all the condition list
        Vector GetNonProjectedResidualsFromConditionList(Vector ConditionsList)
        {
            // Getting the number of conditions from the condition list
            const int nconditions = static_cast<int>(ConditionsList.size());

            const auto& CurrentProcessInfo = mpModelPart.GetProcessInfo();
            const auto cond_begin = mpModelPart.ConditionsBegin();

            int cond_node = cond_begin->GetGeometry().size();
            int TotalDofs = mNodalDofs*cond_node;

            //contributions to the system
            Matrix LHS_Contribution = ZeroMatrix(0, 0);
            Vector RHS_Contribution = ZeroVector(0);

            Element::EquationIdVectorType EquationId;
            Vector VectorResiduals((nconditions)*TotalDofs); // Vector of residuals.
            #pragma omp parallel firstprivate( nconditions, LHS_Contribution, RHS_Contribution, EquationId)
            {
                #pragma omp for nowait
                for (int i=0; i<nconditions;i++) {
                    //detect if the element is active or not. If the user did not make any choice the element is active by default
                    //calculate elemental contribution
                    int cond_Id = ConditionsList(i);
                    mpScheme->CalculateSystemContributions(mpModelPart.GetCondition(cond_Id), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    unsigned int local_size = RHS_Contribution.size();
                    unsigned int i_global = i*TotalDofs;
                    for (unsigned int i_local=0; i_local<local_size; i_local++){
                        VectorResiduals(i_global+i_local)=RHS_Contribution(i_local);
                    }
                }
            }
        return VectorResiduals;
        }

        // Returns a column vector containing the residuals of a given element list (vector), the residuals vector contains all the DoFs residual values within all the elements list
        Vector GetNonProjectedResidualsFromElementList(Vector ElementList)
        {
            // Getting the number of elements from the element list
            const int nelements = static_cast<int>(ElementList.size());

            const auto& CurrentProcessInfo = mpModelPart.GetProcessInfo();
            const auto elem_begin = mpModelPart.ElementsBegin();

            int elem_node = elem_begin->GetGeometry().size();
            int TotalDofs = mNodalDofs*elem_node;

            //contributions to the system
            Matrix LHS_Contribution = ZeroMatrix(0, 0);
            Vector RHS_Contribution = ZeroVector(0);

            Element::EquationIdVectorType EquationId;
            Vector VectorResiduals((nelements)*TotalDofs); // Vector of residuals.
            #pragma omp parallel firstprivate( ElementList, nelements, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo)
            {
                #pragma omp for nowait
                for (int i=0; i<nelements;i++) {
                    //detect if the element is active or not. If the user did not make any choice the element is active by default
                    //calculate elemental contribution
                    int elem_Id = ElementList(i);
                    mpScheme->CalculateSystemContributions(mpModelPart.GetElement(elem_Id), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    unsigned int local_size = RHS_Contribution.size();
                    unsigned int i_global = i*TotalDofs;
                    for (unsigned int i_local=0; i_local<local_size; i_local++){
                        VectorResiduals(i_global+i_local)=RHS_Contribution(i_local); 
                    }
                }
            }
        return VectorResiduals;
        }

        // This function returns the Id's of the condition within the Initializing Model Part
        Vector GetConditionsList()
        {
            // Getting the number of conditions from the model part
            const int nconditions = static_cast<int>(mpModelPart.Conditions().size());

            const auto cond_begin = mpModelPart.ConditionsBegin();

            Vector VectorConditionList(nconditions); 
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

        // This function returns the list of elements (within an input Model Part) that are connected to the nodes conteined in the Initializing Model Part of the utility.
        std::vector<IndexType> GetElementListFromNode(ModelPart& mrMainModelPart)
        {
            const int nNodes = static_cast<int>(mpModelPart.NumberOfNodes()); // Total nodes of the Initializing Model Part

            std::vector<IndexType> ElementList; // Initialize the output Elements list
            Vector NodeList(nNodes); // Initialize a Node list 
            int i = 0;
            // Loop on Initializing Model Part Nodes to create list of nodes.
            for (auto& node : mpModelPart.Nodes()){
                NodeList(i) = node.Id();
                i+=1;
            }
            // Loop on input Model Part Elements
            for (auto& r_elem : mrMainModelPart.Elements()) {
  	            bool salir = 0;
                const auto& geom = r_elem.GetGeometry(); // Get elemental nodes
                // Loop on elemental nodes
                for (auto& node : geom) {
                    // Loop on Initializng Model Part Node list
                    for (int k = 0; k < nNodes; k++){
                        if (NodeList(k)==node.Id()) { //Check if one node coincide
                            ElementList.push_back(r_elem.Id()); // Add element Id to Element List if one node coincide
                            salir = 1; // Flag to break out of Nodal's loop. (Jumps to the next Element whenever a node coincide)
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

        std::vector<IndexType> GetNodeList()
        {
            std::vector<IndexType> NodeList; // Initialize the output Elements list
            // Loop on Initializing Model Part Nodes to create list of nodes.
            for (auto& node : mpModelPart.Nodes()){
                NodeList.push_back(node.Id());
            }
        return NodeList;
        }

        Matrix GetDofsFromElementList(Vector ElementList)
        {
            // Getting the number of elements and conditions from the model
            const int nelements = static_cast<int>(ElementList.size());

            const auto& CurrentProcessInfo = mpModelPart.GetProcessInfo();

            const auto elem_begin = mpModelPart.ElementsBegin();
            int elem_node = elem_begin->GetGeometry().size();
            
            //vector containing the localization in the system of the different terms
            int TotalDofs = mNodalDofs*elem_node;
            Matrix dof_list(nelements, TotalDofs); // Matrix of reduced residuals.
            #pragma omp parallel firstprivate(nelements)
            {
                #pragma omp for nowait
                for (int i=0; i<nelements;i++) {
                    
                    //detect if the element is active or not. If the user did not make any choice the element is active by default
                    //calculate elemental contribution
                    int elem_Id = ElementList(i);
                    Element::Pointer i_elem = mpModelPart.pGetElement(elem_Id);
                    Element::DofsVectorType dofs;
                    i_elem->GetDofList(dofs, CurrentProcessInfo);
                    KRATOS_WATCH(dofs.size())
                    int counter = 0;
                    for (auto& p_dof : dofs){
                        const IndexType dof_id = p_dof->GetId();
                        KRATOS_WATCH(dof_id)
                        dof_list(i,counter) = dof_id;
                        counter += 1;
                    }
                }
            }
        return dof_list;
        }

        // This function assembles the reactions in the restricted model part with by assembling the residual obtained in the process reconstruction of reactions.
        void AssembleReactions(Vector ResidualReconstruction, Vector RestrictedResidualNodes, std::vector<std::vector<IndexType>> ListNodes, std::vector<std::vector<IndexType>> ListElements, int DofElem, int DofNode)
        {
            int nnodes = RestrictedResidualNodes.size();
            #pragma omp parallel firstprivate(nnodes)
            {
                #pragma omp for nowait
                // Loop in nodes
                for (int i=0;i<nnodes;i++){
                    auto& node = mpModelPart.GetNode(RestrictedResidualNodes(i));
                    double reaction_x = 0;
                    double reaction_y = 0;
                    double reaction_z = 0;
                    // List of Nodes is a list of lists of the neighbouring elements of each node
                    // Get the elements list from the node
                    std::vector<IndexType> current_node_list = ListNodes[i];
                    for (IndexType j : current_node_list){
                        int counter = 0;
                        // List of Nodes is a list of lists of the neighbouring nodes of each element.
                        // Loop in all nodes inside each neighboring element to check the corresponding index of contribution and its value.
                        std::vector<IndexType> current_element_list = ListElements[j];
                        for (IndexType k : current_element_list){
                            if (node.Id()==k){
                                int aux_index = (j*DofElem)+(counter*DofNode);
                                reaction_x += ResidualReconstruction(aux_index);
                                reaction_y += ResidualReconstruction(aux_index+1);
                                reaction_z += ResidualReconstruction(aux_index+2);
                            }
                            counter+=1;
                        }
                    }
                    auto& v = node.FastGetSolutionStepValue(REACTION);
                    v[0] = reaction_x;
                    v[1] = reaction_y;
                    v[2] = reaction_z;
                }
            }
        }

        void AssembleStresses(Matrix StressMatrix){
            ModelPart& hrom_model_part = mpModelPart.GetSubModelPart("VISUALIZE_HROM");
            const int nnodes = static_cast<int>(hrom_model_part.NumberOfNodes());
            #pragma omp parallel for
                for (int i_node = 0; i_node < nnodes; ++i_node)
                {
                    auto it_node = hrom_model_part.NodesBegin() + i_node;
                    it_node->SetValue(PK2_STRESS_VECTOR,row(StressMatrix,it_node->Id()-1));
                // node.SetValue(PK2_STRESS_VECTOR,row(StressMatrix,node.Id()-1));
                }     
        }

        // SparseSpaceType::MatrixType BuildExtrapolationOperator(){
        //     // Build a sparse extrapolation operator for stresses (from GP to nodes)
        //     const int nnodes = static_cast<int>(mpModelPart.NumberOfNodes());
        //     const int nelements = static_cast<int>(mpModelPart.NumberOfElements());
        //     const auto elem_begin = mpModelPart.ElementsBegin();

        //     SparseSpaceType::MatrixType P;
        //     P.resize(nnodes,nelements,false);
        //     SparseSpaceType::SetToZero(P);
        //     Vector P_diag (nnodes);
        //     #pragma omp parallel firstprivate(nelements)
        //     {
        //         #pragma omp for nowait
        //         for (int k=0; k<nelements;k++) {
        //             Element::Pointer i_elem = mpModelPart.pGetElement(k);
        //             const auto& elem_geom = i_elem->pGetGeometry();
        //             const Matrix shape_functions_evaluated_in_gp = elem_geom->ShapeFunctionsValues();
        //             double elem_area = elem_geom.Area();
        //             Vector sum_shape_functions_evaluated_in_gp = LumpMatrix(shape_functions_evaluated_in_gp);
        //             for (auto& node : elem_geom){
        //                 P_diag(node.Id()-1) += elem_area;
        //             }
        //             for (int i ; i< shape_functions_evaluated_in_gp.size1(),i++){
        //                 for (int j; j < shape_functions_evaluated_in_gp.size2(),j++){
        //                     P((j).Id()-1,elem.Id-1] = shape_functions_evaluated_in_gp[i,j]*elem_area/sum_shape_functions_evaluated_in_gp[j]
        //                 } 
        //             }
                        
        //         P_diag = np.diag(1/P_diag)
        //         P_diag = lil_matrix(P_diag)
        //         P = P_diag@P
        //         save_npz("shape.npz", P)
        //         }
        // }

        // Vector LumpMatrix(Matrix MatrixToBeLumped){
        //     const auto size_1 = MatrixToBeLumped.size1();
        //     const auto size_2 = MatrixToBeLumped.size2();
        //     Vector LumpedMatrixToVector(size_1,0);
        //     #pragma omp parallel for
        //         for (int i = 0; i < size_1; i++){
        //             for (int j = 0; j < size_2; j++){
        //                 LumpedMatrixToVector(i)+=MatrixToBeLumped(i,j);
        //             }
        //         }
        // return LumpedMatrixToVector;
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
;


#endif // ROM_MODEL_PART_UTILITY_H_INCLUDED  defined
