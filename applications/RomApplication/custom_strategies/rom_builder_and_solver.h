//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Riccardo Rossi
//                  Raul Bravo
//
#if !defined(KRATOS_ROM_BUILDER_AND_SOLVER)
#define KRATOS_ROM_BUILDER_AND_SOLVER

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "utilities/builtin_timer.h"

/* Application includes */
#include "rom_application_variables.h"

namespace Kratos
{

template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class ROMBuilderAndSolver : public BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    /**
     * This struct is used in the component wise calculation only
     * is defined here and is used to declare a member variable in the component wise builder and solver
     * private pointers can only be accessed by means of set and get functions
     * this allows to set and not copy the Element_Variables and Condition_Variables
     * which will be asked and set by another strategy object
     */

    //pointer definition

    KRATOS_CLASS_POINTER_DEFINITION(ROMBuilderAndSolver);

    // The size_t types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// Definition of the classes from the base class
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::TDataType TDataType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    /// Additional definitions
    typedef PointerVectorSet<Element, IndexedObject> ElementsContainerType;
    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;
    typedef boost::numeric::ublas::compressed_matrix<double> CompressedMatrixType;

    /// DoF types definition
    typedef Node<3> NodeType;
    typedef typename NodeType::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /**
     * @brief Default constructor. (with parameters)
     */
    explicit ROMBuilderAndSolver(typename TLinearSolver::Pointer pNewLinearSystemSolver, Parameters ThisParameters)
        : BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSystemSolver)
    {
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "nodal_unknowns" : [],
            "number_of_rom_dofs" : 10
        })");

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // We set the other member variables
        mpLinearSystemSolver = pNewLinearSystemSolver;

        mNodalVariablesNames = ThisParameters["nodal_unknowns"].GetStringArray();

        mNodalDofs = mNodalVariablesNames.size();
        mRomDofs = ThisParameters["number_of_rom_dofs"].GetInt();

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

    /** Destructor.
     */
    ~ROMBuilderAndSolver() = default;

    virtual void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart) override
    {
        KRATOS_TRY;

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Setting up the dofs" << std::endl;

        //Gets the array of elements from the modeler
        auto &r_elements_array = rModelPart.Elements();
        const int number_of_elements = static_cast<int>(r_elements_array.size());

        DofsVectorType dof_list, second_dof_list; // NOTE: The second dof list is only used on constraints to include master/slave relations

        unsigned int nthreads = ParallelUtilities::GetNumThreads();

        typedef std::unordered_set<NodeType::DofType::Pointer, DofPointerHasher> set_type;

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of threads" << nthreads << "\n" << std::endl;

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Initializing element loop" << std::endl;

        /**
         * Here we declare three sets.
         * - The global set: Contains all the DoF of the system
         * - The slave set: The DoF that are not going to be solved, due to MPC formulation
         */
        set_type dof_global_set;
        dof_global_set.reserve(number_of_elements * 20);

        if (mHromWeightsInitialized == false){
            int number_of_hrom_elements=0;
            #pragma omp parallel firstprivate(dof_list, second_dof_list) reduction(+:number_of_hrom_elements)
            {
                const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

                // We create the temporal set and we reserve some space on them
                set_type dofs_tmp_set;
                dofs_tmp_set.reserve(20000);
                // Gets the array of elements from the modeler
                ModelPart::ElementsContainerType selected_elements_private;
                #pragma omp for schedule(guided, 512) nowait
                for (int i = 0; i < number_of_elements; ++i)
                {
                    auto it_elem = r_elements_array.begin() + i;
                    //detect whether the element has a Hyperreduced Weight (H-ROM simulation) or not (ROM simulation)
                    if ((it_elem)->Has(HROM_WEIGHT)){
                        selected_elements_private.push_back(*it_elem.base());
                        number_of_hrom_elements++;
                    }
                    else
                        it_elem->SetValue(HROM_WEIGHT, 1.0);
                    // Gets list of Dof involved on every element
                    pScheme->GetDofList(*it_elem, dof_list, r_current_process_info);
                    dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                }

                // Gets the array of conditions from the modeler
                ConditionsArrayType &r_conditions_array = rModelPart.Conditions();
                const int number_of_conditions = static_cast<int>(r_conditions_array.size());

                ModelPart::ConditionsContainerType selected_conditions_private;
                #pragma omp for schedule(guided, 512) nowait
                for (int i = 0; i < number_of_conditions; ++i)
                {
                    auto it_cond = r_conditions_array.begin() + i;
                    // Gather the H-reduced conditions that are to be considered for assembling. Ignoring those for displaying results only
                    if (it_cond->Has(HROM_WEIGHT)){
                        selected_conditions_private.push_back(*it_cond.base());
                        number_of_hrom_elements++;
                    }
                    else
                        it_cond->SetValue(HROM_WEIGHT, 1.0);
                    // Gets list of Dof involved on every element
                    pScheme->GetDofList(*it_cond, dof_list, r_current_process_info);
                    dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                }
                #pragma omp critical
                {
                    for (auto &cond : selected_conditions_private){
                        mSelectedConditions.push_back(&cond);
                    }
                    for (auto &elem : selected_elements_private){
                        mSelectedElements.push_back(&elem);
                    }

                }

                // Gets the array of constraints from the modeler
                auto &r_constraints_array = rModelPart.MasterSlaveConstraints();
                const int number_of_constraints = static_cast<int>(r_constraints_array.size());
                #pragma omp for schedule(guided, 512) nowait
                for (int i = 0; i < number_of_constraints; ++i)
                {
                    auto it_const = r_constraints_array.begin() + i;

                    // Gets list of Dof involved on every element
                    it_const->GetDofList(dof_list, second_dof_list, r_current_process_info);
                    dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                    dofs_tmp_set.insert(second_dof_list.begin(), second_dof_list.end());
                }

                // We merge all the sets in one thread
                #pragma omp critical
                {
                    dof_global_set.insert(dofs_tmp_set.begin(), dofs_tmp_set.end());
                }
            }
            if (number_of_hrom_elements>0){
                mHromSimulation = true;
            }
            mHromWeightsInitialized = true;
        }
        else{
            #pragma omp parallel firstprivate(dof_list, second_dof_list)
            {
                const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

                // We cleate the temporal set and we reserve some space on them
                set_type dofs_tmp_set;
                dofs_tmp_set.reserve(20000);

                // Gets the array of elements from the modeler
                #pragma omp for schedule(guided, 512) nowait
                for (int i = 0; i < number_of_elements; ++i) {
                    auto it_elem = r_elements_array.begin() + i;

                    // Gets list of Dof involved on every element
                    pScheme->GetDofList(*it_elem, dof_list, r_current_process_info);
                    dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                }

                // Gets the array of conditions from the modeler
                ConditionsArrayType& r_conditions_array = rModelPart.Conditions();
                const int number_of_conditions = static_cast<int>(r_conditions_array.size());
                #pragma omp for  schedule(guided, 512) nowait
                for (int i = 0; i < number_of_conditions; ++i) {
                    auto it_cond = r_conditions_array.begin() + i;

                    // Gets list of Dof involved on every element
                    pScheme->GetDofList(*it_cond, dof_list, r_current_process_info);
                    dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                }

                // Gets the array of constraints from the modeler
                auto& r_constraints_array = rModelPart.MasterSlaveConstraints();
                const int number_of_constraints = static_cast<int>(r_constraints_array.size());
                #pragma omp for  schedule(guided, 512) nowait
                for (int i = 0; i < number_of_constraints; ++i) {
                    auto it_const = r_constraints_array.begin() + i;

                    // Gets list of Dof involved on every element
                    it_const->GetDofList(dof_list, second_dof_list, r_current_process_info);
                    dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                    dofs_tmp_set.insert(second_dof_list.begin(), second_dof_list.end());
                }

                // We merge all the sets in one thread
                #pragma omp critical
                {
                    dof_global_set.insert(dofs_tmp_set.begin(), dofs_tmp_set.end());
                }
            }
        }

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Initializing ordered array filling\n" << std::endl;

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();

        Doftemp.reserve(dof_global_set.size());
        for (auto it = dof_global_set.begin(); it != dof_global_set.end(); it++)
        {
            Doftemp.push_back(*it);
        }
        Doftemp.Sort();

        BaseType::mDofSet = Doftemp;

        //Throws an exception if there are no Degrees Of Freedom involved in the analysis
        KRATOS_ERROR_IF(BaseType::mDofSet.size() == 0) << "No degrees of freedom!" << std::endl;

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "Number of degrees of freedom:" << BaseType::mDofSet.size() << std::endl;

        BaseType::mDofSetIsInitialized = true;

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished setting up the dofs" << std::endl;

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2)) << "End of setup dof set\n"
                                                                          << std::endl;

#ifdef KRATOS_DEBUG
        // If reactions are to be calculated, we check if all the dofs have reactions defined
        // This is tobe done only in debug mode
        if (BaseType::GetCalculateReactionsFlag())
        {
            for (auto dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
            {
                KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction()) << "Reaction variable not set for the following : " << std::endl
                                                                 << "Node : " << dof_iterator->Id() << std::endl
                                                                 << "Dof : " << (*dof_iterator) << std::endl
                                                                 << "Not possible to calculate reactions." << std::endl;
            }
        }
#endif
        KRATOS_CATCH("");
    }

    /**
            organises the dofset in order to speed up the building phase
     */
    virtual void SetUpSystem(
        ModelPart &r_model_part
    ) override
    {
        //int free_id = 0;
        BaseType::mEquationSystemSize = BaseType::mDofSet.size();
        int ndofs = static_cast<int>(BaseType::mDofSet.size());

        #pragma omp parallel for firstprivate(ndofs)
        for (int i = 0; i < static_cast<int>(ndofs); i++){
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + i;
            dof_iterator->SetEquationId(i);
        }
    }


    // Vector ProjectToReducedBasis(
	// 	const TSystemVectorType& rX,
	// 	ModelPart::NodesContainerType& rNodes
	// )
    // {
    //     Vector rom_unknowns = ZeroVector(mRomDofs);
    //     for(const auto& node : rNodes)
    //     {
    //         unsigned int node_aux_id = node.GetValue(AUX_ID);
    //         const auto& nodal_rom_basis = node.GetValue(ROM_BASIS);
	// 			for (int i = 0; i < mRomDofs; ++i) {
	// 				for (int j = 0; j < mNodalDofs; ++j) {
	// 					rom_unknowns[i] += nodal_rom_basis(j, i)*rX(node_aux_id*mNodalDofs + j);
	// 				}
	// 			}
    //     }
    //     return rom_unknowns;
	// }

    void ProjectToFineBasis(
        const TSystemVectorType &rRomUnkowns,
        ModelPart &rModelPart,
        TSystemVectorType &Dx)
    {
        const auto dofs_begin = BaseType::mDofSet.begin();
        const auto dofs_number = BaseType::mDofSet.size();

        #pragma omp parallel firstprivate(dofs_begin, dofs_number)
        {
            const Matrix *pcurrent_rom_nodal_basis = nullptr;
            unsigned int old_dof_id;
            #pragma omp for nowait
            for (int k = 0; k < static_cast<int>(dofs_number); k++){
                auto dof = dofs_begin + k;
                if(pcurrent_rom_nodal_basis == nullptr){
                    pcurrent_rom_nodal_basis = &(rModelPart.pGetNode(dof->Id())->GetValue(ROM_BASIS));
                    old_dof_id = dof->Id();
                }
                else if(dof->Id() != old_dof_id ){
                    pcurrent_rom_nodal_basis = &(rModelPart.pGetNode(dof->Id())->GetValue(ROM_BASIS));
                    old_dof_id = dof->Id();
                }
                Dx[dof->EquationId()] = inner_prod(  row(  *pcurrent_rom_nodal_basis    , mMapPhi[dof->GetVariable().Key()]   )     , rRomUnkowns);
            }
        }
    }

    void GetPhiElemental(
        Matrix &PhiElemental,
        const Element::DofsVectorType &dofs,
        const Element::GeometryType &geom)
    {
        const Matrix *pcurrent_rom_nodal_basis = nullptr;
        int counter = 0;
        for(int k = 0; k < static_cast<int>(dofs.size()); ++k){
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
                noalias(row(PhiElemental, k)) = row(*pcurrent_rom_nodal_basis, mMapPhi[variable_key]);
        }
    }



    /*@{ */

    /**
            Function to perform the building and solving phase at the same time.
            It is ideally the fastest and safer function to use when it is possible to solve
            just after building
     */
    virtual void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        //define a dense matrix to hold the reduced problem
        Matrix Arom = ZeroMatrix(mRomDofs, mRomDofs);
        Vector brom = ZeroVector(mRomDofs);
        TSystemVectorType x(Dx.size());

        const auto forward_projection_timer = BuiltinTimer();
        Vector xrom = ZeroVector(mRomDofs);
        //this->ProjectToReducedBasis(x, rModelPart.Nodes(),xrom);
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Project to reduced basis time: " << forward_projection_timer.ElapsedSeconds() << std::endl;

        //build the system matrix by looping over elements and conditions and assembling to A
        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        auto help_nelements = static_cast<int>(rModelPart.Elements().size());
        auto help_el_begin = rModelPart.ElementsBegin();

        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        auto help_cond_begin = rModelPart.ConditionsBegin();
        auto help_nconditions = static_cast<int>(rModelPart.Conditions().size());

        if ( mHromSimulation == true){
            // In case using the full modelpart, but only a set of selected elemets
            help_el_begin = mSelectedElements.begin();
            help_nelements = static_cast<int>(mSelectedElements.size());

            // Only selected conditions are considered for the calculation on an H-ROM simualtion.
            help_cond_begin = mSelectedConditions.begin();
            help_nconditions = static_cast<int>(mSelectedConditions.size());
        }

        // Getting the array of elements
        const auto nelements = help_nelements;
        const auto el_begin = help_el_begin;

        // Getting the array of the conditions
        const auto cond_begin = help_cond_begin;
        const auto nconditions = help_nconditions;


        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        const auto assembling_timer = BuiltinTimer();


        #pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, EquationId, el_begin, cond_begin)
        {
            Matrix PhiElemental;
            Matrix tempA = ZeroMatrix(mRomDofs,mRomDofs);
            Vector tempb = ZeroVector(mRomDofs);
            Matrix aux;

            #pragma omp for nowait
            for (int k = 0; k < static_cast<int>(nelements); k++)
            {
                auto it_el = el_begin + k;
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if ((it_el)->IsDefined(ACTIVE))
                    element_is_active = (it_el)->Is(ACTIVE);

                if (element_is_active){
                    //calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it_el, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    Element::DofsVectorType dofs;
                    it_el->GetDofList(dofs, CurrentProcessInfo);
                    const auto &geom = it_el->GetGeometry();
                    if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                        PhiElemental.resize(dofs.size(), mRomDofs,false);
                    if(aux.size1() != dofs.size() || aux.size2() != mRomDofs)
                        aux.resize(dofs.size(), mRomDofs,false);
                    GetPhiElemental(PhiElemental, dofs, geom);
                    noalias(aux) = prod(LHS_Contribution, PhiElemental);
                    double h_rom_weight = it_el->GetValue(HROM_WEIGHT);
                    noalias(tempA) += prod(trans(PhiElemental), aux) * h_rom_weight;
                    noalias(tempb) += prod(trans(PhiElemental), RHS_Contribution) * h_rom_weight;
                }
            }

            #pragma omp for nowait
            for (int k = 0; k < static_cast<int>(nconditions); k++){
                auto it = cond_begin + k;

                //detect if the element is active or not. If the user did not make any choice the condition
                //is active by default
                bool condition_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    condition_is_active = (it)->Is(ACTIVE);
                if (condition_is_active){
                    Condition::DofsVectorType dofs;
                    it->GetDofList(dofs, CurrentProcessInfo);
                    //calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    const auto &geom = it->GetGeometry();
                    if(PhiElemental.size1() != dofs.size() || PhiElemental.size2() != mRomDofs)
                        PhiElemental.resize(dofs.size(), mRomDofs,false);
                    if(aux.size1() != dofs.size() || aux.size2() != mRomDofs)
                        aux.resize(dofs.size(), mRomDofs,false);
                    GetPhiElemental(PhiElemental, dofs, geom);
                    noalias(aux) = prod(LHS_Contribution, PhiElemental);
                    double h_rom_weight = it->GetValue(HROM_WEIGHT);
                    noalias(tempA) += prod(trans(PhiElemental), aux) * h_rom_weight;
                    noalias(tempb) += prod(trans(PhiElemental), RHS_Contribution) * h_rom_weight;
                }
            }

            #pragma omp critical
            {
                noalias(Arom) +=tempA;
                noalias(brom) +=tempb;
            }

        }

        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Build time: " << assembling_timer.ElapsedSeconds() << std::endl;
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished parallel building" << std::endl;


        //solve for the rom unkowns dunk = Arom^-1 * brom
        Vector dxrom(xrom.size());
        const auto solving_timer = BuiltinTimer();
        MathUtils<double>::Solve(Arom, dxrom, brom);
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Solve reduced system time: " << solving_timer.ElapsedSeconds() << std::endl;

        // //update database
        // noalias(xrom) += dxrom;

        // project reduced solution back to full order model
        const auto backward_projection_timer = BuiltinTimer();
        ProjectToFineBasis(dxrom, rModelPart, Dx);
        KRATOS_INFO_IF("ROMBuilderAndSolver", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Project to fine basis time: " << backward_projection_timer.ElapsedSeconds() << std::endl;
    }

    void ResizeAndInitializeVectors(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType &pA,
        TSystemVectorPointerType &pDx,
        TSystemVectorPointerType &pb,
        ModelPart &rModelPart) override
    {
        KRATOS_TRY
        if (pA == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
            pA.swap(pNewA);
        }
        if (pDx == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0));
            pDx.swap(pNewDx);
        }
        if (pb == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0));
            pb.swap(pNewb);
        }

        TSystemVectorType &Dx = *pDx;
        TSystemVectorType &b = *pb;

        if (Dx.size() != BaseType::mEquationSystemSize)
            Dx.resize(BaseType::mEquationSystemSize, false);
        if (b.size() != BaseType::mEquationSystemSize)
            b.resize(BaseType::mEquationSystemSize, false);

        KRATOS_CATCH("")
    }

    /*@} */
    /**@name Operations */
    /*@{ */

    /*@} */
    /**@name Access */
    /*@{ */

    /*@} */
    /**@name Inquiry */
    /*@{ */

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "ROMBuilderAndSolver";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    /*@} */
    /**@name Friends */
    /*@{ */

    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */

    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /** Pointer to the Model.
     */
    typename TLinearSolver::Pointer mpLinearSystemSolver;

    //DofsArrayType mDofSet;
    std::vector<DofPointerType> mDofList;

    bool mReshapeMatrixFlag = false;

    /// flag taking care if the dof set was initialized ot not
    bool mDofSetIsInitialized = false;

    /// flag taking in account if it is needed or not to calculate the reactions
    bool mCalculateReactionsFlag = false;

    /// number of degrees of freedom of the problem to be solve
    unsigned int mEquationSystemSize;
    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    int mEchoLevel = 0;

    TSystemVectorPointerType mpReactionsVector;

    std::vector<std::string> mNodalVariablesNames;
    int mNodalDofs;
    unsigned int mRomDofs;
    std::unordered_map<Kratos::VariableData::KeyType,int> mMapPhi;
    ModelPart::ElementsContainerType mSelectedElements;
    ModelPart::ConditionsContainerType mSelectedConditions;
    bool mHromSimulation = false;
    bool mHromWeightsInitialized = false;

    /*@} */
    /**@name Protected Operations*/
    /*@{ */

    /*@} */
    /**@name Protected  Access */
    /*@{ */

    /*@} */
    /**@name Protected Inquiry */
    /*@{ */

    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */

    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */

    /*@} */
    /**@name Member Variables */
    /*@{ */

    /*@} */
    /**@name Private Operators*/
    /*@{ */

    /*@} */
    /**@name Private Operations*/
    /*@{ */

    /*@} */
    /**@name Private  Access */
    /*@{ */

    /*@} */
    /**@name Private Inquiry */
    /*@{ */

    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    /*@} */

}; /* Class ROMBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */

/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_ROM_BUILDER_AND_SOLVER  defined */
