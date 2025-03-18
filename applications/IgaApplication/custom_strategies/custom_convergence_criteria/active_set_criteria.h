// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ /
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Andrea Gorgi
//

#pragma once

// System includes

// External includes

// Project includes
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "iga_application_variables.h"


namespace Kratos
{
///@addtogroup ContactStructuralMechanicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class ActiveSetCriteria
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Custom convergence criteria for the active set of the contact
 * @author Andrea Gorgi
 */
template<class TSparseSpace, class TDenseSpace>
class ActiveSetCriteria
    : public  ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ActiveSetCriteria
    KRATOS_CLASS_POINTER_DEFINITION( ActiveSetCriteria );

    /// The base class definition
    using BaseType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;

    /// The definition of the current class
    using ClassType = ActiveSetCriteria<TSparseSpace, TDenseSpace>;

    /// The dofs array type
    using DofsArrayType = typename BaseType::DofsArrayType;

    /// The sparse matrix type
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    /// The dense vector type
    using TSystemVectorType = typename BaseType::TSystemVectorType;

    // /// The GidIO type
    // using GidIOBaseType = GidIO<>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    explicit ActiveSetCriteria(
        )
        : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param ThisParameters The configuration parameters
     */
    explicit ActiveSetCriteria(Kratos::Parameters ThisParameters)
        : BaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());

        mParameters = ThisParameters;
        // this->AssignSettings(ThisParameters);
    }

    ///Copy constructor
    ActiveSetCriteria( ActiveSetCriteria const& rOther )
      :BaseType(rOther)
    {
    }

    /// Destructor
    ~ActiveSetCriteria() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param ThisParameters The configuration parameters
     */
    typename BaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    /**
     * @brief Criterias that need to be called before getting the solution
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    bool PreCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        // KRATOS_INFO_IF("PRE-CRITERIA:: \n ACTIVE SET CRITERION: Convergence achieved", this->GetEchoLevel()>=0);
        PostCriteria(rModelPart, rDofSet, rA, rDx, rb);

        return true;

        // return(this->PostCriteria(rModelPart, rDofSet, rA, rDx, rb));
    //     KRATOS_ERROR_IF_NOT(mParameters.Has("contact_model_part_name"))
    //         << "Missing \"contact_model_part_name\" section" << std::endl;

    //     ModelPart& contact_model_part = mpModel->HasModelPart(mParameters["contact_model_part_name"].GetString())
    //                                 ? mpModel->GetModelPart(mParameters["contact_model_part_name"].GetString())
    //                                 : mpModel->CreateModelPart(mParameters["contact_model_part_name"].GetString());

    //     // contact_model_part.RemoveConditionsFromAllLevels();

    //     // PHYSICS HERE OR IN IGA_MODELER???????????
    //     // const std::string& rDataFileName = mParameters.Has("physics_file_name")
    //     //     ? mParameters["physics_file_name"].GetString()
    //     //     : "physics.iga.json";

    //     //-----------------------------------------------------------------------------------------------
    //     // Obtain SLAVE interface b_reps
    //     const std::string slave_model_part_name = mParameters["contact_parameters"]["slave_model_part"]["sub_model_part_name"].GetString();

    //     KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(slave_model_part_name)) << "ERROR: SLAVE MODEL PART " 
    //                                             << slave_model_part_name << "NOT CREATED" << std::endl; 

    //     ModelPart& slave_model_part = mpModel->GetModelPart(slave_model_part_name);
    //     GeometriesArrayType geometry_list_slave;
    //     GetCadGeometryList(geometry_list_slave, slave_model_part, mParameters["contact_parameters"]["slave_model_part"]);

    //     const IndexType slave_property_id = mParameters["contact_parameters"]["slave_model_part"]["property_id"].GetInt();
    //     Properties::Pointer p_prop_slave = slave_model_part.pGetProperties(slave_property_id);

    //     // Obtain MASTER interface b_reps
    //     const std::string master_model_part_name = mParameters["contact_parameters"]["master_model_part"]["sub_model_part_name"].GetString();

    //     KRATOS_ERROR_IF_NOT(mpModel->HasModelPart(master_model_part_name)) << "ERROR: MASTER MODEL PART " 
    //                                             << master_model_part_name << "NOT CREATED" << std::endl; 

    //     ModelPart& master_model_part = mpModel->GetModelPart(master_model_part_name);
    //     GeometriesArrayType geometry_list_master;
    //     GetCadGeometryList(geometry_list_master, master_model_part, mParameters["contact_parameters"]["master_model_part"]);

    //     const IndexType master_property_id = mParameters["contact_parameters"]["master_model_part"]["property_id"].GetInt();
    //     Properties::Pointer p_prop_master = master_model_part.pGetProperties(master_property_id);

    //     auto couplingGeometry = Kratos::make_shared<NurbsCouplingGeometry2D<PointType, PointerVector<NodeType>>>(geometry_list_slave, geometry_list_master);

    //     //---------------------------------------------------------------
    //     // ÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇÇ
    //     GeometriesArrayType geometries;
    //     SizeType shape_function_derivatives_order = 1;
    //     if (mParameters.Has("shape_function_derivatives_order")) {
    //         shape_function_derivatives_order = mParameters["shape_function_derivatives_order"].GetInt();
    //     }
    //     else {
    //         KRATOS_INFO_IF("CreateQuadraturePointGeometries", mEchoLevel > 4)
    //             << "shape_function_derivatives_order is not provided and thus being considered as 1. " << std::endl;
    //     }

    //     std::string quadrature_method = mParameters.Has("quadrature_method")
    //         ? mParameters["integration_rule"].GetString()
    //         : "GAUSS";
    //     IntegrationInfo integration_info = couplingGeometry->GetDefaultIntegrationInfo();
    //     for (IndexType i = 0; i < integration_info.LocalSpaceDimension(); ++i) {
    //         if (quadrature_method == "GAUSS") {
    //             integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GAUSS);
    //         }
    //         else if (quadrature_method == "EXTENDED_GAUSS") {
    //             integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::EXTENDED_GAUSS);
    //         }
    //         else if (quadrature_method == "GRID") {
    //             integration_info.SetQuadratureMethod(0, IntegrationInfo::QuadratureMethod::GRID);
    //         }
    //         else {
    //             KRATOS_INFO("CreateQuadraturePointGeometries") << "Quadrature method: " << quadrature_method
    //                 << " is not available. Available options are \"GAUSS\" and \"GRID\". Default quadrature method is being considered." << std::endl;
    //         }
    //     }

    //     couplingGeometry->CreateQuadraturePointGeometries(geometries, shape_function_derivatives_order, integration_info);

    //     SizeType id = 1;
    //     // if (contact_model_part.GetRootModelPart().Conditions().size() > 0)
    //     //     id = contact_model_part.GetRootModelPart().Conditions().back().Id() + 1;
    //     if (contact_model_part.GetRootModelPart().Conditions().size() > 0)
    //         id = contact_model_part.GetRootModelPart().Conditions().back().Id() + 1;
        
    //     KRATOS_ERROR_IF_NOT(mParameters.Has("name"))
    //         << "\"name\" need to be specified." << std::endl;
    //     std::string name = mParameters["name"].GetString();


    //     this->CreateConditions(
    //                     geometries.ptr_begin(), geometries.ptr_end(),
    //                     contact_model_part, name, id, p_prop_master, p_prop_slave);


        // return true;
    }

    /**
     * @brief Compute relative and absolute error.
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    bool PostCriteria(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {

        if (!rModelPart.HasSubModelPart("ContactInterface")) return true;
        ModelPart& r_contact_model_part = rModelPart.GetSubModelPart("ContactInterface");

        for (auto& contact_sub_model_part : r_contact_model_part.SubModelParts())
        {
            KRATOS_INFO_IF("::[ActiveSetCriteria]:: Starting Post Criteria of Contact sub model part ", this->GetEchoLevel() >= 1) 
                            << contact_sub_model_part.Name();
        
            if (contact_sub_model_part.NumberOfConditions() == 0) return true;

            int count_cond = 0;
            int n_CP = (contact_sub_model_part.Conditions().begin())->GetGeometry().GetGeometryPart(0).size();
            int p = (int) sqrt(n_CP);

            int n_GP_per_segment = 2*p+1;
            double toll_tangent_distance = 1e-1;
            double toll = 0;//1e-9;
            double toll_stress = 1e-4;
            int n_cond = contact_sub_model_part.Conditions().size(); 

            Vector length = ZeroVector(n_cond);
            Vector check_per_segment = ZeroVector(n_cond);

            Vector check_per_segment_stress = ZeroVector(n_cond);
            Vector check_per_segment_gap = ZeroVector(n_cond);
            int n_changes = 0;

            int n_active = 0;

            for (auto i_cond(contact_sub_model_part.Conditions().begin()); i_cond != contact_sub_model_part.Conditions().end(); ++i_cond)
            {
                
                // double normal_gap = i_cond->GetValue(NORMAL_GAP);
                Vector normal_stress_slave = i_cond->GetValue(STRESS_SLAVE);
                Vector normal_slave = i_cond->GetValue(NORMAL_SLAVE);
                Vector normal_stress_master = i_cond->GetValue(STRESS_MASTER);
                Vector normal_master = i_cond->GetValue(NORMAL_MASTER);

                // //FIXME:
                // normal_master = (normal_master - normal_slave)/2;

                Vector gap = i_cond->GetValue(GAP);
                double normal_gap_master = inner_prod(gap, normal_master);
                double normal_gap_slave = -inner_prod(gap, normal_slave);

                double check_value_gap = -(0.5*normal_gap_master + 0.5* normal_gap_slave);

                double tangent_gap_master = norm_2(gap - 0.5*normal_master*normal_gap_master + 0.5 * normal_slave*normal_gap_slave);

                double weight = i_cond->GetValue(INTEGRATION_WEIGHT);

                double young_modulus_master = i_cond->GetValue(YOUNG_MODULUS_MASTER); 
                double young_modulus_slave = i_cond->GetValue(YOUNG_MODULUS_SLAVE); 
                const double gamma = young_modulus_slave/(young_modulus_master + young_modulus_slave);

                double true_normal_stress_master = (normal_stress_master[0]* normal_master[0] + normal_stress_master[2]* normal_master[1])*normal_master[0] +
                                                (normal_stress_master[2]* normal_master[0] + normal_stress_master[1]* normal_master[1])*normal_master[1];

                double true_normal_stress_slave = (normal_stress_slave[0]* normal_slave[0] + normal_stress_slave[2]* normal_slave[1])*normal_slave[0] +
                                                (normal_stress_slave[2]* normal_slave[0] + normal_stress_slave[1]* normal_slave[1])*normal_slave[1];


                // int segment_index = (int) count_cond/n_GP_per_segment;
                // double check_value = -(true_normal_stress_master+young_modulus*normal_gap_master);

                
                double check_value_stress = -(gamma*true_normal_stress_master + (1-gamma) *true_normal_stress_slave)/std::min(young_modulus_master, young_modulus_slave);
                // // double check_value = -(yound_modulus*normal_gap);

                // length[segment_index] += weight;
                // check_per_segment[segment_index] += weight*check_value;

                // check_per_segment_stress[segment_index] += check_value_stress;
                // check_per_segment_gap[segment_index] += check_value_gap;

                // // FIXME:
                // if (i_cond->GetValue(SKIN_MASTER_COORDINATES)[0] < 0.01)
                // {
                //     if (i_cond->GetValue(ACTIVATION_LEVEL) == 0)
                //     {
                //         i_cond->SetValue(ACTIVATION_LEVEL, 1);
                //         n_changes++;
                //     }
                //     n_active ++;
                // }
                // else 

                if (check_value_stress< -toll_stress)
                {   
                    if (i_cond->GetValue(ACTIVATION_LEVEL) == 1)
                    {
                    i_cond->SetValue(ACTIVATION_LEVEL, 0);
                    n_changes++;
                    }

                } else if ((check_value_gap > toll) && std::abs(tangent_gap_master) < toll_tangent_distance) //|| check_value_gap + check_value_stress > toll
                {
                    if (i_cond->GetValue(ACTIVATION_LEVEL) == 0)
                    {
                        i_cond->SetValue(ACTIVATION_LEVEL, 1);
                        n_changes++;
                    }
                    n_active ++;
                }
                else if ( i_cond->GetValue(ACTIVATION_LEVEL) == 1)
                    n_active ++;

                count_cond++;

            }
            if (n_active == 0) 
                KRATOS_WATCH("[Warning]:: zero active contact conditions")
                                                                                
            // count_cond = 0;
            // for (auto i_cond(contact_sub_model_part->Conditions().begin()); i_cond != contact_sub_model_part->Conditions().end(); ++i_cond)
            // {
            //     int segment_index = (int) count_cond/n_GP_per_segment;

            //     // // OLD VERSION 
            //     // if (check_per_segment[segment_index]/length[segment_index] > toll)
            //     // {
            //     //     if (i_cond->GetValue(ACTIVATION_LEVEL) == 0)
            //     //     {
            //     //         i_cond->SetValue(ACTIVATION_LEVEL, 1);
            //     //         n_changes++;
            //     //     }
            //     // } else {
            //     //     if (i_cond->GetValue(ACTIVATION_LEVEL) == 1)
            //     //     {
            //     //         i_cond->SetValue(ACTIVATION_LEVEL, 0);
            //     //         n_changes++;
            //     //     }
            //     // }

            //     // NEW VERSION 
            //     if (i_cond->GetValue(ACTIVATION_LEVEL) == 1 &&
            //         check_per_segment_stress[segment_index]/length[segment_index] < -toll)
            //     {
            //         i_cond->SetValue(ACTIVATION_LEVEL, 0);
            //         n_changes++;
            //     }
            //     else if (i_cond->GetValue(ACTIVATION_LEVEL) == 0 &&
            //         check_per_segment_gap[segment_index]/length[segment_index] > toll)
            //     {
            //         i_cond->SetValue(ACTIVATION_LEVEL, 1);
            //         n_changes++;
            //     }

            //     // if (i_cond->GetValue(ACTIVATION_LEVEL) == 1) n_active++;


            //     count_cond++;

            // }
            double min_percentage_change = 0.05; //8/n_cond;
            double rel_change = (double) n_changes/n_active; //n_changes/n_cond;
            if (n_active == 0) rel_change = 0;
            if (rel_change <= min_percentage_change){
                KRATOS_INFO_IF("ACTIVE SET CRITERION: Convergence achieved", this->GetEchoLevel()>=0)
                << n_changes << " changes over " << n_active << " conditions active-> CHANGE: " << rel_change << std::endl;
                return true;
            } 
            else{
                KRATOS_INFO_IF("ACTIVE SET CRITERION: Convergence NOT achieved. -> ", this->GetEchoLevel()>=0) 
                << n_changes << " changes over " << n_active << " conditions active-> CHANGE: " << rel_change*100 << "%" << std::endl;
                return false;
            } 
        }
    }

    /**
     * @brief This function initialize the convergence criteria
     * @param rModelPart The model part of interest
     */
    void Initialize(ModelPart& rModelPart) override
    {
        // // Calling base criteria
        // BaseType::Initialize(rModelPart);

        // // The current process info
        // ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        // r_process_info.SetValue(ACTIVE_SET_COMPUTED, false);
    }

    /**
     * @brief This function initializes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        // Update normal of the conditions
        // ModelPart& r_contact_model_part = rModelPart.GetSubModelPart("Contact");
        // NormalCalculationUtils().CalculateUnitNormals<ModelPart::ConditionsContainerType>(r_contact_model_part, true);
        // const bool frictional_problem = rModelPart.IsDefined(SLIP) ? rModelPart.Is(SLIP) : false;
        // if (frictional_problem) {
        //     const bool has_lm = rModelPart.HasNodalSolutionStepVariable(VECTOR_LAGRANGE_MULTIPLIER);
        //     if (has_lm && mOptions.IsNot(ActiveSetCriteria::PURE_SLIP)) {
        //         MortarUtilities::ComputeNodesTangentModelPart(r_contact_model_part);
        //     } else {
        //         MortarUtilities::ComputeNodesTangentModelPart(r_contact_model_part, &WEIGHTED_SLIP, 1.0, true);
        //     }
        // }

        // IO for debugging
        // if (true || mOptions.Is(ActiveSetCriteria::IO_DEBUG)) {
        //     mpIO->CloseResultFile();
        //     std::ostringstream new_name ;
        //     new_name << "POST_LINEAR_ITER_STEP=""POST_LINEAR_ITER_STEP=" << rModelPart.GetProcessInfo()[STEP];
        //     mpIO->ChangeOutputName(new_name.str());
        // }
    }

    /**
     * @brief This function finalizes the solution step
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual)
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        // // IO for debugging
        // if (mOptions.Is(ActiveSetCriteria::IO_DEBUG)) {
        //     mpIO->FinalizeResults();
        // }
    }

    /**
     * @brief This function finalizes the non-linear iteration
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param rA System matrix (unused)
     * @param rDx Vector of results (variations on nodal variables)
     * @param rb RHS vector (residual + reactions)
     */
    void FinalizeNonLinearIteration(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        const TSystemMatrixType& rA,
        const TSystemVectorType& rDx,
        const TSystemVectorType& rb
        ) override
    {
        // Calling base criteria
        BaseType::FinalizeNonLinearIteration(rModelPart, rDofSet, rA, rDx, rb);

        // // The current process info
        // ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        // r_process_info.SetValue(ACTIVE_SET_COMPUTED, false);
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "echo_level": 4,
            "contact_model_part_name": "IgaModelPart.ContactInterface",
            //"quadrature_method": "GAUSS",
            "shape_function_derivatives_order": 3,
            "name": "SupportContact2DCondition",
            "contact_parameters" : {
                "slave_model_part" : {
                "sub_model_part_name": "IgaModelPart.Body2",
                "property_id" : 3,
                "brep_ids": [5]
                },
                "master_model_part" : {
                "sub_model_part_name": "IgaModelPart.Body1",
                "property_id" : 2,
                "brep_ids": [3]
                },
                "is_SBM" : false
            }
        })" );

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "base_mortar_criteria";
    }

    ///@}
    ///@name Acces
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ActiveSetCriteria";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    Flags mOptions; /// Local flags

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        // BaseType::AssignSettings(ThisParameters);

        // // Set local flags
        // mOptions.Set(ActiveSetCriteria::COMPUTE_DYNAMIC_FACTOR, ThisParameters["compute_dynamic_factor"].GetBool());
        // mOptions.Set(ActiveSetCriteria::IO_DEBUG, ThisParameters["gidio_debug"].GetBool());
        // mOptions.Set(ActiveSetCriteria::PURE_SLIP, ThisParameters["pure_slip"].GetBool());

        // if (mOptions.Is(ActiveSetCriteria::IO_DEBUG)) {
        //     mpIO = Kratos::make_shared<GidIOBaseType>("POST_LINEAR_ITER", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
        // }
    }

    /**
     * @brief This method resets the weighted gap in the nodes of the problem
     * @param rModelPart Reference to the ModelPart containing the contact problem.
     */
    virtual void ResetWeightedGap(ModelPart& rModelPart)
    {
        // auto& r_nodes_array = rModelPart.GetSubModelPart("Contact").Nodes();
        // VariableUtils().SetVariable(WEIGHTED_GAP, 0.0, r_nodes_array);
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    Parameters mParameters;

    // GidIOBaseType::Pointer mpIO; /// The pointer to the debugging IO

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief It computes the mean of the normal in the condition in all the nodes
     * @param rModelPart The model part to compute
     */
    inline void ComputeNodesMeanNormalModelPartWithPairedNormal(ModelPart& rModelPart)
    {
        // // Compute normal and tangent
        // ModelPart& r_contact_model_part = rModelPart.GetSubModelPart("Contact");
        // NormalCalculationUtils().CalculateUnitNormals<ModelPart::ConditionsContainerType>(r_contact_model_part, true);

        // // Iterate over the computing conditions
        // ModelPart& r_computing_contact_model_part = rModelPart.GetSubModelPart("ComputingContact");
        // auto& r_conditions_array = r_computing_contact_model_part.Conditions();
        // block_for_each(r_conditions_array, [&](Condition& rCond) {
        //     // Aux coordinates
        //     Point::CoordinatesArrayType aux_coords;

        //     // We update the paired normal
        //     GeometryType& r_parent_geometry = rCond.GetGeometry().GetGeometryPart(0);
        //     aux_coords = r_parent_geometry.PointLocalCoordinates(aux_coords, r_parent_geometry.Center());
        //     rCond.SetValue(NORMAL, r_parent_geometry.UnitNormal(aux_coords));
        // });
    }

    ///@}
}; // Class ActiveSetCriteria


}  // namespace Kratos
