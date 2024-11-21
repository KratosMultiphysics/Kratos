//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Ruben Zorrilla
//                   Mohammad Reza Hashemi
//

#if !defined(KRATOS_VECTORIAL_CONVECTION_PROCESS_INCLUDED)
#define KRATOS_VECTORIAL_CONVECTION_PROCESS_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "includes/convection_diffusion_settings.h"
#include "fluid_dynamics_application_variables.h"
#include "includes/define.h"
#include "includes/global_pointer_variables.h"
#include "includes/kratos_flags.h"
#include "custom_elements/vectorial_convection_fractional_element.h"
#include "geometries/geometry_data.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/pointer_communicator.h"
#include "utilities/pointer_map_communicator.h"
#include "processes/find_nodal_h_process.h"
#include "utilities/time_discretization.h"
#include "processes/calculate_nodal_area_process.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

namespace Kratos
{
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

    /// Short class definition.
    /**takes a model part full of SIMPLICIAL ELEMENTS (triangles and tetras) and convects a level set distance
     * on the top of it
     */
    template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
    class VectorialConvectionProcess
        : public Process
    {
    private:
        ///@name Type Definitions
        ///@{

        // class ProcessInfoDataContainer
        // {
        // public:
        //     ProcessInfoDataContainer(const ProcessInfo &rInputProcessInfo)
        //         : DeltaTime(rInputProcessInfo.GetValue(DELTA_TIME)), pUnknownVariable(FRACTIONAL_VELOCITY), pConvectionVariable(FRACTIONAL_VELOCITY)
        //     {
        //     }



        // private:
        //     const double DeltaTime;
        //     const Variable<array_1d<double, 3>> pUnknownVariable;
        //     const Variable<array_1d<double, 3>> pConvectionVariable;
        //     // const Variable<double> *pVolumeSourceVariable;
        // };

        ///@}

    public:
        ///@name Type Definitions
        ///@{

        typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> SolvingStrategyType;
        typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderAndSolverPointerType;


        ///@}
        ///@name Pointer Definitions
        ///@{

        /// Pointer definition of VectorialConvectionProcess
        KRATOS_CLASS_POINTER_DEFINITION(VectorialConvectionProcess);

        ///@}
        ///@name Life Cycle
        ///@{

        /**
         * @brief Construct a new Level Set Convection Process object
         * Level set convection proces model constructor
         * @param rModel Model container
         * @param pLinearSolver Linear solver to be used in the level set convection problem
         * @param ThisParameters Json settings encapsulating the process configuration (see also GetDefaultParameters)
         */
        VectorialConvectionProcess(
            Model &rModel,
            typename TLinearSolver::Pointer pLinearSolver,
            Parameters ThisParameters)
            : VectorialConvectionProcess(
                  rModel.GetModelPart(ThisParameters["model_part_name"].GetString()),
                  pLinearSolver,
                  ThisParameters)
        {
        }

        /**
         * @brief Construct a new Level Set Convection Process object
         * Level set convection proces model part constructor
         * @param rBaseModelPart Origin model part
         * @param pLinearSolver Linear solver to be used in the level set convection problem
         * @param ThisParameters Json settings encapsulating the process configuration (see also GetDefaultParameters)
         */
        VectorialConvectionProcess(
            ModelPart &rBaseModelPart,
            typename TLinearSolver::Pointer pLinearSolver,
            Parameters ThisParameters)
            : VectorialConvectionProcess(
                  rBaseModelPart,
                  ThisParameters)
        {
            KRATOS_TRY

            auto p_builder_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>>(pLinearSolver);
            InitializeConvectionStrategy(p_builder_solver);

            KRATOS_CATCH("")
        }

        /// Copy constructor.
        VectorialConvectionProcess(VectorialConvectionProcess const &rOther) = delete;

        /// Destructor.
        ~VectorialConvectionProcess() override
        {
            mrModel.DeleteModelPart(mAuxModelPartName);
        }

        ///@}
        ///@name Operators
        ///@{

        void operator()()
        {
            Execute();
        }

        ///@}
        ///@name Operations
        ///@{

        /**
         * @brief Perform the level-set convection
         * This solver provides a stabilized convection solver based on [Codina, R., 1993. Comput. Methods Appl. Mech. Engrg., 110(3-4), pp.325-342.]
         * It uses the sub-stepping approach to comply with the user defined maximum CFL number.
         * The error compensation is done according to the BFECC algorithm, which requires forward, backward, and the final forward solution steps (that triplicates the computational cost).
         * The error compensation severely disturbs the monotonicity of the results that is compensated for by implementing a limited BFECC algorithm.
         * The limiter relies on the nodal gradient of LevelSetVar (non-historical variable LevelSetGradientVar). For more info see [Kuzmin et al., Comput. Methods Appl. Mech. Engrg., 322 (2017) 23–41].
         */
        void ExecuteInitialize() override
        {

            ComputeNodalArea();
        }

        void Execute() override
        {
            KRATOS_TRY;

            // Fill the auxiliary convection model part if not done yet
            if (mDistancePartIsInitialized == false)
            {
                ReGenerateConvectionModelPart(mrBaseModelPart);
            }

            // If required, calculate nodal element size
            // KRATOS_WATCH("PRIMERO")
            // Note that this is done once assuming no mesh deformation
            if (mElementTauNodal || mCalculateNodalH)
            {
                ComputeNodalH();
                mCalculateNodalH = false;
            }

            // Evaluate steps needed to achieve target max_cfl
            // const auto n_substep = EvaluateNumberOfSubsteps();

            // Save the variables to be employed so that they can be restored after the solution
            // const auto process_info_data = ProcessInfoDataContainer(mpDistanceModelPart->GetProcessInfo());

            // We set these values at every time step as other processes/solvers also use them
            // Note that this function sets the element related data (e.g. stabilization parameters)
            auto fill_process_info_function = GetFillProcessInfoFormulationDataFunction();
            fill_process_info_function(mrBaseModelPart);

            // Set convection problem data
            auto &r_conv_process_info = mpDistanceModelPart->GetProcessInfo();
            const double previous_delta_time = r_conv_process_info.GetValue(DELTA_TIME);
            // const double dt =  previous_delta_time / static_cast<double>(n_substep);
            const double dt = previous_delta_time;



            // If the nodal stabilization tau is to be used, it is also computed in here
            IndexPartition<int>(mpDistanceModelPart->NumberOfNodes()).for_each([&](int i_node)
                                                                               {
            const auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
            mVelocity[i_node] = it_node->FastGetSolutionStepValue(FRACTIONAL_VELOCITY);
            // mVelocityOld[i_node] = it_node->FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 1);
            // mOldDistance[i_node] = it_node->FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 1);

            if (mElementTauNodal) {
                double velocity_norm = norm_2(mVelocity[i_node]);
                const double nodal_h = it_node-> GetValue(NODAL_H);
                const double dynamic_tau = r_conv_process_info.GetValue(DYNAMIC_TAU);
                const double tau = 1.0 / (dynamic_tau / dt + velocity_norm / std::pow(nodal_h,2));

                it_node->GetValue(TAU) = tau;
            } });

            ComputeNodalArea();
            mpSolvingStrategy->InitializeSolutionStep();
            // if (mIsBDFElement)
            // {
            //     NodalOSSProjection();
            // }
            mpSolvingStrategy->Predict();
            // block_for_each(mpDistanceModelPart->Nodes(), [&](Node &rNode){ rNode.FastGetSolutionStepValue(*mpVolumeSourceVar, 0.0); });
            mpSolvingStrategy->SolveSolutionStep(); // forward convection to reach phi_n+1
            mpSolvingStrategy->FinalizeSolutionStep();
            // NodalAccelerationProjection();

            // Reset the processinfo to the original settings
            // process_info_data.RestoreProcessInfoData(mpDistanceModelPart->GetProcessInfo());

            // Reset the velocities and levelset values to the one saved before the solution process
            // IndexPartition<int>(mpDistanceModelPart->NumberOfNodes()).for_each([&](int i_node)
            //                                                                    {
            //                                                                        auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
            //                                                                        it_node->FastGetSolutionStepValue(FRACTIONAL_VELOCITY) = mVelocity[i_node];
            //                                                                        it_node->FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 1) = mVelocityOld[i_node];
            //                                                                        it_node->FastGetSolutionStepValue(FRACTIONAL_VELOCITY, 1) = mOldDistance[i_node];
            //                                                                    });

            KRATOS_CATCH("")
        }

        void Clear() override
        {
            mpDistanceModelPart->Nodes().clear();
            mpDistanceModelPart->Conditions().clear();
            mpDistanceModelPart->Elements().clear();
            // mpDistanceModelPart->GetProcessInfo().clear();
            mDistancePartIsInitialized = false;

            mpSolvingStrategy->Clear();

            mVelocity.clear();
            // mVelocityOld.clear();
            // mOldDistance.clear();
            // mSigmaPlus.clear();
            // mSigmaMinus.clear();
            // mLimiter.clear();

            // mError.clear();
        }

        const Parameters GetDefaultParameters() const override
        {
            Parameters default_parameters = Parameters(R"({
            "model_part_name" : "",
            "echo_level" : 0,
            "convection_model_part_name" : "",
            "max_CFL" : 1.0,
            "element_type" : "levelset_convection_bdf",
            "element_settings" : {}
        })");

            return default_parameters;
        }

        // "levelset_volume_source_variable_name": "HEAT_FLUX",
        ///@}
        ///@name Access
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
            return "VectorialConvectionProcess";
        }

        /// Print information about this object.
        void PrintInfo(std::ostream &rOStream) const override
        {
            rOStream << "VectorialConvectionProcess";
        }

        /// Print object's data.
        void PrintData(std::ostream &rOStream) const override
        {
        }

        ///@}
        ///@name Friends
        ///@{

        ///@}
    protected:
        ///@name Protected static Member Variables
        ///@{

        ///@}
        ///@name Protected member Variables
        ///@{

        ModelPart &mrBaseModelPart;

        Model &mrModel;

        ModelPart *mpDistanceModelPart = nullptr;

        const Variable<array_1d<double, 3>> *mpLevelSetVar = nullptr;

        // const Variable<double>* mpVolumeSourceVar = nullptr;

        const Variable<array_1d<double, 3>> *mpConvectVar = nullptr;

        const Variable<array_1d<double, 3>> *mpLevelSetGradientVar = nullptr;

        double mMaxAllowedCFL = 1.0;

        unsigned int mMaxSubsteps = 0;

        bool mIsBfecc;

        bool mElementRequiresLimiter;

        bool mElementTauNodal;

        bool mCalculateNodalH = true;

        bool mElementRequiresLevelSetGradient;

        bool mEvaluateLimiter;

        double mPowerBfeccLimiter = 2.0;

        double mPowerElementalLimiter = 4.0;

        Vector mError;

        std::vector<array_1d<double, 3>> mOldDistance;

        Vector mSigmaPlus;

        Vector mSigmaMinus;

        Vector mLimiter;

        std::vector<array_1d<double, 3>> mVelocity;

        std::vector<array_1d<double, 3>> mVelocityOld;

        bool mDistancePartIsInitialized = false;

        typename SolvingStrategyType::UniquePointer mpSolvingStrategy;

        std::string mAuxModelPartName;

        std::string mConvectionElementType;

        const Element *mpConvectionFactoryElement = nullptr;

        Parameters mLevelSetConvectionSettings;


        ///@}
        ///@name Protected Operators
        ///@{

        ///@}
        ///@name Protected Operations
        ///@{

        VectorialConvectionProcess(
            ModelPart &rModelPart,
            Parameters ThisParameters)
            : mrBaseModelPart(rModelPart), mrModel(rModelPart.GetModel())
        {
            // Validate the common settings as well as the element formulation specific ones
            ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
            ThisParameters["element_settings"].ValidateAndAssignDefaults(GetConvectionElementDefaultParameters(ThisParameters["element_type"].GetString()));

            std::string element_type = ThisParameters["element_type"].GetString();
            KRATOS_WATCH(element_type)
            if (element_type == "levelset_convection_bdf")
            {
                mIsBDFElement = true;
            }
            // Checks and assign all the required member variables
            CheckAndAssignSettings(ThisParameters);

            // Sets the convection diffusion problem settings
            SetConvectionProblemSettings();
        }

        /**
         * @brief Set the level set convection formulation settings
         * This method sets the convection diffusion settings specifying the variable to be convect, its gradient, and the convection variable
         * Additionally, it also sets the required ProcessInfo variables
         */
        void SetConvectionProblemSettings()
        {
            // Get the base model part process info
            // Note that this will be shared with the auxiliary model part used in the convection resolution
            auto &r_process_info = mrBaseModelPart.GetProcessInfo();
            // Allocate if needed the variable CONVECTION_DIFFUSION_SETTINGS of the process info, and create it if it does not exist
            if (!r_process_info.Has(CONVECTION_DIFFUSION_SETTINGS))
            {
                auto p_conv_diff_settings = Kratos::make_shared<ConvectionDiffusionSettings>();
                r_process_info.SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_diff_settings);
                // p_conv_diff_settings->SetUnknownVariable(*mpLevelSetVar);
                // KRATOS_WATCH(*mpConvectVar)
                // p_conv_diff_settings->SetConvectionVariable(*mpConvectVar);
                // KRATOS_WATCH("SetConvectionProblemSettings")
                // // p_conv_diff_settings->SetVolumeSourceVariable(*mpVolumeSourceVar);
                // // KRATOS_WATCH(*mpVolumeSourceVar)
                // if (mpLevelSetGradientVar)
                // {
                //     p_conv_diff_settings->SetGradientVariable(*mpLevelSetGradientVar);
                // }
            }

            // This call returns a function pointer with the ProcessInfo filling directives
            // If the user-defined level set convection requires nothing to be set, the function does nothing
            auto fill_process_info_function = GetFillProcessInfoFormulationDataFunction();

            fill_process_info_function(mrBaseModelPart);
        }

        virtual void ReGenerateConvectionModelPart(ModelPart &rBaseModelPart)
        {

            KRATOS_TRY

            KRATOS_ERROR_IF(mrModel.HasModelPart(mAuxModelPartName)) << "A process or operation using an auxiliar model_part with the name '" << mAuxModelPartName << "' already exists. Please choose another." << std::endl;

            mpDistanceModelPart = &(mrModel.CreateModelPart(mAuxModelPartName));

            // Check buffer size
            const auto base_buffer_size = rBaseModelPart.GetBufferSize();
            KRATOS_ERROR_IF(base_buffer_size < 2) << "Base model part buffer size is " << base_buffer_size << ". Set it to a minimum value of 2." << std::endl;

            // Generate
            mpDistanceModelPart->Nodes().clear();
            mpDistanceModelPart->Conditions().clear();
            mpDistanceModelPart->Elements().clear();

            mpDistanceModelPart->SetProcessInfo(rBaseModelPart.pGetProcessInfo());
            mpDistanceModelPart->SetBufferSize(base_buffer_size);
            for (auto it_properties = rBaseModelPart.PropertiesBegin(); it_properties != rBaseModelPart.PropertiesEnd(); ++it_properties)
            {
                mpDistanceModelPart->AddProperties(*(it_properties).base());
            }
            mpDistanceModelPart->Tables() = rBaseModelPart.Tables();

            // Assigning the nodes to the new model part
            mpDistanceModelPart->Nodes() = rBaseModelPart.Nodes();

            // // Ensure that the nodes have distance as a DOF
            // VariableUtils().AddDof<Variable<array_1d<double,3>>>(FRACTIONAL_VELOCITY, rBaseModelPart);

            // Generating the elements
            mpDistanceModelPart->Elements().reserve(rBaseModelPart.NumberOfElements());
            KRATOS_ERROR_IF(mpConvectionFactoryElement == nullptr) << "Convection factory element has not been set yet." << std::endl;
            for (auto it_elem = rBaseModelPart.ElementsBegin(); it_elem != rBaseModelPart.ElementsEnd(); ++it_elem)
            {
                // Create the new element from the factory registered one
                auto p_element = mpConvectionFactoryElement->Create(
                    it_elem->Id(),
                    it_elem->pGetGeometry(),
                    it_elem->pGetProperties());

                mpDistanceModelPart->Elements().push_back(p_element);
            }

            // Initialize the nodal and elemental databases
            InitializeDistanceModelPartDatabases();

            // Resize the arrays
            const auto n_nodes = mpDistanceModelPart->NumberOfNodes();
            mVelocity.resize(n_nodes);
            mVelocityOld.resize(n_nodes);
            mOldDistance.resize(n_nodes);

            mDistancePartIsInitialized = true;

            KRATOS_CATCH("")
        }

        /**
         * @brief Initializes the databases values
         * This function initializes is intended to collect all the database initializations
         */
        void InitializeDistanceModelPartDatabases()
        {

            const array_1d<double, 3> aux_zero_vector = ZeroVector(3);
            if (mElementTauNodal)
            {
                block_for_each(mpDistanceModelPart->Nodes(), [&](Node &rNode)
                               { rNode.SetValue(TAU, 0.0); });
            }
        }


        void ComputeNodalArea()
        {
            // Calculate the NODAL_AREA
            CalculateNodalAreaProcess<false> nodal_area_process(mrBaseModelPart);
            nodal_area_process.Execute();
        }
        /**
         * @brief Nodal H calculation
         * This function calculates the nodal h  by executing a process where the nodal h calculaiton is implemented.
         */
        void ComputeNodalH()
        {
            auto nodal_h_process = FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable>(mrBaseModelPart);
            nodal_h_process.Execute();
        }

        ///@}
        ///@name Protected  Access
        ///@{

        ///@}
        ///@name Protected Inquiry
        ///@{

        ///@}
        ///@name Protected LifeCycle
        ///@{

        ///@}
    private:
        ///@name Static Member Variables
        ///@{

        ///@}
        ///@name Member Variables
        ///@{

        ///@}
        ///@name Private Operators
        ///@{
        bool mIsBDFElement = false;
        ///@}
        ///@name Private Operations
        ///@{

        /**
         * @brief Checks and assign the required member variables
         * This function checks the provided parameters, which need to have been already validated and sets the member variables
         * @param ThisParameters Json string containing the already validated process and formulation settings
         */
        void CheckAndAssignSettings(const Parameters ThisParameters)
        {
            mLevelSetConvectionSettings = ThisParameters;

            std::string element_register_name = GetConvectionElementRegisteredName(ThisParameters);


            mpConvectionFactoryElement = &KratosComponents<Element>::Get(element_register_name);

            mElementTauNodal = ThisParameters["element_settings"].Has("tau_nodal") ? ThisParameters["element_settings"]["tau_nodal"].GetBool() : false;

            // Convection related settings
            mMaxAllowedCFL = ThisParameters["max_CFL"].GetDouble();
            // mMaxSubsteps = ThisParameters["max_substeps"].GetInt();


            mMaxAllowedCFL = ThisParameters["max_CFL"].GetDouble();
            mpLevelSetVar = &KratosComponents<Variable<array_1d<double, 3>>>::Get("FRACTIONAL_VELOCITY");
            mpConvectVar = &KratosComponents<Variable<array_1d<double, 3>>>::Get("FRACTIONAL_VELOCITY");
            // mpVolumeSourceVar = &KratosComponents<Variable<double>>::Get(ThisParameters["levelset_volume_source_variable_name"].GetString());
            if (ThisParameters["convection_model_part_name"].GetString() == "")
            {
                mAuxModelPartName = mrBaseModelPart.Name() + "_DistanceConvectionPart";
            }
            else
            {
                mAuxModelPartName = ThisParameters["convection_model_part_name"].GetString();
            }

            // Limiter related settings
            KRATOS_WATCH("")
        }

        std::string GetConvectionElementRegisteredName(Parameters ThisParameters)
        {
            // Convection element formulation settings
            std::string element_type = ThisParameters["element_type"].GetString();
            const auto element_list = GetConvectionElementsList();
            if (std::find(element_list.begin(), element_list.end(), element_type) == element_list.end())
            {
                KRATOS_INFO("") << "Specified \'" << element_type << "\' is not in the available elements list. " << "Attempting to use custom specified element." << std::endl;
                mConvectionElementType = element_type;
            }
            else
            {
                mConvectionElementType = GetConvectionElementName(element_type);
            }
            std::string element_register_name = mConvectionElementType + std::to_string(TDim) + "D" + std::to_string(TDim + 1) + "N";
            if (!KratosComponents<Element>::Has(element_register_name))
            {
                KRATOS_ERROR << "Specified \'" << element_type << "\' is not in the available elements list: " << element_list
                             << " and it is nor registered as a kratos element either. Please check your settings\n";
            }
            return element_register_name;
        }

        /**
         * @brief Get the Convection Elements List object
         * This method returns a list with the available formulations for the level set convection
         * @return const std::vector<std::string> List containing the available formulations
         */
        const virtual inline std::vector<std::string> GetConvectionElementsList()
        {
            std::vector<std::string> elements_list = {
                "levelset_convection_bdf"};
            return elements_list;
        }

        /**
         * @brief Get the Convection Element Name object
         * This method maps the user-defined element name to the Kratos class name
         * @param InputName User-defined element name
         * @return const std::string Kratos convection element class name
         */
        const virtual std::string GetConvectionElementName(std::string InputName)
        {
            const std::map<std::string, std::string> elements_name_map{

                {"levelset_convection_bdf", "VectorialConvectionFractionalElement"}

            };
            return elements_name_map.at(InputName);
        }

        /**
         * @brief Get the Convection Element Default Parameters object
         * For each of the available formulations, this method returns the corresponding settings
         * @param ElementType User-defined element type
         * @return const Parameters Json string encapsulating the input element settings
         */
        const virtual Parameters GetConvectionElementDefaultParameters(const std::string ElementType)
        {
            Parameters default_parameters;

            if (ElementType == "levelset_convection_bdf"){

                default_parameters = Parameters(R"({
                "dynamic_tau" : 0.1,
                "tau_nodal": true
            })");
            }

            return default_parameters;
        }

        /**
         * @brief Get the Fill Process Info Function object
         * This method returns a lambda function with the required operations to be performed in the process info
         * It has to be particularised for all the formulations. If not particularised a do nothing instruction is returned
         * @return const std::function<void(ModelPart&)> A function pointer to be called when setting up the distance model part
         */
        const virtual std::function<void(ModelPart &)> GetFillProcessInfoFormulationDataFunction()
        {
            std::function<void(ModelPart &)> fill_process_info_function;

            if (mConvectionElementType == "VectorialConvectionFractionalElement")
            {
                fill_process_info_function = [this](ModelPart &rModelPart)
                {
                    auto &r_process_info = rModelPart.GetProcessInfo();
                    r_process_info.SetValue(DYNAMIC_TAU, mLevelSetConvectionSettings["element_settings"]["dynamic_tau"].GetDouble());
                };
            }
            else
            {
                fill_process_info_function = [](ModelPart &rModelPart) {};
            }

            return fill_process_info_function;
        }

        void InitializeConvectionStrategy(BuilderAndSolverPointerType pBuilderAndSolver)
        {
            // Check that there is at least one element and node in the model
            KRATOS_ERROR_IF(mrBaseModelPart.NumberOfNodes() == 0) << "The model has no nodes." << std::endl;
            KRATOS_ERROR_IF(mrBaseModelPart.NumberOfElements() == 0) << "The model has no elements." << std::endl;

            // Check that the level set and convection variables are in the nodal database
            // TODO: ESTE CHEKC TIRA ERROR DE
            // VariableUtils().CheckVariableExists<Variable<array_1d<double, 3>>>(FRACTIONAL_STEP, mrBaseModelPart.Nodes());
            // VariableUtils().CheckVariableExists<Variable<double>>(*mpVolumeSourceVar, mrBaseModelPart.Nodes());
            // Check the base model part element family (only simplex elements are supported)
            if constexpr (TDim == 2)
            {
                KRATOS_ERROR_IF(mrBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Triangle) << "In 2D the element type is expected to be a triangle" << std::endl;
            }
            else if constexpr (TDim == 3)
            {
                KRATOS_ERROR_IF(mrBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Tetrahedra) << "In 3D the element type is expected to be a tetrahedra" << std::endl;
            }

            // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
            ReGenerateConvectionModelPart(mrBaseModelPart);

            // Generate  strategy
            bool CalculateReactions = false;
            bool ReformDofAtEachIteration = false;
            bool CalculateNormDxFlag = false;
            auto p_conv_criteria = Kratos::make_shared<DisplacementCriteria<TSparseSpace, TDenseSpace>>(1e-10, 1e-9);
            auto p_scheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>>();
            const std::size_t max_it = 10;
            mpSolvingStrategy = Kratos::make_unique<ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>>(
                *mpDistanceModelPart,
                p_scheme,
                p_conv_criteria,
                pBuilderAndSolver,
                max_it,
                CalculateReactions,
                ReformDofAtEachIteration,
                CalculateNormDxFlag);
            mpSolvingStrategy->SetEchoLevel(1);
            mpSolvingStrategy->Check();
            mpSolvingStrategy->Solve();
        }

        ///@}
        ///@name Private  Access
        ///@{

        ///@}
        ///@name Private Inquiry
        ///@{

        ///@}
        ///@name Un accessible methods
        ///@{

        /// Assignment operator.
        VectorialConvectionProcess &operator=(VectorialConvectionProcess const &rOther);

        ///@}
    }; // Class VectorialConvectionProcess

    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Input stream function
    template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
    inline std::istream &operator>>(
        std::istream &rIStream,
        VectorialConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver> &rThis);

    /// Output stream function
    template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
    inline std::ostream &operator<<(
        std::ostream &rOStream,
        const VectorialConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver> &rThis)
    {

        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

} // namespace Kratos.

#endif // KRATOS_LEVELSET_CONVECTION_BDF_PROCESS_INCLUDED  defined
