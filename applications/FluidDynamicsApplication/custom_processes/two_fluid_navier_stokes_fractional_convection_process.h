//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Uxue Chasco
//                   
//                   
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "fluid_dynamics_application_variables.h"
#include "custom_elements/two_fluid_navier_stokes_fractional_convection.h"
#include "geometries/geometry_data.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "processes/find_nodal_h_process.h"
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
    /**Takes a model part full of SIMPLICIAL ELEMENTS (triangles and tetras) and convects the fractional velocity variable 
     * on the top of it as the first part of the fractinal splitting of two fluid navier stokes approach. 
     */
    template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
    class TwoFluidNavierStokesFractionalConvectionProcess
        : public Process
    {
    private:
        ///@name Type Definitions

        ///@}

    public:
        ///@name Type Definitions
        ///@{

        typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> SolvingStrategyType;
        typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderAndSolverPointerType;


        ///@}
        ///@name Pointer Definitions
        ///@{

        /// Pointer definition of TwoFluidNavierStokesFractionalConvectionProcess
        KRATOS_CLASS_POINTER_DEFINITION(TwoFluidNavierStokesFractionalConvectionProcess);
        ///@}
        ///@name Life Cycle
        ///@{

        /*
         * @brief Construct a new Navier Stokes Fractional Convection Process object
         * NS Fractional convection proces model constructor
         * @param rModel Model container
         * @param pLinearSolver Linear solver to be used in the level set convection problem
         * @param ThisParameters Json settings encapsulating the process configuration (see also GetDefaultParameters)
         */

        TwoFluidNavierStokesFractionalConvectionProcess(
            Model &rModel,
            typename TLinearSolver::Pointer pLinearSolver,
            Parameters ThisParameters)
            :mrModel(rModel),
            mrBaseModelPart(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()))
            {

            KRATOS_TRY
            auto p_builder_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>>(pLinearSolver);
            InitializeConvectionStrategy(p_builder_solver);
            // Validate the common settings as well as the element formulation specific ones
            // Checks and assign all the required member variables
            CheckAndAssignSettings(ThisParameters);
            ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
            ThisParameters["element_settings"].ValidateAndAssignDefaults(GetNSFractionalConvectionElementDefaultParameters(ThisParameters["element_type"].GetString()));
            
            std::string element_type = ThisParameters["element_type"].GetString();

            KRATOS_CATCH("")
        }
        
        /// Copy constructor.
        TwoFluidNavierStokesFractionalConvectionProcess(TwoFluidNavierStokesFractionalConvectionProcess const &rOther) = delete;

        /// Destructor.
        ~TwoFluidNavierStokesFractionalConvectionProcess() override
        {
            mrModel.DeleteModelPart(mAuxModelPartName);
        }

        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        /**
         * @brief Perform Navier Stokes Fractional Vectorial Convection.
         * This solver provides a stabilized convection problem for the fractional velocity unknown.
         * This process is always called within the framework of the two-fluid Navier-Stokes fractional solver, where the Navier-Stokes equations for two fluids are solved using a splitting approach.
         * The first part of the splitting involves solving the convection of the fractional velocity 
         */


        void Execute() override
        {
            KRATOS_TRY;

            // Fill the auxiliary convection model part if not done yet
            if (mFractionalVelocityPartIsInitialized == false)
            {
                ReGenerateNSFractionalConvectionModelPart(mrBaseModelPart);
            }

            // If required, calculate nodal element size
            // Note that this is done once assuming no mesh deformation
            if (mElementTauNodal || mCalculateNodalH)
            {
                ComputeNodalH();
                mCalculateNodalH = false;
            }

            // We set these values at every time step as other processes/solvers also use them
            // Note that this function sets the element related data (e.g. stabilization parameters)
            auto fill_process_info_function = GetFillProcessInfoFormulationDataFunction();
            fill_process_info_function(mrBaseModelPart);

            // Set convection problem data
            auto &r_conv_process_info = mpFractionalVelocityModelPart->GetProcessInfo();
            const double previous_delta_time = r_conv_process_info.GetValue(DELTA_TIME);
            const double dt = previous_delta_time;

            // If the nodal stabilization tau is to be used, it is also computed in here
            IndexPartition<int>(mpFractionalVelocityModelPart->NumberOfNodes()).for_each([&](int i_node){
                const auto it_node = mpFractionalVelocityModelPart->NodesBegin() + i_node;
                mVelocity[i_node] = it_node->FastGetSolutionStepValue(FRACTIONAL_VELOCITY);

                if (mElementTauNodal) {
                    double velocity_norm = norm_2(mVelocity[i_node]);
                    const double nodal_h = it_node-> GetValue(NODAL_H);
                    const double dynamic_tau = r_conv_process_info.GetValue(DYNAMIC_TAU);
                    const double tau = 1.0 / (dynamic_tau / dt + velocity_norm / std::pow(nodal_h,2));

                    it_node->GetValue(TAU) = tau;
                }
            });

            mpSolvingStrategy->InitializeSolutionStep();
            mpSolvingStrategy->Predict();
            mpSolvingStrategy->SolveSolutionStep(); // forward convection to reach phi_n+1
            mpSolvingStrategy->FinalizeSolutionStep();

            KRATOS_CATCH("")
        }

        void Clear() override
        {
            mpFractionalVelocityModelPart->Nodes().clear();
            mpFractionalVelocityModelPart->Conditions().clear();
            mpFractionalVelocityModelPart->Elements().clear();
            mFractionalVelocityPartIsInitialized = false;

            mpSolvingStrategy->Clear();

            mVelocity.clear();
        }

        const Parameters GetDefaultParameters() const override
        {
            Parameters default_parameters = Parameters(R"({
            "model_part_name"     :" ",
            "echo_level" : 0,
            "element_type" : "ns_fractional_velocity_convection",
            "element_settings" : {}
        })");

            return default_parameters;
        }

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
            return "TwoFluidNavierStokesFractionalConvectionProcess";
        }

        /// Print information about this object.
        void PrintInfo(std::ostream &rOStream) const override
        {
            rOStream << "TwoFluidNavierStokesFractionalConvectionProcess";
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

        Model &mrModel;
        ModelPart &mrBaseModelPart;

        ModelPart *mpFractionalVelocityModelPart = nullptr;

        bool mElementTauNodal;

        bool mCalculateNodalH = true;

        std::vector<array_1d<double, 3>> mVelocity;

        bool mFractionalVelocityPartIsInitialized = false;

        typename SolvingStrategyType::UniquePointer mpSolvingStrategy;

        std::string mAuxModelPartName;

        std::string mNSFractionalConvectionElementType;

        const Element *mpNSFractionalConvectionFactoryElement = nullptr;

        Parameters mNSFractionalConvectionSettings;

        ///@}
        ///@name Protected Operators
        ///@{

        ///@}
        ///@name Protected Operations
        ///@{


        virtual void ReGenerateNSFractionalConvectionModelPart(ModelPart &rBaseModelPart)
        {

            KRATOS_TRY

            KRATOS_ERROR_IF(mrModel.HasModelPart(mAuxModelPartName)) << "A process or operation using an auxiliar model_part with the name '" << mAuxModelPartName << "' already exists. Please choose another." << std::endl;

            mpFractionalVelocityModelPart = &(mrModel.CreateModelPart(mAuxModelPartName));

            // Check buffer size
            const auto base_buffer_size = rBaseModelPart.GetBufferSize();
            KRATOS_ERROR_IF(base_buffer_size < 2) << "Base model part buffer size is " << base_buffer_size << ". Set it to a minimum value of 2." << std::endl;

            // Generate
            mpFractionalVelocityModelPart->Nodes().clear();
            mpFractionalVelocityModelPart->Conditions().clear();
            mpFractionalVelocityModelPart->Elements().clear();

            mpFractionalVelocityModelPart->SetProcessInfo(rBaseModelPart.pGetProcessInfo());
            mpFractionalVelocityModelPart->SetBufferSize(base_buffer_size);
            for (auto it_properties = rBaseModelPart.PropertiesBegin(); it_properties != rBaseModelPart.PropertiesEnd(); ++it_properties)
            {
                mpFractionalVelocityModelPart->AddProperties(*(it_properties).base());
            }
            mpFractionalVelocityModelPart->Tables() = rBaseModelPart.Tables();

            // Assigning the nodes to the new model part
            mpFractionalVelocityModelPart->Nodes() = rBaseModelPart.Nodes();

            // Generating the elements
            mpFractionalVelocityModelPart->Elements().reserve(rBaseModelPart.NumberOfElements());
            KRATOS_ERROR_IF(mpNSFractionalConvectionFactoryElement == nullptr) << "Convection factory element has not been set yet." << std::endl;
            for (auto it_elem = rBaseModelPart.ElementsBegin(); it_elem != rBaseModelPart.ElementsEnd(); ++it_elem)
            {
                // Create the new element from the factory registered one
                auto p_element = mpNSFractionalConvectionFactoryElement->Create(
                    it_elem->Id(),
                    it_elem->pGetGeometry(),
                    it_elem->pGetProperties());

                mpFractionalVelocityModelPart->Elements().push_back(p_element);
            }

            // Initialize the nodal and elemental databases
            InitializeFractionalVelocityModelPartDatabases();

            // Resize the arrays
            const auto n_nodes = mpFractionalVelocityModelPart->NumberOfNodes();
            mVelocity.resize(n_nodes);
        
            mFractionalVelocityPartIsInitialized = true;

            KRATOS_CATCH("")
        }

        /**
         * @brief Initializes the databases values
         * This function initializes is intended to collect all the database initializations
         */
        void InitializeFractionalVelocityModelPartDatabases()
        {
            const array_1d<double, 3> aux_zero_vector = ZeroVector(3);
            if (mElementTauNodal)
            {
                block_for_each(mpFractionalVelocityModelPart->Nodes(), [&](Node &rNode)
                               { rNode.SetValue(TAU, 0.0); });
            }
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
            mNSFractionalConvectionSettings = ThisParameters;

            std::string element_register_name = GetNSFractionalConvectionElementRegisteredName(ThisParameters);


            mpNSFractionalConvectionFactoryElement = &KratosComponents<Element>::Get(element_register_name);

            mElementTauNodal = ThisParameters["element_settings"].Has("tau_nodal") ? ThisParameters["element_settings"]["tau_nodal"].GetBool() : false;
            KRATOS_WATCH(ThisParameters["model_part_name"].GetString())
            // Convection related settings
            if (ThisParameters["model_part_name"].GetString() == "")
            {
                
                mAuxModelPartName = mrBaseModelPart.Name() + "_FractionalVelocityConvectionPart";
            }
            else
            {
                mAuxModelPartName = ThisParameters["model_part_name"].GetString();
            }

            // Limiter related settings
            KRATOS_WATCH("")
        }

        std::string GetNSFractionalConvectionElementRegisteredName(Parameters ThisParameters)
        {
            // Convection element formulation settings
            std::string element_type = ThisParameters["element_type"].GetString();
            const auto element_list = GetNSFractionalConvectionElementsList();
            if (std::find(element_list.begin(), element_list.end(), element_type) == element_list.end())
            {
                KRATOS_INFO("") << "Specified \'" << element_type << "\' is not in the available elements list. " << "Attempting to use custom specified element." << std::endl;
                mNSFractionalConvectionElementType = element_type;
            }
            else
            {
                mNSFractionalConvectionElementType = GetNSFractionalConvectionElementName(element_type);
            }
            std::string element_register_name = mNSFractionalConvectionElementType + std::to_string(TDim) + "D" + std::to_string(TDim + 1) + "N";
            if (!KratosComponents<Element>::Has(element_register_name))
            {
                KRATOS_ERROR << "Specified \'" << element_type << "\' is not in the available elements list: " << element_list
                             << " and it is nor registered as a kratos element either. Please check your settings\n";
            }
            return element_register_name;
        }

        /**
         * @brief Get the Fractional Velocity Convection Elements List object
         * This method returns a list with the available formulations for the level set convection
         * @return const std::vector<std::string> List containing the available formulations
         */
        const virtual inline std::vector<std::string> GetNSFractionalConvectionElementsList()
        {
            std::vector<std::string> elements_list = {
                "ns_fractional_velocity_convection"};
            return elements_list;
        }

        /**
         * @brief Get the  Fractional Velocity Convection Element Name object
         * This method maps the user-defined element name to the Kratos class name
         * @param InputName User-defined element name
         * @return const std::string Kratos convection element class name
         */
        const virtual std::string GetNSFractionalConvectionElementName(std::string InputName)
        {
            const std::map<std::string, std::string> elements_name_map{

                {"ns_fractional_velocity_convection", "TwoFluidNavierStokesFractionalConvection"}

            };
            return elements_name_map.at(InputName);
        }

        /**
         * @brief Get the  Fractional Velocity Convection Element Default Parameters object
         * For each of the available formulations, this method returns the corresponding settings
         * @param ElementType User-defined element type
         * @return const Parameters Json string encapsulating the input element settings
         */
        const virtual Parameters GetNSFractionalConvectionElementDefaultParameters(const std::string ElementType)
        {
            Parameters default_parameters;

            if (ElementType == "ns_fractional_velocity_convection"){

                default_parameters = Parameters(R"({
                "dynamic_tau" : 0.1,
                "tau_nodal": false
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

            if (mNSFractionalConvectionElementType == "TwoFluidNavierStokesFractionalConvection")
            {
                fill_process_info_function = [this](ModelPart &rModelPart)
                {
                    auto &r_process_info = rModelPart.GetProcessInfo();
                    r_process_info.SetValue(DYNAMIC_TAU, mNSFractionalConvectionSettings["element_settings"]["dynamic_tau"].GetDouble());
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

            // Check the base model part element family (only simplex elements are supported)
            if constexpr (TDim == 2)
            {
                KRATOS_ERROR_IF(mrBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Triangle) << "In 2D the element type is expected to be a triangle" << std::endl;
            }
            else if constexpr (TDim == 3)
            {
                KRATOS_ERROR_IF(mrBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Tetrahedra) << "In 3D the element type is expected to be a tetrahedra" << std::endl;
            }

            // Generate an auxilary model part and populate it by elements of type TwoFluidNavierStokesFractionalConvection
            ReGenerateNSFractionalConvectionModelPart(mrBaseModelPart);

            // Generate  strategy
            bool CalculateReactions = false;
            bool ReformDofAtEachIteration = false;
            bool CalculateNormDxFlag = false;
            auto p_conv_criteria = Kratos::make_shared<DisplacementCriteria<TSparseSpace, TDenseSpace>>(1e-4, 1e-3);
            auto p_scheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>>();
            const std::size_t max_it = 10;
            mpSolvingStrategy = Kratos::make_unique<ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>>(
                *mpFractionalVelocityModelPart,
                p_scheme,
                p_conv_criteria,
                pBuilderAndSolver,
                max_it,
                CalculateReactions,
                ReformDofAtEachIteration,
                CalculateNormDxFlag);
            mpSolvingStrategy->SetEchoLevel(0);
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
        TwoFluidNavierStokesFractionalConvectionProcess &operator=(TwoFluidNavierStokesFractionalConvectionProcess const &rOther);

        ///@}
    }; // Class TwoFluidNavierStokesFractionalConvectionProcess

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
        TwoFluidNavierStokesFractionalConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver> &rThis);

    /// Output stream function
    template <unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
    inline std::ostream &operator<<(
        std::ostream &rOStream,
        const TwoFluidNavierStokesFractionalConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver> &rThis)
    {

        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

} // namespace Kratos.

