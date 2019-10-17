//    _____  _____ _                  __          __                                                 _ _           _   _
//   / ____|/ ____| |                 \ \        / /                               /\               | (_)         | | (_)
//  | |    | (___ | |__   __ _ _ __ _ _\ \  /\  / / __ __ _ _ __  _ __   ___ _ __ /  \   _ __  _ __ | |_  ___ __ _| |_ _  ___  _ __
//  | |     \___ \| '_ \ / _` | '__| '_ \ \/  \/ / '__/ _` | '_ \| '_ \ / _ \ '__/ /\ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//  | |____ ____) | | | | (_| | |  | |_) \  /\  /| | | (_| | |_) | |_) |  __/ | / ____ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//   \_____|_____/|_| |_|\__,_|_|  | .__/ \/  \/ |_|  \__,_| .__/| .__/ \___|_|/_/    \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                 | |                     | |   | |                    | |   | |
//                                 |_|                     |_|   |_|                    |_|   |_|
//
//
//  License: BSD License
//   license: CSharpWrapperApplication/license.txt
//
//  Main authors:    Hubert Balcerzak
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "custom_includes/kratos_internals.h"
#include "utilities/variable_utils.h"
#include "utilities/read_materials_utility.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"

//#define SKIN_SUBMODEL_PART_NAME "skin_model_part"

namespace CSharpKratosWrapper
{

void KratosInternals::initInternals() {
    mApplication.Register();
    mKernel.Initialize();
}

void KratosInternals::initModelPart() {
    mModel.Reset();

    // Accessing parameters
    mModelpartName = mSettingsParameters["solver_settings"]["model_part_name"].GetString();
    const std::size_t buffer_size = mSettingsParameters["solver_settings"]["buffer_size"].GetInt();
    const std::size_t domain_size = mSettingsParameters["solver_settings"]["domain_size"].GetInt();

    // Creating model part
    auto& r_model_part = mModel.CreateModelPart(mModelpartName, buffer_size);
    r_model_part.GetProcessInfo()[Kratos::DOMAIN_SIZE] = domain_size;

    // Basic variables
    r_model_part.AddNodalSolutionStepVariable(Kratos::DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(Kratos::REACTION);
    r_model_part.AddNodalSolutionStepVariable(Kratos::VOLUME_ACCELERATION);

    // Adding custom variables
    const std::size_t n_variables = mSettingsParameters["solver_settings"]["auxiliary_variables_list"].size();
    for (std::size_t i_var = 0; i_var < n_variables; ++i_var) {
        const std::string& r_variable_name = mSettingsParameters["solver_settings"]["auxiliary_variables_list"].GetArrayItem(i_var).GetString();
        if (Kratos::KratosComponents<Kratos::Variable<double>>::Has(r_variable_name)){
            const auto& r_variable = Kratos::KratosComponents<Kratos::Variable<double>>::Get(r_variable_name);
            r_model_part.AddNodalSolutionStepVariable(r_variable);
        } else if (Kratos::KratosComponents<Kratos::Variable<Kratos::array_1d<double, 3>>>::Has(r_variable_name)) {
            const auto& r_variable = Kratos::KratosComponents<Kratos::Variable<Kratos::array_1d<double, 3>>>::Get(r_variable_name);
            r_model_part.AddNodalSolutionStepVariable(r_variable);
        }
    }

    // NOTE: This is hardcoded
//    r_model_part.CreateSubModelPart(SKIN_SUBMODEL_PART_NAME);
}

void KratosInternals::loadMDPA(const std::string& rMDPAFilePath) {
    // Getting model part
    auto& r_model_part = GetMainModelPart();

    // Reading file
    Kratos::shared_ptr<std::fstream> pFile = Kratos::make_shared<std::fstream>();
    pFile->open(rMDPAFilePath, std::fstream::in);
    Kratos::ModelPartIO(pFile).ReadModelPart(r_model_part);
    pFile->close();
}

void KratosInternals::loadSettingsParameters(const std::string& rJSONFilePath) {
    if (rJSONFilePath != "") {
        std::ifstream infile(rJSONFilePath);
        if (!infile.good()) std::cout << "JSON file: " << rJSONFilePath  << " cannot be found" << std::endl;
        std::stringstream buffer;
        buffer << infile.rdbuf();
        mSettingsParameters = Kratos::Parameters(buffer.str());
    }

    Kratos::Parameters default_parameters = GetDefaultSettings();
    mSettingsParameters.RecursivelyAddMissingParameters(default_parameters);
}

void KratosInternals::initDofs() {
    // Getting model part
    auto& r_model_part = GetMainModelPart();

    // Basic variables
    Kratos::VariableUtils().AddDofWithReaction(Kratos::DISPLACEMENT_X, Kratos::REACTION_X, r_model_part);
    Kratos::VariableUtils().AddDofWithReaction(Kratos::DISPLACEMENT_Y, Kratos::REACTION_Y, r_model_part);
    Kratos::VariableUtils().AddDofWithReaction(Kratos::DISPLACEMENT_Z, Kratos::REACTION_Z, r_model_part);

    // Adding custom dofs
    typedef Kratos::VectorComponentAdaptor<Kratos::array_1d<double,3>> ComponentType;
    const std::size_t n_variables = mSettingsParameters["solver_settings"]["auxiliary_dofs_list"].size();
    for (std::size_t i_var = 0; i_var < n_variables; ++i_var) {
        const std::string& r_variable_name = mSettingsParameters["solver_settings"]["auxiliary_dofs_list"].GetArrayItem(i_var).GetString();
        const std::string& r_reaction_variable_name = mSettingsParameters["solver_settings"]["auxiliary_reaction_list"].GetArrayItem(i_var).GetString();
        if (Kratos::KratosComponents<Kratos::Variable<double>>::Has(r_variable_name)){
            const auto& r_variable = Kratos::KratosComponents<Kratos::Variable<double>>::Get(r_variable_name);
            const auto& r_reaction_variable = Kratos::KratosComponents<Kratos::Variable<double>>::Get(r_reaction_variable_name);
            Kratos::VariableUtils().AddDofWithReaction(r_variable, r_reaction_variable, r_model_part);
        } else if (Kratos::KratosComponents<Kratos::Variable<Kratos::array_1d<double, 3>>>::Has(r_variable_name)) {
            const auto& r_component_x = Kratos::KratosComponents<Kratos::VariableComponent<ComponentType>>::Get(r_variable_name + "_X");
            const auto& r_component_y = Kratos::KratosComponents<Kratos::VariableComponent<ComponentType>>::Get(r_variable_name + "_Y");
            const auto& r_component_z = Kratos::KratosComponents<Kratos::VariableComponent<ComponentType>>::Get(r_variable_name + "_Z");

            const auto& r_reaction_component_x = Kratos::KratosComponents<Kratos::VariableComponent<ComponentType>>::Get(r_reaction_variable_name + "_X");
            const auto& r_reaction_component_y = Kratos::KratosComponents<Kratos::VariableComponent<ComponentType>>::Get(r_reaction_variable_name + "_Y");
            const auto& r_reaction_component_z = Kratos::KratosComponents<Kratos::VariableComponent<ComponentType>>::Get(r_reaction_variable_name + "_Z");

            Kratos::VariableUtils().AddDofWithReaction(r_component_x, r_reaction_component_x, r_model_part);
            Kratos::VariableUtils().AddDofWithReaction(r_component_y, r_reaction_component_y, r_model_part);
            Kratos::VariableUtils().AddDofWithReaction(r_component_z, r_reaction_component_z, r_model_part);
        }
    }
}

void KratosInternals::initProperties() {
    // Getting model part
    auto& r_model_part = GetMainModelPart();

    const std::string& r_materials_filename = mSettingsParameters["solver_settings"]["material_import_settings"]["materials_filename"].GetString();

    // We will hardcode the CL
    if (r_materials_filename == "") {
        Kratos::ConstitutiveLaw::Pointer pCl = Kratos::make_shared<Kratos::HyperElasticIsotropicNeoHookean3D>();
        r_model_part.GetProperties(0).SetValue(Kratos::CONSTITUTIVE_LAW, pCl);
    } else { // We read the materials json
        Kratos::Parameters material_settings = Kratos::Parameters(R"({"Parameters": {"materials_filename": ""}})" );
        material_settings["Parameters"]["materials_filename"].SetString(r_materials_filename);
        Kratos::ReadMaterialsUtility(material_settings, mModel);
    }
}

void KratosInternals::initSolver() {
    LinearSolverType::Pointer pSolver = nullptr;
    // User specified a linear solver
    if (mSettingsParameters["solver_settings"]["linear_solver_settings"].Has("solver_type")) {
        pSolver = LinearSolverFactoryType().Create(mSettingsParameters["solver_settings"]["linear_solver_settings"]);
    } else { // Automatic
        std::array<std::string, 5> linear_solvers_by_speed = {"pardiso_lu", "sparse_lu", "pastix", "super_lu", "skyline_lu_factorization"};
        for (auto& r_name_solver :  linear_solvers_by_speed) {
            if (LinearSolverFactoryType().Has(r_name_solver)) {
                mSettingsParameters["solver_settings"]["linear_solver_settings"].AddEmptyValue("solver_type").SetString(r_name_solver);
                pSolver = LinearSolverFactoryType().Create(mSettingsParameters["solver_settings"]["linear_solver_settings"]);
                break;
            }
        }
    }
    ResidualBasedEliminationBuilderAndSolverType::Pointer pBuilderAndSolver = Kratos::make_shared<ResidualBasedEliminationBuilderAndSolverType>(pSolver);
    ResidualBasedIncrementalUpdateStaticSchemeType::Pointer pScheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticSchemeType>();
    const double residual_relative_tolerance = mSettingsParameters["solver_settings"]["residual_relative_tolerance"].GetDouble();
    const double residual_absolute_tolerance = mSettingsParameters["solver_settings"]["residual_absolute_tolerance"].GetDouble();
    ResidualCriteriaType::Pointer pConvergenceCriterion = Kratos::make_shared<ResidualCriteriaType>(residual_relative_tolerance, residual_absolute_tolerance);
    pConvergenceCriterion->SetEchoLevel(mSettingsParameters["solver_settings"]["echo_level"].GetInt());

    const int maxIters = mSettingsParameters["solver_settings"]["max_iteration"].GetInt();
    const bool computerReactions = mSettingsParameters["solver_settings"]["compute_reactions"].GetBool();
    const bool reformStepDofs = mSettingsParameters["solver_settings"]["reform_dofs_at_each_step"].GetBool();
    const bool moveMeshFlag = mSettingsParameters["solver_settings"]["move_mesh_flag"].GetBool();

    pmStrategy = Kratos::make_shared < ResidualBasedNewtonRaphsonStrategyType >(
        GetMainModelPart(),
        pScheme,
        pSolver,
        pConvergenceCriterion,
        pBuilderAndSolver,
        maxIters,
        computerReactions,
        reformStepDofs,
        moveMeshFlag);

    pmStrategy->SetEchoLevel(mSettingsParameters["problem_data"]["echo_level"].GetInt());
    pmStrategy->Check();
}

void KratosInternals::solve() {
    // TODO: Add everything related with the delta time, time, steps... (this is important in many solvers)
    pmStrategy->Solve();
}

Kratos::ModelPart& KratosInternals::GetMainModelPart() {
    return mModel.GetModelPart(mModelpartName);
}

Kratos::Parameters KratosInternals::GetSettings() {
    return mSettingsParameters;
}

Kratos::Parameters KratosInternals::GetDefaultSettings() {
    Kratos::Parameters default_parameters = Kratos::Parameters(R"(
    {
        "problem_data"    : {
            "problem_name"  : "Structure",
            "parallel_type" : "OpenMP",
            "start_time"    : 0.0,
            "end_time"      : 1.0,
            "echo_level"    : 0
        },
        "solver_settings" : {
            "model_part_name"                   : "Structure",
            "domain_size"                       : 3,
            "echo_level"                        : 0,
            "buffer_size"                       : 2,
            "analysis_type"                     : "non_linear",
            "model_import_settings"             : {
                "input_type"                        : "mdpa",
                "input_filename"                    : "unknown_name"
            },
            "computing_model_part_name"         : "computing_domain",
            "material_import_settings"          :{
                "materials_filename"                : ""
            },
            "time_stepping"                     : { },
            "rotation_dofs"                     : false,
            "reform_dofs_at_each_step"          : true,
            "line_search"                       : false,
            "compute_reactions"                 : true,
            "block_builder"                     : true,
            "clear_storage"                     : false,
            "move_mesh_flag"                    : true,
            "multi_point_constraints_used"      : true,
            "convergence_criterion"             : "residual_criterion",
            "displacement_relative_tolerance"   : 1.0e-4,
            "displacement_absolute_tolerance"   : 1.0e-9,
            "residual_relative_tolerance"       : 1.0e-4,
            "residual_absolute_tolerance"       : 1.0e-9,
            "max_iteration"                     : 10,
            "linear_solver_settings"            : { },
            "problem_domain_sub_model_part_list": [],
            "processes_sub_model_part_list"     : [],
            "auxiliary_variables_list"          : [],
            "auxiliary_dofs_list"               : [],
            "auxiliary_reaction_list"           : []
        },
        "processes"        : {},
        "output_processes" : {}
    })" );

    return default_parameters;
}
}
