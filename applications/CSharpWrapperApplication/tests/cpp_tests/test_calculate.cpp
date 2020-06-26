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
//  Main authors:  Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "custom_includes/kratos_wrapper.h"
#include "includes/kratos_filesystem.h"

namespace Kratos {
    namespace Testing {
        void CreateMDPAFile() {
            std::filebuf buffer;
            buffer.open(FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "file.mdpa"}), std::ios::out);
            std::ostream os(&buffer);
            os
                    << "Begin ModelPartData\nEnd ModelPartData\n\nBegin Properties  0\n    DENSITY 2700.000000\n    YOUNG_MODULUS 7000000.000000\n    POISSON_RATIO 0.300000\n    BODY_FORCE [3] (0.000000,0.000000,0.000000)\n    THICKNESS 1.000000\nEnd Properties\n\nBegin Nodes\n        1        0.0        0.0         0.0\n        2        0.0        0.0         1.0\n        3        1.0        0.0         0.0\n        4        1.0        1.0         0.0\nEnd Nodes\n\nBegin Elements SmallDisplacementElement3D4N\n    1 0 1 2 3 4\nEnd Elements\n\nBegin SubModelPart BasePart // Note that this would be a sub sub modelpart\n    Begin SubModelPartNodes\n        1\n        2\n    End SubModelPartNodes\n    Begin SubModelPart inner_part\n        Begin SubModelPartNodes\n            1\n        End SubModelPartNodes\n    End SubModelPart\nEnd SubModelPart";
            buffer.close();
        }

        void CreateJSONFile() {
            Parameters json_parameters = Parameters(R"(
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
            })");
            const std::string &r_json_text = json_parameters.PrettyPrintJsonString();
            std::filebuf buffer;
            buffer.open(FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "file.json"}), std::ios::out);
            std::ostream os(&buffer);
            os << r_json_text;
            buffer.close();
        }

        /**
        * Checks the correct work of update node position
        */
        KRATOS_TEST_CASE_IN_SUITE(CSharpWrapperUpdateNodePosition, KratosCSharpWrapperApplicationFastSuite) {
            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<Element>::Has("SmallDisplacementElement3D4N"))
                return void();

            // Import mdpa
            CreateMDPAFile();
            const std::string file_name = FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "file.mdpa"});
            CSharpKratosWrapper::KratosWrapper *wrapperInstance = new CSharpKratosWrapper::KratosWrapper();
            wrapperInstance->init(file_name.c_str());
            CSharpKratosWrapper::ModelPartWrapper *mainModelPart = wrapperInstance->getRootModelPartWrapper();
            // Get some API info
            mainModelPart->retrieveResults();
            float *x = mainModelPart->getXCoordinates();
            float *y = mainModelPart->getYCoordinates();
            float *z = mainModelPart->getZCoordinates();
//            int n = mainModelPart->getNodesCount();

//            const float float_epsilon = std::numeric_limits<float>::epsilon();

//             // Non-initialized (calculated)
//             for (int i = 0; i < n; i++) {
//                 KRATOS_CHECK_NEAR(x[i], 0.0, float_epsilon);
//                 KRATOS_CHECK_NEAR(y[i], 0.0, float_epsilon);
//                 KRATOS_CHECK_NEAR(z[i], 0.0, float_epsilon);
//             }

            // None fixed
            wrapperInstance->calculate();
            mainModelPart->retrieveResults();
            x = mainModelPart->getXCoordinates();
            y = mainModelPart->getYCoordinates();
            z = mainModelPart->getZCoordinates();

//             KRATOS_CHECK_NEAR(x[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[1], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[1], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[1], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[2], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[2], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[2], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[3], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[3], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[3], 0.0, float_epsilon);

            // All fixed
            mainModelPart->updateNodePos(0, x[0], y[0], z[0]);
            mainModelPart->updateNodePos(1, x[1], y[1], z[1]);
            mainModelPart->updateNodePos(2, x[2], y[2] + 1.0e-8, z[2]);
            mainModelPart->updateNodePos(3, x[3], y[3], z[3]);

            wrapperInstance->calculate();
            mainModelPart->retrieveResults();
            x = mainModelPart->getXCoordinates();
            y = mainModelPart->getYCoordinates();
            z = mainModelPart->getZCoordinates();

//             KRATOS_CHECK_NEAR(x[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[1], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[1], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[1], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[2], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[2], 1.0e-8, float_epsilon);
//             KRATOS_CHECK_NEAR(z[2], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[3], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[3], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[3], 0.0, float_epsilon);

            // Partially fixed
            mainModelPart->updateNodePos(0, x[0], y[0], z[0]);
            mainModelPart->updateNodePos(1, x[1], y[1], z[1]);
            mainModelPart->updateNodePos(2, x[2], y[2] + 1.0e-8, z[2]);

            wrapperInstance->calculate();
            mainModelPart->retrieveResults();
            x = mainModelPart->getXCoordinates();
            y = mainModelPart->getYCoordinates();
            z = mainModelPart->getZCoordinates();

//             KRATOS_CHECK_NEAR(x[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[1], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[1], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[1], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[2], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[2], 2.0e-8, float_epsilon);
//             KRATOS_CHECK_NEAR(z[2], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[3], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[3], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[3], 0.0, float_epsilon);

            Kratos::filesystem::remove((FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "file.mdpa"})).c_str());
        }

        /**
        * Checks the correct work of update node position with json input
        */
        KRATOS_TEST_CASE_IN_SUITE(CSharpWrapperUpdateNodePositionWithJSON, KratosCSharpWrapperApplicationFastSuite) {
            // In case the StructuralMechanicsApplciation is not compiled we skip the test
            if (!KratosComponents<Element>::Has("SmallDisplacementElement3D4N"))
                return void();

            // Import mdpa
            CreateMDPAFile();
            CreateJSONFile();
            const std::string file_name = FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "file.mdpa"});
            const std::string file_name_json = FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "file.json"});
            CSharpKratosWrapper::KratosWrapper *wrapperInstance = new CSharpKratosWrapper::KratosWrapper();

            wrapperInstance->init(file_name.c_str(), file_name_json.c_str());
            CSharpKratosWrapper::ModelPartWrapper *mainModelPart = wrapperInstance->getRootModelPartWrapper();

            // Get some API info
            mainModelPart->retrieveResults();
            float *x = mainModelPart->getXCoordinates();
            float *y = mainModelPart->getYCoordinates();
            float *z = mainModelPart->getZCoordinates();
//            int n = mainModelPart->getNodesCount();

//            const float float_epsilon = std::numeric_limits<float>::epsilon();

//             // Non-initialized (calculated)
//             for (int i = 0; i < n; i++) {
//                 KRATOS_CHECK_NEAR(x[i], 0.0, float_epsilon);
//                 KRATOS_CHECK_NEAR(y[i], 0.0, float_epsilon);
//                 KRATOS_CHECK_NEAR(z[i], 0.0, float_epsilon);
//             }

            // None fixed
            wrapperInstance->calculate();
            mainModelPart->retrieveResults();
            x = mainModelPart->getXCoordinates();
            y = mainModelPart->getYCoordinates();
            z = mainModelPart->getZCoordinates();

//             KRATOS_CHECK_NEAR(x[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[1], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[1], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[1], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[2], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[2], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[2], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[3], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[3], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[3], 0.0, float_epsilon);

            // All fixed
            mainModelPart->updateNodePos(0, x[0], y[0], z[0]);
            mainModelPart->updateNodePos(1, x[1], y[1], z[1]);
            mainModelPart->updateNodePos(2, x[2], y[2] + 1.0e-8, z[2]);
            mainModelPart->updateNodePos(3, x[3], y[3], z[3]);

            wrapperInstance->calculate();
            mainModelPart->retrieveResults();
            x = mainModelPart->getXCoordinates();
            y = mainModelPart->getYCoordinates();
            z = mainModelPart->getZCoordinates();

//             KRATOS_CHECK_NEAR(x[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[1], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[1], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[1], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[2], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[2], 1.0e-8, float_epsilon);
//             KRATOS_CHECK_NEAR(z[2], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[3], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[3], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[3], 0.0, float_epsilon);

            // Partially fixed
            mainModelPart->updateNodePos(0, x[0], y[0], z[0]);
            mainModelPart->updateNodePos(1, x[1], y[1], z[1]);
            mainModelPart->updateNodePos(2, x[2], y[2] + 1.0e-8, z[2]);

            wrapperInstance->calculate();
            mainModelPart->retrieveResults();
            x = mainModelPart->getXCoordinates();
            y = mainModelPart->getYCoordinates();
            z = mainModelPart->getZCoordinates();

//             KRATOS_CHECK_NEAR(x[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[0], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[1], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[1], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[1], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[2], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[2], 2.0e-8, float_epsilon);
//             KRATOS_CHECK_NEAR(z[2], 0.0, float_epsilon);
//             KRATOS_CHECK_NEAR(x[3], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(y[3], 1.0, float_epsilon);
//             KRATOS_CHECK_NEAR(z[3], 0.0, float_epsilon);

            Kratos::filesystem::remove((FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "file.mdpa"})).c_str());
            Kratos::filesystem::remove((FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "file.json"})).c_str());
        }

    } // namespace Testing
}  // namespace Kratos.
