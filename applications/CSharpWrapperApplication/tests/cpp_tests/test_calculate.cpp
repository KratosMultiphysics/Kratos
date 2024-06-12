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
//                 KRATOS_EXPECT_NEAR(x[i], 0.0, float_epsilon);
//                 KRATOS_EXPECT_NEAR(y[i], 0.0, float_epsilon);
//                 KRATOS_EXPECT_NEAR(z[i], 0.0, float_epsilon);
//             }

            // None fixed
            wrapperInstance->calculate();
            mainModelPart->retrieveResults();
            x = mainModelPart->getXCoordinates();
            y = mainModelPart->getYCoordinates();
            z = mainModelPart->getZCoordinates();

//             KRATOS_EXPECT_NEAR(x[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[1], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[1], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[1], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[2], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[2], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[2], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[3], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[3], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[3], 0.0, float_epsilon);

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

//             KRATOS_EXPECT_NEAR(x[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[1], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[1], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[1], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[2], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[2], 1.0e-8, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[2], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[3], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[3], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[3], 0.0, float_epsilon);

            // Partially fixed
            mainModelPart->updateNodePos(0, x[0], y[0], z[0]);
            mainModelPart->updateNodePos(1, x[1], y[1], z[1]);
            mainModelPart->updateNodePos(2, x[2], y[2] + 1.0e-8, z[2]);

            wrapperInstance->calculate();
            mainModelPart->retrieveResults();
            x = mainModelPart->getXCoordinates();
            y = mainModelPart->getYCoordinates();
            z = mainModelPart->getZCoordinates();

//             KRATOS_EXPECT_NEAR(x[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[1], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[1], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[1], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[2], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[2], 2.0e-8, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[2], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[3], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[3], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[3], 0.0, float_epsilon);

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
//                 KRATOS_EXPECT_NEAR(x[i], 0.0, float_epsilon);
//                 KRATOS_EXPECT_NEAR(y[i], 0.0, float_epsilon);
//                 KRATOS_EXPECT_NEAR(z[i], 0.0, float_epsilon);
//             }

            // None fixed
            wrapperInstance->calculate();
            mainModelPart->retrieveResults();
            x = mainModelPart->getXCoordinates();
            y = mainModelPart->getYCoordinates();
            z = mainModelPart->getZCoordinates();

//             KRATOS_EXPECT_NEAR(x[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[1], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[1], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[1], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[2], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[2], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[2], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[3], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[3], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[3], 0.0, float_epsilon);

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

//             KRATOS_EXPECT_NEAR(x[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[1], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[1], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[1], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[2], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[2], 1.0e-8, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[2], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[3], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[3], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[3], 0.0, float_epsilon);

            // Partially fixed
            mainModelPart->updateNodePos(0, x[0], y[0], z[0]);
            mainModelPart->updateNodePos(1, x[1], y[1], z[1]);
            mainModelPart->updateNodePos(2, x[2], y[2] + 1.0e-8, z[2]);

            wrapperInstance->calculate();
            mainModelPart->retrieveResults();
            x = mainModelPart->getXCoordinates();
            y = mainModelPart->getYCoordinates();
            z = mainModelPart->getZCoordinates();

//             KRATOS_EXPECT_NEAR(x[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[0], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[1], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[1], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[1], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[2], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[2], 2.0e-8, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[2], 0.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(x[3], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(y[3], 1.0, float_epsilon);
//             KRATOS_EXPECT_NEAR(z[3], 0.0, float_epsilon);

            Kratos::filesystem::remove((FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "file.mdpa"})).c_str());
            Kratos::filesystem::remove((FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "file.json"})).c_str());
        }

        KRATOS_TEST_CASE_IN_SUITE(CSharpWrapperInit, KratosCSharpWrapperApplicationFastSuite)
        {
            std::string mdpa_file = "D:\\ProyectosUnity\\ROMtrial\\Assets\\KratosMultiphysics\\Resources\\MoldeV2\\MoldeV2Kratos.mdpa";
            std::string json_file = "D:\\ProyectosUnity\\ROMtrial\\Assets\\KratosMultiphysics\\Resources\\MoldeV2\\ProjectParameters.json";

            CSharpKratosWrapper::KratosWrapper* wrapperInstance = new CSharpKratosWrapper::KratosWrapper();

            wrapperInstance->init(mdpa_file.c_str(), json_file.c_str());

            CSharpKratosWrapper::ModelPartWrapper* modelPartInstance = wrapperInstance->getRootModelPartWrapper();
            int* truanglesRoot = modelPartInstance->getTriangles();
            int triangleSize = modelPartInstance->getTrianglesCount();
            
            int nodeSize = modelPartInstance->getNodesCount();
            std::cout << " Tamano de triangulos: " << triangleSize << std::endl;
            std::cout << " Tamano de nodos: " << nodeSize << std::endl;

            int maxTriangle = 0;
            int maxNodeid = 0;

            for (int i = 0; i < triangleSize; i++) {
                if (maxTriangle < truanglesRoot[i]) {
                    maxTriangle = truanglesRoot[i];
                }
               // std::cout << "   " << truanglesRoot[i] << std::endl;
            }

       /*     for (int i = 0; i < nodeSize; i++) {
                if (maxNodeid < modelPartInstance->getNodes()[i]->Id()) {
                    maxNodeid = modelPartInstance->getNodes()[i]->Id();
                }
            }*/

            std::cout << " MaxTriangle : " << maxTriangle << std::endl;
            std::cout << "Max Node id : " << maxNodeid << std::endl;

            auto sk = modelPartInstance->getSubmodelPart("Parts_Solid_Cross");
            auto nodes = sk->getNodesCount();

            auto& mdpa_debug = modelPartInstance->getKratosModelPart();

            auto ss = modelPartInstance->getSubmodelPart("Parts_Solid_Sphere");
            int* triangles = ss->getTriangles();
            int trianglesLenght = ss->getTrianglesCount();
            auto nodess = ss->getNodesCount();

            bool hassub = modelPartInstance->hasSubmodelPart("DISPLACEMENT_MoveBC");

            int* kratosIdsSMP = modelPartInstance->getKratosNodesIdSubModelPart("DISPLACEMENT_MoveBC");
            
            int sizeIds = modelPartInstance->getKratosIdsCountSubModelPart("DISPLACEMENT_MoveBC");
            //int countids = sizeof(kratosIdsSMP)/sizeof(int);
            //std::cout << " countID:  " << sizeIds << " || " << countids << std::endl;
            //

            //for (int i = 0; i < sizeIds; i++)
            //{
            //    std::cout << kratosIdsSMP[i] << std::endl;
            //}

            //for (int i = 0; i < trianglesLenght; i++) {
            //    std::cout << triangles[i];
            //}

            system("pause");


#pragma region nodes

               int kratosIndex[] =
                {
                    393,
                    415,
                    426,
                    430,
                    435,
                    445,
                    447,
                    469,
                    479,
                    496,
                    520,
                    522,
                    538,
                    546,
                    553,
                    561,
                    581
                };

               int StaticNodesKratos[] =
               {
                   1,
                   3,
                   5,
                   6,
                   9,
                  10,
                  13,
                  19,
                  22,
                  23,
                  28,
                  29,
                  31,
                  34,
                  35,
                  38,
                  39,
                  41,
                  47,
                  49,
                  50,
                  51,
                  57,
                  58,
                  59,
                  60,
                  65,
                  67,
                  69,
                  75,
                  77,
                  80,
                  83,
                  86,
                  87,
                  95,
                  97,
                  98,
                 102,
                 103,
                 106,
                 107,
                 117,
                 122,
                 125,
                 126,
                 127,
                 131,
                 134,
                 138,
                 140,
                 141,
                 142,
                 143,
                 146,
                 147,
                 152,
                 155,
                 165,
                 168,
                 171,
                 173,
                 179,
                 186,
                 187,
                 198,
                 201,
                 202,
                 203,
                 206,
                 210,
                 212,
                 223,
                 224,
                 228,
                 238,
                 244,
                 248,
                 249,
                 265,
                 269,
                 275,
                 283,
                 286,
                 294,
                 295,
                 301,
                 307,
                 319,
                 325,
                 335,
                 345,
                 357,
                 366,
                 384,
                 389,
                 390,
                 392,
                 399,
                 410,
                 412,
                 417,
                 423,
                 433,
                 436,
                 449,
                 458,
                 459,
                 463,
                 471,
                 485,
                 501,
                 508,
                 526,
                 530,
                 535,
                 544,
                 548,
                 549,
                 552,
                 556,
                 558,
                 565,
                 583,
                 601,
                 612,
                 633,
                 634,
                 645,
                 648,
                 652,
                 664,
                 668,
                 675,
                 681,
                 698,
                 706,
                 727,
                 733,
                 740,
                 758,
                 772,
                 787,
                 791,
                 792,
                 808,
                 810,
                 819,
                 821,
                 827,
                 858,
                 864,
                 865,
                 866,
                 876,
                 903,
                 912,
                 927,
                 929,
                 956,
                 963,
                 964,
                 968,
                 985,
                 986,
                 987,
                 989,
                1000,
                1002,
                1009,
                1024,
                1036,
                1038,
                1039,
                1040,
                1044,
                1045,
                1055,
                1056,
                1060,
                1069,
                1070,
                1074,
                1077,
                1078,
                1082,
                1088,
                1091,
                1095,
                1096,
                1103,
                1106,
                1113,
                1114,
                1117,
                1121,
                1126,
                1127,
                1134,
                1139,
                1143,
                1144,
                1148,
                1149,
                1150,
                1151
               };
#pragma endregion

            int sizeMove = sizeof(kratosIndex) / sizeof(kratosIndex[0]);

            int sizeStatic = sizeof(StaticNodesKratos) / sizeof(StaticNodesKratos[0]);

            float displacement = 0;

           

            for (int i = 0; i < 3; i++) {
               
                displacement += .1f;

                for (int j = 0; j < sizeIds; j++) {

            /*        auto node = modelPartInstance->getNode(kratosIdsSMP[j]);
                    float nX = node->X();
                    float nY = node->Y();
                    float nZ = node->Z();*/
                    //std::cout << "Nodo: " << kratosIdsSMP[j] << "  " << node->Info() << "  " << nX << "  " << nY << " " << nZ << std::endl;


                    modelPartInstance->updateNodePosKratosId(kratosIdsSMP[j],0, - displacement, 0);
                }

                for (int j = 0; j < sizeStatic; j++) {
                /*    auto node = modelPartInstance->getNode(StaticNodesKratos[j]);
                    float nX = node->X();
                    float nY = node->Y();
                    float nZ = node->Z();*/
                   // std::cout << "Nodo: " << node->Info() << "  " << nX << "  " << nY << " " << nZ << std::endl;

                    modelPartInstance->updateNodePosKratosId(StaticNodesKratos[j], 0, 0, 0);
                }
         
                modelPartInstance->enableSurfaceStressResults();
                wrapperInstance->calculate();
                modelPartInstance->retrieveResults();
                float *stress = modelPartInstance->getSurfaceStress();

                std::cout << " Tamano de stress: " << stress << std::endl;
                
                //for (int i = 0; i < triangleSize; i++) {
                //    std::cout << "  Stress " << i << ": " << stress[i] << std::endl;
                //}
                //system("pause");
                //std::cout << "Calculation OK" << std::endl;
                //system("pause");

                int nodeindex = 581;
                auto noderef = modelPartInstance->getNode(nodeindex);
                Kratos::array_1d<double, 3> displacement = noderef->GetSolutionStepValue(Kratos::DISPLACEMENT);
                std::cout << i << " .:  Simulation Step :. " << " displacement: " << displacement << " Node 581:   " << noderef->Coordinates() << std::endl;
            }

            std::cout << "Simulation Finalized..." << std::endl;
            system("pause");

        }


        KRATOS_TEST_CASE_IN_SUITE(CSharpWrapperInitROM, KratosCSharpWrapperApplicationFastSuite)
        {
            std::string mdpa_file = "D:\\ProyectosUnity\\ROMtrial\\Assets\\KratosMultiphysics\\Resources\\Bunny\\StanfordBunnyKratos.mdpa";
            std::string json_file = "D:\\ProyectosUnity\\ROMtrial\\Assets\\KratosMultiphysics\\Resources\\Bunny\\ProjectParameters.json";
            std::string rom_json_file = "D:\\ProyectosUnity\\ROMtrial\\Assets\\KratosMultiphysics\\Resources\\Bunny\\RomParameters.json";

            CSharpKratosWrapper::KratosWrapper* wrapperInstance = new CSharpKratosWrapper::KratosWrapper();

            wrapperInstance->initRom(mdpa_file.c_str(), json_file.c_str(), rom_json_file.c_str());


            CSharpKratosWrapper::ModelPartWrapper* modelPartInstance = wrapperInstance->getRootModelPartWrapper();

            std::ofstream myfile;
            myfile.open("nodeRomResults.txt");

            int kratosIndex[] = {
              673,
              683,
              689,
              697,
              723,
              776,
              783,
              787,
              853,
              855,
              887,
              889,
              912,
             1042,
             1048,
             1072,
             1075,
             1096,
             1143,
             1178,
             1193,
             1200,
             1216,
             1217,
             1230,
             1249,
             1298,
             1300,
             1350,
             1354,
             1378,
             1381,
             1383,
             1413,
             1440,
             1471,
             1479,
             1502,
             1529,
             1541,
             1551,
             1564,
             1574,
             1575,
             1613,
             1625,
             1643,
             1660,
             1680,
             1682,
             1686,
             1696,
             1740,
             1747,
             1763,
             1765,
             1776,
             1800,
             1817,
             1830,
             1832,
             1858,
             1864,
             1899,
             1900,
             1922,
             1926,
             1943,
             1948,
             1949,
             1955,
             1958,
             1963,
             2049,
             2051,
             2070,
             2075,
             2082,
             2083,
             2088,
             2093,
             2129,
             2133,
             2153,
             2198,
             2200,
             2208,
             2219,
             2228,
             2255,
             2266,
             2273,
             2298,
             2307,
             2340,
             2348,
             2368,
             2372,
             2386,
             2403,
             2412,
             2414,
             2433,
             2438,
             2488,
             2519,
             2525,
             2537,
             2551,
             2566,
             2617,
             2621,
             2625,
             2632,
             2646,
             2663,
             2665,
             2679,
             2695,
             2731,
             2732,
             2750,
             2765,
             2776,
             2818,
             2822,
             2827,
             2849,
             2893,
             2911,
             2926,
             2989,
             2991,
             2992,
             3027,
             3033,
             3113,
             3168,
             3194,
             3227,
             3268,
             3344,
             3366,
             3440,
             3456,
             3554,
             3598,
             3607,
             3624
            };

            int sizeMove = sizeof(kratosIndex) / sizeof(kratosIndex[0]);

            int StaticNodesKratos[] = {
                    5987,
                    5988,
                    5997,
                    6006,
                    6014,
                    6026,
                    6032,
                    6033,
                    6044,
                    6050,
                    6055,
                    6071,
                    6078,
                    6079,
                    6080,
                    6095,
                    6097,
                    6099,
                    6100,
                    6106,
                    6110,
                    6111,
                    6113,
                    6117,
                    6119,
                    6123,
                    6136,
                    6140,
                    6144,
                    6145,
                    6147,
                    6151,
                    6153,
                    6156,
                    6160,
                    6161,
                    6162,
                    6165,
                    6174,
                    6175,
                    6176,
                    6179,
                    6180,
                    6191,
                    6192,
                    6202,
                    6203,
                    6204,
                    6206,
                    6211,
                    6217,
                    6223,
                    6226,
                    6227,
                    6228,
                    6233,
                    6239,
                    6242,
                    6245,
                    6246,
                    6248,
                    6249,
                    6251,
                    6259,
                    6271,
                    6272,
                    6276,
                    6289,
                    6296,
                    6300,
                    6301,
                    6305,
                    6306,
                    6309,
                    6311,
                    6316,
                    6318,
                    6323,
                    6329,
                    6331,
                    6336,
                    6337,
                    6343,
                    6353,
                    6355,
                    6362,
                    6370,
                    6377,
                    6380,
                    6381,
                    6383,
                    6386,
                    6387,
                    6395,
                    6397,
                    6400,
                    6407,
                    6413,
                    6416,
                    6418,
                    6419,
                    6423,
                    6427,
                    6431,
                    6433,
                    6437,
                    6439,
                    6441,
                    6444,
                    6453,
                    6455,
                    6460,
                    6463,
                    6466,
                    6471,
                    6494,
                    6495,
                    6496,
                    6498,
                    6517,
                    6520,
                    6525,
                    6531,
                    6533,
                    6539,
                    6548,
                    6553,
                    6560,
                    6561,
                    6566,
                    6567,
                    6568,
                    6575,
                    6578,
                    6591,
                    6596,
                    6602,
                    6604,
                    6606,
                    6608,
                    6619,
                    6631,
                    6633,
                    6643,
                    6647,
                    6648,
                    6660,
                    6669,
                    6671,
                    6674,
                    6679,
                    6682,
                    6683,
                    6693,
                    6702,
                    6703,
                    6705,
                    6720,
                    6722,
                    6723,
                    6724,
                    6736,
                    6740,
                    6741,
                    6742,
                    6743,
                    6762,
                    6766,
                    6769,
                    6772,
                    6773,
                    6782,
                    6789,
                    6802,
                    6806,
                    6809,
                    6835,
                    6838,
                    6843,
                    6845,
                    6851,
                    6858,
                    6862,
                    6863,
                    6864,
                    6870,
                    6879,
                    6880,
                    6886,
                    6889,
                    6890,
                    6903,
                    6909,
                    6913,
                    6935,
                    6944,
                    6956,
                    6957,
                    6961,
                    6963,
                    6981,
                    6983,
                    6984,
                    6986,
                    6988,
                    6993,
                    6994,
                    7004,
                    7005,
                    7016,
                    7022,
                    7026,
                    7030,
                    7031,
                    7043,
                    7045,
                    7053,
                    7055,
                    7069,
                    7072,
                    7076,
                    7079,
                    7085,
                    7088,
                    7090,
                    7092,
                    7093,
                    7097,
                    7100,
                    7104,
                    7113,
                    7117,
                    7119,
                    7130,
                    7136,
                    7141,
                    7144,
                    7150,
                    7151,
                    7155,
                    7175,
                    7179,
                    7188,
                    7194,
                    7201,
                    7203,
                    7213,
                    7217,
                    7219,
                    7231,
                    7235,
                    7247,
                    7258,
                    7262,
                    7263,
                    7271,
                    7275,
                    7280,
                    7281,
                    7283,
                    7288,
                    7292,
                    7297,
                    7301,
                    7302,
                    7310,
                    7317,
                    7328,
                    7329,
                    7330,
                    7332,
                    7333,
                    7336,
                    7341,
                    7344,
                    7349,
                    7359,
                    7362,
                    7371,
                    7372,
                    7387,
                    7398,
                    7400,
                    7401,
                    7413,
                    7423,
                    7426,
                    7433,
                    7437,
                    7443,
                    7444,
                    7448,
                    7456,
                    7463,
                    7468,
                    7474,
                    7479,
                    7481,
                    7482,
                    7488,
                    7501,
                    7505,
                    7510,
                    7512,
                    7515,
                    7516,
                    7520,
                    7523,
                    7525,
                    7529,
                    7530,
                    7532,
                    7533,
                    7544,
                    7551,
                    7556,
                    7565,
                    7589,
                    7594,
                    7606,
                    7609,
                    7610,
                    7618,
                    7620,
                    7623,
                    7625,
                    7630,
                    7632,
                    7637,
                    7640,
                    7649,
                    7661,
                    7663,
                    7664,
                    7669,
                    7670,
                    7679,
                    7689,
                    7697,
                    7711,
                    7712,
                    7716,
                    7723,
                    7724,
                    7725,
                    7726,
                    7731,
                    7735,
                    7740,
                    7744,
                    7746,
                    7750,
                    7753,
                    7767,
                    7773,
                    7780,
                    7787,
                    7790,
                    7796,
                    7801,
                    7803,
                    7805,
                    7810,
                    7825,
                    7829,
                    7833,
                    7834,
                    7836,
                    7851,
                    7852,
                    7859,
                    7865,
                    7866,
                    7869,
                    7872,
                    7874,
                    7876,
                    7879,
                    7881,
                    7882,
                    7885,
                    7887,
                    7889,
                    7896,
                    7897,
                    7901,
                    7904,
                    7907,
                    7910,
                    7917,
                    7920,
                    7924,
                    7925,
                    7929,
                    7933,
                    7935,
                    7940,
                    7951,
                    7953,
                    7954,
                    7959,
                    7961,
                    7967,
                    7968,
                    7970,
                    7971,
                    7972,
                    7974,
                    7979,
                    7985,
                    7987,
                    7988,
                    7989,
                    7993,
                    7995,
                    7997,
                    8002,
                    8003,
                    8005,
                    8006,
                    8007,
                    8008,
                    8009,
                    8019,
                    8024,
                    8028,
                    8031,
                    8040,
                    8046,
                    8049,
                    8055,
                    8060,
                    8061,
                    8066,
                    8069,
                    8071,
                    8073,
                    8079,
                    8093,
                    8097,
                    8103,
                    8110,
                    8111,
                    8115,
                    8116,
                    8124,
                    8128,
                    8131,
                    8132,
                    8134,
                    8135,
                    8141,
                    8145,
                    8153,
                    8156,
                    8162,
                    8165,
                    8173,
                    8175,
                    8180,
                    8186,
                    8187,
                    8189,
                    8204,
                    8206,
                    8208,
                    8218,
                    8230,
                    8231,
                    8232,
                    8247,
                    8265,
                    8272,
                    8283,
                    8284,
                    8286,
                    8288,
                    8305,
                    8311,
                    8325,
                    8340,
                    8348,
                    8358,
                    8373,
                    8376,
                    8377,
                    8390,
                    8408,
                    8442,
                    8452,
                    8457,
                    8459,
                    8463,
                    8497,
                    8518,
                    8541,
                    8553,
                    8555,
                    8556,
                    8567,
                    8593,
                    8603,
                    8620,
                    8631,
                    8639,
                    8648,
                    8664,
                    8671,
                    8673,
                    8687,
                    8692,
                    8693,
                    8707,
                    8713,
                    8723,
                    8731,
                    8733,
                    8746,
                    8762,
                    8777,
                    8781,
                    8783,
                    8803,
                    8811,
                    8812,
                    8830,
                    8836,
                    8838,
                    8854,
                    8864
            };


            int sizeStatic = sizeof(StaticNodesKratos) / sizeof(StaticNodesKratos[0]);

            for (int i = 0; i < 30; i++) {

                float nodeFacePressure = (-15000000 * 13) + (i * 15000000);


                for (int i = 0; i < sizeStatic; i++) {

                    auto snode = modelPartInstance->getNode(StaticNodesKratos[i]);

                    float x = snode->X();
                    float y = snode->Y();
                    float z = snode->Z();

                    modelPartInstance->updateNodePos(modelPartInstance->getIdTranslator()->getSurfaceId(StaticNodesKratos[i]), x, y, z);

                }

                for (int i = 0; i < sizeMove; ++i) {
                    int index = modelPartInstance->getIdTranslator()->getSurfaceId(kratosIndex[i]);
                    modelPartInstance->updateNodePressure(index, nodeFacePressure);
                }

                wrapperInstance->calculate();
                system("pause");
                int nodeindex = 673;
                auto noderef = modelPartInstance->getNode(nodeindex);
                Kratos::array_1d<double, 3> displacement = noderef->GetSolutionStepValue(Kratos::DISPLACEMENT);
                auto force = noderef->GetSolutionStepValue(Kratos::POSITIVE_FACE_PRESSURE);
                std::cout << "Node: " << nodeindex << "  " << noderef->Coordinates() << " displacement: " << displacement;
                myfile << i << "  Node: " << nodeindex << "  " << noderef->Coordinates() << " displacement: " << displacement << " Pressure: " << nodeFacePressure << " FacePressure " << force << "\n";

                if (i == 29) {
                    myfile.close();
                }

            }
        }

        KRATOS_TEST_CASE_IN_SUITE(CSharpWrapperInitConvDiff, KratosCSharpWrapperApplicationFastSuite)
        {
            std::string mdpa_file = "D:\\Descargas\\ThermalProblem_Sliders\\Siemens.mdpa";
            std::string json_file = "D:\\Descargas\\ThermalProblem_Sliders\\ProjectParameters.json";
            std::string rom_json_file = "D:\\Descargas\\ThermalProblem_Sliders\\RomParameters.json";

            CSharpKratosWrapper::KratosWrapper* wrapperInstance = new CSharpKratosWrapper::KratosWrapper();

            wrapperInstance->initConvDiff(mdpa_file.c_str(), json_file.c_str());

            CSharpKratosWrapper::ModelPartWrapper* modelPartInstance = wrapperInstance->getRootModelPartWrapper();

        }
    } // namespace Testing
}  // namespace Kratos.
