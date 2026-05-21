//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License: BSD License
//           Kratos default license: kratos/license.txt
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

#include "test_creation_utility.h"

#include "custom_elements/shell_5p_element.h"
#include "custom_utilities/director_utilities.h"

namespace Kratos
{
namespace Testing
{
    ///@name Type Definitions
    ///@{

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    ///@}
    ///@name Operations
    ///@{

    typename Shell5pElement::Pointer GetShell5pElement(
        ModelPart& rModelPart, SizeType PolynomialDegree, IntegrationPoint<3> IntegrationPoint)
    {
        // Set the element properties
        auto p_elem_prop = rModelPart.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 200000000);
        p_elem_prop->SetValue(POISSON_RATIO, 0.0);
        p_elem_prop->SetValue(THICKNESS, 0.01);

        auto p_quadrature_point = TestCreationUtility::GetQuadraturePointGeometry(
            rModelPart, PolynomialDegree, IntegrationPoint);


        return Kratos::make_intrusive<Shell5pElement>(1, p_quadrature_point, p_elem_prop);
    }

    Parameters GetDirectorParametersSimpleTest()
    {
        return Parameters(R"(
        {
            "model_part_name": "ModelPart",
            "brep_ids" : [1] ,
            "linear_solver_settings" : {
                "solver_type": "skyline_lu_factorization"
            }
        })");
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=3.
    KRATOS_TEST_CASE_IN_SUITE(IgaDirectorUtility, KratosIgaFast5PSuite)
    {
        Model current_model;
        auto& r_model_part = current_model.CreateModelPart("ModelPart");

        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(DIRECTORINC);

        IntegrationPoint<3> integration_point(0.0694318442029737, 0.211324865405187, 0.0, 0.086963711284364);
        GetShell5pElement(r_model_part, 3, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);
        TestCreationUtility::AddDirectorInc2DDofs(r_model_part);

        DirectorUtilities(r_model_part, GetDirectorParametersSimpleTest()).ComputeDirectors();
        const double tolerance = 1.0e-8;
        array_1d<double, 3> director = ZeroVector(3);
        director[2] = 1.0;

        KRATOS_EXPECT_VECTOR_NEAR(r_model_part.GetNode(4).GetValue(DIRECTOR), director, tolerance);
        KRATOS_EXPECT_VECTOR_NEAR(r_model_part.GetNode(4).GetValue(DIRECTOR), director, tolerance);
        KRATOS_EXPECT_VECTOR_NEAR(r_model_part.GetNode(8).GetValue(DIRECTOR), director, tolerance);
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=3.
    KRATOS_TEST_CASE_IN_SUITE(IgaShell5pElementP3, KratosIgaFast5PSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");

        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
        const auto& r_process_info = r_model_part.GetProcessInfo();

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(DIRECTORINC);

        IntegrationPoint<3> integration_point(0.0694318442029737, 0.211324865405187, 0.0, 0.086963711284364);
        auto p_shell_5p_element = GetShell5pElement(r_model_part, 3, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);
        TestCreationUtility::AddDirectorInc2DDofs(r_model_part);

        DirectorUtilities(r_model_part, GetDirectorParametersSimpleTest()).ComputeDirectors();

        p_shell_5p_element->Check(r_process_info);
        p_shell_5p_element->Initialize(r_process_info);

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_5p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        const double tolerance = 1.0e-8;

        const std::array<double, 40> expected_LHS_row_0{637725.1237693131, 143581.4724979368, 0, 0, 0, 64285.87864500261, -122155.5804409024, 0, 0, 0, -1057.606121067305, -20626.57435529137, 0, 0, 0, -171.9000696345915, -799.3177017429664, 0, 0, 0, -545148.4726534408, 38472.53960389353, 0, 0, 0, -143047.8132845628, -32731.48913009484, 0, 0, 0, -12241.73636779629, -5526.873941120825, 0, 0, 0, -343.4739178139834, -214.1765326779287, 0, 0, 0};
        const std::array<double, 40> expected_LHS_row_1{ 143581.4724979368,1165930.893369905,0,0,0,32138.83808555151,221748.1263942097,0,0,0,2397.9531052289,13618.07884900591,0,0,0,59.63894400465794,265.8938341073225,0,0,0,-143581.4724979368,-1119642.567811969,0,0,0,-32138.83808555151,-261129.0937139899,0,0,0,-2397.953105228899,-20267.75009343771,0,0,0,-59.63894400465793,-523.5808278316099,0,0,0 };
        const std::array<double, 40> expected_LHS_row_2{ 0,0,601218.6723797393,-11323.91371476317,-44537.44868967543,0,0,95344.66834640411,-2534.710245284323,-9969.126428910287,0,0,4186.824242646201,-189.1205988018452,-743.8195996068691,0,0,31.33125482424369,-4.703575218996611,-18.49936737869886,0,0,-554930.3468218031,-3034.233535030504,-11933.77340934109,0,0,-134725.6356661842,-679.1735632708269,-2671.219375870272,0,0,-10836.495487078,-50.67471172104472,-199.3058610291014,0,0,-289.0182485485311,-1.260319181469182,-4.956890549609012 };
        const std::array<double, 40> expected_RHS{0,0,0,0,0,0,0,0,0,0
                                                 ,0,0,0,0,0,0,0,0,0,0
                                                 ,0,0,0,0,0,0,0,0,0,0
                                                 ,0,0,0,0,0,0,0,0,0,0};

        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
            KRATOS_EXPECT_NEAR(left_hand_side_matrix(0,i), expected_LHS_row_0[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
            KRATOS_EXPECT_NEAR(left_hand_side_matrix(1,i), expected_LHS_row_1[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
            KRATOS_EXPECT_NEAR(left_hand_side_matrix(2,i), expected_LHS_row_2[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
            KRATOS_EXPECT_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=4.
    KRATOS_TEST_CASE_IN_SUITE(IgaShell5pElementP4, KratosIgaFast5PSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
        const auto& r_process_info = r_model_part.GetProcessInfo();

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(DIRECTORINC);

        IntegrationPoint<3> integration_point(0.0469100770306680, 0.211324865405187, 0, 0.0592317212640473);
        auto p_shell_5p_element = GetShell5pElement(r_model_part, 4, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);
        TestCreationUtility::AddDirectorInc2DDofs(r_model_part);

        DirectorUtilities(r_model_part, GetDirectorParametersSimpleTest()).ComputeDirectors();

        p_shell_5p_element->Initialize(r_process_info);

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_5p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        const double tolerance = 1.0e-8;

        const std::array<double, 50> expected_LHS_row_0{ 491667.1614141508,133490.278045005,0,0,0,4078.075758025352,-113779.5275428051,0,0,0,-6544.238833106851,-18740.608229045,0,0,0,-439.346498780415,-954.2258146033695,0,0,0,-8.169843782843573,-15.91645855150754,0,0,0,-379618.7388442269,35768.61219956507,0,0,0,-99581.7841076655,-30487.13252028928,0,0,0,-9186.162068129439,-5021.530840640658,0,0,0,-361.6069404337981,-255.684036419903,0,0,0,-5.190036050463355,-4.264802215239878,0,0,0 };
        const std::array<double, 50> expected_LHS_row_1{ 133490.278045005,850779.6157010947,0,0,0,26281.00066959985,121138.251571485,0,0,0,1940.284546309754,5520.786405805647,0,0,0,63.66583420603015,68.84509465384181,0,0,0,0.7833912401261195,-0.534780027131341,0,0,0,-133490.278045005,-794755.4044161328,0,0,0,-26281.00066959985,-168890.1057463051,0,0,0,-1940.284546309754,-13385.98685642379,0,0,0,-63.66583420603014,-469.3218142609483,0,0,0,-0.7833912401261196,-6.145159889522122,0,0,0 };
        const std::array<double, 50> expected_LHS_row_2{ 0,0,447482.2590384154,-10528.04630042432,-31807.0597047671,0,0,41738.77577650348,-2072.717174038302,-6262.039225936466,0,0,-341.1508091004015,-153.0254175713081,-462.3164121952369,0,0,-123.5004680421911,-5.021166036153184,-15.16981625480087,0,0,-2.901541269991639,-0.06178411917468645,-0.186660574176663,0,0,-391458.0477534534,-2820.981504076165,-8522.675961500841,0,0,-89490.62995132357,-555.3828929216814,-1677.908353561687,0,0,-7524.049641517747,-41.0030370596673,-123.8773092953677,0,0,-276.9762515649156,-1.345417384449826,-4.064740014802409,0,0,-3.778398646661828,-0.01655500483792547,-0.05001555010936651 };
        const std::array<double, 50> expected_RHS{ 0,0,0,0,0,0,0,0,0,0
                                                  ,0,0,0,0,0,0,0,0,0,0
                                                  ,0,0,0,0,0,0,0,0,0,0
                                                  ,0,0,0,0,0,0,0,0,0,0
                                                  ,0,0,0,0,0,0,0,0,0,0 };

        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
            KRATOS_EXPECT_NEAR(left_hand_side_matrix(0,i), expected_LHS_row_0[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
            KRATOS_EXPECT_NEAR(left_hand_side_matrix(1,i), expected_LHS_row_1[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
            KRATOS_EXPECT_NEAR(left_hand_side_matrix(2,i), expected_LHS_row_2[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
            KRATOS_EXPECT_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=5.
    KRATOS_TEST_CASE_IN_SUITE(IgaShell5pElementP5, KratosIgaFast5PSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
        const auto& r_process_info = r_model_part.GetProcessInfo();

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(DIRECTORINC);

        IntegrationPoint<3> integration_point(0.0337652428984240, 0.211324865405187, 0, 0.0428311230947926);
        auto p_shell_5p_element = GetShell5pElement(r_model_part, 5, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);
        TestCreationUtility::AddDirectorInc2DDofs(r_model_part);

        DirectorUtilities(r_model_part, GetDirectorParametersSimpleTest()).ComputeDirectors();

        p_shell_5p_element->Initialize(r_process_info);

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_5p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        const double tolerance = 1.0e-8;

        const std::array<double, 60> expected_LHS_row_0{ 405000.2900314076,123985.6767727672,0,0,0,-33973.98814693153,-106654.8719703566,0,0,0,-9694.630768219549,-16422.36278641331,0,0,0,-594.5854330008631,-887.2782390413756,0,0,0,-14.85851893613687,-20.97888398219133,0,0,0,-0.1350847109561106,-0.1848929735518754,0,0,0,-276681.9120580379,33221.86196428906,0,0,0,-76407.95701336297,-28578.08681330184,0,0,0,-7301.614363424889,-4400.358846430367,0,0,0,-323.6988995799197,-237.7454876128451,0,0,0,-6.853475822902734,-5.621275021134379,0,0,0,-0.05626937995498505,-0.04954192294941391,0,0,0 };
        const std::array<double, 60> expected_LHS_row_1{ 123985.6767727671,658198.3090341121,0,0,0,21663.50600301304,62635.26689033574,0,0,0,1514.070026661864,717.512338770068,0,0,0,52.90944238935797,-102.828836610041,0,0,0,0.9244648677593761,-4.031472288912596,0,0,0,0.006461117359220246,-0.04379510180535657,0,0,0,-123985.6767727671,-594039.1200474274,0,0,0,-21663.50600301305,-117826.239470483,0,0,0,-1514.070026661865,-9215.634904592293,0,0,0,-52.90944238935805,-356.3133296803508,0,0,0,-0.9244648677593776,-6.824525090607215,0,0,0,-0.00646111735922026,-0.05188194365019132,0,0,0 };
        const std::array<double, 60> expected_LHS_row_2{ 0,0,354399.5330218398,-9778.442031659108,-23959.85405612183,0,0,9553.759581134736,-1708.546851272185,-4186.406492157968,0,0,-2992.372809816492,-119.4109382063518,-292.5894168892806,0,0,-232.4714232036346,-4.172836159776343,-10.22458844309223,0,0,-6.296663741683155,-0.07291024540083026,-0.1786500173896845,0,0,-0.05962660425382239,-0.0005095722602915917,-0.001248591232438192,0,0,-290240.344035155,-2620.125645617602,-6420.0235451054,0,0,-64744.732161282,-457.8037490291192,-1121.744238762136,0,0,-5505.749756005726,-31.99606445983464,-78.3990979693758,0,0,-226.6707430867568,-1.118108079159458,-2.739670216267152,0,0,-4.559333637836649,-0.01953624137510743,-0.04786912788737195,0,0,-0.03605044120172546,-0.0001365394756304337,-0.0003345590124083934 };
        const std::array<double, 60> expected_RHS{ 0,0,0,0,0,0,0,0,0,0
                                                  ,0,0,0,0,0,0,0,0,0,0
                                                  ,0,0,0,0,0,0,0,0,0,0
                                                  ,0,0,0,0,0,0,0,0,0,0
                                                  ,0,0,0,0,0,0,0,0,0,0
                                                  ,0,0,0,0,0,0,0,0,0,0 };

        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
            KRATOS_EXPECT_NEAR(left_hand_side_matrix(0,i), expected_LHS_row_0[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
            KRATOS_EXPECT_NEAR(left_hand_side_matrix(1,i), expected_LHS_row_1[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
            KRATOS_EXPECT_NEAR(left_hand_side_matrix(2,i), expected_LHS_row_2[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
            KRATOS_EXPECT_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=3 (with initial displacement)
    KRATOS_TEST_CASE_IN_SUITE(IgaShell5pElementP3Disp, KratosIgaFast5PSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
        const auto& r_process_info = r_model_part.GetProcessInfo();

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(DIRECTORINC);

        IntegrationPoint<3> integration_point(0.0694318442029737, 0.211324865405187, 0.0, 0.086963711284364);
        auto p_shell_5p_element = GetShell5pElement(r_model_part, 3, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);
        TestCreationUtility::AddDirectorInc2DDofs(r_model_part);

        DirectorUtilities(r_model_part, GetDirectorParametersSimpleTest()).ComputeDirectors();

        p_shell_5p_element->Initialize(r_process_info);

        array_1d<double, 3> delta = ZeroVector(3);
        delta[2] = 0.002;

        for (auto& node : p_shell_5p_element->GetGeometry()) // p_shell_3p_element->GetGeometry()[i].
        {
            if(node.Id()==3 || node.Id()==7)
            {
                node.FastGetSolutionStepValue(DISPLACEMENT) = delta/2;
                node.Coordinates() += node.FastGetSolutionStepValue(DISPLACEMENT);
            }
            if(node.Id()==4 || node.Id()==8)
            {
                node.FastGetSolutionStepValue(DISPLACEMENT) = delta;
                node.Coordinates() += node.FastGetSolutionStepValue(DISPLACEMENT);
            }
        }

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_5p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        const double tolerance = 1.0e-8;

        const std::array<double, 40> expected_LHS_row_0{ 637725.1296726804,143581.4724979368,256.4475891677992,-4.553670952982031,0,64285.87362256309,-122155.5804409024,25.85119823820547,-1.019279792209065,0,-1057.606969131077,-20626.57435529137,-0.4252937981081139,-0.07605082474725849,0,-171.9001024986228,-799.3177017429664,-0.06912595535676944,-0.001891442692819822,0,-545148.4710716383,38472.53960389353,-219.2198587444297,-1.220152454448596,0,-143047.8146303215,-32731.48913009484,-57.52363437671553,-0.2731151971837814,0,-12241.73659503429,-5526.873941120825,-4.922753803697947,-0.02037775707474875,0,-343.4739266198741,-214.1765326779287,-0.1381207276965891,-0.000506810542070819,0 };
        const std::array<double, 40> expected_LHS_row_1{ 143581.4724979368,1165930.899273272,57.7382340743063,0,-4.553670952982031,32138.83808555151,221748.1213717702,12.92394989392847,0,-1.019279792209065,2397.9531052289,13618.07800094214,0.9642858182200725,0,-0.07605082474725849,59.63894400465794,265.8938012432913,0.02398253234890642,0,-0.001891442692819822,-143581.4724979368,-1119642.566230166,-57.7382340743063,0,-1.220152454448596,-32138.83808555151,-261129.0950597485,-12.92394989392847,0,-0.2731151971837814,-2397.953105228899,-20267.75032067571,-0.9642858182200726,0,-0.02037775707474875,-59.63894400465793,-523.5808366375006,-0.02398253234890642,0,-0.000506810542070819 };
        const std::array<double, 40> expected_LHS_row_2{ 256.4475891677992,57.73823407430628,601218.7814080479,-11323.91371476317,-44537.44868967543,25.85119823820546,-49.12226747835398,95344.67371947391,-2534.710245284323,-9969.126428910287,-0.4252937981081139,-8.294537989878952,4186.823223559578,-189.1205988018452,-743.8195996068691,-0.06912595535676944,-0.3214286060733575,31.33119416268029,-4.703575218996611,-18.49936737869886,-219.2198587444297,15.47091319260948,-554930.4333946023,-3034.233535030504,-11933.77340934109,-57.52363437671553,-13.16227190121054,-134725.6601438485,-679.1735632708269,-2671.219375870272,-4.922753803697947,-2.222514755977332,-10836.4976938967,-50.67471172104472,-199.3058610291014,-0.1381207276965891,-0.08612653542161754,-289.0183128967335,-1.260319181469182,-4.956890549609012 };
        const std::array<double, 40> expected_RHS{ 0.00288127048596413,5.104804942162921e-13,7.165046390521992,-2.22252097587773,-3.668470065476827e-17,-0.002451313962638624,-4.343042315271281e-13,-6.095844996524769,-0.4974823042471772,-8.211391303132552e-18,-0.0004139164951875832,-7.333441904232834e-14,-1.029313598489092,-0.03711830630249979,-6.12670913029941e-19,-1.604002564640936e-05,-2.841843646835727e-15,-0.03988779550812917,-0.0009231609184914462,-1.523759834827438e-20,0.0007720340973440014,1.367828361770958e-13,1.919868394071895,-0.595522700647668,-9.829635915022635e-18,-0.0006568275972539793,-1.163714681071132e-13,-1.633376744004111,-0.1332999816718051,-2.200235668410309e-18,-0.0001109085906619495,-1.964989835979741e-14,-0.2758037474735249,-0.009945820198165862,-1.641646763724112e-19,-4.297911919587286e-06,-7.614697101851492e-16,-0.01068790259426097,-0.0002473602225937563,-4.082902172009926e-21};

        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_EXPECT_NEAR(left_hand_side_matrix(0,i), expected_LHS_row_0[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_EXPECT_NEAR(left_hand_side_matrix(1,i), expected_LHS_row_1[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_EXPECT_NEAR(left_hand_side_matrix(2,i), expected_LHS_row_2[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          KRATOS_EXPECT_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=4 (with initial displacement)
    KRATOS_TEST_CASE_IN_SUITE(IgaShell5pElementP4Disp, KratosIgaFast5PSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
        const auto& r_process_info = r_model_part.GetProcessInfo();

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(DIRECTORINC);

        IntegrationPoint<3> integration_point(0.0469100770306680, 0.211324865405187, 0, 0.0592317212640473);
        auto p_shell_5p_element = GetShell5pElement(r_model_part, 4, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);
        TestCreationUtility::AddDirectorInc2DDofs(r_model_part);

        DirectorUtilities(r_model_part, GetDirectorParametersSimpleTest()).ComputeDirectors();

        p_shell_5p_element->Initialize(r_process_info);

        array_1d<double, 3> delta = ZeroVector(3);
        delta[2] = 0.002;

        for (auto& node : p_shell_5p_element->GetGeometry()) // p_shell_3p_element->GetGeometry()[i].
        {
            if(node.Id()==3 || node.Id()==8)
            {
                node.FastGetSolutionStepValue(DISPLACEMENT) = delta/4;
                node.Coordinates() += node.FastGetSolutionStepValue(DISPLACEMENT);
            }
            if(node.Id()==4 || node.Id()==9)
            {
                node.FastGetSolutionStepValue(DISPLACEMENT) = delta/1.6;
                node.Coordinates() += node.FastGetSolutionStepValue(DISPLACEMENT);
            }
            if(node.Id()==5 || node.Id()==10)
            {
                node.FastGetSolutionStepValue(DISPLACEMENT) = delta;
                node.Coordinates() += node.FastGetSolutionStepValue(DISPLACEMENT);
            }
        }

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_5p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        const double tolerance = 1.0e-8;

        const std::array<double, 50> expected_LHS_row_0{ 491667.1647522034,133490.278045005,135.1390441188093,-2.893726132503429,0,4078.072912858971,-113779.5275428051,1.120894993675235,-0.5697045473253241,0,-6544.239301733763,-18740.608229045,-1.798741607732272,-0.04206038206210171,0,-439.3465226417482,-954.2258146033695,-0.1207582497707606,-0.001380111652232144,0,-8.169844180849909,-15.91645855150754,-0.002245553427318444,-1.69819086208132e-05,0,-379618.7379497983,35768.61219956507,-104.3415495748018,-0.7753715803211292,0,-99581.78487002554,-30487.13252028928,-27.37092930357384,-0.1526518733801588,0,-9186.162193697643,-5021.530840640658,-2.524897447771153,-0.01127004540688462,0,-361.606946827423,-255.684036419903,-0.0993908483462614,-0.0003697998026803852,0,-5.190036157108832,-4.264802215239878,-0.001426527061080211,-4.550288700886015e-06,0 };
        const std::array<double, 50> expected_LHS_row_1{ 133490.278045005,850779.6190391475,36.69097712824559,0,-2.893726132503429,26281.00066959985,121138.2487263186,7.223564207054831,0,-0.5697045473253241,1940.284546309754,5520.785937178736,0.5333042746898634,0,-0.04206038206210171,63.66583420603015,68.84507079250869,0.0174991145491253,0,-0.001380111652232144,0.7833912401261195,-0.5347804251376762,0.0002153219732169923,0,-1.69819086208132e-05,-133490.278045005,-794755.4035217044,-36.69097712824559,0,-0.7753715803211292,-26281.00066959985,-168890.1065086651,-7.223564207054831,0,-0.1526518733801588,-1940.284546309754,-13385.98698199199,-0.5333042746898634,0,-0.01127004540688462,-63.66583420603014,-469.3218206545732,-0.0174991145491253,0,-0.0003697998026803852,-0.7833912401261196,-6.145159996167598,-0.0002153219732169923,0,-4.550288700886015e-06 };
        const std::array<double, 50> expected_LHS_row_2{ 135.1390441188093,36.69097712824558,447482.299520623,-10528.04630042432,-31807.0597047671,1.120894993675232,-31.27330397295447,41738.77323942495,-2072.717174038302,-6262.039225936466,-1.798741607732271,-5.151021017946189,-341.1517721272649,-153.0254175713081,-462.3164121952369,-0.1207582497707606,-0.2622773587076504,-123.5005250949914,-5.021166036153184,-15.16981625480087,-0.002245553427318444,-0.004374778637281324,-2.901542285208086,-0.06178411917468645,-0.186660574176663,-104.3415495748018,9.831317691022198,-391458.0755382176,-2820.981504076165,-8522.675961500841,-27.37092930357384,-8.379656544206162,-89490.63823682428,-555.3828929216814,-1677.908353561687,-2.524897447771152,-1.38021192195442,-7524.050461076223,-41.0030370596673,-123.8773092953677,-0.09939084834626137,-0.07027700645868278,-276.9762852769895,-1.345417384449826,-4.064740014802409,-0.00142652706108021,-0.001172218402924458,-3.778399145400814,-0.01655500483792547,-0.05001555010936651 };
        const std::array<double, 50> expected_RHS{ 0.001222175032142404,4.816042666670712e-13,4.446556340576824,-1.059491929987973,8.632006757436571e-18,-0.001041712546895015,-4.104921101857641e-13,-3.789991953217198,-0.2085882846993699,1.699432937698264e-18,-0.0001715803110911822,-6.761209142137841e-14,-0.6242489832782576,-0.01539974182990685,1.254664351620343e-19,-8.73644868549722e-06,-3.442641894280185e-15,-0.03178522742185168,-0.0005053059934035036,4.116883410033861e-21,-1.457237073918238e-07,-5.742316564844739e-17,-0.0005301766595157847,-6.217656514709396e-06,5.065716078743543e-23,0.0003274808110271593,1.290454743248247e-13,1.191451180557047,-0.2838900070275687,2.312939239715121e-18,-0.0002791260360490991,-1.099910294236229e-13,-1.015525283185001,-0.05589106243578924,4.553616832471e-19,-4.597480582077843e-05,-1.811660529493758e-14,-0.1672670109453585,-0.004126348386971319,3.361862997887889e-20,-2.340924370874805e-06,-9.224531154019257e-16,-0.008516826018924731,-0.0001353963328630749,1.103115585051659e-21,-3.904654972453512e-08,-1.538649086234006e-17,-0.0001420604077630848,-1.666016041930491e-06,1.357354532384685e-23 };

        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_EXPECT_NEAR(left_hand_side_matrix(0,i), expected_LHS_row_0[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_EXPECT_NEAR(left_hand_side_matrix(1,i), expected_LHS_row_1[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_EXPECT_NEAR(left_hand_side_matrix(2,i), expected_LHS_row_2[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          KRATOS_EXPECT_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=5 (with initial displacement)
    KRATOS_TEST_CASE_IN_SUITE(IgaShell5pElementP5Disp, KratosIgaFast5PSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
        const auto& r_process_info = r_model_part.GetProcessInfo();

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(DIRECTORINC);

        IntegrationPoint<3> integration_point(0.0337652428984240, 0.211324865405187, 0, 0.0428311230947926);
        auto p_shell_5p_element = GetShell5pElement(r_model_part, 5, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);
        TestCreationUtility::AddDirectorInc2DDofs(r_model_part);

        DirectorUtilities(r_model_part, GetDirectorParametersSimpleTest()).ComputeDirectors();

        p_shell_5p_element->Initialize(r_process_info);

        array_1d<double, 3> delta = ZeroVector(3);
        delta[2] = 0.002;

        for (auto& node : p_shell_5p_element->GetGeometry()) // p_shell_3p_element->GetGeometry()[i].
        {
            if(node.Id()==3 || node.Id()==9)
            {
                node.FastGetSolutionStepValue(DISPLACEMENT) = delta*3/20;
                node.Coordinates() += node.FastGetSolutionStepValue(DISPLACEMENT);
            }
            if(node.Id()==4 || node.Id()==10)
            {
                node.FastGetSolutionStepValue(DISPLACEMENT) = delta/2.5;
                node.Coordinates() += node.FastGetSolutionStepValue(DISPLACEMENT);
            }
            if(node.Id()==5 || node.Id()==11)
            {
                node.FastGetSolutionStepValue(DISPLACEMENT) = delta*7/10;
                node.Coordinates() += node.FastGetSolutionStepValue(DISPLACEMENT);
            }
            if(node.Id()==6 || node.Id()==12)
            {
                node.FastGetSolutionStepValue(DISPLACEMENT) = delta;
                node.Coordinates() += node.FastGetSolutionStepValue(DISPLACEMENT);
            }
        }

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_5p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        const double tolerance = 1.0e-8;

        const std::array<double, 60> expected_LHS_row_0{ 405000.2920386971,123985.6767727672,80.66438668106139,-1.947583862517543,0,-33973.98987364071,-106654.8719703566,-6.766639393688726,-0.3402932967357669,0,-9694.631034092492,-16422.36278641331,-1.930891074070916,-0.02378321776677868,0,-594.5854473656216,-887.2782390413756,-0.118424283791966,-0.0008311087123488509,0,-14.85851927577848,-20.97888398219133,-0.002959388786806722,-1.452161979332796e-05,0,-0.1350847139494703,-0.1848929735518754,-2.690498161968747e-05,-1.014921096548247e-07,0,-276681.9115201863,33221.86196428906,-55.10706360277968,-0.5218535231534619,0,-76407.95747603332,-28578.08681330184,-15.21826315126308,-0.09118131405007307,0,-7301.614434665329,-4400.358846430367,-1.45427116697033,-0.006372693994021872,0,-323.6989034289452,-237.7454876128451,-0.06447149260540978,-0.0002226949082963447,0,-6.853475913909427,-5.621275021134379,-0.001365014883927762,-3.89105629641403e-06,0,-0.05626938075705337,-0.04954192294941391,-1.120723894454616e-05,-2.719472882014121e-08,0 };
        const std::array<double, 60> expected_LHS_row_1{ 123985.6767727671,658198.3110414017,24.69437385671005,0,-1.947583862517543,21663.50600301304,62635.26516362656,4.314746107858398,0,-0.3402932967357669,1514.070026661864,717.5120728971251,0.3015591176080054,0,-0.02378321776677868,52.90944238935797,-102.8288509747996,0.01053803620645195,0,-0.0008311087123488509,0.9244648677593761,-4.031472628554199,0.0001841267608974197,0,-1.452161979332796e-05,0.006461117359220246,-0.04379510479871627,1.286868384749658e-06,0,-1.014921096548247e-07,-123985.6767727671,-594039.1195095758,-24.69437385671005,0,-0.5218535231534619,-21663.50600301305,-117826.2399331534,-4.3147461078584,0,-0.09118131405007307,-1514.070026661865,-9215.634975832734,-0.3015591176080056,0,-0.006372693994021872,-52.90944238935805,-356.3133335293762,-0.01053803620645196,0,-0.0002226949082963447,-0.9244648677593776,-6.824525181613907,-0.00018412676089742,0,-3.89105629641403e-06,-0.00646111735922026,-0.05188194445225964,-1.28686838474966e-06,0,-2.719472882014121e-08 };
        const std::array<double, 60> expected_LHS_row_2{ 80.66438668106139,24.69437385671005,354399.5510951506,-9778.442031659108,-23959.85405612183,-6.766639393688727,-21.24257697042332,9553.75650670598,-1708.546851272185,-4186.406492157968,-1.930891074070917,-3.270861415721915,-2992.373460267294,-119.4109382063518,-292.5894168892806,-0.118424283791966,-0.1767202560822225,-232.4714611550978,-4.172836159776343,-10.22458844309223,-0.002959388786806723,-0.004178389130401298,-6.296664670749723,-0.07291024540083026,-0.1786500173896845,-2.690498161968748e-05,-3.682535217948397e-05,-0.05962661260587909,-0.0005095722602915917,-0.001248591232438192,-55.10706360277968,6.616837532497675,-290240.3544730424,-2620.125645617602,-6420.0235451054,-15.21826315126308,-5.691931344380892,-64744.7356549918,-457.8037490291192,-1121.744238762136,-1.45427116697033,-0.8764246748968054,-5505.750116895074,-31.99606445983464,-78.3990979693758,-0.06447149260540978,-0.04735204990345269,-226.6707597766457,-1.118108079159458,-2.739670216267152,-0.001365014883927762,-0.001119595993154008,-4.559334000714971,-0.01953624137510743,-0.04786912788737195,-1.120723894454616e-05,-9.867323377484408e-06,-0.0360504442359528,-0.0001365394756304337,-0.0003345590124083934 };
        const std::array<double, 60> expected_RHS{ 0.0005839978454105505,1.593836517207118e-12,2.932140369516759,-0.5666271650679264,-3.199913269787143e-17,-0.0005023662132058553,-1.371048931691938e-12,-2.522283733491812,-0.09900442786158557,-5.591076496376472e-18,-7.73526802724566e-05,-2.111095586934594e-13,-0.38837286807849,-0.006919454159967685,-3.907622957595076e-19,-4.17926158646332e-06,-1.14059663593173e-14,-0.02098326525615573,-0.0002418015381030824,-1.365525689786227e-20,-9.881482505137843e-08,-2.696836622697073e-16,-0.0004961301517459051,-4.22489855972176e-06,-2.385926725404626e-22,-8.708836395387572e-10,-2.376800132825633e-18,-4.372538553667234e-06,-2.952796409811002e-08,-1.667532550014852e-24,0.0001564817461075064,4.270672076528801e-13,0.7856646441067073,-0.1518272912894872,-8.574141764890973e-18,-0.0001346086219983917,-3.673714540304082e-13,-0.6758438894712882,-0.02652815649261718,-1.498124432024705e-18,-2.072658827185682e-05,-5.656663576640311e-14,-0.1040641963637903,-0.001854062154227513,-1.047044415812915e-19,-1.119829769164367e-06,-3.056219474875622e-15,-0.005622448979954957,-6.479052686332424e-05,-3.65891505822171e-21,-2.647735260966098e-08,-7.226151951703566e-17,-0.0001329376735010455,-1.132058157180858e-06,-6.393071392720023e-23,-2.333525681739323e-10,-6.368616761608131e-19,-1.171618174329084e-06,-7.911994134223759e-09,-4.468140001290899e-25 };

        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_EXPECT_NEAR(left_hand_side_matrix(0,i), expected_LHS_row_0[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_EXPECT_NEAR(left_hand_side_matrix(1,i), expected_LHS_row_1[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_EXPECT_NEAR(left_hand_side_matrix(2,i), expected_LHS_row_2[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          KRATOS_EXPECT_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }
}
}
