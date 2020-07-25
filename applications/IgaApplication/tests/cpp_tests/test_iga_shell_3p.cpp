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

#include "custom_elements/shell_3p_element.h"

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

    typename Shell3pElement::Pointer GetShell3pElement(
        ModelPart& rModelPart, SizeType PolynomialDegree, IntegrationPoint<3> IntegrationPoint)
    {
        // Set the element properties
        auto p_elem_prop = rModelPart.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 200000000);
        p_elem_prop->SetValue(POISSON_RATIO, 0.0);
        p_elem_prop->SetValue(THICKNESS, 0.01);
        const auto& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticPlaneStress2DLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        auto p_quadrature_point = TestCreationUtility::GetQuadraturePointGeometry(
            rModelPart, PolynomialDegree, IntegrationPoint);

        return Kratos::make_intrusive<Shell3pElement>(1, p_quadrature_point, p_elem_prop);
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=3.
    KRATOS_TEST_CASE_IN_SUITE(IgaShell3pElementP3, KratosIgaFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

        IntegrationPoint<3> integration_point(0.0694318442029737, 0.211324865405187, 0.0, 0.086963711284364);
        auto p_shell_3p_element = GetShell3pElement(r_model_part, 3, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);

        p_shell_3p_element->Initialize();

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_3p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        const double tolerance = 1.0e-8;

        const std::array<double, 24> expected_LHS_row_0{637725.123769309,143581.472497936,0,64285.8786450023,-122155.580440902,0,-1057.60612106728,-20626.5743552912,0,-171.900069634590,-799.317701742961,0,-545148.472653437,38472.5396038933,0,-143047.813284562,-32731.4891300946,0,-12241.7363677962,-5526.87394112079,0,-343.473917813981,-214.176532677927,0};
        const std::array<double, 24> expected_LHS_row_1{143581.472497936,1165930.89336990,0,32138.8380855513,221748.126394209,0,2397.95310522888,13618.0788490058,0,59.6389440046576,265.893834107322,0,-143581.472497936,-1119642.56781196,0,-32138.8380855513,-261129.093713988,0,-2397.95310522888,-20267.7500934376,0,-59.6389440046576,-523.580827831607,0};
        const std::array<double, 24> expected_LHS_row_2{0,0,9784.68465798909,0,0,-8327.58890049534,0,0,-1402.84978536515,0,0,-54.2459721286143,0,0,-9781.12109993644,0,0,8320.72766970534,0,0,1405.88157278718,0,0,54.5118574439228};
        const std::array<double, 24> expected_RHS{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(0,i), expected_LHS_row_0[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(1,i), expected_LHS_row_1[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(2,i), expected_LHS_row_2[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          KRATOS_CHECK_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=4.
    KRATOS_TEST_CASE_IN_SUITE(IgaShell3pElementP4, KratosIgaFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

        IntegrationPoint<3> integration_point(0.0469100770306680, 0.211324865405187, 0, 0.0592317212640473);
        auto p_shell_3p_element = GetShell3pElement(r_model_part, 4, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);

        p_shell_3p_element->Initialize();

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_3p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        const double tolerance = 1.0e-8;

        const std::array<double, 30> expected_LHS_row_0{491667.161414150,133490.278045005,0,4078.07575802542,-113779.527542805,0,-6544.23883310685,-18740.6082290450,0,-439.346498780415,-954.225814603370,0,-8.16984378284357,-15.9164585515075,0,-379618.738844226,35768.6121995651,0,-99581.7841076655,-30487.1325202893,0,-9186.16206812944,-5021.53084064066,0,-361.606940433798,-255.684036419903,0,-5.19003605046336,-4.26480221523988,0};
        const std::array<double, 30> expected_LHS_row_1{133490.278045005,850779.615701094,0,26281.0006695998,121138.251571485,0,1940.28454630975,5520.78640580566,0,63.6658342060302,68.8450946538420,0,0.783391240126120,-0.534780027131337,0,-133490.278045005,-794755.404416132,0,-26281.0006695998,-168890.105746305,0,-1940.28454630975,-13385.9868564238,0,-63.6658342060302,-469.321814260949,0,-0.783391240126120,-6.14515988952213,0};
        const std::array<double, 30> expected_LHS_row_2{0,0,11846.6051205351,0,0,-10105.0283553663,0,0,-1656.23498738675,0,0,-83.9478154806674,0,0,-1.39396230134425,0,0,-11837.3538952984,0,0,10087.4365759178,0,0,1663.68728170558,0,0,84.8136642417877,0,0,1.41637343323110};
        const std::array<double, 30> expected_RHS{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(0,i), expected_LHS_row_0[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(1,i), expected_LHS_row_1[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(2,i), expected_LHS_row_2[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          KRATOS_CHECK_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=5.
    KRATOS_TEST_CASE_IN_SUITE(IgaShell3pElementP5, KratosIgaFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

        IntegrationPoint<3> integration_point(0.0337652428984240, 0.211324865405187, 0, 0.0428311230947926);
        auto p_shell_3p_element = GetShell3pElement(r_model_part, 5, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);

        p_shell_3p_element->Initialize();

        Matrix left_hand_side_matrix;
        Vector right_hand_side_vector;
        p_shell_3p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        const double tolerance = 1.0e-8;

        const std::array<double, 36> expected_LHS_row_0{405000.290031407,123985.676772767,0,-33973.9881469315,-106654.871970357,0,-9694.63076821954,-16422.3627864133,0,-594.585433000862,-887.278239041373,0,-14.8585189361368,-20.9788839821913,0,-0.135084710956110,-0.184892973551875,0,-276681.912058038,33221.8619642891,0,-76407.9570133629,-28578.0868133018,0,-7301.61436342488,-4400.35884643036,0,-323.698899579919,-237.745487612845,0,-6.85347582290271,-5.62127502113436,0,-0.0562693799549848,-0.0495419229494137,0};
        const std::array<double, 36> expected_LHS_row_1{123985.676772767,658198.309034111,0,21663.5060030130,62635.2668903357,0,1514.07002666186,717.512338770064,0,52.9094423893578,-102.828836610041,0,0.924464867759373,-4.03147228891258,0,0.00646111735922022,-0.0437951018053564,0,-123985.676772767,-594039.120047427,0,-21663.5060030130,-117826.239470483,0,-1514.07002666186,-9215.63490459227,0,-52.9094423893578,-356.313329680350,0,-0.924464867759373,-6.82452509060719,0,-0.00646111735922022,-0.0518819436501911,0};
        const std::array<double, 36> expected_LHS_row_2{0,0,13572.8850580605,0,0,-11690.6158176292,0,0,-1784.38897035118,0,0,-95.6182409894311,0,0,-2.24242692027938,0,0,-0.0196021703704716,0,0,-13554.5592857502,0,0,11655.8854649824,0,0,1798.93949510416,0,0,97.4059422137838,0,0,2.30799924850204,0,0,0.0203842013433446};
        const std::array<double, 36> expected_RHS{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0};
        
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(0,i), expected_LHS_row_0[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(1,i), expected_LHS_row_1[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(2,i), expected_LHS_row_2[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          KRATOS_CHECK_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=3 (with initial displacement)
    KRATOS_TEST_CASE_IN_SUITE(IgaShell3pElementP3Disp, KratosIgaFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

        IntegrationPoint<3> integration_point(0.0694318442029737, 0.211324865405187, 0.0, 0.086963711284364);
        auto p_shell_3p_element = GetShell3pElement(r_model_part, 3, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);

        p_shell_3p_element->Initialize();

        array_1d<double, 3> delta = ZeroVector(3);
        delta[2] = 0.002;

        for (auto& node : p_shell_3p_element->GetGeometry()) // p_shell_3p_element->GetGeometry()[i].
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
        p_shell_3p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        const double tolerance = 1.0e-8;

        const std::array<double, 24> expected_LHS_row_0{637725.131249061,143581.472486388,252.520188482307,64285.8722840803,-122155.580418668,29.1937486436839,-1057.60719805849,-20626.5743651157,0.137783414571448,-171.900111473322,-799.317702604558,-0.0473527418821970,-545148.472654889,38472.5396007991,-215.284633390820,-143047.813282613,-32731.4891241371,-60.8713018849727,-12241.7363682489,-5526.87394375324,-5.48838012390693,-343.473917859216,-214.176532908792,-0.160052398980594};
        const std::array<double, 24> expected_LHS_row_1{143581.472486388,1165930.89927327,57.7669503941045, 32138.8380829665,221748.121371769,12.9303776672083,2397.95310503603,13618.0780009421,0.964765409263631,59.6389439998610,265.893801243290,0.0239944601482156,-143581.472486388,-1119642.56623016,-57.7669503941045,-32138.8380829665,-261129.095059747,-12.9303776672083,-2397.95310503603,-20267.7503206756,-0.964765409263630,-59.6389439998610,-523.580836637498,-0.0239944601482156};
        const std::array<double, 24> expected_LHS_row_2{252.520188482307,57.7669503941045,9784.79180340765,29.1859025885060,-49.1775574946808,-8327.58224663268,0.145044057410362,-8.27010688204294,-1402.85057562355,-0.0467673295431543,-0.319286017380818,-54.2460236290349,-215.284633390820,15.4786077073091,-9781.20579958086,-60.8734042291214,-13.1770868164348,8320.70191214905,-5.48643464052171,-2.21596846036247,1405.87914536987,-0.159895538217108,-0.0855524305117396,54.5117845395494};
        const std::array<double, 24> expected_RHS{0.00288270349253658,-4.86436346513047e-19,-0.00356240255755126,-0.00245407306041304,4.14107915251829e-19,0.00686024710397255,-0.000412697326582429,6.95144184378740e-20,-0.00303195284467833,-1.59331055411182e-05,2.68664219960341e-21,-0.000265891701742958,0.000772418072843556,-1.30340226217317e-19,-0.000954542888410425,-0.000657566894704646,1.10959881471063e-19,0.00183819767138739,-0.000110581915376245,1.87853471548814e-20,-0.000812409316220805,-4.26926276266246e-06,7.22368215112201e-22,-7.12454667561625e-05};
        
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(0,i), expected_LHS_row_0[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(1,i), expected_LHS_row_1[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(2,i), expected_LHS_row_2[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          KRATOS_CHECK_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=4 (with initial displacement)
    KRATOS_TEST_CASE_IN_SUITE(IgaShell3pElementP4Disp, KratosIgaFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

        IntegrationPoint<3> integration_point(0.0469100770306680, 0.211324865405187, 0, 0.0592317212640473);
        auto p_shell_3p_element = GetShell3pElement(r_model_part, 4, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);

        p_shell_3p_element->Initialize();

        array_1d<double, 3> delta = ZeroVector(3);
        delta[2] = 0.002;

        KRATOS_WATCH(delta)

        for (auto& node : p_shell_3p_element->GetGeometry()) // p_shell_3p_element->GetGeometry()[i].
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
        p_shell_3p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        const double tolerance = 1.0e-8;

        const std::array<double, 30> expected_LHS_row_0{491667.165639896,133490.278033997,131.896156282325,4078.07215948445,-113779.527521874,3.88705251019307,-6544.23942928134,-18740.6082379119,-1.34537181584760,-439.346529298734,-954.225815633574,-0.0977792100485633,-8.16984429455177,-15.9164585781728,-0.00186399115264519,-379618.738846031,35768.6121966157,-101.084397244775,-99581.7841052579,-30487.1325146808,-30.1465771277546,-9186.16206865983,-5021.53084301653,-2.98267513217156,-361.606940504366,-255.684036695945,-0.122728017546118,-5.19003605235412,-4.26480222238480,-0.00181625322171328};
        const std::array<double, 30> expected_LHS_row_1{133490.278033997,850779.619039147,36.7310242288219,26281.0006674328,121138.248726319,7.23144851063466,1940.28454614976,5520.78593717874,0.533886360303219,63.6658342007804,68.8450707925089,0.0175182143075726,0.783391240061524,-0.534780425137672,0.000215556990689751,-133490.278033997,-794755.403521703,-36.7310242288219,-26281.0006674328,-168890.106508665,-7.23144851063466,-1940.28454614976,-13385.9869819920,-0.533886360303219,-63.6658342007804,-469.321820654574,-0.0175182143075726,-0.783391240061524,-6.14515999616760,-0.000215556990689751};
        const std::array<double, 30> expected_LHS_row_2{131.896156282325,36.7310242288219,11846.6444784462,3.87314462029527,-31.3494560011355,-10105.0301718526,-1.33283298918446,-5.11876122577742,-1656.23582025666,-0.0964438387354175,-0.258529237489293,-83.9478654942200,-0.00183029923111139,-0.00427776441970159,-1.39396319288572,-101.084397244775,9.84204827928082,-11837.3805685212,-30.1503037356212,-8.40006141865926,10087.4275717359,-2.97931536369314,-1.37156793669480,1663.68634185947,-0.122370205881164,-0.0692727004050899,84.8136244420074,-0.00180722549854685,-0.00114622352166963,1.41637283404528};
        const std::array<double, 30> expected_RHS{0.00122350899926591,3.18192233730121e-20,-0.00485294650590367,-0.00104424916932570,-2.15460595295816e-20,0.00922852914947831,-0.000170505738850679,-3.63509955927068e-21,-0.00390959070274879,-8.61159892182199e-06,-1.83595107447661e-22,-0.000454234809913011,-1.42492167704367e-07,-3.01598804061074e-24,-1.17571309128358e-05,0.000327838248285511,1.25474665262836e-21,-0.00130034309716833,-0.000279805721617692,-6.68214782026876e-21,0.00247277693292979,-4.56868750299114e-05,-9.74021991313307e-22,-0.00104757167153776,-2.30747097664293e-06,-4.91941607749060e-23,-0.000121711850490296,-3.81806612641452e-08,-8.13993158612733e-25,-3.15031373340135e-06};
        
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(0,i), expected_LHS_row_0[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(1,i), expected_LHS_row_1[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(2,i), expected_LHS_row_2[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          KRATOS_CHECK_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }

    // Tests the stiffness matrix of the Shell3pElement with a polynomial degree of p=5 (with initial displacement)
    KRATOS_TEST_CASE_IN_SUITE(IgaShell3pElementP5Disp, KratosIgaFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart");
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

        IntegrationPoint<3> integration_point(0.0337652428984240, 0.211324865405187, 0, 0.0428311230947926);
        auto p_shell_3p_element = GetShell3pElement(r_model_part, 5, integration_point);

        TestCreationUtility::AddDisplacementDofs(r_model_part);

        p_shell_3p_element->Initialize();

        array_1d<double, 3> delta = ZeroVector(3);
        delta[2] = 0.002;

        for (auto& node : p_shell_3p_element->GetGeometry()) // p_shell_3p_element->GetGeometry()[i].
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
        p_shell_3p_element->CalculateLocalSystem(left_hand_side_matrix, right_hand_side_vector, r_model_part.GetProcessInfo());

        //Check RHS and LHS results
        const double tolerance = 1.0e-8;

        const std::array<double, 36> expected_LHS_row_0{405000.292569059,123985.676762889,77.9812995652517,-33973.9903262899,-106654.871951637,-4.45561680719571,-9694.63110754454,-16422.3627942561,-1.57817313211404,-594.585451523123,-887.278240004960,-0.0995247322300679,-14.8585193784760,-20.9788840175353,-0.00251618673220180,-0.135084714893090,-0.184892973973396,-2.30309776448731e-05,-276681.912060044,33221.8619616424,-52.4019627677404,-76407.9570106778,-28578.0868082858,-17.5444448474909,-7301.61436401748,-4400.35884853185,-1.81328640616881,-323.698899662575,-237.745487871037,-0.0839107598657819,-6.85347582603529,-5.62127503060477,-0.00182561946431025,-0.0562693799929136,-0.0495419230623600,-1.52752718686793e-05};
        const std::array<double, 36> expected_LHS_row_1{123985.676762889,658198.311041401,24.7439681392901,21663.5060012871,62635.2651636265,4.32341151233375,1514.07002654124,717.512072897120,0.302164745763633,52.9094423851426,-102.828850974799,0.0105591999884734,0.924464867685722,-4.03147262855419,0.000184496546933035,0.00646111735870547,-0.0437951047987161,1.28945283231195e-06,-123985.676762889,-594039.119509575,-24.7439681392901,-21663.5060012871,-117826.239933153,-4.32341151233376,-1514.07002654124,-9215.63497583271,-0.302164745763634,-52.9094423851426,-356.313333529375,-0.0105591999884735,-0.924464867685722,-6.82452518161388,-0.000184496546933036,-0.00646111735870547,-0.0518819444522594,-1.28945283231196e-06};
        const std::array<double, 36> expected_LHS_row_2{77.9812995652517,24.7439681392901,13572.9024105490,-4.47656441009672,-21.3365662785762,-11690.6184558682,-1.55942158217499,-3.23148394288604,-1784.38954300685,-0.0974054207236786,-0.171882275361150,-95.6182742887614,-0.00244033919117337,-0.00400093349197397,-2.24242772937873,-2.21370631075372e-05,-3.47089745885002e-05,-0.0196021775794885,-52.4019627677404,6.63012628046413,-13554.5690190363,-17.5500577407716,-5.71711570359764,11655.8815392086,-1.80826194350581,-0.865873512850455,1798.93906736569,-0.0833428920591349,-0.0460557168762441,97.4059219821418,-0.00180529617694380,-0.00107204689814506,2.30799880398585,-1.50357481902977e-05,-9.30024171110098e-06,0.0203841976096565};
        const std::array<double, 36> expected_RHS{0.000585170698120203,2.76630098664650e-19,-0.00588856882017852,-0.000504588969498597,-2.38003632711840e-19,0.0111599251128912,-7.64214415484534e-05,-3.61603118801395e-20,-0.00467558539023629,-4.06484806731845e-06,-1.92055147242158e-21,-0.000574449008272163,-9.46181724570131e-08,-4.47292293772705e-23,-2.10706013290177e-05,-8.20833375518636e-10,-3.87908848972582e-25,-2.51292875247124e-07,0.000156796015995664,7.41228115393349e-20,-0.00157783725994192,-0.000135204206886802,-6.44483698874426e-20,0.00299029292159101,-2.04770635473303e-05,-9.64690832217642e-21,-0.00125281932945657,-1.08917275699319e-06,-5.15929529937906e-22,-0.000153923147859386,-2.53528628991654e-08,-1.19851608897061e-23,-5.64585061014844e-06,-2.19941640090731e-10,-1.04100911876724e-25,-6.73337229861618e-08};

        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(0,i), expected_LHS_row_0[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(1,i), expected_LHS_row_1[i], tolerance);
        }
        for (unsigned int i = 0; i < left_hand_side_matrix.size1(); i++) {
          KRATOS_CHECK_NEAR(left_hand_side_matrix(2,i), expected_LHS_row_2[i], tolerance);
        }
        for (unsigned int i = 0; i < right_hand_side_vector.size(); i++) {
          KRATOS_CHECK_NEAR(right_hand_side_vector(i), expected_RHS[i], tolerance);
        }
    }
}
}
