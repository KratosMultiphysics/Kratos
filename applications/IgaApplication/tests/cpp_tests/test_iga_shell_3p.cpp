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
        const std::array<double, 24> expected_LHS_row_2{0,0,198.447972994059,0,0,-171.85483448436,0,0,-25.7137222612275,0,0,-0.879416248471314,0,0,-194.884414941405,0,0,164.99360369436,0,0,28.7455096832649,0,0,1.14530156377983};
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
        const std::array<double, 30> expected_LHS_row_2{0,0,244.082389493276,0,0,-215.697282151207,0,0,-27.3648093072877,0,0,-1.00974038916239,0,0,-0.0105576456187532,0,0,-234.831164256541,0,0,198.105502702631,0,0,34.8171036261213,0,0,1.87558915028264,0,0,0.0329687775055955};
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
        const std::array<double, 36> expected_LHS_row_2{0,0,285.621720485676,0,0,-260.655462589834,0,0,-24.4416550804041,0,0,-0.530647626131597,0,0,0.00583242108537799,0,0,0.00021238960772272,0,0,-267.295948175435,0,0,225.925109943067,0,0,38.9921798333809,0,0,2.31834885048427,0,0,0.0597399071372794,0,0,0.000569641365150261};
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

        const std::array<double, 24> expected_LHS_row_0{637725.129698895,143581.472486388,256.375088968514,64285.8736029236,-122155.580418668,25.9140942550857,-1057.60697536529,-20626.5743651157,-0.416002461296868,-171.90010284355,-799.317702604558,-0.0688129633389082,-545148.471104723,38472.5396007991,-219.139533876989,-143047.814601456,-32731.4891241371,-57.5916474964064,-12241.7365909422,-5526.87394375324,-4.93459424804399,-343.473926488987,-214.176532908792,-0.138592177524092};
        const std::array<double, 24> expected_LHS_row_1{143581.472486388,1165930.89927327,57.7669503941626,32138.8380829665,221748.121371769,12.9303776672213,2397.95310503603,13618.0780009421,0.964765409264601,59.638943999861,265.89380124329,0.0239944601482397,-143581.472486388,-1119642.56623016,-57.7669503941626,-32138.8380829665,-261129.095059747,-12.9303776672213,-2397.95310503603,-20267.7503206756,-0.964765409264601,-59.638943999861,-523.580836637498,-0.0239944601482397};
        const std::array<double, 24> expected_LHS_row_2{256.375088968514,57.7669503941626,198.556668578815,25.9062481999077,-49.1775574947303,-171.84949946493,-0.408741818457946,-8.27010688205126,-25.7147352128392,-0.0682275509998648,-0.319286017381139,-0.879476378663411,-219.139533876989,15.4786077073247,-194.97066475203,-57.5937498405551,-13.1770868164481,164.969164981308,-4.93264876465877,-2.2159684603647,28.7433049591622,-0.138435316760606,-0.0855524305118256,1.14523728917786};
        const std::array<double, 24> expected_RHS{0.00288270349253658,-5.01616466158904E-019,-0.00356240255755486,-0.00245407306041304,4.27030896326048E-019,0.0068602471039795,-0.000412697326582426,7.16876543642762E-020,-0.00303195284468141,-0.000015933105541118,2.77054484483851E-021,-0.000265891701743228,0.000772418072843557,-1.34407727017432E-019,-0.000954542888411407,-0.000657566894704647,1.14422583813703E-019,0.00183819767138927,-0.000110581915376244,1.93676639663235E-020,-0.000812409316221627,-0.0000042692627626624,7.44849861145786E-022,-0.0000712454667562347};
        
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

        const std::array<double, 30> expected_LHS_row_0{491667.164763356,133490.278033997,135.085211471484,4078.07290659722,-113779.527521874,1.16888308148715,-6544.23930622452,-18740.6082379119,-1.7930810783351,-439.346523032982,-954.225815633574,-0.120575467804896,-8.16984419003921,-15.9164585781728,-0.00224423206224352,-379618.737969491,35768.6121966157,-104.273452434095,-99581.7848523707,-30487.1325146808,-27.4284076989117,-9186.16219171666,-5021.53084301653,-2.53496586966148,-361.606946770118,-255.684036695945,-0.0999317597886351,-5.19003615686668,-4.2648022223848,-0.00143601231209576};
        const std::array<double, 30> expected_LHS_row_1{133490.278033997,850779.619039147,36.7310242286304,26281.0006674328,121138.248726319,7.23144851059694,1940.28454614976,5520.78593717874,0.533886360300434,63.6658342007804,68.8450707925089,0.0175182143074813,0.783391240061524,-0.534780425137672,0.000215556990688626,-133490.278033997,-794755.403521703,-36.7310242286304,-26281.0006674328,-168890.106508665,-7.23144851059695,-1940.28454614976,-13385.986981992,-0.533886360300435,-63.6658342007804,-469.321820654574,-0.0175182143074813,-0.783391240061524,-6.1451599961676,-0.000215556990688627};
        const std::array<double, 30> expected_LHS_row_2{135.085211471484,36.7310242286304,244.122623944223,1.15497519158943,-31.3494560009719,-215.699845750257,-1.78054225167203,-5.11876122575071,-27.3657652340317,-0.119240096491757,-0.258529237487944,-1.0097966684669,-0.00221054014070989,-0.00427776441967927,-0.0105586416727839,-104.273452434095,9.8420482792294,-234.858714019245,-27.4321343067783,-8.40006141861546,198.09724563353,-2.53160610118308,-1.37156793668765,34.8162868368328,-0.0995739481236837,-0.0692727004047286,1.87555561625438,-0.00142698458892938,-0.00114622352166365,0.0329682828323349};
        const std::array<double, 30> expected_RHS{0.00122350899926589,4.38914350736193E-019,-0.00485294650587845,-0.00104424916932567,-3.68996506433975E-019,0.00922852914943032,-0.00017050573885069,-6.03670529118336E-020,-0.00390959070272845,-0.00000861159892182329,-3.0489111467248E-021,-0.000454234809910648,-0.0000001424921677044,-5.04270550557E-023,-0.0000117571309127747,0.000327838248285506,1.10335557272313E-019,-0.00130034309716154,-0.000279805721617685,-9.97812144781241E-020,0.0024727769329169,-0.0000456868750299144,-1.61753030771727E-020,-0.00104757167153232,-0.00000230747097664328,-8.16953279559158E-022,-0.000121711850489663,-0.0000000381806612641542,-1.35177502776037E-023,-0.00000315031373338496};
        
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

        const std::array<double, 36> expected_LHS_row_0{405000.292041965,123985.676762889,80.6277394114706,-33973.989872873,-106654.871951637,-6.73213543680337,-9694.63103772891,-16422.3627942561,-1.92870390614012,-594.585447751078,-887.27824000496,-0.118463439890632,-14.8585192892895,-20.9788840175353,-0.00296397518603882,-0.135084714107063,-0.184892973973396,-0.0000269774668300175,-276681.91153295,33221.8619616424,-55.0484026139658,-76407.9574640947,-28578.0868082858,-15.2679262178776,-7301.61443383311,-4400.35884853185,-1.46275563214187,-323.69890343462,-237.745487871037,-0.0649720522051711,-6.85347591522185,-5.62127503060477,-0.00137783101047213,-0.0562693807789405,-0.04954192306236,-0.000011328782683525};
        const std::array<double, 36> expected_LHS_row_1{123985.676762889,658198.311041401,24.7439681392838,21663.5060012871,62635.2651636265,4.32341151233265,1514.07002654124,717.51207289712,0.302164745763555,52.9094423851426,-102.828850974799,0.0105591999884707,0.924464867685722,-4.03147262855419,0.000184496546932986,0.00646111735870547,-0.0437951047987161,0.00000128945283231161,-123985.676762889,-594039.119509575,-24.7439681392838,-21663.5060012871,-117826.239933153,-4.32341151233266,-1514.07002654124,-9215.63497583271,-0.302164745763557,-52.9094423851426,-356.313333529375,-0.0105591999884708,-0.924464867685722,-6.82452518161388,-0.000184496546932989,-0.00646111735870547,-0.0518819444522594,-0.00000128945283231163};
        const std::array<double, 36> expected_LHS_row_2{80.6277394114706,24.7439681392838,285.639600068782,-6.75308303970435,-21.3365662785708,-260.658554245773,-1.90995235620106,-3.23148394288521,-24.4422975516996,-0.116344128384243,-0.171882275361106,-0.530684697506728,-0.0028881276450104,-0.00400093349197295,0.00583152279946594,-0.0000260835522926819,-0.0000347089745884912,0.000212381612678957,-55.0484026139658,6.63012628046235,-267.306208556034,-15.2735391111583,-5.71711570359619,225.921637586169,-1.45773116947887,-0.865873512850234,38.9918219105388,-0.0644041843985242,-0.0460557168762323,2.31833239088707,-0.00135750772310567,-0.00107204689814478,0.0597395518076507,-0.0000110892590051435,-0.0000093002417110986,0.000569638417489115};
        const std::array<double, 36> expected_RHS{0.000585170698120203,6.72831526014287E-019,-0.0058885688201771,-0.000504588969498596,-5.79645602453129E-019,0.0111599251128885,-0.0000764214415484539,-8.79029646111829E-020,-0.00467558539023508,-0.00000406484806731851,-4.67273747790602E-021,-0.000574449008272016,-0.0000000946181724570153,-1.08792339679098E-022,-0.0000210706013290123,-0.000000000820833375518662,-9.43670366120384E-025,-0.00000025129287524706,0.000156796015995663,1.80284664037728E-019,-0.00157783725994144,-0.000135204206886802,-1.55991059780199E-019,0.00299029292159018,-0.0000204770635473305,-2.35113103357035E-020,-0.00125281932945626,-0.0000010891727569932,-1.2533755475277E-021,-0.000153923147859347,-0.0000000253528628991659,-2.91508195597068E-023,-0.00000564585061014701,-0.000000000219941640090738,-2.53016761580773E-025,-0.0000000673337229861447};

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
