// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:
//


// System includes
#include <string>
#include <vector>
#include <random>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kernel.h"
#include "containers/model.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/properties.h"
#include "includes/process_info.h"
#include "includes/model_part.h"
#include "includes/kratos_components.h"

// Application includes
#include "custom_elements/total_lagrangian.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_3d.h"

namespace Kratos
{
namespace Testing
{

void AssignRandomNodalData(Variable<array_1d<double,3>> const& rVariable, ModelPart& rModelPart, double dmin=-1.0e-2, double dmax=1.0e-2)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(dmin, dmax);
    for (auto& r_node : rModelPart.Nodes())
    {
        array_1d<double, 3>& r_var = r_node.FastGetSolutionStepValue(rVariable);
        auto& r_initial_coords = r_node.GetInitialPosition().Coordinates();
        auto& r_coords = r_node.Coordinates();
        for (std::size_t d = 0; d < static_cast<std::size_t>(rModelPart.GetProcessInfo()[DOMAIN_SIZE]); ++d)
        {
            r_var[d] = dis(gen);
            if (rVariable == DISPLACEMENT)
                r_coords[d] = r_initial_coords[d] + r_var[d];
        }
    }
}

void AssignNodalData3(ModelPart& rModelPart)
{
    KRATOS_CHECK(rModelPart.NumberOfNodes() == 3);
    std::vector<double> dx = {0.00946560399690433464, 0.00662748804282238015, 0.00659328835078053753};
    std::vector<double> dy = {0.00445411379289947958, -0.00237153525872378125, 0.00322058577500765331};
    std::vector<double> dz = {-0.00820947461981770196, 0.00638188396019644284, -0.00772847850915277159};
    for (std::size_t i = 0; i < 3; ++i)
    {
        auto& r_node = *(rModelPart.Nodes().begin() + i);
        r_node.FastGetSolutionStepValue(DISPLACEMENT_X) = dx.at(i);
        r_node.X() = r_node.X0() + r_node.FastGetSolutionStepValue(DISPLACEMENT_X);
        r_node.FastGetSolutionStepValue(DISPLACEMENT_Y) = dy.at(i);
        r_node.Y() = r_node.Y0() + r_node.FastGetSolutionStepValue(DISPLACEMENT_Y);
        if (rModelPart.GetProcessInfo()[DOMAIN_SIZE] == 3)
        {
            r_node.FastGetSolutionStepValue(DISPLACEMENT_Z) = dz.at(i);
            r_node.Z() = r_node.Z0() + r_node.FastGetSolutionStepValue(DISPLACEMENT_Z);
        }
    }
}

void AssignNodalData4(ModelPart& rModelPart)
{
    KRATOS_CHECK(rModelPart.NumberOfNodes() == 4);
    std::vector<double> dx = {0.00946560399690433464, 0.00662748804282238015,
                              0.00659328835078053753, -0.00433392971520204368};
    std::vector<double> dy = {0.00445411379289947958, -0.00237153525872378125,
                              0.00322058577500765331, -0.00818728957827923216};
    std::vector<double> dz = {-0.00820947461981770196, 0.00638188396019644284,
                              -0.00772847850915277159, 0.00205521693631388787};
    for (std::size_t i = 0; i < 4; ++i)
    {
        auto& r_node = *(rModelPart.Nodes().begin() + i);
        r_node.FastGetSolutionStepValue(DISPLACEMENT_X) = dx.at(i);
        r_node.X() = r_node.X0() + r_node.FastGetSolutionStepValue(DISPLACEMENT_X);
        r_node.FastGetSolutionStepValue(DISPLACEMENT_Y) = dy.at(i);
        r_node.Y() = r_node.Y0() + r_node.FastGetSolutionStepValue(DISPLACEMENT_Y);
        if (rModelPart.GetProcessInfo()[DOMAIN_SIZE] == 3)
        {
            r_node.FastGetSolutionStepValue(DISPLACEMENT_Z) = dz.at(i);
            r_node.Z() = r_node.Z0() + r_node.FastGetSolutionStepValue(DISPLACEMENT_Z);
        }
    }
}

void CreateTotalLagrangianTestModelPart(std::string const& rElementName, ModelPart& rModelPart)
{
    KRATOS_TRY;
    ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    const Element& r_elem = KratosComponents<Element>::Get(rElementName);
    r_process_info[DOMAIN_SIZE] = r_elem.WorkingSpaceDimension();
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(DENSITY);
    rModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(THICKNESS);
    Matrix coordinates;
    r_elem.GetGeometry().PointsLocalCoordinates(coordinates);
    if (r_process_info[DOMAIN_SIZE] == 2)
        for (std::size_t i = 0; i < r_elem.GetGeometry().PointsNumber(); ++i)
            rModelPart.CreateNewNode(i + 1, coordinates(i, 0), coordinates(i, 1), 0.0);
    else
        for (std::size_t i = 0; i < r_elem.GetGeometry().PointsNumber(); ++i)
            rModelPart.CreateNewNode(i + 1, coordinates(i, 0), coordinates(i, 1), coordinates(i, 2));
    std::vector<ModelPart::IndexType> node_ids(r_elem.GetGeometry().PointsNumber());
    for (std::size_t i = 0; i < r_elem.GetGeometry().PointsNumber(); ++i)
        node_ids.at(i) = i + 1;
    auto p_prop = rModelPart.pGetProperties(1);
    rModelPart.CreateNewElement(rElementName, 1, node_ids, p_prop);
    rModelPart.SetBufferSize(2);
    for (auto& r_node : rModelPart.Nodes())
    {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
    }
    if (r_process_info[DOMAIN_SIZE] == 2)
        (*p_prop)[CONSTITUTIVE_LAW] = LinearPlaneStrain::Pointer(new LinearPlaneStrain());
    else
        (*p_prop)[CONSTITUTIVE_LAW] = HyperElasticIsotropicKirchhoff3D::Pointer(new HyperElasticIsotropicKirchhoff3D());
    (*p_prop)[DENSITY] = 1000.0;
    (*p_prop)[YOUNG_MODULUS] = 1400000.0;
    (*p_prop)[POISSON_RATIO] = 0.4;
    (*p_prop)[RAYLEIGH_ALPHA] = 0.02;
    (*p_prop)[RAYLEIGH_BETA] = 0.03;
    rModelPart.GetElement(1).Check(r_process_info);
    rModelPart.GetElement(1).Initialize();
    rModelPart.GetElement(1).InitializeSolutionStep(r_process_info);
    rModelPart.GetElement(1).InitializeNonLinearIteration(r_process_info);
    KRATOS_CATCH("");
}

KRATOS_TEST_CASE_IN_SUITE(TotalLagrangian2D3_CalculateLocalSystem, KratosStructuralMechanicsFastSuite)
{
    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("test");
    CreateTotalLagrangianTestModelPart("TotalLagrangianElement2D3N", test_model_part);
    AssignNodalData3(test_model_part);
    auto p_elem = test_model_part.pGetElement(1);
    Vector rhs(6), rhs_ref(6);
    rhs_ref(0) =-7822.53968136757794;
    rhs_ref(1) =-7011.69414932534164;
    rhs_ref(2) = 5422.42769342759675;
    rhs_ref(3) = 2380.19819087945461;
    rhs_ref(4) = 2400.11198793998165;
    rhs_ref(5) = 4631.49595844588657;
    Matrix lhs(6, 6), lhs_ref(6, 6);
    lhs_ref(0, 0) = 1717995.43261930253;
    lhs_ref(0, 1) = 1228007.99107088544;
    lhs_ref(0, 2) = -1480054.41139845038;
    lhs_ref(0, 3) = -238075.960906155466;
    lhs_ref(0, 4) = -237941.021220852126;
    lhs_ref(0, 5) = -989932.030164730037;
    lhs_ref(1, 0) = 1228007.99107088568;
    lhs_ref(1, 1) = 1713784.48650319432;
    lhs_ref(1, 2) = -985010.149534911616;
    lhs_ref(1, 3) = -233066.782798242755;
    lhs_ref(1, 4) = -242997.841535974003;
    lhs_ref(1, 5) = -1480717.70370495156;
    lhs_ref(2, 0) = -1480054.41139845038;
    lhs_ref(2, 1) = -985010.149534911499;
    lhs_ref(2, 2) = 1486054.96427591494;
    lhs_ref(2, 3) = -10926.6087433038465;
    lhs_ref(2, 4) = -6000.55287746450813;
    lhs_ref(2, 5) = 995936.75827821542;
    lhs_ref(3, 0) = -238075.960906155466;
    lhs_ref(3, 1) = -233066.782798242755;
    lhs_ref(3, 2) = -10926.6087433038483;
    lhs_ref(3, 3) = 244008.667851975275;
    lhs_ref(3, 4) = 249002.569649459329;
    lhs_ref(3, 5) = -10941.8850537325252;
    lhs_ref(4, 0) = -237941.021220852097;
    lhs_ref(4, 1) = -242997.841535973974;
    lhs_ref(4, 2) = -6000.55287746450813;
    lhs_ref(4, 3) = 249002.569649459299;
    lhs_ref(4, 4) = 243941.574098316603;
    lhs_ref(4, 5) = -6004.72811348533378;
    lhs_ref(5, 0) = -989932.030164730153;
    lhs_ref(5, 1) = -1480717.70370495156;
    lhs_ref(5, 2) = 995936.758278215537;
    lhs_ref(5, 3) = -10941.885053732527;
    lhs_ref(5, 4) = -6004.72811348533378;
    lhs_ref(5, 5) = 1491659.588758684;
    p_elem->CalculateLocalSystem(lhs, rhs, test_model_part.GetProcessInfo());
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            KRATOS_CHECK_NEAR(lhs(i, j), lhs_ref(i, j), 1e-5);
    for (std::size_t i = 0; i < 6; ++i)
        KRATOS_CHECK_NEAR(rhs(i), rhs_ref(i), 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(TotalLagrangian3D4_CalculateLocalSystem, KratosStructuralMechanicsFastSuite)
{
    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("test");
    CreateTotalLagrangianTestModelPart("TotalLagrangianElement3D4N", test_model_part);
    AssignNodalData4(test_model_part);
    auto p_elem = test_model_part.pGetElement(1);
    Vector rhs(12), rhs_ref(12);
    rhs_ref(0) = 973.552300837932876;
    rhs_ref(1) = 121.952314606513113;
    rhs_ref(2) = 3061.75614150151387;
    rhs_ref(3) = -1731.97991223004692;
    rhs_ref(4) = 818.175593646390325;
    rhs_ref(5) = -114.77383435126508;
    rhs_ref(6) = 795.687321229903091;
    rhs_ref(7) = -1997.59191955930328;
    rhs_ref(8) = 1029.51404683175747;
    rhs_ref(9) = -37.2597098377889395;
    rhs_ref(10) = 1057.46401130639947;
    rhs_ref(11) = -3976.49635398200644;
    Matrix lhs(12, 12), lhs_ref(12, 12);
    lhs_ref(0,0) = 653401.220045543392;
    lhs_ref(0,1) = 397705.399732108228;
    lhs_ref(0,2) = 419040.494193086983;
    lhs_ref(0,3) = -491274.120119777392;
    lhs_ref(0,4) = -78354.1272159001091;
    lhs_ref(0,5) = -90021.7609090902552;
    lhs_ref(0,6) = -81871.3284694776812;
    lhs_ref(0,7) = -325400.735819087189;
    lhs_ref(0,8) = 37.6050580883811563;
    lhs_ref(0,9) = -80255.771456288363;
    lhs_ref(0,10) = 6049.4633028790995;
    lhs_ref(0,11) = -329056.338342085015;
    lhs_ref(1,0) = 397705.399732108228;
    lhs_ref(1,1) = 653231.578756408417;
    lhs_ref(1,2) = 415282.400222067779;
    lhs_ref(1,3) = -324157.260997956153;
    lhs_ref(1,4) = -81378.4649230227806;
    lhs_ref(1,5) = -3147.35741105511261;
    lhs_ref(1,6) = -79877.3783681678469;
    lhs_ref(1,7) = -490850.401433835388;
    lhs_ref(1,8) = -84463.8387448693829;
    lhs_ref(1,9) = 6329.23963401579567;
    lhs_ref(1,10) = -81002.7123995502334;
    lhs_ref(1,11) = -327671.204066143313;
    lhs_ref(2,0) = 419040.494193086983;
    lhs_ref(2,1) = 415282.400222067779;
    lhs_ref(2,2) = 697461.962443730212;
    lhs_ref(2,3) = -342051.845693412353;
    lhs_ref(2,4) = 2174.61754526946333;
    lhs_ref(2,5) = -92322.9731976368494;
    lhs_ref(2,6) = 891.78094151727521;
    lhs_ref(2,7) = -340364.145009994099;
    lhs_ref(2,8) = -85443.3856178080168;
    lhs_ref(2,9) = -77880.4294411918818;
    lhs_ref(2,10) = -77092.8727573431388;
    lhs_ref(2,11) = -519695.603628285171;
    lhs_ref(3,0) = -491274.120119777392;
    lhs_ref(3,1) = -324157.260997956153;
    lhs_ref(3,2) = -342051.845693412353;
    lhs_ref(3,3) = 498918.285717851715;
    lhs_ref(3,4) = -3627.66579176887217;
    lhs_ref(3,5) = 6113.09303605707464;
    lhs_ref(3,6) = -1999.59942692764412;
    lhs_ref(3,7) = 331978.919426071807;
    lhs_ref(3,8) = 156.384413680835564;
    lhs_ref(3,9) = -5644.56617114669825;
    lhs_ref(3,10) = -4193.99263634675663;
    lhs_ref(3,11) = 335782.368243674457;
    lhs_ref(4,0) = -78354.1272159001091;
    lhs_ref(4,1) = -81378.4649230227806;
    lhs_ref(4,2) = 2174.61754526946333;
    lhs_ref(4,3) = -3627.66579176887262;
    lhs_ref(4,4) = 84900.301738915281;
    lhs_ref(4,5) = -1074.02780385466576;
    lhs_ref(4,6) = 83000.8565498197713;
    lhs_ref(4,7) = -3646.71015235031655;
    lhs_ref(4,8) = 1213.3522739833511;
    lhs_ref(4,9) = -1019.06354215078795;
    lhs_ref(4,10) = 124.873336457815611;
    lhs_ref(4,11) = -2313.94201539814867;
    lhs_ref(5,0) = -90021.7609090902843;
    lhs_ref(5,1) = -3147.35741105511261;
    lhs_ref(5,2) = -92322.9731976368494;
    lhs_ref(5,3) = 6113.09303605707464;
    lhs_ref(5,4) = -1074.02780385466576;
    lhs_ref(5,5) = 86895.1865129896323;
    lhs_ref(5,6) = 25.9989198115909517;
    lhs_ref(5,7) = 4857.51298424164906;
    lhs_ref(5,8) = -803.273565265277966;
    lhs_ref(5,9) = 83882.6689532216114;
    lhs_ref(5,10) = -636.127769331870809;
    lhs_ref(5,11) = 6231.06024991249797;
    lhs_ref(6,0) = -81871.3284694776667;
    lhs_ref(6,1) = -79877.3783681678469;
    lhs_ref(6,2) = 891.78094151727521;
    lhs_ref(6,3) = -1999.59942692764412;
    lhs_ref(6,4) = 83000.8565498197568;
    lhs_ref(6,5) = 25.9989198115909517;
    lhs_ref(6,6) = 84862.7659596562735;
    lhs_ref(6,7) = -1987.03891516270164;
    lhs_ref(6,8) = 50.0396256025466428;
    lhs_ref(6,9) = -991.838063250942582;
    lhs_ref(6,10) = -1136.43926648922138;
    lhs_ref(6,11) = -967.819486931412939;
    lhs_ref(7,0) = -325400.735819087247;
    lhs_ref(7,1) = -490850.401433835388;
    lhs_ref(7,2) = -340364.145009994099;
    lhs_ref(7,3) = 331978.919426071807;
    lhs_ref(7,4) = -3646.710152350317;
    lhs_ref(7,5) = 4857.51298424164906;
    lhs_ref(7,6) = -1987.03891516270164;
    lhs_ref(7,7) = 500766.219008822169;
    lhs_ref(7,8) = -832.361853296542108;
    lhs_ref(7,9) = -4591.14469182188623;
    lhs_ref(7,10) = -6269.10742263652173;
    lhs_ref(7,11) = 336338.993879048969;
    lhs_ref(8,0) = 37.6050580883814547;
    lhs_ref(8,1) = -84463.8387448693829;
    lhs_ref(8,2) = -85443.3856178080168;
    lhs_ref(8,3) = 156.384413680835564;
    lhs_ref(8,4) = 1213.3522739833511;
    lhs_ref(8,5) = -803.273565265277966;
    lhs_ref(8,6) = 50.0396256025464936;
    lhs_ref(8,7) = -832.361853296541994;
    lhs_ref(8,8) = 87052.5402584950498;
    lhs_ref(8,9) = -244.029097371763328;
    lhs_ref(8,10) = 84082.8483241825743;
    lhs_ref(8,11) = -805.881075421730429;
    lhs_ref(9,0) = -80255.771456288363;
    lhs_ref(9,1) = 6329.23963401579567;
    lhs_ref(9,2) = -77880.4294411918818;
    lhs_ref(9,3) = -5644.56617114669916;
    lhs_ref(9,4) = -1019.06354215078795;
    lhs_ref(9,5) = 83882.6689532216114;
    lhs_ref(9,6) = -991.838063250942582;
    lhs_ref(9,7) = -4591.14469182188623;
    lhs_ref(9,8) = -244.029097371763356;
    lhs_ref(9,9) = 86892.1756906860101;
    lhs_ref(9,10) = -719.031400043121948;
    lhs_ref(9,11) = -5758.21041465795588;
    lhs_ref(10,0) = 6049.4633028790995;
    lhs_ref(10,1) = -81002.7123995502334;
    lhs_ref(10,2) = -77092.8727573431388;
    lhs_ref(10,3) = -4193.99263634675663;
    lhs_ref(10,4) = 124.873336457815611;
    lhs_ref(10,5) = -636.127769331870809;
    lhs_ref(10,6) = -1136.43926648922161;
    lhs_ref(10,7) = -6269.10742263652173;
    lhs_ref(10,8) = 84082.8483241825743;
    lhs_ref(10,9) = -719.031400043121948;
    lhs_ref(10,10) = 87146.9464857289277;
    lhs_ref(10,11) = -6353.84779750755843;
    lhs_ref(11,0) = -329056.338342085073;
    lhs_ref(11,1) = -327671.204066143313;
    lhs_ref(11,2) = -519695.603628285287;
    lhs_ref(11,3) = 335782.368243674457;
    lhs_ref(11,4) = -2313.94201539814912;
    lhs_ref(11,5) = 6231.06024991249797;
    lhs_ref(11,6) = -967.819486931412939;
    lhs_ref(11,7) = 336338.993879048969;
    lhs_ref(11,8) = -805.881075421730429;
    lhs_ref(11,9) = -5758.21041465795861;
    lhs_ref(11,10) = -6353.84779750755661;
    lhs_ref(11,11) = 514270.424453794491;
    p_elem->CalculateLocalSystem(lhs, rhs, test_model_part.GetProcessInfo());
    for (std::size_t i = 0; i < 12; ++i)
        for (std::size_t j = 0; j < 12; ++j)
        KRATOS_CHECK_NEAR(lhs(i, j), lhs_ref(i, j), 1e-5);
    for (std::size_t i = 0; i < 12; ++i)
        KRATOS_CHECK_NEAR(rhs(i), rhs_ref(i), 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(TotalLagrangian2D3_MassMatrix, KratosStructuralMechanicsFastSuite)
{
    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("test");
    CreateTotalLagrangianTestModelPart("TotalLagrangianElement2D3N", test_model_part);
    AssignNodalData3(test_model_part);
    auto p_elem = test_model_part.pGetElement(1);
    Matrix lhs(6, 6), lhs_ref(6, 6);
    lhs_ref(0, 0) = 83.3333333333332291;
    lhs_ref(0, 1) = 0;
    lhs_ref(0, 2) = 41.6666666666665719;
    lhs_ref(0, 3) = 0;
    lhs_ref(0, 4) = 41.6666666666665719;
    lhs_ref(0, 5) = 0;
    lhs_ref(1, 0) = 0;
    lhs_ref(1, 1) = 83.3333333333332291;
    lhs_ref(1, 2) = 0;
    lhs_ref(1, 3) = 41.6666666666665719;
    lhs_ref(1, 4) = 0;
    lhs_ref(1, 5) = 41.6666666666665719;
    lhs_ref(2, 0) = 41.6666666666665719;
    lhs_ref(2, 1) = 0;
    lhs_ref(2, 2) = 83.3333333333333428;
    lhs_ref(2, 3) = 0;
    lhs_ref(2, 4) = 41.6666666666666217;
    lhs_ref(2, 5) = 0;
    lhs_ref(3, 0) = 0;
    lhs_ref(3, 1) = 41.6666666666665719;
    lhs_ref(3, 2) = 0;
    lhs_ref(3, 3) = 83.3333333333333428;
    lhs_ref(3, 4) = 0;
    lhs_ref(3, 5) = 41.6666666666666217;
    lhs_ref(4, 0) = 41.6666666666665719;
    lhs_ref(4, 1) = 0;
    lhs_ref(4, 2) = 41.6666666666666217;
    lhs_ref(4, 3) = 0;
    lhs_ref(4, 4) = 83.3333333333333428;
    lhs_ref(4, 5) = 0;
    lhs_ref(5, 0) = 0;
    lhs_ref(5, 1) = 41.6666666666665719;
    lhs_ref(5, 2) = 0;
    lhs_ref(5, 3) = 41.6666666666666217;
    lhs_ref(5, 4) = 0;
    lhs_ref(5, 5) = 83.3333333333333428;
    p_elem->CalculateMassMatrix(lhs, test_model_part.GetProcessInfo());
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            KRATOS_CHECK_NEAR(lhs(i, j), lhs_ref(i, j), 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(TotalLagrangian2D3_DampingMatrix, KratosStructuralMechanicsFastSuite)
{
    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("test");
    CreateTotalLagrangianTestModelPart("TotalLagrangianElement2D3N", test_model_part);
    AssignNodalData3(test_model_part);
    auto p_elem = test_model_part.pGetElement(1);
    Matrix lhs(6, 6), lhs_ref(6, 6);
    lhs_ref(0, 0) = 51541.5296452457405;
    lhs_ref(0, 1) = 36840.2397321265598;
    lhs_ref(0, 2) = -44400.7990086201753;
    lhs_ref(0, 3) = -7142.2788271846639;
    lhs_ref(0, 4) = -7137.39730329223039;
    lhs_ref(0, 5) = -29697.9609049419014;
    lhs_ref(1, 0) = 36840.2397321265671;
    lhs_ref(1, 1) = 51415.2012617624932;
    lhs_ref(1, 2) = -29550.3044860473456;
    lhs_ref(1, 3) = -6991.17015061394977;
    lhs_ref(1, 4) = -7289.93524607921972;
    lhs_ref(1, 5) = -44420.6977778152068;
    lhs_ref(2, 0) = -44400.7990086201753;
    lhs_ref(2, 1) = -29550.3044860473456;
    lhs_ref(2, 2) = 44583.3155949441134;
    lhs_ref(2, 3) = -327.798262299115379;
    lhs_ref(2, 4) = -179.183252990601886;
    lhs_ref(2, 5) = 29878.1027483464604;
    lhs_ref(3, 0) = -7142.2788271846639;
    lhs_ref(3, 1) = -6991.17015061394977;
    lhs_ref(3, 2) = -327.798262299115436;
    lhs_ref(3, 3) = 7321.9267022259246;
    lhs_ref(3, 4) = 7470.07708948377967;
    lhs_ref(3, 5) = -327.423218278642423;
    lhs_ref(4, 0) = -7137.39730329222948;
    lhs_ref(4, 1) = -7289.93524607921881;
    lhs_ref(4, 2) = -179.183252990601886;
    lhs_ref(4, 3) = 7470.07708948377876;
    lhs_ref(4, 4) = 7319.91388961616485;
    lhs_ref(4, 5) = -180.141843404560007;
    lhs_ref(5, 0) = -29697.960904941905;
    lhs_ref(5, 1) = -44420.6977778152068;
    lhs_ref(5, 2) = 29878.1027483464641;
    lhs_ref(5, 3) = -327.42321827864248;
    lhs_ref(5, 4) = -180.141843404560007;
    lhs_ref(5, 5) = 44751.4543294271789;
    p_elem->CalculateDampingMatrix(lhs, test_model_part.GetProcessInfo());
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            KRATOS_CHECK_NEAR(lhs(i, j), lhs_ref(i, j), 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(TotalLagrangian3D10_StrainEnergy, KratosStructuralMechanicsFastSuite)
{
    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("test");
    CreateTotalLagrangianTestModelPart("TotalLagrangianElement3D10N", test_model_part);
    KRATOS_CHECK(test_model_part.NumberOfNodes() == 10);
    std::vector<double> dx = {0.00946, 0.00662, 0.00659, 0.00618, 0.00530, 0.00851, 0.00445, -0.00237, 0.00322, 0.00202};
    std::vector<double> dy = {0.00445, -0.00237, 0.00322, 0.00872, -0.00506, 0.00505, 0.00946, 0.00662, 0.00659, 0.00354};
    std::vector<double> dz = {0.00603, -0.00535, 0.00328, 0.00542, 0.00732, 0.00515, 0.00113, -0.00258, 0.00577, 0.00836};
    // Apply small deformation.
    for (std::size_t i = 0; i < test_model_part.NumberOfNodes(); ++i)
    {
        auto& r_node = *(test_model_part.Nodes().begin() + i);
        r_node.FastGetSolutionStepValue(DISPLACEMENT_X) = dx.at(i);
        r_node.X() = r_node.X0() + r_node.FastGetSolutionStepValue(DISPLACEMENT_X);
        r_node.FastGetSolutionStepValue(DISPLACEMENT_Y) = dy.at(i);
        r_node.Y() = r_node.Y0() + r_node.FastGetSolutionStepValue(DISPLACEMENT_Y);
        r_node.FastGetSolutionStepValue(DISPLACEMENT_Z) = dz.at(i);
        r_node.Z() = r_node.Z0() + r_node.FastGetSolutionStepValue(DISPLACEMENT_Z);
    }
    // Calculate strain energy.
    auto p_elem = test_model_part.pGetElement(1);
    std::vector<double> weights;
    std::vector<double> strain_energies;
    p_elem->CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, weights,
                                         test_model_part.GetProcessInfo());
    p_elem->CalculateOnIntegrationPoints(STRAIN_ENERGY, strain_energies,
                                         test_model_part.GetProcessInfo());
    double element_strain_energy = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i)
        element_strain_energy += weights[i] * strain_energies[i];
    // Apply large rigid body rotation (alpha=pi/2).
    for (auto& r_node : test_model_part.Nodes())
    {
        const double x = r_node.X();
        r_node.X() = -r_node.Y();
        r_node.FastGetSolutionStepValue(DISPLACEMENT_X) = r_node.X() - r_node.X0();
        r_node.Y() = x;
        r_node.FastGetSolutionStepValue(DISPLACEMENT_Y) = r_node.Y() - r_node.Y0();
    }
    // Calculate strain energy on rotated element.
    p_elem->CalculateOnIntegrationPoints(INTEGRATION_WEIGHT, weights,
                                         test_model_part.GetProcessInfo());
    p_elem->CalculateOnIntegrationPoints(STRAIN_ENERGY, strain_energies,
                                         test_model_part.GetProcessInfo());
    double rotated_element_strain_energy = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i)
        rotated_element_strain_energy += weights[i] * strain_energies[i];
    // Check that strain energy didn't change.
    KRATOS_CHECK_NEAR(rotated_element_strain_energy, element_strain_energy, 1e-7);
}

KRATOS_TEST_CASE_IN_SUITE(TotalLagrangian3D4_SensitivityMatrix, KratosStructuralMechanicsFastSuite)
{
    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("test");
    CreateTotalLagrangianTestModelPart("TotalLagrangianElement3D4N", test_model_part);
    AssignRandomNodalData(DISPLACEMENT, test_model_part, -1.0e-2, 1.0e-2);
    AssignRandomNodalData(ACCELERATION, test_model_part, -10.0, 10.0);
    AssignRandomNodalData(VOLUME_ACCELERATION, test_model_part, -10.0, 10.0);
    auto p_elem = test_model_part.pGetElement(1);
    Matrix sensitivity_matrix, semi_analytic_sensitivity_matrix;
    p_elem->CalculateSensitivityMatrix(SHAPE_SENSITIVITY, sensitivity_matrix,
                                       test_model_part.GetProcessInfo());
    semi_analytic_sensitivity_matrix.resize(sensitivity_matrix.size1(),
                                            sensitivity_matrix.size2(), false);
    Vector R, R_perturb, semi_analytic_sensitivity_vector, acceleration;
    Matrix mass_matrix;
    auto get_res = [&p_elem, &test_model_part, &mass_matrix, &acceleration](Vector& rRes) {
        p_elem->CalculateRightHandSide(rRes, test_model_part.GetProcessInfo());
        p_elem->CalculateMassMatrix(mass_matrix, test_model_part.GetProcessInfo());
        p_elem->GetSecondDerivativesVector(acceleration);
        if (rRes.size() != acceleration.size())
            rRes.resize(acceleration.size(), false);
        noalias(rRes) -= prod(mass_matrix, acceleration);
    };
    get_res(R);
    std::size_t ws_dim = p_elem->GetGeometry().WorkingSpaceDimension();
    for (std::size_t i = 0; i < p_elem->GetGeometry().PointsNumber(); ++i)
        for (std::size_t d = 0; d < ws_dim; ++d)
        {
            const double delta = 0.00000001;
            p_elem->GetGeometry()[i].GetInitialPosition()[d] += delta;
            p_elem->GetGeometry()[i].Coordinates()[d] += delta;
            get_res(R_perturb);
            p_elem->GetGeometry()[i].GetInitialPosition()[d] -= delta;
            p_elem->GetGeometry()[i].Coordinates()[d] -= delta;
            semi_analytic_sensitivity_vector = (1.0 / delta) * (R_perturb - R);
            for (std::size_t k = 0; k < semi_analytic_sensitivity_vector.size(); ++k)
                semi_analytic_sensitivity_matrix(i * ws_dim + d, k) =
                    semi_analytic_sensitivity_vector(k);
        }
    for (std::size_t i = 0; i < sensitivity_matrix.size1(); ++i)
        for (std::size_t j = 0; j < sensitivity_matrix.size2(); ++j)
            KRATOS_CHECK_NEAR(sensitivity_matrix(i, j), semi_analytic_sensitivity_matrix(i,j), 1e-1);
}
}
}
