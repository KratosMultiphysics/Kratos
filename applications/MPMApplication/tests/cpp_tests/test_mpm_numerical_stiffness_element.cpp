// // KRATOS  ___|  |                   |                   |
// //       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
// //             | |   |    |   | (    |   |   | |   (   | |
// //       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
// //
// //  License:         BSD License
// //                   license: StructuralMechanicsApplication/license.txt
// //
// //  Main authors:    Andi Makarim Katili
// //

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "mpm_application_variables.h"
#include "containers/model.h"

#include "custom_utilities/material_point_search_utility.h"
#include "utilities/quadrature_points_utility.h"

// Application includes
#include "custom_elements/mpm_soft_stiffness.hpp"

namespace Kratos
{
namespace Testing
{

void AddVariablesAndDofs(ModelPart& rModelPart)
{
    // Add required variables
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    rModelPart.AddNodalSolutionStepVariable(TOTAL_MP_VOLUME);
    for (auto& r_node : rModelPart.Nodes()){
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(DISPLACEMENT_Z);
    }
}

double CalculateGridVolume(Element::Pointer pElement)
{
    double volume = pElement->GetGeometry().DomainSize();
    if (pElement->GetGeometry().WorkingSpaceDimension() == 2 && pElement->GetProperties().Has(THICKNESS)) {
        volume *= pElement->GetProperties()[THICKNESS];
    }
    return volume;
}

void SetTotalMPVolumeOnGrid(Element::Pointer pElement, const double rVolumeRatio)
{
    pElement->GetGeometry().SetValue(TOTAL_MP_VOLUME, rVolumeRatio * CalculateGridVolume(pElement));
}

void InitializeAllElement(ModelPart& rModelPart)
{
    for (auto& r_element : rModelPart.Elements()) {
        r_element.Initialize(rModelPart.GetProcessInfo());
    }
}

void InitializeSolutionStepAllElement(ModelPart& rModelPart)
{
    for (auto& r_element : rModelPart.Elements()) {
        r_element.InitializeSolutionStep(rModelPart.GetProcessInfo());
    }
}

ModelPart::GeometryType::Pointer Generate2D4NGeometry(ModelPart& rModelpart)
{
    
    rModelpart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    AddVariablesAndDofs(rModelpart);

    // Create Geometry
    rModelpart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    rModelpart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    rModelpart.CreateNewNode(3, 1.0 , 1.0 , 0.0);
    rModelpart.CreateNewNode(4, 0.0 , 1.0 , 0.0);

    std::vector<ModelPart::IndexType> geometry_node_ids {1, 2, 3, 4};

    ModelPart::GeometryType::Pointer p_geometry = rModelpart.CreateNewGeometry("Quadrilateral2D4", geometry_node_ids);

    
    return p_geometry;
}

ModelPart::GeometryType::Pointer Generate2D3NGeometry(ModelPart& rModelpart)
{
        
    rModelpart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    AddVariablesAndDofs(rModelpart);

    // Create Geometry
    rModelpart.CreateNewNode(1, 0.0 , 0.0 , 0.0);
    rModelpart.CreateNewNode(2, 1.0 , 0.0 , 0.0);
    rModelpart.CreateNewNode(3, 0.0 , 1.0 , 0.0);

    std::vector<ModelPart::IndexType> geometry_node_ids {1, 2, 3};

    ModelPart::GeometryType::Pointer p_geometry = rModelpart.CreateNewGeometry("Triangle2D3", geometry_node_ids);
    
    return p_geometry;
}


ModelPart::PropertiesType::Pointer CreateCommonProperties(ModelPart& rModelpart)
{
    auto p_elem_prop = rModelpart.CreateNewProperties(0);
    ConstitutiveLaw const& r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("LinearElasticIsotropicPlaneStrain2DLaw");
    auto p_this_law = r_clone_cl.Clone();
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_this_law);
    p_elem_prop->SetValue(DENSITY, 1000);
    p_elem_prop->SetValue(THICKNESS, 0.1);
    p_elem_prop->SetValue(YOUNG_MODULUS, 1.0e6);
    p_elem_prop->SetValue(POISSON_RATIO, 0.0);
    p_elem_prop->SetValue(PENALTY_FACTOR, 1e-4);
    p_elem_prop->SetValue(VOLUME_RATIO_THRESHOLD, 0.9);
    return p_elem_prop;
}

void AssignPredefinedDisplacement(Element::Pointer pElement)
{
    const unsigned int number_of_nodes = pElement->GetGeometry().size();
    const unsigned int dimension = pElement->GetGeometry().WorkingSpaceDimension();

    double displacement = 0.0;
    for(unsigned int i = 0; i < number_of_nodes; i++){
        array_1d<double, 3>& disp = pElement->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        for(unsigned int j = 0; j < dimension; j++){
            disp[j] = displacement;
            displacement += 0.1;
        }
    }
}

void ConductMPMNumericalStiffnessMatrixTest(Element::Pointer pElement, const bool rRefIsActive, const Matrix rRefLHS, const Vector rRefRHS, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    pElement->InitializeSolutionStep(rCurrentProcessInfo);
    
    if (pElement->IsActive() == true)
    {
        KRATOS_EXPECT_EQ(pElement->IsActive(), rRefIsActive);
        
        double tolerance = 1e-8;
        Matrix lhs = ZeroMatrix(8,8);
        Vector rhs = ZeroVector(8);

        pElement->CalculateRightHandSide(rhs, rCurrentProcessInfo);
        KRATOS_EXPECT_VECTOR_NEAR(rhs, rRefRHS, tolerance);

        pElement->CalculateLeftHandSide(lhs, rCurrentProcessInfo);
        KRATOS_EXPECT_MATRIX_NEAR(lhs, rRefLHS, tolerance);
        
        pElement->CalculateLocalSystem(lhs, rhs, rCurrentProcessInfo);
        KRATOS_EXPECT_VECTOR_NEAR(rhs, rRefRHS, tolerance);
        KRATOS_EXPECT_MATRIX_NEAR(lhs, rRefLHS, tolerance);
    }
    else if (pElement->IsActive() == false)
    {
        KRATOS_EXPECT_EQ(pElement->IsActive(), rRefIsActive);
    }
    KRATOS_CATCH("ConductMPMNumericalStiffnessMatrixTest");
}

KRATOS_TEST_CASE_IN_SUITE(TestNumericalStiffnessElement2D4N, KratosMPMFastSuite)
{
    // Create Modelpart
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart",1);
    
    auto p_geometry = Generate2D4NGeometry(r_model_part);
    
    // Create material properties
    ModelPart::PropertiesType::Pointer p_properties = CreateCommonProperties(r_model_part);
    
    // Create test element and assign TOTAL_MP_VOLUME by volume ratio
    const Element& new_element = KratosComponents<Element>::Get("MPMSoftStiffness2D4N");
    Element::Pointer p_element = new_element.Create(1, p_geometry, p_properties);
    r_model_part.AddElement(p_element);
    
    
    InitializeAllElement(r_model_part);
    KRATOS_EXPECT_DOUBLE_EQ(p_element->GetValue(GRID_VOLUME), 0.1);
    
    bool is_active;
    Matrix ref_lhs = ZeroMatrix(8,8);
    Vector ref_rhs = ZeroVector(8);
    
    const ProcessInfo& process_info = r_model_part.GetProcessInfo();

    // Empty grid
    SetTotalMPVolumeOnGrid(p_element, 0.0);
    is_active = false;
    ConductMPMNumericalStiffnessMatrixTest(p_element, is_active, ref_lhs, ref_rhs, process_info);
    
    // Almost Full Grid
    SetTotalMPVolumeOnGrid(p_element, 0.99);
    is_active = false;
    ConductMPMNumericalStiffnessMatrixTest(p_element, is_active, ref_lhs, ref_rhs, process_info);
    
    // Half Filled Grid 
    SetTotalMPVolumeOnGrid(p_element, 0.5);
    is_active = true;

    ref_lhs(0,0)=   10 ; ref_lhs(0,1)=  2.5; ref_lhs(0,2)=   -5; ref_lhs(0,3)= -2.5; ref_lhs(0,4)=   -5; ref_lhs(0,5)= -2.5; ref_lhs(0,6)=    0; ref_lhs(0,7)=  2.5;
    ref_lhs(1,0)=  2.5 ; ref_lhs(1,1)=   10; ref_lhs(1,2)=  2.5; ref_lhs(1,3)=    0; ref_lhs(1,4)= -2.5; ref_lhs(1,5)=   -5; ref_lhs(1,6)= -2.5; ref_lhs(1,7)=   -5;
    ref_lhs(2,0)=   -5 ; ref_lhs(2,1)=  2.5; ref_lhs(2,2)=   10; ref_lhs(2,3)= -2.5; ref_lhs(2,4)=    0; ref_lhs(2,5)= -2.5; ref_lhs(2,6)=   -5; ref_lhs(2,7)=  2.5;
    ref_lhs(3,0)= -2.5 ; ref_lhs(3,1)=    0; ref_lhs(3,2)= -2.5; ref_lhs(3,3)=   10; ref_lhs(3,4)=  2.5; ref_lhs(3,5)=   -5; ref_lhs(3,6)=  2.5; ref_lhs(3,7)=   -5;
    ref_lhs(4,0)=   -5 ; ref_lhs(4,1)= -2.5; ref_lhs(4,2)=    0; ref_lhs(4,3)=  2.5; ref_lhs(4,4)=   10; ref_lhs(4,5)=  2.5; ref_lhs(4,6)=   -5; ref_lhs(4,7)= -2.5;
    ref_lhs(5,0)= -2.5 ; ref_lhs(5,1)=   -5; ref_lhs(5,2)= -2.5; ref_lhs(5,3)=   -5; ref_lhs(5,4)=  2.5; ref_lhs(5,5)=   10; ref_lhs(5,6)=  2.5; ref_lhs(5,7)=    0;
    ref_lhs(6,0)=    0 ; ref_lhs(6,1)= -2.5; ref_lhs(6,2)=   -5; ref_lhs(6,3)=  2.5; ref_lhs(6,4)=   -5; ref_lhs(6,5)=  2.5; ref_lhs(6,6)=   10; ref_lhs(6,7)= -2.5;
    ref_lhs(7,0)=  2.5 ; ref_lhs(7,1)=   -5; ref_lhs(7,2)=  2.5; ref_lhs(7,3)=   -5; ref_lhs(7,4)= -2.5; ref_lhs(7,5)=    0; ref_lhs(7,6)= -2.5; ref_lhs(7,7)=   10;
    
    ref_rhs(0)= 3.0;
    ref_rhs(1)= 7.0;
    ref_rhs(2)= 1.0;
    ref_rhs(3)= 1.0;
    ref_rhs(4)=-1.0;
    ref_rhs(5)=-5.0;
    ref_rhs(6)=-3.0;
    ref_rhs(7)=-3.0;

    AssignPredefinedDisplacement(p_element);
 
    ConductMPMNumericalStiffnessMatrixTest(p_element, is_active, ref_lhs, ref_rhs, process_info);

}

KRATOS_TEST_CASE_IN_SUITE(TestNumericalStiffnessElement2D3N, KratosMPMFastSuite)
{
    // Create Modelpart
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart",1);
    
    auto p_geometry = Generate2D3NGeometry(r_model_part);
    
    // Create material properties
    ModelPart::PropertiesType::Pointer p_properties = CreateCommonProperties(r_model_part);
    
    // Create test element and assign TOTAL_MP_VOLUME by volume ratio
    const Element& new_element = KratosComponents<Element>::Get("MPMSoftStiffness2D3N");
    Element::Pointer p_element = new_element.Create(1, p_geometry, p_properties);
    r_model_part.AddElement(p_element);
    
    
    InitializeAllElement(r_model_part);
    KRATOS_EXPECT_DOUBLE_EQ(p_element->GetValue(GRID_VOLUME), 0.05);
    
    bool is_active;
    Matrix ref_lhs = ZeroMatrix(6,6);
    Vector ref_rhs = ZeroVector(6);
    
    const ProcessInfo& process_info = r_model_part.GetProcessInfo();

    // Empty grid
    SetTotalMPVolumeOnGrid(p_element, 0.0);
    is_active = false;
    ConductMPMNumericalStiffnessMatrixTest(p_element, is_active, ref_lhs, ref_rhs, process_info);
    
    // Almost Full Grid
    SetTotalMPVolumeOnGrid(p_element, 0.99);
    is_active = false;
    ConductMPMNumericalStiffnessMatrixTest(p_element, is_active, ref_lhs, ref_rhs, process_info);
    
    // Half Filled Grid 
    SetTotalMPVolumeOnGrid(p_element, 0.5);
    is_active = true;

    ref_lhs(0,0)=  3.75; ref_lhs(0,1)=  1.25; ref_lhs(0,2)= -2.5; ref_lhs(0,3)= -1.25; ref_lhs(0,4)= -1.25; ref_lhs(0,5)=    0;
    ref_lhs(1,0)=  1.25; ref_lhs(1,1)=  3.75; ref_lhs(1,2)=    0; ref_lhs(1,3)= -1.25; ref_lhs(1,4)= -1.25; ref_lhs(1,5)= -2.5;
    ref_lhs(2,0)=  -2.5; ref_lhs(2,1)=     0; ref_lhs(2,2)=  2.5; ref_lhs(2,3)=     0; ref_lhs(2,4)=     0; ref_lhs(2,5)=    0;
    ref_lhs(3,0)= -1.25; ref_lhs(3,1)= -1.25; ref_lhs(3,2)=    0; ref_lhs(3,3)=  1.25; ref_lhs(3,4)=  1.25; ref_lhs(3,5)=    0;
    ref_lhs(4,0)= -1.25; ref_lhs(4,1)= -1.25; ref_lhs(4,2)=    0; ref_lhs(4,3)=  1.25; ref_lhs(4,4)=  1.25; ref_lhs(4,5)=    0;
    ref_lhs(5,0)=     0; ref_lhs(5,1)=  -2.5; ref_lhs(5,2)=    0; ref_lhs(5,3)=     0; ref_lhs(5,4)=     0; ref_lhs(5,5)=  2.5;
    
    ref_rhs(0)=  1.25;
    ref_rhs(1)=  1.75;
    ref_rhs(2)=  -0.5;
    ref_rhs(3)= -0.75;
    ref_rhs(4)= -0.75;
    ref_rhs(5)=    -1;

    AssignPredefinedDisplacement(p_element);
 
    ConductMPMNumericalStiffnessMatrixTest(p_element, is_active, ref_lhs, ref_rhs, process_info);

}

// TODO:
// KRATOS_TEST_CASE_IN_SUITE(TestNumericalStiffnessElement3D4N, KratosMPMFastSuite)
// KRATOS_TEST_CASE_IN_SUITE(TestNumericalStiffnessElement3D8N, KratosMPMFastSuite)


} // namespace Testing
} // namespace Kratos

