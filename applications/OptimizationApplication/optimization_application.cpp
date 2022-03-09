// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_27.h"
#include "optimization_application.h"
#include "optimization_application_variables.h"

// ==============================================================================

namespace Kratos
{
    KratosOptimizationApplication::KratosOptimizationApplication() :
        KratosApplication("OptimizationApplication"),  
        /* ELEMENTS */
        mHelmholtzSurfShape3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
        mHelmholtzBulkShape3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
        /* CONDITIONS */            
        mHelmholtzSurfShapeCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3))))              
    {}

 	void KratosOptimizationApplication::Register()
 	{
        KRATOS_INFO("") << "Initializing KratosOptimizationApplication..." << std::endl;

        // Register variables 

        //strain energy
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_STRAIN_ENERGY_D_X);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_STRAIN_ENERGY_D_CX);
        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_D_T);
        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_D_CT);       
        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_D_P);
        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_D_CP);    
        
        //mass
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_MASS_D_X);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_MASS_D_CX);
        KRATOS_REGISTER_VARIABLE(D_MASS_D_T);
        KRATOS_REGISTER_VARIABLE(D_MASS_D_CT);       
        KRATOS_REGISTER_VARIABLE(D_MASS_D_P);
        KRATOS_REGISTER_VARIABLE(D_MASS_D_CP);      

        //eigenfrequency
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_EIGEN_FREQ_D_X);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_EIGEN_FREQ_D_CX);
        KRATOS_REGISTER_VARIABLE(D_EIGEN_FREQ_D_T);
        KRATOS_REGISTER_VARIABLE(D_EIGEN_FREQ_D_CT);       
        KRATOS_REGISTER_VARIABLE(D_EIGEN_FREQ_D_P);
        KRATOS_REGISTER_VARIABLE(D_EIGEN_FREQ_D_CP);            

        //local_stress
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_LOCAL_STRESS_D_X);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_LOCAL_STRESS_D_CX);
        KRATOS_REGISTER_VARIABLE(D_LOCAL_STRESS_D_T);
        KRATOS_REGISTER_VARIABLE(D_LOCAL_STRESS_D_CT);       
        KRATOS_REGISTER_VARIABLE(D_LOCAL_STRESS_D_P);
        KRATOS_REGISTER_VARIABLE(D_LOCAL_STRESS_D_CP);                  

        //max_stress
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_MAX_STRESS_D_X);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_MAX_STRESS_D_CX);
        KRATOS_REGISTER_VARIABLE(D_MAX_STRESS_D_T);
        KRATOS_REGISTER_VARIABLE(D_MAX_STRESS_D_CT);       
        KRATOS_REGISTER_VARIABLE(D_MAX_STRESS_D_P);
        KRATOS_REGISTER_VARIABLE(D_MAX_STRESS_D_CP); 

        // shape control
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CX);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_CX);  
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_X);   

        // For implicit vertex-morphing with Helmholtz PDE
        KRATOS_REGISTER_VARIABLE( HELMHOLTZ_MASS_MATRIX_SHAPE );
        KRATOS_REGISTER_VARIABLE( HELMHOLTZ_SURF_RADIUS_SHAPE );
        KRATOS_REGISTER_VARIABLE( HELMHOLTZ_BULK_RADIUS_SHAPE );
        KRATOS_REGISTER_VARIABLE( COMPUTE_CONTROL_POINTS_SHAPE );
        KRATOS_REGISTER_VARIABLE( HELMHOLTZ_SURF_POISSON_RATIO_SHAPE );
        KRATOS_REGISTER_VARIABLE( HELMHOLTZ_BULK_POISSON_RATIO_SHAPE );
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( HELMHOLTZ_VARS_SHAPE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( HELMHOLTZ_SOURCE_SHAPE);  

        // Shape optimization elements
        KRATOS_REGISTER_ELEMENT("HelmholtzSurfShape3D3N", mHelmholtzSurfShape3D3N); 
        KRATOS_REGISTER_ELEMENT("HelmholtzBulkShape3D4N", mHelmholtzBulkShape3D4N); 

        // Shape optimization conditions
        KRATOS_REGISTER_CONDITION("HelmholtzSurfShapeCondition3D3N", mHelmholtzSurfShapeCondition3D3N); 

        // KRATOS_REGISTER_VARIABLE(TEST_MAP);
 	}

}  // namespace Kratos.


