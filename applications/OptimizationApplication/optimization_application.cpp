//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl,
//                   Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_27.h"

// Application includes
#include "optimization_application.h"
#include "optimization_application_variables.h"

namespace Kratos
{
    KratosOptimizationApplication::KratosOptimizationApplication() :
        KratosApplication("OptimizationApplication"),
        /* ELEMENTS */
        mHelmholtzSurfShape3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
        mHelmholtzSurfThickness3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
        mHelmholtzBulkShape3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
        mHelmholtzBulkTopology3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
        mAdjointSmallDisplacementElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4))),Element::Pointer() ),
        /* CONDITIONS */
        mHelmholtzSurfShapeCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3))))
    {}

 	void KratosOptimizationApplication::Register()
 	{

        KRATOS_INFO("") << "_______        ____________             \n"       
                        << "__  __ \\_________  /___    |_______________             \n"
                        << "_  / / /__  __ \\  __/_  /| |__  __ \\__  __ \\            \n"
                        << "/ /_/ /__  /_/ / /_ _  ___ |_  /_/ /_  /_/ /            \n"
                        << "\\____/ _  .___/\\__/ /_/  |_|  .___/_  .___/             \n"
                        << "       /_/                 /_/     /_/                  \n"
                        << "Initializing KratosOptimizationApplication... " << std::endl;


        // Register variables

        //Auxilary field
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(AUXILIARY_FIELD);

        //linear function
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_LINEAR_D_X);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_LINEAR_D_CX);

        //strain energy
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_STRAIN_ENERGY_1_D_X);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_STRAIN_ENERGY_1_D_CX);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_STRAIN_ENERGY_2_D_X);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_STRAIN_ENERGY_2_D_CX);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_STRAIN_ENERGY_3_D_X);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_STRAIN_ENERGY_3_D_CX);

        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_1_D_FD);
        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_1_D_CD);

        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_2_D_FD);
        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_2_D_CD);

        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_3_D_FD);
        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_3_D_CD);

        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_1_D_FT);
        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_1_D_CT);

        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_2_D_FT);
        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_2_D_CT);

        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_3_D_FT);
        KRATOS_REGISTER_VARIABLE(D_STRAIN_ENERGY_3_D_CT);

        // partitioning
        KRATOS_REGISTER_VARIABLE(D_INTERFACE_D_FD);
        KRATOS_REGISTER_VARIABLE(D_INTERFACE_D_CD);
        KRATOS_REGISTER_VARIABLE(D_PARTITION_MASS_D_FD);
        KRATOS_REGISTER_VARIABLE(D_PARTITION_MASS_D_CD);


        //symmetry plane
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_PLANE_SYMMETRY_D_X);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_PLANE_SYMMETRY_D_CX);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NEAREST_NEIGHBOUR_POINT);
        KRATOS_REGISTER_VARIABLE(NEAREST_NEIGHBOUR_COND_ID);
        KRATOS_REGISTER_VARIABLE(NEAREST_NEIGHBOUR_DIST);

        //mass
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_MASS_D_X);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_MASS_D_CX);
        KRATOS_REGISTER_VARIABLE(D_MASS_D_FT);
        KRATOS_REGISTER_VARIABLE(D_MASS_D_CT);
        KRATOS_REGISTER_VARIABLE(D_MASS_D_PD);
        KRATOS_REGISTER_VARIABLE(D_MASS_D_FD);
        KRATOS_REGISTER_VARIABLE(D_MASS_D_CD);

        //max overhang angle
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_MAX_OVERHANG_ANGLE_D_X);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_MAX_OVERHANG_ANGLE_D_CX);

        //stress
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_STRESS_D_X);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_STRESS_D_CX);
        KRATOS_REGISTER_VARIABLE(D_STRESS_D_FD);
        KRATOS_REGISTER_VARIABLE(D_STRESS_D_CD);

        // shape control
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SX);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CX);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_CX);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(D_X);

        // thickness control
        KRATOS_REGISTER_VARIABLE(PT);
        KRATOS_REGISTER_VARIABLE(PPT);
        KRATOS_REGISTER_VARIABLE(FT);
        KRATOS_REGISTER_VARIABLE(CT);
        KRATOS_REGISTER_VARIABLE(D_CT);
        KRATOS_REGISTER_VARIABLE(D_PT);
        KRATOS_REGISTER_VARIABLE(D_PT_D_FT);
        KRATOS_REGISTER_VARIABLE(D_PPT_D_FT);

        // density control
        KRATOS_REGISTER_VARIABLE(PD);
        KRATOS_REGISTER_VARIABLE(PE);
        KRATOS_REGISTER_VARIABLE(D_PD_D_FD);
        KRATOS_REGISTER_VARIABLE(FD);
        KRATOS_REGISTER_VARIABLE(CD);
        KRATOS_REGISTER_VARIABLE(D_PE_D_FD);
        KRATOS_REGISTER_VARIABLE(D_CD);
        KRATOS_REGISTER_VARIABLE(D_PD);

        // For implicit vertex-morphing with Helmholtz PDE
        KRATOS_REGISTER_VARIABLE( HELMHOLTZ_MASS_MATRIX );
        KRATOS_REGISTER_VARIABLE( HELMHOLTZ_SURF_RADIUS_SHAPE );
        KRATOS_REGISTER_VARIABLE( HELMHOLTZ_BULK_RADIUS_SHAPE );
        KRATOS_REGISTER_VARIABLE( COMPUTE_CONTROL_POINTS_SHAPE );
        KRATOS_REGISTER_VARIABLE( ELEMENT_STRAIN_ENERGY );
        KRATOS_REGISTER_VARIABLE( HELMHOLTZ_SURF_POISSON_RATIO_SHAPE );
        KRATOS_REGISTER_VARIABLE( HELMHOLTZ_BULK_POISSON_RATIO_SHAPE );
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( HELMHOLTZ_VARS_SHAPE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( HELMHOLTZ_SOURCE_SHAPE);

        // For thickness optimization
        KRATOS_REGISTER_VARIABLE(HELMHOLTZ_VAR_THICKNESS);
        KRATOS_REGISTER_VARIABLE(HELMHOLTZ_SOURCE_THICKNESS);
        KRATOS_REGISTER_VARIABLE(HELMHOLTZ_RADIUS_THICKNESS);

        // For topology optimization
        KRATOS_REGISTER_VARIABLE(HELMHOLTZ_VAR_DENSITY);
        KRATOS_REGISTER_VARIABLE(HELMHOLTZ_SOURCE_DENSITY);
        KRATOS_REGISTER_VARIABLE(HELMHOLTZ_RADIUS_DENSITY);
        KRATOS_REGISTER_VARIABLE(COMPUTE_CONTROL_DENSITIES);

        KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SHAPE);
        KRATOS_REGISTER_VARIABLE(CROSS_AREA);
        KRATOS_REGISTER_VARIABLE(DENSITY_SENSITIVITY);
        KRATOS_REGISTER_VARIABLE(THICKNESS_SENSITIVITY);
        KRATOS_REGISTER_VARIABLE(CROSS_AREA_SENSITIVITY);
        KRATOS_REGISTER_VARIABLE(YOUNG_MODULUS_SENSITIVITY);
        KRATOS_REGISTER_VARIABLE(POISSON_RATIO_SENSITIVITY);

        // do not expose the following variables to python. They are used
        // as temporary data holders. They can be changed
        // at any point of time in an analysis.
        // Hence, not recommended to be used for calculations
        // unless existing values on those variables are not of interest
        KRATOS_CREATE_VARIABLE(double, TEMPORARY_SCALAR_VARIABLE_1);
        KRATOS_CREATE_VARIABLE(double, TEMPORARY_SCALAR_VARIABLE_2);
        KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(TEMPORARY_ARRAY3_VARIABLE_1);
        KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(TEMPORARY_ARRAY3_VARIABLE_2);

        // Shape optimization elements
        KRATOS_REGISTER_ELEMENT("HelmholtzSurfShape3D3N", mHelmholtzSurfShape3D3N);
        KRATOS_REGISTER_ELEMENT("HelmholtzBulkShape3D4N", mHelmholtzBulkShape3D4N);

        // Topology optimization elements
        KRATOS_REGISTER_ELEMENT("HelmholtzBulkTopology3D4N", mHelmholtzBulkTopology3D4N);

        // Adjoint elements
        KRATOS_REGISTER_ELEMENT("AdjointSmallDisplacementElement3D4N", mAdjointSmallDisplacementElement3D4N);

        // Adjoint RHS
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_RHS);

        // Thickness optimization elements
        KRATOS_REGISTER_ELEMENT("HelmholtzSurfThickness3D3N", mHelmholtzSurfThickness3D3N);

        // Shape optimization conditions
        KRATOS_REGISTER_CONDITION("HelmholtzSurfShapeCondition3D3N", mHelmholtzSurfShapeCondition3D3N);
 	}

}  // namespace Kratos.


