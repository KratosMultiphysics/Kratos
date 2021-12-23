//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//


// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"

#include "fem_to_dem_application.h"
#include "fem_to_dem_application_variables.h"


namespace Kratos {

KratosFemToDemApplication::KratosFemToDemApplication(): KratosApplication("FemToDemApplication"),
mSmallStrainModifiedMohrCoulombFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mSmallStrainModifiedMohrCoulombFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mSmallStrainRankineFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mSmallStrainRankineFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mSmallStrainSimoJuFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mSmallStrainSimoJuFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mSmallStrainDruckerPragerFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mSmallStrainDruckerPragerFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mSmallStrainVonMisesFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mSmallStrainVonMisesFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mSmallStrainTrescaFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mSmallStrainTrescaFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mSmallStrainMohrCoulombFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mSmallStrainMohrCoulombFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mTotalLagrangianModifiedMohrCoulombFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mTotalLagrangianModifiedMohrCoulombFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mTotalLagrangianRankineFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mTotalLagrangianRankineFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mTotalLagrangianSimoJuFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mTotalLagrangianSimoJuFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mTotalLagrangianDruckerPragerFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mTotalLagrangianDruckerPragerFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mTotalLagrangianVonMisesFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mTotalLagrangianVonMisesFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mTotalLagrangianTrescaFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mTotalLagrangianTrescaFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mTotalLagrangianMohrCoulombFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mTotalLagrangianMohrCoulombFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mTotalLagrangianMixturesModifiedMohrCoulombFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mTotalLagrangianMixturesModifiedMohrCoulombFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mTotalLagrangianMixturesRankineFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mTotalLagrangianMixturesRankineFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mTotalLagrangianMixturesSimoJuFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mTotalLagrangianMixturesSimoJuFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mTotalLagrangianMixturesDruckerPragerFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mTotalLagrangianMixturesDruckerPragerFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mTotalLagrangianMixturesVonMisesFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mTotalLagrangianMixturesVonMisesFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mTotalLagrangianMixturesTrescaFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mTotalLagrangianMixturesTrescaFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
mTotalLagrangianMixturesMohrCoulombFemDemElement2D(0, Element::GeometryType::Pointer(new Triangle2D3 <Node<3> >(Element::GeometryType::PointsArrayType(3)))),
mTotalLagrangianMixturesMohrCoulombFemDemElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4))))
{}

void KratosFemToDemApplication::Register()
{
    //REGISTER VARIABLES FEM2DEM
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BACKUP_LAST_STRUCTURAL_VELOCITY)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BACKUP_LAST_STRUCTURAL_DISPLACEMENT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SMOOTHED_STRUCTURAL_VELOCITY)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ACCELERATION_BACKUP)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_BACKUP)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(OLD_RELAXED_VELOCITY)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RELAXED_VELOCITY)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FSI_INTERFACE_RESIDUAL)

    KRATOS_REGISTER_VARIABLE(HARDENING_MODULUS)
    KRATOS_REGISTER_VARIABLE(PRESSURE_VOLUME)
    KRATOS_REGISTER_VARIABLE(PRESSURE_INITIAL_VOLUME)
    KRATOS_REGISTER_VARIABLE(DEM_PARTICLE_POINTER)
    KRATOS_REGISTER_VARIABLE(FRAGILE)
    KRATOS_REGISTER_VARIABLE(VOLUME_COUNTED)
    KRATOS_REGISTER_VARIABLE(ERASED_VOLUME)
    KRATOS_REGISTER_VARIABLE(COHESION_MC)
    KRATOS_REGISTER_VARIABLE(RECOMPUTE_NEIGHBOURS)
    KRATOS_REGISTER_VARIABLE(GENERATE_DEM)
    KRATOS_REGISTER_VARIABLE(DISPLACEMENT_INCREMENT)
    KRATOS_REGISTER_VARIABLE(DAMAGE_ELEMENT)
    KRATOS_REGISTER_VARIABLE(TIME_UNIT_CONVERTER)
    KRATOS_REGISTER_VARIABLE(STRESS_VECTOR)
    KRATOS_REGISTER_VARIABLE(YIELD_STRESS_C)
    KRATOS_REGISTER_VARIABLE(YIELD_STRESS_T)
    KRATOS_REGISTER_VARIABLE(FRAC_ENERGY_T)
    KRATOS_REGISTER_VARIABLE(FRAC_ENERGY_C)
    KRATOS_REGISTER_VARIABLE(INTERNAL_PRESSURE_ITERATION)
    KRATOS_REGISTER_VARIABLE(PFEM_PRESSURE_ITERATION)
    KRATOS_REGISTER_VARIABLE(STRESS_VECTOR_INTEGRATED)
    // KRATOS_REGISTER_VARIABLE(THRESHOLD)
    KRATOS_REGISTER_VARIABLE(SMOOTHED_STRESS_VECTOR)
    KRATOS_REGISTER_VARIABLE(YIELD_SURFACE)
    KRATOS_REGISTER_VARIABLE(STRAIN_VECTOR)
    KRATOS_REGISTER_VARIABLE(SMOOTHING)
    KRATOS_REGISTER_VARIABLE(IS_DAMAGED)
    KRATOS_REGISTER_VARIABLE(TANGENT_CONSTITUTIVE_TENSOR)
    KRATOS_REGISTER_VARIABLE(RECONSTRUCT_PRESSURE_LOAD)
    KRATOS_REGISTER_VARIABLE(IS_DYNAMIC)
    KRATOS_REGISTER_VARIABLE(STRESS_THRESHOLD)
    KRATOS_REGISTER_VARIABLE(DEMFEM_CONTACT)
    KRATOS_REGISTER_VARIABLE(INTEGRATION_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(MAPPING_PROCEDURE)
    KRATOS_REGISTER_VARIABLE(INITIAL_THRESHOLD)
    KRATOS_REGISTER_VARIABLE(IS_DEM)
    KRATOS_REGISTER_VARIABLE(DEM_RADIUS)
    KRATOS_REGISTER_VARIABLE(DEM_GENERATED)
    KRATOS_REGISTER_VARIABLE(INACTIVE_NODE)
    KRATOS_REGISTER_VARIABLE(NUMBER_OF_ACTIVE_ELEMENTS)
    KRATOS_REGISTER_VARIABLE(NODAL_FORCE_APPLIED)
    KRATOS_REGISTER_VARIABLE(NODAL_FORCE_X)
    KRATOS_REGISTER_VARIABLE(NODAL_FORCE_Y)
    KRATOS_REGISTER_VARIABLE(NODAL_FORCE_Z)
    KRATOS_REGISTER_VARIABLE(NODAL_STRESS_VECTOR)
    KRATOS_REGISTER_VARIABLE(EQUIVALENT_NODAL_STRESS)
    KRATOS_REGISTER_VARIABLE(PRESSURE_EXPANDED)
    KRATOS_REGISTER_VARIABLE(IS_SKIN)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(EQUIVALENT_NODAL_STRESS_GRADIENT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(AUXILIAR_GRADIENT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FORCE_LOAD)

    KRATOS_REGISTER_VARIABLE(STRAIN_TENSOR);
    KRATOS_REGISTER_VARIABLE(STRESS_TENSOR);
    KRATOS_REGISTER_VARIABLE(STRESS_TENSOR_INTEGRATED);
    // KRATOS_REGISTER_VARIABLE(HENCKY_STRAIN_VECTOR);

    // KRATOS_REGISTER_VARIABLE(GREEN_LAGRANGE_STRAIN_VECTOR);
    // KRATOS_REGISTER_VARIABLE(BIOT_STRAIN_VECTOR);
    // KRATOS_REGISTER_VARIABLE(ALMANSI_STRAIN_VECTOR);
    // KRATOS_REGISTER_VARIABLE(KIRCHHOFF_STRESS_VECTOR);

    // KRATOS_REGISTER_VARIABLE(CONSTITUTIVE_MATRIX_PK2);
    // KRATOS_REGISTER_VARIABLE(CONSTITUTIVE_MATRIX_KIRCHHOFF);
    // KRATOS_REGISTER_VARIABLE(ALMANSI_STRAIN_TENSOR);

    // KRATOS_REGISTER_VARIABLE(RAYLEIGH_ALPHA);
    // KRATOS_REGISTER_VARIABLE(RAYLEIGH_BETA);
    // KRATOS_REGISTER_VARIABLE(VON_MISES_STRESS);

    // Composite
    KRATOS_REGISTER_VARIABLE(PLASTIC_UNIAXIAL_STRESS);
    KRATOS_REGISTER_VARIABLE(MATRIX_STRESS_TENSOR);
    KRATOS_REGISTER_VARIABLE(FIBER_STRESS_TENSOR);
    KRATOS_REGISTER_VARIABLE(MATRIX_STRESS_VECTOR);
    KRATOS_REGISTER_VARIABLE(FIBER_STRESS_VECTOR);
    KRATOS_REGISTER_VARIABLE(YOUNG_MODULUS_FIBER);
    KRATOS_REGISTER_VARIABLE(DENSITY_FIBER);
    KRATOS_REGISTER_VARIABLE(POISSON_RATIO_FIBER);
    KRATOS_REGISTER_VARIABLE(FIBER_VOLUMETRIC_PART);
    KRATOS_REGISTER_VARIABLE(MATRIX_STRESS_TENSOR_INTEGRATED);
    // KRATOS_REGISTER_VARIABLE(PLASTIC_STRAIN_VECTOR);
    KRATOS_REGISTER_VARIABLE(MAX_PLASTIC_STRAIN);

    KRATOS_REGISTER_VARIABLE(YIELD_STRESS_C_FIBER);
    KRATOS_REGISTER_VARIABLE(YIELD_STRESS_T_FIBER);
    KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_FIBER);
    KRATOS_REGISTER_VARIABLE(ACUMULATED_PLASTIC_STRAIN);
    KRATOS_REGISTER_VARIABLE(EQUIVALENT_STRESS_VM);
    KRATOS_REGISTER_VARIABLE(NODAL_DAMAGE);
    KRATOS_REGISTER_VARIABLE(IS_TAKEN);
    KRATOS_REGISTER_VARIABLE(PRESSURE_ID);

    // Hardening variables plasticity
    KRATOS_REGISTER_VARIABLE(HARDENING_LAW);
    // KRATOS_REGISTER_VARIABLE(MAXIMUM_STRESS);
    // KRATOS_REGISTER_VARIABLE(MAXIMUM_STRESS_POSITION);

    // KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD);
    // KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD);

    //Register element
    KRATOS_REGISTER_ELEMENT("SmallStrainModifiedMohrCoulombFemDemElement2D", mSmallStrainModifiedMohrCoulombFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("SmallStrainModifiedMohrCoulombFemDemElement3D", mSmallStrainModifiedMohrCoulombFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("SmallStrainRankineFemDemElement2D", mSmallStrainRankineFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("SmallStrainRankineFemDemElement3D", mSmallStrainRankineFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("SmallStrainSimoJuFemDemElement2D", mSmallStrainSimoJuFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("SmallStrainSimoJuFemDemElement3D", mSmallStrainSimoJuFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("SmallStrainDruckerPragerFemDemElement2D", mSmallStrainDruckerPragerFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("SmallStrainDruckerPragerFemDemElement3D", mSmallStrainDruckerPragerFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("SmallStrainVonMisesFemDemElement2D", mSmallStrainVonMisesFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("SmallStrainVonMisesFemDemElement3D", mSmallStrainVonMisesFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("SmallStrainTrescaFemDemElement2D", mSmallStrainTrescaFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("SmallStrainTrescaFemDemElement3D", mSmallStrainTrescaFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("SmallStrainMohrCoulombFemDemElement3D", mSmallStrainMohrCoulombFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("SmallStrainMohrCoulombFemDemElement2D", mSmallStrainMohrCoulombFemDemElement2D)

    KRATOS_REGISTER_ELEMENT("TotalLagrangianModifiedMohrCoulombFemDemElement2D", mTotalLagrangianModifiedMohrCoulombFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianModifiedMohrCoulombFemDemElement3D", mTotalLagrangianModifiedMohrCoulombFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianRankineFemDemElement2D", mTotalLagrangianRankineFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianRankineFemDemElement3D", mTotalLagrangianRankineFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianSimoJuFemDemElement2D", mTotalLagrangianSimoJuFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianSimoJuFemDemElement3D", mTotalLagrangianSimoJuFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianDruckerPragerFemDemElement2D", mTotalLagrangianDruckerPragerFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianDruckerPragerFemDemElement3D", mTotalLagrangianDruckerPragerFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianVonMisesFemDemElement2D", mTotalLagrangianVonMisesFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianVonMisesFemDemElement3D", mTotalLagrangianVonMisesFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianTrescaFemDemElement2D", mTotalLagrangianTrescaFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianTrescaFemDemElement3D", mTotalLagrangianTrescaFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMohrCoulombFemDemElement2D", mTotalLagrangianMohrCoulombFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMohrCoulombFemDemElement3D", mTotalLagrangianMohrCoulombFemDemElement3D)

    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixturesModifiedMohrCoulombFemDemElement2D", mTotalLagrangianMixturesModifiedMohrCoulombFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixturesModifiedMohrCoulombFemDemElement3D", mTotalLagrangianMixturesModifiedMohrCoulombFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixturesRankineFemDemElement2D", mTotalLagrangianMixturesRankineFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixturesRankineFemDemElement3D", mTotalLagrangianMixturesRankineFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixturesSimoJuFemDemElement2D", mTotalLagrangianMixturesSimoJuFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixturesSimoJuFemDemElement3D", mTotalLagrangianMixturesSimoJuFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixturesDruckerPragerFemDemElement2D", mTotalLagrangianMixturesDruckerPragerFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixturesDruckerPragerFemDemElement3D", mTotalLagrangianMixturesDruckerPragerFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixturesVonMisesFemDemElement2D", mTotalLagrangianMixturesVonMisesFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixturesVonMisesFemDemElement3D", mTotalLagrangianMixturesVonMisesFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixturesTrescaFemDemElement2D", mTotalLagrangianMixturesTrescaFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixturesTrescaFemDemElement3D", mTotalLagrangianMixturesTrescaFemDemElement3D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixturesMohrCoulombFemDemElement2D", mTotalLagrangianMixturesMohrCoulombFemDemElement2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixturesMohrCoulombFemDemElement3D", mTotalLagrangianMixturesMohrCoulombFemDemElement3D)

    //Register Constitutive Laws
    // Serializer::Register("ElasticIsotropic3D", mElasticIsotropic3D);
    // Serializer::Register("LinearPlaneStress", mLinearPlaneStressFEMDEM);
    // Serializer::Register("LinearPlaneStrain", mLinearPlaneStrain);
    // Serializer::Register("HyperElasticIsotropicNeoHookean3D", mHyperElasticIsotropicNeoHookean3D);
    // Serializer::Register("HyperElasticIsotropicNeoHookeanPlaneStrain2D", mHyperElasticIsotropicNeoHookeanPlaneStrain2D);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("ElasticIsotropic3DFEMDEM", mElasticIsotropic3DFEMDEM)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearPlaneStressFEMDEM", mLinearPlaneStressFEMDEM)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearPlaneStrainFEMDEM", mLinearPlaneStrainFEMDEM)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookean3DFEMDEM", mHyperElasticIsotropicNeoHookean3DFEMDEM)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlaneStrain2DFEMDEM", mHyperElasticIsotropicNeoHookeanPlaneStrain2DFEMDEM)

}

}  // namespace Kratos.
