from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# path for applications
import sys
import os.path

Import_ExternalSolversApplication = False
Import_SolidMechanicsApplication = False
Import_PfemApplication = False
Import_PfemSolidMechanicsApplication = False
Import_PfemFluidDynamicsApplication = False
Import_ContactMechanicsApplication = False
Import_ConstitutiveModelsApplication = False
Import_UmatApplication = False
Import_MachiningApplication = False
Import_StringDynamicsApplication = False
Import_ConvectionDiffusionApplication = False
Import_MeshMovingApplication = False
Import_IncompressibleFluidApplication = False
Import_StructuralApplication = False
Import_StructuralMechanicsApplication = False
Import_FSIApplication = False
Import_ConstitutiveLawsApplication = False
Import_ULFApplication = False
Import_MeshingApplication = False
Import_KratosMPISearchApplication = False
Import_KratosTrilinosApplication = False
Import_KratosMetisApplication = False
Import_PoissonApplication = False
Import_ElectrostaticApplication = False
Import_MagnetostaticApplication = False
Import_DamApplication = False
Import_TestApplication = False
Import_OpenCLApplication = False
Import_PodApplication = False
Import_LevelSetApplication = False
Import_FluidDynamicsApplication = False
Import_KratosDEMApplication = False
Import_KratosSwimmingDEMApplication = False
Import_KratosMixedElementApplication = False
Import_ThermoMechanicalApplication = False
Import_DEM_FEM_Application = False
Import_MultiScaleApplication = False
Import_ContactStructuralMechanicsApplication = False
Import_KratosMappingApplication = False
Import_ConstitutiveModelsApplication = False
Import_ShallowWaterApplication = False
Import_DelaunayMeshingApplication = False
Import_FluidRveLagrangeMultipliersApplication=False
Import_PoromechanicsApplication = False
Import_FluidTransportApplication = False

print("Applications Available:")
print("Import_FluidRveLagrangeMultipliersApplication: False")
print("Import_ExternalSolversApplication: False")
print("Import_SolidMechanicsApplication: False")
print("Import_PfemApplication: False")
print("Import_PfemSolidMechanicsApplication: False")
print("Import_PfemFluidDynamicsApplication: False")
print("Import_ContactMechanicsApplication: False")
print("Import_ConstitutiveModelsApplication: False")
print("Import_UmatApplication: False")
print("Import_MachiningApplication: False")
print("Import_StringDynamicsApplication: False")
print("Import_ConvectionDiffusionApplication: False")
print("Import_MeshingApplication: False")
print("Import_MeshMovingApplication: False")
print("Import_IncompressibleFluidApplication: False")
print("Import_StructuralApplication: False")
print("Import_StructuralMechanicsApplication: False")
print("Import_FSIApplication: False")
print("Import_ConstitutiveLawsApplication: False")
print("Import_ULFApplication: False")
print("Import_KratosMPISearchApplication: False")
print("Import_KratosTrilinosApplication: False")
print("Import_KratosMetisApplication: False")
print("Import_PoissonApplication: False")
print("Import_ElectrostaticApplication: False")
print("Import_MagnetostaticApplication: False")
print("Import_DamApplication: False")
print("Import_TestApplication: False")
print("Import_OpenCLApplication: False")
print("Import_PodApplication: False")
print("Import_LevelSetApplication: False")
print("Import_FluidDynamicsApplication: False")
print("Import_KratosDEMApplication: False")
print("Import_KratosSwimmingDEMApplication: False")
print("Import_KratosMixedElementApplication: False")
print("Import_ThermoMechanicalApplication: False")
print("Import_DEM_FEM_Application: False")
print("Import_MultiScaleApplication: False")
print("Import_ContactStructuralMechanicsApplication: False")
print("Import_KratosMappingApplication: False")
print("Import_ConstitutiveModelsApplication: False")
print("Import_ShallowWaterApplication: False")
print("Import_DelaunayMeshingApplication: False")
print("Import_PoromechanicsApplication: False")
print("Import_FluidTransportApplication: False")

application_directory = os.path.dirname(os.path.realpath(__file__))

def ImportApplications(kernel, applications_path=application_directory):
    # importing the applications
    print("Applications Available:")
    print("Import_FluidRveLagrangeMultipliersApplication: " + str(Import_FluidRveLagrangeMultipliersApplication))
    print("Import_ExternalSolversApplication: " + str(Import_ExternalSolversApplication))
    print("Import_SolidMechanicsApplication: " + str(Import_SolidMechanicsApplication))
    print("Import_PfemApplication: " + str(Import_PfemApplication))
    print("Import_PfemSolidMechanicsApplication: " + str(Import_PfemSolidMechanicsApplication))
    print("Import_PfemFluidDynamicsApplication: " + str(Import_PfemFluidDynamicsApplication))
    print("Import_ContactMechanicsApplication: " + str(Import_ContactMechanicsApplication))
    print("Import_ConstitutiveModelsApplication: " + str(Import_ConstitutiveModelsApplication))
    print("Import_UmatApplication: " + str(Import_UmatApplication))
    print("Import_MachiningApplication: " + str(Import_MachiningApplication))
    print("Import_StringDynamicsApplication: " + str(Import_StringDynamicsApplication))
    print("Import_ConvectionDiffusionApplication: " + str(Import_ConvectionDiffusionApplication))
    print("Import_MeshingApplication: " + str(Import_MeshingApplication))
    print("Import_MeshMovingApplication: " + str(Import_MeshMovingApplication))
    print("Import_IncompressibleFluidApplication: " + str(Import_IncompressibleFluidApplication))
    print("Import_StructuralApplication: " + str(Import_StructuralApplication))
    print("Import_StructuralMechanicsApplication: " + str(Import_StructuralMechanicsApplication))
    print("Import_FSIApplication: " + str(Import_FSIApplication))
    print("Import_ConstitutiveLawsApplication: " + str(Import_ConstitutiveLawsApplication))
    print("Import_ULFApplication: " + str(Import_ULFApplication))
    print("Import_KratosMPISearchApplication:  " + str(Import_KratosMPISearchApplication))
    print("Import_KratosTrilinosApplication: " + str(Import_KratosTrilinosApplication))
    print("Import_KratosMetisApplication: " + str(Import_KratosMetisApplication))
    print("Import_PoissonApplication: " + str(Import_PoissonApplication))
    print("Import_ElectrostaticApplication: " + str(Import_ElectrostaticApplication))
    print("Import_MagnetostaticApplication: " + str(Import_MagnetostaticApplication))
    print("Import_DamApplication: " + str(Import_DamApplication))
    print("Import_TestApplication: " + str(Import_TestApplication))
    print("Import_OpenCLApplication: " + str(Import_OpenCLApplication))
    print("Import_PodApplication: " + str(Import_PodApplication))
    print("Import_LevelSetApplication:" + str(Import_LevelSetApplication))
    print("Import_FluidDynamicsApplication: " + str(Import_FluidDynamicsApplication))
    print("Import_KratosMixedElementApplication: " + str(Import_KratosMixedElementApplication))
    print("Import_KratosDEMApplication:  " + str(Import_KratosDEMApplication))
    print("Import_KratosSwimmingDEMApplication:  " + str(Import_KratosSwimmingDEMApplication))
    print("Import_ThermoMechanicalApplication: " + str(Import_ThermoMechanicalApplication))
    print("Import_DEM_FEM_Application: " + str(Import_DEM_FEM_Application))
    print("Import_MultiScaleApplication: " + str(Import_MultiScaleApplication))
    print("Import_ContactStructuralMechanicsApplication: " + str(Import_ContactStructuralMechanicsApplication))
    print("Import_KratosMappingApplication: " + str(Import_KratosMappingApplication))
    print("Import_ConstitutiveModelsApplication: " + str(Import_ConstitutiveModelsApplication))
    print("Import_ShallowWaterApplication: " + str(Import_ShallowWaterApplication))
    print("Import_DelaunayMeshingApplication: " + str(Import_DelaunayMeshingApplication))
    print("Import_PoromechanicsApplication: " + str(Import_PoromechanicsApplication))
    print("Import_FluidTransportApplication: " + str(Import_FluidTransportApplication))

    if(Import_FluidRveLagrangeMultipliersApplication):
        print("importing KratosFluidRveLagrangeMultipliersApplication ...")
        sys.path.append(applications_path + '/fluid_rve_lagrange_multipliers_application/python_scripts')
        sys.path.append(applications_path + '/fluid_rve_lagrange_multipliers_application/Linux')
        from KratosFluidRveLagrangeMultipliersApplication import *
        fluid_rve_lagrange_multipliers_application = KratosFluidRveLagrangeMultipliersApplication()
        kernel.ImportApplication(fluid_rve_lagrange_multipliers_application)
        print("KratosFluidRveLagrangeMultipliersApplication sucessfully imported")

    if(Import_ExternalSolversApplication):
        print("importing KratosExternalSolversApplication ...")
        sys.path.append(applications_path + '/ExternalSolversApplication/python_scripts')
        sys.path.append(applications_path + '/ExternalSolversApplication/Linux')
        from KratosExternalSolversApplication import *
        external_solvers_application = KratosExternalSolversApplication()
        kernel.ImportApplication(external_solvers_application)
        print("KratosExternalSolversApplication sucessfully imported")


    if(Import_SolidMechanicsApplication):
        print("importing KratosSolidMechanicsApplication ...")
        sys.path.append(applications_path + '/SolidMechanicsApplication/python_scripts')
        sys.path.append(applications_path + '/SolidMechanicsApplication/python_scripts/constitutive_laws')
        sys.path.append(applications_path + '/SolidMechanicsApplication/Linux')
        from KratosSolidMechanicsApplication import *
        solid_mechanics_application = KratosSolidMechanicsApplication()
        kernel.ImportApplication(solid_mechanics_application)
        print("KratosSolidMechanicsApplication Succesfully imported")

    if(Import_PfemApplication):
        print("importing KratosPfemApplication ...")
        sys.path.append(applications_path + '/PfemApplication/python_scripts')
        sys.path.append(applications_path + '/PfemApplication/Linux')
        from KratosPfemApplication import *
        pfem_application = KratosPfemApplication()
        kernel.ImportApplication(pfem_application)
        print("KratosPfemApplication Succesfully imported")

    if(Import_PfemSolidMechanicsApplication):
        print("importing KratosPfemSolidMechanicsApplication ...")
        sys.path.append(applications_path + '/PfemSolidMechanicsApplication/python_scripts')
        sys.path.append(applications_path + '/PfemSolidMechanicsApplication/Linux')
        from KratosPfemSolidMechanicsApplication import *
        pfem_solid_mechanics_application = KratosPfemSolidMechanicsApplication()
        kernel.ImportApplication(pfem_solid_mechanics_application)
        print("KratosPfemSolidMechanicsApplication Succesfully imported")

    if(Import_PfemFluidDynamicsApplication):
        print("importing KratosPfemFluidDynamicsApplication ...")
        sys.path.append(applications_path + '/PfemFluidDynamicsApplication/python_scripts')
        sys.path.append(applications_path + '/PfemFluidDynamicsApplication/Linux')
        from KratosPfemFluidDynamicsApplication import *
        pfem_fluid_dynamics_application = KratosPfemFluidDynamicsApplication()
        kernel.ImportApplication(pfem_fluid_dynamics_application)
        print("KratosPfemFluidDynamicsApplication Succesfully imported")

    if(Import_ContactMechanicsApplication):
        print("importing KratosContactMechanicsApplication ...")
        sys.path.append(applications_path + '/ContactMechanicsApplication/python_scripts')
        sys.path.append(applications_path + '/ContactMechanicsApplication/Linux')
        from KratosContactMechanicsApplication import *
        contact_mechanics_application = KratosContactMechanicsApplication()
        kernel.ImportApplication(contact_mechanics_application)
        print("KratosContactMechanicsApplication Succesfully imported")

    if(Import_ConstitutiveModelsApplication):
        print("importing KratosConstitutiveModelsApplication ...")
        sys.path.append(applications_path + '/ConstitutiveModelsApplication/python_scripts')
        sys.path.append(applications_path + '/ConstitutiveModelsApplication/Linux')
        from KratosConstitutiveModelsApplication import *
        constitutive_models_application = KratosConstitutiveModelsApplication()
        kernel.ImportApplication(constitutive_models_application)
        print("KratosConstitutiveModelsApplication Succesfully imported")

    if(Import_UmatApplication):
        print("importing KratosUmatApplication ...")
        sys.path.append(applications_path + '/UmatApplication/python_scripts')
        sys.path.append(applications_path + '/UmatApplication/Linux')
        from KratosUmatApplication import *
        umat_application = KratosUmatApplication()
        kernel.ImportApplication(umat_application)
        print("KratosUmatApplication Succesfully imported")

    if(Import_MachiningApplication):
        print("importing KratosMachiningApplication ...")
        sys.path.append(applications_path + '/MachiningApplication/python_scripts')
        sys.path.append(applications_path + '/MachiningApplication/Linux')
        from KratosMachiningApplication import *
        machining_application = KratosMachiningApplication()
        kernel.ImportApplication(machining_application)
        print("KratosMachiningApplication Succesfully imported")

    if(Import_StringDynamicsApplication):
        print("importing KratosStringDynamicsApplication ...")
        sys.path.append(applications_path + '/StringDynamicsApplication/python_scripts')
        sys.path.append(applications_path + '/StringDynamicsApplication/Linux')
        from KratosStringDynamicsApplication import *
        pfem_solid_mechanics_application = KratosStringDynamicsApplication()
        kernel.ImportApplication(string_dynamics_application)
        print("KratosStringDynamicsApplication Succesfully imported")

    if(Import_ConvectionDiffusionApplication):
        print("importing KratosConvectionDiffusionApplication ...")
        sys.path.append(applications_path + '/convection_diffusion_application/python_scripts')
        sys.path.append(applications_path + '/convection_diffusion_application/Linux')
        from KratosConvectionDiffusionApplication import *
        convection_diffusion_application = KratosConvectionDiffusionApplication()
        kernel.ImportApplication(convection_diffusion_application)
        print("KratosConvectionDiffusionApplication Succesfully imported")

    if(Import_MeshingApplication):
        print("importing KratosMeshingApplication ...")
        sys.path.append(applications_path + '/MeshingApplication/python_scripts')
        from KratosMeshingApplication import *
        meshing_application = KratosMeshingApplication()
        kernel.ImportApplication(meshing_application)
        print("KratosMeshingApplication sucessfully imported")

    if(Import_MeshMovingApplication):
        print("importing KratosMeshMovingApplication ...")
        sys.path.append(applications_path + '/MeshMovingApplication/python_scripts')
        sys.path.append(applications_path + '/MeshMovingApplication/Linux')
        from KratosMeshMovingApplication import *
        mesh_moving_app = KratosMeshMovingApplication()
        kernel.ImportApplication(mesh_moving_app)
        print("KratosMeshMovingApplication Succesfully imported")

    if(Import_IncompressibleFluidApplication):
        print("importing KratosIncompressibleFluidApplication ...")
        sys.path.append(applications_path + '/incompressible_fluid_application/python_scripts')
        sys.path.append(applications_path + '/incompressible_fluid_application/Linux')
        from KratosIncompressibleFluidApplication import *
        print("KratosIncompressibleFluidApplication lib loaded")
        incompressible_fluid_application = KratosIncompressibleFluidApplication()
        print("KratosIncompressibleFluidApplication application created")
        kernel.ImportApplication(incompressible_fluid_application)
        print("KratosIncompressibleFluidApplication Succesfully imported")

    if(Import_StructuralApplication):
        print("importing KratosStructuralApplication ...")
        sys.path.append(applications_path + '/structural_application/python_scripts')
        sys.path.append(applications_path + '/structural_application/Linux')
        from KratosStructuralApplication import *
        structural_application = KratosStructuralApplication()
        kernel.ImportApplication(structural_application)
        print("KratosStructuralApplication Succesfully imported")

    if(Import_StructuralMechanicsApplication):
        print("importing StructuralMehcanicsApplication ...")
        sys.path.append(applications_path + '/StructuralMechanicsApplication/python_scripts')
        sys.path.append(applications_path + '/StructuralMechanicsApplication/Linux')
        from StructuralMechanicsApplication import *
        structural_mechanics_application = StructuralMechanicsApplication()
        kernel.ImportApplication(structural_mechanics_application)
        print("StructuralMechanicsApplication Succesfully imported")

    if(Import_FSIApplication):
        print("importing FSIapplication ...")
        sys.path.append(applications_path + '/FSIapplication/python_scripts')
        sys.path.append(applications_path + '/FSIapplication/Linux')
        from KratosFSIApplication import *
        fsi_application = KratosFSIApplication()
        kernel.ImportApplication(fsi_application)
        print("FSIapplication Succesfully imported")

    if(Import_ConstitutiveLawsApplication):
        print("importing KratosConstitutiveLawsApplication ...")
        sys.path.append(applications_path + '/constitutive_laws_application/python_scripts')
        sys.path.append(applications_path + '/constitutive_laws_application/Linux')
        from KratosConstitutiveLawsApplication import *
        constitutive_laws_application = KratosConstitutiveLawsApplication()
        kernel.ImportApplication(constitutive_laws_application)
        print("KratosConstitutiveLawsApplication successfully imported")

    if(Import_ULFApplication):
        print("importing KratosULFApplication ...")
        sys.path.append(applications_path + '/ULFapplication/python_scripts')
        sys.path.append(applications_path + '/ULFapplication/Linux')
        from KratosULFApplication import *
        ulf_application = KratosULFApplication()
        kernel.ImportApplication(ulf_application)
        print("KratosULFApplication sucessfully imported")

    if(Import_KratosMPISearchApplication):
        print("importing KratosMPISearchApplication ...")
        sys.path.append(applications_path + '/mpi_search_application/python_scripts')
        from KratosMPISearchApplication import *
        mpi_search_application = KratosMPISearchApplication()
        kernel.ImportApplication(mpi_search_application)
        print("KratosMPISearchApplication sucessfully imported")

    if(Import_KratosTrilinosApplication):
        print("importing KratosTrilinosApplication ...")
        sys.path.append(applications_path + '/trilinos_application/python_scripts')
        from KratosTrilinosApplication import *
        trilinos_application = KratosTrilinosApplication()
        kernel.ImportApplication(trilinos_application)
        print("KratosTrilinosApplication sucessfully imported")

    if(Import_KratosMetisApplication):
        print("importing KratosMetisApplication ...")
        sys.path.append(applications_path + '/metis_application/python_scripts')
        from KratosMetisApplication import *
        metis_application = KratosMetisApplication()
        kernel.ImportApplication(metis_application)
        print("KratosMetisApplication sucessfully imported")

    if(Import_PoissonApplication):
        print("importing KratosR1PoissonApplication ...")
        sys.path.append(applications_path + '/kPoisson/python_scripts')
        from KratosR1PoissonApplication import *
        kPoisson = KratosR1PoissonApplication()
        kernel.ImportApplication(kPoisson)
        print("Kratos PoissonApplication1 sucessfully imported")

    if(Import_ElectrostaticApplication):
        print("importing KratosR1ElectrostaticApplication ...")
        sys.path.append(applications_path + '/kElectrostatic/python_scripts')
        from KratosR1ElectrostaticApplication import *
        kElectrostatic = KratosR1ElectrostaticApplication()
        kernel.ImportApplication(kElectrostatic)
        print("Kratos ElectromagnecticApplication1 sucessfully imported")

    if(Import_MagnetostaticApplication):
        print("importing KratosR1MagnetostaticApplication ...")
        sys.path.append(applications_path + '/kMagnetostatic/python_scripts')
        from KratosR1MagnetostaticApplication import *
        kMagnetostatic = KratosR1MagnetostaticApplication()
        kernel.ImportApplication(kMagnetostatic)
        print("Kratos ElectromagneticApplication2 sucessfully imported")

    if(Import_DamApplication):
        print("importing KratosDamApplication ...")
        sys.path.append(applications_path + '/DamApplication/python_scripts')
        from KratosDamApplication import *
        dam_application = KratosDamApplication()
        kernel.ImportApplication(dam_application)
        print("Kratos DamApplication sucessfully imported")

    if(Import_TestApplication):
        print("importing KratosTestApplication ...")
        sys.path.append(applications_path + '/TestApplication/python_scripts')
        from KratosTestApplication import *
        test_application = KratosTestApplication()
        kernel.ImportApplication(test_application)
        print("Kratos TestApplication sucessfully imported")

    if(Import_OpenCLApplication):
        print("importing KratosOpenCLApplication ...")
        sys.path.append(applications_path + '/OpenCLapplication/python_scripts')
        from KratosOpenCLApplication import *
        opencl_application = KratosOpenCLApplication()
        kernel.ImportApplication(opencl_application)
        print("KratosOpenCLApplication sucessfully imported")

    if(Import_PodApplication):
        print("importing KratosPodApplication ...")
        sys.path.append(applications_path + '/PODApplication/python_scripts')
        from KratosPodApplication import *
        pod_application = KratosPodApplication()
        kernel.ImportApplication(pod_application)
        print("KratosPodApplication sucessfully imported")

    if(Import_LevelSetApplication):
        print("importing KratosPodApplication ...")
        sys.path.append(applications_path + '/LevelSetApplication/python_scripts')
        from KratosLevelSetApplication import *
        levelset_application = KratosLevelSetApplication()
        kernel.ImportApplication(levelset_application)
        print("KratosLevelSetApplication sucessfully imported")

    if(Import_FluidDynamicsApplication):
        print("importing KratosFluidDynamicsApplication ...")
        sys.path.append(applications_path + '/FluidDynamicsApplication/python_scripts')
        from KratosFluidDynamicsApplication import *
        fluid_dynamics_application = KratosFluidDynamicsApplication()
        kernel.ImportApplication(fluid_dynamics_application)
        print("KratosFluidDynamicsApplication sucessfully imported")

    if(Import_KratosDEMApplication):
        print("importing KratosDEMApplication ...")
        sys.path.append(applications_path + '/DEM_application/python_scripts')
        from KratosDEMApplication import *
        DEM_application = KratosDEMApplication()
        kernel.ImportApplication(DEM_application)
        print("KratosDEMApplication sucessfully imported")

    if(Import_KratosSwimmingDEMApplication):
        print("importing KratosSwimmingDEMApplication ...")
        sys.path.append(applications_path + '/swimming_DEM_application/python_scripts')
        from KratosSwimmingDEMApplication import *
        swimming_DEM_application = KratosSwimmingDEMApplication()
        kernel.ImportApplication(swimming_DEM_application)
        print("KratosSwimmingDEMApplication sucessfully imported")

    if(Import_KratosMixedElementApplication):
        print("importing KratosMixedElementApplication ...")
        sys.path.append(applications_path + '/MixedElementApplication/python_scripts')
        from KratosMixedElementApplication import *
        mixed_element_application = KratosMixedElementApplication()
        kernel.ImportApplication(mixed_element_application)
        print("KratosMixedElementApplication sucessfully imported")

    if(Import_ThermoMechanicalApplication):
        print("importing ThermoMechanicalApplication ...")
        sys.path.append(applications_path + '/ThermoMechanicalApplication/python_scripts')
        from KratosThermoMechanicalApplication import *
        thermo_mechanical_application = KratosThermoMechanicalApplication()
        kernel.ImportApplication(thermo_mechanical_application)
        print("KratosThermoMechanicalApplication sucessfully imported")

    if(Import_DEM_FEM_Application):
        print("importing DEM_FEM_Application ...")
        sys.path.append(applications_path + '/DEM_FEM_Application/python_scripts')
        from KratosDEM_FEM_Application import *
        dem_fem_application = KratosDEM_FEM_Application()
        kernel.ImportApplication(dem_fem_application)
        print("KratosDem_Fem_Application sucessfully imported")

    if(Import_MultiScaleApplication):
        print("importing KratosMultiscaleApplication ...")
        sys.path.append(applications_path + '/MultiScaleApplication/python_scripts')
        from KratosMultiscaleApplication import *
        multi_scale_application = KratosMultiScaleApplication()
        kernel.ImportApplication(multi_scale_application)
        print("KratosMultiScaleApplication sucessfully imported")

    if(Import_ContactStructuralMechanicsApplication):
        print("importing KratosStructuralContactMechanicsApplication ...")
        sys.path.append(applications_path + '/ContactStructuralMechanics/python_scripts')
        sys.path.append(applications_path + '/ContactStructuralMechanics/Linux')
        from KratosContactStructuralMechanicsApplication import *
        contact_structural_mechanics_application = KratosContactStructuralMechanicsApplication()
        kernel.ImportApplication(contact_structural_mechanics_application)
        print("KratosContactStructuralMechanicsApplication Succesfully imported")

    if(Import_KratosMappingApplication):
        print("importing KratosMappingApplication ...")
        sys.path.append(applications_path + '/MappingApplication/python_scripts')
        sys.path.append(applications_path + '/MappingApplication/Linux')
        from MappingApplication import *
        mapping_application = KratosMappingApplication()
        kernel.ImportApplication(mapping_application)
        print("KratosMappingApplication Succesfully imported")

    if(Import_ConstitutiveModelsApplication):
        print("importing KratosConstitutiveModelsApplication ...")
        sys.path.append(applications_path + '/ConstitutiveModels/python_scripts')
        sys.path.append(applications_path + '/ConstitutiveModels/Linux')
        from KratosConstitutiveModelsApplication import *
        constitutive_models_application = KratosConstitutiveModelsApplication()
        kernel.AddApplication(constitutive_models_application)
        print("KratosConstitutiveModelsApplication Succesfully imported")

    if(Import_ShallowWaterApplication):
        print("importing KratosShallowWaterApplication ...")
        sys.path.append(applications_path + '/ConstitutiveModels/python_scripts')
        sys.path.append(applications_path + '/ConstitutiveModels/Linux')
        from KratosShallowWaterApplication import *
        shallow_water_application = KratosShallowWaterApplication()
        kernel.AddApplication(shallow_water_application)
        print("KratosShallowWaterApplication Succesfully imported")

    if(Import_DelaunayMeshingApplication):
        print("importing KratosDelaunayMeshingApplication ...")
        sys.path.append(applications_path + '/DelaunayMeshing/python_scripts')
        sys.path.append(applications_path + '/DelaunayMeshing/Linux')
        from KratosDelaunayMeshingApplication import *
        delaunay_meshing_application = KratosDelaunayMeshingApplication()
        kernel.ImportApplication(delaunay_meshing_application)
        print("KratosDelaunayMeshingApplication Succesfully imported")

    if(Import_PoromechanicsApplication):
        print("importing KratosPoromechanicsApplication ...")
        sys.path.append(applications_path + '/Poromechanics/python_scripts')
        sys.path.append(applications_path + '/Poromechanics/Linux')
        from KratosPoromechanicsApplication import *
        poromechanics_application = KratosPoromechanicsApplication()
        kernel.AddApplication(poromechanics_application)
        print("KratosPoromechanicsApplication Succesfully imported")

    if(Import_FluidTransportApplication):
        print("importing KratosFluidTransportApplication ...")
        sys.path.append(applications_path + '/FluidTransport/python_scripts')
        sys.path.append(applications_path + '/FluidTransport/Linux')
        from KratosFluidTransportApplication import *
        fluid_transport_application = KratosFluidTransportApplication()
        kernel.AddApplication(fluid_transport_application)
        print("KratosFluidTransportApplication Succesfully imported")

    # dynamic renumbering of variables to ensure the consistency
    kernel.Initialize()

    if(Import_FluidRveLagrangeMultipliersApplication):
        kernel.InitializeApplication(fluid_rve_lagrange_multipliers_application)
    if(Import_SolidMechanicsApplication):
        kernel.InitializeApplication(solid_mechanics_application)
    if(Import_PfemApplication):
        kernel.InitializeApplication(pfem_application)
    if(Import_PfemSolidMechanicsApplication):
        kernel.InitializeApplication(pfem_solid_mechanics_application)
    if(Import_MachiningApplication):
        kernel.InitializeApplication(machining_application)
    if(Import_PfemFluidDynamicsApplication):
        kernel.InitializeApplication(pfem_fluid_dynamics_application)
    if(Import_StringDynamicsApplication):
        kernel.InitializeApplication(string_dynamics_application)
    if(Import_MeshMovingApplication):
        kernel.InitializeApplication(mesh_moving_app)
    if(Import_IncompressibleFluidApplication):
        kernel.InitializeApplication(incompressible_fluid_application)
    if(Import_StructuralApplication):
        kernel.InitializeApplication(structural_application)
    if(Import_StructuralMechanicsApplication):
        kernel.InitializeApplication(structural_mechanics_application)
    if(Import_ConvectionDiffusionApplication):
        kernel.InitializeApplication(convection_diffusion_application)
    if(Import_FSIApplication):
        kernel.InitializeApplication(fsi_application)
    if(Import_ExternalSolversApplication):
        kernel.InitializeApplication(external_solvers_application)
    if(Import_ConstitutiveLawsApplication):
        kernel.InitializeApplication(constitutive_laws_application)
    if(Import_ULFApplication):
        kernel.InitializeApplication(ulf_application)
    if(Import_MeshingApplication):
        kernel.InitializeApplication(meshing_application)
    if(Import_KratosMPISearchApplication):
        kernel.InitializeApplication(mpi_search_application)
    if(Import_KratosTrilinosApplication):
        kernel.InitializeApplication(trilinos_application)
    if(Import_KratosMetisApplication):
        kernel.InitializeApplication(metis_application)
    if(Import_PoissonApplication):
        kernel.InitializeApplication(kPoisson)
    if(Import_ElectrostaticApplication):
        kernel.InitializeApplication(kElectrostatic)
    if(Import_MagnetostaticApplication):
        kernel.InitializeApplication(kMagnetostatic)
    if(Import_DamApplication):
        kernel.InitializeApplication(dam_application)
    if(Import_TestApplication):
        kernel.InitializeApplication(test_application)
    if(Import_OpenCLApplication):
        kernel.InitializeApplication(opencl_application)
    if(Import_PodApplication):
        kernel.InitializeApplication(pod_application)
    if(Import_LevelSetApplication):
        kernel.InitializeApplication(levelset_application)
    if(Import_FluidDynamicsApplication):
        kernel.InitializeApplication(fluid_dynamics_application)
    if(Import_KratosDEMApplication):
        kernel.InitializeApplication(DEM_application)
    if(Import_KratosSwimmingDEMApplication):
        kernel.InitializeApplication(swimming_DEM_application)
    if(Import_KratosMixedElementApplication):
        kernel.InitializeApplication(mixed_element_application)
    if(Import_ThermoMechanicalApplication):
        kernel.InitializeApplication(thermo_mechanical_application)
    if(Import_DEM_FEM_Application):
        kernel.InitializeApplication(dem_fem_application)
    if(Import_MultiScaleApplication):
        kernel.InitializeApplication(MultiScaleApplication)
    if(Import_ContactMechanicsApplication):
        kernel.InitializeApplication(contact_mechanics_application)
    if(Import_ContactStructuralMechanicsApplication):
        kernel.InitializeApplication(contact_structural_mechanics_application)
    if(Import_KratosMappingApplication):
        kernel.InitializeApplication(mapping_application)
    if(Import_ConstitutiveModelsApplication):
        kernel.InitializeApplication(constitutive_models_application)
    if(Import_ShallowWaterApplication):
        kernel.InitializeApplication(shallow_water_application)
    if(Import_DelaunayMeshingApplication):
        kernel.InitializeApplication(delaunay_meshing_application)
    if(Import_PoromechanicsApplication):
        kernel.InitializeApplication(poromechanics_application)
    if(Import_FluidTransportApplication):
        kernel.InitializeApplication(fluid_transport_application)

# def ImportApplications(kernel  ):
    # import os.path
    # application_directory = os.path.dirname( os.path.realpath(__file__)  )
    # ImportApplications(kernel, application_directory )
