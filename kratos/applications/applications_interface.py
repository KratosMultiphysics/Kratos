#path for applications
import sys

Import_ALEApplication = False
Import_IncompressibleFluidApplication = False
Import_StructuralApplication = False
Import_ConvectionDiffusionApplication = False
Import_FSIApplication = False
Import_PFEMApplication = False
Import_ExternalSolversApplication = False
Import_ULFApplication = False
Import_MeshingApplication = False
Import_KratosMKLSolversApplication = False
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
Import_KratosMixedElementApplication = False

print "Applications Available:"
print "Import_ALEApplication: False"
print "Import_IncompressibleFluidApplication: False"
print "Import_StructuralApplication: False"
print "Import_ConvectionDiffusionApplication: False"
print "Import_FSIApplication: False"
print "Import_ExternalSolversApplication: False"
print "Import_ULFApplication: False"
print "Import_MeshingApplication: False"
print "Import_KratosMKLSolversApplication: False"
print "Import_KratosTrilinosApplication: False"
print "Import_KratosMetisApplication: False"
print "Import_PoissonApplication: False"  
print "Import_ElectrostaticApplication: False"
print "Import_MagnetostaticApplication: False"
print "Import_DamApplication: False"
print "Import_TestApplication: False"
print "Import_OpenCLApplication: False"
print "Import_PodApplication: False"
print "Import_LevelSetApplication: False"
print "Import_FluidDynamicsApplication: False"
print "Import_KratosMixedElementApplication: False"


def ImportApplications(kernel, applications_path):
    ##########################################################################
    ##importing the applications
    print "Applications Available:"
    print "Import_ALEApplication: "+str(Import_ALEApplication)
    print "Import_IncompressibleFluidApplication: "+str(Import_IncompressibleFluidApplication)
    print "Import_StructuralApplication: "+str(Import_StructuralApplication)
    print "Import_ConvectionDiffusionApplication: "+str(Import_ConvectionDiffusionApplication)
    print "Import_FSIApplication: "+str(Import_FSIApplication)
    print "Import_ExternalSolversApplication: "+str(Import_ExternalSolversApplication)
    print "Import_ULFApplication: "+str(Import_ULFApplication)
    print "Import_MeshingApplication: "+str(Import_MeshingApplication)
    print "Import_KratosTrilinosApplication: "+str(Import_KratosTrilinosApplication)
    print "Import_KratosMetisApplication: "+str(Import_KratosMetisApplication)
    print "Import_PoissonApplication: "+str(Import_PoissonApplication)
    print "Import_ElectrostaticApplication: "+str(Import_ElectrostaticApplication)
    print "Import_MagnetostaticApplication: "+str(Import_MagnetostaticApplication)
    print "Import_DamApplication: "+str(Import_DamApplication)
    print "Import_TestApplication: "+str(Import_TestApplication)
    print "Import_OpenCLApplication: "+str(Import_OpenCLApplication)
    print "Import_PodApplication: "+str(Import_PodApplication)
    print "Import_LevelSetApplication:"+str(Import_LevelSetApplication)
    print "Import_FluidDynamicsApplication: "+str(Import_FluidDynamicsApplication)
    print "Import_KratosMixedElementApplication: "+str(Import_KratosMixedElementApplication)

    if(Import_ALEApplication == True):
        print "importing KratosALEApplication ..."
        sys.path.append(applications_path + '/ALEapplication/python_scripts' )
        sys.path.append(applications_path + '/ALEapplication/Linux')
        from KratosALEApplication import *
        ale_app = KratosALEApplication()
        kernel.AddApplication(ale_app)
        print "KratosALEApplication Succesfully imported"
        
    if(Import_IncompressibleFluidApplication == True):
        print "importing KratosIncompressibleFluidApplication ..."
        sys.path.append(applications_path + '/incompressible_fluid_application/python_scripts') 
        sys.path.append(applications_path + '/incompressible_fluid_application/Linux') 
        from KratosIncompressibleFluidApplication import *
        print "KratosIncompressibleFluidApplication lib loaded"
        incompressible_fluid_application = KratosIncompressibleFluidApplication()
        print "KratosIncompressibleFluidApplication application created"
        kernel.AddApplication(incompressible_fluid_application)
        print "KratosIncompressibleFluidApplication Succesfully imported"

    if(Import_StructuralApplication == True):
        print "importing KratosStructuralApplication ..."
        sys.path.append(applications_path + '/structural_application/python_scripts') 
        sys.path.append(applications_path + '/structural_application/Linux') 
        from KratosStructuralApplication import *
        structural_application = KratosStructuralApplication()
        kernel.AddApplication(structural_application)
        print "KratosStructuralApplication Succesfully imported"

    if(Import_ConvectionDiffusionApplication == True):
        print "importing KratosConvectionDiffusionApplication ..."
        sys.path.append(applications_path + '/convection_diffusion_application/python_scripts') 
        sys.path.append(applications_path + '/convection_diffusion_application/Linux') 
        from KratosConvectionDiffusionApplication import *
        convection_diffusion_application = KratosConvectionDiffusionApplication()
        kernel.AddApplication(convection_diffusion_application)
        print "KratosConvectionDiffusionApplication Succesfully imported"

    if(Import_FSIApplication == True):
        print "importing FSIapplication ..."
        sys.path.append(applications_path + '/FSIapplication/python_scripts') 
        sys.path.append(applications_path + '/FSIapplication/Linux') 
        from KratosFSIApplication import *
        fsi_application = KratosFSIApplication()
        kernel.AddApplication(fsi_application)
        print "FSIapplication Succesfully imported"

    if(Import_PFEMApplication == True):
        print "importing KratosPFEMApplication ..."
        sys.path.append(applications_path + '/PFEMapplication/python_scripts') 
        sys.path.append(applications_path + '/PFEMapplication/Linux') 
        from KratosPFEMApplication import *
        pfem_application = KratosPFEMApplication()
        kernel.AddApplication(pfem_application)
        print "KratosPFEMApplication Succesfully imported"
	
    if(Import_ExternalSolversApplication == True):
        print "importing KratosExternalSolversApplication ..."
        sys.path.append(applications_path + '/ExternalSolversApplication/python_scripts')
        sys.path.append(applications_path + '/ExternalSolversApplication/Linux')
        from KratosExternalSolversApplication import *
        external_solvers_application = KratosExternalSolversApplication()
        kernel.AddApplication(external_solvers_application)
        print "KratosExternalSolversApplication sucessfully imported"
 	
    if(Import_ULFApplication == True):
        print "importing KratosULFApplication ..."
        sys.path.append(applications_path + '/ULFapplication/python_scripts')
        sys.path.append(applications_path + '/ULFapplication/Linux')
        from KratosULFApplication import *
        ulf_application = KratosULFApplication()
        kernel.AddApplication(ulf_application)
        print "KratosULFApplication sucessfully imported"

    if(Import_MeshingApplication == True):
        print "importing KratosMeshingApplication ..."
        sys.path.append(applications_path + '/MeshingApplication/python_scripts')
        from KratosMeshingApplication import *
        meshing_application = KratosMeshingApplication()
        kernel.AddApplication(meshing_application)
        print "KratosMeshingApplication sucessfully imported"

    if(Import_KratosMKLSolversApplication == True):
        print "importing KratosMKLSolversApplication ..."
        sys.path.append(applications_path + '/mkl_solvers_application/python_scripts')
        from KratosMKLSolversApplication import *
        mkl_solvers_application = KratosMKLSolversApplication()
        kernel.AddApplication(mkl_solvers_application)
        print "KratosMKLSolversApplication sucessfully imported"

    if(Import_KratosTrilinosApplication == True):
        print "importing KratosTrilinosApplication ..."
        sys.path.append(applications_path + '/trilinos_application/python_scripts')
        from KratosTrilinosApplication import *
        trilinos_application = KratosTrilinosApplication()
        kernel.AddApplication(trilinos_application)
        print "KratosTrilinosApplication sucessfully imported"
        
    if(Import_KratosMetisApplication == True):
        print "importing KratosMetisApplication ..."
        sys.path.append(applications_path + '/metis_application/python_scripts')
        from KratosMetisApplication import *
        metis_application = KratosMetisApplication()
        kernel.AddApplication(metis_application)
        print "KratosMetisApplication sucessfully imported"

    if(Import_PoissonApplication == True):
        print "importing KratosR1PoissonApplication ..."
        sys.path.append(applications_path + '/kPoisson/python_scripts')
        from KratosR1PoissonApplication import *
        kPoisson = KratosR1PoissonApplication()
        kernel.AddApplication(kPoisson)
        print "Kratos PoissonApplication1 sucessfully imported"
        
    if(Import_ElectrostaticApplication == True):
        print "importing KratosR1ElectrostaticApplication ..."
        sys.path.append(applications_path + '/kElectrostatic/python_scripts')
        from KratosR1ElectrostaticApplication import *
        kElectrostatic = KratosR1ElectrostaticApplication()
        kernel.AddApplication(kElectrostatic)
        print "Kratos ElectromagnecticApplication1 sucessfully imported"

    if(Import_MagnetostaticApplication == True):
        print "importing KratosR1MagnetostaticApplication ..."
        sys.path.append(applications_path + '/kMagnetostatic/python_scripts')
        from KratosR1MagnetostaticApplication import *
        kMagnetostatic = KratosR1MagnetostaticApplication()
        kernel.AddApplication(kMagnetostatic)
        print "Kratos ElectromagneticApplication2 sucessfully imported"       
        
    if(Import_DamApplication == True):
        print "importing KratosDamApplication ..."
        sys.path.append(applications_path + '/DamApplication/python_scripts')
        from KratosDamApplication import *
        dam_application = KratosDamApplication()
        kernel.AddApplication(dam_application)
        print "Kratos DamApplication sucessfully imported"

    if(Import_TestApplication == True):
        print "importing KratosTestApplication ..."
        sys.path.append(applications_path + '/TestApplication/python_scripts')
        from KratosTestApplication import *
        test_application = KratosTestApplication()
        kernel.AddApplication(test_application)
        print "Kratos TestApplication sucessfully imported"

    if(Import_OpenCLApplication == True):
        print "importing KratosOpenCLApplication ..."
        sys.path.append(applications_path + '/OpenCLapplication/python_scripts')
        from KratosOpenCLApplication import *
        opencl_application = KratosOpenCLApplication()
        kernel.AddApplication(opencl_application)
        print "KratosOpenCLApplication sucessfully imported"

    if(Import_PodApplication == True):
        print "importing KratosPodApplication ..."
        sys.path.append(applications_path + '/PODApplication/python_scripts')
        from KratosPodApplication import *
        pod_application = KratosPodApplication()
        kernel.AddApplication(pod_application)
        print "KratosPodApplication sucessfully imported"
    if(Import_LevelSetApplication == True):
        print "importing KratosPodApplication ..."
        sys.path.append(applications_path + '/LevelSetApplication/python_scripts')
        from KratosLevelSetApplication import *
        levelset_application = KratosLevelSetApplication()
        kernel.AddApplication(levelset_application)
        print "KratosLevelSetApplication sucessfully imported"
    if(Import_FluidDynamicsApplication == True):
        print "importing KratosFluidDynamicsApplication ..."
        sys.path.append(applications_path + '/FluidDynamicsApplication/python_scripts')
        from KratosFluidDynamicsApplication import *
        fluid_dynamics_application = KratosFluidDynamicsApplication()
        kernel.AddApplication(fluid_dynamics_application)
        print "KratosFluidDynamicsApplication sucessfully imported"
    if(Import_KratosMixedElementApplication == True):
        print "importing KratosMixedElementApplication ..."
        sys.path.append(applications_path + '/MixedElementApplication/python_scripts')
        from KratosMixedElementApplication import *
        mixed_element_application = KratosMixedElementApplication()
        kernel.AddApplication(mixed_element_application)
        print "KratosMixedElementApplication sucessfully imported"
        
    ##########################################################################
    ##dynamic renumbering of variables to ensure the consistency
    kernel.Initialize()
    if(Import_ALEApplication == True):
        kernel.InitializeApplication(ale_app);
    if(Import_IncompressibleFluidApplication == True):
        kernel.InitializeApplication(incompressible_fluid_application);    
    if(Import_StructuralApplication == True):
        kernel.InitializeApplication(structural_application);    
    if(Import_ConvectionDiffusionApplication == True):
        kernel.InitializeApplication(convection_diffusion_application); 
    if(Import_FSIApplication == True):
        kernel.InitializeApplication(fsi_application); 
    if(Import_PFEMApplication == True):
        kernel.InitializeApplication(pfem_application); 
    if(Import_ExternalSolversApplication == True):
        kernel.InitializeApplication(external_solvers_application);
    if(Import_ULFApplication == True):
        kernel.InitializeApplication(ulf_application);
    if(Import_MeshingApplication == True):
        kernel.InitializeApplication(meshing_application);
    if(Import_KratosMKLSolversApplication == True):
        kernel.InitializeApplication(mkl_solvers_application);
    if(Import_KratosTrilinosApplication == True):
        kernel.InitializeApplication(trilinos_application);
    if(Import_KratosMetisApplication == True):
        kernel.InitializeApplication(metis_application);
    if(Import_PoissonApplication == True):
        kernel.InitializeApplication(kPoisson);
    if(Import_ElectrostaticApplication == True):
        kernel.InitializeApplication(kElectrostatic);
    if(Import_MagnetostaticApplication == True):
        kernel.InitializeApplication(kMagnetostatic);
    if(Import_DamApplication == True):
        kernel.InitializeApplication(dam_application);
    if(Import_TestApplication == True):
        kernel.InitializeApplication(test_application);
    if(Import_OpenCLApplication == True):
        kernel.InitializeApplication(opencl_application);
    if(Import_PodApplication == True):
        kernel.InitializeApplication(pod_application);
    if(Import_LevelSetApplication == True):
        kernel.InitializeApplication(levelset_application);
    if(Import_FluidDynamicsApplication == True):
        kernel.InitializeApplication(fluid_dynamics_application);
    if(Import_KratosMixedElementApplication == True):
        kernel.InitializeApplication(mixed_element_application);
