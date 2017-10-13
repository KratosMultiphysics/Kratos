from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# path for applications
import sys
import os.path

Import_SolidMechanicsApplication = False
Import_PfemApplication = False
Import_PfemSolidMechanicsApplication = False
Import_MachiningApplication = False
Import_PfemFluidDynamicsApplication = False
Import_StringDynamicsApplication = False
Import_ExternalSolversApplication = False
Import_ContactMechanicsApplication = False
Import_ConstitutiveModelsApplication = False

print("Applications Available:")
print("Import_SolidMechanicsApplication: False")
print("Import_PfemApplication: False")
print("Import_PfemSolidMechanicsApplication: False")
print("Import_MachiningApplication: False")
print("Import_PfemFluidDynamicsApplication: False")
print("Import_StringDynamicsApplication: False")
print("Import_ExternalSolversApplication: False")
print("Import_ContactMechanicsApplication: False")
print("Import_ConstitutiveModelsApplication: False")

application_directory = os.path.dirname(os.path.realpath(__file__))

def ImportApplications(kernel, applications_path=application_directory):
    # importing the applications
    print("Applications Available:")
    print("Import_SolidMechanicsApplication: " + str(Import_SolidMechanicsApplication))
    print("Import_PfemApplication: " + str(Import_PfemApplication))
    print("Import_PfemSolidMechanicsApplication: " + str(Import_PfemSolidMechanicsApplication))
    print("Import_MachiningApplication: " + str(Import_MachiningApplication))
    print("Import_PfemFluidDynamicsApplication: " + str(Import_PfemFluidDynamicsApplication))
    print("Import_StringDynamicsApplication: " + str(Import_StringDynamicsApplication))
    print("Import_ExternalSolversApplication: " + str(Import_ExternalSolversApplication))
    print("Import_ContactMechanicsApplication: " + str(Import_ContactMechanicsApplication))   
    print("Import_ConstitutiveModelsApplication: " + str(Import_ConstitutiveModelsApplication))

    if(Import_SolidMechanicsApplication):
        print("importing KratosSolidMechanicsApplication ...")
        sys.path.append(applications_path + '/SolidMechanicsApplication/python_scripts')
        sys.path.append(applications_path + '/SolidMechanicsApplication/python_scripts/constitutive_laws')
        sys.path.append(applications_path + '/SolidMechanicsApplication/Linux')
        from KratosSolidMechanicsApplication import *
        solid_mechanics_application = KratosSolidMechanicsApplication()
        kernel.AddApplication(solid_mechanics_application)
        print("KratosSolidMechanicsApplication Succesfully imported")

    if(Import_PfemApplication):
        print("importing KratosPfemApplication ...")
        sys.path.append(applications_path + '/PfemApplication/python_scripts')
        sys.path.append(applications_path + '/PfemApplication/Linux')
        from KratosPfemApplication import *
        pfem_application = KratosPfemApplication()
        kernel.AddApplication(pfem_application)
        print("KratosPfemApplication Succesfully imported")

    if(Import_PfemSolidMechanicsApplication):
        print("importing KratosPfemSolidMechanicsApplication ...")
        sys.path.append(applications_path + '/PfemSolidMechanicsApplication/python_scripts')
        sys.path.append(applications_path + '/PfemSolidMechanicsApplication/Linux')
        from KratosPfemSolidMechanicsApplication import *
        pfem_solid_mechanics_application = KratosPfemSolidMechanicsApplication()
        kernel.AddApplication(pfem_solid_mechanics_application)
        print("KratosPfemSolidMechanicsApplication Succesfully imported")

    if(Import_MachiningApplication):
        print("importing KratosMachiningApplication ...")
        sys.path.append(applications_path + '/MachiningApplication/python_scripts')
        sys.path.append(applications_path + '/MachiningApplication/Linux')
        from KratosMachiningApplication import *
        machining_application = KratosMachiningApplication()
        kernel.AddApplication(machining_application)
        print("KratosMachiningApplication Succesfully imported")

    if(Import_PfemFluidDynamicsApplication):
        print("importing KratosPfemFluidDynamicsApplication ...")
        sys.path.append(applications_path + '/PfemFluidDynamicsApplication/python_scripts')
        sys.path.append(applications_path + '/PfemFluidDynamicsApplication/Linux')
        from KratosPfemFluidDynamicsApplication import *
        pfem_fluid_dynamics_application = KratosPfemFluidDynamicsApplication()
        kernel.AddApplication(pfem_fluid_dynamics_application)
        print("KratosPfemFluidDynamicsApplication Succesfully imported")

    if(Import_StringDynamicsApplication):
        print("importing KratosStringDynamicsApplication ...")
        sys.path.append(applications_path + '/StringDynamicsApplication/python_scripts')
        sys.path.append(applications_path + '/StringDynamicsApplication/Linux')
        from KratosStringDynamicsApplication import *
        pfem_solid_mechanics_application = KratosStringDynamicsApplication()
        kernel.AddApplication(string_dynamics_application)
        print("KratosStringDynamicsApplication Succesfully imported")

    if(Import_ExternalSolversApplication):
        print("importing KratosExternalSolversApplication ...")
        sys.path.append(applications_path + '/ExternalSolversApplication/python_scripts')
        sys.path.append(applications_path + '/ExternalSolversApplication/Linux')
        from KratosExternalSolversApplication import *
        external_solvers_application = KratosExternalSolversApplication()
        kernel.AddApplication(external_solvers_application)
        print("KratosExternalSolversApplication sucessfully imported")
		
    if(Import_ContactMechanicsApplication):
        print("importing KratosContactMechanicsApplication ...")
        sys.path.append(applications_path + '/ContactMechanics/python_scripts')
        sys.path.append(applications_path + '/ContactMechanics/Linux')
        from KratosContactMechanicsApplication import *
        contact_mechanics_application = KratosContactMechanicsApplication()
        kernel.AddApplication(contact_mechanics_application)
        print("KratosContactMechanicsApplication Succesfully imported")
        
    if(Import_ConstitutiveModelsApplication):
        print("importing KratosConstitutiveModelsApplication ...")
        sys.path.append(applications_path + '/ConstitutiveModels/python_scripts')
        sys.path.append(applications_path + '/ConstitutiveModels/Linux')
        from KratosConstitutiveModelsApplication import *
        constitutive_models_application = KratosConstitutiveModelsApplication()
        kernel.AddApplication(constitutive_models_application)
        print("KratosConstitutiveModelsApplication Succesfully imported")


    # dynamic renumbering of variables to ensure the consistency
    kernel.Initialize()
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
    if(Import_ExternalSolversApplication):
        kernel.InitializeApplication(external_solvers_application)
    if(Import_ContactMechanicsApplication):
        kernel.InitializeApplication(contact_mechanics_application)
    if(Import_ConstitutiveModelsApplication):
        kernel.InitializeApplication(constitutive_models_application)

# def ImportApplications(kernel  ):
    # import os.path
    # application_directory = os.path.dirname( os.path.realpath(__file__)  )
    # ImportApplications(kernel, application_directory )
