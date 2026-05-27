import os
import sys 

from hatchling.builders.hooks.plugin.interface import BuildHookInterface
from hatchling.metadata.plugin.interface import MetadataHookInterface

class CustomHook(BuildHookInterface):
    def initialize(self, version, build_data):
        # Generate Wheel tags. This is necessary because Kratos is platform specific.
        build_data["pure_pyton"] = False
        py_version = f"{sys.version_info.major}{sys.version_info.minor}"
        
        if sys.platform.startswith("linux"):
            # Linux: Force the PyPI-compliant manylinux tag
            build_data["infer_tag"] = False
            build_data["tag"] = f"cp{py_version}-cp{py_version}-manylinux_2_34_x86_64"
        elif sys.platform == "darwin":
            build_data["infer_tag"] = True
        else:
            build_data["infer_tag"] = True

class CustomMeta(MetadataHookInterface):
    def update(self, metadata):
        # Get the version
        version = os.environ["KRATOS_VERSION"]

        # Fill version globaly for all projects
        metadata["version"] = version

        # Add [Extras] info
        if self.metadata.name == "KratosMultiphysics":
            metadata["optional-dependencies"] = {
                "FluidDynamicsApplication": [f"KratosFluidDynamicsApplication >= {version}"],
                "StructuralMechanicsApplication" : [f"KratosStructuralMechanicsApplication >= {version}"],
                "DEMApplication": [f"KratosDEMApplication >= {version}"],
                "ContactStructuralMechanicsApplication": [f"KratosContactStructuralMechanicsApplication >= {version}"],
                "MPMApplication": [f"KratosMPMApplication >= {version}"],
                "ConvectionDiffusionApplication" : [f"KratosConvectionDiffusionApplication >= {version}"],
                "DamApplication": [f"KratosDamApplication >= {version}"],
                "PoromechanicsApplication": [f"KratosPoromechanicsApplication >= {version}"],
                "FSIApplication": [f"KratosFSIApplication >= {version}"],
                "LinearSolversApplication": [f"KratosLinearSolversApplication >= {version}"],
                "ConstitutiveLawsApplication": [f"KratosConstitutiveLawsApplication >= {version}"],
                "MeshingApplication": [f"KratosMeshingApplication >= {version}"],
                "MetisApplication": [f"KratosMetisApplication >= {version}"],
                "DemStructuresCouplingApplication": [f"KratosDemStructuresCouplingApplication >= {version}"],
                "MeshMovingApplication": [f"KratosMeshMovingApplication >= {version}"],
                "ShapeOptimizationApplication": [f"KratosShapeOptimizationApplication >= {version}"],
                "CoSimulationApplication": [f"KratosCoSimulationApplication >= {version}"],
                "CableNetApplication": [f"KratosCableNetApplication >= {version}"],
                "RANSApplication": [f"KratosRANSApplication >= {version}"],
                "MappingApplication": [f"KratosMappingApplication >= {version}"],
                "HDF5Application": [f"KratosHDF5Application >= {version}"],
                "MedApplication": [f"KratosMedApplication >= {version}"],
                "IgaApplication": [f"KratosIgaApplication >= {version}"],
                "ChimeraApplication": [f"KratosChimeraApplication >= {version}"],
                "StatisticsApplication": [f"KratosStatisticsApplication >= {version}"],
                "RomApplication": [f"KratosRomApplication >= {version}"],
                "ShallowWaterApplication": [f"KratosShallowWaterApplication >= {version}"],
                "OptimizationApplication": [f"KratosOptimizationApplication >= {version}"],
                "GeoMechanicsApplication": [f"KratosGeoMechanicsApplication >= {version}"],
                "SystemIdentificationApplication": [f"KratosSystemIdentificationApplication >= {version}"],
                "all": [
                    f"KratosContactStructuralMechanicsApplication >= {version}",
                    f"KratosConvectionDiffusionApplication >= {version}",
                    f"KratosConstitutiveLawsApplication >= {version}",
                    f"KratosCoSimulationApplication >= {version}",
                    f"KratosDEMApplication >= {version}",
                    f"KratosDamApplication >= {version}",
                    f"KratosFluidDynamicsApplication >= {version}",
                    f"KratosFSIApplication >= {version}",
                    f"KratosLinearSolversApplication >= {version}",
                    f"KratosMappingApplication >= {version}",
                    f"KratosMeshingApplication >= {version}",
                    f"KratosMPMApplication >= {version}",
                    f"KratosPoroMechanicsApplication >= {version}",
                    f"KratosShallowWaterApplication >= {version}",
                    f"KratosStructuralMechanicsApplication >= {version}",
                    f"KratosGeoMechanicsApplication >= {version}"
                ]
            }