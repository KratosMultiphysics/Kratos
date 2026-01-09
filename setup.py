from setuptools import setup
from setuptools.command.install import install

import sys
import platform

def set_install_requirements():
    """
    Creates kratos requirements list
    """
    release_tag = "test_release0.1.2"
    kratos_version = "10.3.1.3"
    python_version_part = ""
    platform_part = ""

    # get platform part of wheel the name
    if (sys.platform == "win32"):
        platform_part = "-win_amd64.whl"
    elif (sys.platform == "linux"):
        platform_part = "-manylinux_2_34_x86_64.whl"

    # get python version part of the wheel name

    kratos_multiphysics = "KratosMultiphysics"
    kratos_linear_solvers_application = "KratosLinearSolversApplication"
    kratos_structural_mechanics_application = "KratosStructuralMechanicsApplication"
    kratos_geo_mechanics_application = "KratosGeoMechanicsApplication"
    kratos_railway_application = "KratosRailwayApplication"

    if (platform.python_version().startswith("3.10.")):
        python_version_part = '-cp310-cp310'
    elif (platform.python_version().startswith("3.11.")):
        python_version_part = '-cp311-cp311'
    elif (platform.python_version().startswith("3.12.")):
        python_version_part = '-cp312-cp312'

        kratos_multiphysics = "kratosmultiphysics"
        kratos_linear_solvers_application = "kratoslinearsolversapplication"
        kratos_structural_mechanics_application = "kratosstructuralmechanicsapplication"
        kratos_geo_mechanics_application = "kratosgeomechanicsapplication"
        kratos_railway_application = "kratosrailwayapplication"


    requirements = [
        f"KratosMultiphysics @ https://github.com/StemVibrations/Kratos/releases/download/{release_tag}/{kratos_multiphysics}-{kratos_version}{python_version_part}{platform_part}",
        f"KratosLinearSolversApplication @ https://github.com/StemVibrations/Kratos/releases/download/{release_tag}/{kratos_linear_solvers_application}-{kratos_version}{python_version_part}{platform_part}",
        f"KratosStructuralMechanicsApplication @ https://github.com/StemVibrations/Kratos/releases/download/{release_tag}/{kratos_structural_mechanics_application}-{kratos_version}{python_version_part}{platform_part}",
        f"KratosGeoMechanicsApplication @ https://github.com/StemVibrations/Kratos/releases/download/{release_tag}/{kratos_geo_mechanics_application}-{kratos_version}{python_version_part}{platform_part}",
        f"KratosRailwayApplication @ https://github.com/StemVibrations/Kratos/releases/download/{release_tag}/{kratos_railway_application}-{kratos_version}{python_version_part}{platform_part}",
        f"orjson==3.11.0",
        f"numpy>=2.0.2"
                    ]

    return requirements



if __name__ == '__main__':
    setup(
        install_requires=set_install_requirements()
    )
