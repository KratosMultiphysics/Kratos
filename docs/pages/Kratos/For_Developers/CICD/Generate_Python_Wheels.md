---
title: How to generate python wheels
keywords: 
tags: [How-to-generate-python-wheels.md]
sidebar: kratos_for_developers
summary: 
---

## Package metadata
It is possible go define package-specific metadata, like author, author email, description, README file. 
To do that, edit wheel definitions in `scripts/wheels/<platform>/applications/<application_name>.json` in the Kratos repository.
It's important to edit both linux and windows definitions, since they are treated as separate packages.

### Kratos version
Before each new release Kratos version must be changed. This is done by modifying KRATOS_VERSION variables at the top of [`scripts/wheels/windows/build.ps1`](https://github.com/KratosMultiphysics/Kratos/blob/master/scripts/wheels/windows/build.ps1) and [`scripts/wheels/linux/build.sh`](https://github.com/KratosMultiphysics/Kratos/blob/master/scripts/wheels/linux/build.sh). 
Allowed version strings are defined by the [PEP440](https://www.python.org/dev/peps/pep-0440/).

### Special wheels
Build process for KratosCore and Kratos-all wheels is slighty different. Beacuse of that, their definition files are located outside of the `applications` directory.
You can find them in `scripts/wheels/windows/KratosMultiphysics.json`, `scripts/wheels/linux/KratosMultiphysics.json` and `scrips/wheels/linux/KratosMultiphysics-all.json`
Json structure is exactly the same, as for other applications.

### Package definition example
```json
{
    "wheel_name": "KratosContactStructuralMechanicsApplication",
    "included_modules": ["ContactStructuralMechanicsApplication"],
    "included_binaries": ["KratosContactStructuralMechanicsApplication.cpython-${PYTHON}m-x86_64-linux-gnu.so"],
    "excluded_binaries": ["libKratosStructuralMechanicsCore"],
    "dependencies": ["KratosMultiphysics==${KRATOS_VERSION}", "KratosStructuralMechanicsApplication==${KRATOS_VERSION}"],
    "author": "Vicente Mataix Ferrándiz",
    "author_email": "vmataix@cimne.upc.edu",
    "description": "KRATOS Multiphysics (\"Kratos\") is a framework for building parallel, multi-disciplinary simulation software, aiming at modularity, extensibility, and high performance. Kratos is written in C++, and counts with an extensive Python interface.",
    "readme": "scripts/wheels/README.md"
}
```

`wheel_name` - this is the wheel name. Will be passed to `pip install` command when instaling this package.

`included_modules` - included python modules. Those are copied from `KratosMultiphysics/<module_name>`

`included_binaries` - binary files copied from `libs` directory. On windows you must include all required binaries and their dependencies. On linux binary dependencies are managed automatically and in most cases it's enough to specify the binary python module. Supports keywords.

`excluded_binaries` - On linux only, optional. Used to specify all binary dependencies, that are **guaranteed** to be already provided, for example by another Kratos application, which is included in `dependencies`. **Skip file extension.** You don't have to specify here binaries provided by Kratos Core wheel - those are filtered automatically. Supports keywords.

`dependecies` - names and versions of required python packages, as they are available on the PyPi. Supports keywords. For more information visit [PEP440 - version specifiers](https://www.python.org/dev/peps/pep-0440/#version-specifiers).

`author` - application/package author name

`author_email` - application/package author contact email. Will be public

`description` - short package description

`readme` - path to the readme file in the Kratos repository

### Keywords
In some places in the json file you may use keywords, which will be replaced with correct information during build. Available keywords:
* `${PYTHON}` - two digits and optional 'm' representing currently used python version. For example `35m` or `38`.
* `${KRATOS_VERSION}` - String representing currently built Kratos version. Wheels will be tagged using that version. Useful when definig dependencies on other Kratos wheels.

## Adding new application
1. Make sure, that every application, that will be included in the new wheel package is being built by the wheel generation scripts. 
This can be checked and eventually modified in `scripts/wheels/windows/configure.bat` and `scripts/wheels/linux/configure.sh` files.
2. Define json files. Place them in `scripts/wheels/<platform>/applications` directory for each platform, that you want to support. For more information on json structure visit "Package metadata" section.
3. Push your changes to any branch in the Kratos repository.
4. That's it! Now you can build your wheels as described in the "Wheel generation section". 
Remember, that due to binary dependecies between wheels, all of them must be built without any code changes in the Kratos Core or any other dependencies.

### Depending on another Kratos application
On windows this process is quite easy. All you have to do is add dependent wheel name in the `dependencies` array in the json file. Remember to set dependency version using `${KRATOS_VERSION}` keyword.

On linux you have to also blacklist all binaries provided by the dependency app in `excluded_binaries` array, so they don't get packaged twice. You may use `${PYTHON}` and `${KRATOS_VERSION}` keywords.


For a working example check `ContactStructuralMechanicsApplication.json`

## Wheel generation

### Linux
1. Install docker.
2. Run command `sudo docker run --rm -v local_output_dir:/workspace/out kratosmultiphysics/kratos-wheelbuilder-linux [OPTION]...`

Available options:

`--branch branch_name` - branch to compile. By default `master`.

`--cpus number_of_cores` - number of cores assigned to the build. By default `4`.

`--cotire ON/OFF` - enable/disable cotire build. By default `OFF`.

`--pythons python_versions` - coma separated python versions to target. By default `"35,36,37,38"`.

`--repository git_repository_url` - repository to use during build. By default `https://github.com/KratosMultiphysics/Kratos.git`.

3. Wait. Wheel files will appear in `local_output_dir`.

### Windows
1. Install docker and enable windows containers.
2. Run command `docker run --rm --memory=<memory> --cpus=<number_of_cores> -v local_output_dir:c:\out kratosmultiphysics/kratos-wheelbuilder-windows [OPTION]...`

`<memory>` - maximum amount of RAM assigned to the container. Recomended value is `8g`.

`<number_of_cores>` - number of cores assigned to the container. By default `1`.

Available options:

`-branch` - branch to compile. By default `master`.

`-cotire ON/OFF` - enable/disable cotire build. By default `OFF`.

`-pythons python_versions` - coma separated python versions to target. By default `"35,36,37,38"`.

`-cpus number_of_cores` - number of cores assigned to the build. This is not the same as `--cpus` option provided by docker. By default the same as <number_of_cores>.

3. Wait. Wheel files will appear in `local_output_dir`.

### MacOS

As Apple Inc. Does not provide a reliable virtualization solution capable of our scaling needs, currently MacOS releases are not supported in our pipelines.

## Build environment

Docker images used in wheel generation process are available on the docker hub of Kratos:
* [Linux](https://hub.docker.com/r/kratosmultiphysics/kratos-wheelbuilder-linux)
* [Windows](https://hub.docker.com/r/kratosmultiphysics/kratos-wheelbuilder-windows)

The corresponding Dockerfiles are available here: 
* [Linux](https://github.com/KratosMultiphysics/Kratos/tree/master/scripts/docker_files/docker_file_wheelbuilder_linux)
* [Windows](https://github.com/KratosMultiphysics/Kratos/tree/master/scripts/docker_files/docker_file_wheelbuilder_windows)