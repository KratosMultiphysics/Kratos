import setuptools
import os
import json


with open("wheel.json", "r") as conf_file:
    conf = json.loads(conf_file.read())

with open(os.path.join(os.environ["KRATOS_ROOT"], conf["readme"]), "r") as fh:
    long_description = fh.read()


import shutil

for module in conf["included_modules"]:
    shutil.copytree(os.path.join(os.environ["KRATOS_ROOT"], "KratosMultiphysics", module), os.path.join("KratosMultiphysics", module))

for binary in conf["included_binaries"]:
    shutil.copy(os.path.join(os.environ["KRATOS_ROOT"], "libs", binary.replace("${PYTHON}", os.environ["PYTHON"])), os.path.join("KratosMultiphysics", ".libs"))

class BinaryDistribution(setuptools.Distribution):
    def has_ext_modules(foo):
        return True

class EmptyListWithLength(list):
    def __len__(self):
        return 1

setuptools.setup(
    name=conf["wheel_name"],
    version="7.0+" + os.environ['HASH'],
    author=conf["author"],
    author_email=conf["author_email"],
    description=conf["description"],
    url="https://github.com/KratosMultiphysics/",
    packages=setuptools.find_namespace_packages(),
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=conf["dependencies"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Mathematics",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Other Audience",
        "Intended Audience :: Developers"
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console"
    ],
    package_data={
        'KratosMultiphysics': list(map(lambda x: ".libs/" + x, os.listdir("KratosMultiphysics/.libs")))
    },
    python_requires='>=3.5',
    ext_modules=EmptyListWithLength(),
    distclass=BinaryDistribution
)
