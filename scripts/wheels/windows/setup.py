import setuptools
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

class BinaryDistribution(setuptools.Distribution):
    def has_ext_modules(foo):
        return True


setuptools.setup(
    name="KratosMultiphysics",
    version="7.0-" + os.environ['HASH'],
    author="Kratos Team",
    author_email="kratos@listas.cimne.upc.edu",
    description="KRATOS Multiphysics (\"Kratos\") is a framework for building parallel, multi-disciplinary simulation software, aiming at modularity, extensibility, and high performance. Kratos is written in C++, and counts with an extensive Python interface.",
    url="https://github.com/KratosMultiphysics/",
    packages=setuptools.find_packages(),
    long_description=long_description,
    long_description_content_type="text/markdown",
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
    distclass=BinaryDistribution
)
