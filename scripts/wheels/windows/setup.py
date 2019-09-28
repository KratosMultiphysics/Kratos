import setuptools
import os

#
# with open("README.md", "r") as fh:
#     long_description = fh.read()

class BinaryDistribution(setuptools.Distribution):
    def has_ext_modules(foo):
        return True


setuptools.setup(
    name="KratosMultiphysicsWheelTest",
    version="0.0.10",

    author="Example Author",
    author_email="author@example.com",
    description="Kratos Multiphisics wheel test",
    url="https://github.com/KratosMultiphysics/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    package_data={
        'KratosMultiphysics': list(map(lambda x: "libs/" + x, os.listdir("KratosMultiphysics/libs")))
    },
    python_requires='>=3.5',
    distclass=BinaryDistribution
)
