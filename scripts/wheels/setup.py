import setuptools
import os
import json
import glob
import shutil

kratos_version = os.environ["KRATOS_VERSION"]

def replaceKeyword(str):
    return str.replace("${KRATOS_VERSION}", kratos_version).replace("${PYTHON}", os.environ["PYTHON"])

def replaceKeywords(stringArray):
    return list(map(lambda str: replaceKeyword(str), stringArray))


with open("wheel.json", "r") as conf_file:
    conf = json.loads(conf_file.read())

with open(os.path.join(os.environ["KRATOS_ROOT"], conf["readme"]), "r") as fh:
    long_description = fh.read()

for module in conf["included_modules"]:
    src = os.path.join(os.environ["KRATOS_ROOT"], "bin", "Release", replaceKeyword("python_${PYTHON_VERSION}"), "KratosMultiphysics", module)
    dst = os.path.join("KratosMultiphysics", module)
    if os.path.exists(dst):
        shutil.rmtree(dst)
    try:
        shutil.copytree(src, dst)
    except Exception as e:
        print("Warning copying {}: {}".format(src,e))

for binary in conf["included_binaries"]:
    for file in glob.glob(os.path.join(os.environ["KRATOS_ROOT"], "bin", "Release", replaceKeyword("python_${PYTHON_VERSION}"), "libs", replaceKeyword(binary))):
        print("Adding {} matching binary: {}".format(file, binary))
        shutil.copy(file, os.path.join("KratosMultiphysics", ".libs"))

if "excluded_binaries" in conf:
    f = open("excluded.txt", "w")
    f.writelines(list(map(lambda x: replaceKeyword(x) + "\n", conf["excluded_binaries"])))
    f.close()

package_classifiers = [
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Other Audience",
        "Intended Audience :: Developers",
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console"
    ]

if "licence" not in conf.keys():
    package_classifiers.append("License :: OSI Approved :: BSD License")
else:
    package_classifiers.append(conf["licence"])

class BinaryDistribution(setuptools.Distribution):
    def has_ext_modules(foo):
        return True

class EmptyListWithLength(list):
    def __len__(self):
        return 1

setuptools.setup(
    name=conf["wheel_name"],
    version=kratos_version,
    author=conf["author"],
    author_email=conf["author_email"],
    description=conf["description"],
    url="https://github.com/KratosMultiphysics/",
    packages=setuptools.find_namespace_packages(),
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=replaceKeywords(conf["dependencies"]),
    classifiers=package_classifiers,
    package_data={
        'KratosMultiphysics': list(map(lambda x: ".libs/" + x, os.listdir("KratosMultiphysics/.libs")))
    },
    python_requires='>=3.8',
    ext_modules=EmptyListWithLength(),
    distclass=BinaryDistribution
)
