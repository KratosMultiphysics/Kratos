---
title: How to run multiple cases using PyCOMPSs
keywords: 
tags: [How-to-run-multiple-cases-using-PyCOMPSs.md]
sidebar: kratos_for_developers
summary: 
---

# Overview
This tutorial gives a brief but exhaustive overview about the integration of PyCOMPSs inside the Kratos environment in order to launch multiple cases of a problem concurrently.

# Content
* [What is PyCOMPSs?][what]
* [How to install PyCOMPSs?][how-to-install]
* [How can I use this library?][how-to-use]
	* [Problem definition][problem-def]
	* [Integration of PyCOMPSs][integration-pycompss]
		* [Task definition][task]
		* [Use and abuse of compss_wait_on][compss_wait_on]
		* [Running with PyCOMPSs][runcompss]
* [Important observations][observations]

[what]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-run-multiple-cases-using-PyCOMPSs#what-is-pycompss
[how-to-install]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-run-multiple-cases-using-PyCOMPSs#how-to-install-pycompss
[how-to-use]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-run-multiple-cases-using-PyCOMPSs#how-can-i-use-this-library
[problem-def]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-run-multiple-cases-using-PyCOMPSs#problem-definition
[integration-pycompss]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-run-multiple-cases-using-PyCOMPSs#integration-of-pycompss
[task]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-run-multiple-cases-using-PyCOMPSs#task-definition
[compss_wait_on]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-run-multiple-cases-using-PyCOMPSs#use-and-abuse-of-compss_wait_on
[runcompss]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-run-multiple-cases-using-PyCOMPSs#running-with-pycompss
[observations]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-run-multiple-cases-using-PyCOMPSs#important-observations

# What is PyCOMPSs?
COMP Superscalar [(COMPSs)](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar) is a framework developed at Barcelona Supercomputing Center ([BSC](https://www.bsc.es/)) whose aim is to ease the running of applications in distributed environments.
Exploiting this structure, the developer can program following the sequential programming paradigm and avoid considering parallelization, distribution, data distribution, etc.
PyCOMPSs is the python library required in order to use COMPSs in a python environment.

# How to install PyCOMPSs?
The following section guides the user toward the installation of the PyCOMPSs library.
Check the presence of `apt-get` and all the usual Ubuntu installers. In general these commands should not be needed.
```console
# clean APT
sudo -E apt-get update ; sudo apt-get install -y --no-install-recommends apt-utils
```
COMPSs needs Java 8, with other versions does not work. So we uninstall any existing Java version and we install Java 8 to be sure that the used Java is the one with version 8.
```console
# remove any Java version
dpkg-query -W -f='${binary:Package}\n' | grep -E -e '^(ia32-)?(sun|oracle)-java' -e '^openjdk-' -e '^icedtea' -e '^(default|gcj)-j(re|dk)' -e '^gcj-(.*)-j(re|dk)' | xargs sudo apt-get -y remove
```
Install all the dependencies.
```console
# install basic COMPSs dependencies
sudo apt-get -y --no-install-recommends install openjdk-8-jre openjdk-8-jdk
sudo apt-get -y --no-install-recommends install python
sudo apt-get -y --no-install-recommends install maven subversion
sudo apt-get -y --no-install-recommends install graphviz xdg-utils
sudo apt-get -y --no-install-recommends install libtool automake build-essential
sudo bash -c "export DEBIAN_FRONTEND=noninteractive && apt-get install -y --no-install-recommends openssh-server openssh-client"
sudo apt-get -y --no-install-recommends install libxml2 libxml2-dev gfortran libpapi-dev papi-tools
sudo apt-get -y --no-install-recommends install openmpi-bin openmpi-doc libopenmpi-dev uuid-runtime curl bc git
sudo apt-get install libboost-dev
```
Download and install COMPSs.
```console
# setup COMPSs version and path
compss_folder_name="compss"
compss_path="$HOME"/"${compss_folder_name}"
# setup bash environment
grep -v "JAVA_HOME" "$HOME"/.bashrc > "$HOME"/newbashrc
echo "export JAVA_HOME=\"/usr/lib/jvm/java-8-openjdk-amd64/\"" >> "$HOME"/newbashrc
echo "export PYTHONPATH=$PYTHONPATH:\"/opt/COMPSs/Bindings/python/3/\"" >> "$HOME"/newbashrc
echo "source /etc/profile.d/compss.sh" >> "$HOME"/newbashrc
mv "$HOME"/newbashrc "$HOME"/.bashrc
export JAVA_HOME="/usr/lib/jvm/java-8-openjdk-amd64/"
# download COMPSs
sudo rm -rf "${compss_path}"
cd "$HOME"
git clone https://github.com/bsc-wdc/compss.git -b 2.8 "${compss_folder_name}"
cd "${compss_path}"
./submodules_get.sh
./submodules_patch.sh
cd -
# install COMPSs without monitor (-M) and autoparallel (-A)
cd "${compss_path}"/builders
sudo -E ./buildlocal -M -A
cd -
```

In case of problem during the installation, try to add the flag `-K (stream backend)` to the `buildlocal`. In other words, use:
```console
sudo -E ./buildlocal -M -A -K
```

In addition to the previous steps, we need to be able to `ssh` to localhost without password in the local machine we run compss. We already have this in a cluster but in local it is not always the case. The steps to follow are the followings.
Check if we already have a ssh key:
```console
ls ~/.ssh
```
In case the file "id_rsa.pub" exists, that means that we already have ssh keys. Otherwise, we should run the following command:
```console
ssh-keygen
```
Once the "id_rsa.pub" file has been generated, we should run this command in order to be able to have passwordless ssh to localhost:
```console
cat ~/.ssh/id_rsa.pub >> ~/.ssh/authorized_keys
```
These final steps are needed because COMPSs needs passwordless ssh connection between all the nodes used in a computation.

A small observation now: in case you are compiling again, so it is not the first time you compile, and you get the error
```console
cp: cannot stat '"${compss_path}"/builders/tmp/compss/runtime/adaptors/agent/master/*.jar': No such file or directory
```
then you need to checkout the whole software on a different folder and run again the compilation.

# How can I use this library?
The PyCOMPSs library can be used inside a python environment to run a distributed simulation, thus it allows to launch multiple simulations, as we will see.
Since the goal of this section is to describe how to implement COMPSs inside Kratos, we choose to solve a very simple [problem](https://github.com/KratosMultiphysics/Documentation/tree/master/Wiki_files/How-to-run-multiple-cases-using-PyCOMPSs): a two-dimensional Poisson equation.

Suppose we have to run many times the same simulation changing only few parameters between different scenarios. To avoid launching many times the same simulation, that can be highly time-consuming, we can launch all the simulations required together exploiting PyCOMPSs.

## Problem definition
Let's consider the stationary heat equation with a varying heat flux, a square two-dimensional domain and Dirichlet boundary conditions. The problem is

<img src="https://render.githubusercontent.com/render/math?math=\nabla\cdot(K\nabla u)=\varepsilon f \ , u\in\Omega \\">

with boundary condition

<img src="https://render.githubusercontent.com/render/math?math=u=0 \ , u\in\partial(\Omega) \,,">

where <img src="https://render.githubusercontent.com/render/math?math=\Omega=[0,1]^{2}">, <img src="https://render.githubusercontent.com/render/math?math=f=-432(x^2+y^2-x-y)"> and <img src="https://render.githubusercontent.com/render/math?math=\varepsilon\sim~Beta(2,6)">, i.e. <img src="https://render.githubusercontent.com/render/math?math=\varepsilon"> follows a beta distribution. The thermal diffusivity is <img src="https://render.githubusercontent.com/render/math?math=K=1"> for simplicity. The Quantity of Interest (QoI) we are interested about is the integral over the whole domain of the temperature, meaning:

<img src="https://render.githubusercontent.com/render/math?math=QoI=\int_{\Omega}u(x,y)dxdy\,.">

In this problem we have a heat flux that varies following a beta distribution, but it is also possible to define a list of values for the heat flux and at each iteration to use one of them.

We report [here](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/How-to-run-multiple-cases-using-PyCOMPSs/launch-multiple-simulations-pycompss.py) the python file containing the `SimulationScenario` (the `AnalysisStage` of the problem), the main and the functions needed to run the tutorial.

We want to highlight that in the `SimulationScenario` class we added the `EvaluateQuantityOfInterest(self)` function with respect to the `AnalysisStage` base class. This function computes the QoI, given the results of the analysis. In addition, we see `ModifyInitialProperties(self)` modifies the property `KratosMultiphysics.HEAT_FLUX` of our PDE.

## Integration of PyCOMPSs
We start now analyzing how PyCOMPSs is defined inside our code. First of all we observe the import statements:

```python
from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import *
```

We import:
* `task`: the attribute we give to each function we want to launch in distributed environment,
* `compss_wait_on`: the command that brings the return of a task back to our local machine (explained deeper later),
* `parameter`: all the parameter attributes, required for example to read the [project_parameters.json](https://github.com/KratosMultiphysics/Documentation/tree/master/Wiki_files/How-to-run-multiple-cases-using-PyCOMPSs/problem_settings/project_parameters.json) file inside a task.

Observe that these three lines in the code are commented, since we import them with the single command
```python
from exaqute.ExaquteTaskPyCOMPSs import *
```


### Task definition
Each function that we want to run in distributed environment should be defined as a task. In our example we have two tasks:
```python
@ExaquteTask(parameter_file_name=FILE_IN,returns=2)
def SerializeModelParameters_Task(parameter_file_name):
```
and
```python
@ExaquteTask(returns=1)
def ExecuteInstance_Task(pickled_model,pickled_parameters,heat_flux_list,instance):
```

Many different attributes can be given to customize each task, in this tutorial we only describe the ones we used here. In a task definition we always have to define the number of returns of the task, e.g.
```python
returns=2
```
denotes that the number of returns of the task will be two. And these two will be `pycompss.runtime.binding.Future` objects, always!
On the other hand,
```python
parameter_file_name=FILE_IN
```
says that the string `parameter_file_name` is read-only.
For more parameter options we leave as a reference the [COMPSs user manual](http://compss.bsc.es/releases/compss/latest/docs/COMPSs_User_Manual_App_Development.pdf), chapter 3 "Python Binding".

### Use and abuse of compss_wait_on
The `compss_wait_on` command brings the `pycompss.runtime.binding.Future` back to its "real" nature, i.e. if in the task we compute `qoi` that is a float, and then we return this variable, `qoi` is a `pycompss.runtime.binding.Future`. Only the execution of the `compss_wait_on` command convert it from `pycompss.runtime.binding.Future` to `float`.
To clarify what we have just said, we report the following example taken from our problem.
From each `ExecuteInstance_Task` task we return a single float number, whose type is `pycompss.runtime.binding.Future`. We append each of these values in the `qoi` list, and we print type and value of the first member of `qoi` before and after the call of `compss_wait_on`.

```python
print(type(qoi[0]),qoi[0])
qoi = compss_wait_on(qoi)
print(type(qoi[0]),qoi[0])
```
gives
```python
<class 'pycompss.runtime.binding.Future'> <pycompss.runtime.binding.Future object at 0x7f7f1ff766a0>
<class 'float'> 1.14
```
where `1.14` is just the value of the variable for the first instance.
Therefore it is clear the importance of this command for the visualization of the results.

On the other hand, we should not abuse of the `compss_wait_on` command. In fact, what this command does is synchronizing `qoi` with our local machine, thus it is a synchronization point and we should minimize its usage and locate it as in the and as possible in our code to fully exploit the distribution potential!

### Running with PyCOMPSs
The problem of this tutorial can be directly run without difficulties using python3:
```console
python3 launch-multiple-simulations-pycompss.py
```

On the  other hand, to run with PyCOMPSs the command is the following:
```console
runcompss \
    --lang=python \
    --python_interpreter=python3 \
    --pythonpath=$TEST_DIR \
    --graph=false \
    --trace=false \
    $TEST_DIR/launch-multiple-simulations-pycompss.py
```

The graph of our simulation is
![graph]
and shows the parallelism of the execution. The blue circles are the `ExecuteInstance_Task` tasks, while the red hexagon is the synchronization point. To obtain the graph and the trace, one should run with `--graph=true` and `--trace=true`, respectively. The graph and the trace can be found in `~/.COMPSs/execution_id/monitor` and `~/.COMPSs/execution_id/trace`, respectively. In order to visualize the graph connections, we should execute `compss_gengraph complete_graph.dot`, which generates the pdf of the graph. In old PyCOMPSs versions the command was `gengraph complete_graph.dot`
To see the trace, let's open the `*.prv` file using the `wxparaver` software, which can be downloaded from the [BSC page](https://tools.bsc.es/downloads).

# Important observations
To launch simulation using PyCOMPSs we must always use the ABSOLUTE path referring to any file, that in our case are [project_parameters.json](https://github.com/KratosMultiphysics/Documentation/tree/master/Wiki_files/How-to-run-multiple-cases-using-PyCOMPSs/launch-multiple-simulations-pycompss.py#L134) and [model_part.mdpa](hhttps://github.com/KratosMultiphysics/Documentation/tree/master/Wiki_files/How-to-run-multiple-cases-using-PyCOMPSs/problem_settings/project_parameters.json#L18) files.
Thus we should have
```console
parameter_file_name = "/ABSOLUTE/PATH/TO/project_parameters.json"
```
and
```console
"input_filename": "/ABSOLUTE/PATH/TO/model_part"
```

Another important observation is that not all Kratos objects can be used as PyCOMPSs inputs, e.g., `KratosMultiphysics.Model` and `KratosMultiphysics.Parameters` classes cannot. To overcome this issue, we exploit the `KratosMultiphysics.StreamSerializer` class, which can serialize all the Kratos classes, and this class can be used as input in a task. In our problem, the function taking care of the serialization is `SerializeModelParameters_Task`.

[graph]: https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/How-to-run-multiple-cases-using-PyCOMPSs/complete_graph.png