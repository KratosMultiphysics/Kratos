## CoSimulation Application PYTHON ONLY VERSION

This is the description for the python-only version of the CoSimulationApplication.
Kratos is only used as a **container for data** within the CoSimulationApplication, the CoSimApp itself is coded in python.
For the python-only version, Kratos is replaced by pyKratos, which is a python version of Kratos which implements the basic features (of KratosMultiphysics) required for CoSimulation. The version of pyKratos within this folder is based on the [work of Riccardo Rossi](https://github.com/RiccardoRossi/pyKratos)

In order to install the python-version, please follow the instructions in [configure_py_co_sim.sh](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/CoSimulationApplication/custom_data_structure/configure_py_co_sim.sh). For this **no compilation** is required.

The tests can be executed by running `python test_CoSimulationApplication.py` in the folder `tests`
