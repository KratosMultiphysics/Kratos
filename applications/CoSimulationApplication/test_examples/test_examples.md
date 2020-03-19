# Test Examples

This documentation describes the different test examples.
Currently all these examples calculate the flow in a flexible tube.


## Folder and file structure

This section describes the different folders and files which are provided.

- MainKratos.py: Main file which has to be run with a parameterfile as argument
- project_parameters_X.json: Parameter file in json format
- setup_X: Setup folder containing all files for setting up solver X
- setup_X.sh: Bash file which has to be run to set up solver X
- readme.md: A description of the specific example

When the setup files are run, working directories are created which have to match the ones specified in the parameter file.
These folder are expandable and are deleted when the setup files are (re)run.

## Running a case

In order to run a test example, the case have to be made ready by running the setup files for both solvers.
Then, the calculation is started by running MainKratos.py with the parameter file as argument.

## Debug files
The folder `test_examples` also contains a folder `debug_files` with files for debug purposes. In order to use those,
the debug boolean `self.debug` has to be `True` in the corresponding solver wrappers.
