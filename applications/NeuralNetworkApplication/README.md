## Neural Network Application

The Neural Network Application contains the interface for implementing surrogate models based on neural networks in Kratos Multiphysics.

### Requirements

The following main external packages are required:

  - [Tensorflow](https://www.tensorflow.org/): at least version 2.0, tested on version 2.7.0.
  - [Keras](https://keras.io/): tested on version 2.7.0.
 

### Features

  - Dataset generation utilities integrated with other Applications.
  - Generation of neural network models based on Keras API.
  - Training, testing and hyperparameter tuning of neural network models.
  - Use of pre-trained models as subtitutes of the solver in an analysis stage.
  - Implementation of a pre-trained surrogate model to be used in a multiphysics simulation (see CoSimulationApplication)

### Workflow

#### Data generation
To generate a suitable dataset, insert a [data_generator_process](python_scripts/data_generator_process.py) in the JSON file for the project parameters and 
adjust the variables as needed. The data generated is in a suitable format for training a neural network. Neuarl Network application supports the perturbation
of variables by a random noise.

#### Neural network modeling
This application supports neural network models created externally as long as they are created with Keras. It also gives support for creating a neural network
model sequentially layer by layer through the application itself. Most of the basic layers and functionalities present in Keras are currently available
with identical parameters as the original counterparts. Consult the documentation on Keras layers for more information [here](https://keras.io/api/layers/).

#### Neural network training, testing and hyperaparameter tuning
Each of the different modes can be chosen in the project parameters, and allow to generate a working neural network model. The model can be used stand-alone or
inside Kratos Multiphysics. For hyperparameter tuning, the alternatives present in Keras are available.

#### Use in multiphysics simulation
For the use of pre-trained neural networks as surrogate models in multiphysics simulations, it is recommended the use of CoSimulationApplication. The integration
is achieved through a wrapper adapted to that application.

### Examples
There are examples in the [Examples](https://github.com/KratosMultiphysics/Examples) repository from KratosMultiphysics project in the corresponding branch. 
Currently, there are examples available for:

  - Static structural mechanics
  - Dynamic structural mechanics
  - Fluid dynamics
  - FSI

The "structural_full_example" is significatively representative of the full workflow for a simple model. Additionally, the proposed solution for Mok's benchmark 
shows the use of a surrogate model based on neural networks in a co-simulation case.
