---
title: Postprocessing
keywords: 
tags: [Postprocessing.md]
sidebar: kratos_workshops
summary: 
---
# Postprocessing
Upon finishing the simulation, the results should be processed in order to be comprehensibly displayed in the final report. This part of the guide will help you with the technicalities of the postprocessing.

___
## 1. Postprocessing in Paraview
Paraview can't read h5 files by default. In order to visualize the result, we need to generate an ".xdmf" file. Follow these steps to generate it:

- Load the default Kratos version at the cluster:
```shell
$ module load kratos/Master
```
- Navigate to the folder where you have the h5 files with the "cd" command.
- Create the ".xdmf" file by running:
```shell
$ convertH5toXdmf <name_of_files_until_dash>
```
- An ".xdmf" file will be created. Copy all the h5 files and the xdmf file to your computer to visualize. 

**Note that the xdmf file only serves as a link between the Paraview and the h5 files. Paraview need both the h5 and xdmf files.* 

After the files are ready to be read in Paraview, here are some helpful commands when using the software:
- **Paraview → Edit → Settings → General**, search for “Cache” and tick “Cache Geometry For Animation” 
to speed up picture creation for animation.
- **Import flow domain**: Velocity in model part “FluidModelPart.fluid_computational_model_part”:
  - Create **slices** &rarr; define a plane by its origin (coordinates) and the normal (coordinates or select normal).
  - Use **stream tracer** to visualize velocity (one of the symbols above the pipeline browser).
  - **Set opacity of the flow domain** so that slices, stream tracer and structure can be seen. Adjust lighting according to preference.
  - Vorticity visualization using the **Q-criterion** (later use iso-surface for Q-values in range of 0.2 - 0.01 1/s): follow [this link](https://discourse.paraview.org/t/qcriterion-in-paraview/2355){:target="_blank"} for a guide in using Q-criterion in Paraview.
- **Import pressure on Structure**: Pressure in model part  “FluidModelPart.NoSlip3D_structure” (or similar ModelPart name instead of structure):
  -  Same principles regarding creating slices and visualization as with the velocity.
- Animations:
  - Using save animation will create standalone JPG or PNG pictures for the results above.
  - Merge pictures into: GIF with the software of your choice **or/also** any video format with the software of your choice.

____
## 2. Point and line ASCII data
The ["point_output_process"](Preprocessing.html#21-point-output-process){:target="_blank"} and ["line_output_process"](Preprocessing.html#22-line-output-process){:target="_blank"} that you have defined in ProjectParametersCustom.json generates an ascii output (.dat files), with time series of respective pressure and velocities. We recommend (and support) you to create your own pythons scripts with numpy and matplotlib to visualize the data. Here's an example of a python script to generate a simple plot of the pressure of a certain point:

```python
import numpy as np

import matplotlib.pyplot as plt

path = <path_to_the_file_including_the_me>

data = np.loadtxt(path)

time = data[:,0]

pressure = data[:,1]

plt.plot(time, pressure)
```

This script reads a certain ascii output (in this case: pressure output) and then splits the file in two variables (time and pressure). Time is the first column of the ascii file, while pressure the second column. Then you can plot both of them to create a plot with the time series of the pressure output of the simulation. Similarly, you can do the same for the global forces (base forces and moments). 

With Numpy, you also have useful commands, which can help you find mean, maximum values etc. Combine these different commands to extract the results of interest. Visit this page on [Matplotlib](https://matplotlib.org/){:target="_blank"} for more details of the features available to visualize your data.

Other softwares are also available, such as Excel or Matlab, however we recommend and support Python.

____
## 3. Converting level forces to ParOptBeam

For the level forces, you need to convert them to a format that ParOptBeam understands, using the python script (**convert_kratos_to_paroptbeam.py**) that is available in the sample files we provided you. The script reads the level force file from the respective folder and converts them accordingly. For this, make sure that under `"n_level_files"` in the script, you input the [number of level forces received from the simulation](Preprocessing.html#23-force-output-process){:target="_blank"}. You can also sample the level forces for the ParOptbeam format in more intervals than the current level forces, under `"number_of_sampling_interval_cases"`. Run the script to receive the level force output, which will then serve as an input for the CSD in ParOptBeam.

More details on setting up the beam model and the CSD parameters for ParOptBeam will be explained in [the next chapter](ParOptBeam_Guide.html){:target="_blank"}.
