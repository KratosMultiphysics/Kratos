## ExxonMobil Tutorial for the Poromechanics application in Kratos 

Authors: Ignasi de Pouplana Guillermo Casas Danilo B. Cavalcanti Karen Lorena Casallas ~~OO~~ May 2023 

## Contents 

1. Installing 

2. Example 1: Terzaghi’s consolidation problem 

3. Example 2: Adding a discontinuity 

4. Example 3: 3D problem 

2 

Kratos – Poromechanics application 

Installing 

## Contents 

## 1. Kratos 

## 2. GUI (GiD + Poromechanics problemtype ) 

4 

## Kratos 

## https://github.com/KratosMultiphysics/Kratos 

5 

## Installing dependencies 

## Download and install: 

 Git; 

- To clone Kratos repository, open a command prompt window, navigate to the location where you will store the program, and type the following: 

git clone https://github.com/KratosMultiphysics/Kratos Kratos 

- Visual Studio 2019 or higher; 

   - Install it with Desktop development with C++ workload. 

- CMake 3.20 or higher; 

- Python 3.8 or higher; 

   - In the installation, add Python to the system PATH. 

- Boost 1.78.0 

   - extract it outside the Kratos folder. 

6 

## Setting configuration file to compile 

## Kratos 

1. Go to the folder where Kratos was installed; 

2. Open the scripts folder 

7 

## Creating the configure.bat file 

- ⮚ Make a copy of the standad_configure.bat file 

- ⮚ Here, we named the copy 

kratos_poromechanics_debug.bat 

8 

@ echo off 

rem Please do not modify this script rem For any question please contact with us in: rem - https://github.com/KratosMultiphysics/Kratos rem Optional parameters: rem You can find a list will all the compiation options in INSTALL.md or here: rem - https://github.com/KratosMultiphysics/Kratos/wiki/Compilation-options rem Set compiler set CC=cl.exe set CXX=cl.exe 

## Copy the text of this slide to your configure.bat file a 

## rem Set variables 

## Example of Release configuration file 

if not defined KRATOS_SOURCE set KRATOS_SOURCE= %~dp0 .. if not defined KRATOS_BUILD set KRATOS_BUILD= %KRATOS_SOURCE% /build set KRATOS_APP_DIR=applications 

rem Warning: In windows this option only works if you run through a terminal with admin privileges set KRATOS_INSTALL_PYTHON_USING_LINKS=ON 

rem Set basic configuration set KRATOS_BUILD_TYPE=Release set BOOST_ROOT= path_to_boost \boost_1_78_0 set PYTHON_EXECUTABLE= path_to_python \python3.9\python.exe 

Edit these lines. Change the **path_to_boost** and **path_to_python** to the folders where each one of them were installed. 

rem Set applications to compile set KRATOS_APPLICATIONS= 

CALL :add_app %KRATOS_APP_DIR% \ PoromechanicsApplication; 

## rem Clean 

del /F /Q " %KRATOS_BUILD% \ %KRATOS_BUILD_TYPE% \cmake_install.cmake" del /F /Q " %KRATOS_BUILD% \ %KRATOS_BUILD_TYPE% \CMakeCache.txt" del /F /Q " %KRATOS_BUILD% \ %KRATOS_BUILD_TYPE% \CMakeFiles" 

## rem Configure 

## @ echo on 

cmake -G"Visual Studio 16 2019" -H" %KRATOS_SOURCE% " -B" %KRATOS_BUILD% \ %KRATOS_BUILD_TYPE% " ^ -DUSE_EIGEN_MKL=OFF ^ 

-DFORCE_LOCAL_ZLIB_COMPILATION=ON ^ 

- -DCMAKE_CXX_FLAGS=" /MP3 " 

## rem Build 

cmake --build " %KRATOS_BUILD% / %KRATOS_BUILD_TYPE% " --target install -- /property:configuration= %KRATOS_BUILD_TYPE% /p:Platform=x64 goto:eof 

rem Function to add apps :add_app 

set KRATOS_APPLICATIONS=% KRATOS_APPLICATIONS% % 1; goto:eof 

9 

@ echo off 

rem Please do not modify this script rem For any question please contact with us in: rem - https://github.com/KratosMultiphysics/Kratos rem Optional parameters: rem You can find a list will all the compiation options in INSTALL.md or here: rem - https://github.com/KratosMultiphysics/Kratos/wiki/Compilation-options rem Set compiler set CC=cl.exe set CXX=cl.exe 

## Copy the text of this slide to your configure.bat file Le 

## rem Set variables 

## Example of Debug configuration file 

if not defined KRATOS_SOURCE set KRATOS_SOURCE= %~dp0 .. if not defined KRATOS_BUILD set KRATOS_BUILD= %KRATOS_SOURCE% /build set KRATOS_APP_DIR=applications 

rem Warning: In windows this option only works if you run through a terminal with admin privileges set KRATOS_INSTALL_PYTHON_USING_LINKS=ON 

rem Set basic configuration set KRATOS_BUILD_TYPE=FullDebug set BOOST_ROOT= path_to_boost \boost_1_78_0 set PYTHON_EXECUTABLE= path_to_python \python3.9\python.exe 

Edit these lines. Change the **path_to_boost** and **path_to_python** to the folders where each one of them were installed. 

rem Set applications to compile 

## set KRATOS_APPLICATIONS= 

CALL :add_app %KRATOS_APP_DIR% \ PoromechanicsApplication; 

## rem Clean 

del /F /Q " %KRATOS_BUILD% \ %KRATOS_BUILD_TYPE% \cmake_install.cmake" del /F /Q " %KRATOS_BUILD% \ %KRATOS_BUILD_TYPE% \CMakeCache.txt" del /F /Q " %KRATOS_BUILD% \ %KRATOS_BUILD_TYPE% \CMakeFiles" 

## rem Configure 

## @ echo on 

cmake -G"Visual Studio 16 2019" -H" %KRATOS_SOURCE% " -B" %KRATOS_BUILD% \ %KRATOS_BUILD_TYPE% " ^ -DINSTALL_RUNKRATOS=OFF ^ 

-DUSE_EIGEN_MKL=OFF ^ 

-DFORCE_LOCAL_ZLIB_COMPILATION=ON ^ 

-DCMAKE_CXX_FLAGS=" /MP3 " 

## rem Build 

cmake --build " %KRATOS_BUILD% / %KRATOS_BUILD_TYPE% " --target install -- /property:configuration= %KRATOS_BUILD_TYPE% /p:Platform=x64 goto:eof 

## rem Function to add apps 

## :add_app 

set KRATOS_APPLICATIONS=% KRATOS_APPLICATIONS% % 1; goto:eof 

10 

## Run the command prompt as administrator 

11 

## Set up the path 

1. Open the command prompt; 

2. Change the current path to the one where the gid model is: 

   - Write “cd” and the PATH to the scripts folder 

   - Press Enter 

12 

## Run the bat file 

## 1. Write the name of the file.bat; 

## 2. Press **ENTER** . 

13 

## Compilation output 

## Check if Python and Boost libraries were found 

14 

## Compilation output 

Check the applications 

15 

## Compilation output 

16 

## Compilation output 

## End of the compilation 

17 

## Environment variables 

## 1. Open System Properties; 2. Click in Environment variable 

18 

## Environment variables 

## In the system variables : 1. Edit the Path; 

## 2. Add the folder Kratos\bin\FullDebug\libs here 

19 

## Environment variables 

## In the system variables : 

## 1. Create the variable PYTHONPATH 2. Add the folder Kratos\bin\FullDebug here 

20 

## Testing Kratos 

## 1. Open a CMD window, and type the following: 

## python 

## 2. Press ENTER ; 

## 3. Once inside python, type the following: import KratosMultiphysics 

## 4. Press ENTER ; 5. You should see the KRATOS logo. 

21 

## GUI (GiD + Poromechanics problemtype ) 

## The last release version of GiD14 is still compatible with the poromechanics problemtype 

https://downloads.gidsimulation.com/#GiD_Official_Versions/Windows/win-x64/Old/ 

22 

## KratosPoroGiDInterface 

To clone KratosPoroGiDInterface repository, open a command prompt window, navigate to the location where you will store the program, and type the following: 

_git clone https://github.com/ipouplana/KratosPoroGiDInterface KratosPoroGiDInterface_ 

https://github.com/ipouplana/KratosP ———— oroGiDInterface.git 

23 

## Adding the Poromechanics 

## problemtype 

1. Go the folder which the repository was cloned; 

2. Create a shortcut for the “Poromechanics_Application.gid ” folder; 

3. Move the Poromechanics_Application.gid shortcut to the problemtypes GiD folder 

24 

## To execute Kratos from GiD 

## 3. Edit the file Poromechanics_application.win.bat 

rem USED BY GiD Substitute rem OutputFile: %1.info rem ErrorFile: %1.err rem Information rem basename = %1 rem currentdirectory = %2 rem problemtypedirectory = %3 rem Setting paths. WARNING: one should check them before running this file set PATH= Kratos_PATH \\Kratos\\bin\\Release\\libs ;%PATH% set PYTHONPATH= Kratos_PATH \\Kratos\\bin\\Release ;%PYTHONPATH% rem Execute the program %1. info 2> %1 .err Python_PATH \\python3.9\\python.exe MainKratos.py > 

Substitute **Kratos_PATH** and **Python_PATH** according to where those were installed. 

25 

## Using the Poromechanics problemtype 

## 4. Open GiD and load the Poromechanics Application problemtype 

26 

## Changing the branch 

27 

## Updating a branch 

28 

## Updating a branch with the master 

Remember that everytime you update the branch or make a change in the C++ code , you must recompile. 

29 

## Updating the KratosPoroGiDInterface 

30 

## Updating the problemtype 

31 

## Updating the problemtype 

32 

Kratos – Poromechanics application 

Example 01 

## Problem description 

## Terzaghi’s consolidation problem 

**==> picture [736 x 333] intentionally omitted <==**

**----- Start of picture text -----**<br>
w<br>0<br>p =<br>E  = 1000 kPa K K<br>x  =   = 1 m/d<br>1.0 y<br>ν  = 0.3<br>γ f  = 10 kN/m³<br>0<br>q =<br>eee<br>1.0<br>0 0<br>q =  q =<br>**----- End of picture text -----**<br>


34 

## Open GiD 

35 

## Creating the geometry 

36 

## Creating the geometry 

37 

## Creating the surface 

38 

## Creating the surface 

39 

## Creating the parts 

40 

## Creating the material 

41 

## Assign the material 

**==> picture [9 x 11] intentionally omitted <==**

**----- Start of picture text -----**<br>
4<br>**----- End of picture text -----**<br>


42 

## Assign the boundary conditions 

43 

## Boundary condition at the bottom edge 

44 

## Boundary condition at the lateral edges 

45 

## Boundary condition at the top edge 

46 

## Boundary condition at the top edge 

47 

## Assigning the loads 

48 

## Select face load 

49 

## Assign the load 

50 

## Solver 

51 

## Setting up the solver 

**==> picture [195 x 62] intentionally omitted <==**

**----- Start of picture text -----**<br>
<br>FIC Stabilization;<br> Non-local Damage;<br><br>.<br>Fracture propagation<br>**----- End of picture text -----**<br>


52 

## Setting up the solver 

53 

## Parallel configuration 

54 

## Select the output parameters 

55 

## Configuration of the mesh generation 

56 

## Structured mesh 

7 subdivisions on each 

57 

## Select the element type 

58 

## Mesh generation 

59 

## Generated mesh 

60 

## Running Kratos with GiD 

61 

## Running Kratos with GiD 

62 

## Visualizing the results 

67 

## Post -process visualization 

68 

## Selecting the output field to visualize 

69 

## Vertical displacement field 

70 

## Water pore -pressure profile 

71 

Kratos – Poromechanics application 

Example 2 

## Problem description 

## Terzaghi’s consolidation problem 

_E_ = 1000 kPa ν = 0.3 

_K K x_ = = 1 m/d _y_ γ _f_ = 10 kN/m³ 

**==> picture [737 x 333] intentionally omitted <==**

**----- Start of picture text -----**<br>
w<br>0<br>p =<br>1.0<br>0<br>q =<br>ee<br>0.5 0.5<br>0 0<br>q =  q =<br>**----- End of picture text -----**<br>


73 

## Create the domain 

74 

## Creating the interface 

75 

## Creating the interface 

76 

## Creating interface 

77 

## Intersect lines 

78 

## Intersect lines 

79 

## Contact surface 

80 

## Contact surface 

- < * 

- 1. Select contact surface creation; 

81 

## Creating the domain surface 

82 

## Creating the domain surface 

83 

## Creating the domain surface 

84 

## Creating the domain surface 

85 

## Creating the material 

86 

## Creating the material of the interface 

87 

## Interface material 

88 

## Boundary condition at the bottom edge 

89 

## Boundary conditions at the lateral edges 

90 

## Boundary condition at the top edge 

91 

## Load at the top edge 

92 

## Generate the mesh 

93 

## Run Kratos with GiD 

94 

## Run Kratos with GiD 

95 

## Water field pressure 

96 

Kratos – Poromechanics application 

Example 3 

**==> picture [960 x 492] intentionally omitted <==**

**----- Start of picture text -----**<br>
3D problem with a fracture<br>eee<br>w<br>V0 p46<br>p = 0<br>1.0<br>z<br>1.0<br>y<br>x<br>ea u<br>1.0<br>**----- End of picture text -----**<br>


98 

## Top view 

## Terzaghi’s consolidation problem 

_E_ = 1000 kPa 

ν = 0.3 

_K K x_ = = 1 m/d _y_ γ _f_ = 10 kN/m³ 

0 _q =_ 

~~e~~ @ 0.5 0.5 ~~.~~ @ e ~~pod~~ . _y_ 

0 _q =_ 

99 

## Creating the model 

Open the poromechanics problemtype; ’ 

100 

## Creating the model 

101 

## Creating the model 

102 

## Creating the model 

103 

## Creating the model 

104 

## Creating the model 

105 

## Creating the model 

106 

## Creating the model 

107 

## Creating the model 

108 

## Creating the model 

109 

## Viewing the 3D model 

110 

## Creating the contact volume 

111 

## Create and assign groups 

Create the following groups: 

112 

## Checking the assigned groups 

113 

## Checking the assigned groups 

114 

## Checking the assigned groups 

115 

## Checking the assigned groups 

116 

## Checking the assigned groups 

117 

## Checking the assigned groups 

118 

## Checking the assigned groups 

119 

## Creating the material 

120 

## Assigning the material 

**==> picture [216 x 15] intentionally omitted <==**

**----- Start of picture text -----**<br>
Select the BulkDomain group and  ;-S<br>**----- End of picture text -----**<br>


121 

## Creating the interface material 

|Select the 3D elastic cohesive<br>law;<br>Fill the material properties;<br>Select the InterfaceVol group<br>and assign it<br>i<br>a<br>/<br>/<br>||
|---|



122 

## Assigning the boundary conditions 

Select the BottomFace group and ; 

123 

## Assigning the boundary conditions 

**==> picture [219 x 15] intentionally omitted <==**

**----- Start of picture text -----**<br>
Select the LateralXDirFaces group  ,<br>**----- End of picture text -----**<br>


124 

## Assigning the boundary conditions 

**==> picture [219 x 15] intentionally omitted <==**

**----- Start of picture text -----**<br>
Select the LateralYDirFaces group  .<br>**----- End of picture text -----**<br>


125 

## Assigning the boundary conditions 

126 

## Assignin the load 

127 

## Problem data definition 

128 

## Creating a structured mesh 

129 

## Mesh generated 

130 

## Run Kratos from GiD 

131 

## Pore-pressure field 

132 

## Contact 

Ignasi de Pouplana – ipouplana@cimne.upc.edu - ignasi.de.pouplana@upc.edu Danilo Cavalcanti    - dborges@cimne.upc.edu 

—_ 

Tutorial for Kratos Poromechanics application 

