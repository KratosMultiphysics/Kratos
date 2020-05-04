
# mmg - Surface and volume remeshers 
[![Build Status](https://ci.inria.fr/mmg/buildStatus/icon?job=Mmg-multiConf)](https://ci.inria.fr/mmg/job/Mmg-multiConf/)[![Quality Gate](https://sonarqube.inria.fr/sonarqube/api/badges/gate?key=cardamom%3Ammg%3Agithub%3Adevelop)](https://sonarqube.inria.fr/sonarqube/dashboard?id=cardamom%3Ammg%3Agithub%3Adevelop)

[![Vulnerabilities](https://sonarqube.inria.fr/sonarqube/api/badges/measure?key=cardamom%3Ammg%3Agithub%3Adevelop&metric=vulnerabilities)](https://sonarqube.inria.fr/sonarqube/component_measures?id=cardamom%3Ammg%3Agithub%3Adevelop&metric=vulnerabilities)[![Potential Bugs](https://sonarqube.inria.fr/sonarqube/api/badges/measure?key=cardamom%3Ammg%3Agithub%3Adevelop&metric=bugs)](https://sonarqube.inria.fr/sonarqube/component_measures?id=cardamom%3Ammg%3Agithub%3Adevelop&metric=bugs)[![Coverage](https://sonarqube.inria.fr/sonarqube/api/badges/measure?key=cardamom%3Ammg%3Agithub%3Adevelop&metric=coverage)](https://sonarqube.inria.fr/sonarqube/component_measures?id=cardamom%3Ammg%3Agithub%3Adevelop&metric=coverage)

[![New Vulnerabilities](https://sonarqube.inria.fr/sonarqube/api/badges/measure?key=cardamom%3Ammg%3Agithub%3Adevelop&metric=new_vulnerabilities)](https://sonarqube.inria.fr/sonarqube/component_measures?id=cardamom%3Ammg%3Agithub%3Adevelop&metric=new_vulnerabilities)[![New Potential Bugs](https://sonarqube.inria.fr/sonarqube/api/badges/measure?key=cardamom%3Ammg%3Agithub%3Adevelop&metric=new_bugs)](https://sonarqube.inria.fr/sonarqube/component_measures?id=cardamom%3Ammg%3Agithub%3Adevelop&metric=new_bugs)








mmg is an open source software for bidimensional and tridimensional surface and volume remeshing.

It provides 3 applications and 4 libraries:
  * the **mmg2d** application and library: mesh generation from a set of edges, adaptation and optimization of a bidimensionnal triangulation and isovalue discretization;
  * the **mmgs** application and library: adaptation and optimization of a surface triangulation and isovalue discretization;
  * the **mmg3d** application and library: adaptation and optimization of a tetrahedral mesh, isovalue discretization and lagrangian movement;
  * the **mmg** library gathering the **mmg2d**, **mmgs** and **mmg3d** libraries.

[//]: # ( comment )

## Get and compile the mmg project
### Needed tools
To get and build Mmg, you will need:
  * **Git**: to download the code you will have to use a git manager. You can install a git manager from the link below but there are many other git clients that you can use:
    * [Official Git client](https://git-scm.com/download) (command line program)
    * [GitKraken](https://www.gitkraken.com/)
    * [SourceTree](https://www.sourcetreeapp.com/)  

    Note that if you uses Microsoft Visual Studio (Windows OS), you can simply activate the Git Module of the application.

  * **CMake** : Mmg uses the CMake building system that can be downloaded on the
    following web page:
    [https://cmake.org/download/](https://cmake.org/download/). On Windows OS,
    once CMake is installed, please <span style="color:red"> do not forget to
    mark the option: 
    ```
    "Add CMake to the system PATH for all users"
    ```
    </span>  

### Mmg download and compilation
#### Unix-like OS (Linux, MacOS...)

  1. Get the repository:  
      ```Shell
      wget https://github.com/MmgTools/mmg/archive/master.zip
      ```
     or
      ```Shell
      git clone https://github.com/MmgTools/mmg.git
      ```

  The project sources are available under the **_src/_** directory, see:
   * **_src/mmg2d/_**   for files related to the mmg2d application;
   * **_src/mmgs/_**   for files related to the mmgs application;
   * **_src/mmg3d/_**  for files related to the mmg3d application;
   * **_src/common/_** for files related to the both.

  2. Fast compilation (build both **mmg2d**, **mmgs**, **mmg3d**, the mmg2d static library (**libmmg3d.a**), the mmgs static library (**libmmgs.a**), the mmg3d static library (**libmmg3d.a**) and the mmg static library (**libmmg.a**)):  
```Shell
      cd mmg  
      mkdir build  
      cd build  
      cmake ..  
      make  
      make install
```

  If the `make install` command fail, try to run the `sudo make install` command.
  If you don't have root access, please refers to the [Installation section](https://github.com/MmgTools/Mmg/wiki/Setup-guide#iii-installation) of the [setup guide](https://github.com/MmgTools/Mmg/wiki/Setup-guide#setup-guide).

  The **mmg2d**, **mmgs** and **mmg3d** applications are available under the `mmg2d_O3`, `mmgs_O3` and `mmg3d_O3` commands. 

Note that if you use some specific options and want to set it easily, you can use a shell script to execute the previous commands. An example is provided [here](https://github.com/MmgTools/mmg/wiki/Configure-script-for-CMake-(UNIX-like-OS)).

#### Windows OS
The following compilation can be performed in any modern version of *Windows*
(AKA 7, 8, 8.1 and 10). A basic knowledge of Windows is assumed (execute
commands in cmd, create directories, etc...).

##### Compile with VisualStudio

Universal windows platform development
  1. Get the **Visual Studio** software: it can be downloaded [here](https://www.visualstudio.com/downloads/);

  2. if not done during the previous step, download **C/C++** compilers: in the Visual Studio searching zone, search **C compiler** and install the **"Visual C++ compilers and libraries"** (individual componant) and the MSBuild componant;

  3. in the Visual Studio searching zone, search the **git** word and select the installation of the **"GitHub extension for VisualStudio"**;

  4. stay in VisualStudio and clone the Mmg repository from the following url: https://github.com/MmgTools/mmg.git;

  5. Use **CMake** to configure and generate your project. It can be done either with the graphic mode of CMake (you have to select the "VisualStudio" generator) or with a command line. In this case, it is highly recommended to specify that you intent to build a VisualStudio project.
       For example, if you are using VisualStudio 2017:
  ```Shell
    cmake -G "Visual Studio 15 2017 Win64" ^
    configure
  ```
  
  Note that you can use a script to make this step easier (an example of script is provided [here](https://github.com/MmgTools/mmg/wiki/Configure-script-for-CMake-(Windows-OS))).  
  
   Once the configuration script has finished without errors a `mmg.sln` file will be generated in the cmake_build directory.

  6. Double click this file and the visual studio project will open. Then choose the project configuration (Release, Debug...).  
  Please, make sure that the project is set to Win32 or x64 and change it if is not. 
  Finally, in order to compile Mmg, right click the `INSTALL` project and select the option `BUILD`.

##### Compile with MinGW

  1. Get a **C Compiler**: 
      * **MinGW** can be downloaded [here](http://mingw.org/). We recommand to install the *mingw-developer-tools*, *mingw32-base*, *mingw32-gcc-fortran*, *mingw32-gcc-g++* and *msys-base* packages;
      * Edit the environment variables and add MinGW in your **PATH** variable. It can be done in the **advanced system settings** panel. (note that you must modify the **PATH** variable, not the **Path** one);
      * **MinGW** binaries are probably in **C:\MinGW\bin**
      * the MinGW terminal is in **C:\MinGW\msys\1.0\msys**

  2. Clone the Mmg repository from the following url: https://github.com/MmgTools/mmg.git;

  3. Quit and restart the *CMake* application to take the PATH modification into account then use CMake to configure and generate your project (select the MinGW Makefiles generator of CMake). If you have installed the scotch libraries, you will need to set explicitely the libraries paths;
  4. Build the Mmg applications: in the minGW prompt (**C:\MinGW\msys\1.0\msys**) run:
```Shell
       mingw32-make
```

Again, if you use some specific options and want to make the CMake configuration step easier, you can use a batch script. An example of script is provided [here](https://github.com/MmgTools/mmg/wiki/Configure-script-for-CMake-(Windows-OS)).

## Documentation
### Project's web page
Project's actualities and software tutorials can be found on the [mmgtools](http://www.mmgtools.org) web page.

### **Mmg**'s forum
Share your comments and issues with other members of the Mmg community on the [Mmg forum](https://forum.mmgtools.org/).

### GitHub's Wiki
More detailed informations about the compilation and configuration of the mmg's applications are available on the project [wiki](https://github.com/MmgTools/mmg/wiki).

### Man-pages
Man pages are available inside the **_doc/man_** directory:
  * To see the **mmg2d** man page, just tap `man ./doc/man/mmg2d.1.gz`
  * To see the **mmgs** man page, just tap `man ./doc/man/mmgs.1.gz`
  * To see the **mmg3d** man page, just tap `man ./doc/man/mmg3d.1.gz`

### Code documentation
Run the `make doc` command to build the Doxygen documentation.
  * To see the **mmg2d** documentation, open up the **_mmg/doc/mmg2d/html/index.html_** file;
  * To see the **mmgs** documentation, open up the **_mmg/doc/mmgs/html/index.html_** file;
  * To see the **mmg3d** documentation, open up the **_mmg/doc/mmg3d/html/index.html_** file.

## Platforms
The **mmg** applications are validated on OS X and on most of the Linux platforms. 

## Contributing
Your contributions to the **mmg** project are welcomed. You can help us to improve
our code by many means:
  * pull requests: please follow the [wiki's guideline](https://github.com/MmgTools/Mmg/wiki/Developers-wiki#pull-requests);
  * feature requests: please use the [Mmg forum](https://forum.mmgtools.org/);
  * bug reports: please use the [GitHub issue tracker](https://github.com/MmgTools/mmg/issues/new);

## About the team
mmg's current developers and maintainers are Charles Dapogny, Cécile Dobrzynski, Pascal Frey and Algiane Froehly.

Contact: contact@mmgtools.org

## License and copyright
Code is under the [terms of the GNU Lesser General Public License](https://raw.githubusercontent.com/MmgTools/mmg/master/LICENSE).

Copyright © Bx INP/Inria/UBordeaux/UPMC, 2004- .

## Reference
[Three-dimensional adaptive domain remeshing, implicit domain meshing, and applications to free and moving boundary problems - _C. Dapogny, C. Dobrzynski and P. Frey_ - April 1, 2014 - _JCP_](http://www.sciencedirect.com/science/article/pii/S0021999114000266)
