---
title: Installation Guides
keywords: 
tags: [Installation_Guides.md]
sidebar: kratos_workshops
summary: 
---
# Installation Guides
The following part of the guide will show you all the necessary tools to be prepared for the Structural Wind Engineering (SWE) course and the project work. The guide contains the necessary steps for the installation of the tools and a couple of tips to get started. It is advised to go through all the installation steps, without neglecting any parts, in order to avoid potential errors when using the tools during the course.

It is important to note that, given the annual organization of the Structural Wind Engineering (SWE) course, the mandatory program versions may be subject to modification. As such, participants will be clarified regarding the specific versions required during the course.

___
## 1. Python
Most of the tools necessary for the course and project require to use python. The version necessary for the course is Python 3.10, together with the following modules:
- Numpy
- Sympy
- Matplotlib
- SciPy
- etc.

An Integrated Development Environment (**IDE**) is highly recommended to carry on the Python processes. The project course will mostly focus on using **Visual Studio Code**, however, other IDE-s such as **Spyder** or **PyCharm** are also available.

Another python-based program you will need is **Jupyter Notebook**.
During the semester you will be supplied with Jupyter Notebook Scripts, which are helpful to consolidate your knowledge but also to use as helpful tools for the project work. 

___
### 1.1. Anaconda 

An efficient way to install Python together with all its modules is installing the latest **Anaconda package** (Anaconda 22.10), as it comes with some necessary modules and a web application of **Jupyter Notebook**. Visual Studio Code can also be installed together with Anaconda. 

Anaconda is available for [download under this page](https://www.anaconda.com/distribution/){:target="_blank"}. It is available for Windows, Linux, Mac OS. You can choose the Operating System of your choice by clicking on its respective symbol.

Choose the Python 3.10 version with 64-bit, if possible. After the download is finished, run the .exe file and go through each step of the setup. Typically, your install location should be in your "Users" file in the local disk. 

Example: `C:\Users\maxmustermann\anaconda3`

___
### 1.2. Visual Studio Code
Visual Studio Code is a helpful IDE and code editor. It comes with built-in support for JavaScript, TypeScript and Node.js and has a rich ecosystem of extensions for other languages and runtimes (such as C++, C#, Java, Python, PHP, Go, .NET). 

#### 1. Download and Install Visual Studio Code
Visual Studio Code should come together with the Anaconda package. You can also [download it separately here](https://code.visualstudio.com/Download){:target="_blank"} for the operating system of your choice.

Here are some [helpful videos to get started](https://code.visualstudio.com/docs/getstarted/introvideos){:target="_blank"}  using Visual Studio Code.

#### 2. Add Python Path to Visual Studio Code
Make sure to **add Python to your path** if not already available in visual studio code integrated terminal. For Windows, the process will be as the following:

**Start menu &rarr; search for "Edit the system environment variables"**. Click on **"Environment Variables"**. Under **"System variables"** find **"Path"** and click on **Edit**. Press **"New"** to add the python location to the path. In case you installed python via the anaconda package, go to the install location of anaconda (in the users file by default) and copy its path. 

Also add the anaconda script path, which has the same path as anaconda, but with \Scripts added to it. For example:

- Anaconda Path: ```C:\Users\maxmustermann\anaconda3```
- Script Path: ```C:\Users\maxmustermann\anaconda3\Scripts```

#### 3. Visual Studio Code Interface
The window of Visual Studio Code is split into:
- File explorer (*top side bar*): Where all the opened folders and their respective files are located. It is located on the left of the window. 

- Code editor (*center*): Located in the center of the window. The code of your file will be displayed here. In case you install extensions, you can also open different files, such as pdf, excel, word, etc.

- Terminal (*below code editor*): Located in the bottom part of the window, the terminal (also "command prompt" on windows), is used to run your code. 
  In case the terminal is not displayed, you can go to View &rarr; Terminal to open it.


#### 4. Useful tips
To open the terminal, you can go to View &rarr; Terminal command. You can have multiple terminal simultaneously. To add more terminals, go to **Menu bar &rarr; Terminal &rarr; New Terminal** or press **(Ctrl + Shift + `)** on your keyboard.

To open a File/Folder, go to **File &rarr; Open File / Open Folder**. This will add it to the File Explorer, while removing the old folder from your File Explorer. In order to add a folder you should go to **File &rarr; Add Folder to Workspace"**.

You can run the python scripts directly by **right-click on script &rarr; "Run Python File in Terminal"** or type **"python.exe file_name.py"** in the terminal.

Save the updated scripts by clicking on script and then going to **"File &rarr; Save"**. You can also save all files by going to **"File &rarr; Save all"**

#### 5. Useful Extensions
In order to unlock some functions on Visual Studio Code and make the user experience easier, you can install extensions. 

To search for extensions, you can either click on the **Extensions** button on the left sidebar, or **(Ctrl + Shift + X)**.

There are plenty of them in the Extensions: Marketplace. Sometimes, when you have a file which you cannot open due to the missing extension, Visual Studio Code will recommend an extension which you can install in order to open the file. 

Some useful extensions for the course are the following:

For Python:
- Python
- Python Extended
- Python for VSCode
- Python-autopep8

For JSON:
- JSON Tools
- Prettify JSON

For Remote Computing:
- Remote SSH

Jupyter Notebook:
- Jupyter

___
### 1.3. Jupyter Notebook
Jupyter Notebook is a web-based interactive computing platform. Its biggest advantage is its interactive block output, making it easy for the user to understand certain processes in the code easier. 

#### 1. Getting started with Jupyter Notebook
Jupyter Notebook is installed together with the anaconda package. In case it is not installed due to reasons, you can write and run **"pip install notebook"** in the terminal. 

After installing Jupyter Notebook, you can open it either from the Anaconda Navigator, or directly from the start menu **"Start &rarr; Jupyter Notebook"**. In Linux you need to navigate to the required folder in the terminal window and then type **"jupyter notebook"**. This is also applicable in windows. The drive of which the jupyter notebook will operate depends on the directory it is launched from the terminal. 

#### 2. Adding Jupyter shortcut to the desired directory
In order to add the Jupyter shortcut to the desired directory, you can open the file location **"Start &rarr; right-click Jupyter Notebook &rarr; Open file location"**. When in the file location, go to **"right-click Jupyter Notebook Shortcut &rarr; Properties"**. Under the properties you can go **to "Shortcut &rarr; Target"**. Replace **%USERPROFILE%** with the desired directory (```D:\..```). Now you have a shortcut of Jupyter Notebook in your desired directory.

#### 3. Jupyter Notebook Interface
The Jupyter Notebook Interface consists of the **Code Cell**, where the code is written and will run, and the **Markdown Cell** in which comments and explanations are written in markdown, which is a simple markup language. 

You can run the full script by going to **"Kernel &rarr; Restart & Run All"**. In case you want to run only one block of script, left-click to select the script, then click on **"Run"**. 

#### 4. Keyboard Shortcuts for Jupyter Notebook
For an efficient way of working, some keyboard shortcuts are presented:

To enter the **COMMAND** mode press **ESC** or click anywhere outside the cell. You will see grey border around the cell with blue left margin. When you are in Command mode, you can edit your notebook but you can't type in the cells.

Once in command mode, you can (**case sensitive**):

| Command | Function |
| ------- | ---------|
| **&uarr;**/**&darr;** | Scroll up or down your cells | 
| **A** / **B** | Insert a new cell above or below the active cell | 
| **M** | Transform the active cell to a Markdown cell | 
| **Y** | Set the active cell to a code cell | 
| **D + D** | Delete the active cell | 
| **Z** | Undo cell deletion | 
| **Shift + &uarr;**/**&darr;** | Select multiple cells at once | 
| **Shift + M** | Merge the selected cells | 
| **Ctrl + Shift + -** | In edit mode, split the active cell at the cursor | 

You can also click and Shift + Click in the margin to the left of your cells to select them.

For more keyboard shortcuts, go to **Help &rarr; Keyboard Shortcuts**

#### 5. Using Jupyter in the browser

Alternatively, you can [open Jupyter in the browser](https://jupyter.org/try-jupyter/lab?path=notebooks%2FIntro.ipynb){:target="_blank"} without having to install the application. Using the web app, you can upload your Jupyter Notebooks (.ipynb files) and then view, edit and run your scripts via the browser. 

___
### 1.4. FFMPEG 
FFmpeg is a multimedia framework that can be used for various operations in numerous media files. FFMPEG is necessary for the Beam Model in the Project Work when using ParOptBeam as well as for some of the Jupyter Notebooks offered during the course.

#### Download and Install FFMPEG
FFmpeg is available for [download here](https://www.ffmpeg.org/download.html){:target="_blank"}.

#### 1. For Windows:

Hover your mouse on the Windows symbol. Click on **"Windows builds from gyan.dev"**.

When located in the new website, select **"ffmpeg-git-full.7z"**. Make sure to have a program to extract the contents of ffmpeg. After extracting the content of ffmpeg, rename it to FFmpeg. Save ffmpeg in a directory of your choice, for example under local disk C:\ .

Next step is to add ffmpeg to your path. Go to **Start menu &rarr; search for "Edit the system environment variables"**. Click on **"Environment Variables"**. Under **"User variables for ..."** find **"Path"** and click on **Edit**. Add a new path which would be the path in which your FFmpeg is, followed by ```\bin```. For example, if your FFmpeg folder is saved as FFMPEG under local disk C:\ the path would be ```C:\FFMPEG\bin```.

You can also visit [this website](https://www.wikihow.com/Install-FFmpeg-on-Windows){:target="_blank"}, which shows the installation steps of FFmpeg.

#### 2. For Linux:
In ubuntu, start by updating the packages list:

```shell
$ sudo apt update
```

Next, install FFmpeg by typing the following command:

```shell
$ sudo apt-get install ffmpeg
```

Verify your installation using:

```shell
$ ffmpeg -version
```

#### 3. For Mac OS:
Download and install FFmpeg from [the website](https://www.ffmpeg.org/download.html){:target="_blank"} and follow the installation steps.

___
### 1.5. More useful links for python (optional)

#### 1. Windows
- [WinPython](https://winpython.github.io/){:target="_blank"} comes together with Spyder as IDE.
- [Notepad++](https://notepad-plus-plus.org/){:target="_blank"}, Lightweight text editor for editing text files.
- [PyCharm](https://www.jetbrains.com/pycharm/){:target="_blank"}, another IDE, just like Spyder or VS Code.
- [Eclipse](https://www.eclipse.org/){:target="_blank"} together with the [Python extension](http://www.pydev.org/){:target="_blank"}.

#### 2. Mac OS
- ~~Get the free version of [Canopy](https://store.enthought.com/downloads/#default){:target="_blank"} here.~~

#### 3. Linux
The sources for the installation of Python can be found here:
 - [Python](https://www.python.org/downloads/source/){:target="_blank"}
 - [Sci-Py packages](https://www.scipy.org/install.html){:target="_blank"}

___
## 2. Installation Guide for GiD & Kratos module
GiD and the Kratos module for GiD will be used during the project work for the CFD simulation of the structure. GiD supports the geometrical modelling of the structure, meshing, boundary conditions, system properties, etc. Kratos is used as a solver, from which we receive our results. It is a multi-physics simulation code. Aside from CFD, GiD + Kratos might also be used for other purposes, such as Structural/Dynamic Analysis, Fluid-Structure-Interaction etc. Kratos module in GiD add some features to support Kratos in GiD, which will generate the input for the solver.

___
### 2.1. Download and Install GiD

[Download GiD at their website](https://www.gidsimulation.com/gid-for-science/downloads/){:target="_blank"}. Choose your operating system **(Windows, Linux, Mac OS)**. For Mac, only GiD is available, without the precompiled Kratos problemtype. 

Download the recommended version during the course. The scripts of the course are frequently tested with the newest GiD versions. In this case you are recommended to use GiD 16.1.3.

**For Windows** &rarr; Download the file and go through the installation. 

**For Linux** &rarr; Install using the graphical way or through console:
- Enable executable: chmod +x gid16.0.1-linux-x64-Install, then double click to use the graphical install.

**or**
- ./gid16.0.1-linux-x64-Install --mode console.

In order to use GiD, you will need to activate a professional license, as the free license is limited. You will be provided with access to this license during the course.

After registering, open GiD and go to **" Help &rarr; Register GiD &rarr; Named user &rarr; sign in"**. To sign in, put your TUM email adress and password you used during registration.

___
### 2.2. Download and Install Kratos module

#### 1. Option: Install Kratos module from GiD
You can install Kratos module directly from GiD, by going to **"Menu bar &rarr; Data &rarr; Problemtype &rarr; Internet Retrieve"**. A window will then open. In the offered modules, select **Kratos (9.2.2) &rarr; Retrieve Module**. 

#### 2. Option: Install Kratos module from GiD archive files
You can also install [Kratos module](https://www.gidsimulation.com/downloads/kratos/){:target="_blank"} manually from [GiD archive files.](https://downloads.gidsimulation.com/#gidmodules/){:target="_blank"}. Download the latest version of Kratos (9.2.2) for your respective operating system. The x64 version is recommended as further exercises have been created and tested using this one. After downloading the folder, unzip it and move the **"kratos.gid"** folder into the **“problemtypes“** folder inside the installation directory of GiD (overwriting the existing folder). By default, the installation directory is:

```C:\Program Files\GiD\GiD 16.1.3\problemtypes```

The next time GiD is started, the kratos problemtype can be found under **Data &rarr; Problem type &rarr; kratos**.

___
## 3. Paraview
Paraview is an open-source, visualization programm which will be used for the postprocessing of the CFD-Simulation during project work.

### 3.1. Download and Install Paraview

It is available for download [under this link](https://www.paraview.org/download/){:target="_blank"}. Download the .exe installable file for **Paraview version 5.11** according to your operating system (**Windows, Linux, Mac OS**). Paraview can then be simply installed from the .exe file.

On Windows, you can open paraview from the Start Menu. In Linux, open the terminal and write **"paraview"**.

___
## 4. Remote Computing
Some additional tools, which are presented here, are required to [execute the CFD simulation remotely](Remote_Computing.html){:target="_blank"} from the cluster of the Chair. 

___
### 4.1. VPN
It is important to establish a vpn connection to the LRZ server of TUM. To connect to this server, the **eduVPN** client is recommended. [Download the client](https://www.eduvpn.org/client-apps/){:target="_blank"} according to your operating system (Windows, macOS and Linux all available). Download the **.exe** file, run it and follow the installation steps. After opening eduVPN, you can log in with your TUM credentials in order to establish a connection.

___
### 4.2. Secure Shell (SSH) Client 
After being logged in to the TUM LRZ server, it is important to establish an SSH  connection  with the Chair's cluster. A recommended SSH client is [Git Bash](https://gitforwindows.org/){:target="_blank"}. The way to establish connection and commands in Git are more similar to Linux. Download the **.exe** file, run it and follow the installation steps. 

Other popular SSH clients are:
- [PuTTY](https://www.putty.org/){:target="_blank"}
- [Open SSH](https://learn.microsoft.com/en-us/windows-server/administration/openssh/openssh_install_firstuse?tabs=gui){:target="_blank"} 
- [VSCode](https://code.visualstudio.com/docs/remote/ssh){:target="_blank"} 

The other steps regarding [Remote Computing](Remote_Computing.html){:target="_blank"} are explained in a seperate part of the Hitchhiker Guide.