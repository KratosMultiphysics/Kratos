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

## Content
### [1. Python](https://github.com/enisalite/Hitchhiker-Guide-SWE/blob/main/1_Installation_Guides.md#1-python-1) 
- [1.1. Anaconda](https://github.com/enisalite/Hitchhiker-Guide-SWE/blob/main/1_Installation_Guides.md#11-anaconda) 
- [1.2. Visual Studio Code](https://github.com/enisalite/Hitchhiker-Guide-SWE/blob/main/1_Installation_Guides.md#12-visual-studio-code) 
- [1.3. Jupyter Notebook](https://github.com/enisalite/Hitchhiker-Guide-SWE/blob/main/1_Installation_Guides.md#13-jupyter-notebook) 
- [1.4. FFmpeg](https://github.com/enisalite/Hitchhiker-Guide-SWE/blob/main/1_Installation_Guides.md#14-ffmpeg) 
- [1.5. Useful links](https://github.com/enisalite/Hitchhiker-Guide-SWE/blob/main/1_Installation_Guides.md#15-more-useful-links-for-python-optional) 

### [2. GiD & Kratos](https://github.com/enisalite/Hitchhiker-Guide-SWE/blob/main/1_Installation_Guides.md#2-installation-guide-for-gid--kratos) 
- [2.1. GiD](https://github.com/enisalite/Hitchhiker-Guide-SWE/blob/main/1_Installation_Guides.md#21-download-and-install-gid) 
- [2.2. Kratos](https://github.com/enisalite/Hitchhiker-Guide-SWE/blob/main/1_Installation_Guides.md#22-download-and-install-kratos) 

### [3. Paraview](https://github.com/enisalite/Hitchhiker-Guide-SWE/blob/main/1_Installation_Guides.md#3-paraview-1)

- [3.1. Download and Install Paraview](https://github.com/enisalite/Hitchhiker-Guide-SWE/blob/main/1_Installation_Guides.md#31-download-and-install-paraview)

### [4. Remote Computing](https://github.com/enisalite/Hitchhiker-Guide-SWE/blob/main/1_Installation_Guides.md#4-remote-computing-1) 
- [4.1. VPN](https://github.com/enisalite/Hitchhiker-Guide-SWE/blob/main/1_Installation_Guides.md#41-vpn) 
- [4.2. Secure Shell (SSH) Client](https://github.com/enisalite/Hitchhiker-Guide-SWE/blob/main/1_Installation_Guides.md#42-secure-shell-ssh-client) 

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

Anaconda is available under https://www.anaconda.com/distribution/. It is available for Windows, Linux, Mac OS). You can choose the Operating System of your choice by clicking on its respective symbol.

Choose the Python 3.10 version with 64-bit, if possible. After the download is finished, run the .exe file and go through each step of the setup. Typically, your install location should be in your "Users" file in the local disk. 

Example: C:\Users\maxmustermann\anaconda3

___
### 1.2. Visual Studio Code
Visual Studio Code is a helpful IDE and code editor. It comes with built-in support for JavaScript, TypeScript and Node.js and has a rich ecosystem of extensions for other languages and runtimes (such as C++, C#, Java, Python, PHP, Go, .NET). 

#### 1. Download and Install Visual Studio Code
Visual Studio Code should come together with the Anaconda package. You can also download it separately under https://code.visualstudio.com/Download for the operating system of your choice.

The videos under https://code.visualstudio.com/docs/getstarted/introvideos are helpful to get started using Visual Studio Code.

#### 2. Add Python Path to Visual Studio Code
Make sure to **add Python to your path** if not already available in visual studio code integrated terminal. For Windows, the process will be as the following:

"**Start -> Edit the system environment variables**". Click on **"Environment Variables"**. Under **"System variables"** find **"Path"** and click on **Edit**. Press **"New"** to add the python location to the path. In case you installed python via the anaconda package, go to the install location of anaconda (typically in the users file) and copy its path. 

Also add the anaconda script path, which has the same path as anaconda, but with \Scripts added to it. For example:

- Anaconda Path: C:\Users\maxmustermann\anaconda3
- Script Path: C:\Users\maxmustermann\anaconda3\Scripts

#### 3. Visual Studio Code Interface
The window of Visual Studio Code is split into:
- File explorer (*left*): Where all the opened folders and their respective files are located. It is located on the left of the window. 

- Code editor (*center*): Located in the center of the window. The code of your file will be displayed here. In case you install extensions, you can also open different files, such as pdf, excel, word, etc.

- Terminal (*below code editor*): Located in the bottom part of the window, the terminal (also "command prompt" on windows), is used to run your code. 
  In case the terminal is not displayed, you can go to View -> Terminal to open it.


#### 4. Useful tips
To open the terminal, you can go to View -> Terminal command.

To open a File/Folder, go to **File -> Open File / Open Folder**. This will add it to the File Explorer, while removing the old folder from your File Explorer. In order to add a folder you should go to **File -> Add Folder to Workspace"**.

You can run the python scripts directly by **right-click on script -> "Run Python File in Terminal"** or type **"python.exe file_name.py"** in the terminal.

Save the updated scripts by clicking on script and then going to **"File -> Save"**. You can also save all files by going to **"File -> Save all"**

#### 5. Useful Extensions
In order to unlock some functions on Visual Studio Code and make the user experience easier, you can install extensions. 

To search for extensions, you can either click on the **Extensions** button, or **(Ctrl + Shift + X)**.

There are plenty of them in the Extensions: Marketplace. Sometimes, when you have a file which you cannot open due to the missing extension, Visual Studio Code will recommend an extension which you can install in order to open the file. 

Some important extensions for the course are the following:

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

After installing Jupyter Notebook, you can open it either from the Anaconda Navigator, or by using the shortcut **"Start -> Jupyter Notebook"**. In Linux you need to navigate to the required folder in the terminal window and then type **"jupyter notebook"**.

#### 2. Adding Jupyter shortcut to the desired directory
In order to add the Jupyter shortcut to the desired directory, you can open the file location **"Start -> right-click Jupyter Notebook -> Open file location"**. When in the file location, go to **"right-click Jupyter Notebook Shortcut -> Properties"**. Under the properties you can go to "Shortcut -> Target". Replace **%USERPROFILE%** with the desired directory (C:\..). Now you have a shortcut of Jupyter Notebook in your desired directory.

#### 3. Jupyter Notebook Interface
The Jupyter Notebook Interface consists of the **Code Cell**, where the code is written and will run, and the **Markdown Cell** in which comments and explanations are written in markdown, which is a simple markup language. 

You can run the full script by going to **"Kernel -> Restart & Run All"**. In case you want to run only one block of script, left-click to select the script, then click on **"Run"**. 

#### 4. Keyboard Shortcuts for Jupyter Notebook
For an efficient way of working, some keyboard shortcuts are presented:
- Toggle between edit and command mode by selecting block, then using Esc and Enter, respectively.

Once in command mode, you can (**case sensitive**):
- Scroll up and down your cells with your Up and Down keys.
- Press A or B to insert a new cell above or below the active cell.
- M will transform the active cell to a Markdown cell.
- Y will set the active cell to a code cell.
- D + D (D twice) will delete the active cell.
- Z will undo cell deletion.
- Hold Shift and press Up or Down to select multiple cells at once.
- With multiple cells selected, Shift + M will merge your selection.
- Ctrl + Shift + -, in edit mode, will split the active cell at the cursor.
- You can also click and Shift + Click in the margin to the left of your cells to select them.

#### 5. Using Jupyter in the browser

Under **"https://jupyter.org/try -> JupyterLab"** you can open Jupyter in the browser, without having to install the application. You can upload your Jupyter Notebooks (.ipynb files) and then view, edit and run your scripts via the browser. 

___
### 1.4. FFMPEG 
FFmpeg is a multimedia framework that can be used for various operations in numerous media files. FFMPEG is necessary for the Beam Model in the Project Work when using ParOptBeam as well as for some of the Jupyter Notebooks offered during the course.

#### Download and Install FFMPEG
You can download FFMPEG under https://ffmpeg.org/download.html#build-windows. 

#### 1. For Windows:

Hover your mouse on the Windows symbol. Click on **"Windows builds from gyan.dev"**.

When located in the new website, select **"ffmpeg-git-full.7z"**. Make sure to have a program to extract the contents of ffmpeg. After extracting the content of ffmpeg, rename it to FFmpeg. Save ffmpeg in a directory of your choice, for example under local disk C:\ .

Next step is to add ffmpeg to your path. Go to Start -> Edit the system environment variables. Click on "Environment Variables". Under "User variables for ..." find "Path" and click on Edit. Add a new path which would be the path in which your FFmpeg is plus \bin. For example, if your FFmpeg folder is saved as FFMPEG under local disk C:\ the path would be **C:\FFMPEG\bin**.

#### 2. For Linux:
In ubuntu, use $ sudo apt-get install ffmpeg.

#### 3. For Mac OS:
Download and install FFmpeg from (https://www.ffmpeg.org/download.html). 

You can also visit this website, which shows the installation steps of FFmpeg: https://www.wikihow.com/Install-FFmpeg-on-Windows.

___
### 1.5. More useful links for python (optional)

#### 1. Windows
- WinPython (https://winpython.github.io/ ) comes together with Spyder as IDE.
- Edit Python files (e.g. helloWorld.py) with a text editor (e.g Notepad++ (https://notepad-plus-plus.org/).
- PyCharm (https://www.jetbrains.com/pycharm/), which is another IDE, as Spyder or VS Code.
- Eclipse (https://www.eclipse.org/ ) together with the Python extension (http://www.pydev.org/).

#### 2. Mac OS
- Get the free version of Canopy from https://store.enthought.com/downloads/#default.

#### 3. Linux
The sources for the installation of Python can be found here:
 - Python: https://www.python.org/downloads/source/.
 - Sci-Py packages: https://www.scipy.org/install.html.

___
## 2. Installation Guide for GiD & Kratos
GiD and Kratos will be used during the project work for the CFD simulation of the structure. GiD supports the geometrical modelling of the structure, meshing, boundary conditions, system properties, etc. Kratos is used as a solver, from which we receive our results. It is a multi-physics simulation code. Aside from CFD, GiD + Kratos might also be used for other purposes, such as Structural/Dynamic Analysis, Fluid-Structure-Interaction etc.

___
### 2.1. Download and Install GiD

To download GiD, visit https://www.gidsimulation.com/gid-for-science/downloads/. Choose your operating system **(Windows, Linux, Mac OS)**. For Mac, only GiD is available, without the precompiled Kratos problemtype. 

Download the recommended version during the course. The scripts of the course are frequently tested with the newest GiD versions. In this case you are recommended to use GiD 16.1.3.

**For Windows** -> Download the file and go through the installation. 

**For Linux** -> Install using the graphical way or through console:
- Enable executable: chmod +x gid16.0.1-linux-x64-Install, then double click to use the graphical install.

**or**
- ./gid16.0.1-linux-x64-Install --mode console.

In order to use GiD, you will need to activate a professional license, as the free license is limited. You will be provided with access to this license during the course.

After registering, open GiD and go to **" Help -> Register GiD -> Named user -> sign in"**. To sign in, put your TUM email adress and password you used during registration.

___
### 2.2. Download and Install Kratos

#### 1. Option: Install Kratos from GiD
You can install Kratos directly from GiD, by going to **"Data -> Problemtype -> Internet Retrieve"**. A window will then open. In the offered modules, select **Kratos (9.2.2) -> Retrieve Module**. 

#### 2. Option: Install Kratos from GiD archive files
In a more manual way, you can also install Kratos from GiD archive files.

Visit ftp://www.gidhome.com/pub/gidmodules/ *or* www.gidhome.com. Go to Download section -> Other downloads -> "Download a GiD Problem Type". Download the latest version of Kratos (9.2.2) for your respective operating system. The x64 version is recommended as further exercises have been created and tested using this one.

After downloading the folder, unzip it and move the "kratos.gid" folder into the “problemtypes“ directory in the installation directory of GiD (overwriting the existing folder):

e.g. move to “C:\Program Files\GiD\GiD 16.1.3\problemtypes“.

The next time GiD is started, the kratos problemtype can be found under *Data -> Problem type -> kratos*.

___
## 3. Paraview
Paraview is an open-source, visualization programm which will be used for the postprocessing of the CFD-Simulation during project work.

### 3.1. Download and Install Paraview

It is available to download under https://www.paraview.org/download/. Choose your operating system (**Windows, Linux, Mac OS**) and download the .exe installable file for **Paraview version 5.11**. Paraview can then be simply installed from the .exe file.

To open paraview on Windows, a shortcut will be created in the Start Menu. In Linux, open the terminal and write "paraview".

___
## 4. Remote Computing
Some additional tools, which are presented here, are necessary to execute the CFD simulation remotely from the cluster of the Chair. 

___
### 4.1. VPN
It is important to establish a vpn connection to the LRZ server of TUM. To connect to this server, the **eduVPN** client is recommended. Under *https://www.eduvpn.org/client-apps/* you can download the client depending on your operating system (Windows, macOS and Linux all available). Download the **.exe** file, run it and follow the installation steps. After opening eduVPN, you can log in with your TUM credentials in order to establish a connection.

___
### 4.2. Secure Shell (SSH) Client 
After being logged in to the TUM LRZ server, it is important to establish an SSH  connection  with the Chair's cluster. A recommended SSH client is  **Git Bash**, under *https://gitforwindows.org/*. The way to establish connection and commands in Git are more similar to Linux. Download the **.exe** file, run it and follow the installation steps. 

Other popular SSH clients are PuTTY, Open SSH, VSCode .etc.

The other steps regarding Remote Computing are explained in its seperate part of the Hitchhiker Guide.
