FROM mcr.microsoft.com/windows/servercore:ltsc2019

# Restore the default Windows shell for correct batch processing below.
SHELL ["cmd", "/S", "/C"]

# Download the Build Tools bootstrapper.
RUN powershell.exe -Command \
	mkdir c:\TEMP; \
    wget https://aka.ms/vs/17/release/vs_buildtools.exe -OutFile c:\TEMP\vs_buildtools.exe

# Install Build Tools
RUN C:\TEMP\vs_buildtools.exe --quiet --wait --norestart --nocache \
    --add Microsoft.VisualStudio.Product.BuildTools \
    --add Microsoft.VisualStudio.Workload.VCTools \
    --add Microsoft.VisualStudio.Component.VC.CMake.Project \
    --add Microsoft.VisualStudio.Component.Windows10SDK.19041 \
    --add Microsoft.VisualStudio.Component.TestTools.BuildTools \
    --add Microsoft.VisualStudio.Component.VC.140 \
 || IF "%ERRORLEVEL%"=="3010" EXIT 0 ; \
 	$ErrorActionPreference = 'Stop'

#install chocolatey
RUN powershell.exe -Command \
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
    iex ((New-Object System.Net.WebClient).DownloadString('http://chocolatey.org/install.ps1'))

#install git
RUN powershell.exe -Command \
    choco install -y 7zip.install

#install git
RUN powershell.exe -Command \
    choco install -y git

#install cmake
RUN powershell.exe -Command \
    choco install -y cmake

#download and extract boost
RUN powershell.exe -Command \
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
    wget https://altushost-swe.dl.sourceforge.net/project/boost/boost/1.74.0/boost_1_74_0.tar.gz -OutFile c:\TEMP\boost.tar.gz; \
    mkdir c:\boost;

#Decompress and clean
RUN powershell.exe -Command \
    ls c:\TEMP; \
    7z x "c:\TEMP\boost.tar.gz" -o"c:\boost"; \
    rm "c:\TEMP\boost.tar.gz"

#Install python 3.8
RUN powershell.exe -Command \
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
    wget https://www.python.org/ftp/python/3.8.2/python-3.8.2-amd64.exe -OutFile c:\temp\python38.exe; \
    mkdir c:\python\38; \
    Start-Process c:\temp\python38.exe -ArgumentList '/quiet InstallAllUsers=1 PrependPath=0 TargetDir=c:\\python\\38' -Wait; \
    c:\python\38\python.exe -m pip install --upgrade pip; \
    c:\python\38\python.exe -m pip install --upgrade numpy sympy h5py

#set env variables
ENV CMAKE       "C:\\cmake\\cmake-3.14.1-win64-x64\\bin\\cmake.exe"
ENV SEVEN_ZIP   "C:\\7zip\\7z.exe"
ENV VCVARS      "C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\BuildTools\\VC\\Auxiliary\\Build\\vcvars64.bat\\"
ENV BOOST       "C:\\boost\\boost_1_71_0"
ENV PYTHONROOT  "C:\\python"

# Default to PowerShell if no other command specified.
CMD ["powershell.exe"]

