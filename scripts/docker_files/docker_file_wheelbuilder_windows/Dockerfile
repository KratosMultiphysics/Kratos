FROM mcr.microsoft.com/windows/servercore:ltsc2019


# Restore the default Windows shell for correct batch processing below.
SHELL ["cmd", "/S", "/C"]

# Download the Build Tools bootstrapper.
RUN powershell.exe -Command \
	mkdir c:\TEMP; \
    wget https://aka.ms/vs/16/release/vs_buildtools.exe -OutFile c:\TEMP\vs_buildtools.exe

# Install Build Tools
RUN C:\TEMP\vs_buildtools.exe --quiet --wait --norestart --nocache \
    --add Microsoft.VisualStudio.Product.BuildTools \
    --add Microsoft.VisualStudio.Workload.VCTools \
    --add Microsoft.VisualStudio.Component.VC.CMake.Project \
    --add Microsoft.VisualStudio.Component.Windows10SDK.17763 \
    --add Microsoft.VisualStudio.Component.TestTools.BuildTools \
    --add Microsoft.VisualStudio.Component.VC.140 \
 || IF "%ERRORLEVEL%"=="3010" EXIT 0 ; \
 	$ErrorActionPreference = 'Stop'

#install chocolatey
ENV chocolateyVersion '1.4.0'
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
    wget https://archives.boost.io/release/1.74.0/source/boost_1_74_0.zip -OutFile c:\TEMP\boost.zip; \
    mkdir c:\boost; \
    7z x c:\TEMP\boost.zip -o"c:\boost"

#Install python 3.8
RUN powershell.exe -Command \
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
    wget https://www.python.org/ftp/python/3.8.2/python-3.8.2-amd64.exe -OutFile c:\temp\python38.exe; \
    mkdir c:\python\38; \
    Start-Process c:\temp\python38.exe -ArgumentList '/quiet InstallAllUsers=1 PrependPath=0 TargetDir=c:\\python\\38' -Wait; \
    c:\python\38\python.exe -m pip install --upgrade pip; \
    c:\python\38\python.exe -m pip install mypy; \
    c:\python\38\python.exe -m pip install --upgrade setuptools wheel

#Install python 3.9
RUN powershell.exe -Command \
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
    wget https://www.python.org/ftp/python/3.9.0/python-3.9.0-amd64.exe -OutFile c:\temp\python39.exe; \
    mkdir c:\python\39; \
    Start-Process c:\temp\python39.exe -ArgumentList '/quiet InstallAllUsers=1 PrependPath=0 TargetDir=c:\\python\\39' -Wait; \
    c:\python\39\python.exe -m pip install --upgrade pip; \
    c:\python\39\python.exe -m pip install mypy; \
    c:\python\39\python.exe -m pip install --upgrade setuptools wheel

#Install python 3.10
RUN powershell.exe -Command \
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
    wget https://www.python.org/ftp/python/3.10.0/python-3.10.0-amd64.exe -OutFile c:\temp\python310.exe; \
    mkdir c:\python\310; \
    Start-Process c:\temp\python310.exe -ArgumentList '/quiet InstallAllUsers=1 PrependPath=0 TargetDir=c:\\python\\310' -Wait; \
    c:\python\310\python.exe -m pip install --upgrade pip; \
    c:\python\310\python.exe -m pip install --upgrade setuptools wheel

#Install python 3.11
RUN powershell.exe -Command \
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
    wget https://www.python.org/ftp/python/3.11.0/python-3.11.0-amd64.exe -OutFile c:\temp\python311.exe; \
    mkdir c:\python\311; \
    Start-Process c:\temp\python311.exe -ArgumentList '/quiet InstallAllUsers=1 PrependPath=0 TargetDir=c:\\python\\311' -Wait; \
    c:\python\311\python.exe -m pip install --upgrade pip; \
    c:\python\311\python.exe -m pip install --upgrade setuptools wheel

# Download and extract hdf5
RUN powershell.exe -Command \
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
    wget https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_12_2/source/hdf5-1.12.2.zip -OutFile c:\TEMP\hdf5.zip; \
    mkdir c:\hdf5; \
    mkdir c:\hdf5\build; \
    mkdir c:\hdf5\bin; \
    mkdir c:\hdf5\source; \
    7z x c:\TEMP\hdf5.zip -o"c:\hdf5\source"

# Download and extract med
RUN powershell.exe -Command \
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
    wget 'https://www.salome-platform.org/?sdm_process_download=1"&"download_id=2540' -OutFile c:\TEMP\med.tar.bz; \
    mkdir c:\med\build; \
    mkdir c:\med\bin; \
    mkdir c:\med\source; \
    7z x c:\TEMP\med.tar.bz -o"c:\TEMP\med.tar"; \
    7z x c:\TEMP\med.tar -o"c:\med\source"; \
    dir "c:\med\source"

# Build Hdf5 lib
RUN powershell.exe -Command \
    cd 'c:/Program Files/CMake/bin/'; \
    ./cmake.exe -H"c:\hdf5\source\hdf5-1.12.2" -B"c:\hdf5\build" -DCMAKE_INSTALL_PREFIX="c:\hdf5\bin" -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=ON; \
    ./cmake.exe --build "c:\hdf5\build" --target install -- /property:configuration=Release

# Build med lib
Run dir c:\med\source
RUN powershell.exe -Command \
    cd 'c:/Program Files/CMake/bin/'; \
    ./cmake.exe -H"c:\med\source\med-5.0.0" -B"c:\med\build" -DCMAKE_INSTALL_PREFIX="c:\med\bin" -DMEDFILE_BUILD_TESTS=OFF -DMEDFILE_BUILD_SHARED_LIBS=ON -DHDF5_ROOT="c:\hdf5\bin"; \
    ./cmake.exe --build "c:\med\build" --target install -- /property:configuration=Release

# #Downloand blas/lapack
# RUN powershell.exe -Command \
# 	[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; \
# 	mkdir c:\libs; \
# 	wget https://icl.cs.utk.edu/lapack-for-windows/libraries/VisualStudio/3.7.0/Dynamic-MINGW/Win64/libblas.lib -OutFile c:\libs\libblas.lib; \
# 	wget https://icl.cs.utk.edu/lapack-for-windows/libraries/VisualStudio/3.7.0/Dynamic-MINGW/Win64/liblapack.lib -OutFile c:\libs\liblapack.lib

#set env variables
ENV CMAKE       "C:\\cmake\\cmake-3.14.1-win64-x64\\bin\\cmake.exe"
ENV SEVEN_ZIP   "C:\\7zip\\7z.exe"
ENV VCVARS      "C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\BuildTools\\VC\\Auxiliary\\Build\\vcvars64.bat\\"
ENV BOOST       "C:\\boost\\boost_1_71_0"
ENV PYTHONROOT  "C:\\python"

COPY start.ps1 c:\\scripts\\start.ps1

# Default to PowerShell if no other command specified.
ENTRYPOINT ["powershell.exe", "-NoLogo", "-ExecutionPolicy", "Bypass", "-Command", "c:\\scripts\\start.ps1"]

