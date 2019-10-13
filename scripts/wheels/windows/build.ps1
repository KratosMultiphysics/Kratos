$pythons = "37","36","35"

foreach ($python in $pythons){
    echo "Begining build for python $($python)"

    #env cleanup
    mkdir c:\wheel
    cd c:\kratos\kratos
    git clean -ffxd
    $env:hash=$(git show -s --format=%h) #used in version number
    cd cmake_build
    cp c:\kratos\Kratos\scripts\wheels\windows\configure.bat .\configure.bat


    $pythonPath = "$($env:python)\$($python)\python.exe"
    cmd.exe /c "call configure.bat $($pythonPath)"
    MSBuild.exe /m INSTALL.vcxproj /p:Configuration=Custom /p:Platform="x64"

    echo "Finished build"
    echo "Begining wheel construction for python $($python)"

    cd c:\kratos\kratos
    cp -r KratosMultiphysics c:\wheel
    mkdir c:\wheel\KratosMultiphysics\.libs
    cp -r libs\* c:\wheel\KratosMultiphysics\.libs
    cp scripts\wheels\windows\__init__.py c:\wheel\KratosMultiphysics\__init__.py
    cp scripts\wheels\windows\setup.py c:\wheel\setup.py
    cp scripts\wheels\windows\README.md c:\wheel\README.md
    cd c:\wheel

    & $pythonPath setup.py bdist_wheel
    cp c:\wheel\dist\* c:\out\
    cd c:\
    rm -r c:\wheel

    echo "Finished wheel construction for python $($python)"
}
