param([System.String]$cotire="OFF")

$pythons = "39"
$env:kratos_version = "8.2.0"

$kratosRoot = "c:\kratos\kratos"
$env:kratos_root = $kratosRoot
$wheelRoot = "c:\wheel"
$wheelOutDir = "c:\out"

function exec_build($python, $pythonPath) {
    cmd.exe /c "call configure.bat $($pythonPath) $($kratosRoot) OFF"
    cmake --build "$($kratosRoot)/build/Release" --target install -- /property:configuration=Release /p:Platform=x64
}

function exec_build_cotire($python, $pythonPath) {
    cmd.exe /c "call configure.bat $($pythonPath) $($kratosRoot) ON"
    cmake --build "$($kratosRoot)/build/Release" --target install -- /property:configuration=Release /p:Platform=x64
}

function build ($python, $pythonPath) {
    cd $kratosRoot
    cp "$($kratosRoot)\scripts\wheels\windows\configure.bat" .\configure.bat

    if($cotire -eq "ON"){
        exec_build_cotire $python $pythonPath
    }else {
        exec_build $python $pythonPath
    }
}

function  setup_wheel_dir {
    cd $kratosRoot
    mkdir c:\wheel
    cp scripts\wheels\setup.py c:\wheel\setup.py
    mkdir c:\wheel\KratosMultiphysics
    mkdir c:\wheel\KratosMultiphysics\.libs
}

function create_core_wheel ($pythonPath) {
    setup_wheel_dir
    cd $kratosRoot
    cp bin\release\KratosMultiphysics\* "$($wheelRoot)\KratosMultiphysics"
    cp scripts\wheels\windows\KratosMultiphysics.json "$($wheelRoot)\wheel.json"


    cp scripts\wheels\__init__.py "$($wheelRoot)\KratosMultiphysics\__init__.py"
    cd $wheelRoot
    & $pythonPath setup.py bdist_wheel
    cp "$($wheelRoot)\dist\*" $wheelOutDir
    cd c:\
    rm -r $wheelRoot
}

function create_application_wheel ($pythonPath, $app) {
    setup_wheel_dir
    cp "$($kratosRoot)\scripts\wheels\windows\applications\$($app)" c:\wheel\wheel.json
    cd $wheelRoot
    & $pythonPath setup.py bdist_wheel #pythonpath
    cp "$($wheelRoot)\dist\*" $wheelOutDir
    cd c:\
    rm -r $wheelRoot
}

foreach ($python in $pythons){
    Write-Host "Begining build for python $($python)"
    $env:python = $python

    cd $kratosRoot
    git clean -ffxd
    $pythonPath = "$($env:pythonRoot)\$($python)\python.exe"

    build $python $pythonPath

    Write-Host "Finished build"
    Write-Host "Begining wheel construction for python $($python)"

    create_core_wheel $pythonPath

    $applications = Get-ChildItem "$($kratosRoot)\scripts\wheels\windows\applications"

    foreach($app in $applications) {
        create_application_wheel $pythonPath $app
    }

    Write-Host "Finished wheel construction for python $($python)"
}
