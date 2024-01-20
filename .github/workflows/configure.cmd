call "%ProgramFiles%\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" x64 || goto :error

set CC=cl.exe
set CXX=cl.exe

set KRATOS_SOURCE=%cd%
set KRATOS_BUILD=%cd%\build
set KRATOS_APP_DIR=applications

REM Set applications to compile .. see "ci_apps_windows.json"
set KRATOS_APPLICATIONS_FILE="%KRATOS_SOURCE%\ci_compiled_apps.txt"

del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\cmake_install.cmake"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeCache.txt"
del /F /Q "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%\CMakeFiles"

cmake                                                 ^
  -G"Visual Studio 17 2022"                           ^
  -H"%KRATOS_SOURCE%"                                 ^
  -B"%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%"              ^
  -DBOOST_ROOT="%TEMP%\boost"                         ^
  -DCMAKE_CXX_FLAGS="/Od /we4661 /we4804 /WX /wd4996" ^
  -DFORCE_LOCAL_ZLIB_COMPILATION=ON                   ^
  -DCMAKE_UNITY_BUILD=ON                                    || goto :error

cmake --build "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%" --target all_build -- /property:configuration=%KRATOS_BUILD_TYPE% /p:Platform=x64 || goto :error
cmake --build "%KRATOS_BUILD%\%KRATOS_BUILD_TYPE%" --target install -- /property:configuration=%KRATOS_BUILD_TYPE% /p:Platform=x64 || goto :error

goto :EOF

:error
echo Failed with error #%errorlevel%.
exit /b %errorlevel%
