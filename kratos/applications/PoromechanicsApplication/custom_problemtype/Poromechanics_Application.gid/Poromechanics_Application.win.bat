REM @ECHO OFF
REM Identification for arguments
REM basename                          = %1
REM Project directory                 = %2
REM Problem directory                 = %3
 
REM OutputFile: "%2\%1.info"
REM ErrorFile: "%2\%1.err"

del "%2\%1.info"
del "%2\%1.err"
del "%2\%1.mdpa"
del "%2\ProjectParameters.py"

move "%2\%1.dat" "%2\%1.mdpa"
move "%2\%1-1.dat" "%2\ProjectParameters.py"
copy "%3\\..\\..\\python_scripts\\main_script.py" "%2\\"

REM Updathe PATH
set PATH=C:\\KratosMultiphysics\\KratosMultiphysics\\libs;%PATH%

C:\\KratosMultiphysics\\runkratos.exe %2\\main_script.py > "%2\%1.info" 2> "%2\%1.err"
