REM @ECHO OFF

REM Identification for arguments
REM basename               = %1
REM Project directory      = %2
REM Problem Type directory = %3
 
REM OutputFile: "%2\%1.info"
REM ErrorFile: "%2\%1.err"


del "%2\\%1.info"
del "%2\\%1.err"

move "%2\\%1.dat" "%2\\%1.mdpa"
move "%2\\%1-1.dat" "%2\\ProjectParameters.py"
copy "%3\\..\\..\\python_scripts\\thermo_mechanic_script.py" "%2\\"

REM WARNING: one should check the following paths before running this file

set PATH=C:\\KratosInstall;C:\\KratosInstall\\libs;%PATH%
"C:\\KratosInstall\\runkratos" "%2\\thermo_mechanic_script.py" > "%2\\%1.info" 2> "%2\\%1.err"
