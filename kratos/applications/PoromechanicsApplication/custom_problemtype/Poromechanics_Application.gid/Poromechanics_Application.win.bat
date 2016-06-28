rem USED BY GiD
rem OutputFile: %1.log
rem ErrorFile: %1.err

rem Information
rem basename = %1
rem currentdirectory = %2
rem problemtypedirectory = %3
 
rem Remove old files
del %1.post.bin
del %1.post.res
del %1.post.msh
del %1.log
del %1.err

rem Update input files
move %1.dat %1.mdpa
move %1-1.dat ProjectParameters.json
copy %3\\..\\..\\python_scripts\\poromechanics_main.py %2\\

rem Setting paths. WARNING: one should check them before running this file
set PATH=C:\\KratosInstall;C:\\KratosInstall\\libs;%PATH%

rem Execute the program
C:\\KratosInstall\\runkratos poromechanics_main.py > %1.log 2> %1.err
