rem USED BY GiD
rem OutputFile: %1.info
rem ErrorFile: %1.err

rem Information
rem basename = %1
rem currentdirectory = %2
rem problemtypedirectory = %3

rem Setting paths. WARNING: one should check them before running this file
set PYTHONPATH=%PYTHONPATH%;C:\\Kratos\\bin\\Release
set PATH=%PATH%;C:\\Kratos\\bin\\Release\\libs

rem Execute the program
C:\\Kratos\\runkratos MainKratos.py > %1.info 2> %1.err
