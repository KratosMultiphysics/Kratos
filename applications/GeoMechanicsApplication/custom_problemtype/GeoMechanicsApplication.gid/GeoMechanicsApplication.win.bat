rem USED BY GiD
rem OutputFile: %1.info
rem ErrorFile: %1.err

rem Information
rem basename = %1
rem currentdirectory = %2
rem problemtypedirectory = %3

rem Setting paths. WARNING: one should check them before running this file
set PATH=D:\\src\\Kratos;D:\\src\\Kratos\\libs;%PATH%

rem Execute the program
runkratos AddConditionsToExcavation.py
runkratos split_linear_beam_element_3n.py

runkratos MainKratos.py > %1.info 2> %1.err
