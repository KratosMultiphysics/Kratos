@ECHO OFF

DEL %2\%1.info
DEL %2\%1.flavia.res
DEL %2\%1.flavia.dat
DEL %2\%1.err

DEL %2\%1.mdpa
DEL %2\%1_aux.win.bat

REN %2\%1.dat %2\%1.mdpa
DEL %2\%1-1.dat
REN %2\%1-2.dat %2\%1_aux.win.bat

CALL %1_aux.win.bat %1 %2

REM D:\kratos\kratos\sources\Release\Kratos %1 > %2\%1.info

