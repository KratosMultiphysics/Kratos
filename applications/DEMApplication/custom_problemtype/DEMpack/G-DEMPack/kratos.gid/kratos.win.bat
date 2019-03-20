REM @ECHO OFF
REM Identification for arguments
REM basename                          = %1
REM Project directory                 = %2
REM Problem directory                 = %3
 
REM OutputFile: "%2\%1.info"
REM ErrorFile: "%2\%1.err"
 
DEL "%2\%1.info"
DEL "%2\%1.err"

REM Updathe PATH
set PATH=%3\\kratos;%3\\kratos\\libs;%3\\kratos\\bin;%PATH%
set PYTHONPATH=%3\\kratos

REM Set the number of threads for OpenMP
REM export OMP_NUM_THREADS=%5
set OMP_NUM_THREADS=%5

REM Run Python using the script KratosOpenMP.py
"%3\\kratos\\runkratos" KratosDEM.py > "%2\\%1.info" 2> "%2\\%1.err"
