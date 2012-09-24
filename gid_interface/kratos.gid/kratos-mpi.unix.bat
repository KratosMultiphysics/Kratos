REM @ECHO OFF
REM Identification for arguments
REM basename                          = %1
REM Project directory                 = %2
REM Problem directory                 = %3
REM GiD installation path             = %4
REM Number of processor               = %5

REM ECHO 1: %1; 2: %2; 3: %3; 4: %4 5: %5 > test.txt
 
REM OutputFile: %2\%1.info
REM ErrorFile: %2\%1.err
 
DEL %2\%1.info
DEL %2\%1.post.bin
DEL %2\%1.err

REM Run the python script
mpirun -np %5 python kratosMPI.py > %2\%1.info 2> %2\%1.err
