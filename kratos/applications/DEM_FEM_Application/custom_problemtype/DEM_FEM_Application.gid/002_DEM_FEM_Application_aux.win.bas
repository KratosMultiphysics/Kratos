REM ECHO continuamos en otro.bat >> c:\temp\pp.txt
SET write_python_file=*GenData(Python_script_file)
SET problemtype_name=*Tcl(GiD_Info Project ProblemType)

DEL %2\%problemtype_name%_var.py

RENAME %2\%1-3.dat %problemtype_name%_var.py

ECHO problem_name="%1" >> %problemtype_name%_var.py
ECHO problem_path="%2" >> %problemtype_name%_var.py
ECHO kratos_path="%KRATOS_PATH%" >> %problemtype_name%_var.py

REM ECHO problem_name="%1" >> c:\temp\pp.txt
REM ECHO problem_path="%2" >> c:\temp\pp.txt
REM ECHO kratos_path="%KRATOS_PATH%" >> c:\temp\pp.txt
REM echo pepito=%3 >> c:\temp\pp.txt

REM IF %write_python_file%==1 COPY /Y %3\**.py %2\

REM IF EXIST %1.py (
REM  python $1.py %1 %2 %KRATOS_PATH% >& %1.info
REM )

copy %3\script.py %2\

python %2\script.py > %2\%1.info
