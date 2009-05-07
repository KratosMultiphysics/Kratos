
SET write_python_file=*GenData(Write_a_new_Python_script_file_for_Kratos,Int)
SET problemtype_name=*Tcl(GiD_Info Project ProblemType)

DEL %2\%problemtype_name%_var.py

RENAME %2\%1-3.dat %problemtype_name%_var.py

ECHO problem_name="%1" >> %problemtype_name%_var.py
ECHO problem_path="%2" >> %problemtype_name%_var.py
ECHO kratos_path="%KRATOS_PATH%" >> %problemtype_name%_var.py

IF %write_python_file%==1 COPY /Y %3\**.py %2\

REM IF EXIST %1.py (
REM  python $1.py %1 %2 %KRATOS_PATH% > %1.log
REM )
