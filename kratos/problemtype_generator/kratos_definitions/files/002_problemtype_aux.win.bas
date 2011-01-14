SET write_python_file=*GenData(Python_script_file)
SET file_location=*Tcl(getWinPyScript)
SET problemtype_name=*Tcl(GiD_Info Project ProblemType)
DEL %2\%problemtype_name%_var.py

RENAME %2\%1-3.dat %problemtype_name%_var.py

ECHO problem_name="%1" >> %problemtype_name%_var.py
ECHO problem_path="%2" >> %problemtype_name%_var.py
ECHO kratos_path="%KRATOS_PATH%" >> %problemtype_name%_var.py

COPY /Y "%3\analysis.py" %2\

IF %write_python_file%==Use_Default (
 COPY /Y "%3\script.py" %2\
)

IF %write_python_file%==Copy_From (
 REM ECHO in_2nd_if > %2\%1.info
 REM ECHO "%file_location%" > %2\%1.info
 DEL %2\script.py
 COPY /Y "%file_location%" %2\script.py
)

IF EXIST script.py (
REM CALL C:\python25\python.exe script.py %1 %2 %KRATOS_PATH% > %2\%1.info 2> %2\%1.err
REM python script.py %1 %2 %KRATOS_PATH% > %2\%1.info 2> %2\%1.err
CMD /c %2\script.py %1 %2 %KRATOS_PATH% > %2\%1.info 2> %2\%1.err
)
