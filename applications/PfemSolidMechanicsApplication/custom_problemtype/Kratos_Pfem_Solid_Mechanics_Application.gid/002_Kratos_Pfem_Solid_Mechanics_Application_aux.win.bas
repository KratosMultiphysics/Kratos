SET write_python_file=*GenData(Python_script_file)
SET file_location=*Tcl(getWinPyScript)
SET problemtype_name=*Tcl(GiD_Info Project ProblemType)
DEL problem_settings.py

RENAME %1-3.dat problem_settings.py


IF %write_python_file%==Use_Default (
 COPY /Y "%3\script.py" %2\
)

IF %write_python_file%==Copy_From (
 DEL %2\script.py
 COPY /Y "%file_location%" %2\script.py
)

IF EXIST script.py (
 echo "This script was using runkratos. Please contact the developer so he can update it to the new mechanism"
)
