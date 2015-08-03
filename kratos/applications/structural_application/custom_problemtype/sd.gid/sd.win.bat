REM    %1: name of the current project
REM    %2: path of the current project
REM    %3: path of the problem type
REM    delete previous result file 
DEL "%2\*.post.res" 
DEL "%2\*.post.msh"
DEL "%2\*.post.bin"
COPY "%3\sed.exe" "%2"
COPY "%3\regex2.dll" "%2"
COPY "%3\libintl3.dll" "%2"
COPY "%3\libiconv2.dll" "%2"
REM renaming Kratos input files
python.exe "%3/clean_mdpa.py" "%2/%1.dat" "%2/%1.bak"
REN "%2/%1.bak" "%2/%1.mdpa"
DEL "%2/%1.dat"
REN "%2/%1-1.dat" "%2/%1.py"
REN "%2/%1-2.dat" "%2/%1_distributed_include.py"
REN "%2/%1-3.dat" "%2/%1_layers.py"
REN "%2/%1-4.dat" "%2/%1_shared_include.py"
REM append ess to python script
TYPE "%2\%1.ess" >> "%2\%1.py"
SED "s\rEpLaCeMeNtStRiNg\%1\g" < "%2\%1.py" > "%2\%1.py_changed"
REN "%2\%1.py_changed" "%2\%1.py"
CD "%2"
CD ..

