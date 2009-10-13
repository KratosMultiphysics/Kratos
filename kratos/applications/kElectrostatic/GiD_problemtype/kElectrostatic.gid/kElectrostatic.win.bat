@ECHO OFF

rem OutputFile: %1.log

del %2\%1.info
del %2\%1.flavia.res
del %2\%1.flavia.dat
del %2\%1.err

del %2\%1.py

del %2\%1.py
del %2\kElectrostatic.py
del %2\%1.mdpa

ren %2\%1.dat %2\%1.py
ren %2\%1-1.dat %2\kElectrostatic.py
ren %2\%1-2.dat %2\%1.mdpa

rem D:\kratosR1\Python25\python %2\%1.py >%2\%1.log
D:\Python26\python %2\%1.py >%2\%1.log

del %2\%1.post.res
ren %2\%1_0.post.bin %2\%1.post.res


