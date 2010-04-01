@ECHO OFF

rem OutputFile: %1.log

del %2\%1.info
del %2\%1.flavia.res
del %2\%1.flavia.dat
del %2\%1.err

del %2\%1.py

del %2\%1.prop
del %2\%1.elem
del %2\%1.cond
del %2\%1.init
del %2\%1.node

ren %2\%1.dat %2\%1.py
ren %2\%1-1.dat %2\%1.prop
ren %2\%1-2.dat %2\%1.elem
ren %2\%1-3.dat %2\%1.cond
ren %2\%1-4.dat %2\%1.init
ren %2\%1-5.dat %2\%1.node

D:\kratosR1\Python25\python %2\%1.py >%2\%1.log

del %2\%1.post.res
ren %2\%1_0.post.bin %2\%1.post.res


