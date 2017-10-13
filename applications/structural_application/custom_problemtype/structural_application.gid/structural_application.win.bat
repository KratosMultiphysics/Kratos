@ECHO OFF

del %2\%1.info
del %2\%1.flavia.res
del %2\%1.flavia.dat
del %2\%1.err

del %2\%1.txt
                    
del %2\%1_structure.prop 
del %2\%1_structure.elem 
del %2\%1_structure.cond 
del %2\%1_structure.init 
del %2\%1_structure.node


ren %2\%1.dat %2\%1.txt


ren %2\%1-1.dat %2\%1_structure.prop
ren %2\%1-2.dat %2\%1_structure.elem
ren %2\%1-3.dat %2\%1_structure.cond
ren %2\%1-4.dat %2\%1_structure.init
ren %2\%1-5.dat %2\%1_structure.node


D:\kratos\kratos\sources\Release\Kratos %1 > %2\%1.info

