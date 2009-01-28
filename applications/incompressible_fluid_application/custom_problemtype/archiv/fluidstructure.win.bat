@ECHO OFF

del %2\%1.info
del %2\%1.flavia.res
del %2\%1.flavia.dat
del %2\%1.err

del %2\%1.txt

del %2\%1_fluid.prop     
del %2\%1_fluid.elem     
del %2\%1_fluid.cond     
del %2\%1_fluid.init     
del %2\%1_fluid.node     
                    
del %2\%1_structure.prop 
del %2\%1_structure.elem 
del %2\%1_structure.cond 
del %2\%1_structure.init 
del %2\%1_structure.node

del coupling.dat


ren %2\%1.dat %2\%1.txt

ren %2\%1-1.dat %2\%1_fluid.prop
ren %2\%1-2.dat %2\%1_fluid.elem
ren %2\%1-3.dat %2\%1_fluid.cond
ren %2\%1-4.dat %2\%1_fluid.init
ren %2\%1-5.dat %2\%1_fluid.node


ren %2\%1-6.dat %2\%1_structure.prop
ren %2\%1-7.dat %2\%1_structure.elem
ren %2\%1-8.dat %2\%1_structure.cond
ren %2\%1-9.dat %2\%1_structure.init
ren %2\%1-10.dat %2\%1_structure.node

ren %2\%1-11.dat coupling.dat

D:\kratos\kratos\sources\Release\Kratos %1 > %2\%1.info

