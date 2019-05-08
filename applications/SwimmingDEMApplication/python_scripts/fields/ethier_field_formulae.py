import sympy as sp
from sympy.printing import ccode
import math
import time as timer

time = sp.symbols('time')
coor = list(sp.symbols('coor:3'))
mA = sp.symbols('mA')
mD = sp.symbols('mD')

exp_d2t  = sp.exp(- mD * mD * time)
exp_ax   = sp.exp(mA * coor[0])
exp_ay   = sp.exp(mA * coor[1])
exp_az   = sp.exp(mA * coor[2])
sin_axdy = sp.sin(mA * coor[0] + mD * coor[1])
cos_axdy = sp.cos(mA * coor[0] + mD * coor[1])
sin_aydz = sp.sin(mA * coor[1] + mD * coor[2])
cos_aydz = sp.cos(mA * coor[1] + mD * coor[2])
sin_azdx = sp.sin(mA * coor[2] + mD * coor[0])
cos_azdx = sp.cos(mA * coor[2] + mD * coor[0])

p = - 0.5 * mA**2 * (exp_ax**2
                   + exp_ay**2
                   + exp_az**2
                   + 2 * sin_axdy * cos_azdx * exp_ay * exp_az
                   + 2 * sin_aydz * cos_axdy * exp_az * exp_ax
                   + 2 * sin_azdx * cos_aydz * exp_ax * exp_ay) * exp_d2t ** 2

def GetCppCode(title, formula, elapsed_time):
    code = str(ccode(formula))
    code = code.replace('coor0', 'coors[0]')
    code = code.replace('coor1', 'coors[1]')
    code = code.replace('coor2', 'coors[2]')
    code = code.replace('exp(', 'std::exp(')
    code = code.replace('pow(', 'std::pow(')
    code = code.replace('sin(', 'std::sin(')
    code = code.replace('cos(', 'std::cos(')
    print('Finished formula for ' + title + ' (elapsed time: ', '{:1f}'.format(elapsed_time), 's)', flush=True)
    return '\n' + title + ' =\n' + code + '\n'

start = timer.time()

with open('ethier_formulae.flae', 'w') as f:
    f.write(GetCppCode('p', sp.simplify(p), timer.time() - start))
    f.write(GetCppCode('dpdt', sp.simplify(sp.diff(p,time)), timer.time() - start))
    f.write(GetCppCode('dpd0', sp.simplify(sp.diff(p,coor[0])), timer.time() - start))
    f.write(GetCppCode('dpd1', sp.simplify(sp.diff(p,coor[1])), timer.time() - start))
    f.write(GetCppCode('dpd2', sp.simplify(sp.diff(p,coor[2])), timer.time() - start))
    f.write(GetCppCode('dpdtdt', sp.simplify(sp.diff(p,time,2)), timer.time() - start))
    f.write(GetCppCode('dpdtd0', sp.simplify(sp.diff(p,time,coor[0])), timer.time() - start))
    f.write(GetCppCode('dpdtd1', sp.simplify(sp.diff(p,time,coor[1])), timer.time() - start))
    f.write(GetCppCode('dpdtd2', sp.simplify(sp.diff(p,time,coor[2])), timer.time() - start))
    f.write(GetCppCode('dpd0d0', sp.simplify(sp.diff(p,coor[0], 2)), timer.time() - start))
    f.write(GetCppCode('dpd0d1', sp.simplify(sp.diff(p,coor[0], coor[1])), timer.time() - start))
    f.write(GetCppCode('dpd0d2', sp.simplify(sp.diff(p,coor[0], coor[2])), timer.time() - start))
    f.write(GetCppCode('dpd1d1', sp.simplify(sp.diff(p,coor[1], 2)), timer.time() - start))
    f.write(GetCppCode('dpd1d2', sp.simplify(sp.diff(p,coor[1], coor[2])), timer.time() - start))
    f.write(GetCppCode('dpd2d2', sp.simplify(sp.diff(p,coor[2], 2)), timer.time() - start))

