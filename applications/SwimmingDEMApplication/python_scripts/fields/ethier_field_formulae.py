import sympy as sp
from sympy.printing import ccode
import math

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

def GetCppCode(formula):
    code = str(ccode(formula))
    code = code.replace('exp(', 'std::exp(')
    code = code.replace('pow(', 'std::pow(')
    code = code.replace('sin(', 'std::sin(')
    code = code.replace('cos(', 'std::cos(')
    return code

# print('p = ', ccode(sp.simplify(p)))
# print('dpdt = ', ccode(sp.simplify(sp.diff(p,time))))
# print('dpd0 = ', GetCppCode(sp.simplify(sp.diff(p,coor[0]))))
# print('dpd1 = ', GetCppCode(sp.simplify(sp.diff(p,coor[1]))))
# print('dpd2 = ', GetCppCode(sp.simplify(sp.diff(p,coor[2]))))
# print('dpdtdt = ', GetCppCode(sp.simplify(sp.diff(p,time,2))))
# print('dpdtd0 = ', GetCppCode(sp.simplify(sp.diff(p,time,coor[0]))))
# print('dpdtd1 = ', GetCppCode(sp.simplify(sp.diff(p,time,coor[1]))))
# print('dpdtd2 = ', GetCppCode(sp.simplify(sp.diff(p,time,coor[2]))))
# print('dpd0d0 = ', GetCppCode(sp.simplify(sp.diff(p,coor[0], 2))))
# print('dpd0d1 = ', GetCppCode(sp.simplify(sp.diff(p,coor[0], coor[1]))))
# print('dpd0d2 = ', GetCppCode(sp.simplify(sp.diff(p,coor[0], coor[2]))))
# print('dpd1d1 = ', GetCppCode(sp.simplify(sp.diff(p,coor[1], 2))))
# print('dpd1d2 = ', GetCppCode(sp.simplify(sp.diff(p,coor[1], coor[2]))))
print('dpd2d2 = ', GetCppCode(sp.simplify(sp.diff(p,coor[2], 2))))

