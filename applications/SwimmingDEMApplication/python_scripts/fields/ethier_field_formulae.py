import sympy as sp
from sympy.printing import ccode
import math

time = sp.symbols('time')
coor = list(sp.symbols('coor:3'))
mA = sp.symbols('mA')
mD = sp.symbols('mD')

mExpD2T  = sp.exp(- mD * mD * time)
mExpAX   = sp.exp(mA * coor[0])
mExpAY   = sp.exp(mA * coor[1])
mExpAZ   = sp.exp(mA * coor[2])
mSinAXDY = sp.sin(mA * coor[0] + mD * coor[1])
mCosAXDY = sp.cos(mA * coor[0] + mD * coor[1])
mSinAYDZ = sp.sin(mA * coor[1] + mD * coor[2])
mCosAYDZ = sp.cos(mA * coor[1] + mD * coor[2])
mSinAZDX = sp.sin(mA * coor[2] + mD * coor[0])
mCosAZDX = sp.cos(mA * coor[2] + mD * coor[0])

p = - 0.5 * mA**2 * (mExpAX**2
                   + mExpAY**2
                   + mExpAZ**2
                   + 2 * mSinAXDY * mCosAZDX * mExpAY * mExpAZ
                   + 2 * mSinAYDZ * mCosAXDY * mExpAZ * mExpAX
                   + 2 * mSinAZDX * mCosAYDZ * mExpAX * mExpAY) * mExpD2T ** 2

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

