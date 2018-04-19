from KratosMultiphysics import *
from KratosMultiphysics.ThermoMechanicalApplication import *

def InitializeTables(model_part):
    tb1 = PiecewiseLinearTable(TEMPERATURE, DENSITY)
    tb1.AddRow(0.0000000000e+00,      3.000000000e+03)
    model_part.AddTable(1, tb1)

    tb2 = PiecewiseLinearTable(TEMPERATURE, SPECIFIC_HEAT)
    tb2.AddRow(0.0000000000e+00,      10.0000000000e+02)
    tb2.AddRow(5000.0000000000e+00,      10.0000000000e+02)
    model_part.AddTable(2, tb2)

    tb3 = PiecewiseLinearTable(TEMPERATURE, SOLIDFRACTION)
    tb3.AddRow(0.0000000,      1.0000000000e+00)
    tb3.AddRow(600.0000000,      1.0000000000e+00)
    tb3.AddRow(610.0000000,      0.0000000000e+00)
    tb3.AddRow(5000.0000000,      0.0000000000e+00)
    model_part.AddTable(3, tb3)

    tb4 = PiecewiseLinearTable(TEMPERATURE, SOLIDFRACTION_RATE)
    tb4.AddRow(0.000,      0.0000000000e+00)
    tb4.AddRow(600.0,      0.0000000000e+00)
    tb4.AddRow(600.0,      -0.1000000000e+00)
    tb4.AddRow(610.0,      -0.1000000000e+00)
    tb4.AddRow(610.0,      0.0000000000e+00)
    tb4.AddRow(1000.0,     0.0000000000e+00)
    model_part.AddTable(4, tb4)

    tb5 = PiecewiseLinearTable(TEMPERATURE, CONDUCTIVITY)
    tb5.AddRow(0.0000000000e+00,      1.5000000000e+02)
    tb4.AddRow(5000.0000000000e+00,      1.5000000000e+02)
    model_part.AddTable(5, tb5)

    tb6 = PiecewiseLinearTable(TEMPERATURE, VISCOSITY)
    tb6.AddRow(0.0,      1.0000000000e-03)
    tb6.AddRow(5.000000000e+03,      1.0000000000e-03)
    model_part.AddTable(6, tb6)

    tb7 = PiecewiseLinearTable(TEMPERATURE, HTC)
    tb7.AddRow(0.0000000000e+00,      2.0000000000e+03)
    tb7.AddRow(5000.0000000000e+00,      2.0000000000e+03)
    model_part.AddTable(7, tb7)
