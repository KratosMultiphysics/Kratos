from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *


# PRESSURE CALCULATION


def ApplyPressure(Pressure, model_part, solver, SKIN, BOT, TOP, LAT, XLAT, XBOT, XTOP, XBOTCORNER, XTOPCORNER, alpha_top, alpha_bot, alpha_lat):

    print("Applying Pressure", "\n")

    skin_list = list()
    top_nodes_list = list()
    bot_nodes_list = list()
    total_cross_section = 0.0

    # Cylinder dimensions

    h = 0.3
    d = 0.15

    surface = 2 * (3.141592 * d * d * 0.25) +(3.141592*d*h)

    top_pressure = 0.0
    bot_pressure = 0.0

    for node in XLAT:

        r = node.GetSolutionStepValue(RADIUS, 0)
        x = node.X
        y = node.Y
        z = node.Z

        values = Array3()
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0

        cross_section = 3.141592 * r * r

        vect = zeros(3, double)

        # vector normal al centre:
        vect_moduli = sqrt(x * x + z * z)

        if(vect_moduli > 0.0):
            vect[0] = -x / vect_moduli
            vect[1] = 0
            vect[2] = -z / vect_moduli

        values[0] = cross_section * alpha_lat * Pressure * vect[0]
        values[1] = 0.0
        values[2] = cross_section * alpha_lat * Pressure * vect[2]

        node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)

    for node in XTOP:

        r = node.GetSolutionStepValue(RADIUS, 0)
        x = node.X
        y = node.Y
        z = node.Z

        values = Array3()
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0

        cross_section = 3.141592 * r * r

        values[0] = 0.0
        values[1] = -cross_section * alpha_top * Pressure
        values[2] = 0.0

        node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)

    for node in XBOT:

        r = node.GetSolutionStepValue(RADIUS, 0)
        x = node.X
        y = node.Y
        z = node.Z

        values = Array3()
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0

        cross_section = 3.141592 * r * r

        values[0] = 0.0
        values[1] = cross_section * alpha_bot * Pressure
        values[2] = 0.0

        node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)

    for node in XTOPCORNER:

        r = node.GetSolutionStepValue(RADIUS, 0)
        x = node.X
        y = node.Y
        z = node.Z

        values = Array3()
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0

        cross_section = 3.141592 * r * r

        vect = zeros(3, double)

        # vector normal al centre:
        vect_moduli = sqrt(x * x + z * z)

        if(vect_moduli > 0.0):
            vect[0] = -x / vect_moduli
            vect[1] = 0
            vect[2] = -z / vect_moduli

        values[0] = cross_section * alpha_lat * Pressure * vect[0] * 0.70710678
        values[1] = -cross_section * alpha_top * Pressure * 0.70710678
        values[2] = cross_section * alpha_lat * Pressure * vect[2] * 0.70710678

        node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)

    for node in XBOTCORNER:

        r = node.GetSolutionStepValue(RADIUS, 0)
        x = node.X
        y = node.Y
        z = node.Z

        values = Array3()
        values[0] = 0.0
        values[1] = 0.0
        values[2] = 0.0

        cross_section = 3.141592 * r * r

        vect = zeros(3, double)

        # vector normal al centre:
        vect_moduli = sqrt(x * x + z * z)

        if(vect_moduli > 0.0):
            vect[0] = -x / vect_moduli
            vect[1] = 0
            vect[2] = -z / vect_moduli

        values[0] = cross_section * alpha_lat * Pressure * vect[0] * 0.70710678
        values[1] = cross_section * alpha_bot * Pressure * 0.70710678
        values[2] = cross_section * alpha_lat * Pressure * vect[2] * 0.70710678

        node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE, values)
