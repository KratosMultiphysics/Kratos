def calculate_1D_consolidation(y_coord, H, T_v):
    """
    Calculates the analytical solution for 1d consolidation on linear elastic soil [@@ref@@]

    :param y_coord: vertical coordinate
    :param H: sample height
    :param T_v:  dimensionless time factor
    :return: relative excess pore pressure
    """
    from math import fabs, cos, pi, exp

    convergence_criterion = 1e-10

    j = 1
    rel_p_old = 1
    rel_p = 0
    max_iterations = 1001
    min_iterations=20
    min_iterations_reached = False
    while fabs(rel_p_old - rel_p) > convergence_criterion and j < max_iterations or not min_iterations_reached:

        rel_p_old = rel_p
        rel_p = (-1) ** (j - 1) / (2 * j - 1) * cos((2 * j - 1) * pi / 2 * y_coord / H) * exp(
            -1 * (2 * j - 1) ** 2 * pi ** 2 / 4 * T_v) + rel_p_old
        j += 1

        if (j > min_iterations):
            min_iterations_reached = True

    rel_p = 4.0 / pi * rel_p
    return rel_p


def rigid_footing(x, B, delta, G, nu, settlement):
    """
    Calculates analytical solution for reaction pressure of settlement controlled rigid footing on linear elastic soil
    [@@ref@@]

    :param x: x-coordinate
    :param B: width footing
    :param delta: geometry dependent factor
    :param G:  shear modulus
    :param nu: poison ratio
    :param settlement:  settlement

    :return: vertical reaction pressure
    """
    from math import pi, sqrt

    reaction_force = settlement * 2.0 * (1.0 + nu) * G / delta
    sigma_v = -2.0 / pi * reaction_force / 2.0 / (B * sqrt(1.0 - (x / B) ** 2.0))

    return sigma_v


def calculate_max_deflections_ring(force, r, young, m_inertia):
    """
    todo Extend description
    ref: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.912.904&rep=rep1&type=pdf
    :param force: Point load
    :param r: radius
    :param young: Young's modulus
    :param m_inertia: area moment of inertia
    :return: relative increase in horizontal and vertical diameter
    """
    from math import pi
    eps_h = (1/2 - 2/pi) * force * r ** 3 / (young*m_inertia)

    eps_v = (pi/4 - 2/pi) * force * r ** 3 / (young*m_inertia)
    return eps_h, eps_v

def calculate_bending_moments_ring(force, r, theta):
    """
    todo extend description
    ref http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.912.904&rep=rep1&type=pdf
    :param force: point load
    :param r: radius
    :param theta: angle
    :return: bending moment
    """
    from math import cos, pi
    moment = force * r * (1/pi - cos(theta)/2)
    return moment
