import math


def calculate_relative_water_pressure(y_coord, height, t_v):
    """
    Calculate the relative water pressure based on the analytical solution for one dimensional consolidation.
    See Program 16.1 in Section 16.2 "Solution" of the book "Soil Mechanics" by Arnold Verruijt (available for
    download from https://ocw.tudelft.nl/wp-content/uploads/SoilMechBook.pdf).

    The local variable names used here deliberately match the ones used by Program 16.1, to make it easy to
    verify the correctness of the code, despite the fact that the chosen names are not particularly descriptive.
    """
    # Variables defined at line with label 160
    tt = t_v
    a = 4.0 / math.pi
    pp = math.pi * math.pi / 4.0

    # Variables defined at line with label 170
    z = y_coord / height  # relative vertical position
    p = 0.0
    c = -1.0
    j = 0

    jt = 0.0
    while jt < 20.0:
        j += 1
        c *= -1.0
        jj = 2 * j - 1
        jt = jj * jj * pp * tt
        p += (a * c / jj) * math.cos(jj * math.pi * z / 2.0) * math.exp(-1.0 * jt)

    return p
    
def calculate_degree_of_1d_consolidation(t_v):
    """
    Calculates the analytical solution for degree of consolidation in 1d consolidation on linear elastic soil [@@ref@@]

    :param t_v:  dimensionless time factor
    :return   : degree of consolidation
    """
    convergence_criterion = 1e-10

    j = 1
    rel_d_old = 1
    rel_d = 0
    max_iterations = 1001
    min_iterations = 20
    while math.fabs(rel_d_old - rel_d) > convergence_criterion and j < max_iterations or j > min_iterations:

        rel_d_old = rel_d
        rel_d = 1 / (2 * j - 1) ** 2 * math.exp(-1 * (2 * j - 1) ** 2 * (math.pi ** 2) / 4 * t_v) + rel_d_old
        j += 1

    return 1.0 - 8.0 * rel_d / (math.pi ** 2)

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
    eps_h = (1/2 - 2/math.pi) * force * r ** 3 / (young*m_inertia)

    eps_v = (math.pi/4 - 2/math.pi) * force * r ** 3 / (young*m_inertia)
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
    return force * r * (1/math.pi - math.cos(theta)/2)
