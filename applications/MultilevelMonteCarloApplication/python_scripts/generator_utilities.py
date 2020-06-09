import numpy as np

# TODO: move in an external file and create a generator class, e.g. MonteCarlo.generator.Generate()

'''
function generatng the random variable
'''
def GenerateSample(problem_name):

    # Poisson equation
    if (problem_name == "poisson_square_2d"):
        alpha = 2.0
        beta = 6.0
        number_samples = 1
        sample = np.random.beta(alpha,beta,number_samples)

    # Compressible potential flow problem
    if (problem_name == "body_fitted_ellipse"):
        sample = []
        mean_Mach = 0.3
        std_deviation_Mach = 0.01
        number_samples = 1
        sample.append(np.random.normal(mean_Mach,std_deviation_Mach,number_samples))
        mean_angle_attack = 0.0 # [rad] = 0 [degrees] airfoil already has 5 degrees
        std_deviation_angle_attack = np.deg2rad(0.1)
        sample.append(np.random.normal(mean_angle_attack,std_deviation_angle_attack,number_samples))
        if sample[0] >= 1.0 or sample[0] <= 0.0 :
            raise Exception ("stochastic Mach number computed > 1 or < 0")

    # # Problem zero problem
    if (problem_name == "ProblemZero"):
        sample = []
        uref = 10 # wind speed of 10 m/s
        su = 1.0 #  a standard deviation of 10 % of the wind speed is considered
        sample.append(np.random.normal(uref,su))
        alpha0 = 0.12 # corresponds to open terrain from annexure to Euro Code DIN EN 1991-1-4 NA
        salpha = 0.012 # a standard deviation of 10 % is considered
        sample.append(np.random.normal(alpha0,salpha))
        y_0_mean = 0.02 # corresponds to "open country terrain," M. Andre's dissertation, p30 laso needs to be considered as uncertain ?
        sy_0 = 0.002 # a standard deviation of 10 % is considered
        sample.append(np.random.normal(y_0_mean, sy_0))
        u0_bar = 10 # wind speed of 10 m/s
        su = 1.0 # a standard deviation of 10 % of the wind speed is considered
        sample.append(np.random.normal(u0_bar,su))

    return sample
