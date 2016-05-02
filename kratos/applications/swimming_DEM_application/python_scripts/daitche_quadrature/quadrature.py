import math 
import cmath
import mpmath
import matplotlib.pyplot as plt
from bigfloat import *

def ExactIntegrationOfSinusKernel(t, a = None, b = None):
    with precision(300):
        if a == None and b == None:
            return 0.5 * math.pi * math.sqrt(t) * (mpmath.angerj(0.5, t) - mpmath.angerj(-0.5, t))
        if a == None and b != None:
            a = t
        elif b == None:
            b = 0.
        mpmath.mp.dps = 50
        mpmath.mp.pretty = True
        pi = mpmath.mp.pi
        pi = +pi
        fcos = mpmath.fresnelc
        fsin = mpmath.fresnels
        
        arg_a = mpmath.sqrt(2 * (t - a) / pi)
        arg_b = mpmath.sqrt(2 * (t - b) / mpmath.mp.pi) 
        return mpmath.sqrt(2 * mpmath.mp.pi) * ((fsin(arg_b) - fsin(arg_a)) * mpmath.cos(t) + (fcos(arg_a) - fcos(arg_b)) * mpmath.sin(t))
       

def ApproximateQuadrature(times, f):
    values = [0.0 for t in times]
    acc_sum = 2 * math.sqrt(times[-1] - times[-2]) * f(times[-1])
    
    for i in range(len(values) - 1):
        if i == 0:
            delta_t = times[1] - times[0]
        else:
            delta_t = times[i] - times[i - 1]
                    
        acc_sum += 0.5 * delta_t * (f(times[i]) + f(times[i - 1])) / math.sqrt(times[-1] - times[i])
        
    return acc_sum

def Alpha(n, j):
    with precision(300):
        four_thirds = BigFloat(4.) / BigFloat(3.)
        exponent = BigFloat(1.5)
        if 0 < j and j < n:
            return four_thirds * ((j - 1) ** exponent + (j + 1) ** exponent - 2 * j ** exponent)
        elif j == 0:
            return four_thirds
        else:
            return four_thirds * ((n - 1) ** exponent - n ** exponent + exponent * sqrt(n))

def Beta(n, j):
    sqrt_2 = math.sqrt(2)
    sqrt_3 = math.sqrt(3)
    sqrt_n = math.sqrt(n) 
    
    if n >= 4:
        if 2 < j and j < n - 1:
            return 8. / 15 * (  (j + 2) ** 2.5 - 3 * (j + 1) ** 2.5 + 3 * j ** 2.5 - (j - 1) ** 2.5)\
                   + 2. / 3 * (- (j + 2) ** 1.5 + 3 * (j + 1) ** 1.5 - 3 * j ** 1.5 + (j - 1) ** 1.5)
        elif j == 0:
            return 4. / 5 * sqrt_2
        elif j == 1:
            return  14. / 5 * sqrt_3 - 12. / 5 * sqrt_2
        elif j == 2:
            return 176. / 15 - 42. / 5 * sqrt_3 + 12. / 5 * sqrt_2
        elif j == n - 1:
            return 8. / 15 * (- 2 * n ** 2.5 + 3 * (n - 1) ** 2.5 - (n - 2) ** 2.5)\
                  + 2. / 3 * (  4 * n ** 1.5 - 3 * (n - 1) ** 1.5 + (n - 2) ** 1.5)
        else:
            return 8. / 15 * (n ** 2.5 - (n - 1) ** 2.5) + 2. / 3 * (- 3 * n ** 1.5 + (n - 1) ** 1.5) + 2 * sqrt_n
        
    elif n == 2:
        if j == 0:
            return 12. / 15 * sqrt_2
        elif j == 1:
            return 16. / 15 * sqrt_2
        else:
            return 2. / 15 * sqrt_2
        
    else:
        if j == 0:
            return 4. / 5 * sqrt_2
        elif j == 1:
            return 14. / 5 * sqrt_3 - 12. / 5 * sqrt_2
        elif j == 2:
            return - 8. / 5 * sqrt_3 + 12. / 5 * sqrt_2
        else:
            return 4. / 5 * sqrt_3 - 4. / 5 * sqrt_2

def Gamma(n, j):

    with precision(200):
        sqrt_2 = sqrt(2)
        sqrt_3 = sqrt(3)
        sqrt_5 = sqrt(5)
        sqrt_6 = sqrt(6) 
        sqrt_n = sqrt(n) 
        one = BigFloat(1)
        if n >= 7:
            #if 3 < j and j < n - 3:
                #answer = 16. / 105 * (    (j + 2) ** 3.5     + (j - 2) ** 3.5 - 4 * (j + 1) ** 3.5 - 4 * (j - 1) ** 3.5 + 6 * j ** 3.5)\
                            #+ 2. / 9 * (4 * (j + 1) ** 1.5 + 4 * (j - 1) ** 1.5     - (j + 2) ** 1.5     - (j - 2) ** 1.5 - 6 * j ** 1.5)
                #return float(answer)
            #elif j == 0:
                #answer = 244. / 315 * sqrt_2
                #return float(answer)
            #elif j == 1:
                #answer = 362. / 105 * sqrt_3 - 976. / 315 * sqrt_2
                #return float(answer)            
            #elif j == 2:
                #answer = 5584. / 315 - 1448. / 105 * sqrt_3 + 488. / 105 * sqrt_2
                #return float(answer)                        
            #elif j == 3:
                #answer = 1130. / 63 * sqrt_5 - 22336. / 315 + 724. / 35 * sqrt_3 - 976. / 315 * sqrt_2   
                #return float(answer)                        
            #elif j == n - 3:
                #answer = 16. / 105 * (n ** 3.5 - 4 * (n - 2) ** 3.5     + 6 * (n - 3) ** 3.5      - 4 * (n - 4) ** 3.5          + (n - 5) ** 3.5)\
                        #- 8. / 15 * n ** 2.5 + 4. / 9 * n ** 1.5 + 8. / 9 * (n - 2) ** 1.5 - 4. / 3 * (n - 3) ** 1.5 + 8. / 9 * (n - 4) ** 1.5 - 2. / 9 * (n - 5) ** 1.5 
                #return float(answer)                            
            #elif j == n - 2:
                #answer = 16. / 105 * ((n - 4) ** 3.5 - 4 * (n - 3) ** 3.5      + 6 * (n - 2) ** 3.5            - 3 * n ** 3.5)\
                            #+ 32. / 15 * n ** 2.5       - 2 * n ** 1.5 - 4. / 3 * (n - 2) ** 1.5 + 8. / 9 * (n - 3) ** 1.5 - 2. / 9 * (n - 4) ** 1.5
                #return float(answer)                                 
            #elif j == n - 1:
                #answer = 16. / 105 * (3 * n ** 3.5 - 4 * (n - 2) ** 3.5 + (n - 3) ** 3.5) - 8. / 3 * n ** 2.5 + 4 * n ** 1.5 + 8. / 9 * (n - 2) ** 1.5 - 2. / 9 * (n - 3) ** 1.5
                #return float(answer)                        
            #else:
                #answer = 16. / 105 * ((n - 2) ** 3.5 - n ** 3.5) + 16. / 15 * n ** 2.5 - 22. / 9 * n ** 1.5 - 2. / 9 * (n - 2) ** 1.5 + 2 * sqrt_n
                #return float(answer)     
            if 3 < j and j < n - 3:
                answer = 16. / (one*105) * (    (j + 2) ** (one*3.5)     + (j - 2) ** (one*3.5) - 4 * (j + 1) ** (one*3.5) - 4 * (j - 1) ** (one*3.5) + 6 * j ** (one*3.5))\
                            + 2. / (one*9) * (4 * (j + 1) ** (one*1.5) + 4 * (j - 1) ** (one*1.5)     - (j + 2) ** (one*1.5)     - (j - 2) ** (one*1.5) - 6 * j ** (one*1.5))
                return answer        
            elif j == 0:
                answer = 244. / (one*315) * sqrt_2
                return float(answer)
            elif j == 1:
                answer = 362. / (one*105) * sqrt_3 - 976. / (one*315) * sqrt_2
                return float(answer)            
            elif j == 2:
                answer = 5584. / (one*315) - 1448. / (one*105) * sqrt_3 + 488. / (one*105) * sqrt_2
                return float(answer)                        
            elif j == 3:
                answer = 1130. / (one*63) * sqrt_5 - 22336. / (one*315) + 724. / (one*35) * sqrt_3 - 976. / (one*315) * sqrt_2   
                return float(answer)                        
            elif j == n - 3:
                answer = 16. / (one*105) * (n ** (one*3.5) - 4 * (n - 2) ** (one*3.5)     + 6 * (n - 3) ** (one*3.5)      - 4 * (n - 4) ** (one*3.5)          + (n - 5) ** (one*3.5))\
                        - 8. / (one*15) * n ** (one*2.5) + 4. / (one*9) * n ** (one*1.5) + 8. / (one*9) * (n - 2) ** (one*1.5) - 4. / (one*3) * (n - 3) ** (one*1.5) + 8. / (one*9) * (n - 4) ** (one*1.5) - 2. / (one*9) * (n - 5) ** (one*1.5) 
                return float(answer)                            
            elif j == n - 2:
                answer = 16. / (one*105) * ((n - 4) ** (one*3.5) - 4 * (n - 3) ** (one*3.5)      + 6 * (n - 2) ** (one*3.5)            - 3 * n ** (one*3.5))\
                            + 32. / (one*15) * n ** (one*2.5)       - 2 * n ** (one*1.5) - 4. / (one*3) * (n - 2) ** (one*1.5) + 8. / (one*9) * (n - 3) ** (one*1.5) - 2. / (one*9) * (n - 4) ** (one*1.5)
                return float(answer)                                 
            elif j == n - 1:
                answer = 16. / (one*105) * (3 * n ** (one*3.5) - 4 * (n - 2) ** (one*3.5) + (n - 3) ** (one*3.5)) - 8. / (one*3) * n ** (one*2.5) + 4 * n ** (one*1.5) + 8. / (one*9) * (n - 2) ** (one*1.5) - 2. / (one*9) * (n - 3) ** (one*1.5)
                return float(answer)                        
            else:
                answer = 16. / (one*105) * ((n - 2) ** (one*3.5) - n ** (one*3.5)) + 16. / (one*15) * n ** (one*2.5) - 22. / (one*9) * n ** (one*1.5) - 2. / (one*9) * (n - 2) ** (one*1.5) + 2 * sqrt_n
                return float(answer)                        
            
        elif n == 3:
            if j == 0:
                return 68. / 105 * sqrt_3
            elif j == 1:
                return 6. / 7 * sqrt_3
            elif j == 2:
                return 12. / 35 * sqrt_3  
            else:
                return 16. / 105 * sqrt_3
        elif n == 4:
            if j == 0:
                return 244. / 315 * sqrt_2
            elif j == 1:
                return 1888. / 315 - 976. / 315 * sqrt_2
            elif j == 2:
                return - 656. / 105 + 488. / 105 * sqrt_2
            elif j == 3:
                return 544. / 105 - 976. / 315 * sqrt_2
            else:
                return - 292. / 315 + 244. / 315 * sqrt_2
        elif n == 5:
            if j == 0:
                return 244. / 315 * sqrt_2
            elif j == 1:
                return 362. / 105 * sqrt_3 - 976. / 315 * sqrt_2
            elif j == 2:
                return 500. / 63 * sqrt_5 - 1448. / 105 * sqrt_3 + 488. / 105 * sqrt_2
            elif j == 3:
                return - 290. / 21 * sqrt_5 + 724. / 35 * sqrt_3 - 976. / 315 * sqrt_2
            elif j == 4:
                return 220. / 21 * sqrt_5 - 1448. / 105 * sqrt_3 + 244. / 315 * sqrt_2            
            else:
                return - 164. / 63 * sqrt_5 + 362. / 105 * sqrt_3
        else:
            if j == 0:
                return 244. / 315 * sqrt_2
            elif j == 1:
                return 362. / 105 * sqrt_3 - 976. / 315 * sqrt_2
            elif j == 2:
                return 5584. / 315 - 1448. / 105 * sqrt_3 + 488. / 105 * sqrt_2
            elif j == 3:
                return 344. / 21 * sqrt_6 - 22336. / 315 + 724. / 35 * sqrt_3 - 976. / 315 * sqrt_2
            elif j == 4:
                return - 1188. / 35 * sqrt_6 + 11168. / 105 - 1448. / 105 * sqrt_3 + 244. / 315 * sqrt_2
            elif j == 5:
                return 936. / 35 * sqrt_6 - 22336. / 315 + 362. / 105 * sqrt_3
            else:
                return - 754. / 105 * sqrt_6 + 5584. / 315
        
def Coefficient(order, n, j):
    if order == 1:
        return Alpha(n, j)
    elif order == 2:
        return Beta(n, j)
    else: 
        return Gamma(n, j)

def Daitche(order, times, f):
    t = times[- 1]
    t0 = times[0]
    sqrt_of_h = math.sqrt(times[-1] - times[-2])    
    n = len(times) - 1
    total = 0.0   
    
    for j in range(0 , n + 1):
        coefficient = Coefficient(order, n, j)
        total += coefficient * f(times[-j - 1])
    
    return sqrt_of_h * total

def Phi(t):
    if t > 1e-10:
        answer = (exp(t) - 1) / t
    else:
        answer = 1 + 0.5 * t + 1. / 6 * t ** 2
    return answer

def Hinsberg(m, times, f):
    if len(times) < 4:
        return Bombardelli(times, f, 1)
    else:
        import hinsberg_optimization as op
        t_win = 0.4 * times[- 1]
        for i in range(len(times)):
            if times[i] >= t_win:
                break
        old_times = [time * times[i] / times[-1] for time in times]
        recent_times = [times[i] + time * (times[-1] - times[0] - times[i]) / (times[-1] - times[0]) for time in times]
        print(recent_times)
        F_win = Daitche(1, recent_times, f)
        F_win2 = Daitche(2, recent_times, f)
        F_win3 = Daitche(3, recent_times, f) 
        tis = [0.1, 0.3, 1., 3., 10., 40., 190., 1000., 6500., 50000.]
        tis = [0.1, 0.3, 0.5, 0.7, 0.8, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 20, 190, 500, 1000, 5000, 50000]        
        a0 = [0.2 for ti in tis]
        functional = op.Functional(tis)
        op.GetExponentialsCoefficients(functional, a0)
        F_tail = 0.
        contributions = [0.0 for coefficient in a0]

        for i in range(len(a0)):
            ti = tis[i]
            for k in range(2, len(old_times)):
                delta_t = old_times[k] - old_times[k - 1]                
                normalized_dt = delta_t / (2 * ti)
                normalized_t = old_times[k - 1] / (2 * ti)
                fn = f(old_times[k - 1])
                fn_plus_1 = f(old_times[k])
                Fdi = 2 * sqrt(exp(1.) * ti) * exp(- normalized_t) * (fn * (1 - Phi(- normalized_dt)) + fn_plus_1 * exp(- normalized_dt) * (Phi(normalized_dt) - 1.))
            
                if len(old_times) > 1:
                    Fre = exp(- normalized_dt) * contributions[i]
                else:
                    Fre = 0.
                contributions[i] = Fdi + Fre
                
            F_tail += a0[i] * contributions[i]
        final_time = times[- 1]
        t_win = recent_times[0]
        print("times",times)
        print("old_times", old_times)
        print("recent_times",recent_times)   
        print("exact over the whole history", ExactIntegrationOfSinusKernel(final_time))
        print("Daitche(1, times, f)", Daitche(1, times, f))        
        print("exact over recent times", ExactIntegrationOfSinusKernel(final_time, final_time, t_win))
        print("Daitche(1, recent_times, f)", Daitche(1, recent_times, f))
        print("Bombardelli(recent_times, f)", Bombardelli(recent_times, f, 2))            
        print("Daitche - OldDaitche", Daitche(1, times, f) -  Daitche(1, old_times, f))            
        print("Bombardelli - OldBombardelli", Bombardelli(times, f) - Bombardelli(old_times, f))                    
        print("exact over old times", ExactIntegrationOfSinusKernel(final_time, 0., t_win))
        print(F_tail)
        print(F_win)        
        print(F_win2) 
        print(F_win3) 
        print("F_win + F_tail", F_win + F_tail)
        #para
        return F_win + F_tail

def Bombardelli(times, f, order = 1):
    with precision(200):
        q = - 0.5
        t = times[- 1]
        a = times[0]
        if a == 0.5:
            print("TIMES", times)
        N = len(times) - 1
        h = (t - a) / N
        #initial_approx_deriv = 0.5 / h * (- f(a + 2 * h) + 4 * f(a + h) - 3 * f(a))          
        initial_approx_deriv = cos(a)
        constant_initial_correction = 2 * sqrt(t - a) * f(a)#t ** (- q) / gamma(1 - q) * f(a)
        linear_initial_correction = 2. / 3 * sqrt(t - a) * (a + 2 * t) * initial_approx_deriv 
        #linear_initial_correction = t ** (1 - q) / gamma(2 - q) * initial_approx_deriv
        linear_correction_option = 0

        if order == 1:
            coeff = h ** (- q)
            values = [gamma(k - q) / gamma(k + 1) * (f(t - k * h) - f(a)) for k in range(N)]
            initial_value_correction = constant_initial_correction
        else:
            coeff = h ** (- q) * gamma(- q)
            values = [(- 1) ** k * gamma(q + 1) / (gamma(k + 1) * gamma(q - k + 1)) * (f(t - (k  - 0.5 * q) * h) - f(a) - linear_correction_option * (t - (k  - 0.5 * q) * h - a) * initial_approx_deriv) for k in range(N)]
            initial_value_correction =  gamma(- q) * (constant_initial_correction + linear_correction_option * linear_initial_correction)
        return coeff * sum(values) + initial_value_correction

def SubstituteRichardsons(approx_successive_values, k, order, level = - 1):
    with precision(200):
        one = BigFloat(1)
        n = len(approx_successive_values)
        if level > n or level < 0:
            max_n = n
        else:
            max_n = level
        
        richardsons = [value for value in approx_successive_values]
        while max_n:
            max_n -= 1
            n -= 1
            new_richardsons = []
            
            for i in range(n):
                new_richardsons.append((k ** order * richardsons[i + 1] - richardsons[i]) / (one*(k ** order - 1)))
                
            richardsons = [value for value in new_richardsons]        
            
            for i in range(n):
                approx_successive_values[- i - 1] = richardsons[- i - 1]
            order += 1
            
def SubstituteEmpiricalRichardsons(approx_successive_values, k, order, level = - 1):
    with precision(200):
        one = BigFloat(1)
        n = len(approx_successive_values)
        if level > n or level < 0:
            max_n = n
        else:
            max_n = level
        
        richardsons = [value for value in approx_successive_values]
        while max_n:
            max_n -= 1
            n -= 1
            new_richardsons = []

            for i in range(n):
                approx_order = order
                
                if i > 0:
                    eiminus1 = abs((richardsons[i] - richardsons[i-1]) / (one*richardsons[i]))
                    ei = abs((richardsons[i+1] - richardsons[i]) / (one*richardsons[i+1]))
                    approx_order = max(log(eiminus1 / (one*ei)) / log(one*k), order)
                
                new_richardsons.append((k ** approx_order * richardsons[i + 1] - richardsons[i]) / (one*(k ** approx_order - 1)))

            richardsons = [value for value in new_richardsons]        
            
            for i in range(n):
                approx_successive_values[- i - 1] = richardsons[- i - 1]
            order += 1            

def SubstituteShanks(approx_sequence):
    with precision(200):
        one = BigFloat(1)
        my_list = approx_sequence
        shanks = [(my_list[i + 1] * my_list[i - 1] - my_list[i] ** 2) / (one*(my_list[i + 1] - 2 * my_list[i] + my_list[i - 1])) for i in range(1, len(my_list) - 1)]
        
        while len(shanks) > 2:        
            for i in range(len(shanks)):
                my_list[- i - 1] = shanks[- i - 1]
            temp_shanks = [(shanks[i + 1] * shanks[i - 1] - shanks[i] ** 2) / (one*(shanks[i + 1] - 2 * shanks[i] + shanks[i - 1])) for i in range(1, len(shanks) - 1)]
            shanks = temp_shanks
t = 1.0
f = math.sin
n_discretizations = 6
min_exp = 2
k = 2
n_div = [k ** (min_exp + i) for i in range(n_discretizations)]
m = 10
order_bomb = 2

print("Sizes: ", n_div)
exact_values = []
approx_values_naive = []
approx_values_1 = []
approx_values_2 = []
approx_values_3 = []
approx_values_bomb = []
approx_values_hins = []
errors_naive = []
errors_1 = []
errors_2 = []
errors_3 = []
errors_bomb = []
errors_hins = []

for n_divisions in n_div:
    h = t / n_divisions 
    times = [h * delta for delta in range(n_divisions)]
    times.append(t)
    values = [float(ExactIntegrationOfSinusKernel(t)) for t in times]
    exact_value = values[-1]
    exact_values.append(exact_value)
    approx_value_naive = ApproximateQuadrature(times, f)
    approx_value_1 = Daitche(1, times, f)
    approx_value_2 = Daitche(2, times, f)
    approx_value_3 = Daitche(3, times, f)
    approx_value_bomb = Bombardelli(times, f, order_bomb)
    approx_value_hins = Hinsberg(m, times, f)    
    approx_values_naive.append(approx_value_naive)
    approx_values_1.append(approx_value_1)
    approx_values_2.append(approx_value_2)
    approx_values_3.append(approx_value_3)
    approx_values_bomb.append(approx_value_bomb)
    approx_values_hins.append(approx_value_hins)    
    errors_naive.append(abs((approx_value_naive - exact_value) / exact_value))
    errors_1.append(abs((approx_value_1 - exact_value) / exact_value))
    errors_2.append(abs((approx_value_2 - exact_value) / exact_value))
    errors_3.append(abs((approx_value_3 - exact_value) / exact_value))
    errors_bomb.append(abs((approx_value_bomb - exact_value) / exact_value))
    errors_hins.append(abs((approx_value_hins - exact_value) / exact_value))    

approx_values_naive_rich = [value for value in approx_values_naive] 
approx_values_1_rich = [value for value in approx_values_1] 
approx_values_2_rich = [value for value in approx_values_2] 
approx_values_3_rich = [value for value in approx_values_3]
approx_values_bomb_rich = [value for value in approx_values_bomb]
approx_values_hins_rich = [value for value in approx_values_hins]
approx_values_naive_rich_emp = [value for value in approx_values_naive] 
approx_values_1_rich_emp = [value for value in approx_values_1] 
approx_values_2_rich_emp = [value for value in approx_values_2] 
approx_values_3_rich_emp = [value for value in approx_values_3]
approx_values_bomb_rich_emp = [value for value in approx_values_bomb]
approx_values_hins_rich_emp = [value for value in approx_values_hins] 

#approx_values_naive_rich[-1] = mpmath.richardson(approx_values_naive_rich)[0]
#approx_values_1_rich[-1] = mpmath.richardson(approx_values_1_rich)[0]
#approx_values_2_rich[-1] = mpmath.richardson(approx_values_2_rich)[0]
#approx_values_3_rich[-1] = mpmath.richardson(approx_values_3_rich)[0]
SubstituteRichardsons(approx_values_naive_rich, k, 0.5)
SubstituteRichardsons(approx_values_1_rich, k, 2)
SubstituteRichardsons(approx_values_2_rich, k, 3)
SubstituteRichardsons(approx_values_3_rich, k, 4)
SubstituteRichardsons(approx_values_bomb_rich, k, order_bomb)
SubstituteRichardsons(approx_values_hins_rich, k, 1)
SubstituteEmpiricalRichardsons(approx_values_naive_rich_emp, k, 0.5)
SubstituteEmpiricalRichardsons(approx_values_1_rich_emp, k, 2)
SubstituteEmpiricalRichardsons(approx_values_2_rich_emp, k, 3)
SubstituteEmpiricalRichardsons(approx_values_3_rich_emp, k, 4)
SubstituteEmpiricalRichardsons(approx_values_bomb_rich_emp, k, 2)
SubstituteEmpiricalRichardsons(approx_values_hins_rich_emp, k, 1)
approx_values_naive_shank = [value for value in approx_values_naive]
approx_values_1_shank = [value for value in approx_values_1] 
approx_values_2_shank = [value for value in approx_values_2] 
approx_values_3_shank = [value for value in approx_values_3]
approx_values_bomb_shank = [value for value in approx_values_bomb]
approx_values_hins_shank = [value for value in approx_values_hins] 
SubstituteShanks(approx_values_naive_shank)
SubstituteShanks(approx_values_1_shank)
SubstituteShanks(approx_values_2_shank)
SubstituteShanks(approx_values_3_shank)
SubstituteShanks(approx_values_bomb_shank)
SubstituteShanks(approx_values_hins_shank)
errors_naive_rich = [abs((approx_values_naive_rich[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_1_rich = [abs((approx_values_1_rich[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_2_rich = [abs((approx_values_2_rich[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_3_rich = [abs((approx_values_3_rich[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_bomb_rich = [abs((approx_values_bomb_rich[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_hins_rich = [abs((approx_values_hins_rich[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_naive_rich_emp = [abs((approx_values_naive_rich_emp[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_1_rich_emp = [abs((approx_values_1_rich_emp[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_2_rich_emp = [abs((approx_values_2_rich_emp[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_3_rich_emp = [abs((approx_values_3_rich_emp[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_bomb_rich_emp = [abs((approx_values_bomb_rich_emp[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_hins_rich_emp = [abs((approx_values_hins_rich_emp[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_naive_shank = [abs((approx_values_naive_shank[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_1_shank = [abs((approx_values_1_shank[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_2_shank = [abs((approx_values_2_shank[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_3_shank = [abs((approx_values_3_shank[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_bomb_shank = [abs((approx_values_bomb_shank[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_hins_shank = [abs((approx_values_hins_shank[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
theoretical_slope_naive = []
theoretical_slope_naive = [errors_naive[0] * 0.5 ** (i / 2) for i in range(len(n_div))]
theoretical_slope_1 = [errors_1[0] * 0.5 ** (2 * i) for i in range(len(n_div))]
theoretical_slope_2 = [errors_2[0] * 0.5 ** (3 * i) for i in range(len(n_div))]
theoretical_slope_3 = [errors_3[0] * 0.5 ** (4 * i) for i in range(len(n_div))]
theoretical_slope_bomb = [errors_bomb[0] * 0.5 ** (order_bomb * i) for i in range(len(n_div))]
theoretical_slope_hins = [errors_hins[0] * 0.5 ** (2 * i) for i in range(len(n_div))]


plt.plot(n_div, errors_naive, color='r')
plt.plot(n_div, errors_1, color='b')
plt.plot(n_div, errors_2, color='g')
plt.plot(n_div, errors_3, color='k')
plt.plot(n_div, errors_bomb, color='c')
plt.plot(n_div, errors_hins, color='m')
plt.plot(n_div, errors_naive_rich, color='r', linestyle = '--')
plt.plot(n_div, errors_1_rich, color='b', linestyle = '--')
plt.plot(n_div, errors_2_rich, color='g', linestyle = '--')
plt.plot(n_div, errors_3_rich, color='k', linestyle = '--')
plt.plot(n_div, errors_bomb_rich, color='c', linestyle = '--')
plt.plot(n_div, errors_hins_rich, color='m', linestyle = '--')
#plt.plot(n_div, errors_naive_rich_emp, color='r', linestyle = '-.')
#plt.plot(n_div, errors_1_rich_emp, color='b', linestyle = '-.')
#plt.plot(n_div, errors_2_rich_emp, color='g', linestyle = '-.')
#plt.plot(n_div, errors_3_rich_emp, color='k', linestyle = '-.')
#plt.plot(n_div, errors_bomb_rich_emp, color='c', linestyle = '-.')
#plt.plot(n_div, errors_hins_rich_emp, color='m', linestyle = '-.')
#plt.plot(n_div, errors_naive_shank, color='r', linestyle = '-.')
#plt.plot(n_div, errors_1_shank, color='b', linestyle = '-.')
#plt.plot(n_div, errors_2_shank, color='g', linestyle = '-.')
#plt.plot(n_div, errors_3_shank, color='k', linestyle = '-.')
#plt.plot(n_div, errors_bomb_shank, color='c', linestyle = '-.')
plt.plot(n_div, errors_hins_shank, color='m', linestyle = '-.')
plt.plot(n_div, theoretical_slope_naive, color='r', linestyle = ':')
plt.plot(n_div, theoretical_slope_1, color='b', linestyle = ':')
plt.plot(n_div, theoretical_slope_2, color='g', linestyle = ':')
plt.plot(n_div, theoretical_slope_3, color='k', linestyle = ':')
plt.plot(n_div, theoretical_slope_bomb, color='c', linestyle = ':')
plt.plot(n_div, theoretical_slope_hins, color='m', linestyle = ':')
plt.loglog()
#plt.show()
#print("Naive: ", ApproximateQuadrature(times, f))
#print("First Order Daitche: ", Daitche(1, times, f))
#print("Second Order Daitche: ", Daitche(2, times, f))
#print("Third Order Daitche: ", Daitche(3, times, f))
#plt.plot(times, values)
#final_time = 1.0
#h = 0.001
#order = 1
#optimal_a =  1 / (2 * (order + 1))
#best_a = optimal_a
#best_i = 0
#best_error = 999999999999999
#optimal_error = 0
#n_exponents = 20
#max_a = 3 * optimal_a
#delta_a = max_a / n_exponents
#def Id(x):
    #return 1.
#best_times = []

#for i in range(n_exponents):
    #a = delta_a * i
    #times = [h * i ** (a) for i in range(int(final_time / h))]
    #norm = sum(times)
    #times = [sum(times[:i]) / norm for i in range(len(times) + 1)]
    #new_error = ApproximateQuadrature(times, f) - exact_values[-1]
    #print(a)
    #print(new_error)
    #if best_error > new_error:
        #best_error = new_error
        #best_i = i
        #best_a = a
        #best_times = times
    #if int(a) == int(optimal_a):
        #optimal_error = new_error
        #optimal_times = times
#print("best a", best_a)
#print("optimal a ", optimal_a)
#print("best error", best_error)
#print("optimal_error", optimal_error)
#print("best times", best_times)
#print("optimal times", optimal_times)
#plt.plot(best_times, optimal_times)
plt.show()