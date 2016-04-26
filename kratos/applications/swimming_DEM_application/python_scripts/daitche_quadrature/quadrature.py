import math 
import cmath
import mpmath
import matplotlib.pyplot as plt
from bigfloat import *

def ExactIntegrationOfSinusKernel(t):
    return 0.5 * math.pi * math.sqrt(t) * (mpmath.angerj(0.5, t) - mpmath.angerj(-0.5, t))
       

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
    if 0 < j and j < n:
        return 4. / 3 * ((j - 1) ** 1.5 + (j + 1) ** 1.5 - 2 * j ** 1.5)
    elif j == 0:
        return 4. / 3
    else:
        return 4. / 3 * ((n - 1) ** 1.5 - n ** 1.5 + 1.5 * math.sqrt(n))

def Beta(n, j):
    sqrt_2 = math.sqrt(2)
    sqrt_3 = math.sqrt(3)
    sqrt_n = math.sqrt(n) 
    
    if n >= 4:
        if 2 < j and j < n - 1:
            return 8. / 15 * (  (j + 2) ** 2.5 - 3 * (j + 1) ** 2.5 + 3 * j ** 2.5 - (j - 1) ** 2.5) +\
                    2. / 3 * (- (j + 2) ** 1.5 + 3 * (j + 1) ** 1.5 - 3 * j ** 1.5 + (j - 1) ** 1.5)
        elif j == 0:
            return 4. / 5 * sqrt_2
        elif j == 1:
            return  14. / 5 * sqrt_3 - 12. / 5 * sqrt_2
        elif j == 2:
            return 176. / 15 - 42. / 5 * sqrt_3 + 12. / 5 * sqrt_2
        elif j == n - 1:
            return 8. / 15 * (- 2 * n ** 2.5 + 3 * (n - 1) ** 2.5 - (n - 2) ** 2.5) +\
                    2. / 3 * (  4 * n ** 1.5 - 3 * (n - 1) ** 1.5 + (n - 2) ** 1.5)
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
    sqrt_of_h = math.sqrt(times[-1] - times[-2])    
    n = len(times) - 1
    total = 0.0
    
    for j in range(0 , n + 1):
        coefficient = Coefficient(order, n, j)
        total += coefficient * f(times[-j - 1])
    
    return sqrt_of_h * total

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
            print(richardsons)
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
print("Sizes: ", n_div)
exact_values = []
approx_values_naive = []
approx_values_1 = []
approx_values_2 = []
approx_values_3 = []
errors_naive = []
errors_1 = []
errors_2 = []
errors_3 = []

for n_divisions in n_div:
    h = t / n_divisions 
    times = [h * delta for delta in range(n_divisions + 1)]
    values = [float(ExactIntegrationOfSinusKernel(t)) for t in times]
    exact_value = values[-1]
    exact_values.append(exact_value)
    approx_value_naive = ApproximateQuadrature(times, f)
    approx_value_1 = Daitche(1, times, f)
    approx_value_2 = Daitche(2, times, f)
    approx_value_3 = Daitche(3, times, f)
    approx_values_naive.append(approx_value_naive)
    approx_values_1.append(approx_value_1)
    approx_values_2.append(approx_value_2)
    approx_values_3.append(approx_value_3)
    errors_naive.append(abs((approx_value_naive - exact_value) / exact_value))
    errors_1.append(abs((approx_value_1 - exact_value) / exact_value))
    errors_2.append(abs((approx_value_2 - exact_value) / exact_value))
    errors_3.append(abs((approx_value_3 - exact_value) / exact_value))

approx_values_naive_rich = [value for value in approx_values_naive] 
approx_values_1_rich = [value for value in approx_values_1] 
approx_values_2_rich = [value for value in approx_values_2] 
approx_values_3_rich = [value for value in approx_values_3] 
#approx_values_naive_rich[-1] = mpmath.richardson(approx_values_naive_rich)[0]
#approx_values_1_rich[-1] = mpmath.richardson(approx_values_1_rich)[0]
#approx_values_2_rich[-1] = mpmath.richardson(approx_values_2_rich)[0]
#approx_values_3_rich[-1] = mpmath.richardson(approx_values_3_rich)[0]
SubstituteRichardsons(approx_values_naive_rich, k, 0.5, 1)
SubstituteRichardsons(approx_values_1_rich, k, 2)
SubstituteRichardsons(approx_values_2_rich, k, 3)
SubstituteRichardsons(approx_values_3_rich, k, 4)
approx_values_naive_shank = [value for value in approx_values_naive]
approx_values_1_shank = [value for value in approx_values_1] 
approx_values_2_shank = [value for value in approx_values_2] 
approx_values_3_shank = [value for value in approx_values_3] 
SubstituteShanks(approx_values_naive_shank)
SubstituteShanks(approx_values_1_shank)
SubstituteShanks(approx_values_2_shank)
SubstituteShanks(approx_values_3_shank)
errors_naive_rich = [abs((approx_values_naive_rich[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_1_rich = [abs((approx_values_1_rich[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_2_rich = [abs((approx_values_2_rich[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_3_rich = [abs((approx_values_3_rich[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_naive_shank = [abs((approx_values_naive_shank[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_1_shank = [abs((approx_values_1_shank[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_2_shank = [abs((approx_values_2_shank[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
errors_3_shank = [abs((approx_values_3_shank[i] - exact_values[i]) / exact_values[i]) for i in range(len(exact_values))]
theoretical_slope_naive = []
theoretical_slope_naive = [errors_naive[0] * 0.5 ** math.sqrt(i) for i in range(len(n_div))]
theoretical_slope_1 = [errors_1[0] * 0.5 ** (2 * i) for i in range(len(n_div))]
theoretical_slope_2 = [errors_2[0] * 0.5 ** (3 * i) for i in range(len(n_div))]
theoretical_slope_3 = [errors_3[0] * 0.5 ** (4 * i) for i in range(len(n_div))]


plt.plot(n_div, errors_naive, color='r')
plt.plot(n_div, errors_1, color='b')
plt.plot(n_div, errors_2, color='g')
plt.plot(n_div, errors_3, color='k')
plt.plot(n_div, errors_naive_rich, color='r', linestyle = '--')
plt.plot(n_div, errors_1_rich, color='b', linestyle = '--')
plt.plot(n_div, errors_2_rich, color='g', linestyle = '--')
plt.plot(n_div, errors_3_rich, color='k', linestyle = '--')
plt.plot(n_div, errors_naive_shank, color='r', linestyle = '-.')
plt.plot(n_div, errors_1_shank, color='b', linestyle = '-.')
plt.plot(n_div, errors_2_shank, color='g', linestyle = '-.')
plt.plot(n_div, errors_3_shank, color='k', linestyle = '-.')
plt.plot(n_div, theoretical_slope_naive, color='r', linestyle = ':')
plt.plot(n_div, theoretical_slope_1, color='b', linestyle = ':')
plt.plot(n_div, theoretical_slope_2, color='g', linestyle = ':')
plt.plot(n_div, theoretical_slope_3, color='k', linestyle = ':')
plt.loglog()
#plt.show()
#print("Naive: ", ApproximateQuadrature(times, f))
#print("First Order Daitche: ", Daitche(1, times, f))
#print("Second Order Daitche: ", Daitche(2, times, f))
#print("Third Order Daitche: ", Daitche(3, times, f))
#plt.plot(times, values)
plt.show()