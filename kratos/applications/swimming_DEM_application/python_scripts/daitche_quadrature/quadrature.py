import math 
import cmath
import mpmath
import matplotlib.pyplot as plt

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
            return 8. / 15 * (n ** 2.5 - (n - 1) ** 2.5) + 2. / 3 * (-3 * n ** 1.5 + (n - 1) ** 1.5) + 2 * sqrt_n
        
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
    return 0.0

def Coefficient(order, n, j):
    if order == 1:
        return Alpha(n, j)
    elif order == 2:
        return Beta(n, j)
    else: 
        return Gamma(n, j)

def Daitche(order, times, f):
    sqrt_of_h = math.sqrt(times[-1] - times[-2])    
    n = len(times) + 1
    total = 0.0
    
    for j in range(1 , n):
        coefficient = Coefficient(order, n, j)
        total += coefficient * f(times[-j])
    
    return sqrt_of_h * total
    

t = 1.0
f = math.sin
n_div = [10, 40, 160, 640, 2560]
errors_naive = []
errors_1 = []
errors_2 = []

for n_divisions in n_div:
    h = t / n_divisions 
    times = [h * delta for delta in range(n_divisions + 1)]
    values = [ExactIntegrationOfSinusKernel(t) for t in times]
    errors_naive.append(abs((ApproximateQuadrature(times, f) - values[-1]) / values[-1]))
    errors_1.append(abs((Daitche(1, times, f) - values[-1]) / values[-1]))
    errors_2.append(abs((Daitche(2, times, f) - values[-1]) / values[-1]))
    
plt.plot(n_div, errors_naive, color='r')
plt.plot(n_div, errors_1, color='b')
plt.plot(n_div, errors_2, color='g')
plt.loglog()
plt.show()
#print("Naive: ", ApproximateQuadrature(times, f))
#print("First Order Daitche: ", Daitche(1, times, f))
#print("Second Order Daitche: ", Daitche(2, times, f))
#print("Third Order Daitche: ", Daitche(3, times, f))
#plt.plot(times, values)
#plt.show()