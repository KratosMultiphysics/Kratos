import matplotlib.pyplot as plt
import math
        
def ExtractVariableHistory(file_path, var_name, results):
    with open(file_path, mode = 'r') as f:
        indices = None
        first_line_already_read = False        
        for line in f:     
            if first_line_already_read:
                values = line.split()
                entries = [values[i] for i in indices]
                results.extend([float(entry) for entry in entries])
                
            else:
                first_line_already_read = True
                headers = line.split()
                i = 0
                for entry in headers:
                    if var_name in entry:
                        if entry[-2:] == '_X':
                            indices = [i, i + 1, i + 2]                            
                        else:
                            indices = [i]     
                        break
                    i += 1

def ComputeNormOfErrorInTime(reference_path, my_path, var_name):
    time = []
    reference_result = []
    my_result = []
    ExtractVariableHistory(reference_path, 'Time', time)
    ExtractVariableHistory(reference_path, var_name, reference_result)
    ExtractVariableHistory(my_path, var_name, my_result)
    n_time = len(time)
    n_var = len(reference_result) 
    
    is_scalar = False
    
    if n_time == n_var:
        is_scalar = True

    if not n_var == len(my_result):
        raise 'Both data sets must be of the same length!'
    
    if is_scalar:
        difference = [my_result[i] - reference_result[i] for i in range(n_var)]
        plt.plot(time, difference, label='Difference norm')
        plt.plot(time, reference_result, label='Reference solution norm')
        plt.plot(time, my_result, label='Approximation norm')     
        
    else:
        exact_norm = [math.sqrt(reference_result[3 * i] ** 2 + reference_result[3 * i + 1] ** 2 + reference_result[3 * i + 2] ** 2) for i in range(n_time)]
        my_norm = [math.sqrt(my_result[3 * i] ** 2 + my_result[3 * i + 1] ** 2 + my_result[3 * i + 2] ** 2) for i in range(n_time)]
        difference = [my_result[i] - reference_result[i] for i in range(n_var)]
        difference_norm = [math.sqrt(difference[3 * i] ** 2 + difference[3 * i + 1] ** 2 + difference[3 * i + 2] ** 2) for i in range(n_time)]
        plt.plot(time, difference_norm, label='Difference')
        plt.plot(time, exact_norm, label='Reference solution')
        plt.plot(time, my_norm, label='Approximation')      
    
    plt.xlabel('time')
    plt.ylabel(var_name)
    plt.legend()    
    plt.show()    

if __name__ == "__main__":
    node_id = 604
    var_name = 'DISPLACEMENT'
    my_path        = './fluid_imposed_Post_Files/results_node_' + str(node_id) + '.txt'
    reference_path = './Reference_Post_Files/results_node_' + str(node_id) + '.txt'    
    ComputeNormOfErrorInTime(reference_path, my_path, var_name)
