from exaqute.ExaquteTaskLocal import *

@ExaquteTask(returns=1)
def check_vector(*collection_in):
    base_string = ""
    for elem in from_args_to_vector(collection_in):
        base_string += str(len(elem))
    return base_string

if __name__ == "__main__":
    vec = [[1,2,3],[4,5]]
    result = check_vector(*(from_vector_to_args(vec)))
