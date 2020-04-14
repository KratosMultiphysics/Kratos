import CoSimIO

# print(CoSimIO.__str__)
# help(CoSimIO)



u = 15

def func_print_u():
    print(u)


# CoSimIO.FunctionCapture(func_print_u)

def ProceedToNextTimeStep(new_time):


    return new_time + 0.1

def SolveSolutionStep():
    print("SolveSolutionStep!", u)

def SolveSolutionStep2(some_var):
    print("SolveSolutionStep!")

connection_name = "Python_solver"

CoSimIO.Connect(connection_name, "py_solver_settings.txt")

# CoSimIO.FunctionCapture(SolveSolutionStep2)
CoSimIO.Register_SolveSolutionStep(connection_name, SolveSolutionStep)


# CoSimIO.Run(connection_name)


CoSimIO.Disconnect(connection_name)