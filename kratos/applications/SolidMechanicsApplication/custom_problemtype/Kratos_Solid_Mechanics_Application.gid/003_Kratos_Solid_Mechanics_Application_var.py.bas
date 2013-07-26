domain_size = *GenData(DOMAIN_SIZE)

#Problem Data
#####################################
ProblemType = "*GenData(Problem_Type)"
SolverType  = "*GenData(Solver_Type)"
*format "%10.5e"
time_step   = *GenData(Time_Step)
*format "%8i"
nsteps = *GenData(Number_of_steps,INT)
*format "%8i"
NumberofThreads = *GenData(Number_of_threads,INT)

#Solver Data
#####################################
LinearSolver = "*GenData(Linear_Solver)"
*format "%10.5e"
Linear_Solver_Tolerance = *GenData(Linear_Solver_Tolerance)
*format "%8i"
Linear_Solver_Max_Iteration = *GenData(Linear_Solver_Max_Iteration,INT)
Convergence_Criteria = "*GenData(Convergence_Criteria)"
*format "%10.5e"
Convergence_Tolerance = *GenData(Convergence_Tolerance)
*format "%10.5e"
Absolute_Tolerance = *GenData(Absolute_Tolerance)
*format "%4i"
Max_Iter = *GenData(Max_Iter,INT)

#Constraints Data
#####################################
Incremental_Load = "*GenData(Incremental_Load)"
Incremental_Displacement = "*GenData(Incremental_Displacement)"

#PostProcess Data
#####################################
echo_level = *GenData(Echo_Level)
WriteResults = "*GenData(WriteResults)"
WriteMesh = "*GenData(WriteMesh)"
WriteParticles = "*GenData(WriteParticles)"
WriteConditions = "*GenData(WriteConditions)"
WriteFrequency = *GenData(WriteFrequency)
PlotGraphs = "*GenData(Plot_Graphs)"
PlotFrequency = *GenData(PlotFrequency)
PrintLists = "*GenData(Print_List_Files)"
number_of_lists = *GenData(List_Files)
file_list = []
*for(i=1;i<=GenData(List_Files,INT);i=i+1)
file_list.append(*GenData(List_Files,*i))
*end
SaveRestart = "*GenData(Print_Restart)"
Restart_Interval = *GenData(Restart_Interval)
LoadRestart = "*GenData(Load_Restart)"
Restart_Step = *GenData(Load_Step)

# Declare Python Variables

