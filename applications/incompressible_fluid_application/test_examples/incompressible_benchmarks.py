import os
import sys

kratos_benchmarking_path = '../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

def Run():
	Msg = ""
	Text = "== Incompressible_Fluid ==========\n"

        ################################################################
	# column

	#Text += "column: "
	#os.chdir("column.gid")
	#sys.path.append(os.getcwd())

	##import column_benchmark
	##Msg = column_benchmark.Run()

	#print "Running column.py..."
	#Msg = benchmarking.RunBenchmark("column.py", "column_ref.txt")	
	
	#if (Msg == True):
	#	Text += "OK\n"
	#	print "colum example succesful"
	#else:
	#	Text += "FAILED\n"
	#	Text += Msg
	#	Text += "\n\n"
	#	print "colum example FAILED"

	#os.chdir("..")
        ################################################################
	# naca3d test for edge based level set solver
	try:
	    Text += "naca3d test: "
	    os.chdir("naca3d.gid")
	    sys.path.append(os.getcwd())
	except OSError:
	    print "Skipping naca3d benchmark: directory does not exist"
	    Text += "FAILED: directory not found\n"
	else:
	    print "running the benchmark for naca3d test..."
	    Msg = benchmarking.RunBenchmark("run_benchmark.py", "benchmark_reference_solution.txt")
	    
	    if (Msg == True):
		    Text += "OK\n"
		    print "naca3d test example succesful"
	    else:
		    Text += "FAILED\n"
		    Text += Msg
		    Text += "\n\n"
		    print "naca3d test example FAILED"

	    os.chdir("..")
	   
	
        ################################################################
	# mass conservation test for edge based level set solver

	try:
	    Text += "mass conservation test: "
	    os.chdir("mass_conservation.gid")
	    sys.path.append(os.getcwd())
	except OSError:
	    print "mass conservation test: directory does not exist"
	    Text += "FAILED: directory not found\n"
	else:
	    print "running the benchmark for mass_conservation test..."
	    Msg = benchmarking.RunBenchmark("run_benchmark.py", "benchmark_reference_solution.txt")
	    
	    if (Msg == True):
		    Text += "OK\n"
		    print "mass_conservation test example succesful"
	    else:
		    Text += "FAILED\n"
		    Text += Msg
		    Text += "\n\n"
		    print "mass_conservation test example FAILED"

	    os.chdir("..")
	
        # Add other examples here	

##        ################################################################
##	# cavity2D
##
##	Text += "cavity2d: "
##	os.chdir("cavity2d.gid")
##	sys.path.append(os.getcwd())
##
##	#import column_benchmark
##	#Msg = column_benchmark.Run()
##
##	print "Running cavity.py..."
##	Msg = benchmarking.RunBenchmark("cavity2d.py", "cavity_ref.txt")	
##	
##	if (Msg == True):
##		Text += "OK\n"
##		print "cavity2d example succesful"
##	else:
##		Text += "FAILED\n"
##		Text += Msg
##		Text += "\n\n"
##		print "cavity2d example FAILED"
##
##	os.chdir("..")
##	
##        # Add other examples here


        ###############################################################################
	# cavity2D
	try:
	    Text += "cavity2D: "
	    os.chdir("cavity2D.gid")
	    sys.path.append(os.getcwd())
	except OSError:
	    print "cavity2D test: directory does not exist"
	    Text += "FAILED: directory not found\n"
	else:
	    print "verifying  test_fractstep_cont_laplacian.py..."
	    Msg = benchmarking.RunBenchmark("test_fractstep_cont_laplacian.py", "fractstep_cont_laplacian_benchmarking_ref.txt")

	    if (Msg == True):
		Text += "OK\n"
		print "test_fractstep_cont_laplacian example succesful"
	    else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "test_fractstep_cont_laplacian example FAILED"

	    print "verifying  test_fractstep_discrete_laplacian.py..."
	    Msg = benchmarking.RunBenchmark("test_fractstep_discrete_laplacian.py", "fractstep_discrete_laplacian_benchmarking_ref.txt")


	    if (Msg == True):
		Text += "OK\n"
		print "test_fractstep_discrete_laplacian example succesful"
	    else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "test_fractstep_discrete_laplacian example FAILED"

	    os.chdir("..")
	
        ###############################################################################
	# cavity2D
	try:
	    Text += "cavity3D: "
	    os.chdir("cavity3D.gid")
	    sys.path.append(os.getcwd())
	except OSError:
	    print "cavity3D test: directory does not exist"
	    Text += "FAILED: directory not found\n"
	else:
	    print "verifying  test_fractstep_cont_laplacian.py...3D"
	    Msg = benchmarking.RunBenchmark("test_fractstep_cont_laplacian.py", "fractstep_cont_laplacian_benchmarking_ref.txt")

	    if (Msg == True):
		Text += "OK\n"
		print "test_fractstep_cont_laplacian 3D example succesful"
	    else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "test_fractstep_cont_laplacian 3D example FAILED"
	    os.chdir("..")
        
        ################################################################
	# cylinder
	try:
	    Text += "cylinder: "
	    os.chdir("cylinder.gid")
	    sys.path.append(os.getcwd())
	except OSError:
	    print "cylinder test: directory does not exist"
	    Text += "FAILED: directory not found\n"
	else:
	    #import column_benchmark
	    #Msg = column_benchmark.Run()

	    print "Running cylinder.py..."
	    Msg = benchmarking.RunBenchmark("run_example.py", "cylinder_ref.txt")	
	    
	    if (Msg == True):
		    Text += "OK\n"
		    print "cylinder example succesful"
	    else:
		    Text += "FAILED\n"
		    Text += Msg
		    Text += "\n\n"
		    print "cylinder example FAILED"

	    os.chdir("..")
	################################################################
	# cilinderGLS
	try:
	    Text += "cilinderGLS: "
	    os.chdir("cilinderGLS.gid")
	    sys.path.append(os.getcwd())
	except OSError:
	    print "cilinderGLS test: directory does not exist"
	    Text += "FAILED: directory not found\n"
	else:
	    #import column_benchmark
	    #Msg = column_benchmark.Run()

	    print "Running Cilinder GLS example with explicit Runge-Kutta and FRAC STEP..."
	    Msg = benchmarking.RunBenchmark("cil_gls.py", "cil_gls_ref.txt")	
	    
	    if (Msg == True):
		    Text += "OK\n"
		    print "cilinderGLS example succesful"
	    else:
		    Text += "FAILED\n"
		    Text += Msg
		    Text += "\n\n"
		    print "cilinderGLS example FAILED"

	    os.chdir("..")
        # Add other examples here
        ################################################################
	# dam2d
	try:
	    Text += "dam2d: "
	    os.chdir("dam2d.gid")
	    sys.path.append(os.getcwd())
	except OSError:
	    print "dam2d test: directory does not exist"
	    Text += "FAILED: directory not found\n"
	else:

	    print "Running dam2d.py..."
	    Msg = benchmarking.RunBenchmark("run_example.py", "dam2d_ref.txt")	
	    
	    if (Msg == True):
		    Text += "OK\n"
		    print "dam2d example succesful"
	    else:
		    Text += "FAILED\n"
		    Text += Msg
		    Text += "\n\n"
		    print "dam2d example FAILED"

	    os.chdir("..")


	################################################################
	# StillWater_Edgebased
	try:
	    Text += "StillWater_Edgebased: "
	    os.chdir("StillWater_Edgebased.gid")
	    sys.path.append(os.getcwd())
	except OSError:
	    print "StillWater_Edgebased test: directory does not exist"
	    Text += "FAILED: directory not found\n"
	else:	

	    print "Running StillWater_Edgebased..."
	    Msg = benchmarking.RunBenchmark("StillWater_Edgebased_script.py", "StillWater_Edgebased_ref.txt")	
	    
	    if (Msg == True):
		    Text += "OK\n"
		    print "StillWater_Edgebased example succesful"
	    else:
		    Text += "FAILED\n"
		    Text += Msg
		    Text += "\n\n"
		    print "StillWater_Edgebased example FAILED"

	    os.chdir("..")

	################################################################
	# StillWater_Edgebased_halfPorous
	try:
	    Text += "StillWater_Edgebased_halfPorous: "
	    os.chdir("StillWater_Edgebased_halfPorous.gid")
	    sys.path.append(os.getcwd())
	except OSError:
	    print "StillWater_Edgebased_halfPorous test: directory does not exist"
	    Text += "FAILED: directory not found\n"
	else:	
	    print "Running StillWater_Edgebased_halfPorous..."
	    Msg = benchmarking.RunBenchmark("StillWater_Edgebased_halfPorous_script.py", "StillWater_Edgebased_halfPorous_ref.txt")	
	    
	    if (Msg == True):
		    Text += "OK\n"
		    print "StillWater_Edgebased_halfPorous example succesful"
	    else:
		    Text += "FAILED\n"
		    Text += Msg
		    Text += "\n\n"
		    print "StillWater_Edgebased_halfPorous example FAILED"

	    os.chdir("..")


        ################################################################
	# cylinder

	#Text += "cylinder_3d: "
	#os.chdir("cylinder_3d.gid")
	#sys.path.append(os.getcwd())

	##import column_benchmark
	##Msg = column_benchmark.Run()

	#print "Running cylinder_3d.py..."
	#Msg = benchmarking.RunBenchmark("run_example.py", "cylinder_3d_ref.txt")	
	
	#if (Msg == True):
		#Text += "OK\n"
		#print "cylinder_3d example succesful"
	#else:
		#Text += "FAILED\n"
		#Text += Msg
		#Text += "\n\n"
		#print "cylinder_3d example FAILED"

	#os.chdir("..")


	
        ################################################################
	# cavityMonolithic2D
	try:
	    Text += "cavity_monolithic_3d: "
	    os.chdir("CavityMonolithic3D.gid")
	    sys.path.append(os.getcwd())
	except OSError:
	    print "cavity_monolithic_3d test: directory does not exist"
	    Text += "FAILED: directory not found\n"
	else:	
	    #import column_benchmark
	    #Msg = column_benchmark.Run()

	    print "Running script.py..."
	    Msg = benchmarking.RunBenchmark("script.py", "cavity_monolithic_3d_ref.txt")	
	    
	    if (Msg == True):
		    Text += "OK\n"
		    print "cavity_monolithic_3d example succesful"
	    else:
		    Text += "FAILED\n"
		    Text += Msg
		    Text += "\n\n"
		    print "cavity_monolithic_3d example FAILED"

	    os.chdir("..")

	# Add other examples here
        ################################################################
        print "resume of all of the examples for the fluid application :"
        print Text
	return Text

if __name__=='__main__':
  Run();
