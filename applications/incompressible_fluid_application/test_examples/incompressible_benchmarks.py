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
	# mass conservation test for edge based level set solver

	Text += "mass conservation test: "
	os.chdir("mass_conservation.gid")
	sys.path.append(os.getcwd())

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

        ################################################################
	# cavity2D

	Text += "cavity2d: "
	os.chdir("cavity2d.gid")
	sys.path.append(os.getcwd())

	#import column_benchmark
	#Msg = column_benchmark.Run()

	print "Running cavity.py..."
	Msg = benchmarking.RunBenchmark("cavity2d.py", "cavity_ref.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "cavity2d example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "cavity2d example FAILED"

	os.chdir("..")
	
        # Add other examples here

        ################################################################
	# cylinder

	Text += "cylinder: "
	os.chdir("cylinder.gid")
	sys.path.append(os.getcwd())

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

	Text += "cilinderGLS: "
	os.chdir("cilinderGLS.gid")
	sys.path.append(os.getcwd())

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

	Text += "dam2d: "
	os.chdir("dam2d.gid")
	sys.path.append(os.getcwd())

	#import column_benchmark
	#Msg = column_benchmark.Run()

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
	# QuietWater_LevelSet

	Text += "QuietWater_LevelSet: "
	os.chdir("QuietWater_LevelSet.gid")
	sys.path.append(os.getcwd())

	print "Running QuietWater_LevelSet..."
	Msg = benchmarking.RunBenchmark("split_level_set_QuietWater.py", "QuietWater_LS_ref.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "QuietWater_LevelSet example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "QuietWater_LevelSet example FAILED"

	os.chdir("..")

        ################################################################
	# cylinder

	Text += "cylinder_3d: "
	os.chdir("cylinder_3d.gid")
	sys.path.append(os.getcwd())

	#import column_benchmark
	#Msg = column_benchmark.Run()

	print "Running cylinder_3d.py..."
	Msg = benchmarking.RunBenchmark("run_example.py", "cylinder_3d_ref.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "cylinder_3d example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "cylinder_3d example FAILED"

	os.chdir("..")


	
        ################################################################


	# Add other examples here
        ################################################################
        print "resume of all of the examples for the fluid application :"
        print Text
	return Text

Run();
