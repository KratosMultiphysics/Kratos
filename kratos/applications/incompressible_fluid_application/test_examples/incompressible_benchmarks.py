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

	Text += "column: "
	os.chdir("column.gid")
	sys.path.append(os.getcwd())

	#import column_benchmark
	#Msg = column_benchmark.Run()

	print "Running column.py..."
	Msg = benchmarking.RunBenchmark("column.py", "column_ref.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "colum example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "colum example FAILED"

	os.chdir("..")

	

        ################################################################
	# cavity2D

	Text += "cavity2d: "
	os.chdir("cavity2d.gid")
	sys.path.append(os.getcwd())

	#import column_benchmark
	#Msg = column_benchmark.Run()

	print "Running column.py..."
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
        print "resume of all of the examples for the fluid application :"
        print Text
	return Text

Run();
