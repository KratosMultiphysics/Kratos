import os
import sys

kratos_benchmarking_path = '../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

def Run():
	Msg = ""
	Text = "== Convection_Diffusion ==========\n"
	
        ################################################################
	# rotatingcone_PureConvection

	Text += "rotationgcone: "
	os.chdir("square.gid")
	sys.path.append(os.getcwd())

	#import column_benchmark
	#Msg = column_benchmark.Run()

	print "Running rotatingcone_PureConvection.py..."
	Msg = benchmarking.RunBenchmark("rotatingcone_PureConvectionBenchmarking.py", "rotatingcone_PureConvection_ref.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "square example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "square example FAILED"

	os.chdir("..")
	
        # Add other examples here


        
        ################################################################
        print "resume of all of the examples for the fluid application :"
        print Text
	return Text

Run();
