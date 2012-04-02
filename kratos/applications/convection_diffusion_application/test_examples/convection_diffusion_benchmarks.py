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

##	Text += "rotationgcone: "
##	os.chdir("square.gid")
##	sys.path.append(os.getcwd())
##
##	print "Running rotatingcone_PureConvection.py..."
##	Msg = benchmarking.RunBenchmark("rotatingcone_PureConvectionBenchmarking.py", "rotatingcone_PureConvection_ref.txt")	
##	
##	if (Msg == True):
##		Text += "OK\n"
##		print "square example succesful"
##	else:
##		Text += "FAILED\n"
##		Text += Msg
##		Text += "\n\n"
##		print "square example FAILED"
##
##	os.chdir("..")


	################################################################
	# testConvection Edgebased

	Text += "square_edgebased: "
	os.chdir("square_edgebased.gid")
	sys.path.append(os.getcwd())

	print "Running edgebased_PureConvection.py..."
	Msg = benchmarking.RunBenchmark("edgebased_PureConvection.py", "test_pureconvectionsolver_benchmarking_ref.txt")
	
	if (Msg == True):
		Text += "OK\n"
		print "testConvectionEdgebased example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "testConvectionEdgebased example FAILED"

	os.chdir("..")
        # Add other examples here


        
        ################################################################
        print "resume of all of the examples for the fluid application :"
        print Text
	return Text

Run();
