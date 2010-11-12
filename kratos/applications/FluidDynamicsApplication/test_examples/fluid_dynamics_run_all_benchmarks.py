import os
import sys

kratos_benchmarking_path = '../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
	Msg = ""
	Text = "== Fluid Dynamics ==========\n"

        ###############################################################################
	#VMS2D element test
	
	Text += "VMS2D element test: "
	os.chdir("vms2d_test")	
	sys.path.append(os.getcwd())

	print "running the vms2d_test benchmark test..."
        Msg = benchmarking.RunBenchmark("script_elemtest.py", "vms2d_test_ref.txt")

        if (Msg == True):
            Text += "OK\n"
            print "VMS2D element test succesful"
        else:
            Text += "FAILED\n"
            Text += Msg
            Text += "\n\n"
            print "VMS2D element test FAILED"


	os.chdir("..")
        ###############################################################################


	# Add other examples here



        ################################################################
        print "resume of all of the examples for the Fluid Dynamics application :"
        print Text
	return Text

Run();

