import os
import sys

kratos_benchmarking_path = '../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
	Msg = ""
	Text = "========== DEM ==========\n"

        ###############################################################################
	
	Text += "DEM element test: "
	os.chdir("two_balls_no_damp.gid")	
	sys.path.append(os.getcwd())

	print "running the DEM two_balls_no_damp test..."
        Msg = benchmarking.RunBenchmark("two_balls_no_damp.py", "two_balls_no_damp_ref.txt")

        if (Msg == True):
            Text += "OK\n"
            print "two_balls_no_damp test succesful"
        else:
            Text += "FAILED\n"
            Text += Msg
            Text += "\n\n"
            print "two_balls_no_damp example test FAILED"


	os.chdir("..")
        ###############################################################################


	# Add other examples here



        ################################################################
        print "resume of all of the examples for the DEM application :"
        print Text
	return Text

if __name__=='__main__':
  Run();


