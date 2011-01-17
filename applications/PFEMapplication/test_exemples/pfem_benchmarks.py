import os
import sys

kratos_benchmarking_path = '../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
	Msg = ""
	Text = "== PFEM ==========\n"
        
##	# dam2d
##
##	Text += "dam_2d: "
##	os.chdir("dam2d.gid")	
##	sys.path.append(os.getcwd())
##
##	print "running the pfem dam2d benchmark test..."
##        Msg = benchmarking.RunBenchmark("dam2d.py","dam2d_ref.txt")
##	
##	if (Msg == True):
##		Text += "OK\n"
##        	print "pfem dam2d test example succesful"
##	else:
##		Text += "FAILED\n"
##		Text += Msg
##		Text += "\n\n"
##        	print "pfem dam2d test example FAILED"
##
##
##	os.chdir("..")       


        ###############################################################################
	#dam 2D NON-NEWTONIAN
	
	Text += "dam_2d_Non-Newtonian: "
	os.chdir("dam2dNonNewt_Y1500.gid")	
	sys.path.append(os.getcwd())

	print "running the dam2dNonNewt benchmark test..."
        Msg = benchmarking.RunBenchmark("dam2dNonNewt.py", "dam2dNonNewt_ref.txt")

        if (Msg == True):
            Text += "OK\n"
            print "dam2dNonNewt example succesful"
        else:
            Text += "FAILED\n"
            Text += Msg
            Text += "\n\n"
            print "dam2dNonNewt example FAILED"


	os.chdir("..")
        ###############################################################################
	#slope 2D NON-NEWTONIAN COUPLED
	
	Text += "coupled_slope_2d_Non-Newtonian: "
	os.chdir("Coupled_slope2dNonNewt.gid")	
	sys.path.append(os.getcwd())

	print "running the CoupledNonNewtSlope benchmark test..."
        Msg = benchmarking.RunBenchmark("CoupledNonNewtSlope.py", "CoupledNonNewtSlope_ref.txt")

        if (Msg == True):
            Text += "OK\n"
            print "CoupledNonNewtSlope example succesful"
        else:
            Text += "FAILED\n"
            Text += Msg
            Text += "\n\n"
            print "CoupledNonNewtSlope example FAILED"




	os.chdir("..")
        ###############################################################################
	#slope 2D COUETTE NON-NEWTONIAN
	
	Text += "COUETTE_2d_Non-Newtonian: "
	os.chdir("CouetteNonNewtonian2d.gid")	
	sys.path.append(os.getcwd())

	print "running the Couette2dNonNewt benchmark test..."
        Msg = benchmarking.RunBenchmark("Couette2dNonNewt.py", "Couette2dNonNewt_ref.txt")

        if (Msg == True):
            Text += "OK\n"
            print "Couette2dNonNewt example succesful"
        else:
            Text += "FAILED\n"
            Text += Msg
            Text += "\n\n"
            print "Couette2dNonNewt example FAILED"




	os.chdir("..")
        ###############################################################################


	# Add other examples here



        ################################################################
        print "resume of all of the examples for the pfem application :"
        print Text
	return Text

Run();
