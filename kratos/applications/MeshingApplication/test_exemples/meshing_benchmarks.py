import os
import sys

kratos_benchmarking_path = '../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

def Run():
	Msg = ""
	Text = "========== Meshing Aplication ==========\n"

        ################################################################

	Text += "adaptive_mesher2d: "
	os.chdir("adaptive_mesher2d.gid")
	sys.path.append(os.getcwd())

	print "Running Adaptive Mesher 2d benchmark..."
	Msg = benchmarking.RunBenchmark("remesh.py", "adaptive_mesher2d_ref.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "Adaptive_mesher2d benchmarking example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "Adaptive_mesher2d benchmarking example FAILED"

	os.chdir("..")

        ################################################################

	Text += "adaptive_mesher3d: "
	os.chdir("adaptive_mesher3d.gid")
	sys.path.append(os.getcwd())

	print "Running Adaptive Mesher 3d benchmark..."
	Msg = benchmarking.RunBenchmark("remesh.py", "adaptive_mesher3d_ref.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "Adaptive_mesher2d benchmarking example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "Adaptive_mesher2d benchmarking example FAILED"

	os.chdir("..")
	

        ################################################################

	Text += "ProjectionTest_2D: "
	os.chdir("Moving_Circle.gid")
	sys.path.append(os.getcwd())

	print "Running Moving_Circle benchmark..."
	Msg = benchmarking.RunBenchmark("ProjectionTest_2D_benchmarking.py", "ProjectionTest_2D_benchmarking_ref.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "Moving_Circle benchmarking example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "Moving_Circle benchmarking example FAILED"

	os.chdir("..")
	

        ################################################################

	Text += "ProjectionTest_3D: "
	os.chdir("Moving_Sphere.gid")
	sys.path.append(os.getcwd())

	print "Running Moving_Sphere benchmark..."
	Msg = benchmarking.RunBenchmark("ProjectionTest_3D_benchmarking.py", "ProjectionTest_3D_benchmarking_ref.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "Moving_Sphere benchmarking example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "Moving_Sphere benchmarking example FAILED"

	os.chdir("..")
	

        ################################################################
        ################################################################



        print "resume of all of the examples for the meshing application :"
        print Text
	return Text

Run();
