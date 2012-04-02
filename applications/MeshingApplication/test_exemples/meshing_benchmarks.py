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

	Text += "Mapping_2d: "
	os.chdir("Mapping_2d.gid")
	sys.path.append(os.getcwd())

	print "Running Mapping_2d benchmark..."
	Msg = benchmarking.RunBenchmark("ProjectionTest_2D_benchmarking.py", "ProjectionTest_2D_benchmarking_ref.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "Mapping_2d benchmarking example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "Mapping_2d benchmarking example FAILED"

	os.chdir("..")
	

        ################################################################

	Text += "Mapping_3d: "
	os.chdir("Mapping_3d.gid")
	sys.path.append(os.getcwd())

	print "Running Mapping_3d benchmark..."
	Msg = benchmarking.RunBenchmark("ProjectionTest_3D_benchmarking.py", "ProjectionTest_3D_benchmarking_ref.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "Mapping_3d benchmarking example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "Mapping_3d benchmarking example FAILED"

	os.chdir("..")
	
        ################################################################

	Text += "Mapping_2d_BinBased: "
	os.chdir("Mapping_2d_BinBased.gid")
	sys.path.append(os.getcwd())

	print "Running Mapping_2d_BinBased benchmark..."
	Msg = benchmarking.RunBenchmark("ProjectionTestBinBased_2D_benchmarking.py", "ProjectionTestBinBased_2D_benchmarking_ref.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "Mapping_2d_BinBased benchmarking example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "Mapping_2d_BinBased benchmarking example FAILED"

	os.chdir("..")
	
       ################################################################

	Text += "Mapping_3d_BinBased: "
	os.chdir("Mapping_3d_BinBased.gid")
	sys.path.append(os.getcwd())

	print "Running Mapping_3d_BinBased benchmark..."
	Msg = benchmarking.RunBenchmark("ProjectionTestBinBased_3D_benchmarking.py", "ProjectionTestBinBased_3D_benchmarking_ref.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "Mapping_3d_BinBased benchmarking example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "Mapping_3d_BinBased benchmarking example FAILED"

	os.chdir("..")
        ################################################################
        ################################################################



        print "resume of all of the examples for the meshing application :"
        print Text
	return Text

Run();
