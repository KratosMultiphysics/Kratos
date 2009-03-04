import os
import sys

kratos_benchmarking_path = '../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

def Run():
	Msg = ""
	Text = "== Structural Aplication ==========\n"

        ################################################################

	Text += "Patch_Test_Total_Lagrangian_3D_8N: "
	os.chdir("Patch_Test_Total_Lagrangian_3D_8N.gid")
	sys.path.append(os.getcwd())

	print "Running Patch_Test_Total_Lagrangian_3D_8N..."
	Msg = benchmarking.RunBenchmark("Patch_Test_Total_Lagrangian_3D_8N_benchmarking.py", "Patch_Test_Total_Lagrangian_3D_8N.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "Patch_Test_Total_Lagrangian_3D_8N_benchmarking example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "Patch_Test_Total_Lagrangian_3D_8N_benchmarking.py example FAILED"

	os.chdir("..")

	

        ################################################################


	Text += "Patch_Test_Total_Lagrangian_4N: "
	os.chdir("Patch_Test_Total_Lagrangian_4N.gid")
	sys.path.append(os.getcwd())

	print "Running Patch_Test_Total_Lagrangian_4N..."
	Msg = benchmarking.RunBenchmark("Patch_Test_Total_Lagrangian_4N_benchmarking.py", "Patch_Test_Total_Lagrangian_4N.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "Patch_Test_Total_Lagrangian_4N_benchmarking example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "Patch_Test_Total_Lagrangian_4N_benchmarking.py example FAILED"

	os.chdir("..")



        
        ################################################################
        print "resume of all of the examples for the fluid application :"
        print Text
	return Text

Run();
