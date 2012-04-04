import os
import sys

kratos_benchmarking_path = '../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

def Run():
	Msg = ""
	Text = "========== Structural Aplication ==========\n"

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
        Text += "cantilever2dDynamic:"
	os.chdir("cantilever2d.gid")
	sys.path.append(os.getcwd())

	print "Running cantilever2ddynamic"
	Msg = benchmarking.RunBenchmark("cantilever2ddynamic_benchmarking.py", "cantilever2ddynamic.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "cantilever2dynamic succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "cantilever2dynamic example FAILED"

	os.chdir("..")

        ################################################################

       	Text += "cantilever2dstatic using MKLPardisoSolver: "
	os.chdir("cantilever2d.gid")
	sys.path.append(os.getcwd())

	print "Running cantilever2dstatic using MKLPardisoSolver"
	Msg = benchmarking.RunBenchmark("cantilever2dstatic_Mkl_benchmarking.py", "cantilever2dstatic_superlu.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "cantilever2dstatic using MKLPardisoSolver succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "cantilever2dstatic using MKLPardisoSolver example FAILED"

	os.chdir("..")

        ################################################################ 

	Text += "cantilever2dstatic_superlu: "
	os.chdir("cantilever2d.gid")
	sys.path.append(os.getcwd())

	print "Running cantilever2dstatic witdh SuperLUSolver"
	Msg = benchmarking.RunBenchmark("cantilever2dstatic_superlu_benchmarking.py", "cantilever2dstatic_superlu.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "cantilever2dstatic witdh SuperLUSolver example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "cantilever2dstatic witdh SuperLUSolver example FAILED"

	os.chdir("..")

        ################################################################
	Text += "cantilever2dstatic: "
	os.chdir("cantilever2d.gid")
	sys.path.append(os.getcwd())

	print "Running cantilever2dstatic with classical solver"
	Msg = benchmarking.RunBenchmark("cantilever2dstatic_benchmarking.py", "cantilever2dstatic_superlu.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "cantilever2dstatic with classical solver example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "cantilever2dstatic with classical solver example FAILED"

	os.chdir("..")
        ################################################################
	
	Text += "cantilever2dstatic_parallel:"
	os.chdir("cantilever2d.gid")
	sys.path.append(os.getcwd())

	print "Running cantilever2dstatic_parallel using ParallelSkylineLUFactorizationSolver"
	Msg = benchmarking.RunBenchmark("cantilever2dstatic_parallel_benchmarking.py", "cantilever2dstatic_superlu.txt")

	if (Msg == True):
		Text += "OK\n"
		print "cantilever2dstatic_parallel using ParallelSkylineLUFactorizationSolver  example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "cantilever2dstatic_parallel using ParallelSkylineLUFactorizationSolver example FAILED"

	os.chdir("..")


        ################################################################

	Text += "cantilever2dstatic using ParallelMKLPardisoSolver: "
	os.chdir("cantilever2d.gid")
	sys.path.append(os.getcwd())

	print "Running cantilever2dstatic using ParallelMKLPardisoSolver"
	Msg = benchmarking.RunBenchmark("cantilever2dstatic_parallel_pardiso_benchmarking.py", "cantilever2dstatic_superlu.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "cantilever2dstatic using ParallelMKLPardisoSolver  example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "cantilever2dstatic using ParallelMKLPardisoSolver example FAILED"

	os.chdir("..")

        ################################################################
        
        
	Text += "cantilever3dstatic: "
	os.chdir("cantilever3d.gid")
	sys.path.append(os.getcwd())

	print "Running cantilever3dstatic"
	Msg = benchmarking.RunBenchmark("cantilever3dstatic_superlu_benchmarking.py", "cantilever3dstatic.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "cantilever3dstatic example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "cantilever3dstatic example FAILED"

	os.chdir("..")


        ################################################################
        
        
	Text += "cantilever3ddynamic: "
	os.chdir("cantilever3d.gid")
	sys.path.append(os.getcwd())

	print "Running cantilever3dynamic"
	Msg = benchmarking.RunBenchmark("cantilever3ddynamic_benchmarking.py", "cantilever3ddynamic.txt")	
	
	if (Msg == True):
		Text += "OK\n"
		print "cantilever3ddynamic example succesful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "cantilever2dynamic  example FAILED"

	os.chdir("..")

        ################################################################



	Text += "balken contact benchmark: "
	os.chdir("balken.gid")
	sys.path.append(os.getcwd())

	print "Running balken contact benchmark..."
	Msg = benchmarking.RunBenchmark("balken_benchmarking.py","balken_ref.txt")

	if (Msg == True):
		Text += "OK\n"
		print "balken contact benchmark example successful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "balken contact benchmark example FAILED"

	os.chdir("..")


        ################################################################


	Text += "arc length desplacement benchmark: "
	os.chdir("arc_length_des.gid")
	sys.path.append(os.getcwd())

	print "Running arc length desplacement  benchmark..."
	Msg = benchmarking.RunBenchmark("arc_length_des_benchmarking.py","arc_length_des_benchmarking.txt")

	if (Msg == True):
		Text += "OK\n"
		print "arc_length_des benchmark example successful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "arc_length_des benchmark example FAILED"

	os.chdir("..")


        ################################################################


	Text += "arc length benchmark: "
	os.chdir("arc_length.gid")
	sys.path.append(os.getcwd())

	print "Running arc length  benchmark..."
	Msg = benchmarking.RunBenchmark("arc_length_benchmarking.py","arc_length_benchmarking.txt")

	if (Msg == True):
		Text += "OK\n"
		print "arc_length benchmark example successful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "arc_length benchmark example FAILED"

	os.chdir("..")


        ################################################################

	Text += "Pendulo Kratos Length benchmark: "
	os.chdir("Pendulo_Kratos.gid")
	sys.path.append(os.getcwd())

	print "Pendulo Kratos Length benchmark..."
	Msg = benchmarking.RunBenchmark("Pendulo_Kratos_benchmarking.py","Pendulo_Kratos_benchmarking.txt")

	if (Msg == True):
		Text += "OK\n"
		print "Pendulo Kratos benchmark example successful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "Pendulo Kratos benchmark example FAILED"

	os.chdir("..")

        ################################################################

	Text += "PlasticitJ2 with EBST collection benchmarks: "
	os.chdir("plasticJ2.gid")
	sys.path.append(os.getcwd())

	print "PlasticitJ2 with EBST TENSION benchmarks..."
	Msg = benchmarking.RunBenchmark("tension.py","tension_ref.txt")

	if (Msg == True):
		Text += "OK\n"
		print "PlasticitJ2 with EBST TENSION benchmarks successful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "PlasticitJ2 with EBST TENSION benchmarks FAILED"

	print "PlasticitJ2 with EBST TORSION benchmarks..."
	Msg = benchmarking.RunBenchmark("torsion.py","torsion_ref.txt")

	if (Msg == True):
		Text += "OK\n"
		print "PlasticitJ2 with EBST TORSION benchmarks successful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "PlasticitJ2 with EBST TORSION benchmarks FAILED"

	print "PlasticitJ2 with EBST VERTICAL benchmarks..."
	Msg = benchmarking.RunBenchmark("vertical.py","vertical_ref.txt")

	if (Msg == True):
		Text += "OK\n"
		print "PlasticitJ2 with EBST VERTICAL benchmarks successful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "PlasticitJ2 with EBST VERTICAL benchmarks FAILED"

	print "PlasticitJ2 with EBST FORCE benchmarks..."
	Msg = benchmarking.RunBenchmark("force.py","force_ref.txt")

	if (Msg == True):
		Text += "OK\n"
		print "PlasticitJ2 with EBST FORCE benchmarks successful"
	else:
		Text += "FAILED\n"
		Text += Msg
		Text += "\n\n"
		print "PlasticitJ2 with EBST FORCE benchmarks FAILED"

	os.chdir("..")

        ################################################################



        print "Resume of all of the examples for the structural application :"
        print Text
	return Text

if __name__=='__main__':
  Run();

