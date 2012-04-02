import os
import sys

kratos_benchmarking_path = '../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

def Run():
	Msg = ""
	Text = "========== Structural Aplication ==========\n"

        ################################################################

	#Text += "Patch_Test_Total_Lagrangian_3D_8N: "
	#os.chdir("Patch_Test_Total_Lagrangian_3D_8N.gid")
	#sys.path.append(os.getcwd())

	#print "Running Patch_Test_Total_Lagrangian_3D_8N..."
	#Msg = benchmarking.RunBenchmark("Patch_Test_Total_Lagrangian_3D_8N_benchmarking.py", "Patch_Test_Total_Lagrangian_3D_8N.txt")	
	
	#if (Msg == True):
		#Text += "OK\n"
		#print "Patch_Test_Total_Lagrangian_3D_8N_benchmarking example succesful"
	#else:
		#Text += "FAILED\n"
		#Text += Msg
		#Text += "\n\n"
		#print "Patch_Test_Total_Lagrangian_3D_8N_benchmarking.py example FAILED"

	#os.chdir("..")

	

        ################################################################


	#Text += "Patch_Test_Total_Lagrangian_4N: "
	#os.chdir("Patch_Test_Total_Lagrangian_4N.gid")
	#sys.path.append(os.getcwd())

	#print "Running Patch_Test_Total_Lagrangian_4N..."
	#Msg = benchmarking.RunBenchmark("Patch_Test_Total_Lagrangian_4N_benchmarking.py", "Patch_Test_Total_Lagrangian_4N.txt")	
	
	#if (Msg == True):
		#Text += "OK\n"
		#print "Patch_Test_Total_Lagrangian_4N_benchmarking example succesful"
	#else:
		#Text += "FAILED\n"
		#Text += Msg
		#Text += "\n\n"
		#print "Patch_Test_Total_Lagrangian_4N_benchmarking.py example FAILED"

	#os.chdir("..")


        
        ################################################################


	#Text += "cantilever2dstatic_superlu: "
	#os.chdir("cantilever2d.gid")
	#sys.path.append(os.getcwd())

	#print "Running cantilever2dstatic_superlu..."
	#Msg = benchmarking.RunBenchmark("cantilever2dstatic_superlu_benchmarking.py", "cantilever2dstatic_superlu.txt")	
	
	#if (Msg == True):
		#Text += "OK\n"
		#print "cantilever2dstatic_superlu example succesful"
	#else:
		#Text += "FAILED\n"
		#Text += Msg
		#Text += "\n\n"
		#print "cantilever2dstatic_superlu example FAILED"

	#os.chdir("..")


        ################################################################
	
#	Text += "cantilever2dstatic_parallel: "
#	os.chdir("cantilever2d.gid")
#	sys.path.append(os.getcwd())
#
#	print "Running cantilever2dstatic_parallel..."
#	Msg = benchmarking.RunBenchmark("cantilever2dstatic_parallel_benchmarking.py", "cantilever2dstatic_superlu.txt")
#
#	if (Msg == True):
#		Text += "OK\n"
#		print "cantilever2dstatic_parallel example succesful"
#	else:
#		Text += "FAILED\n"
#		Text += Msg
#		Text += "\n\n"
#		print "cantilever2dstatic_parallel example FAILED"
#
#	os.chdir("..")


        ################################################################

	#Text += "cantilever2dstatic_parallel_pardiso: "
	#os.chdir("cantilever2d.gid")
	#sys.path.append(os.getcwd())

	#print "Running cantilever2dstatic_parallel_pardiso..."
	#Msg = benchmarking.RunBenchmark("cantilever2dstatic_parallel_pardiso_benchmarking.py", "cantilever2dstatic_superlu.txt")	
	
	#if (Msg == True):
		#Text += "OK\n"
		#print "cantilever2dstatic_parallel_pardiso example succesful"
	#else:
		#Text += "FAILED\n"
		#Text += Msg
		#Text += "\n\n"
		#print "cantilever2dstatic_parallel example FAILED"

	#os.chdir("..")


        ################################################################

	#Text += "balken contact benchmark: "
	#os.chdir("balken.gid")
	#sys.path.append(os.getcwd())

	#print "Running balken contact benchmark..."
	#Msg = benchmarking.RunBenchmark("balken_benchmarking.py","balken_ref.txt")

	#if (Msg == True):
		#Text += "OK\n"
		#print "balken contact benchmark example successful"
	#else:
		#Text += "FAILED\n"
		#Text += Msg
		#Text += "\n\n"
		#print "balken contact benchmark example FAILED"

	#os.chdir("..")


        ################################################################


	#Text += "arc length desplacement benchmark: "
	#os.chdir("arc_length_des.gid")
	#sys.path.append(os.getcwd())

	#print "Running arc length desplacement  benchmark..."
	#Msg = benchmarking.RunBenchmark("arc_length_des_benchmarking.py","arc_length_des_benchmarking.txt")

	#if (Msg == True):
		#Text += "OK\n"
		#print "arc_length_des benchmark example successful"
	#else:
		#Text += "FAILED\n"
		#Text += Msg
		#Text += "\n\n"
		#print "arc_length_des benchmark example FAILED"

	#os.chdir("..")


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

	#Text += "Pendulo Kratos Length benchmark: "
	#os.chdir("Pendulo_Kratos.gid")
	#sys.path.append(os.getcwd())

	#print "Pendulo Kratos Length benchmark..."
	#Msg = benchmarking.RunBenchmark("Pendulo_Kratos_benchmarking.py","Pendulo_Kratos_benchmarking.txt")

	#if (Msg == True):
		#Text += "OK\n"
		#print "Pendulo Kratos benchmark example successful"
	#else:
		#Text += "FAILED\n"
		#Text += Msg
		#Text += "\n\n"
		#print "Pendulo Kratos benchmark example FAILED"

	#os.chdir("..")


        ################################################################

	#Text += "PlasticitJ2 with EBST collection benchmarks: "
	#os.chdir("plasticJ2.gid")
	#sys.path.append(os.getcwd())

	#print "PlasticitJ2 with EBST TENSION benchmarks..."
	#Msg = benchmarking.RunBenchmark("tension.py","tension_ref.txt")

	#if (Msg == True):
		#Text += "OK\n"
		#print "PlasticitJ2 with EBST TENSION benchmarks successful"
	#else:
		#Text += "FAILED\n"
		#Text += Msg
		#Text += "\n\n"
		#print "PlasticitJ2 with EBST TENSION benchmarks FAILED"

	#print "PlasticitJ2 with EBST TORSION benchmarks..."
	#Msg = benchmarking.RunBenchmark("torsion.py","torsion_ref.txt")

	#if (Msg == True):
		#Text += "OK\n"
		#print "PlasticitJ2 with EBST TORSION benchmarks successful"
	#else:
		#Text += "FAILED\n"
		#Text += Msg
		#Text += "\n\n"
		#print "PlasticitJ2 with EBST TORSION benchmarks FAILED"

	#print "PlasticitJ2 with EBST VERTICAL benchmarks..."
	#Msg = benchmarking.RunBenchmark("vertical.py","vertical_ref.txt")

	#if (Msg == True):
		#Text += "OK\n"
		#print "PlasticitJ2 with EBST VERTICAL benchmarks successful"
	#else:
		#Text += "FAILED\n"
		#Text += Msg
		#Text += "\n\n"
		#print "PlasticitJ2 with EBST VERTICAL benchmarks FAILED"

	#print "PlasticitJ2 with EBST FORCE benchmarks..."
	#Msg = benchmarking.RunBenchmark("force.py","force_ref.txt")

	#if (Msg == True):
		#Text += "OK\n"
		#print "PlasticitJ2 with EBST FORCE benchmarks successful"
	#else:
		#Text += "FAILED\n"
		#Text += Msg
		#Text += "\n\n"
		#print "PlasticitJ2 with EBST FORCE benchmarks FAILED"

	#os.chdir("..")

        ################################################################



        print "resume of all of the examples for the fluid application :"
        print Text
	return Text

Run();
