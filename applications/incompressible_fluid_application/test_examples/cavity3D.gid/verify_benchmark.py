import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

Text = " "

###############################################################################
print "verifying  test_fractstep_cont_laplacian.py...3D case"
Msg = benchmarking.RunBenchmark("test_fractstep_cont_laplacian.py", "fractstep_cont_laplacian_benchmarking_ref.txt")

if (Msg == True):
    Text += "OK\n"
    print "test_fractstep_cont_laplacian 3D example succesful"
else:
    Text += "FAILED\n"
    Text += Msg
    Text += "\n\n"
    print "test_fractstep_cont_laplacian 3D example FAILED"


print Text
