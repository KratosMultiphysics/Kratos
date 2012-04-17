import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

Text = " "

###############################################################################
print "verifying  connectivity_preserve_modeler benchmark..."
Msg = benchmarking.RunBenchmark("do_test.py", "connectivity_preserve_benchmarking_ref.txt")

if (Msg == True):
    Text += "OK\n"
    print "connectivity_preserve_modeler benchmark succesful"
else:
    Text += "FAILED\n"
    Text += Msg
    Text += "\n\n"
    print "connectivity_preserve_modeler benchmark FAILED"


print Text
