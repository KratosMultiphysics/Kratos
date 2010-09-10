import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

Text = " "

###############################################################################
print "verifying  dam2dNonNewt.py..."
Msg = benchmarking.RunBenchmark("dam2dNonNewt.py", "dam2dNonNewt_ref.txt")

if (Msg == True):
    Text += "OK\n"
    print "dam2dNonNewt example SUCCESFUL"
else:
    Text += "FAILED\n"
    Text += Msg
    Text += "\n\n"
    print "dam2dNonNewt example FAILED"

print Text
