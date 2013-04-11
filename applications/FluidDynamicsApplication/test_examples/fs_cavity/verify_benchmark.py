import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

Text = " "

###############################################################################
print "Verifying fractional step element: cavity flow..."
Msg = benchmarking.RunBenchmark("fs_benchmark.py", "fs_benchmark_ref.txt")

if (Msg == True):
    Text += "OK\n"
    print "Fractional step element: cavity flow example succesful"
else:
    Text += "FAILED\n"
    Text += Msg
    Text += "\n\n"
    print "Fractional step element: cavity flow example FAILED"

print Text
