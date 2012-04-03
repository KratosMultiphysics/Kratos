import sys
kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

Text = " "

################################################################
# StillWater_Elembased_halfPorous

print "Running StillWater_Elembased_halfPorous..."
Msg = benchmarking.RunBenchmark("StillWater_Elembased_halfPorous_script.py", "StillWater_Elembased_halfPorous_ref.txt")	

if (Msg == True):
        Text += "OK\n"
        print "StillWater_Elembased_halfPorous example succesful"
else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print "StillWater_Elembased_halfPorous example FAILED"

print Text
