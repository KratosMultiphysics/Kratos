import sys

kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)

import benchmarking

def Run():
  print "Running VMS test: parabolic flow in a trapezoidal domain solved using OSS"
  return benchmarking.RunBenchmark("test.py", "trapezoid_exact.txt")

if __name__=='__main__':
  if Run():
    print "Test successful"
  else:
    print "Test failed"
