import sys

kratos_benchmarking_path = '../../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)

import benchmarking

def Run():
  print "Running test of cylinder with slip BC, wall law and outlet close to the cylinder"
  return benchmarking.RunBenchmark("run_test.py", "drag_reference.txt")

if __name__=='__main__':
  Msg = Run()
  if( Msg == True ):
    print "Test successful"
  else:
    print "Test failed"
    print Msg
