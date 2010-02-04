import os
import sys

kratos_benchmarking_path = '../../../benchmarking' 
sys.path.append(kratos_benchmarking_path)
import benchmarking

def Run():
    Msg = ""
    Text = "== FSI_application ==========\n"

    # Example benchmark problem:
    ################################################################
    # column

    #Text += "column: "
    #os.chdir("column.gid")
    #sys.path.append(os.getcwd())

    ##import column_benchmark
    ##Msg = column_benchmark.Run()

    ##print "Running column.py..."
    #Msg = benchmarking.RunBenchmark("column.py", "column_ref.txt")

    #if (Msg == True):
    #	Text += "OK\n"
    #	print "colum example succesful"
    #else:
    #	Text += "FAILED\n"
    #	Text += Msg
    #	Text += "\n\n"
    #	print "colum example FAILED"

    #os.chdir("..")

    ################################################################

    # non-conformant mesh test

    Text += "meshtest: "
    os.chdir("meshtest.gid")
    sys.path.append(os.getcwd())

    print "Running mesh.py..."
    Msg = benchmarking.RunBenchmark("mesh.py","mesh_ref.txt")

    if (Msg == True):
        Text += "OK\n"
        print "non-conformant mesh example successful"
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print "non-conformant mesh example FAILED"

    os.chdir("..")

    ################################################################

    # Add other examples here

    ################################################################
    print "resume of all of the examples for the FSI application :"
    print Text
    return Text

Run()
