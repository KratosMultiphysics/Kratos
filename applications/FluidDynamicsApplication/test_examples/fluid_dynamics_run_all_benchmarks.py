from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import sys

kratos_benchmarking_path = '../../../benchmarking'
sys.path.append(kratos_benchmarking_path)
import benchmarking


def Run():
    Msg = ""
    Text = "== Fluid Dynamics ==========\n"

    #
    # VMS2D element test

    Text += "VMS2D element test: "
    os.chdir("vms2d_test")
    sys.path.append(os.getcwd())

    print("running the vms2d_test benchmark test...")
    successful,Msg = benchmarking.RunBenchmark("script_elemtest.py", "vms2d_test_ref.txt")

    if(successful==True):
        Text += "OK\n"
        print("VMS2D element test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("VMS2D element test FAILED")

    os.chdir("..")
    #
    # parabolic flow in a trapezoidal domain: test VMS2D OSS implementation

    Text += "Parabolic flow OSS test: "

    os.chdir("oss_trapezoid")
    sys.path.append(os.getcwd())

    import trapezoid_benchmark
    success, Msg = trapezoid_benchmark.Run()

    if(successful==True):
        Text += "OK\n"
        print("oss_trapezoid test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("oss_trapezoid test FAILED")

    os.chdir("..")
    #
    # slip condition and wall law (MonolithicWallCondition2D + ResidualBasedVelocityBossakSchemeTurbulent)

    Text += "Slip condition and wall law test: "

    os.chdir("slip_test")
    sys.path.append(os.getcwd())

    import slip_test_benchmark
    success, Msg = slip_test_benchmark.Run()

    if(successful==True):
        Text += "OK\n"
        print("slip condition test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("slip condition test FAILED")

    os.chdir("..")
    #
    # fractional step 2D cavity

    Text += "Fractional step element test: 2D cavity flow: "

    os.chdir("fs_cavity")
    sys.path.append(os.getcwd())

    import fs_cavity_verify_benchmark
    success, Msg = fs_cavity_verify_benchmark.Run()

    if(successful==True):
        Text += "OK\n"
        print("Fractional step 2D cavity test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("Fractional step 2D cavity test FAILED")

    os.chdir("..")
    #
    # test of cylinder with slip BC, wall law and outlet close to the cylinder

    Text += "test of cylinder with slip BC, wall law and outlet close to the cylinder: "

    os.chdir("cyl_bench")
    sys.path.append(os.getcwd())

    import cyl_verify_benchmark
    success, Msg = cyl_verify_benchmark.Run()

    if(successful==True):
        Text += "OK\n"
        print("cylinder with slip BC, wall law and outlet test successful")
    else:
        Text += "FAILED\n"
        Text += Msg
        Text += "\n\n"
        print("cylinder with slip BC, wall law and outlet test FAILED")

    os.chdir("..")

    # Add other examples here

    #
    print("resume of all of the examples for the Fluid Dynamics application :")
    print(Text)
    return Text

if __name__ == '__main__':
    Run()
