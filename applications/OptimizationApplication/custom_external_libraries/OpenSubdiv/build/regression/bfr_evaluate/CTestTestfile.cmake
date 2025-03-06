# CMake generated Testfile for 
# Source directory: /home/basti/KratosMultiphysics/OpenSubdiv/regression/bfr_evaluate
# Build directory: /home/basti/KratosMultiphysics/OpenSubdiv/build/regression/bfr_evaluate
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(bfr_evaluate_pos "/home/basti/KratosMultiphysics/OpenSubdiv/build/bin/bfr_evaluate" "-all" "-silent" "-l" "3" "-pass" "2" "-d1")
set_tests_properties(bfr_evaluate_pos PROPERTIES  _BACKTRACE_TRIPLES "/home/basti/KratosMultiphysics/OpenSubdiv/regression/bfr_evaluate/CMakeLists.txt;48;add_test;/home/basti/KratosMultiphysics/OpenSubdiv/regression/bfr_evaluate/CMakeLists.txt;0;")
add_test(bfr_evaluate_uv1 "/home/basti/KratosMultiphysics/OpenSubdiv/build/bin/bfr_evaluate" "-all" "-silent" "-l" "3" "-pass" "0" "-skippos" "-uv" "-uvint" "1")
set_tests_properties(bfr_evaluate_uv1 PROPERTIES  _BACKTRACE_TRIPLES "/home/basti/KratosMultiphysics/OpenSubdiv/regression/bfr_evaluate/CMakeLists.txt;50;add_test;/home/basti/KratosMultiphysics/OpenSubdiv/regression/bfr_evaluate/CMakeLists.txt;0;")
add_test(bfr_evaluate_uv5 "/home/basti/KratosMultiphysics/OpenSubdiv/build/bin/bfr_evaluate" "-all" "-silent" "-l" "3" "-pass" "0" "-skippos" "-uv" "-uvint" "5")
set_tests_properties(bfr_evaluate_uv5 PROPERTIES  _BACKTRACE_TRIPLES "/home/basti/KratosMultiphysics/OpenSubdiv/regression/bfr_evaluate/CMakeLists.txt;52;add_test;/home/basti/KratosMultiphysics/OpenSubdiv/regression/bfr_evaluate/CMakeLists.txt;0;")
