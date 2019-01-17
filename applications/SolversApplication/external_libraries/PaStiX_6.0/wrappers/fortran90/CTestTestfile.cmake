# CMake generated Testfile for 
# Source directory: /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90
# Build directory: /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(fortran_fsimple "./fsimple")
add_test(fortran_flaplacian "./flaplacian")
add_test(fortran_fstep-by-step "./fstep-by-step")
add_test(fortran_fmultilap_seq "/bin/bash" "-c" "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90/fmultilap < /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90/examples/test_seq.in")
add_test(fortran_fmultilap_mt "/bin/bash" "-c" "/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90/fmultilap < /home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/wrappers/fortran90/examples/test_mt.in")
