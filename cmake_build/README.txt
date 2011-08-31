This directory is designed for building Kratos using the cmake tool

to compile first of all copy the example file to your own, say local_compile.sh

    cp example_configure.sh.do_not_touch local_compile.sh

and edit local_compile.sh ... to customarize it following your taste.



REMARK: it is highly advisable that the "local_compile.sh" file is located 
in the directory cmake_build. If another directory is chosen one should change
"cmake .." with "cmake /path/to/kratos/root"



then run (assuming the file you copied to is local_compile.sh)

    sh local_compile.sh

hopefully you should enjoy the kratos compilation with CMAKE
