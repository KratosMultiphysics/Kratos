# Set EMPIRE environment 
# Please change that
export EMPIRE_BASE_DIR=/home/jwolf/kratos/applications/empire_application/custom_external_libraries/empire
export MPICH2_BASE_DIR=/home/jwolf/kratos/applications/empire_application/custom_external_libraries/empire/mpich2
######################################################################################
# Set PATH ENV variables for INCLUDE, LINKING and EXECUTION 
	export EMPIRE_API_LIB_NAME=EMPIRE_API
	export EMPIRE_LD_LIBRARY_PATH=$EMPIRE_BASE_DIR/lib
	export EMPIRE_API_INC_ON_MACHINE=$EMPIRE_BASE_DIR/EMPIRE_API/src/include/
	export EMPIRE_API_LIB_ON_MACHINE=$EMPIRE_LD_LIBRARY_PATH/lib$EMPIRE_API_LIB_NAME.a	
	export EMPIRE_API_LIBSO_ON_MACHINE=$EMPIRE_LD_LIBRARY_PATH/lib$EMPIRE_API_LIB_NAME.so
	export LD_LIBRARY_PATH=$EMPIRE_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
	export EMP_PATH=$EMPIRE_BASE_DIR/bin:$EMP_PATH
	export PATH=$EMP_PATH:$PATH
######################################################################################
# Set PATH ENV variables for MPICH2	
	export LD_LIBRARY_PATH=$MPICH2_BASE_DIR/lib:$LD_LIBRARY_PATH
	export PATH=$MPICH2_BASE_DIR/bin:$PATH
