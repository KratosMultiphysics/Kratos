SHELL = /bin/sh
RANLIB = ranlib
cur_dir := $(shell pwd)
KRATOS_PATH := ../..
EXTLIB_PATH := ${KRATOS_PATH}/external_libraries
APPLICATION_SRC_PATH := ${cur_dir}
APPLICATION_INC :=  ${KRATOS_PATH}/kratos

LIBAPPLICATION_PREFIX := libapplication

EXE_NAME := application

MAIN_PATH := ${cur_dir}

vpath %cpp ${MAIN_PATH} 

# set the external dependencies
#

# BOOST include path

# ifeq (${BOOST_INC},)
BOOST_INC :=  ${EXTLIB_PATH}/boost
# endif

# gidpost path

GIDPOST_INC := ${EXTLIB_PATH}

LIB_GIDPOST := libgidpost.a
LIB_GIDPOST_PATH := ${EXTLIB_PATH}/gidpost/unix/release
EXT_LIBS :=  -L${LIB_GIDPOST_PATH} -lgidpost
CXXFLAGS := -I${GIDPOST_INC} -I${BOOST_INC} -I${APPLICATION_INC}

vpath %.cpp ${APPLICATION_SRC_PATH} ${APPLICATION_INC}/sources

target := release                                                               
ifeq (${target},debug)
OUTPUT_DIR := debug
LIB_KRATOS := -L${KRATOS_PATH}/bin/debug -lkratos
CFLAGS   := -g ${CFLAGS}
CXXFLAGS := -g ${CXXFLAGS}
LDFLAGS  := -g
LIBAPPLICATION := ${LIBAPPLICATION_PREFIX}d.a
else
OUTPUT_DIR := release
LIB_KRATOS := -L${KRATOS_PATH}/bin/release -lkratos
CFLAGS :=   -O3 -funroll-loops ${CFLAGS}
CXXFLAGS := -O3 -funroll-loops ${CXXFLAGS}
LDFLAGS := -O3
LIBAPPLICATION := ${LIBAPPLICATION_PREFIX}.a
endif
status := $(shell mkdir -p $(OUTPUT_DIR))        

libapplication_src = \
main.cpp

libapplication_objs = ${libapplication_src:%.cpp=$(OUTPUT_DIR)/%.o}

# DEP_FILES = ${libapplication_src:%.cpp=$.deps/%.P}

all: ${OUTPUT_DIR}/${EXE_NAME}

${OUTPUT_DIR}/${EXE_NAME}: ${OUTPUT_DIR}/main.o ${OUTPUT_DIR}/${LIBAPPLICATION}
	${CXX} -o $@ $< ${OUTPUT_DIR}/${LIBAPPLICATION} ${EXT_LIBS} ${LIB_KRATOS}

${OUTPUT_DIR}/%.o : %.cpp
	$(CXX) ${CXXFLAGS} -c $< -o $@

# DEPS_MAGIC := $(shell mkdir .deps > /dev/null 2>&1 || :)

# -include $(DEP_FILES)

lib: ${OUTPUT_DIR}/${LIBAPPLICATION}

# .SUFFIXES:
# .SUFFIXES: .S .c .cpp .lo .o .obj .s



${OUTPUT_DIR}/${LIBAPPLICATION}: ${libapplication_objs}
	-rm -f $@
	$(AR) cru $@ ${libapplication_objs}
	$(RANLIB) $@

