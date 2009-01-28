CC = gcc

CFORTRAN_DIR := ..
SOURCE_DIR := ../source
EXA_DIR := ../examples

INC_FLAGS := -I${SOURCE_DIR} -I${CFORTRAN_DIR}

target := release

ifeq (${target},debug)
OUTPUT_DIR := debug
CFLAGS   := -g2 -Wall ${CFLAGS}
CXXFLAGS := -g2 -Wall ${CXXFLAGS}
LDFLAGS  := -g2
LIB_NAME_SUFFIX := d.a
else
OUTPUT_DIR := release
CFLAGS :=   -O2 -fPIC -funroll-loops ${CFLAGS}
CXXFLAGS := -O2 -fPIC -funroll-loops ${CXXFLAGS}
LDFLAGS := -O2
LIB_NAME_SUFFIX := .a
endif
status := $(shell mkdir -p $(OUTPUT_DIR))

LIB_NAME := ${OUTPUT_DIR}/libgidpost${LIB_NAME_SUFFIX}
EXE_NAMEC = ${OUTPUT_DIR}/tbpc
EXE_NAMEF = ${OUTPUT_DIR}/tbpf

LST_INCLUDES = gidpost.h gidpostInt.h testpost.post.msh zlib.h zlibint1.h zlibint2.h

vpath %cpp ${SOURCE_DIR}
vpath %c ${SOURCE_DIR}
vpath %c ${EXA_DIR}
vpath %f ${EXA_DIR}
vpath %o ${OUTPUT_DIR}

CPP_SRCS = gidpost.cpp gidpostInt.cpp
C_SRCS = zlibint1.c zlibint2.c gidpostfor.c
F_SRC = testpostfor.f

OBJS = ${CPP_SRCS:%.cpp=${OUTPUT_DIR}/%.o} ${C_SRCS:%.c=${OUTPUT_DIR}/%.o}

lib: $(LIB_NAME)

$(LIB_NAME): $(OBJS)
	ar cr $(LIB_NAME) $(OBJS)

testf: $(EXE_NAMEF)

$(EXE_NAMEF): ${OUTPUT_DIR}/testpostfor.o $(LIB_NAME)
	${FC} -g -o $(EXE_NAMEF) ${OUTPUT_DIR}/testpostfor.o $(LIB_NAME) -lstdc++

testc: $(EXE_NAMEC)

$(EXE_NAMEC) : ${OUTPUT_DIR}/testpost.o $(LIB_NAME)
	${CXX} $(LDFLAGS) -o $(EXE_NAMEC) ${OUTPUT_DIR}/testpost.o $(LIB_NAME) -lm
clean:
	rm -rf $(OBJS)

clobber:
	rm $(EXE_NAME) $(LIB_NAME) $(OBJS)

depend: 
	${CC} -MM -Df2cFortran ${CPP_SRCS} ${C_SRCS}

.SUFFIXES : .o .cpp .c .f

${OUTPUT_DIR}/%.o: %.f
	${FC} -g -o $@ -c $<

${OUTPUT_DIR}/%.o: %.c
	${CC} -c $(CFLAGS) -Df2cFortran $(INC_FLAGS) -o $@ $<

${OUTPUT_DIR}/%.o: %.cpp
	${CXX} -c $(CXXFLAGS) $(INC_FLAGS) -o $@ $<

gidpost.o: gidpost.cpp gidpostInt.h zlib.h gidpost.h
gidpostInt.o: gidpostInt.cpp gidpostInt.h zlib.h
zlibint1.o: zlibint1.c zlib.h zlibint1.h
zlibint2.o: zlibint2.c zlibint1.h zlib.h zlibint2.h
gidpostfor.o: gidpostfor.c cfortran/cfortran.h gidpost.h
