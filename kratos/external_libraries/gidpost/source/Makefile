
# -*- mode: Makefile;-*-

# to compile debug version, do: make DEBUG=yes

HOST = $(shell hostname)
OS = $(shell if [ -e /usr/include/linux ]; then echo linux; else echo mac; fi)
ifeq ($(WINDOWS),yes)
OS = windows
endif

ifeq ($(HOST),rrg7.local)
  LIBSDIR=/Users/ramsan/gidproject/libs
else ifeq ($(HOST),hoschi)
  LIBSDIR=/Users/miguel
else
  LIBSDIR=
endif

OBJDIR = $(if $(filter yes,$(DEBUG)),debug,release)

ifeq ($(USER),ramsan)
HDF5_DIRECTORY = /home/ramsan/gidproject/gid/hdf5-1.8.3/hdf5
else
TCL_DIRECTORY = /opt/hdf5
endif

CC = gcc
CPPFLAGS=-DUSE_TCL_STUBS -DTCLVERSION
STRIP=strip

INCLUDE_DIRECTORIES = $(HDF5_DIRECTORY)/include
LDFLAGS= -static-libgcc

ifeq ($(OS),mac)
  INCLUDE_DIRECTORIES +=  /opt/local/include
  CPPFLAGS += -fast -arch i386
  LDFLAGS += -dynamiclib -arch i386
  LIB_DIRECTORIES = $(LIBSDIR)/lib /opt/local/lib
  LIBEXT = dylib
  LD = g++
else
  CPPFLAGS += -fPIC
  LDFLAGS += -shared -fPIC
  LIB_DIRECTORIES = /usr/local/ActiveTcl-8.5/lib
  LIBEXT = so
  LD = gcc
endif

CPPFLAGS += $(addprefix -I ,$(INCLUDE_DIRECTORIES))
LDFLAGS += $(addprefix -L,$(LIB_DIRECTORIES)) -lz

LIBS+ = $(HDF5_DIRECTORY)/lib/libhdf5_hl.a $(HDF5_DIRECTORY)/lib/libhdf5.a

EXE = $(OBJDIR)/gidpost.$(LIBEXT)
EXE_INSTALL = 

ifeq ($(OBJDIR),debug)
CPPFLAGS += -g -D_DEBUG -DDEBUG -Wall
LDFLAGS  += -g
endif
ifeq ($(OBJDIR),release)
CPPFLAGS += -O3
LDFLAGS  += -O3
endif

ifeq ($(OS),windows)
  CC = i586-mingw32msvc-gcc
  LD = i586-mingw32msvc-g++
  STRIP = i586-mingw32msvc-strip
  CPPFLAGS += -D_WINDOWS -DNO_SIGACTION -DWIN32 -D_WIN32 -I ~/code/mingw/include
  OBJDIR = release_mingw
  LDFLAGS = -m32 -shared
  LIBS0 =  kernel32 user32 gdi32 win32k ntdll ntoskrnl \
	advapi32 uuid mpr netapi32 winmm ws2_32 rpcrt4 \
	ole32 shell32 comdlg32 odbc32 odbccp32
  LIBS = -L ~/code/mingw/lib $(addprefix -l,$(LIBS0))
  EXE = $(OBJDIR)/reporter.dll
  EXE_INSTALL =
endif

SRCS = \
    gidpost.c     gidpostHash.c  gidpostInt.c  hdf5c.c    recycle.c \
    gidpostHDF5.c  hashtab.c lookupa.c
#    zlibint2.c  zlibint1.c

OBJS = $(addprefix $(OBJDIR)/,$(SRCS:.c=.o))

all: compile $(EXE_INSTALL)
compile: $(OBJDIR) $(EXE)

-include $(addprefix $(OBJDIR)/,$(SRCS:.c=.d))

$(EXE): $(OBJS)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJS) $(LIBS)
ifneq ($(OBJDIR),debug)
	 $(STRIP) $(EXE)
endif

$(OBJDIR):
	mkdir $(OBJDIR)

$(OBJDIR)/%.o: %.c
	$(CC) -MMD -c $(CPPFLAGS) $< -o $@

$(EXE_INSTALL): $(EXE)
ifneq ($(OBJDIR),debug)
	cp $(EXE) $(EXE_INSTALL)
endif

clean:
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d $(EXE)

copy: compile
	cp $(EXE) $(EXE_INSTALL)

