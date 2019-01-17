##
#
# @file old/Makefile.in
#
# @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                      Univ. Bordeaux. All rights reserved.
#
# Makefile template for installed sources.
#
# @version 6.0.1
# @author Mathieu Faverge
# @date 2018-07-16
#
##
PASTIX_DIR=/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/PasTiX

PKG_CONFIG_PATH:=${PKG_CONFIG_PATH}:/home/josep_maria/Kratos/applications/SolversApplication/external_libraries/pastix-lib/pastix/PasTiX/lib/pkgconfig

EXTRALIBS= /usr/lib/x86_64-linux-gnu/liblapacke.so /usr/lib/x86_64-linux-gnu/libopenblas.so /usr/lib/x86_64-linux-gnu/libm.so /home/josep_maria/Kratos/libs/libscotch.a /usr/lib/x86_64-linux-gnu/libscotcherrexit.so /usr/lib/x86_64-linux-gnu/libpthread.so /usr/lib/x86_64-linux-gnu/libz.so /usr/lib/x86_64-linux-gnu/librt.so

CFLAGS=`PKG_CONFIG_PATH=${PKG_CONFIG_PATH} pkg-config --cflags pastix` -Wall -O2
LDFLAGS=`PKG_CONFIG_PATH=${PKG_CONFIG_PATH} pkg-config --libs pastix` ${EXTRALIBS}

EXAMPLES=       simple	        \
		step-by-step

default: $(EXAMPLES)

% : %.c
	$(CC) -o $@ $< ${CFLAGS} ${LDFLAGS}

clean:
	rm -f ${EXAMPLES}
