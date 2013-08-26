m#!/bin/bash

# generate_makefile.sh v20080402

#common
	export COMMON_SRC_DIR_LIST="util edt malla"
	export RELEASE_COMPILE_OPTIONS="-O3 -funroll-loops -DNDEBUG"
	export DEBUG_COMPILE_OPTIONS="-g2 -D_DEBUG"
#for QT
	export QT_GUI_SRC_DIR_LIST="window window/QT"
	export QT_BIN_NAME="msuite_qt"
#for FLTK
	export FLTK_GUI_SRC_DIR_LIST="window window/FLTK"
	export FLTK_BIN_NAME="msuite_fltk"
#for CLI
	export CLI_GUI_SRC_DIR_LIST="script"
	export CLI_BIN_NAME="msuite_cli"




















#---------------------------------------------------------------
#      check script's arguments
#---------------------------------------------------------------

export SHOW_HELP="0"
if [ "$1" = "" ] || [ "$1" = "qt" ] || [ "$1" = "QT" ]; then
	export CONFIGURATION="QT"
elif [ "$1" = "fltk" ] || [ "$1" = "FLTK" ]; then
	export CONFIGURATION="FLTK"
elif [ "$1" = "cli" ] || [ "$1" = "CLI" ]; then
	export CONFIGURATION="CLI"
else
	export SHOW_HELP="1"
fi

if [ "$2" = "" ] || [ "$2" = "release" ] || [ "$2" = "RELEASE" ]; then
	export CONFIGURATION="$CONFIGURATION release"
elif [ "$2" = "debug" ] || [ "$2" = "DEBUG" ]; then
	export CONFIGURATION="$CONFIGURATION debug"
else
	export SHOW_HELP="1"
fi

if [ "$SHOW_HELP" = "1" ]; then
	echo
	echo "Usar: $0 <cual>"
	echo "   o: $0 <cual> <como>"
	echo
	echo "   <cual> puede ser \"qt\", \"fltk\" o \"cli\""
	echo "   <como> puede ser \"debug\" o \"release\""
	echo "   los valores por defecto son \"qt\" y \"release\""
	echo
	exit
fi


#---------------------------------------------------------------
#    basic checks
#---------------------------------------------------------------

echo
echo "Preparando Configuracion: $CONFIGURATION..."
echo

# check if there's a compiler
if [ "$CXX" = "" ]; then
	export CXX="g++"
fi
echo -n "      g++: "
if $CXX --version >&/dev/null; then
	$CXX --version | head -n 1
else
	echo "ERROR: no hay compilador $CXX (modifique la variable CXX si el ejecutable tiene otro nombre)"
	echo
	exit
fi

# check if there's a make tool
echo -n "      make: "
if make --version >&/dev/null; then
	make --version | head -n 1
else
	echo "ERROR: no se encuentra la herramienta make"
	echo
	exit
fi


#---------------------------------------------------------------
#    QT specific stuf
#---------------------------------------------------------------
if [ "$1" = "" ] || [ "$1" = "qt" ] || [ "$1" = "QT" ]; then

	export LIB_LIST="QtOpenGL QtCore QtGui"
	export GUI_SRC_DIR_LIST="$QT_GUI_SRC_DIR_LIST"
	export BIN_NAME="$QT_BIN_NAME"

	# check if libraries are present
	export MISSING=0
	export LIB_CFLAGS=""
	for A in $LIB_LIST; do
		echo -n "      $A: "
		if ! pkg-config --silence-errors --modversion $A; then
			echo "ERROR: falta la libreria $A"
			export MISSING=1
		fi
	done
	if [ "$MISSING" = "1" ]; then
		echo
		exit
	fi

	# build lib includes part of arguments for g++
	export LIB_CFLAGS=""
	for A in $LIB_LIST; do
		export LIB_CFLAGS="$LIB_CFLAGS $(pkg-config --cflags $A)"
	done
	
	# build lib objects part of arguments for g++
	export LIB_LIBS=""
	for A in $LIB_LIST; do
		export LIB_LIBS="$LIB_LIBS $(pkg-config --libs $A)"
	done


#---------------------------------------------------------------
#    FLTK specific stuf
#---------------------------------------------------------------
elif [ "$1" = "FLTK" ] || [ $1 = "fltk" ]; then

	export GUI_SRC_DIR_LIST="$FLTK_GUI_SRC_DIR_LIST"
	export BIN_NAME="$FLTK_BIN_NAME"

	# check if libraries are present
	export MISSING=0
	export LIB_CFLAGS=""
	echo -n "      Fltk: "
#	if fltk2-config --version 2>/dev/null; then
#		export FLTK_CONFIG="fltk2-config"
#	else
		if fltk-config --version 2>/dev/null; then
			export FLTK_CONFIG="fltk-config"
		else
			echo "ERROR: falta la libreria fltk-1.x"
			exit
		fi
#	fi

	# build lib includes part of arguments for g++
	export LIB_CFLAGS="$($FLTK_CONFIG --use-gl --use-glut --cxxflags)"
	
	# build lib objects part of arguments for g++
	export LIB_LIBS="$($FLTK_CONFIG --use-gl --use-glut --ldflags)"


#---------------------------------------------------------------
#     CLI specific stuf
#---------------------------------------------------------------
elif [ "$1" = "cli" ] || [ $1 = "CLI" ]; then

	export LIB_LIST="$CLI_LIB_LIST"
	export GUI_SRC_DIR_LIST="$CLI_GUI_SRC_DIR_LIST"
	export BIN_NAME="$CLI_BIN_NAME"

fi


#---------------------------------------------------------------
#     RELEASE specific stuff
#---------------------------------------------------------------
if [ "$2" = "" ] || [ "$2" = "release" ] || [ "$2" = "RELEASE" ]; then
	export COMPILE_OPTIONS="$RELEASE_COMPILE_OPTIONS"
#---------------------------------------------------------------

#---------------------------------------------------------------
#    DEBUG specific stuff
#---------------------------------------------------------------
elif [ "$2" = "debug" ] || [ $2 = "DEBUG" ]; then
	export COMPILE_OPTIONS="$DEBUG_COMPILE_OPTIONS"

fi



#---------------------------------------------------------------
#     COMMON stuff
#---------------------------------------------------------------
# build common includes part of arguments for g++
export INCLUDE_PARAMS=""
for A in $COMMON_SRC_DIR_LIST; do
	export INCLUDE_PARAMS="$INCLUDE_PARAMS -I$A"
done

# build gui includes part of arguments for g++
export GUI_INCLUDE_PARAMS=""
for A in $GUI_SRC_DIR_LIST; do
	export GUI_INCLUDE_PARAMS="$GUI_INCLUDE_PARAMS -I$A"
done


#---------------------------------------------------------------
#     init Makefile
#---------------------------------------------------------------
echo > Makefile
echo "#Configuracion $CONFIGURATION" >>Makefile
echo >> Makefile
echo "CXX = $CXX" >>Makefile
echo "LIBS = $LIB_LIBS" >>Makefile
echo "OPTS = $COMPILE_OPTIONS $INCLUDE_PARAMS" >>Makefile
echo "GUI_OPTS = $GUI_INCLUDE_PARAMS $LIB_CFLAGS" >>Makefile
echo >> Makefile
echo "all: $BIN_NAME" >> Makefile
echo >> Makefile
echo 


#---------------------------------------------------------------
#    common objects rules
#---------------------------------------------------------------
export COMMON_OBJS=""
for D in $COMMON_SRC_DIR_LIST; do
	if ! test -d $D; then
		echo ERROR: No se pudo abrir el directorio de fuentes $D
		exit
	fi
	echo "Analizando $D:"
	for A in $D/*.cpp; do
		echo "      $A..."
		export LINE_FOR_MAKEFILE="$($CXX $INCLUDE_PARAMS -MM $A)"
		export COMMON_OBJS="$COMMON_OBJS "$D/$(echo $LINE_FOR_MAKEFILE | cut -d \: -f 1)""
		echo "$D/$LINE_FOR_MAKEFILE" >>Makefile
		echo -e "\\t\${CXX} \${OPTS} $A -c -o \$@" >>Makefile
		echo >>Makefile
	done
	echo
done


#---------------------------------------------------------------
#    gui objects rules
#---------------------------------------------------------------
export GUI_OBJS=""
for D in $GUI_SRC_DIR_LIST; do
	if ! test -d $D; then
		echo ERROR: No se pudo abrir el directorio de fuentes $D
		exit
	fi
	echo "Analizando $D:"
	for A in $D/*.cpp; do
		echo "      $A..."
		export LINE_FOR_MAKEFILE="$($CXX $INCLUDE_PARAMS $GUI_INCLUDE_PARAMS $LIB_CFLAGS -MM $A)"
		export GUI_OBJS="$GUI_OBJS "$D/$(echo $LINE_FOR_MAKEFILE | cut -d \: -f 1)""
		echo "$D/$LINE_FOR_MAKEFILE" >>Makefile
		echo -e "\\t\${CXX} \${OPTS} \${GUI_OPTS} $A -c -o \$@" >>Makefile
		echo >>Makefile
	done
	echo
done


#---------------------------------------------------------------
#     linking and cleaning rules
#---------------------------------------------------------------

echo "Finalizando configuracion..."

echo "$BIN_NAME: $COMMON_OBJS $GUI_OBJS" >>Makefile
echo -e "\\t\${CXX} $COMMON_OBJS $GUI_OBJS \${LIBS} -o \$@" >>Makefile
echo >> Makefile

echo "clean: " >>Makefile
echo -e "\\trm -rf $COMMON_OBJS $GUI_OBJS $BIN_NAME" >>Makefile

rm -rf $GUI_OBJS


#---------------------------------------------------------------

echo 
echo "Ok. Ahora escriba \"make\" y presione Enter."
echo

