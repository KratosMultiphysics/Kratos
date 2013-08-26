
#Configuracion CLI release

CXX = g++
LIBS = 
OPTS = -O3 -funroll-loops -DNDEBUG  -Iutil -Iedt -Imalla
GUI_OPTS =  -Iscript 

all: msuite_cli

util/case.o: util/case.cpp util/case.h
	${CXX} ${OPTS} util/case.cpp -c -o $@

util/cone.o: util/cone.cpp util/utiles.h util/punto.h util/cone.h
	${CXX} ${OPTS} util/cone.cpp -c -o $@

util/filename.o: util/filename.cpp util/filename.h
	${CXX} ${OPTS} util/filename.cpp -c -o $@

util/octree.o: util/octree.cpp util/utiles.h util/octree.h util/punto.h \
  util/array1.h util/pline.h malla/nodo.h util/punto.h util/flags.h \
  util/pline.h
	${CXX} ${OPTS} util/octree.cpp -c -o $@

util/pline.o: util/pline.cpp util/pline.h util/utiles.h
	${CXX} ${OPTS} util/pline.cpp -c -o $@

util/punto.o: util/punto.cpp util/punto.h util/utiles.h
	${CXX} ${OPTS} util/punto.cpp -c -o $@

util/tiempo.o: util/tiempo.cpp util/tiempo.h
	${CXX} ${OPTS} util/tiempo.cpp -c -o $@

edt/alpha.o: edt/alpha.cpp util/tiempo.h edt/voronoi.h util/array1.h \
  util/pline.h util/utiles.h util/pline.h malla/malla.h util/utiles.h \
  malla/mensajes.h util/array1.h util/pline.h util/octree.h util/punto.h \
  util/array1.h malla/nodo.h util/punto.h util/flags.h malla/elemento.h \
  util/filename.h malla/nodo.h edt/esfera.h util/utiles.h util/flags.h \
  malla/nodo.h
	${CXX} ${OPTS} edt/alpha.cpp -c -o $@

edt/edt.o: edt/edt.cpp util/tiempo.h malla/malla.h util/utiles.h \
  malla/mensajes.h util/array1.h util/pline.h util/utiles.h util/pline.h \
  util/octree.h util/punto.h util/array1.h malla/nodo.h util/punto.h \
  util/flags.h malla/elemento.h util/filename.h malla/nodo.h
	${CXX} ${OPTS} edt/edt.cpp -c -o $@

edt/esfera.o: edt/esfera.cpp edt/esfera.h util/utiles.h util/array1.h \
  util/pline.h util/utiles.h util/flags.h malla/nodo.h util/punto.h \
  util/flags.h util/pline.h
	${CXX} ${OPTS} edt/esfera.cpp -c -o $@

edt/hydir.o: edt/hydir.cpp util/tiempo.h edt/voronoi.h util/array1.h \
  util/pline.h util/utiles.h util/pline.h malla/malla.h util/utiles.h \
  malla/mensajes.h util/array1.h util/pline.h util/octree.h util/punto.h \
  util/array1.h malla/nodo.h util/punto.h util/flags.h malla/elemento.h \
  util/filename.h malla/nodo.h edt/esfera.h util/utiles.h util/flags.h \
  malla/nodo.h
	${CXX} ${OPTS} edt/hydir.cpp -c -o $@

edt/s_squeeze.o: edt/s_squeeze.cpp util/tiempo.h edt/voronoi.h util/array1.h \
  util/pline.h util/utiles.h util/pline.h malla/malla.h util/utiles.h \
  malla/mensajes.h util/array1.h util/pline.h util/octree.h util/punto.h \
  util/array1.h malla/nodo.h util/punto.h util/flags.h malla/elemento.h \
  util/filename.h malla/nodo.h edt/esfera.h util/utiles.h util/flags.h \
  malla/nodo.h
	${CXX} ${OPTS} edt/s_squeeze.cpp -c -o $@

edt/slivers.o: edt/slivers.cpp util/tiempo.h edt/voronoi.h util/array1.h \
  util/pline.h util/utiles.h util/pline.h malla/malla.h util/utiles.h \
  malla/mensajes.h util/array1.h util/pline.h util/octree.h util/punto.h \
  util/array1.h malla/nodo.h util/punto.h util/flags.h malla/elemento.h \
  util/filename.h malla/nodo.h edt/esfera.h util/utiles.h util/flags.h \
  malla/nodo.h
	${CXX} ${OPTS} edt/slivers.cpp -c -o $@

edt/vor_out.o: edt/vor_out.cpp edt/voronoi.h util/array1.h util/pline.h \
  util/utiles.h util/pline.h malla/malla.h util/utiles.h malla/mensajes.h \
  util/array1.h util/pline.h util/octree.h util/punto.h util/array1.h \
  malla/nodo.h util/punto.h util/flags.h malla/elemento.h util/filename.h \
  malla/nodo.h edt/esfera.h util/utiles.h util/flags.h malla/nodo.h
	${CXX} ${OPTS} edt/vor_out.cpp -c -o $@

edt/voronoi.o: edt/voronoi.cpp util/tiempo.h util/pline.h util/utiles.h \
  edt/voronoi.h util/array1.h util/pline.h malla/malla.h util/utiles.h \
  malla/mensajes.h util/array1.h util/pline.h util/octree.h util/punto.h \
  util/array1.h malla/nodo.h util/punto.h util/flags.h malla/elemento.h \
  util/filename.h malla/nodo.h edt/esfera.h util/utiles.h util/flags.h \
  malla/nodo.h
	${CXX} ${OPTS} edt/voronoi.cpp -c -o $@

malla/corte.o: malla/corte.cpp util/tiempo.h malla/malla.h util/utiles.h \
  malla/mensajes.h util/array1.h util/pline.h util/utiles.h util/pline.h \
  util/octree.h util/punto.h util/array1.h malla/nodo.h util/punto.h \
  util/flags.h malla/elemento.h util/filename.h malla/nodo.h
	${CXX} ${OPTS} malla/corte.cpp -c -o $@

malla/delaunay.o: malla/delaunay.cpp util/tiempo.h malla/malla.h util/utiles.h \
  malla/mensajes.h util/array1.h util/pline.h util/utiles.h util/pline.h \
  util/octree.h util/punto.h util/array1.h malla/nodo.h util/punto.h \
  util/flags.h malla/elemento.h util/filename.h malla/nodo.h \
  edt/voronoi.h util/array1.h util/pline.h malla/malla.h edt/esfera.h \
  util/utiles.h util/flags.h malla/nodo.h
	${CXX} ${OPTS} malla/delaunay.cpp -c -o $@

malla/elemento.o: malla/elemento.cpp util/utiles.h util/case.h malla/elemento.h \
  util/flags.h
	${CXX} ${OPTS} malla/elemento.cpp -c -o $@

malla/fforma.o: malla/fforma.cpp edt/esfera.h util/utiles.h util/array1.h \
  util/pline.h util/utiles.h util/flags.h malla/nodo.h util/punto.h \
  util/flags.h util/pline.h malla/malla.h util/utiles.h malla/mensajes.h \
  util/array1.h util/octree.h util/punto.h util/array1.h malla/nodo.h \
  malla/elemento.h util/filename.h malla/nodo.h
	${CXX} ${OPTS} malla/fforma.cpp -c -o $@

malla/grabamalla.o: malla/grabamalla.cpp util/tiempo.h malla/malla.h \
  util/utiles.h malla/mensajes.h util/array1.h util/pline.h util/utiles.h \
  util/pline.h util/octree.h util/punto.h util/array1.h malla/nodo.h \
  util/punto.h util/flags.h malla/elemento.h util/filename.h malla/nodo.h
	${CXX} ${OPTS} malla/grabamalla.cpp -c -o $@

malla/leemalla.o: malla/leemalla.cpp util/case.h util/tiempo.h malla/malla.h \
  util/utiles.h malla/mensajes.h util/array1.h util/pline.h util/utiles.h \
  util/pline.h util/octree.h util/punto.h util/array1.h malla/nodo.h \
  util/punto.h util/flags.h malla/elemento.h util/filename.h malla/nodo.h
	${CXX} ${OPTS} malla/leemalla.cpp -c -o $@

malla/malla.o: malla/malla.cpp util/tiempo.h malla/malla.h util/utiles.h \
  malla/mensajes.h util/array1.h util/pline.h util/utiles.h util/pline.h \
  util/octree.h util/punto.h util/array1.h malla/nodo.h util/punto.h \
  util/flags.h malla/elemento.h util/filename.h malla/nodo.h util/case.h
	${CXX} ${OPTS} malla/malla.cpp -c -o $@

malla/meshelm.o: malla/meshelm.cpp util/tiempo.h malla/malla.h util/utiles.h \
  malla/mensajes.h util/array1.h util/pline.h util/utiles.h util/pline.h \
  util/octree.h util/punto.h util/array1.h malla/nodo.h util/punto.h \
  util/flags.h malla/elemento.h util/filename.h malla/nodo.h edt/esfera.h \
  util/utiles.h util/array1.h util/flags.h malla/nodo.h
	${CXX} ${OPTS} malla/meshelm.cpp -c -o $@

malla/no_remalla.o: malla/no_remalla.cpp malla/malla.h util/utiles.h \
  malla/mensajes.h util/array1.h util/pline.h util/utiles.h util/pline.h \
  util/octree.h util/punto.h util/array1.h malla/nodo.h util/punto.h \
  util/flags.h malla/elemento.h util/filename.h malla/nodo.h
	${CXX} ${OPTS} malla/no_remalla.cpp -c -o $@

malla/orienta.o: malla/orienta.cpp util/tiempo.h util/cone.h util/punto.h \
  util/utiles.h malla/malla.h util/utiles.h malla/mensajes.h \
  util/array1.h util/pline.h util/pline.h util/octree.h util/array1.h \
  malla/nodo.h util/punto.h util/flags.h malla/elemento.h util/filename.h \
  malla/nodo.h
	${CXX} ${OPTS} malla/orienta.cpp -c -o $@

malla/refina.o: malla/refina.cpp malla/malla.h util/utiles.h malla/mensajes.h \
  util/array1.h util/pline.h util/utiles.h util/pline.h util/octree.h \
  util/punto.h util/array1.h malla/nodo.h util/punto.h util/flags.h \
  malla/elemento.h util/filename.h malla/nodo.h edt/voronoi.h \
  util/array1.h util/pline.h malla/malla.h edt/esfera.h util/utiles.h \
  util/flags.h malla/nodo.h
	${CXX} ${OPTS} malla/refina.cpp -c -o $@

malla/smooth.o: malla/smooth.cpp util/tiempo.h malla/malla.h util/utiles.h \
  malla/mensajes.h util/array1.h util/pline.h util/utiles.h util/pline.h \
  util/octree.h util/punto.h util/array1.h malla/nodo.h util/punto.h \
  util/flags.h malla/elemento.h util/filename.h malla/nodo.h util/cone.h
	${CXX} ${OPTS} malla/smooth.cpp -c -o $@

malla/surface.o: malla/surface.cpp malla/malla.h util/utiles.h malla/mensajes.h \
  util/array1.h util/pline.h util/utiles.h util/pline.h util/octree.h \
  util/punto.h util/array1.h malla/nodo.h util/punto.h util/flags.h \
  malla/elemento.h util/filename.h malla/nodo.h
	${CXX} ${OPTS} malla/surface.cpp -c -o $@

malla/vecinos_y_frontera.o: malla/vecinos_y_frontera.cpp util/tiempo.h \
  malla/malla.h util/utiles.h malla/mensajes.h util/array1.h util/pline.h \
  util/utiles.h util/pline.h util/octree.h util/punto.h util/array1.h \
  malla/nodo.h util/punto.h util/flags.h malla/elemento.h util/filename.h \
  malla/nodo.h
	${CXX} ${OPTS} malla/vecinos_y_frontera.cpp -c -o $@

script/polygen.o: script/polygen.cpp util/tiempo.h malla/malla.h util/utiles.h \
  malla/mensajes.h util/array1.h util/pline.h util/utiles.h util/pline.h \
  util/octree.h util/punto.h util/array1.h malla/nodo.h util/punto.h \
  util/flags.h malla/elemento.h util/filename.h malla/nodo.h
	${CXX} ${OPTS} ${GUI_OPTS} script/polygen.cpp -c -o $@

msuite_cli:  util/case.o util/cone.o util/filename.o util/octree.o util/pline.o util/punto.o util/tiempo.o edt/alpha.o edt/edt.o edt/esfera.o edt/hydir.o edt/s_squeeze.o edt/slivers.o edt/vor_out.o edt/voronoi.o malla/corte.o malla/delaunay.o malla/elemento.o malla/fforma.o malla/grabamalla.o malla/leemalla.o malla/malla.o malla/meshelm.o malla/no_remalla.o malla/orienta.o malla/refina.o malla/smooth.o malla/surface.o malla/vecinos_y_frontera.o  script/polygen.o
	${CXX}  util/case.o util/cone.o util/filename.o util/octree.o util/pline.o util/punto.o util/tiempo.o edt/alpha.o edt/edt.o edt/esfera.o edt/hydir.o edt/s_squeeze.o edt/slivers.o edt/vor_out.o edt/voronoi.o malla/corte.o malla/delaunay.o malla/elemento.o malla/fforma.o malla/grabamalla.o malla/leemalla.o malla/malla.o malla/meshelm.o malla/no_remalla.o malla/orienta.o malla/refina.o malla/smooth.o malla/surface.o malla/vecinos_y_frontera.o  script/polygen.o ${LIBS} -o $@

clean: 
	rm -rf  util/case.o util/cone.o util/filename.o util/octree.o util/pline.o util/punto.o util/tiempo.o edt/alpha.o edt/edt.o edt/esfera.o edt/hydir.o edt/s_squeeze.o edt/slivers.o edt/vor_out.o edt/voronoi.o malla/corte.o malla/delaunay.o malla/elemento.o malla/fforma.o malla/grabamalla.o malla/leemalla.o malla/malla.o malla/meshelm.o malla/no_remalla.o malla/orienta.o malla/refina.o malla/smooth.o malla/surface.o malla/vecinos_y_frontera.o  script/polygen.o msuite_cli
