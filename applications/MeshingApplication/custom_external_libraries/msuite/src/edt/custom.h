///////////////////////////////////////////////////////////////////////////////////////////
//   parametros "tocables"
///////////////////////////////////////////////////////////////////////////////////////////

static bool my_defaults(){  
  perturba=false;//con array regular, al rehacer aparecen esferas no-huecas
  reponer=true; // reponer posicion de nodos (true)
  parcial=true; // pone solo si esta a mas de dnmin*h (false) (dnmin en voronoi.cpp)
  test_close_to_boundary=false; // testear y en caso positivo no poner el nodo

  _pb=1e-8;
  _ph=1e-3;

  rmax=0.8; //si r>rmax*h pone un punto en el centro de la esfera .7
  dmax=100; //si dist entre nodos frontera > dmax*h => agrega pto medio 1.5
  dnmin=.5;// d<dnmin*h => no agrega el punto (parcial) 0.5
  dfmin=.5;// d_frontera<dfmin*h => no agrega el punto (test_close_to_boundary) 0.5

  return true;
}

static bool temp=my_defaults();

///////////////////////////////////////////////////////////////////////////////////////////

static pline IndiceNodoQueda; //para nodos muy cercanos en delaunay
static array1<punto> NodoAcombinar;
