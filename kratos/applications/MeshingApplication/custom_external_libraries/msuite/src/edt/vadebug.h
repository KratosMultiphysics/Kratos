// _DEBUG -> autocad para voronoi

#ifdef _DEBUG

#include <fstream>

// using namespace std; no!!!, va en quien lo incluye

static void debugout(
  const array1<nodo> &n,
  const array1h<esfera> &s,
  const punto &pt,
  const pline *_sn,  // esferas que se graban
  const pline *_tf,  // caras de la frontera de la cavidad
  bool linea=false  // si graba una linea recorriendo las esferas
  )
{
  int i,j,nv=s[0].NV;
  punto g,pm[4];

  ofstream fs("debug.scr");
  ofstream ft("debug.txt");
  ofstream fd("debug.dat");
  if (!fs.is_open()||!ft.is_open()||!fd.is_open()) return;
  fs << "osmode 0\n";
  double h=0;
  if (_sn&&(*_sn).len) h=pt.distancia(n[s[(*_sn)[0]][1]])/8;
  else if (_tf&&(*_tf).len) h=pt.distancia(n[(*_tf)[0]])/8;
  else h=.1;

  punto::o_separator(",");

  // el punto
  fs << "layer t * s 0 \n";
  fs << "point " << pt << '\n';
  fs << "text j mid " << pt << ' ' << h << " 0 P\n";

  ordlist nn;   
  if (_tf){
    const pline &tf=*_tf;
    for (i=0;i<tf.len;i+=3){
      nn+=tf[i]; nn+=tf[i+1]; if (nv==4) nn+=tf[i+2];
    }
  }
  if (_sn) {
    const pline &sn=*_sn;
    for (i=0;i<sn.len;i++){
      for (j=0;j<nv;j++){
        nn+=s[sn[i]][j];
      }
    }
  }

  fd << nn.len+1+((linea)? (*_sn).len : 0) << " nodes # x y z\n";
  fd << n.len << ' ' << pt << '\n';
  fs << "layer m nn s nn \n";
  for (i=0; i<nn.len; i++){
    fd << nn[i] << ' ' << n[nn[i]] << '\n';
    fs << "text j mid " << n[nn[i]] << ' ' << h << " 0 " << nn[i] << '\n';
  }
  if (linea) for (i=0;i<(*_sn).len;i++) fd << n.len+i+1 << ' ' << s[(*_sn)[i]].c << '\n';

  if (_sn){
    const pline &sn=*_sn;
    // esferas que tienen al nodo o esferas entre el nodo y la cavidad
    fs << "layer m ev s ev \n";
    ft << "esferas del nodo: " << sn.len << '\n';
    fd << sn.len << ((nv==3)? " triangle # f\n" : " tetrahedra # f\n");
    for (i=0;i<sn.len;i++) {
      const esfera &si=s[sn[i]];
      fd << sn[i] << ' ';
      for (j=0;j<nv;j++) fd << si[j] << ' '; fd << si.f << '\n';
      g.zero();
      for (j=0;j<nv;j++) g+=n[si[j]]; g/=nv;
      for (j=0;j<nv;j++) pm[j]=n[si[j]]+(g-n[si[j]]).dir()*h/5;
      ft << i << "; indice: " << sn[i] << "; nodos: ";
      if (nv==3){
        fs
          << "3dface "
          << pm[0] << ' '
          << pm[1] << ' '
          << pm[2] << "  \n";
        ft
          << si[0] << ' '
          << si[1] << ' '
          << si[2] << "; vecinos: "
          << si.vecino[0] << ' '
          << si.vecino[1] << ' '
          << si.vecino[2] << "; volumen: "
          << si.vt << '\n';
      }
      else{
        fs
          << "pface "
          << pm[0] << ' '
          << pm[1] << ' '
          << pm[2] << ' '
          << pm[3] << "  1 2 3  2 3 4  3 4 1  4 1 2  \n";
        ft
          << si[0] << ' '
          << si[1] << ' '
          << si[2] << ' '
          << si[3] << "; vecinos: "
          << si.vecino[0] << ' '
          << si.vecino[1] << ' '
          << si.vecino[2] << ' '
          << si.vecino[3] << "; volumen: "
          << si.vt << '\n';
      }
    }
  }

  if (_tf){
    const pline &tf=*_tf;
    // caras
    fs << "layer m tf s tf \n";
    ft << "caras de frontera de la caviad: " << tf.len/3 << endl;
    fd << tf.len/3;
    if (nv==3) fd << " segment #\n";
    else fd << " triangle #\n";
    for (i=0;i<tf.len;i+=3) {
      ft << i << " nodos: ";
      fd << i << ' ';
      if (nv==3){
        ft
          << tf[i] << ' '
          << tf[i+1] << '\n';
        fd
          << tf[i] << ' '
          << tf[i+1] << '\n';
        fs
          << "line "
          << n[tf[i]] << ' '
          << n[tf[i+1]] << " \n";
      }
      else{
        ft
          << tf[i] << ' '
          << tf[i+1] << ' '
          << tf[i+2] << '\n';
        fd
          << tf[i] << ' '
          << tf[i+1] << ' '
          << tf[i+2] << '\n';
        fs
          << "3dface "
          << n[tf[i]] << ' '
          << n[tf[i+1]] << ' '
          << n[tf[i+2]] << "  \n";
      }
    }
  }

  if (linea){
    const pline &sn=*_sn;
    // linea (recorrido por los centros de esferas)
    fs << "layer m linea s linea \nline\n";
    fd << sn.len-1 << " segment #\n";
    fs << s[sn[0]].c << '\n';
    fd << "0 " << n.len+i << ' ';
    for (i=1;i<sn.len-1;i++) {
      fs << s[sn[i]].c << '\n';
      fd << n.len+i+1 << '\n' << i << ' ' << n.len+i+1 << ' ';
    }    
    fs << s[sn[i]].c << "\n\n";
    fd << n.len+i+1 << '\n';
  }
  fs << "layer s 0 \n";
  fs.close();
  ft.close();
  fd.close();
  punto::io_reset();
}

#endif // _DEBUG
