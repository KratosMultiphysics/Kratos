#include "voronoi.h"

#include <ostream>
#include <fstream>
#include <iomanip> // setw
#include <cstring> // strcat strcpy
#include <cstdio> // sprintf remove rename

using namespace std;

//===============================
//    #defines utiles
//===============================
#define _oi(x) setw(10) << x
#define _ol(x) setw(15) << x
#define _or(x) setprecision(16) << setw(25) << x
#define _of(x) setprecision(3) << setw(10) << x

#define _OF(trail)\
  char fn[_max_file_len];strcpy(fn,filename);strcat(fn,trail);\
  back(fn); ofstream f(fn); if (!f) return false;

#define _close f.close(); return true;

#define _makename(trail) \
  char filename[_max_file_len]; \
  if (arch&&arch[0]) strcpy(filename,arch); \
  else {strcpy(filename,m->nombre);strcat(filename,trail);}

//==================================
// crea un backup (si el archivo existia)
static int back(const char* arch)
{
  char name[_max_file_len]; strcpy(name, arch); strcat(name,".bak");
  remove (name); // borra el viejo
  return rename(arch,name); // crea el nuevo
}
//==================================

bool voronoi::graba_dat(const char *arch) {
  _makename("");
  _OF(".sph");
  int i,j;
//  const array1<nodo> &n=m->n;

  punto::io_z(NV==4); punto::o_precision(16); 
  punto::o_width(25); punto::o_separator(" ");
/*
  if (dir){
    f << dir.len << " Normals\n";
    for (i=0;i<dir.len;i++){
      f << _oi(i)
        << dir[i]
        << endl;
    }
  }
*/
  double vt=0;
  for (i=0;i<s.len;i++) if (!s[i].f.es(e_borrado)) vt+=s[i].vt;

  f << s.len-s.borrado.len 
    << " Spheres: # f nodes r c vol\t(Volume: "
    << vt 
    << ")\n";

  for (i=0;i<s.len;i++){
    if (s[i].f.es(e_borrado)) continue;
    f << _oi(i);
    f << _ol(s[i].f);
    for (j=0;j<NV;j++) f << _oi(s[i].n[j]);
    f << _or(s[i].r) << s[i].c << _or(s[i].vt) << endl;
  }
  punto::io_reset();
  _close;
}

