// tratamiento de nombres de un archivo

#ifndef _ARCHIVOS_H
#define _ARCHIVOS_H

#include <cstring>

static const int _max_file_len=256;
static const int _max_ext_len=32;

// crea un backup (si el archivo existia)
int backup(const char* arch);

// devuelve extension
// abcd.xyz devuelve &x
// abcd o error devuelve 0
// abcd. devuelve (&.)+1 (un char* que contiene un 0)
char *ext_begin(char* name);
const char *ext_begin(const char* name);
// devuelve (&/) o (&\) (el mayor) +1; o name si no hay
char *name_begin(char* name);
const char *name_begin(const char* name);

// \\maquina\sharedpath\path1\path2\pathn\nombre.ext
//                drive:path1\path2\pathn\nombre.ext
// |               path                  |nombre|ext
// path tiene la ultima barra
class nombre_archivo{
public:
  char
    path[_max_file_len],
    nombre[_max_file_len],
    ext[_max_ext_len];

  nombre_archivo() {path[0]=0;strcpy(nombre,"noname");ext[0]=0;}
  nombre_archivo(const nombre_archivo &n){
    strcpy(path,n.path); strcpy(nombre,n.nombre); strcpy(ext,n.ext);
  }
  nombre_archivo(const char *n) {makefrom(n);}

  const char *makefrom(const char*);

  const char *completo() const;
  char *completo(char *) const;

  // devuelve completo()
  const char *cambiar_ext(const char *);
  const char *cambiar_nombre(const char *);
  const char *agregar_nombre(const char *);
  const char *cambiar_path(const char *);
};

#endif
