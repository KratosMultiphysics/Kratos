// Desglose de nombre de archivo

#include "filename.h"
#include <cstdio>

char *ext_begin(char* n){
  if (!n) return 0;
  n=strrchr(n,'.');
  // para evitar ".." "..\" ".\" "..//" ".//" 
  if (!n||*(n-1)=='.'||*(n+1)=='\\'||*(n+1)=='/') return 0;
  return n+1;
}

const char *ext_begin(const char* n){
  if (!n) return 0;
  n=strrchr(n,'.');
  // para evitar ".." "..\" ".\" "..//" ".//" 
  if (!n||*(n-1)=='.'||*(n+1)=='\\'||*(n+1)=='/') return 0;
  return n+1;
}

char *name_begin(char* name){
  char *slash=strrchr(name,'/'), *backslash=strrchr(name,'\\');
  if (!slash&&!backslash) return name;
  return ((slash>backslash)? slash : backslash) + 1; //mezcla
}
const char *name_begin(const char* name){
  const char *slash=strrchr(name,'/'), *backslash=strrchr(name,'\\');
  if (!slash&&!backslash) return name;
  return ((slash>backslash)? slash : backslash) + 1; //mezcla
}

const char *nombre_archivo::makefrom(const char* input){
  const char *n=name_begin(input), *e=ext_begin(input);
  if (n!=input) memcpy(path,input,n-input);
  if (e){
    memcpy(nombre,n,e-n-2); // sin '.'
    strcpy(ext,e);
  }
  else strcpy(nombre,n);
  return input;
}

char *nombre_archivo::completo(char *output) const{
  sprintf(output,"%s%s.%s",path,nombre,ext);
  return output;
}
const char *nombre_archivo::completo() const{
  char output[_max_file_len];
  return completo(output);
}

const char *nombre_archivo::cambiar_ext(const char *newext){
  if (newext) strcpy(ext,newext);
  else ext[0]=0;
  return completo();
}
const char *nombre_archivo::cambiar_nombre(const char *newnombre){
  if (newnombre) strcpy(nombre,newnombre);// si no lo deja como esta
  return completo();
}
const char *nombre_archivo::agregar_nombre(const char *cacho){
  if (cacho) strcat(nombre,cacho);// si no lo deja como esta
  return completo();
}
const char *nombre_archivo::cambiar_path(const char *newpath){
  if (newpath) strcpy(path,newpath);
  else path[0]=0;
  return completo();
}

// crea un backup (si el archivo existia)
int backup(const char* arch){
  char name[_max_file_len]; strcpy(name, arch); strcat(name,".bak");
  remove(name); // borra el viejo
  return rename(arch,name); // crea el nuevo
}

