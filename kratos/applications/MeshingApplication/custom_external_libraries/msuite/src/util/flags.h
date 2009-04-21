// los flags son int para usar enums y por portabilidad

#ifndef _FLAGS_H
#define _FLAGS_H

class flagtype{
public:
  int f;

  flagtype(): f(0){};
  flagtype(int flag){f=flag;}

  inline operator int() const {return f;}

  const flagtype operator &(flagtype &o) const {return flagtype(f&o.f);}
  const flagtype operator |(flagtype &o) const {return flagtype(f|o.f);}
  
  inline void reset() {f=0;} 
  inline void ini() {f=0;} 

  inline bool operator==(int flag){return (f==flag);}
  inline int operator=(int flag){return (f=flag);}
  inline int operator|=(int flag){return (f|=flag);}
  inline int operator&=(int flag){return (f&=flag);}

  // (para varios, n entra oreado)
  inline void set(int n) {f|=n;} 
  inline void reset(int n) {f&=(~n);} 
  inline void mask(int n) {f&=n;} 
  inline bool es_todos(int n) const {return ((f&n)==n);}
  inline bool es_alguno(int n) const {return ((f&n)!=0);}
  inline bool es(int n) const {return ((f&n)!=0);} // alguno
  inline bool noes(int n) const {return ((f&n)==0);} // no es ninguno
  inline bool noes_ninguno(int n) const {return ((f&n)==0);} // no es ninguno
  inline bool noes_alguno(int n) const {return ((f&n)!=n);} // no es alguno
};

#endif

