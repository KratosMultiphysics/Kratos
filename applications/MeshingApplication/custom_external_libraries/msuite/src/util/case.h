#ifndef _CASE_H
#define _CASE_H

// convierte
char* mayusc(char *s);
char* minusc(char *s);
int strcmp_nocase(const char *c1, const char *c2);

//Terminal <-> ISO 8859-1 o Unicode (acentos y enies)
char *T2U(char *txt);
char *U2T(char *txt);

#endif

