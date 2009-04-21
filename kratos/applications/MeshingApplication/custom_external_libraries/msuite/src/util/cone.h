#ifndef _CONE_H_
#define _CONE_H_

#include "punto.h"

// cone returns versor vc minimizing maximum angle [0,180] with the given versors
// cc is the maximum cosine (aperture)
// The axis must give positive scalar with mean dir which must be the first point in the list

// puntos iguales joden cuando ninguno entra en el cono del otro
// por ser unitarios, ERRADM es absoluto (asi que no hay problemas)

bool cone(const punto* const l, int nl, punto &vc, double &cc);

#endif
