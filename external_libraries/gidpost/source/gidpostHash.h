/* -*- mode: c++ -*-
 *
 *  gidpostHash.h --
 *
 *    This file declare the interface to the pool of file handlers
 *    which is a mapping from unsigned integers (INT*4) to
 *    CPostFile*. The external API will work with unsigned integers as
 *    file handlers.
 */

#ifndef __GIDPOSTHASH__
#define __GIDPOSTHASH__

#include "gidpost.h"

typedef struct _CPostFile CPostFile; // defined in gipostInt.h

int GiD_HashInit( void );
int GiD_HashDone( void );

GiD_FILE GiD_HashAdd   ( CPostFile *data);
CPostFile *GiD_HashFind  (GiD_FILE fd);
int     GiD_HashRemove(GiD_FILE fd);

#endif
