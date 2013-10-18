/* gidpost 2.1 */
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

int GiD_HashInit();
int GiD_HashDone();

GiD_FILE GiD_HashAdd   (void *data);
void    *GiD_HashFind  (GiD_FILE fd);
void     GiD_HashRemove(GiD_FILE fd);

#endif
