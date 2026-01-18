#include <stdio.h>
#include <stdlib.h>

#include "libmmgcommon_private.h"

int main() {

  uint16_t tag;
  char *tags_name;

  tag = 0;
  printf("%s\n", MMG5_Get_tagName(tag));

  tag = UINT16_MAX;
  printf("%s\n", MMG5_Get_tagName(tag));

  tag &= ~MG_NUL;
  printf("%s\n", MMG5_Get_tagName(tag));


  return 0;

}
