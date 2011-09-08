/*
  class CHARACTER
  ===============
  A minimal class used when passing string arguments from C++
  to FORTRAN 77 (received as FORTRAN 77 CHARACTER strings), and
  subsequently returned back to C++ as properly zero terminated
  strings.

  Method used for zero-termination:
  =================================
  When the CHARACTER destructor is activated the zero-termination
  of the c-string is automatically managed. Zero termination is
  also done each time a string array is subscripted using
  CHARACTER::operator()(size_t index)

  FORTRAN Assumptions:
  ====================
  (1) F77 truncates strings when CHARACTER variable is short
  (2) F77 pads variable with blanks when assigned string is short
  (3) F77 represents a string as a pointer followed by a length
  (4) A string array is stored in contiguous memory

  Author: Carsten A. Arnholm, 20-AUG-1995

  Updates:
      04-MAR-1996 Added features for handling arrays of strings
      16-MAR-1996 Tested array features, explicit padding included
      29-JUL-1996 Tested portability to SGI/Unix, moved decl. of destructor
      04-APR-1997 Using strncpy instead of strcpy in operator=(char* str);
*/

class CHARACTER {
public:
    CHARACTER(char* cstring);
    CHARACTER(char* cstring, const size_t lstr);
   ~CHARACTER();
    CHARACTER operator()(size_t index);
    void  pad(size_t first,size_t howmany=1);
    void  operator=(char* str);
    operator char*();
public:
    char*   rep;  // Actual string
    size_t  len;  // String length
};

inline CHARACTER::CHARACTER(char* cstring)
: rep(cstring), len(strlen(cstring))
{}

inline CHARACTER::CHARACTER(char* cstring, const size_t lstr)
: rep(cstring), len(lstr)
{
   // find position from where to start padding
   size_t slen   = strlen(rep);                // upper limit
   size_t actual = (slen < len)? slen : len;   // actual <= len.
   for(size_t i=actual;i<len;i++) rep[i]=' ';  // Do the padding.
}

inline CHARACTER::~CHARACTER() {
   if(rep[len] == '\0') return;     // catches string constants

   for(int i=len-1;i>=0;i--) {
     if(rep[i] == '\0') break;      // already zero terminated

     if(rep[i] != ' ') {            // non-blank discovered, so
        rep[i+1] = '\0';            // zero-terminate and jump out
        break;
     }
   }
}

inline CHARACTER CHARACTER::operator()(size_t index)
{
    // Construct a temporary CHARACTER object for the array element
    // identified by "index" in order to zero-terminate that element
    size_t pos = index*len;          // start pos of array element
    CHARACTER element(rep+pos,len);  // construct new CHARACTER.
    return element;                  // destructor called here.
}

inline void  CHARACTER::pad(size_t first,size_t howmany)
{

   size_t pos=0,i=0,stop=first+howmany-1;
   for(size_t index=first; index<=stop; index++) {
      pos = index*len;
      size_t slen   = strlen(rep+pos);             // upper limit
      size_t actual = (slen < len)? slen : len;
      for(i=pos+actual;i<pos+len;i++) rep[i]=' ';  // Do the padding.
   }
}

inline void CHARACTER::operator=(char* str)
{
   strncpy(rep,str,len);    // this will copy a zero if str < rep
   rep[len-1] = '\0';       // zero terminate in case strncpy did not
   size_t slen   = strlen(rep);                // upper limit
   size_t actual = (slen < len)? slen : len;   // actual <= len.
   for(size_t i=actual;i<len;i++) rep[i]=' ';  // Do the padding.
}

inline CHARACTER::operator char*()
{
    return rep;
}
