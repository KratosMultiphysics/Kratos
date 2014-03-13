#if !defined(DEM_DEFINITIONS_H_INCLUDED )
#define  DEM_DEFINITIONS_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes


#if defined( DEBUG_MACRO )
  #define KRATOS_DEBUG(variable) \
  std::cout << #variable << " : " << variable << std::endl;

#else

  #define KRATOS_DEBUG(variable) \
  {}

#endif


#endif // DEM_DEFINITIONS_H_INCLUDED defined 


