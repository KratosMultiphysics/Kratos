//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                    
//

#ifndef COLOR_UTILITIES_H_INCLUDED
#define COLOR_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"

/* FOREGROUND */    
#ifdef KRATOS_COMPILED_IN_WINDOWS
    #define BASE_STRING_ANSI   "\033"
#else
    #define BASE_STRING_ANSI   "\x1B"
#endif

#define RST   BASE_STRING_ANSI "[0m"
#define BOLD_FONT  BASE_STRING_ANSI "[1m"
#define KRED  BASE_STRING_ANSI "[31m"
#define KGRN  BASE_STRING_ANSI "[32m"
#define KYEL  BASE_STRING_ANSI "[33m"
#define KBLU  BASE_STRING_ANSI "[34m"
#define KMAG  BASE_STRING_ANSI "[35m"
#define KCYN  BASE_STRING_ANSI "[36m"
#define KWHT  BASE_STRING_ANSI "[37m"
#define UNDERLINE_FONT  BASE_STRING_ANSI "[4m"

#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST

#define BOLDFONT(x) BOLD_FONT x RST
#define UNDL(x) UNDERLINE_FONT x RST

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
/**
 * @namespace ColorUtilities
 * @ingroup KratosCore
 * @brief This namespace includes several utilities necessaries for use ANSI colors
 * @author Vicente Mataix Ferrandiz
 */
namespace ColorUtilities
{
    /**
     * @brief This method inits the standard STD output in each OS
     */
    void KRATOS_API(KRATOS_CORE) InitializeSTDOutput();

}; // namespace ColorUtilities
}  // namespace Kratos

#endif  /* COLOR_UTILITIES_H_INCLUDED */
