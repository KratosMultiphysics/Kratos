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
#ifdef _WIN32
#include <Windows.h>
#endif
// External includes

// Project includes

/* FOREGROUND */
#ifdef _WIN32
    #define PREPARE_STL_OUTPUT SetConsoleMode(GetStdHandle(STD_OUTPUT_HANDLE), ENABLE_VIRTUAL_TERMINAL_PROCESSING);
#else
    #define PREPARE_STL_OUTPUT void();
#endif
    
#define RST   "\x1B[0m"
#define BOLD  "\x1B[1m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST

#define BOLDFONT(x) BOLD x RST
#define UNDL(x) "\x1B[4m" x RST

#endif  /* COLOR_UTILITIES_H_INCLUDED */
