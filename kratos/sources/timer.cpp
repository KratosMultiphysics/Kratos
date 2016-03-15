//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//











#include "utilities/timer.h"


namespace Kratos
{
/// Default constructor.
Timer::Timer() {}

Timer::ContainerType Timer::msTimeTable;
std::ofstream Timer::msOutputFile;

bool Timer::msPrintOnScreen = false;

#ifndef _OPENMP
double Timer::msGlobalStart = std::clock()/static_cast<double>(CLOCKS_PER_SEC);
#else
double Timer::msGlobalStart = omp_get_wtime();
#endif




}

