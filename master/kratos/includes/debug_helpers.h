//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

/* System includes */
#include <sstream>
#include <iomanip>

/* External includes */

/* Project includes */

// Print vector values with precision
#define KRATOS_WATCH_VECTOR_WITH_PRECISION(variable, precision)                           \
    {                                                                                     \
        std::stringstream values;                                                         \
        values << std::scientific << std::setprecision(precision);                        \
        for (std::size_t i = 0; i < variable.size(); i++)                                 \
        {                                                                                 \
            values << #variable << "[" << i << "] = " << variable[i] << ";" << std::endl; \
        }                                                                                 \
        std::cout << #variable << " : " << std::endl                                      \
                  << values.str() << std::endl;                                           \
    }

// Print matrix values with precision
#define KRATOS_WATCH_MATRIX_WITH_PRECISION(variable, precision)         \
    {                                                                   \
        std::stringstream values;                                       \
        values << std::scientific << std::setprecision(precision);      \
        for (std::size_t i = 0; i < variable.size1(); i++)              \
        {                                                               \
            for (std::size_t j = 0; j < variable.size2(); j++)          \
            {                                                           \
                values << #variable << "(" << i << ", " << j            \
                       << ") = " << variable(i, j) << ";" << std::endl; \
            }                                                           \
        }                                                               \
        std::cout << #variable << " : " << std::endl                    \
                  << values.str() << std::endl;                         \
    }

// Cout a message with appended line and function information.
#define KRATOS_WATCH_LINE(...) std::cout << __VA_ARGS__ << KRATOS_CODE_LOCATION << std::endl
