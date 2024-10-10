//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

/* FOREGROUND */
#if !defined(_WIN32)
    #define RST  "\x1B[0m"
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

    #define BOLDFONT(x) "\x1B[1m" x RST
    #define FAINTFONT(x) "\x1B[2m" x RST
    #define ITAFONT(x) "\x1B[3m" x RST
    #define UNDL(x) "\x1B[4m" x RST
    #define INVFONT(x) "\x1B[7m" x RST
    #define STRFONT(x) "\x1B[9m" x RST
#else
    #define RST  ""
    #define KRED ""
    #define KGRN ""
    #define KYEL ""
    #define KBLU ""
    #define KMAG ""
    #define KCYN ""
    #define KWHT ""

    #define FRED(x) KRED x RST
    #define FGRN(x) KGRN x RST
    #define FYEL(x) KYEL x RST
    #define FBLU(x) KBLU x RST
    #define FMAG(x) KMAG x RST
    #define FCYN(x) KCYN x RST
    #define FWHT(x) KWHT x RST

    #define BOLDFONT(x) "" x RST
    #define FAINTFONT(x) "" x RST
    #define ITAFONT(x) "" x RST
    #define UNDL(x) "" x RST
    #define INVFONT(x) "" x RST
    #define STRFONT(x) "" x RST
#endif