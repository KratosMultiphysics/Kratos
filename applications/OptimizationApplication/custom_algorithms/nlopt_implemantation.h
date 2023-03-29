//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//                   Behzad Tabari, https://github.com/behzadtabari
//

#ifndef NLOPT_IMP_H
#define NLOPT_IMP_H

#ifndef ALGORITHM_BASE_H
#define ALGORITHM_BASE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "algorithm_base.h"
#include <nlopt.hpp> // we need to also always add the CMake configuration of NLOPt 
                     // suited to user's directory
// ==============================================================================

namespace Kratos 
{

 // we need to somehow redifine the analytical function of our problem 
    double myfunc (const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){

    if (!grad.empty()) {
        // grad here must be modified 
    }
    return ; // we need the anlytical function of the problem
}

// constraints of the problem ..
double myconstraint(const std::vector<double> &x, std::vector<double> &grad, void *data){
    
    if (!grad.empty()) {
        
           }
    return ;
}

    class KRATOS_API(OPTIMIZATION_APPLICATION) Nlopt_Implementation : public OptimizationAlgorithm
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(Nlopt_Implementation);
// unlike other algorithm which woud give the opt type in the constructor
// here first we need to get the Algorithm type and give to the constructor 
// I think if we make it user defiend it'll be more user friendly so the user needs to
// know which algorithm would suit 
    

    Nlopt_Implementation (std::string OptName, Model& rModel, Parameters& rOptSettings, int Dim)
    : OptimizationAlgorithm(OptName,OptType,rModel,rOptSettings), OptType_(OptType), Dim_(Dim)  {}

    virtual ~Nlopt_Implementation () {};


    void Initialize() {
        // 
    nlopt::opt opt(nlopt::OptType_, Dim_);

    } 


    private : 
    std :: string OptType_; // it is indeed the type of the algorithm
    int Dim_; // the dimensionality of the should be assesd by user
}; // Classs Nlopt_Implementation
}
