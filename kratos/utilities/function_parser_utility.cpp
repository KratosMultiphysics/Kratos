//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes

// External includes
#include "exprtk/exprtk.hpp"

// Project includes
#include "includes/global_variables.h"
#include "utilities/function_parser_utility.h"
#include "utilities/string_utilities.h"

namespace Kratos
{

GenericFunctionUtility::GenericFunctionUtility(
    const std::string& rFunctionBody,
    Parameters LocalSystem
    ) : mFunctionBody(StringUtilities::ReplaceAllSubstrings(rFunctionBody, "**", "^")) // Correcting function from python style to exprtk
{
    // Defining namespace
    mNameSpace.insert(std::pair<std::string, double>("x", 0.0));
    mNameSpace.insert(std::pair<std::string, double>("y", 0.0));
    mNameSpace.insert(std::pair<std::string, double>("z", 0.0));
    mNameSpace.insert(std::pair<std::string, double>("t", 0.0));
    mNameSpace.insert(std::pair<std::string, double>("X", 0.0));
    mNameSpace.insert(std::pair<std::string, double>("Y", 0.0));
    mNameSpace.insert(std::pair<std::string, double>("Z", 0.0));
    mNameSpace.insert(std::pair<std::string, double>("pi", Globals::Pi));

    // Initialize exprtk classes
    InitializeExprtk();

    // Compiling expression
    mpParser = new exprtk::parser<double>();
    mpParser->compile(mFunctionBody, *mpExpression);

    // Here get the local system if it is provided
    if(LocalSystem.Has("origin")) {
        mUseLocalSystem = true;

        for(IndexType i = 0; i<3; ++i) {
            mCenterCoordinates[i] = LocalSystem["origin"][i].GetDouble();
            for(IndexType j = 0; j<3; ++j) {
                mRotationMatrix(i, j) = LocalSystem["axes"][i][j].GetDouble();
            }
        }
    } else {
        mUseLocalSystem = false;
    }

    // Check if it depends on space
    if (mFunctionBody.find(std::string("x")) == std::string::npos &&
        mFunctionBody.find(std::string("y")) == std::string::npos &&
        mFunctionBody.find(std::string("z")) == std::string::npos &&
        mFunctionBody.find(std::string("X")) == std::string::npos &&
        mFunctionBody.find(std::string("Y")) == std::string::npos &&
        mFunctionBody.find(std::string("Z")) == std::string::npos) {
        mDependsOnSpace = false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

GenericFunctionUtility::GenericFunctionUtility(GenericFunctionUtility const& rOther)
    : mNameSpace(rOther.mNameSpace),
      mFunctionBody(rOther.mFunctionBody),
      mDependsOnSpace(rOther.mDependsOnSpace),
      mUseLocalSystem(rOther.mUseLocalSystem),
      mRotationMatrix(rOther.mRotationMatrix),
      mCenterCoordinates(rOther.mCenterCoordinates)
{
    // Initialize exprtk classes
    InitializeExprtk();
}

/***********************************************************************************/
/***********************************************************************************/

GenericFunctionUtility::~GenericFunctionUtility()
{
//     delete [] mpSymbolTable; /// The symbol table of exprtk
//     delete [] mpExpression;  /// The expression of exprtk
//     delete [] mpParser;      /// The parser of exprtk
}

/***********************************************************************************/
/***********************************************************************************/

std::string GenericFunctionUtility::FunctionBody()
{
    return StringUtilities::ReplaceAllSubstrings(mFunctionBody, "^", "**");
}

/***********************************************************************************/
/***********************************************************************************/

double GenericFunctionUtility::RotateAndCallFunction(
    const double x,
    const double y,
    const double z,
    const double t,
    const double X,
    const double Y,
    const double Z
    )
{
    array_1d<double,3> xglobal;
    xglobal[0] = x;
    xglobal[1] = y;
    xglobal[2] = z;
    array_1d<double,3> xlocal = prod(mRotationMatrix, (xglobal - mCenterCoordinates) );
    array_1d<double,3> xglobal_initial;
    xglobal_initial[0] = X;
    xglobal_initial[1] = Y;
    xglobal_initial[2] = Z;
    array_1d<double,3> xlocal_initial = prod(mRotationMatrix, (xglobal_initial - mCenterCoordinates) );
    return CallFunction(xlocal[0],xlocal[1],xlocal[2], t,xlocal_initial[0],xlocal_initial[1],xlocal_initial[2]);
}

/***********************************************************************************/
/***********************************************************************************/

double GenericFunctionUtility::CallFunction(
    const double x,
    const double y,
    const double z,
    const double t,
    const double X,
    const double Y,
    const double Z
    )
{
    mNameSpace["x"] = x;
    mNameSpace["y"] = y;
    mNameSpace["z"] = z;
    mNameSpace["X"] = X;
    mNameSpace["Y"] = Y;
    mNameSpace["Z"] = Z;
    mNameSpace["t"] = t;

    return mpExpression->value();
}

/***********************************************************************************/
/***********************************************************************************/

void GenericFunctionUtility::InitializeExprtk()
{
    // Defining table
    double& x = mNameSpace["x"];
    double& y = mNameSpace["y"];
    double& z = mNameSpace["z"];
    double& t = mNameSpace["t"];
    double& X = mNameSpace["X"];
    double& Y = mNameSpace["Y"];
    double& Z = mNameSpace["Z"];
    double& pi = mNameSpace["pi"];

    mpSymbolTable = new exprtk::symbol_table<double>();
    mpSymbolTable->add_variable("x",x);
    mpSymbolTable->add_variable("y",y);
    mpSymbolTable->add_variable("z",z);
    mpSymbolTable->add_variable("t",t);
    mpSymbolTable->add_variable("X",X);
    mpSymbolTable->add_variable("Y",Y);
    mpSymbolTable->add_variable("Z",Z);
    mpSymbolTable->add_variable("pi", pi);

    // Creating expression
    mpExpression = new exprtk::expression<double>();
    mpExpression->register_symbol_table(*mpSymbolTable);

    // Compiling expression
    mpParser = new exprtk::parser<double>();
    mpParser->compile(mFunctionBody, *mpExpression);
}

} /// namespace Kratos
