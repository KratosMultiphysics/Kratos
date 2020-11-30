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

// Project includes
#include "includes/global_variables.h"
#include "utilities/function_parser_utility.h"

namespace Kratos
{

GenericFunctionUtility::GenericFunctionUtility(
    const std::string& rFunctionBody,
    Parameters LocalSystem
    ) : mFunctionBody(ReplaceAllSubstrings(rFunctionBody, "**", "^")) // Correcting function from python style to exprtk
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

    // Defining table
    double& x = mNameSpace["x"];
    double& y = mNameSpace["y"];
    double& z = mNameSpace["z"];
    double& t = mNameSpace["t"];
    double& X = mNameSpace["X"];
    double& Y = mNameSpace["Y"];
    double& Z = mNameSpace["Z"];
    double& pi = mNameSpace["pi"];

    mSymbolTable.add_variable("x",x);
    mSymbolTable.add_variable("y",y);
    mSymbolTable.add_variable("z",z);
    mSymbolTable.add_variable("t",t);
    mSymbolTable.add_variable("X",X);
    mSymbolTable.add_variable("Y",Y);
    mSymbolTable.add_variable("Z",Z);
    mSymbolTable.add_variable("pi", pi);

    // Creating expression
    mExpression.register_symbol_table(mSymbolTable);

    // Compiling expression
    mParser.compile(mFunctionBody, mExpression);

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
    // Defining table
    double& x = mNameSpace["x"];
    double& y = mNameSpace["y"];
    double& z = mNameSpace["z"];
    double& t = mNameSpace["t"];
    double& X = mNameSpace["X"];
    double& Y = mNameSpace["Y"];
    double& Z = mNameSpace["Z"];
    double& pi = mNameSpace["pi"];

    mSymbolTable.add_variable("x",x);
    mSymbolTable.add_variable("y",y);
    mSymbolTable.add_variable("z",z);
    mSymbolTable.add_variable("t",t);
    mSymbolTable.add_variable("X",X);
    mSymbolTable.add_variable("Y",Y);
    mSymbolTable.add_variable("Z",Z);
    mSymbolTable.add_variable("pi", pi);

    // Creating expression
    mExpression.register_symbol_table(mSymbolTable);

    // Compiling expression
    mParser.compile(mFunctionBody, mExpression);
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

    return mExpression.value();
}

/***********************************************************************************/
/***********************************************************************************/

std::string GenericFunctionUtility::ReplaceAllSubstrings(std::string str, const std::string& from, const std::string& to)
{
    std::size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}

} /// namespace Kratos
