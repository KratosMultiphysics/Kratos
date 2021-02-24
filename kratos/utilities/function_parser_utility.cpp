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
#include "tinyexpr/tinyexpr/tinyexpr.h"

// Project includes
#include "includes/global_variables.h"
#include "utilities/function_parser_utility.h"
#include "utilities/string_utilities.h"

namespace Kratos
{

GenericFunctionUtility::GenericFunctionUtility(
    const std::string& rFunctionBody,
    Parameters LocalSystem
    ) : mOriginalFunctionBody(rFunctionBody),
        mFunctionBody(rFunctionBody)
{
    // Correcting function from python style to tinyexpr
    std::unordered_map<std::string,std::string> aux_replace_letters({{"**","^"},{"X","a"},{"Y","b"},{"Z","c"}});
    for (auto& r_pair : aux_replace_letters) {
        mFunctionBody = StringUtilities::ReplaceAllSubstrings(mFunctionBody, r_pair.first, r_pair.second);
    }

    // Defining namespace
    mNameSpace.insert(std::pair<std::string, double>("x", 0.0));
    mNameSpace.insert(std::pair<std::string, double>("y", 0.0));
    mNameSpace.insert(std::pair<std::string, double>("z", 0.0));
    mNameSpace.insert(std::pair<std::string, double>("t", 0.0));
    mNameSpace.insert(std::pair<std::string, double>("X", 0.0));
    mNameSpace.insert(std::pair<std::string, double>("Y", 0.0));
    mNameSpace.insert(std::pair<std::string, double>("Z", 0.0));
    mNameSpace.insert(std::pair<std::string, double>("pi", Globals::Pi));

    // Initialize parser classes
    InitializeParser();

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
    if (mOriginalFunctionBody.find(std::string("x")) == std::string::npos &&
        mOriginalFunctionBody.find(std::string("y")) == std::string::npos &&
        mOriginalFunctionBody.find(std::string("z")) == std::string::npos &&
        mOriginalFunctionBody.find(std::string("X")) == std::string::npos &&
        mOriginalFunctionBody.find(std::string("Y")) == std::string::npos &&
        mOriginalFunctionBody.find(std::string("Z")) == std::string::npos) {
        mDependsOnSpace = false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

GenericFunctionUtility::GenericFunctionUtility(GenericFunctionUtility const& rOther)
    : mNameSpace(rOther.mNameSpace),
      mOriginalFunctionBody(rOther.mOriginalFunctionBody),
      mFunctionBody(rOther.mFunctionBody),
      mDependsOnSpace(rOther.mDependsOnSpace),
      mUseLocalSystem(rOther.mUseLocalSystem),
      mRotationMatrix(rOther.mRotationMatrix),
      mCenterCoordinates(rOther.mCenterCoordinates)
{
    // Initialize parser classes
    InitializeParser();
}

/***********************************************************************************/
/***********************************************************************************/

GenericFunctionUtility::~GenericFunctionUtility()
{
    te_free(mpTinyExpr);
}

/***********************************************************************************/
/***********************************************************************************/

std::string GenericFunctionUtility::FunctionBody()
{
    return mOriginalFunctionBody;
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

    return te_eval(mpTinyExpr);
}

/***********************************************************************************/
/***********************************************************************************/

void GenericFunctionUtility::InitializeParser()
{
    // Initialize
    if (mpTinyExpr == nullptr) {
        int err;

        /* Defining table */
        double& x = mNameSpace["x"];
        double& y = mNameSpace["y"];
        double& z = mNameSpace["z"];
        double& t = mNameSpace["t"];
        double& X = mNameSpace["X"];
        double& Y = mNameSpace["Y"];
        double& Z = mNameSpace["Z"];
        double& pi = mNameSpace["pi"];

        /* Store variable names and pointers. */
        const te_variable vars[] = {{"x", &x}, {"y", &y}, {"z", &z}, {"t", &t}, {"a", &X}, {"b", &Y}, {"c", &Z}, {"pi", &pi}};

        /* Compile the expression with variables. */
        mpTinyExpr = te_compile(mFunctionBody.c_str(), vars, 8, &err);
    }
}

} /// namespace Kratos
