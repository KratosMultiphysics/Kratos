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
#include <regex>

// External includes
#include "tinyexpr/tinyexpr/tinyexpr.h"

// Project includes
#include "utilities/function_parser_utility.h"
#include "utilities/string_utilities.h"

namespace Kratos
{

// A check only to check that the x is not part of exp
inline bool CheckThereIsNotx(const std::string& rString)
{
    const auto first_x = rString.find(std::string("x"));
    if (first_x != std::string::npos) {
        const std::string exp_string = "exp";
        const char e_char = exp_string[0];
        const char x_char = exp_string[1];
        const char p_char = exp_string[2];
        for(std::size_t i = first_x; i < rString.size(); ++i) {
            if(rString[i] == x_char) {
                if (!(rString[i - 1] == e_char && rString[i + 1] == p_char)) {
                    return false;
                }
            }
        }
    }

    return true;
}

/***********************************************************************************/
/***********************************************************************************/

BasicGenericFunctionUtility::BasicGenericFunctionUtility(const std::string& rFunctionBody)
    : mFunctionBody(rFunctionBody)
{
    // Correcting function from python style to tinyexpr
    std::unordered_map<std::string,std::string> aux_replace_letters({{"**","^"}});
    for (auto& r_pair : aux_replace_letters) {
        mFunctionBody = StringUtilities::ReplaceAllSubstrings(mFunctionBody, r_pair.first, r_pair.second);
    }

    // Initialize parser classes
    InitializeParser();

    // Check if it depends on space
    const std::regex space_regex("\\b[xyz]\\b", std::regex_constants::icase);
    mDependsOnSpace = std::regex_search(mFunctionBody, space_regex);
}

/***********************************************************************************/
/***********************************************************************************/

BasicGenericFunctionUtility::BasicGenericFunctionUtility(BasicGenericFunctionUtility const& rOther)
    : mValues(rOther.mValues),
      mFunctionBody(rOther.mFunctionBody),
      mDependsOnSpace(rOther.mDependsOnSpace)
{
    // Initialize parser classes
    InitializeParser();
}

/***********************************************************************************/
/***********************************************************************************/

BasicGenericFunctionUtility::~BasicGenericFunctionUtility()
{
    for (std::size_t i = 0; i < mpTinyExpr.size(); ++i) {
        te_free(mpTinyExpr[i]);
        mpTinyExpr[i] =  nullptr;
    }
}

/***********************************************************************************/
/***********************************************************************************/

std::string BasicGenericFunctionUtility::FunctionBody()
{
    std::string aux_string = mFunctionBody;
    std::unordered_map<std::string,std::string> aux_replace_letters({{"**","^"}});
    for (auto& r_pair : aux_replace_letters) {
         aux_string = StringUtilities::ReplaceAllSubstrings(aux_string, r_pair.second, r_pair.first);
    }
    return aux_string;
}

/***********************************************************************************/
/***********************************************************************************/

double BasicGenericFunctionUtility::RotateAndCallFunction(
    const double x,
    const double y,
    const double z,
    const double t,
    const double X,
    const double Y,
    const double Z
    )
{
    return CallFunction(x, y, z, t, X, Y, Z);
}

/***********************************************************************************/
/***********************************************************************************/

double BasicGenericFunctionUtility::CallFunction(
    const double x,
    const double y,
    const double z,
    const double t,
    const double X,
    const double Y,
    const double Z
    )
{
    mValues[0] = x;
    mValues[1] = y;
    mValues[2] = z;
    mValues[3] = t;
    mValues[4] = X;
    mValues[5] = Y;
    mValues[6] = Z;

    // Default function
    if (mpTinyExpr.size() == 1) {
      return te_eval(mpTinyExpr[0]);
    } else { // Ternary expression
      if (te_eval(mpTinyExpr[0]) > 0.0) {
          return te_eval(mpTinyExpr[1]);
      } else {
          return te_eval(mpTinyExpr[2]);
      }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BasicGenericFunctionUtility::InitializeParser()
{
    // Initialize
    if (mpTinyExpr[0] == nullptr) {
        int err;

        /* Defining table */
        double& x = mValues[0];
        double& y = mValues[1];
        double& z = mValues[2];
        double& t = mValues[3];
        double& X = mValues[4];
        double& Y = mValues[5];
        double& Z = mValues[6];

        /* Store variable names and pointers. */
        const te_variable vars[] = {{"x", &x}, {"y", &y}, {"z", &z}, {"t", &t}, {"X", &X}, {"Y", &Y}, {"Z", &Z}};

        auto fct_checker = [](const bool IsValid, const std::string& rFunction, const int ErrPos){
            if (IsValid) return;

            std::stringstream ss;
            ss << "\nParsing error in function: " << rFunction << '\n';
            ss <<   "Error occurred near here : ";
            for(int i=0; i<ErrPos-1; ++i) ss << ' ';
            ss << "^ (char ["<< ErrPos-1 << "])\n";
            ss << "Check your locale (e.g. if \".\" or \",\" is used as decimal point)";

            KRATOS_ERROR << ss.str() << std::endl;
        };

        /* Compile the expression with variables. */
        const bool python_like_ternary = StringUtilities::ContainsPartialString(mFunctionBody, "if") ? true : false;
        if (!python_like_ternary) {
            mpTinyExpr[0] = te_compile(mFunctionBody.c_str(), vars, 7, &err);
            fct_checker(mpTinyExpr[0], mFunctionBody, err);
        } else { // Ternary operator
            mpTinyExpr.resize(3, nullptr);
            std::string condition, first_function, second_function;

            // C like ternary operator
            std::vector<std::string> splitted_string;
            KRATOS_ERROR_IF_NOT(StringUtilities::ContainsPartialString(mFunctionBody, "else")) << "Parsing error in function: " << mFunctionBody << " if defined, but not else" << std::endl;
            std::string aux_string = StringUtilities::ReplaceAllSubstrings(mFunctionBody, "if", "$");
            KRATOS_ERROR_IF(splitted_string.size() > 2) << "Nested ternary functions not supported" << std::endl;
            splitted_string = StringUtilities::SplitStringByDelimiter(aux_string, '$');
            first_function = splitted_string[0];
            aux_string = StringUtilities::ReplaceAllSubstrings(splitted_string[1], "else", "$");
            splitted_string = StringUtilities::SplitStringByDelimiter(aux_string, '$');
            condition = splitted_string[0];
            second_function = splitted_string[1];

            // Parsing the functions
            mpTinyExpr[1] = te_compile(first_function.c_str(), vars, 7, &err);
            fct_checker(mpTinyExpr[1], first_function, err);

            mpTinyExpr[2] = te_compile(second_function.c_str(), vars, 7, &err);
            fct_checker(mpTinyExpr[2], second_function, err);

            // Parsing the condition
            if (StringUtilities::ContainsPartialString(condition, "==")) {
                // Auxiliar string function
                const std::string aux_string = StringUtilities::ReplaceAllSubstrings(condition, "==", "=");
                splitted_string = StringUtilities::SplitStringByDelimiter(aux_string, '=');
                const std::string aux_function = "zIsEqual(" + splitted_string[0] + ", " + splitted_string[1] + ")";
                mpTinyExpr[0] = te_compile(aux_function.c_str(), vars, 7, &err);
            } else if (StringUtilities::ContainsPartialString(condition, "<=")) {
                // Auxiliar string function
                const std::string aux_string = StringUtilities::ReplaceAllSubstrings(condition, "<=", "<");
                splitted_string = StringUtilities::SplitStringByDelimiter(aux_string, '<');
                const std::string aux_function = "zIsLessEqual(" + splitted_string[0] + ", " + splitted_string[1] + ")";
                mpTinyExpr[0] = te_compile(aux_function.c_str(), vars, 7, &err);
            } else if (StringUtilities::ContainsPartialString(condition, ">=")) {
                // Auxiliar string function
                const std::string aux_string = StringUtilities::ReplaceAllSubstrings(condition, ">=", ">");
                splitted_string = StringUtilities::SplitStringByDelimiter(aux_string, '>');
                const std::string aux_function = "zIsGreaterEqual(" + splitted_string[0] + ", " + splitted_string[1] + ")";
                mpTinyExpr[0] = te_compile(aux_function.c_str(), vars, 7, &err);
            } else if (StringUtilities::ContainsPartialString(condition, "<")) {
                // Auxiliar string function
                splitted_string = StringUtilities::SplitStringByDelimiter(condition, '<');
                const std::string aux_function = "zIsLess(" + splitted_string[0] + ", " + splitted_string[1] + ")";
                mpTinyExpr[0] = te_compile(aux_function.c_str(), vars, 7, &err);
            } else if (StringUtilities::ContainsPartialString(condition, ">")) {
                // Auxiliar string function
                splitted_string = StringUtilities::SplitStringByDelimiter(condition, '>');
                const std::string aux_function = "zIsGreater(" + splitted_string[0] + ", " + splitted_string[1] + ")";
                mpTinyExpr[0] = te_compile(aux_function.c_str(), vars, 7, &err);
            } else {
                KRATOS_ERROR << "Cannot identify condition: " << condition << std::endl;
            }
            fct_checker(mpTinyExpr[0], condition, err);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

GenericFunctionUtility::GenericFunctionUtility(
    const std::string& rFunctionBody,
    Parameters LocalSystem
    ) : BasicGenericFunctionUtility(rFunctionBody)
{
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
}

/***********************************************************************************/
/***********************************************************************************/

GenericFunctionUtility::GenericFunctionUtility(GenericFunctionUtility const& rOther)
    : BasicGenericFunctionUtility(rOther),
      mUseLocalSystem(rOther.mUseLocalSystem),
      mRotationMatrix(rOther.mRotationMatrix),
      mCenterCoordinates(rOther.mCenterCoordinates)
{
    // Initialize parser classes
    InitializeParser();
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

} /// namespace Kratos
