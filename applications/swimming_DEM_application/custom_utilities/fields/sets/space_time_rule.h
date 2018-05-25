#ifndef KRATOS_SPACE_TIME_RULE_H
#define KRATOS_SPACE_TIME_RULE_H
// /* External includes */

// System includes

// Project includes
#include "includes/variables.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "../real_functions.h"
#include "../real_field.h"


namespace Kratos
{
class SpaceTimeRule
{
public:

KRATOS_CLASS_POINTER_DEFINITION(SpaceTimeRule);

/// Default constructor.

SpaceTimeRule(){}

/// Destructor.

virtual ~SpaceTimeRule(){}

//***************************************************************************************************************
//***************************************************************************************************************

virtual bool CheckIfRuleIsMet(const double time, const double coor_x, const double coor_y, const double coor_z){return true;}

//***************************************************************************************************************
//***************************************************************************************************************

///@}
///@name Inquiry
///@{

///@}
///@name Input and output
///@{

/// Turn back information as a stemplate<class T, std::size_t dim> tring.

virtual std::string Info() const
{
    return "SpaceTimeRule";
}

/// Print information about this object.

virtual void PrintInfo(std::ostream& rOStream) const{}

/// Print object's data.

virtual void PrintData(std::ostream& rOStream) const{}


///@}
///@name Friends
///@{

///@}

protected:
///@name Protected static Member r_variables
///@{


///@}
///@name Protected member r_variables
///@{ template<class T, std::size_t dim>


///@}
///@name Protected Operators
///@{


///@}
///@name Protected Operations
///@{


///@}
///@name Protected  Access
///@{

///@}
///@name Protected Inquiry
///@{


///@}
///@name Protected LifeCycle
///@{


///@}

private:

///@name Static Member r_variables
///@{


///@}
///@name Member r_variables
///@{

///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{


///@}
///@name Private  Access
///@{


///@}
///@name Private Inquiry
///@{


///@}
///@name Un accessible methods
///@{

/// Assignment operator.
SpaceTimeRule & operator=(SpaceTimeRule const& rOther);


///@}

}; // Class SpaceTimeRule

class BoundingBoxRule: public SpaceTimeRule
{
public:

BoundingBoxRule(): SpaceTimeRule()
{
    mLowTime  = - std::numeric_limits<double>::max();
    mHighTime =   std::numeric_limits<double>::max();
    mLowX     = - std::numeric_limits<double>::max();
    mHighX    =   std::numeric_limits<double>::max();
    mLowY     = - std::numeric_limits<double>::max();
    mHighY    =   std::numeric_limits<double>::max();
    mLowZ     = - std::numeric_limits<double>::max();
    mHighZ    =   std::numeric_limits<double>::max();
}

BoundingBoxRule(const double min_time,
                const double max_time,
                const double min_x,
                const double max_x,
                const double min_y,
                const double max_y,
                const double min_z,
                const double max_z): SpaceTimeRule()
{
    mLowTime  = min_time;
    mHighTime = max_time;
    mLowX     = min_x;
    mHighX    = max_x;
    mLowY     = min_y;
    mHighY    = max_y;
    mLowZ     = min_z;
    mHighZ    = max_z;
}

~BoundingBoxRule(){}

void SetTimeBoundingInterval(const double& low, const double& high)
{
    mLowTime  = low;
    mHighTime = high;
    Check();
}

void SetXBoundingInterval(const double& low, const double& high)
{
    mLowX  = low;
    mHighX = high;
    Check();
}

void SetYBoundingInterval(const double& low, const double& high)
{
    mLowY  = low;
    mHighY = high;
    Check();
}

void SetZBoundingInterval(const double& low, const double& high)
{
mLowZ  = low;
mHighZ = high;
Check();
}

void SetSpaceTimeBoundingBox(const array_1d<double, 4>& low, const array_1d<double, 4>& high)
{
    mLowTime  =  low[0];
    mHighTime = high[0];
    mLowX     =  low[1];
    mHighX    = high[1];
    mLowY     =  low[2];
    mHighY    = high[2];
    mLowZ     =  low[3];
    mHighZ    = high[3];
    Check();
}

bool CheckIfRuleIsMet(const double time, const double coor_x, const double coor_y, const double coor_z) override
{
    bool low  = time < mLowTime  || coor_x < mLowX  || coor_y < mLowY  || coor_z < mLowZ;
    bool high = time > mHighTime || coor_x > mHighX || coor_y > mHighY || coor_z > mHighZ;

    return !(high || low);
}

virtual std::string Info() const override
{
    std::ostringstream os;
    os << "Bounding box limits : "  << std::endl
       << "min time: " << mLowTime  << std::endl
       << "max time: " << mHighTime << std::endl
       << "min x : "   << mLowX     << std::endl
       << "max x : "   << mHighX    << std::endl
       << "min y : "   << mLowY     << std::endl
       << "max y : "   << mHighY    << std::endl
       << "min z : "   << mLowZ     << std::endl
       << "max z : "   << mHighZ    << std::endl;

    std::string info = os.str();

    return info;
}

void PrintData( std::ostream& rOStream = std::cout ) const override
{
    rOStream << "Bounding box limits : "  << std::endl
             << "min time: " << mLowTime  << std::endl
             << "max time: " << mHighTime << std::endl
             << "min x : "   << mLowX     << std::endl
             << "max x : "   << mHighX    << std::endl
             << "min y : "   << mLowY     << std::endl
             << "max y : "   << mHighY    << std::endl
             << "min z : "   << mLowZ     << std::endl
             << "max z : "   << mHighZ    << std::endl;
}

protected:
void Check()
{
    KRATOS_TRY

    if (mHighTime < mLowTime || mHighX < mLowX || mHighY < mLowY || mHighZ < mLowZ){
        KRATOS_THROW_ERROR(std::runtime_error, "Entered low bounds greater than corresponding bounds", "");
    }

    KRATOS_CATCH("")
}

private:

double mLowTime;
double mHighTime;
double mLowX;
double mHighX;
double mLowY;
double mHighY;
double mLowZ;
double mHighZ;
};

//***************************************************************************************************************
//***************************************************************************************************************

class MoreThanRule: public SpaceTimeRule
{
public:

MoreThanRule(RealField::Pointer field, const double value)
    : mValue(value),
      mpField(field),
      mpLowField(field),
      mCase(1){}

MoreThanRule(const double value, RealField::Pointer field)
    : mValue(value),
      mpField(field),
      mpLowField(field),
      mCase(2){}

MoreThanRule(RealField::Pointer field_high, RealField::Pointer field_low)
    : mValue(0.0),
      mpField(field_high),
      mpLowField(field_low),
      mCase(3){}

~MoreThanRule(){}

bool CheckIfRuleIsMet(const double time, const double coor_x, const double coor_y, const double coor_z) override
{
    array_1d<double, 3> coor;
    coor[0] = coor_x;
    coor[1] = coor_y;
    coor[2] = coor_z;

    if (mCase == 1){
        return mpField->Evaluate(time, coor) > mValue;
    }

    else if (mCase == 2){
        return mpField->Evaluate(time, coor) < mValue;
    }

    else {
        return mpField->Evaluate(time, coor) > mpLowField->Evaluate(time, coor);
    }
}

protected:

void Check(){}

private:

double mValue;
RealField::Pointer mpField;
RealField::Pointer mpLowField;
int mCase;
};

//***************************************************************************************************************
//***************************************************************************************************************

class EqualToRule: public SpaceTimeRule
{
public:
EqualToRule(RealField::Pointer field, const double value, const double tol)
    : mValue(value),
      mTol(tol),
      mpField1(field),
      mpField2(field),
      mCase(1)
{
    this->Check();
}

EqualToRule(RealField::Pointer field_1, RealField::Pointer field_2, const double tol)
    : mValue(0.0),
      mTol(tol),
      mpField1(field_1),
      mpField2(field_2),
      mCase(2)
{
    this->Check();
}

~EqualToRule(){}

bool CheckIfRuleIsMet(const double time, const double coor_x, const double coor_y, const double coor_z) override
{
    array_1d<double, 3> coor;
    coor[0] = coor_x;
    coor[1] = coor_y;
    coor[2] = coor_z;

    double valuation_1 = mpField1->Evaluate(time, coor);

    if (mCase == 1){
        return valuation_1 < mValue + mTol && valuation_1 > mValue - mTol;
    }

    else {
        double valuation_2 = mpField2->Evaluate(time, coor);
        return valuation_1 < valuation_2 + mTol && valuation_1 > valuation_2 - mTol;
    }
}

protected:

void Check()
{
    KRATOS_TRY

    if (mTol < 0.0){
        KRATOS_THROW_ERROR(std::runtime_error, "Entered tolerance must be a positive number", "");
    }

    KRATOS_CATCH("")
}

private:

double mValue;
double mTol;
RealField::Pointer mpField1;
RealField::Pointer mpField2;
int mCase;
};

//***************************************************************************************************************
//***************************************************************************************************************

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
#endif // KRATOS_SPACE_TIME_RULE_H
