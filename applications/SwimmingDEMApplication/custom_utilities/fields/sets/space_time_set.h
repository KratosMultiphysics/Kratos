#ifndef KRATOS_SPACE_TIME_SET_H
#define KRATOS_SPACE_TIME_SET_H

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
#include "space_time_rule.h"


namespace Kratos
{
class SpaceTimeSet
{
public:

KRATOS_CLASS_POINTER_DEFINITION(SpaceTimeSet);

/// Default constructor.

SpaceTimeSet(const double min_time,
             const double max_time,
             const double min_x,
             const double min_y,
             const double min_z,
             const double max_x,
             const double max_y,
             const double max_z)
{
    BoundingBoxRule::Pointer b_box_rule(new BoundingBoxRule(min_time,
                                                            max_time,
                                                            min_x,
                                                            min_y,
                                                            min_z,
                                                            max_x,
                                                            max_y,
                                                            max_z));
    std::vector<SpaceTimeRule::Pointer> union_rules;
    union_rules.push_back(b_box_rule);
    mUnionOfRules.push_back(union_rules);
}

SpaceTimeSet()
{
    BoundingBoxRule::Pointer b_box_rule(new BoundingBoxRule());
    std::vector<SpaceTimeRule::Pointer> union_rules;
    union_rules.push_back(b_box_rule);
    mUnionOfRules.push_back(union_rules);
}

/// Destructor.

virtual ~SpaceTimeSet(){}

//***************************************************************************************************************
//***************************************************************************************************************

void AddAndRule(SpaceTimeRule::Pointer p_rule)
{
    for (unsigned int i = 0; i < mUnionOfRules.size(); ++i){
        mUnionOfRules[i].push_back(p_rule);
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void AddOrRule(SpaceTimeRule::Pointer p_rule)
{
std::vector<SpaceTimeRule::Pointer> p_rule_array;
p_rule_array.push_back(p_rule);
mUnionOfRules.push_back(p_rule_array);

}

//***************************************************************************************************************
//***************************************************************************************************************

void AddAndRules(std::vector<SpaceTimeRule::Pointer> p_rules)
{
    for (unsigned int i = 0; i < mUnionOfRules.size(); ++i){

        for (unsigned int j = 0; j < p_rules.size(); ++j){
            mUnionOfRules[i].push_back(p_rules[j]);
        }
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

void AddOrRules(std::vector<SpaceTimeRule::Pointer> p_rules)
{
    mUnionOfRules.push_back(p_rules);
}

//***************************************************************************************************************
//***************************************************************************************************************

bool IsIn(const double time, const double coor_x, const double coor_y, const double coor_z)
{
    bool rules_are_met = true;

    for (unsigned int i = 0; i < mUnionOfRules.size(); ++i){

        for (unsigned int j = 0; j < mUnionOfRules[i].size(); ++j){
            rules_are_met = rules_are_met && mUnionOfRules[i][j]->CheckIfRuleIsMet(time, coor_x, coor_y, coor_z);
        }

        if (rules_are_met){
            return true;
        }

    }

    return false;
}

//***************************************************************************************************************
//***************************************************************************************************************

///@}
///@name Inquiry
///@{
double GetLowTime(){return mLowTime;}
double GetHighTime(){return mHighTime;}
double GetLowX(){return mLowX;}
double GetHighX(){return mHighX;}
double GetLowY(){return mLowY;}
double GetHighY(){return mHighY;}
double GetLowZ(){return mLowZ;}
double GetHighZ(){return mHighZ;}
std::vector<std::vector<SpaceTimeRule::Pointer> > GetRules(){return mUnionOfRules;}


///@}
///@name Input and output
///@{

/// Turn back information as a stemplate<class T, std::size_t dim> tring.

virtual std::string Info() const
{
    return "";
}

/// Print information about this object.

virtual void PrintInfo(std::ostream& rOStream) const
{
}

/// Print object's data.

virtual void PrintData(std::ostream& rOStream) const
{
}


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
double mLowTime;
double mHighTime;
double mLowX;
double mHighX;
double mLowY;
double mHighY;
double mLowZ;
double mHighZ;
std::vector<std::vector<SpaceTimeRule::Pointer> > mUnionOfRules;

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
SpaceTimeSet & operator=(SpaceTimeSet const& rOther);


///@}

}; // Class SpaceTimeSet

SpaceTimeSet::Pointer SetUnion(SpaceTimeSet::Pointer set_1, SpaceTimeSet::Pointer set_2);
SpaceTimeSet::Pointer SetIntersection(SpaceTimeSet::Pointer set_1, SpaceTimeSet::Pointer set_2);

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
#endif // KRATOS_SPACE_TIME_SET_H
