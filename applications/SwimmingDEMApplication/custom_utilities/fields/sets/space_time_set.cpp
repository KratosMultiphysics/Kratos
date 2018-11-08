#include "space_time_set.h"

namespace Kratos
{
SpaceTimeSet::Pointer SetUnion(SpaceTimeSet::Pointer set_1, SpaceTimeSet::Pointer set_2)
{
    SpaceTimeSet::Pointer set(new SpaceTimeSet());
    std::vector<std::vector<SpaceTimeRule::Pointer> > rules_1 = set_1->GetRules();
    std::vector<std::vector<SpaceTimeRule::Pointer> > rules_2 = set_2->GetRules();

    for (unsigned int i = 0; i < rules_1.size(); ++i){
        set->AddOrRules(rules_1[i]);
    }

    for (unsigned int i = 0; i < rules_2.size(); ++i){
        set->AddOrRules(rules_2[i]);
    }

    return set;
}

SpaceTimeSet::Pointer SetIntersection(SpaceTimeSet::Pointer set_1, SpaceTimeSet::Pointer set_2)
{
    SpaceTimeSet::Pointer set(new SpaceTimeSet());
    std::vector<std::vector<SpaceTimeRule::Pointer> > rules_1 = set_1->GetRules();
    std::vector<std::vector<SpaceTimeRule::Pointer> > rules_2 = set_2->GetRules();
    std::vector<SpaceTimeRule::Pointer> u_rules_1;

    for (unsigned int i = 0; i < rules_1.size(); ++i){

        for (unsigned int j = 0; j < rules_1[i].size(); ++j){
          u_rules_1.push_back(rules_1[i][j]);
        }
    }

    for (unsigned int i = 0; i < rules_2.size(); ++i){
        set->AddOrRules(rules_2[i]);
    }

    for (unsigned int i = 0; i < rules_2.size(); ++i){
        set->AddAndRules(u_rules_1);
    }

    return set;
}
}
