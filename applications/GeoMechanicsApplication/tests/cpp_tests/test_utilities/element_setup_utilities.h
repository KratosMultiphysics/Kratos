// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

#include "geo_aliases.h"
#include "geometries/point.h"
#include "includes/condition.h"
#include "includes/element.h"

#include <vector>

namespace Kratos::Testing
{

class ElementSetupUtilities
{
public:
    static std::vector<Kratos::Point> CreatePointsFor2D3NElement();
    static std::vector<Kratos::Point> CreatePointsFor2D3NLineEntity();
    static std::vector<Kratos::Point> CreatePointsFor2D6NElement();
    static std::vector<Kratos::Point> CreatePointsFor2D10NElement();
    static std::vector<Kratos::Point> CreatePointsFor2D15NElement();
    static std::vector<Kratos::Point> CreatePointsFor3D10NElement();

    static Element::Pointer Create2D3NElement(const PointerVector<Node>& rNodes,
                                              const Properties::Pointer& rProperties);
    static Element::Pointer Create2D3NElement();

    static Element::Pointer Create2D6NElement(const PointerVector<Node>& rNodes,
                                              const Properties::Pointer& rProperties);
    static Element::Pointer Create2D6NElement();

    static Element::Pointer Create2D6NDiffOrderElement(const PointerVector<Node>& rNodes,
                                                       const Properties::Pointer& rProperties);
    static Element::Pointer Create2D6NDiffOrderElement();

    static Element::Pointer Create2D10NElement(const PointerVector<Node>& rNodes,
                                               const Properties::Pointer& rProperties);
    static Element::Pointer Create2D10NElement();

    static Element::Pointer Create2D15NElement(const PointerVector<Node>& rNodes,
                                               const Properties::Pointer& rProperties);
    static Element::Pointer Create2D15NElement();

    static Element::Pointer Create3D10NElement(const PointerVector<Node>& rNodes,
                                               const Properties::Pointer& rProperties);
    static Element::Pointer Create3D10NElement();

    static Condition::Pointer Create2D3NLineCondition(const PointerVector<Node>& rNodes,
                                                      const Properties::Pointer& rProperties);
    static Condition::Pointer Create2D3NLineCondition();

    template <class EntityPointerType>
    static void AddVariablesToEntity(EntityPointerType& rpEntity,
                                     const Kratos::Geo::ConstVariableDataRefs& rSolutionStepVariables,
                                     const Kratos::Geo::ConstVariableRefs& rDegreesOfFreedom = {})
    {
        auto p_variable_list = make_intrusive<VariablesList>();
        for (const auto& r_variable_ref : rSolutionStepVariables) {
            p_variable_list->Add(r_variable_ref);
        }

        for (auto& r_node : rpEntity->GetGeometry()) {
            r_node.SetSolutionStepVariablesList(p_variable_list);
            for (const auto& r_degree_of_freedom : rDegreesOfFreedom) {
                r_node.AddDof(r_degree_of_freedom.get());
            }
        }
    }
};

} // namespace Kratos::Testing
