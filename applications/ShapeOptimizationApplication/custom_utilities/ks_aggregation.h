// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//
// ==============================================================================

#ifndef KS_AGGREGATION_UTILITY_H
#define KS_AGGREGATION_UTILITY_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/binbased_fast_point_locator.h"
#include "input_output/vtk_output.h"
#include "processes/integration_values_extrapolation_to_nodes_process.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{


///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

*/

class KreisselmeierSteinhauserAggregationUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KreisselmeierSteinhauserAggregationUtility
    KRATOS_CLASS_POINTER_DEFINITION(KreisselmeierSteinhauserAggregationUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KreisselmeierSteinhauserAggregationUtility( ModelPart& modelPart )
        : mrModelPart( modelPart )
    {
        mpPointLocator = Kratos::make_unique<BinBasedFastPointLocator<3>>(mrModelPart);
        mpPointLocator->UpdateSearchDatabase();
    }

    /// Destructor.
    virtual ~KreisselmeierSteinhauserAggregationUtility()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    double CalculateValue(const Variable<double>& rVariable, const double rho) const
    {
        return InternalCalculateAggregationValue(mrModelPart, rVariable, rho);
    }

    double CalculateValue(const ModelPart& rModelPart, const Variable<double>& rVariable, const double rho) const
    {
        return InternalCalculateAggregationValue(rModelPart, rVariable, rho);
    }

    double CalculateValue(const ModelPart& rModelPart, const Variable<Matrix>& rVariable, const double rho) const
    {
        return InternalCalculateAggregationValue(rModelPart, rVariable, rho);
    }

    // --------------------------------------------------------------------------

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "KreisselmeierSteinhauserAggregationUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "KreisselmeierSteinhauserAggregationUtility";
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
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


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
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    // Initialized by class constructor
    ModelPart& mrModelPart;
    BinBasedFastPointLocator<3>::UniquePointer mpPointLocator;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{
    double InternalCalculateAggregationValue(const ModelPart& rModelPart, const Variable<double>& rVariable, const double rho) const
    {
        double aggregate_value = 0.0;
        double max_value = 0.0;
        std::vector<double> stresses_vector;
        const int max_results = 10000;
        Vector weights;
        typename BinBasedFastPointLocator<3>::ResultContainerType results(max_results);
        typename BinBasedFastPointLocator<3>::ResultIteratorType result_begin =
            results.begin();

        stresses_vector.reserve((rModelPart.NumberOfNodes()));
        for (auto& r_node : rModelPart.Nodes())
        {
            Element::Pointer r_host_element;
            bool is_found = false;
            is_found = mpPointLocator->FindPointOnMesh(r_node.Coordinates(), weights,
                                                r_host_element, result_begin, max_results, 1e-1);
            if(is_found){
                const int num_gps = r_host_element->GetGeometry().IntegrationPoints(r_host_element->GetIntegrationMethod()).size();
                std::vector<double> stress_vector(num_gps);
                double value = 0.0;
                const auto& r_geom = r_host_element->GetGeometry();
                for(unsigned int i = 0; i< r_geom.size(); ++i )
                {
                    value += r_geom[i].GetValue(rVariable) * weights[i];
                }
                if(std::abs(value) > max_value)
                    max_value = value;
                stresses_vector.push_back(value);
            }
            else {
                KRATOS_WARNING("KreisselmeierSteinhauserAggregationUtility")<<"Node "<< r_node <<" not found !"<<std::endl;
            }
        }
        for (const auto stress : stresses_vector)
            aggregate_value += std::exp(rho * ((stress)/max_value));

        aggregate_value = (1/rho) * std::log(aggregate_value);
        return aggregate_value;
    }

    double InternalCalculateAggregationValue(const ModelPart& rModelPart, const Variable<Matrix>& rVariable, const double rho) const
    {
        double aggregate_value = 0.0;
        double max_value = 0.0;
        std::vector<double> stresses_vector;
        const int max_results = 10000;
        Vector weights;
        typename BinBasedFastPointLocator<3>::ResultContainerType results(max_results);
        typename BinBasedFastPointLocator<3>::ResultIteratorType result_begin =
            results.begin();

        stresses_vector.reserve((rModelPart.NumberOfNodes()));
        for (auto& r_node : rModelPart.Nodes())
        {
            Element::Pointer r_host_element;
            bool is_found = false;
            is_found = mpPointLocator->FindPointOnMesh(r_node.Coordinates(), weights,
                                                r_host_element, result_begin, max_results, 1e-1);
            if(is_found){
                const int num_gps = r_host_element->GetGeometry().IntegrationPoints(r_host_element->GetIntegrationMethod()).size();
                std::vector<double> stress_vector(num_gps);
                double value = 0.0;
                const auto& r_geom = r_host_element->GetGeometry();
                for(unsigned int i = 0; i< r_geom.size(); ++i )
                {
                    value += r_geom[i].GetValue(rVariable)(0,0) * weights[i];
                }
                if(std::abs(value) > max_value)
                    max_value = value;
                stresses_vector.push_back(value);
            }
            else {
                KRATOS_WARNING("KreisselmeierSteinhauserAggregationUtility")<<"Node "<< r_node <<" not found !"<<std::endl;
            }
        }
        for (const auto stress : stresses_vector)
            aggregate_value += std::exp(rho * ((stress)/max_value));

        aggregate_value = (1/rho) * std::log(aggregate_value);
        return aggregate_value;
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class KreisselmeierSteinhauserAggregationUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif