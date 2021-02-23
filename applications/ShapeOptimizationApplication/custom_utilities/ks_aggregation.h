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

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{
    double InternalCalculateAggregationValue(const ModelPart& rModelPart, const Variable<double>& rVariable, const double rho) const
    {
        double value = 0.0;
        const auto& r_process_info = rModelPart.GetProcessInfo();
        // GetGeometry().IntegrationPoints(mIntegrationMethod).size()
        for (auto& r_elem : rModelPart.Elements())
        {
            const int num_gps = r_elem.GetGeometry().IntegrationPoints(r_elem.GetIntegrationMethod()).size();
            std::vector<double> tmp_stress_vector(num_gps);
            r_elem.CalculateOnIntegrationPoints(rVariable, tmp_stress_vector, r_process_info);          
        }
        return value;
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