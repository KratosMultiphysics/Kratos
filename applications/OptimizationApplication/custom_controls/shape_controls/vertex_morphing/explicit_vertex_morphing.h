// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

#ifndef EXPLICIT_VERTEX_MORPHING_H
#define EXPLICIT_VERTEX_MORPHING_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "mapper_vertex_morphing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "custom_controls/shape_controls/shape_control.h"

// ==============================================================================

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

class KRATOS_API(OPTIMIZATION_APPLICATION) ExplicitVertexMorphing : public ShapeControl
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of ExplicitVertexMorphing
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitVertexMorphing);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ExplicitVertexMorphing( std::string ControlName, Model& rModel, Parameters ControlSettings )
        : ShapeControl(ControlName,rModel,ControlSettings){}

    /// Destructor.
    virtual ~ExplicitVertexMorphing()
    {
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void Initialize() override {
        // KRATOS_ERROR <<" HI REZA YOU CALLED ME"<<std::endl;

    };
    // --------------------------------------------------------------------------
    void Update() override {};  
    // --------------------------------------------------------------------------
    void MapControlUpdate(const Variable<array_3d> &rOriginVariable, const Variable<array_3d> &rDestinationVariable) override{};
    // --------------------------------------------------------------------------
    void MapFirstDerivative(const Variable<array_3d> &rDerivativeVariable, const Variable<array_3d> &rMappedDerivativeVariable) override{};  

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
    virtual std::string Info() const override
    {
        return "ExplicitVertexMorphing";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ExplicitVertexMorphing";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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

    // Initialized by class constructor

    
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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // --------------------------------------------------------------------------

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
//      ExplicitVertexMorphing& operator=(ExplicitVertexMorphing const& rOther);

    /// Copy constructor.
//      ExplicitVertexMorphing(ExplicitVertexMorphing const& rOther);


    ///@}

}; // Class ExplicitVertexMorphing

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // EXPLICIT_VERTEX_MORPHING_H
