// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef UNIVERSAL_FILE_IO_H
#define UNIVERSAL_FILE_IO_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/model_part.h"

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

Follows the universal I-deas file convetion so that results can be imported in SIEMENS NX:
https://docs.plm.automation.siemens.com/tdoc/nx/10/nx_help/#uid:index_advanced:xid602249:id625716:id625821

*/

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) UniversalFileIO
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UniversalFileIO( ModelPart& OutputModelPart, std::string OutputFilenameWithoutExtension, std::string WriteConditionsFlag, Parameters NodalResults );

    /// Destructor.
    virtual ~UniversalFileIO()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void InitializeLogging();

    // --------------------------------------------------------------------------
    void WriteMeshToResultFile();

    // --------------------------------------------------------------------------
    void InitializeOutputFile();

    // --------------------------------------------------------------------------
    void WriteUnits();

    // --------------------------------------------------------------------------
    void WriteNodes();

    // --------------------------------------------------------------------------
    void WriteElements();

    // --------------------------------------------------------------------------
    void WriteAllElementsButNoConditions();

    // --------------------------------------------------------------------------
    void WriteConditionsAsDummyElements();

    // --------------------------------------------------------------------------
    void LogNodalResults( const int optimizationIteration );

    // ==============================================================================

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a std::string.
    virtual std::string Info() const
    {
        return "UniversalFileIO";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "UniversalFileIO";
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

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    ModelPart& mrOutputModelPart;
    std::string mOutputFilenameWithExtension;
    Parameters mrNodalResults;
    std::string mWriteConditionsFlag;

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
//      UniversalFileIO& operator=(UniversalFileIO const& rOther);

    /// Copy constructor.
//      UniversalFileIO(UniversalFileIO const& rOther);


    ///@}

}; // Class UniversalFileIO

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // UNIVERSAL_FILE_IO_H
