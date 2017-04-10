// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef UNIVERSAL_FILE_IO_H
#define UNIVERSAL_FILE_IO_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iomanip>      // for setprecision

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "../../kratos/includes/define.h"
#include "../../kratos/processes/process.h"
#include "../../kratos/includes/node.h"
#include "../../kratos/includes/element.h"
#include "../../kratos/includes/model_part.h"
#include "../../kratos/includes/kratos_flags.h"
#include "shape_optimization_application.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

using namespace std;

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

class UniversalFileIO
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UniversalFileIO( ModelPart& designSurface, Parameters& optimizationSettings )
        : mrDesignSurface( designSurface )
    {
        setPrecisionForOutput();
        mOutputFilename = createCompleteOutputFilenameWithPath( optimizationSettings );
    }

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
    void setPrecisionForOutput()
    {
        cout.precision(12);
    }

    // --------------------------------------------------------------------------
    string createCompleteOutputFilenameWithPath( Parameters& optimizationSettings  )
    {
        string outputDirectory = optimizationSettings["output"]["output_directory"].GetString();
        string outputFilename = outputDirectory + "/" + optimizationSettings["output"]["design_history_filename"].GetString() + ".unv";
        return outputFilename;
    }

    // --------------------------------------------------------------------------
    void initializeLogging()
    {
        KRATOS_TRY;

        ofstream outputFile;
        outputFile.open(mOutputFilename, ios::out | ios::trunc );

        outputFile << "Writing something initial.\n";

        // const unsigned int domain_size = mrDesignSurface.GetProcessInfo().GetValue(DOMAIN_SIZE);

        // for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        // {

        // }

        outputFile.close();

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void logCurrentDesign( const int optimizationIteration )
    {
        KRATOS_TRY;

        ofstream outputFile;
        outputFile.open(mOutputFilename, ios::out | ios::app );

        outputFile << "Writing following optimization step: "<< optimizationIteration << ".\n";

        // const unsigned int domain_size = mrDesignSurface.GetProcessInfo().GetValue(DOMAIN_SIZE);

        // for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
        // {

        // }

        outputFile.close();

        KRATOS_CATCH("");
    }         

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

    /// Turn back information as a string.
    virtual string Info() const
    {
        return "UniversalFileIO";
    }

    /// Print information about this object.
    virtual void PrintInfo(ostream& rOStream) const
    {
        rOStream << "UniversalFileIO";
    }

    /// Print object's data.
    virtual void PrintData(ostream& rOStream) const
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
    ModelPart& mrDesignSurface;
    string mOutputFilename;

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
