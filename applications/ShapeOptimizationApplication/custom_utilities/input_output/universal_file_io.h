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
#include <fstream>
#include <string>
#include <algorithm>
#include <iomanip>      // for std::setprecision

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

class UniversalFileIO
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UniversalFileIO( ModelPart& OutputModelPart, std::string OutputFilenameWithoutExtension, std::string WriteConditionsFlag, Parameters NodalResults )
        : mrOutputModelPart( OutputModelPart ),
          mOutputFilenameWithExtension( OutputFilenameWithoutExtension + ".unv" ),
          mrNodalResults( NodalResults )
    {
        if(WriteConditionsFlag.compare("WriteElementsOnly")==0 || WriteConditionsFlag.compare("WriteConditionsOnly")==0 )
            mWriteConditionsFlag = WriteConditionsFlag;
        else
            KRATOS_ERROR << "Wrong value specified for \"WriteConditionsFlag\" in UniversalFileIO. Options are: \"WriteElementsOnly\" or \"WriteConditionsOnly\"" << std::endl;
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
    void InitializeLogging()
    {
        WriteMeshToResultFile();
    }

    // --------------------------------------------------------------------------
    void WriteMeshToResultFile()
    {
        InitializeOutputFile();
        WriteUnits();
        WriteNodes();
        WriteElements();
    }

    // --------------------------------------------------------------------------
    void InitializeOutputFile()
    {
        std::ofstream outputFile;
        outputFile.open(mOutputFilenameWithExtension, std::ios::out | std::ios::trunc );
        outputFile.close();
    }
    // --------------------------------------------------------------------------
    void WriteUnits()
    {
        std::ofstream outputFile;
        outputFile.open(mOutputFilenameWithExtension, std::ios::out | std::ios::app );
        outputFile << std::scientific;
        outputFile << std::setprecision(15);

        const int dataSetNumberForNodes = 164;
        const int unitCode = 5;
        const int temperatureMode = 2;
        const double unitFactorLenght = 1.0;
        const double unitFactorForce = 1.0;
        const double unitFactorTemperature = 1.0;
        const double temperatureOffset = 273.15;

        outputFile << std::setw(6) << "-1" << "\n";
        outputFile << std::setw(6) << dataSetNumberForNodes << "\n";
        outputFile << std::setw(10) << unitCode << std::setw(30) << temperatureMode << "\n";
        outputFile << std::setw(25) << unitFactorLenght << std::setw(25) <<  unitFactorForce << std::setw(25) << unitFactorTemperature << "\n";
        outputFile << std::setw(25) << temperatureOffset << "\n";
        outputFile << std::setw(6) << "-1" << "\n";

        outputFile.close();
    }

    // --------------------------------------------------------------------------
    void WriteNodes()
    {
        std::ofstream outputFile;
        outputFile.open(mOutputFilenameWithExtension, std::ios::out | std::ios::app );
        outputFile << std::scientific;
        outputFile << std::setprecision(15);

        const int dataSetNumberForNodes = 2411;
        const int exportCoordinateSystemNumber = 0;
        const int displacementCoordinateSystemNumber = 0;
        const int color = 0;

        outputFile << std::setw(6) << "-1" << "\n";
        outputFile << std::setw(6) << dataSetNumberForNodes << "\n";
        for (auto & node_i : mrOutputModelPart.Nodes())
        {
            int node_label = node_i.Id();
            double x_coordinate = node_i.X();
            double y_coordinate = node_i.Y();
            double z_coordinate = node_i.Z();
            outputFile << std::setw(10) << node_label << std::setw(10) << exportCoordinateSystemNumber << std::setw(10) << displacementCoordinateSystemNumber << std::setw(10) << color << "\n";
            outputFile << std::setw(25) << x_coordinate << std::setw(25) <<  y_coordinate << std::setw(25) << z_coordinate << "\n";
        }
        outputFile << std::setw(6) << "-1" << "\n";

        outputFile.close();
    }

    // --------------------------------------------------------------------------
    void WriteElements()
    {
        if(mWriteConditionsFlag.compare("WriteElementsOnly")==0)
            WriteAllElementsButNoConditions();
        else if(mWriteConditionsFlag.compare("WriteConditionsOnly")==0)
            WriteConditionsAsDummyElements();
    }

    // --------------------------------------------------------------------------
    void WriteAllElementsButNoConditions()
    {
        std::ofstream outputFile;
        outputFile.open(mOutputFilenameWithExtension, std::ios::out | std::ios::app );

        const int dataSetNumberForElements = 2412;
        const int physicalPropertyTableNumber = 1;
        const int materialPropertyTableNumber = 1;
        const int color = 0;

        outputFile << std::setw(6) << "-1" << "\n";
        outputFile << std::setw(6) << dataSetNumberForElements << "\n";

        for (auto & element_i : mrOutputModelPart.Elements())
        {
            const int element_label =element_i.Id();
            ModelPart::ConditionType::GeometryType element_geometry = element_i.GetGeometry();

            // Write triangles
            if( element_geometry.size()==3 && element_geometry.Dimension()==2 )
            {
                const int feDescriptorId = 41; // Plane Stress Linear Triangle
                const int numberOfNodes = 3;
                outputFile << std::setw(10) << element_label;
                outputFile << std::setw(10) << feDescriptorId;
                outputFile << std::setw(10) << physicalPropertyTableNumber;
                outputFile << std::setw(10) << materialPropertyTableNumber ;
                outputFile << std::setw(10) << color;
                outputFile << std::setw(10) << numberOfNodes << "\n";
                outputFile << std::setw(10) << element_geometry[0].Id();
                outputFile << std::setw(10) << element_geometry[1].Id();
                outputFile << std::setw(10) << element_geometry[2].Id() << "\n";
            }
            // Write quads
            else if( element_geometry.size()==4 && element_geometry.Dimension()==2 )
            {
                const int feDescriptorId = 44; // Plane Stress Linear Quadrilateral
                const int numberOfNodes = 4;
                outputFile << std::setw(10) << element_label;
                outputFile << std::setw(10) << feDescriptorId;
                outputFile << std::setw(10) << physicalPropertyTableNumber;
                outputFile << std::setw(10) << materialPropertyTableNumber;
                outputFile << std::setw(10) << color;
                outputFile << std::setw(10) << numberOfNodes << "\n";
                outputFile << std::setw(10) << element_geometry[0].Id();
                outputFile << std::setw(10) << element_geometry[1].Id();
                outputFile << std::setw(10) << element_geometry[2].Id();
                outputFile << std::setw(10) << element_geometry[3].Id() << "\n";
            }
            // Write tetrahedras
            else if( element_geometry.size()==4 && element_geometry.Dimension()==3 )
            {
                const int feDescriptorId = 111; // Solid linear tetrahedron
                const int numberOfNodes = 4;
                outputFile << std::setw(10) << element_label;
                outputFile << std::setw(10) << feDescriptorId;
                outputFile << std::setw(10) << physicalPropertyTableNumber;
                outputFile << std::setw(10) << materialPropertyTableNumber;
                outputFile << std::setw(10) << color;
                outputFile << std::setw(10) << numberOfNodes << "\n";
                outputFile << std::setw(10) << element_geometry[0].Id();
                outputFile << std::setw(10) << element_geometry[1].Id();
                outputFile << std::setw(10) << element_geometry[2].Id();
                outputFile << std::setw(10) << element_geometry[3].Id() << "\n";
            }
            else
                KRATOS_ERROR << "Design surface contains elements with geometries for which no UNV-output is implemented!" << std::endl;
        }

        outputFile << std::setw(6) << "-1" << "\n";
        outputFile.close();
    }

    // --------------------------------------------------------------------------
    void WriteConditionsAsDummyElements()
    {
        std::ofstream outputFile;
        outputFile.open(mOutputFilenameWithExtension, std::ios::out | std::ios::app );

        const int dataSetNumberForElements = 2412;
        const int physicalPropertyTableNumber = 1;
        const int materialPropertyTableNumber = 1;
        const int color = 0;

        outputFile << std::setw(6) << "-1" << "\n";
        outputFile << std::setw(6) << dataSetNumberForElements << "\n";

        for (auto & condition_i : mrOutputModelPart.Conditions())
        {
            const int element_label = condition_i.Id();
            ModelPart::ConditionType::GeometryType condition_geometry = condition_i.GetGeometry();

            if( condition_geometry.size()==3 )
            {
                const int feDescriptorId = 41; // Plane Stress Linear Triangle
                const int numberOfNodes = 3;
                outputFile << std::setw(10) << element_label;
                outputFile << std::setw(10) << feDescriptorId;
                outputFile << std::setw(10) << physicalPropertyTableNumber;
                outputFile << std::setw(10) << materialPropertyTableNumber ;
                outputFile << std::setw(10) << color;
                outputFile << std::setw(10) << numberOfNodes << "\n";
                outputFile << std::setw(10) << condition_geometry[0].Id();
                outputFile << std::setw(10) << condition_geometry[1].Id();
                outputFile << std::setw(10) << condition_geometry[2].Id() << "\n";
            }
            else if( condition_geometry.size()==4 )
            {
                const int feDescriptorId = 44; // Plane Stress Linear Quadrilateral
                const int numberOfNodes = 4;
                outputFile << std::setw(10) << element_label;
                outputFile << std::setw(10) << feDescriptorId;
                outputFile << std::setw(10) << physicalPropertyTableNumber;
                outputFile << std::setw(10) << materialPropertyTableNumber;
                outputFile << std::setw(10) << color;
                outputFile << std::setw(10) << numberOfNodes << "\n";
                outputFile << std::setw(10) << condition_geometry[0].Id();
                outputFile << std::setw(10) << condition_geometry[1].Id();
                outputFile << std::setw(10) << condition_geometry[2].Id();
                outputFile << std::setw(10) << condition_geometry[3].Id() << "\n";
            }
            else
                KRATOS_ERROR << "Design surface contains conditions with geometries for which no UNV-output is implemented!" << std::endl;
        }

        outputFile << std::setw(6) << "-1" << "\n";
        outputFile.close();
    }

    // --------------------------------------------------------------------------
    void LogNodalResults( const int optimizationIteration )
    {
        std::ofstream outputFile;
        outputFile.open(mOutputFilenameWithExtension, std::ios::out | std::ios::app );
        outputFile << std::scientific;
        outputFile << std::setprecision(5);

        const int dataSetNumberForResults = 2414;

        const int analysisDataSetLabel = 1;
        const int dataSetLocation = 1; // Data at nodes
        std::string analysisDataSetName = "KratosShapeOptimization";
        const int expressionNameKey = 1;

        const int modelType = 1; // Structural
        const int analysisType = 4; // Transient analysis
        const int resultType = 95; // Unknown 3DOF Vector
        const int displacementResultType = 8; // Displacement
        const int dataType = 4; // Double precision floating point
        const int numberOfDataValuesForTheDataComponent = 3;

        const int DesignSetId = 1;
        const int iterationNumber = 0;
        const int solutionSetId = 1;
        const int boundaryCondition = 0;
        const int loadSet = 1;
        const int modeNumber = 0;
        const int timeStepNumber = optimizationIteration;
        const int frequencyNumber = 0;

        const int creationOption = 0;
        const int numberRetained = 0;

        const double givenTime = 0.0;
        const double frequency = 0.0;
        const double eigenvalue = 0.0;
        const double modalMass = 0.0;
        const double viscousDampingRation = 0.0;
        const double hystereticDampingRatio = 0.0;
        const double realPartEigenvalue = 0.0;
        const double imaginaryPartEigenvalue = 0.0;
        const double realPartOfModalA = 0.0;
        const double realPartOfMass = 0.0;
        const double imaginaryPartOfModalA = 0.0;
        const double imaginaryPartOfMass = 0.0;

        // Loop over all nodal result variables
        for(unsigned int entry = 0; entry<mrNodalResults.size(); entry++)
        {
            std::string nodalResultName = mrNodalResults[entry].GetString();

            unsigned int dataCharacteristic = 0; // 0: unknown, 1: Scalar value, 2: 3 DOF global translation vector
            if( KratosComponents<Variable<double>>::Has(nodalResultName))
                dataCharacteristic = 1;
            else if( KratosComponents<Variable< array_1d<double,3>>>::Has(nodalResultName))
                dataCharacteristic = 2;

            outputFile << std::setw(6) << "-1" << "\n";
            outputFile << std::setw(6) << dataSetNumberForResults << "\n";

            outputFile << std::setw(10) << analysisDataSetLabel << "\n";
            outputFile << "LOADCASE_NAME_KEY " << analysisDataSetName << "\n";
            outputFile << std::setw(10) << dataSetLocation << "\n";
            outputFile << "RESULT_NAME_KEY " << nodalResultName << "\n";
            outputFile << "NONE" << "\n";
            outputFile << "EXPRESSION_NAME_KEY " << expressionNameKey << "\n";
            outputFile << "NONE" << "\n";
            outputFile << "NONE" << "\n";

            outputFile << std::setw(10) << modelType;
            outputFile << std::setw(10) << analysisType;
            outputFile << std::setw(10) << dataCharacteristic;
            if (nodalResultName.compare("SHAPE_CHANGE") == 0 || nodalResultName.compare("MESH_CHANGE") == 0)
                outputFile << std::setw(10) << displacementResultType;
            else
                outputFile << std::setw(10) << resultType;
            outputFile << std::setw(10) << dataType;
            outputFile << std::setw(10) << numberOfDataValuesForTheDataComponent << "\n";

            outputFile << std::setw(10) << DesignSetId;
            outputFile << std::setw(10) << iterationNumber;
            outputFile << std::setw(10) << solutionSetId;
            outputFile << std::setw(10) << boundaryCondition;
            outputFile << std::setw(10) << loadSet;
            outputFile << std::setw(10) << modeNumber;
            outputFile << std::setw(10) << timeStepNumber;
            outputFile << std::setw(10) << frequencyNumber << "\n";

            outputFile << std::setw(10) << creationOption;
            outputFile << std::setw(10) << numberRetained << "\n";

            outputFile << std::setw(13) << givenTime;
            outputFile << std::setw(13) << frequency;
            outputFile << std::setw(13) << eigenvalue;
            outputFile << std::setw(13) << modalMass;
            outputFile << std::setw(13) << viscousDampingRation;
            outputFile << std::setw(13) << hystereticDampingRatio << "\n";

            outputFile << std::setw(13) << realPartEigenvalue;
            outputFile << std::setw(13) << imaginaryPartEigenvalue;
            outputFile << std::setw(13) << realPartOfModalA;
            outputFile << std::setw(13) << realPartOfMass;
            outputFile << std::setw(13) << imaginaryPartOfModalA;
            outputFile << std::setw(13) << imaginaryPartOfMass << "\n";

            for (auto & node_i : mrOutputModelPart.Nodes())
            {
                outputFile << std::setw(10) << node_i.Id() << "\n";
                if(dataCharacteristic==1)
                {
                    const Variable<double>& nodalResultVariable = KratosComponents<Variable<double>>::Get(nodalResultName);
                    double& nodalResult = node_i.FastGetSolutionStepValue(nodalResultVariable);
                    outputFile << std::setw(13) << nodalResult << "\n";
                }
                else if(dataCharacteristic==2)
                {
                    const Variable< array_1d<double,3>>& nodalResultVariable = KratosComponents<Variable< array_1d<double,3>>>::Get(nodalResultName);
                    array_1d<double,3>& nodalResult = node_i.FastGetSolutionStepValue(nodalResultVariable);
                    outputFile << std::setw(13) << nodalResult[0];
                    outputFile << std::setw(13) << nodalResult[1];
                    outputFile << std::setw(13) << nodalResult[2] << "\n";
                }
            }
            outputFile << std::setw(6) << "-1" << "\n";
        }

        outputFile.close();
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
