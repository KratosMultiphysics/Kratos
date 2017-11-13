// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Sferza Massimo, https://github.com/IIIaxS
//
// ==============================================================================

#ifndef VTK_FILE_IO_H
#define VTK_FILE_IO_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
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

*/

class VTKFileIO
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VTKFileIO( ModelPart& OutputModelPart, std::string OutputFilenamePrefix, std::string WriteConditionsFlag, Parameters NodalResults )
        : mrOutputModelPart( OutputModelPart ),
          mOutputFilenameWithoutExtension( OutputFilenamePrefix ),
          mrNodalResults( NodalResults )
    {
        mDefaultPrecision = 15;
        if(WriteConditionsFlag.compare("WriteElementsOnly")==0 || WriteConditionsFlag.compare("WriteConditionsOnly")==0 )
            mWriteConditionsFlag = WriteConditionsFlag;        
        else
            KRATOS_ERROR << "Wrong value specified for \"WriteConditionsFlag\" in VTKFileIO. Options are: \"WriteElementsOnly\" or \"WriteConditionsOnly\"" << std::endl;
    }

    /// Destructor.
    virtual ~VTKFileIO()
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
        mKratosIdToVtkId = CreateMapFromKratosIdToVTKId();
        mVtkCellListSize = DetermineVtkCellListSize();
    }

    // --------------------------------------------------------------------------
    std::map<int,int> CreateMapFromKratosIdToVTKId()
    {
        std::map<int,int> kratos_id_to_vtk;
        int vtk_id = 0;

        for(auto & node_i : mrOutputModelPart.Nodes())
        {
            int KratosId = node_i.Id();
            kratos_id_to_vtk[KratosId] = vtk_id;
            vtk_id++;
        }

        return kratos_id_to_vtk;
    }

    // --------------------------------------------------------------------------
    unsigned int DetermineVtkCellListSize()
    {
         unsigned int vtk_cell_list_size = 0;

         if(mWriteConditionsFlag.compare("WriteElementsOnly")==0)
         {
            for (auto & element_i : mrOutputModelPart.Elements())
            {
                vtk_cell_list_size++;
                vtk_cell_list_size += element_i.GetGeometry().size();
            }
         }
         else if(mWriteConditionsFlag.compare("WriteConditionsOnly")==0)
         {
             for (auto & condition_i : mrOutputModelPart.Conditions())
             {
                 vtk_cell_list_size++;
                 vtk_cell_list_size += condition_i.GetGeometry().size();
             }
         }

        return vtk_cell_list_size;
    }

    // --------------------------------------------------------------------------
    void LogNodalResults( const int optimizationIteration )
    {
        UpdateOutputFilename( optimizationIteration );
        WriteHeader();
        WriteMesh();
        WriteNodalResults();
    }

    // --------------------------------------------------------------------------
    void UpdateOutputFilename( const int optimizationIteration )
    {
        std::string outputFilename = mOutputFilenameWithoutExtension + "_" + std::to_string(optimizationIteration) + ".vtk";
        mOutputFilenameWithExtension = outputFilename;
    }

    // --------------------------------------------------------------------------
    void WriteHeader()
    {
        std::ofstream outputFile;
        outputFile.open(mOutputFilenameWithExtension, std::ios::out | std::ios::trunc );
        outputFile << "# vtk DataFile Version 4.0" << "\n";
        outputFile << "vtk output" << "\n";
        outputFile << "ASCII" << "\n";
        outputFile << "DATASET UNSTRUCTURED_GRID" << "\n";
        outputFile.close();
    }

    // --------------------------------------------------------------------------
    void WriteMesh()
    {
        WriteNodes();
        WriteElements();
        WriteElementTypes();
    }

    // --------------------------------------------------------------------------
    void WriteNodalResults()
    {
        WriteFirstNodalResultsAsPointData();
        WriteOtherNodalResultsAsFieldData();
    }

    // --------------------------------------------------------------------------
    void WriteNodes()
    {
        std::ofstream outputFile;
        outputFile.open(mOutputFilenameWithExtension, std::ios::out | std::ios::app );
        outputFile << std::scientific;
        outputFile << std::setprecision(mDefaultPrecision);

        outputFile << "POINTS " << mrOutputModelPart.NumberOfNodes() << " float" << "\n";
        for (auto & node_i : mrOutputModelPart.Nodes())
        {
            double x_coordinate = node_i.X0();
            double y_coordinate = node_i.Y0();
            double z_coordinate = node_i.Z0();
            outputFile << " " << x_coordinate;
            outputFile << " " << y_coordinate;
            outputFile << " " << z_coordinate << "\n";
        }

        outputFile.close();
    }

    // --------------------------------------------------------------------------
    void WriteElements()
    {
        if(mWriteConditionsFlag.compare("WriteElementsOnly")==0)
            WriteAllElementsButNoConditions();
        else if(mWriteConditionsFlag.compare("WriteConditionsOnly")==0)
            WriteConditionsAsElements();
    }

    // --------------------------------------------------------------------------
    void WriteAllElementsButNoConditions()
    {
        std::ofstream outputFile;
        outputFile.open(mOutputFilenameWithExtension, std::ios::out | std::ios::app );

        outputFile << "CELLS " << mrOutputModelPart.NumberOfElements() << " " << mVtkCellListSize << "\n";
        for (auto & element_i : mrOutputModelPart.Elements())
        {
            ModelPart::ElementType::GeometryType& element_geometry = element_i.GetGeometry();
            const unsigned int numberOfNodes = element_geometry.size();

            outputFile << numberOfNodes;
            for (unsigned int i=0; i<numberOfNodes; i++)
                outputFile << " " << mKratosIdToVtkId[element_geometry[i].Id()];
            outputFile << "\n";
        }

        outputFile.close();               
    }

    // --------------------------------------------------------------------------
    void WriteConditionsAsElements()
    {
        std::ofstream outputFile;
        outputFile.open(mOutputFilenameWithExtension, std::ios::out | std::ios::app );

        outputFile << "CELLS " << mrOutputModelPart.NumberOfConditions() << " " << mVtkCellListSize << "\n";
        for (auto & condition_i : mrOutputModelPart.Conditions())
        {
            ModelPart::ConditionType::GeometryType& condition_geometry = condition_i.GetGeometry();
            const unsigned int numberOfNodes = condition_geometry.size();

            outputFile << numberOfNodes;
            for (unsigned int i=0; i<numberOfNodes; i++)
                outputFile << " " << mKratosIdToVtkId[condition_geometry[i].Id()];
            outputFile << "\n";
        }

        outputFile.close();        
    }

    // --------------------------------------------------------------------------
    void WriteElementTypes()
    {
        if(mWriteConditionsFlag.compare("WriteElementsOnly")==0)
            WriteElementTypesOfElementsOnly();
        else if(mWriteConditionsFlag.compare("WriteConditionsOnly")==0)
            WriteElementTypesOfConditionsOnly();
    }

    // --------------------------------------------------------------------------
    void WriteElementTypesOfElementsOnly()
    {
        std::ofstream outputFile;
        outputFile.open(mOutputFilenameWithExtension, std::ios::out | std::ios::app );

        outputFile << "CELL_TYPES " << mrOutputModelPart.NumberOfElements() << "\n";
        for (auto & element_i : mrOutputModelPart.Elements())
        {
            const unsigned int numberOfNodes =  element_i.GetGeometry().size();
            const unsigned int dimension = element_i.GetGeometry().Dimension();              
            
            unsigned int vtk_cell_type;
            if( numberOfNodes == 3  && dimension == 2) // triangle
                vtk_cell_type = 5;
            else if( numberOfNodes == 4  && dimension == 2) // quad
                vtk_cell_type = 9;
            else if( numberOfNodes == 4  && dimension == 3) // tet
                vtk_cell_type = 10;
            else if( numberOfNodes == 8 && dimension == 3) // hex
                vtk_cell_type = 12;                                        
            else
                KRATOS_ERROR << "Optimization model part contains elements with geometries for which no VTK-output is implemented!" << std::endl;

            outputFile << vtk_cell_type << "\n";
        }

        outputFile.close();
    }

    // --------------------------------------------------------------------------
    void WriteElementTypesOfConditionsOnly()
    { 
        std::ofstream outputFile;
        outputFile.open(mOutputFilenameWithExtension, std::ios::out | std::ios::app );

        outputFile << "CELL_TYPES " << mrOutputModelPart.NumberOfConditions() << "\n";
        for (auto & condition_i : mrOutputModelPart.Conditions())
        {
            const unsigned int numberOfNodes =  condition_i.GetGeometry().size();
            const unsigned int dimension = condition_i.GetGeometry().Dimension();
            
            unsigned int vtk_cell_type;
            if( numberOfNodes == 3 && dimension == 2) // triangle
                vtk_cell_type = 5;
            else if( numberOfNodes == 4  && dimension == 2) // quad
                vtk_cell_type = 9;
            else
                KRATOS_ERROR << "Design surface contains conditions with geometries for which no VTK-output is implemented!" << std::endl;

            outputFile << vtk_cell_type << "\n";
        }

        outputFile.close();
    }

    // --------------------------------------------------------------------------
    void WriteFirstNodalResultsAsPointData()
    {
        std::ofstream outputFile;
        outputFile.open(mOutputFilenameWithExtension, std::ios::out | std::ios::app );

        // Write nodal results header
        outputFile << "POINT_DATA " << mrOutputModelPart.NumberOfNodes() << "\n";

        // Write nodal results variable header
        std::string nodalResultName = mrNodalResults[0].GetString();
        unsigned int dataCharacteristic = 0; // 0: unknown, 1: Scalar value, 2: 3 DOF global translation vector
        if( KratosComponents<Variable<double>>::Has(nodalResultName))
        {
            dataCharacteristic = 1;
            outputFile << "SCALARS " << nodalResultName << " float" << "\n";
        }
        else if( KratosComponents<Variable< array_1d<double,3>>>::Has(nodalResultName))
        {
            dataCharacteristic = 2;
            outputFile << "VECTORS " << nodalResultName << " float" << "\n";
        }

        // Write nodal results
        outputFile << std::scientific;
        outputFile << std::setprecision(mDefaultPrecision);
        for (ModelPart::NodeIterator node_i = mrOutputModelPart.NodesBegin(); node_i != mrOutputModelPart.NodesEnd(); ++node_i)
        {
            if(dataCharacteristic==1)
            {
                Variable<double> nodalResultVariable = KratosComponents<Variable<double>>::Get(nodalResultName);
                double& nodalResult = node_i->FastGetSolutionStepValue(nodalResultVariable);
                outputFile << nodalResult << "\n";
            }
            else if(dataCharacteristic==2)
            {
                Variable< array_1d<double,3>> nodalResultVariable = KratosComponents<Variable< array_1d<double,3>>>::Get(nodalResultName);
                array_1d<double,3>& nodalResult = node_i->FastGetSolutionStepValue(nodalResultVariable);
                outputFile << nodalResult[0] << " ";
                outputFile << nodalResult[1] << " ";
                outputFile << nodalResult[2] << "\n";
            }
        }

        outputFile.close();
    }

    // --------------------------------------------------------------------------
    void WriteOtherNodalResultsAsFieldData()
    {
        std::ofstream outputFile;
        outputFile.open(mOutputFilenameWithExtension, std::ios::out | std::ios::app );

        // Write nodal results header
        outputFile << "FIELD FieldData " << mrNodalResults.size()-1 << "\n";

        for(unsigned int entry = 1; entry < mrNodalResults.size(); entry++)
        {
            // Write nodal results variable header
            std::string nodalResultName = mrNodalResults[entry].GetString();
            unsigned int dataCharacteristic = 0; // 0: unknown, 1: Scalar value, 2: 3 DOF global translation vector
            if( KratosComponents<Variable<double>>::Has(nodalResultName))
            {
                dataCharacteristic = 1;
                outputFile << nodalResultName << " 1 " << mrOutputModelPart.NumberOfNodes() << " float" << "\n";
            }
            else if( KratosComponents<Variable< array_1d<double,3>>>::Has(nodalResultName))
            {
                dataCharacteristic = 2;
                outputFile << nodalResultName << " 3 " << mrOutputModelPart.NumberOfNodes() << " float" << "\n";
            }

            // Write nodal results
            outputFile << std::scientific;
            outputFile << std::setprecision(mDefaultPrecision);
            for (ModelPart::NodeIterator node_i = mrOutputModelPart.NodesBegin(); node_i != mrOutputModelPart.NodesEnd(); ++node_i)
            {
                if(dataCharacteristic==1)
                {
                    Variable<double> nodalResultVariable = KratosComponents<Variable<double>>::Get(nodalResultName);
                    double& nodalResult = node_i->FastGetSolutionStepValue(nodalResultVariable);
                    outputFile << nodalResult << "\n";
                }
                else if(dataCharacteristic==2)
                {
                    Variable< array_1d<double,3>> nodalResultVariable = KratosComponents<Variable< array_1d<double,3>>>::Get(nodalResultName);
                    array_1d<double,3>& nodalResult = node_i->FastGetSolutionStepValue(nodalResultVariable);
                    outputFile << nodalResult[0] << " ";
                    outputFile << nodalResult[1] << " ";
                    outputFile << nodalResult[2] << "\n";
                }
            }
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
        return "VTKFileIO";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "VTKFileIO";
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
    std::string mOutputFilenameWithoutExtension;
    Parameters mrNodalResults; 
    unsigned int mDefaultPrecision;
    std::string mWriteConditionsFlag;   
    std::map<int,int> mKratosIdToVtkId;
    unsigned int mVtkCellListSize;
    std::string mOutputFilenameWithExtension;

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
//      VTKFileIO& operator=(VTKFileIO const& rOther);

    /// Copy constructor.
//      VTKFileIO(VTKFileIO const& rOther);


    ///@}

}; // Class VTKFileIO

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // VTK_FILE_IO_H
