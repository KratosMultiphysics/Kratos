//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_INPUT_PRIMAL_SOLUTION_PROCESS_H_INCLUDED )
#define KRATOS_INPUT_PRIMAL_SOLUTION_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>

// External includes
#include "H5Cpp.h"

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "includes/model_part.h"

namespace Kratos
{

///@addtogroup AdjointFluidApplication
///@{

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

/// Process to input primal solution (e.g., velocity and pressure) needed for adjoint solution.
/**
 * Requires mesh topology and node numbering to remain constant.
 */
class InputPrimalSolutionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(InputPrimalSolutionProcess);

    typedef ModelPart::NodeType NodeType;

    typedef hsize_t SizeType;

    // typedef typename TDenseSpace::VectorType DenseVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    InputPrimalSolutionProcess(ModelPart& rModelPart, Parameters& rParameters)
    : Process(), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        Parameters DefaultParams(R"(
        {
            "model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "file_name": "PLEASE_SPECIFY_H5_FILE_NAME",
            "variable_list": ["VELOCITY", "ACCELERATION", "PRESSURE"],
            "nodal_value_list": [],
            "processinfo_value_list": []
        })");

        // try accessing parameters without defaults so that an error is thrown
        // if they don't exist
        rParameters["file_name"];

        rParameters.ValidateAndAssignDefaults(DefaultParams);

        // input filename
        std::stringstream FilenameStream;
        FilenameStream << rParameters["file_name"].GetString() << '_' << mrModelPart.GetCommunicator().MyPID() << ".h5";
        mFilename = FilenameStream.str();

        // nodal variable names to input
        mVariableNames.resize(rParameters["variable_list"].size());
        for (unsigned int i = 0; i < mVariableNames.size(); i++)
        {
            mVariableNames[i] = rParameters["variable_list"].GetArrayItem(i).GetString();
        }

        // nodal value names to intput
        mNodalValueNames.resize(rParameters["nodal_value_list"].size());
        for( unsigned int i=0; i < mNodalValueNames.size(); ++i )
        {
            mNodalValueNames[i] = rParameters["nodal_value_list"].GetArrayItem(i).GetString();
        }

        // nodal value names to intput
        mProcessInfoValueNames.resize(rParameters["processinfo_value_list"].size());
        for( unsigned int i=0; i < mProcessInfoValueNames.size(); ++i )
        {
            mProcessInfoValueNames[i] = rParameters["processinfo_value_list"].GetArrayItem(i).GetString();
        }

        mNumNodes = 0;

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~InputPrimalSolutionProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    /// This operator simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    void Execute() override {}

    void ExecuteInitialize() override
    {
        KRATOS_TRY

        // test for valid input parameters / model part
        for (unsigned int i = 0; i < mVariableNames.size(); i++)
        {
            if (KratosComponents< Variable<double> >::Has(mVariableNames[i]))
            {
                if (mrModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents< Variable<double> >::Get(mVariableNames[i])) == false)
                {
                    KRATOS_THROW_ERROR(std::runtime_error, "variable is not found in nodal solution steps variable list: ", mVariableNames[i])
                }
            }
            else if (KratosComponents< Variable<array_1d<double,3> > >::Has(mVariableNames[i]))
            {
                if (mrModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents< Variable<array_1d<double,3> > >::Get(mVariableNames[i])) == false)
                {
                    KRATOS_THROW_ERROR(std::runtime_error,"variable is not found in nodal solution steps variable list: ", mVariableNames[i])
                }
            }
            else
            {
                KRATOS_THROW_ERROR(std::invalid_argument, "variable type not supported: ", mVariableNames[i])
            }
        }

        // nodes should be ordered by node id before IO begins
        if (mrModelPart.Nodes().IsSorted() == false)
        {
            mrModelPart.Nodes().Sort();
        }

        // number of nodes must be constant for IO
        mNumNodes = static_cast<SizeType>(mrModelPart.Nodes().size());

        KRATOS_CATCH("")
    }

    void ExecuteBeforeSolutionLoop() override
    {
        KRATOS_TRY

        // read the values from the file before starting the computation
        // open file to read data
        H5::H5File File(mFilename.c_str(), H5F_ACC_RDONLY);
        
        H5::DataSet NodeIdDataset = File.openDataSet("/NodalData/Id");

        hsize_t dims[3];
        NodeIdDataset.getSpace().getSimpleExtentDims(dims);

        if (mNumNodes != dims[0])
            KRATOS_THROW_ERROR(std::runtime_error, "inconsistent dimension in file: ", mFilename)

        // if (mNumNodes != static_cast<SizeType>(mrModelPart.Nodes().size()))
        if (mNumNodes != mrModelPart.Nodes().size())
            KRATOS_THROW_ERROR(std::runtime_error, "detected change in number of nodes.", "")

        std::vector<unsigned int> NodeIdBuffer(mrModelPart.Nodes().size());
        NodeIdDataset.read(NodeIdBuffer.data(), H5::PredType::NATIVE_UINT);

        // check that the node ordering is consistent between the model part and file
        int Delta = 0;
        int itNode = 0;
        for (const auto& rNode : mrModelPart.Nodes())
        {
            Delta += std::abs(rNode.Id() - NodeIdBuffer[itNode]);
            itNode++;
        }

        if (Delta != 0)
            KRATOS_ERROR_IF(Delta != 0) << "detected mismatch of node ids in file: " << mFilename << std::endl;

        std::string values_group_name = "/NodalData/Values";

        for( std::size_t i = 0; i < mNodalValueNames.size(); ++i)
        {
            std::string dataset_name = values_group_name + "/" + mNodalValueNames[i];

            // File.getObjinfo(dataset_name.c_str(),&status);
            H5::DataSet value_dataset = File.openDataSet(dataset_name.c_str());
            value_dataset.getSpace().getSimpleExtentDims(dims);

            if (KratosComponents< Variable<Vector> >::Has(mNodalValueNames[i]))
            {
                std::vector<double> value_data_buffer(dims[0] * dims[1]);
                value_dataset.read(value_data_buffer.data(), H5::PredType::NATIVE_DOUBLE);
                auto it_value = std::begin(value_data_buffer);

                for( auto& rNode : mrModelPart.Nodes() )
                {
                    Vector temp_vector( dims[1] );
                    for( size_t j = 0; j < dims[1]; ++j)
                    {
                        temp_vector[j] = *it_value;
                        it_value++;
                    }
                    rNode.GetValue(KratosComponents< Variable<Vector> >::Get(mNodalValueNames[i])) = temp_vector;
                }
            }

            else if( KratosComponents< Variable<Matrix> >::Has(mNodalValueNames[i]) )
            {
                std::vector<double> value_data_buffer(dims[0] * dims[1] * dims[2]);
                value_dataset.read(value_data_buffer.data(), H5::PredType::NATIVE_DOUBLE);
                auto it_value = std::begin(value_data_buffer);

                std::for_each(std::begin(mrModelPart.Nodes()), std::end(mrModelPart.Nodes()),
                   [&](NodeType& rNode)
                {
                    Matrix temp_matrix( dims[2], dims[1] );
                    for( size_t j = 0; j < dims[1]; ++j )
                    {
                        for( size_t k = 0; k < dims[2]; ++k )
                        {
                            temp_matrix(k,j) = *it_value;
                            it_value++;
                        }
                    }
                    rNode.GetValue(KratosComponents< Variable<Matrix> >::Get(mNodalValueNames[i])) = temp_matrix;
                });
            }

            else
            {
                KRATOS_ERROR << "Reading from value " << mNodalValueNames[i] << " is not implemented right now" << std::endl;
            }
        } //nodal values
        
        std::string procinfo_group_name("/ProcessInfo");
        hsize_t procinfo_dims[1];
        for( std::size_t i = 0; i < mProcessInfoValueNames.size(); ++i)
        {
            std::string dataset_name = procinfo_group_name + "/" + mProcessInfoValueNames[i];
            
            H5::DataSet procinfo_dataset = File.openDataSet(dataset_name.c_str());
            procinfo_dataset.getSpace().getSimpleExtentDims(procinfo_dims);

            if (KratosComponents< Variable<Vector> >::Has(mProcessInfoValueNames[i]))
            {
                const Variable< Vector >& rVariable = KratosComponents< Variable<Vector> >::Get(mProcessInfoValueNames[i]);
                
                std::vector<double> procinfo_data_buffer(procinfo_dims[0]); 
                procinfo_dataset.read(procinfo_data_buffer.data(), H5::PredType::NATIVE_DOUBLE);
                
                Vector temp_vector( procinfo_data_buffer.size() );
                for( size_t j = 0; j < procinfo_data_buffer.size(); ++j )
                {
                    temp_vector[j] = procinfo_data_buffer[j];
                }
                mrModelPart.GetProcessInfo()[rVariable] = temp_vector;
            }
        } //processinfo values

        KRATOS_CATCH("")
    }

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY

        // open file to read data
        H5::H5File File(mFilename.c_str(), H5F_ACC_RDONLY);

        H5::DataSet NodeIdDataset = File.openDataSet("/NodalData/Id");

        hsize_t Dims[2];
        NodeIdDataset.getSpace().getSimpleExtentDims(Dims);

        if (mNumNodes != Dims[0])
            KRATOS_THROW_ERROR(std::runtime_error, "inconsistent dimension in file: ", mFilename)

        if (mNumNodes != mrModelPart.Nodes().size())
            KRATOS_THROW_ERROR(std::runtime_error, "detected change in number of nodes.", "")

        std::vector<unsigned int> NodeIdBuffer(mrModelPart.Nodes().size());
        NodeIdDataset.read(NodeIdBuffer.data(), H5::PredType::NATIVE_UINT);

        // check that the node ordering is consistent between the model part and file
        int Delta = 0;
        int itNode = 0;
        for (const auto& rNode : mrModelPart.Nodes())
        {
            Delta += std::abs(rNode.Id() - NodeIdBuffer[itNode]);
            itNode++;
        }

        if (Delta != 0)
            KRATOS_ERROR_IF(Delta != 0) << "detected mismatch of node ids in file: " << mFilename << std::endl;

        // input time step data
        std::stringstream TimeStepPathStream;
        TimeStepPathStream << "/NodalData/" << std::fixed << std::setprecision(10) << mrModelPart.GetProcessInfo()[TIME];
        std::string TimeStepGroupName = TimeStepPathStream.str();

        for (unsigned int i = 0; i < mVariableNames.size(); i++)
        {
            std::string DatasetName = TimeStepGroupName + "/" + mVariableNames[i];
            H5::DataSet VariableDataset = File.openDataSet(DatasetName.c_str());
            if (KratosComponents< Variable<double> >::Has(mVariableNames[i]))
            {
                std::vector<double> VariableDataBuffer(mrModelPart.Nodes().size());
                VariableDataset.read(VariableDataBuffer.data(), H5::PredType::NATIVE_DOUBLE);
                auto itNode = std::begin(VariableDataBuffer);
                std::for_each(std::begin(mrModelPart.Nodes()), std::end(mrModelPart.Nodes()),
                   [&](NodeType& rNode)
                   {
                        rNode.FastGetSolutionStepValue(KratosComponents< Variable<double> >::Get(mVariableNames[i])) = *itNode;
                        itNode++;
                   });
            }
            else if (KratosComponents< Variable<array_1d<double,3> > >::Has(mVariableNames[i]))
            {
                std::vector<double> VariableDataBuffer(3 * mrModelPart.Nodes().size());
                VariableDataset.read(VariableDataBuffer.data(), H5::PredType::NATIVE_DOUBLE);
                unsigned int BlockBegin = 0;
                for (auto it = std::begin(mrModelPart.Nodes()); it != std::end(mrModelPart.Nodes()); it++)
                {
                    array_1d<double,3>& rData = it->FastGetSolutionStepValue(KratosComponents< Variable<array_1d<double,3> > >::Get(mVariableNames[i]));
                    for (unsigned int k = 0; k < 3; k++)
                        rData[k] = VariableDataBuffer[BlockBegin + k];
                    BlockBegin += 3;
                }
            }
        }

        KRATOS_CATCH("")
    }

    void ExecuteFinalizeSolutionStep() override
    {
    }

    void ExecuteBeforeOutputStep() override
    {
    }

    void ExecuteAfterOutputStep() override
    {
    }

    void ExecuteFinalize() override
    {
    }

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
    std::string Info() const override
    {
        return "InputPrimalSolutionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << this->Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    std::string mFilename;
    SizeType mNumNodes;
    std::vector<std::string> mVariableNames;
    std::vector<std::string> mNodalValueNames;
    std::vector<std::string> mProcessInfoValueNames;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
}; /* Class InputPrimalSolutionProcess */

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  InputPrimalSolutionProcess& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InputPrimalSolutionProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Adjoint Fluid Application group

}  /* namespace Kratos */

#endif /* KRATOS_INPUT_PRIMAL_SOLUTION_PROCESS_H_INCLUDED defined */
