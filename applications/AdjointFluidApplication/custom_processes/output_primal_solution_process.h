//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_OUTPUT_PRIMAL_SOLUTION_PROCESS_H_INCLUDED )
#define KRATOS_OUTPUT_PRIMAL_SOLUTION_PROCESS_H_INCLUDED

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

/// Process to output primal solution (e.g., velocity and pressure) needed for adjoint solution.
/**
 * Requires mesh topology and node numbering to remain constant.
 */
class OutputPrimalSolutionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(OutputPrimalSolutionProcess);

    typedef ModelPart::NodeType NodeType;

    typedef hsize_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    OutputPrimalSolutionProcess(ModelPart& rModelPart, Parameters& rParameters)
    : Process(), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        Parameters DefaultParams(R"(
        {
            "model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "file_name": "PLEASE_SPECIFY_H5_FILE_NAME",
            "variable_list": ["VELOCITY", "ACCELERATION", "PRESSURE"],
            "nodal_value_list": [],
            "processinfo_value_list": [],
            "alpha_bossak": 0.0
        })");

        // try accessing parameters without defaults so that an error is thrown
        // if they don't exist
        rParameters["file_name"];
        rParameters.ValidateAndAssignDefaults(DefaultParams);
        
        // output filename
        std::stringstream FilenameStream;
        FilenameStream << rParameters["file_name"].GetString() << '_' << mrModelPart.GetCommunicator().MyPID() << ".h5";
        mFilename = FilenameStream.str();

        // nodal variable names to output
        mVariableNames.resize(rParameters["variable_list"].size());
        for (unsigned int i = 0; i < mVariableNames.size(); i++)
        {
            mVariableNames[i] = rParameters["variable_list"].GetArrayItem(i).GetString();
        }

        // nodal value names to output
        mNodalValueNames.resize(rParameters["nodal_value_list"].size());
        for( unsigned int i=0; i < mNodalValueNames.size(); ++i )
        {
            mNodalValueNames[i] = rParameters["nodal_value_list"].GetArrayItem(i).GetString();
        }

        // nodal value names to output
        mProcessInfoValueNames.resize(rParameters["processinfo_value_list"].size());
        for( unsigned int i=0; i < mProcessInfoValueNames.size(); ++i )
        {
            mProcessInfoValueNames[i] = rParameters["processinfo_value_list"].GetArrayItem(i).GetString();
        }

        mAlphaBossak = rParameters["alpha_bossak"].GetDouble();

        mNumNodes = 0;

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~OutputPrimalSolutionProcess() override {}

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

        // file initialization
        H5::H5File File(mFilename.c_str(), H5F_ACC_TRUNC);

        File.createGroup("NodalData");

        // create vector of node ids
        std::vector<unsigned int> NodeIdBuffer(mrModelPart.Nodes().size());
        auto itBuffer = std::begin(NodeIdBuffer);
        std::for_each(std::begin(mrModelPart.Nodes()), std::end(mrModelPart.Nodes()),
           [&](NodeType& rNode)
           {
                *itBuffer = static_cast<unsigned int>(rNode.Id());
                itBuffer++;
           });

        // write node ids to file
        hsize_t Dims[1];
        Dims[0] = mNumNodes;
        H5::DataSpace NodeIdDataspace(1, Dims);

        H5::DataSet NodeIdDataset(File.createDataSet("/NodalData/Id", H5::PredType::NATIVE_UINT, NodeIdDataspace));

        NodeIdDataset.write(NodeIdBuffer.data(), H5::PredType::NATIVE_UINT);

        KRATOS_CATCH("")
    }

    void ExecuteBeforeSolutionLoop() override
    {
    }

    void ExecuteInitializeSolutionStep() override
    {
    }

    void ExecuteFinalizeSolutionStep() override
    {
        KRATOS_TRY

        // open file to append data
        H5::H5File File(mFilename.c_str(), H5F_ACC_RDWR);

        H5::DataSet NodeIdDataset = File.openDataSet("/NodalData/Id");

        hsize_t Dims[2];
        NodeIdDataset.getSpace().getSimpleExtentDims(Dims);

        if (mNumNodes != Dims[0])
            KRATOS_THROW_ERROR(std::runtime_error, "inconsistent dimension in file: ", mFilename)

        if (mNumNodes != static_cast<SizeType>(mrModelPart.Nodes().size()))
            KRATOS_THROW_ERROR(std::runtime_error, "detected change in number of nodes.", "")

        std::vector<unsigned int> NodeIdBuffer(mrModelPart.Nodes().size());
        NodeIdDataset.read(NodeIdBuffer.data(), H5::PredType::NATIVE_UINT);

        // check that the node ordering is consistent between the model part and file
        int Delta = 0;
        auto itNode = std::begin(NodeIdBuffer);
        std::for_each(std::begin(mrModelPart.Nodes()), std::end(mrModelPart.Nodes()),
           [&](NodeType& rNode)
           {
                Delta += std::abs(static_cast<int>(rNode.Id()) - static_cast<int>(*itNode));
                itNode++;
           });

        if (Delta != 0)
            KRATOS_THROW_ERROR(std::runtime_error, "detected mismatch of node ids in file: ", mFilename);

        // output time step data
        std::stringstream TimeStepPathStream;
        TimeStepPathStream << "/NodalData/" << std::fixed << std::setprecision(10) << mrModelPart.GetProcessInfo()[TIME];
        std::string TimeStepGroupName = TimeStepPathStream.str();
        File.createGroup(TimeStepGroupName.c_str());

        // write each variable
        for (unsigned int i = 0; i < mVariableNames.size(); i++)
        {
            std::string DatasetName = TimeStepGroupName + "/" + mVariableNames[i];
            if (KratosComponents< Variable<double> >::Has(mVariableNames[i]))
            {
                // double variable
                H5::DataSpace VariableDataspace(1, Dims);
                H5::DataSet VariableDataset(File.createDataSet(DatasetName.c_str(), H5::PredType::NATIVE_DOUBLE, VariableDataspace));
                std::vector<double> VariableDataBuffer(mrModelPart.Nodes().size());
                auto itNode = std::begin(VariableDataBuffer);

                // fill buffer
                std::for_each(std::begin(mrModelPart.Nodes()), std::end(mrModelPart.Nodes()),
                   [&](NodeType& rNode)
                   {
                        *itNode = rNode.FastGetSolutionStepValue(KratosComponents< Variable<double> >::Get(mVariableNames[i]));
                        itNode++;
                   });

                // write buffer
                VariableDataset.write(VariableDataBuffer.data(), H5::PredType::NATIVE_DOUBLE);
            }
            else if (KratosComponents< Variable<array_1d<double,3> > >::Has(mVariableNames[i]))
            {
                // array_1d variable
                Dims[1] = 3;
                H5::DataSpace VariableDataspace(2, Dims);
                H5::DataSet VariableDataset(File.createDataSet(DatasetName.c_str(), H5::PredType::NATIVE_DOUBLE, VariableDataspace));
                std::vector<double> VariableDataBuffer(3 * mrModelPart.Nodes().size());
                const Variable< array_1d<double, 3> >& rVariable = KratosComponents< Variable<array_1d<double,3> > >::Get(mVariableNames[i]);

                // fill the buffer
                if (rVariable == ACCELERATION)
                {
                    // Bossak acceleration
                    unsigned int BlockBegin = 0;
                    for (auto it = std::begin(mrModelPart.Nodes()); it != std::end(mrModelPart.Nodes()); it++)
                    {
                        array_1d<double, 3> accel = (1.0 - mAlphaBossak) * it->FastGetSolutionStepValue(rVariable,0);
                        accel += mAlphaBossak * it->FastGetSolutionStepValue(rVariable,1);
                        for (unsigned int k = 0; k < 3; k++)
                            VariableDataBuffer[BlockBegin + k] = accel[k];
                        BlockBegin += 3;
                    }
                }
                else
                {
                    unsigned int BlockBegin = 0;
                    for (auto it = std::begin(mrModelPart.Nodes()); it != std::end(mrModelPart.Nodes()); it++)
                    {
                        const array_1d<double, 3>& rData = it->FastGetSolutionStepValue(rVariable);
                        for (unsigned int k = 0; k < 3; k++)
                            VariableDataBuffer[BlockBegin + k] = rData[k];
                        BlockBegin += 3;
                    }
                }

                // write buffer
                VariableDataset.write(VariableDataBuffer.data(), H5::PredType::NATIVE_DOUBLE);
            }
        }

        KRATOS_CATCH("")
    }

    void ExecuteBeforeOutputStep() override
    {
    }

    void ExecuteAfterOutputStep() override
    {
    }

    void ExecuteFinalize() override
    {
        // write the values to the file after finishing the computation
        KRATOS_TRY

        if( mNodalValueNames.size() > 0 || mProcessInfoValueNames.size() > 0 )
        {
            // open file to append data
            H5::H5File file(mFilename.c_str(), H5F_ACC_RDWR);
            
            H5::DataSet NodeIdDataset = file.openDataSet("/NodalData/Id");

            hsize_t dims[3];
            NodeIdDataset.getSpace().getSimpleExtentDims(dims);
            
            if (mNumNodes != dims[0])
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
            
            // output of the nodal values
            // std::stringstream nodal_value_stream;
            // nodal_value_stream << "/NodalData/Values";
            // std::string values_group_name = nodal_value_stream.str();
            std::string values_group_name = "/NodalData/Values";
            file.createGroup(values_group_name.c_str());

            //write each value
            for( std::size_t i = 0; i < mNodalValueNames.size(); ++i )
            {
                std::string dataset_name = values_group_name + "/" + mNodalValueNames[i];
                if( KratosComponents< Variable<Vector> >::Has(mNodalValueNames[i]) )
                {
                    const Variable< Vector >& rVariable = KratosComponents< Variable<Vector> >::Get(mNodalValueNames[i]);
                    
                    //this assumes, that all vectors on the nodes have the same length:
                    unsigned int vector_size = 0;
                    if( mrModelPart.Nodes().size() > 0 )
                    {
                        vector_size = std::begin(mrModelPart.Nodes())->GetValue(rVariable).size();
                    }
                    
                    dims[1] = vector_size;
                    dims[2] = 0;

                    H5::DataSpace value_dataspace(2, dims);
                    H5::DataSet value_dataset(file.createDataSet(dataset_name.c_str(), H5::PredType::NATIVE_DOUBLE, value_dataspace));
                    std::vector<double> value_data_buffer(vector_size * mrModelPart.Nodes().size());
                    
                    unsigned int block_begin = 0;
                    for( auto it = std::begin(mrModelPart.Nodes()); it != std::end(mrModelPart.Nodes()); ++it )
                    {
                        auto& rData = it->GetValue(rVariable);
                        if( rData.size() != vector_size )
                            KRATOS_ERROR << "value vector " << mNodalValueNames[i] << " on node " << it->Id() << " has an unexpected length" << std::endl;

                        for( std::size_t j = 0; j < vector_size; ++j)
                            value_data_buffer[block_begin + j] = rData[j];
                        block_begin += vector_size;
                    }

                    //write buffer
                    value_dataset.write(value_data_buffer.data(), H5::PredType::NATIVE_DOUBLE);
               
                }

                else if( KratosComponents< Variable<Matrix> >::Has(mNodalValueNames[i]) )
                {
                    const Variable< Matrix >& rVariable = KratosComponents< Variable<Matrix> >::Get(mNodalValueNames[i]);
                    
                    //this assumes, that all matrices on the nodes have the same size:
                    const unsigned int size_1 = std::begin(mrModelPart.Nodes())->GetValue(rVariable).size1();
                    const unsigned int size_2 = std::begin(mrModelPart.Nodes())->GetValue(rVariable).size2();
                    dims[1] = size_2; //dofs
                    dims[2] = size_1; //modes

                    H5::DataSpace value_dataspace(3, dims);
                    H5::DataSet value_dataset(file.createDataSet(dataset_name.c_str(), H5::PredType::NATIVE_DOUBLE, value_dataspace));
                    std::vector<double> value_data_buffer(size_1 * size_2 * mrModelPart.Nodes().size());

                    unsigned int block_begin = 0;
                    for( auto it = std::begin(mrModelPart.Nodes()); it != std::end(mrModelPart.Nodes()); ++it )
                    {
                        auto& rData = it->GetValue(rVariable);
                        if( rData.size1() != size_1 || rData.size2() != size_2 )
                            KRATOS_ERROR << "value matrix " << mNodalValueNames[i] << " on node " << it->Id() << " has an unexpected length" << std::endl;

                        for( std::size_t j = 0; j < size_2; ++j )
                        {
                            for( std::size_t k = 0; k < size_1; ++k)
                            {
                                value_data_buffer[block_begin + k] = rData(k,j);
                            }
                            block_begin += size_1;
                        }
                    }

                    value_dataset.write(value_data_buffer.data(), H5::PredType::NATIVE_DOUBLE);
                }

                else
                {
                    std::cout << "output_primal_solution_process.h: output of this value type has not yet been implemented!" << std::endl;
                }
            } //nodal values

            std::string procinfo_group_name = "/ProcessInfo";
            file.createGroup(procinfo_group_name.c_str());

            //write each value
            for( std::size_t i = 0; i < mProcessInfoValueNames.size(); ++i )
            {
                hsize_t procinfo_dims[1];
                std::string dataset_name = procinfo_group_name + "/" + mProcessInfoValueNames[i];
                if( KratosComponents< Variable<Vector> >::Has(mProcessInfoValueNames[i]) )
                {
                    const Variable< Vector >& rVariable = KratosComponents< Variable<Vector> >::Get(mProcessInfoValueNames[i]);

                    const auto& rProcessInfo = mrModelPart.GetProcessInfo();
                    const auto& rData = rProcessInfo[rVariable];
                    const unsigned int vector_size = rData.size();
                    
                    procinfo_dims[0] = vector_size;

                    H5::DataSpace procinfo_dataspace(1, procinfo_dims);
                    H5::DataSet procinfo_dataset(file.createDataSet(dataset_name.c_str(), H5::PredType::NATIVE_DOUBLE, procinfo_dataspace));
                    std::vector<double> value_data_buffer(vector_size);

                    for( size_t j = 0; j < vector_size; ++j )
                    {
                        value_data_buffer[j] = rData[j];
                    }

                    //write buffer
                    procinfo_dataset.write(value_data_buffer.data(), H5::PredType::NATIVE_DOUBLE);
               
                }
            }

        }

        KRATOS_CATCH("")
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
        return "OutputPrimalSolutionProcess";
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
    double mAlphaBossak;

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
}; /* Class OutputPrimalSolutionProcess */

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  OutputPrimalSolutionProcess& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const OutputPrimalSolutionProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Adjoint Fluid Application group

}  /* namespace Kratos */

#endif /* KRATOS_OUTPUT_PRIMAL_SOLUTION_PROCESS_H_INCLUDED defined */
