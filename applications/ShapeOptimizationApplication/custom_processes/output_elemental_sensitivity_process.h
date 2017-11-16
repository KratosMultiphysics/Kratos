//  KratosShapeOptimizationApplication
//
//  License:		 BSD License
//					 license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Long Chen
//                   

// fusseder TODO: copied from Adjoint fluid application. Maybe rename the class

#if !defined(KRATOS_OUTPUT_ELEMENTAL_SENSITIVITY_PROCESS_H_INCLUDED )
#define KRATOS_OUTPUT_ELEMENTAL_SENSITIVITY_PROCESS_H_INCLUDED

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

///@addtogroup ShapeOptimizationApplication
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

/// Process to output primal solution (e.g., displacements) needed for adjoint solution.
/**
 * Requires mesh topology and node numbering to remain constant.
 */
class OutputElementalSensitivityProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(OutputElementalSensitivityProcess);

    typedef ModelPart::NodeType NodeType;
    typedef ModelPart::ElementType ElementType;

    typedef hsize_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    OutputElementalSensitivityProcess(ModelPart& rModelPart, Parameters& rParameters)
    : Process(), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        Parameters DefaultParams(R"(
        {
            "model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "file_name": "PLEASE_SPECIFY_H5_FILE_NAME",
            "variable_list": ["ADJOINT_DISPLACEMENT"]
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
        mVariables.resize(rParameters["variable_list"].size());
        for (unsigned int i = 0; i < mVariables.size(); i++)
        {
            mVariables[i] = rParameters["variable_list"].GetArrayItem(i).GetString();
        }

        mAlphaBossak = 0.0; //rParameters["alpha_bossak"].GetDouble();

        mNumElements = 0;

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~OutputElementalSensitivityProcess() override {}

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


            // test not available since there is no global ElementalSolutionStepVariablesList
                /*

            // test for valid input parameters / model part
            for (unsigned int i = 0; i < mVariables.size(); i++)
            {
                if (KratosComponents< Variable<double> >::Has(mVariables[i]))
                {
                    if (mrModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents< Variable<double> >::Get(mVariables[i])) == false)
                    {
                        KRATOS_THROW_ERROR(std::runtime_error, "variable is not found in nodal solution steps variable list: ", mVariables[i])
                    }
                }
                else if (KratosComponents< Variable<array_1d<double,3> > >::Has(mVariables[i]))
                {
                    if (mrModelPart.GetNodalSolutionStepVariablesList().Has(KratosComponents< Variable<array_1d<double,3> > >::Get(mVariables[i])) == false)
                    {
                        KRATOS_THROW_ERROR(std::runtime_error,"variable is not found in nodal solution steps variable list: ", mVariables[i])
                    }
                }
                else
                {
                    KRATOS_THROW_ERROR(std::invalid_argument, "variable type not supported: ", mVariables[i])
                }
            }

            */


        // elements should be ordered by node id before IO begins
        if (mrModelPart.Elements().IsSorted() == false)
        {
            mrModelPart.Elements().Sort();
        }

        // number of Elements must be constant for IO
        mNumElements = static_cast<SizeType>(mrModelPart.Elements().size());

        // file initialization
        H5::H5File File(mFilename.c_str(), H5F_ACC_TRUNC);

        File.createGroup("ElementalData");

        //// create vector of node ids
        //std::vector<unsigned int> NodeIdBuffer(mrModelPart.Nodes().size());
        //auto itBuffer = std::begin(NodeIdBuffer);
        //std::for_each(std::begin(mrModelPart.Nodes()), std::end(mrModelPart.Nodes()),
        //   [&](NodeType& rNode)
        //   {
        //        *itBuffer = static_cast<unsigned int>(rNode.Id());
        //        itBuffer++;
        //   });


        // create vector of element ids
        std::vector<unsigned int> ElementIdBuffer(mrModelPart.Elements().size());
        auto itBuffer = std::begin(ElementIdBuffer);
        std::for_each(std::begin(mrModelPart.Elements()), std::end(mrModelPart.Elements()),
            [&](ElementType& rElement)
        {
            *itBuffer = static_cast<unsigned int>(rElement.Id());
            itBuffer++;
        });

        //// write node ids to file
        //hsize_t Dims[1];
        //Dims[0] = mNumNodes;
        //H5::DataSpace NodeIdDataspace(1, Dims);

        //H5::DataSet NodeIdDataset(File.createDataSet("/NodalData/Id", H5::PredType::NATIVE_UINT, NodeIdDataspace));

        //NodeIdDataset.write(NodeIdBuffer.data(), H5::PredType::NATIVE_UINT);

        // write element ids to file
        hsize_t Dims[1];
        Dims[0] = mNumElements;
        H5::DataSpace ElementIdDataspace(1, Dims);

        H5::DataSet ElementIdDataset(File.createDataSet("/ElementalData/Id", H5::PredType::NATIVE_UINT, ElementIdDataspace));

        ElementIdDataset.write(ElementIdBuffer.data(), H5::PredType::NATIVE_UINT);




        // std::cout << "###############" << std::endl;
        // std::cout << "LCHEN: OutputElementalSensitivityProcess.ExecuteInitialize() finished!" << std::endl;

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
        OutputNodalSensitivity();
        OutputElementalSensitivity();
    }

    void OutputElementalSensitivity()
    {
        KRATOS_TRY;

        // open file to append data
        H5::H5File File(mFilename.c_str(), H5F_ACC_RDWR);

        H5::DataSet ElementIdDataset = File.openDataSet("ElementalData/Id");

        hsize_t Dims[2];
        ElementIdDataset.getSpace().getSimpleExtentDims(Dims);


        if (mNumElements != Dims[0])
        {
            KRATOS_THROW_ERROR(std::runtime_error, "OutputElementalSensitivityProcess: inconsistent dimension in file: ", mFilename);
        }

        if (mNumElements != static_cast<SizeType>(mrModelPart.Elements().size()))
        {
            KRATOS_THROW_ERROR(std::runtime_error, "OutputElementalSensitivityProcess: detected change in number of nodes.", "")
        }


        std::vector<unsigned int> ElementIdBuffer(mrModelPart.Elements().size());
        ElementIdDataset.read(ElementIdBuffer.data(), H5::PredType::NATIVE_UINT);

        // check that the element ordering is consistent between the model part and file
        int Delta = 0;
        auto itElement = std::begin(ElementIdBuffer);
        std::for_each(std::begin(mrModelPart.Elements()), std::end(mrModelPart.Elements()),
            [&](ElementType& rElement)
        {
            Delta += std::abs(static_cast<int>(rElement.Id()) - static_cast<int>(*itElement));
            itElement++;
        });

        if (Delta != 0)
            KRATOS_THROW_ERROR(std::runtime_error, "detected mismatch of node ids in HDF5 file: ", mFilename);

        // output time step data
        std::stringstream TimeStepPatchStream;
        TimeStepPatchStream << "/ElementalData/" << std::fixed << std::setprecision(10) << mrModelPart.GetProcessInfo()[TIME];
        std::string TimeStepGroupName = TimeStepPatchStream.str();
        File.createGroup(TimeStepGroupName.c_str());

        // write each variable
        for (unsigned int i = 0; i < mVariables.size(); i++)
        {
            std::string DatasetName = TimeStepGroupName + "/" + mVariables[i];
            if (KratosComponents< Variable<double> >::Has(mVariables[i]))
            {
                // double variable
                H5::DataSpace VariableDataspace(1, Dims);
                H5::DataSet VariableDataset(File.createDataSet(DatasetName.c_str(), H5::PredType::NATIVE_DOUBLE, VariableDataspace));
                std::vector<double> VariableDataBuffer(mrModelPart.Elements().size());
                auto itElement = std::begin(VariableDataBuffer);

                // fill buffer
                std::for_each(std::begin(mrModelPart.Elements()), std::end(mrModelPart.Elements()),
                    [&](ElementType& rElement)
                {
                    *itElement = rElement.GetValue(KratosComponents< Variable<double> >::Get(mVariables[i]));
                    itElement++;
                });

                // write buffer
                VariableDataset.write(VariableDataBuffer.data(), H5::PredType::NATIVE_DOUBLE);
                
                std::cout << "LCHEN: write elemental data done!" << std::endl;
                std::cout << "LCHEN: write elemental data done!" << std::endl;

            }
            else if (KratosComponents< Variable<array_1d<double, 3> > >::Has(mVariables[i]))
            {
                KRATOS_THROW_ERROR(std::runtime_error, "this type of variable is not yet implemented!", mFilename);
            }

        }



        KRATOS_CATCH("");
    }

    void OutputNodalSensitivity()
    {

        KRATOS_TRY
            /*
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
        for (unsigned int i = 0; i < mVariables.size(); i++)
        {
            std::string DatasetName = TimeStepGroupName + "/" + mVariables[i];
            if (KratosComponents< Variable<double> >::Has(mVariables[i]))
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
                    *itNode = rNode.FastGetSolutionStepValue(KratosComponents< Variable<double> >::Get(mVariables[i]));
                    itNode++;
                });

                // write buffer
                VariableDataset.write(VariableDataBuffer.data(), H5::PredType::NATIVE_DOUBLE);
            }
            else if (KratosComponents< Variable<array_1d<double, 3> > >::Has(mVariables[i]))
            {
                // array_1d variable
                Dims[1] = 3;
                H5::DataSpace VariableDataspace(2, Dims);
                H5::DataSet VariableDataset(File.createDataSet(DatasetName.c_str(), H5::PredType::NATIVE_DOUBLE, VariableDataspace));
                std::vector<double> VariableDataBuffer(3 * mrModelPart.Nodes().size());
                const Variable< array_1d<double, 3> >& rVariable = KratosComponents< Variable<array_1d<double, 3> > >::Get(mVariables[i]);

                // fill the buffer
                if (rVariable == ACCELERATION)
                {
                    // Bossak acceleration
                    unsigned int BlockBegin = 0;
                    for (auto it = std::begin(mrModelPart.Nodes()); it != std::end(mrModelPart.Nodes()); it++)
                    {
                        array_1d<double, 3> accel = (1.0 - mAlphaBossak) * it->FastGetSolutionStepValue(rVariable, 0);
                        accel += mAlphaBossak * it->FastGetSolutionStepValue(rVariable, 1);
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
        }*/

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
        return "OutputElementalSensitivityProcess";
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
    std::vector<std::string> mVariables;
    double mAlphaBossak;

    SizeType mNumElements;

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
}; /* Class OutputElementalSensitivityProcess */

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
    OutputElementalSensitivityProcess& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const OutputElementalSensitivityProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Shape Optimization Application group

}  /* namespace Kratos */

#endif /* KRATOS_OUTPUT_ELEMENTAL_SENSITIVITY_PROCESS_H_INCLUDED defined */
