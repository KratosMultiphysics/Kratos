//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_DRAG_RESPONSE_FUNCTION)
#define KRATOS_DRAG_RESPONSE_FUNCTION

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/process_info.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "custom_utilities/response_function.h"

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// A response function for drag.
template <unsigned int TDim>
class DragResponseFunction : public ResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(DragResponseFunction);

    typedef ResponseFunction BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DragResponseFunction(ModelPart& rModelPart, Parameters& rParameters)
      : ResponseFunction(rModelPart, rParameters)
    {
        KRATOS_TRY;

        Parameters DefaultParams(R"(
        {
            "response_type": "drag",
            "structure_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "sensitivity_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "nodal_sensitivity_variables" : [],
            "drag_direction": [1.0, 0.0, 0.0]
        })");

        rParameters.ValidateAndAssignDefaults(DefaultParams);

        mStructureModelPartName = rParameters["structure_model_part_name"].GetString();

        if (rParameters["drag_direction"].IsArray() == false ||
            rParameters["drag_direction"].size() != 3)
        {
            KRATOS_THROW_ERROR(std::runtime_error,
                               "drag_direction vector is not a vector or does "
                               "not have size 3:",
                               rParameters.PrettyPrintJsonString())
        }

        for (unsigned int d = 0; d < TDim; ++d)
            mDragDirection[d] = rParameters["drag_direction"][d].GetDouble();

        if (std::abs(norm_2(mDragDirection) - 1.0) > 1e-3)
        {
            const double magnitude = norm_2(mDragDirection);
            if (magnitude == 0.0)
                KRATOS_THROW_ERROR(std::runtime_error,
                                   "drag_direction is not properly defined.",
                                   "")

            std::cout
                << "WARNING: non unit vector detected in \"drag_direction\": "
                << rParameters.PrettyPrintJsonString() << std::endl;
            std::cout << "normalizing \"drag_direction\"..." << std::endl;

            for (unsigned int d = 0; d < TDim; d++)
                mDragDirection[d] /= magnitude;
        }

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~DragResponseFunction() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        KRATOS_TRY;

        BaseType::Initialize();

        ModelPart& r_model_part = this->GetModelPart();

        if (r_model_part.HasSubModelPart(mStructureModelPartName) == false)
            KRATOS_ERROR << "Invalid structure_model_part_name: \""
                         << mStructureModelPartName << "\"." << std::endl;

#pragma omp parallel
        {
            ModelPart::NodeIterator nodes_begin;
            ModelPart::NodeIterator nodes_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Nodes(), nodes_begin, nodes_end);

            for (auto it = nodes_begin; it != nodes_end; ++it)
                it->Set(STRUCTURE, false);
        }

        // mark structure
        ModelPart& r_structure_model_part = r_model_part.GetSubModelPart(mStructureModelPartName);

#pragma omp parallel
        {
            ModelPart::NodeIterator nodes_begin;
            ModelPart::NodeIterator nodes_end;
            OpenMPUtils::PartitionedIterators(r_structure_model_part.Nodes(), nodes_begin, nodes_end);
            for (auto it = nodes_begin; it != nodes_end; ++it)
                it->Set(STRUCTURE, true);
        }

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        ModelPart& r_model_part = this->GetModelPart();

        // allocate auxiliary memory. this is done here instead of Initialize()
        // in case of restart.
        int num_threads = OpenMPUtils::GetNumThreads();
        mElementIds.resize(num_threads);
        mDragFlagVector.resize(num_threads);

        // use first element to initialize drag flag vector
        Element& r_elem = *std::begin(r_model_part.Elements());
#pragma omp parallel
        {
            // initialize drag flag and element id vectors
            int k = OpenMPUtils::ThisThread();
            mElementIds[k] = r_elem.Id() + 1; // force initialization
            this->GetDragFlagVector(r_elem);
        }

        KRATOS_CATCH("");
    }

    void CalculateFirstDerivativesGradient(const Element& rElem,
                                           const Matrix& rAdjointMatrix,
                                           Vector& rGradient,
                                           ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;

        if (rGradient.size() != rAdjointMatrix.size1())
            rGradient.resize(rAdjointMatrix.size1(), false);

        Vector& r_drag_flag_vector = this->GetDragFlagVector(rElem);
        noalias(rGradient) = prod(rAdjointMatrix, r_drag_flag_vector);

        KRATOS_CATCH("");
    }

    void CalculateSecondDerivativesGradient(const Element& rElem,
                                            const Matrix& rAdjointMatrix,
                                            Vector& rGradient,
                                            ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;

        if (rGradient.size() != rAdjointMatrix.size1())
            rGradient.resize(rAdjointMatrix.size1(), false);

        Vector& r_drag_flag_vector = this->GetDragFlagVector(rElem);
        noalias(rGradient) = prod(rAdjointMatrix, r_drag_flag_vector);

        KRATOS_CATCH("");
    }

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void CalculateSensitivityGradient(const Element& rElem,
                                      const Variable<array_1d<double,3>>& rVariable,
                                      const Matrix& rDerivativesMatrix,
                                      Vector& rGradient,
                                      ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;

        if (rGradient.size() != rDerivativesMatrix.size1())
            rGradient.resize(rDerivativesMatrix.size1(), false);

        Vector& r_drag_flag_vector = this->GetDragFlagVector(rElem);
        noalias(rGradient) = prod(rDerivativesMatrix, r_drag_flag_vector);

        KRATOS_CATCH("");
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::string mStructureModelPartName;
    array_1d<double, TDim> mDragDirection;
    std::vector<Vector> mDragFlagVector;
    std::vector<unsigned int> mElementIds;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    Vector& GetDragFlagVector(const Element& rElement)
    {
        int k = OpenMPUtils::ThisThread();

        // if needed, compute the drag flag vector for this element
        if (rElement.Id() != mElementIds[k])
        {
            const unsigned int num_nodes = rElement.GetGeometry().PointsNumber();
            const unsigned int local_size = (TDim + 1) * num_nodes;

            if (mDragFlagVector[k].size() != local_size)
                mDragFlagVector[k].resize(local_size, false);

            unsigned int local_index = 0;
            for (unsigned int i_node = 0; i_node < num_nodes; ++i_node)
            {
                if (rElement.GetGeometry()[i_node].Is(STRUCTURE))
                {
                    for (unsigned int d = 0; d < TDim; d++)
                        mDragFlagVector[k][local_index++] = mDragDirection[d];
                }
                else
                {
                    for (unsigned int d = 0; d < TDim; d++)
                        mDragFlagVector[k][local_index++] = 0.0;
                }

                mDragFlagVector[k][local_index++] = 0.0; // pressure dof
            }

            mElementIds[k] = rElement.Id();
        }

        return mDragFlagVector[k];
    }

    ///@}
};

///@} // Kratos Classes

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_DRAG_RESPONSE_FUNCTION defined */
