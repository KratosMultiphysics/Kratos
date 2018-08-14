// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    
//

#if !defined(KRATOS_ADJOINT_SOLID_ELEMENT_H_INCLUDED )
#define  KRATOS_ADJOINT_SOLID_ELEMENT_H_INCLUDED

// System includes

// External include

// Project includes
#include "includes/define.h"
#include "includes/element.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class AdjointSolidElement
 * @ingroup StructuralMechanicsApplication
 * @brief A template class for creating adjoint elements for solids.
 */
template <class TPrimalElement>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) AdjointSolidElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AdjointSolidElement);

    ///@}
    ///@name Life Cycle
    ///@{

    AdjointSolidElement(IndexType NewId = 0)
        : Element(NewId), mPrimalElement(NewId, pGetGeometry())
    {
    }

    AdjointSolidElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry), mPrimalElement(NewId, pGeometry)
    {
    }

    AdjointSolidElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties), mPrimalElement(NewId, pGeometry, pProperties)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared<AdjointSolidElement<TPrimalElement>>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    void Initialize() override
    {
        KRATOS_TRY;
        mPrimalElement.Initialize();
        KRATOS_CATCH("");
    }

    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        mPrimalElement.InitializeSolutionStep(rCurrentProcessInfo);
        KRATOS_CATCH("");
    }

    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        mPrimalElement.InitializeNonLinearIteration(rCurrentProcessInfo);
        KRATOS_CATCH("");
    }

    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        mPrimalElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
        KRATOS_CATCH("");
    }
    
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        mPrimalElement.FinalizeSolutionStep(rCurrentProcessInfo);
        KRATOS_CATCH("");
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        mPrimalElement.CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
        noalias(rLeftHandSideMatrix) = -rLeftHandSideMatrix;
        KRATOS_CATCH("");
    }

    void GetValuesVector(Vector& rValues, int Step = 0) override
    {
        KRATOS_TRY;
        const auto& r_geom = mPrimalElement.GetGeometry();
        const unsigned dimension = r_geom.WorkingSpaceDimension();
        const unsigned mat_size = r_geom.size() * dimension;
        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);
        for (unsigned int i = 0; i < r_geom.size(); ++i)
        {
            const array_1d<double, 3>& adjoint_displacement =
                r_geom[i].FastGetSolutionStepValue(ADJOINT_DISPLACEMENT, Step);
            const unsigned index = i * dimension;
            for (unsigned k = 0; k < dimension; ++k)
                rValues[index + k] = adjoint_displacement[k];
        }
        KRATOS_CATCH("");
    }

    void EquationIdVector(EquationIdVectorType& rResult,
                          ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        auto& r_geom = mPrimalElement.GetGeometry();
        const unsigned number_of_nodes = r_geom.size();
        const unsigned dimension = r_geom.WorkingSpaceDimension();

        if (rResult.size() != dimension * number_of_nodes)
            rResult.resize(dimension * number_of_nodes, false);

        const unsigned pos = r_geom[0].GetDofPosition(ADJOINT_DISPLACEMENT_X);

        if (dimension == 2)
        {
            for (unsigned i = 0; i < number_of_nodes; ++i)
            {
                const unsigned index = i * 2;
                rResult[index] =
                    r_geom[i].GetDof(ADJOINT_DISPLACEMENT_X, pos).EquationId();
                rResult[index + 1] =
                    r_geom[i].GetDof(ADJOINT_DISPLACEMENT_Y, pos + 1).EquationId();
            }
        }
        else
        {
            for (unsigned i = 0; i < number_of_nodes; ++i)
            {
                const unsigned index = i * 3;
                rResult[index] =
                    r_geom[i].GetDof(ADJOINT_DISPLACEMENT_X, pos).EquationId();
                rResult[index + 1] =
                    r_geom[i].GetDof(ADJOINT_DISPLACEMENT_Y, pos + 1).EquationId();
                rResult[index + 2] =
                    r_geom[i].GetDof(ADJOINT_DISPLACEMENT_Z, pos + 2).EquationId();
            }
        }
        KRATOS_CATCH("");
    }

    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        auto& r_geom = mPrimalElement.GetGeometry();
        const unsigned number_of_nodes = r_geom.size();
        const unsigned dimension = r_geom.WorkingSpaceDimension();
        rElementalDofList.resize(0);
        rElementalDofList.reserve(dimension * number_of_nodes);

        if (dimension == 2)
        {
            for (unsigned i = 0; i < number_of_nodes; ++i)
            {
                rElementalDofList.push_back(r_geom[i].pGetDof(ADJOINT_DISPLACEMENT_X));
                rElementalDofList.push_back(r_geom[i].pGetDof(ADJOINT_DISPLACEMENT_Y));
            }
        }
        else
        {
            for (unsigned int i = 0; i < number_of_nodes; ++i)
            {
                rElementalDofList.push_back(r_geom[i].pGetDof(ADJOINT_DISPLACEMENT_X));
                rElementalDofList.push_back(r_geom[i].pGetDof(ADJOINT_DISPLACEMENT_Y));
                rElementalDofList.push_back(r_geom[i].pGetDof(ADJOINT_DISPLACEMENT_Z));
            }
        }
        KRATOS_CATCH("");
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        return mPrimalElement.Check(rCurrentProcessInfo);
        KRATOS_CATCH("");
    }

    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rDesignVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        mPrimalElement.CalculateSensitivityMatrix(rDesignVariable, rOutput, rCurrentProcessInfo);
        KRATOS_CATCH("");
    }

    ///@}

protected:

private:
    ///@name Member Variables
    ///@{

    TPrimalElement mPrimalElement;

    ///@}
    ///@name Private Operations
    ///@{
   
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
        rSerializer.save("mPrimalElement", mPrimalElement);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
        rSerializer.load("mPrimalElement", mPrimalElement);
    }
    ///@}
};

///@}

} // namespace Kratos.
#endif // KRATOS_ADJOINT_SOLID_ELEMENT_H_INCLUDED  defined 
