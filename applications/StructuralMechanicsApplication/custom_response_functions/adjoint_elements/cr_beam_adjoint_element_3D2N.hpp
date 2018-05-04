// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder 
//


#if !defined(KRATOS_CR_BEAM_ADJOINT_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_CR_BEAM_ADJOINT_ELEMENT_3D2N_H_INCLUDED


#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "custom_elements/cr_beam_element_linear_3D2N.hpp"

namespace Kratos
{

    class CrBeamAdjointElement3D2N : public CrBeamElementLinear3D2N
    {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(CrBeamAdjointElement3D2N);

        typedef Element BaseType;
        typedef BaseType::GeometryType GeometryType;
        typedef BaseType::NodesArrayType NodesArrayType;
        typedef BaseType::PropertiesType PropertiesType;
        typedef BaseType::IndexType IndexType;
        typedef BaseType::SizeType SizeType;
        typedef BaseType::MatrixType MatrixType;
        typedef BaseType::VectorType VectorType;
        typedef BaseType::EquationIdVectorType EquationIdVectorType;
        typedef BaseType::DofsVectorType DofsVectorType;

        CrBeamAdjointElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry);
        CrBeamAdjointElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties);

        ~CrBeamAdjointElement3D2N() override;

        BaseType::Pointer Create(
            IndexType NewId,
            NodesArrayType const& rThisNodes,
            PropertiesType::Pointer pProperties) const override;

        BaseType::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
         PropertiesType::Pointer pProperties) const override;

        void EquationIdVector(
            EquationIdVectorType& rResult,
            ProcessInfo& rCurrentProcessInfo) override;

        void GetDofList(
            DofsVectorType& rElementalDofList,
            ProcessInfo& rCurrentProcessInfo) override;

        double GetDisturbanceMeasureCorrectionFactor(const Variable<double>& rVariable);

        double GetDisturbanceMeasureCorrectionFactor(const Variable<array_1d<double,3>>& rDesignVariable);

        void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;
        void Calculate(const Variable<Vector >& rVariable,
                           Vector& rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override;

        void Calculate(const Variable<Matrix >& rVariable,
                           Matrix& rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable, Matrix& rOutput,
                                                   const ProcessInfo& rCurrentProcessInfo);

        void CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable,
                                            const Variable<Vector>& rStressVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo);

        void CalculateStressDesignVariableDerivative(const Variable<array_1d<double,3>>& rDesignVariable,
                                            const Variable<Vector>& rStressVariable,
                                            Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);

        void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                          std::vector<double>& rOutput,
                          const ProcessInfo& rCurrentProcessInfo) override;

        void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                         std::vector<double>& rValues,
                         const ProcessInfo& rCurrentProcessInfo) override;

        void GetValuesVector(Vector& rValues, int Step = 0) override;

        int Check(const ProcessInfo& rCurrentProcessInfo) override;


    protected:
        CrBeamAdjointElement3D2N(): CrBeamElementLinear3D2N()
        {
        }


    private:
        friend class Serializer;
        void save(Serializer& rSerializer) const override;
        void load(Serializer& rSerializer) override;
    };


}

#endif
