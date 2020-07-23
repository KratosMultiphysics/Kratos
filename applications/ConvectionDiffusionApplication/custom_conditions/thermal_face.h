// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Riccardo Rossi
//                   Ruben Zorrilla
//

#if !defined(KRATOS_THERMAL_FACE_H_INCLUDED )
#define  KRATOS_THERMAL_FACE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "includes/condition.h"
#include "geometries/geometry.h"
#include "includes/variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// A basic Neumann condition for convection-diffusion problems.
/** It applies a flux condition based on the nodal values of the
 *  variable defined as SurfaceSourceVariable by the
 *  CONVECTION_DIFFUSION_SETTINGS variable in the given ProcessInfo.
 */
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) ThermalFace: public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ThermalFace
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ThermalFace);

    /**
     * @brief Gauss pt. data structure
     * Auxiliar data structure to pass the Gauss pt. data
     */
    struct ConditionDataStruct
    {
        double Weight;              // Gauss point weight
        array_1d<double, 3> Normal; // Condition normal
        Vector N;                   // Gauss point shape functions values

        double Emissivity;            // Ambient emissivity value
        double AmbientTemperature;    // Ambient temperature value
        double ConvectionCoefficient; // Ambient convection coefficient
        Vector UnknownValues;         // Previous iteration unknown values
        Vector FaceHeatFluxValues;    // Nodal face heat flux values

        double inline GaussPointUnknown() const
        {
            return InterpolateInGaussPoint(UnknownValues);
        }

        double inline GaussPointFaceHeatFlux() const
        {
            return InterpolateInGaussPoint(FaceHeatFluxValues);
        }

        double inline InterpolateInGaussPoint(const Vector &rNodalValues) const
        {
            double gauss_pt_val = 0.0;
            for (unsigned int i = 0; i < N.size(); ++i) {
                gauss_pt_val += N[i] * rNodalValues[i];
            }
            return gauss_pt_val;
        }
    };

    typedef Condition::MatrixType MatrixType;
    typedef Condition::VectorType VectorType;

    /// Stefan Boltzmann constant for radiation in SI units: [W / (m^2 K^4)].
    constexpr static double StefanBoltzmann = 5.67e-8;

    ///@}
    ///@name Life Cycle
    ///@{

    ThermalFace(
        IndexType NewId,
        Geometry< Node<3> >::Pointer pGeometry);

    ThermalFace(
        IndexType NewId,
        Geometry< Node<3> >::Pointer pGeometry,
        Properties::Pointer pProperties);

    /// Destructor.
    ~ThermalFace() override;

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties) const override;

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(
        DofsVectorType& ConditionalDofList,
        ProcessInfo& CurrentProcessInfo) override;

    GeometryData::IntegrationMethod GetIntegrationMethod() override;

    void GetValueOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(
        const Variable<array_1d<double, 6 > >& rVariable,
        std::vector<array_1d<double, 6 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

protected:

    ///@name Protected Life Cycle
    ///@{

    // Internal default constructor for serialization
    ThermalFace();

    ///@}
    ///@name Protected Operations
    ///@{

    void AddIntegrationPointRHSContribution(
        VectorType& rRightHandSideVector,
        const ConditionDataStruct &rData);

    void AddIntegrationPointLHSContribution(
        MatrixType& rLeftHandSideMatrix,
        const ConditionDataStruct &rData);

    void FillConditionDataStructure(
        const ProcessInfo &rCurrentProcessInfo,
        ConditionDataStruct &rData);

    ///@}

private:

    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ThermalFace& operator=(ThermalFace const& rOther);

    /// Copy constructor.
    ThermalFace(ThermalFace const& rOther);

    ///@}

}; // Class ThermalFace

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_THERMAL_FACE_H_INCLUDED  defined
