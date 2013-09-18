//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_BEAM_ELEMENT_H_INCLUDED )
#define  KRATOS_BEAM_ELEMENT_H_INCLUDED


// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"



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

/// Beam Element for 3D space dimension

/**
 * Implements a Small Displacement definition for structural analysis.
 * This works for line geometries in 3D :: it must be extended to 2D and large displacement
 */

class BeamElement
    :public Element
{

    typedef GeometryData::IntegrationMethod IntegrationMethod;


private:
    ///@name Static Member Variables



    double mArea;                            // Area de la seccion tranversal de la viga.
    double mInertia_x;                       // Momento de Inercia alredor del eje Ix local.
    double mInertia_y;                       // Momento de Inercia alrededor del eje Iy local.
    double mInertia_Polar;                   // Momento Polar de Inercia
    double mlength;                          // Longitud del Elemento.


    void CalculateSectionProperties();

    void CalculateLocalMatrix(Matrix& LocalMatrix);

    void CalculateTransformationMatrix(Matrix& Rotation);

    void CalculateBodyForce(Matrix& Rotation,  Vector& LocalBody, Vector& GlobalBody);

    void CalculateLocalNodalStress(Vector& Stress);

    double CalculateInternalAxil(   const double& Ao, const double& Load, const double& X);
    double CalculateInternalShear(  const double& Vo, const double& Load, const double& X);
    double CalculateInternalMoment( const double& Mo, const double& Vo,   const double& Load, const double& X);

    void CalculateDistributedBodyForce(const int Direction, Vector& Load);


public:



    KRATOS_CLASS_POINTER_DEFINITION(BeamElement);

    /// Default constructor.
    BeamElement(IndexType NewId, GeometryType::Pointer pGeometry);
    BeamElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~BeamElement();


    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

    void Initialize();

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);

    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

    void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);

    void GetValuesVector(Vector& values, int Step);

    void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);

    void CalculateRHS(Vector& rRightHandSideVector);

    void CalculateLHS(Matrix& rLeftHandSideMatrix);

    void CalculateElementalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                                  ProcessInfo& rCurrentProcessInfo,
                                  bool CalculateStiffnessMatrixFlag,
                                  bool CalculateResidualVectorFlag);

    void  GetFirstDerivativesVector(Vector& values, int Step);
    void  GetSecondDerivativesVector(Vector& values, int Step);

    void CalculateOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable,
                                       std::vector< array_1d<double,3> >& Output,
                                       const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable,
                                      std::vector<array_1d<double,3> >& rValues,
                                      const ProcessInfo& rCurrentProcessInfo);

    IntegrationMethod GetIntegrationMethod() const;

    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo);

private:
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;


    // A private default constructor necessary for serialization
    BeamElement() {};


    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
    }




}; // Class BeamElement

} // Namespace Kratos.


#endif //  KRATOS_BEAM_ELEMENT_H_INCLUDED defined

