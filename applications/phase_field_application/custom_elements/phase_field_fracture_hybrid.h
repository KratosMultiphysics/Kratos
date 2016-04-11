/*
==============================================================================
KratosPhaseFieldApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
hgbk2008@gmail.com
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 7 Dec 2015 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_PHASE_FIELD_APPLICATION_PHASE_FIELD_FRACTURE_HYBRID_H_INCLUDED )
#define  KRATOS_PHASE_FIELD_APPLICATION_PHASE_FIELD_FRACTURE_HYBRID_H_INCLUDED


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
#include "phase_field_application/custom_utilities/isotropic_tensor_utility.h"


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
REF: D. Schillinger, M. Borden, H. K. Stolarski, Isogeometric collocation for phase-field fracture models
     M. Ambati, L. De Lorenzis, A review on phase-field models of brittle fracture and a new fast hybrid formulation
REMARKS:
    +   for this PhaseFieldFractureHybrid element, we should not enable MoveMeshFlag because Inverse of Jacobian is not computed a priori
    +   this element only works with tensile crack
    +   this element only works with staggerred scheme
    +   this element works for 3D and plane strain case
    +   this element supports the linear elastic material by default
 */
class PhaseFieldFractureHybrid : public Element
{

public:

    ///@name Type Definitions
    ///@{
    typedef GeometryData::IntegrationMethod IntegrationMethod;
    typedef GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef GeometryType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

    KRATOS_CLASS_POINTER_DEFINITION( PhaseFieldFractureHybrid );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PhaseFieldFractureHybrid( IndexType NewId, GeometryType::Pointer pGeometry );
    PhaseFieldFractureHybrid( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Destructor.
    virtual ~PhaseFieldFractureHybrid();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    IntegrationMethod GetIntegrationMethod() const {return mThisIntegrationMethod;}

    virtual Element::Pointer Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const;

    virtual Element::Pointer Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const;

    virtual void Initialize();

    virtual void EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo );

    virtual void GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo );

    virtual void InitializeSolutionStep( ProcessInfo& CurrentProcessInfo );

    virtual void InitializeNonLinearIteration( ProcessInfo& CurrentProcessInfo );

    virtual void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo );

    virtual void CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo); // this is required for calculating the reaction force

    virtual void FinalizeNonLinearIteration( ProcessInfo& CurrentProcessInfo );

    virtual void FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo );

    ///@}
    ///@name Access
    ///@{

    virtual void CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo );

    virtual void CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );

    virtual void CalculateOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo );

    virtual void GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo );

    virtual void GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo );

    virtual void GetValueOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo );

    virtual void SetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo );

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
//      virtual String Info() const;

    /// Print information about this object.
//      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;


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

    IntegrationMethod mThisIntegrationMethod;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    
    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    // A private default constructor necessary for serialization
    PhaseFieldFractureHybrid() {}

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
        rSerializer.save( "mIsInitialized", mIsInitialized );
        rSerializer.save( "mCurrentStresses", mCurrentStresses );
        rSerializer.save( "mReferenceEnergyDensity", mReferenceEnergyDensity );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
        rSerializer.load( "mIsInitialized", mIsInitialized );
        rSerializer.load( "mCurrentStresses", mCurrentStresses );
        rSerializer.load( "mReferenceEnergyDensity", mReferenceEnergyDensity );
    }

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo );


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

    std::vector<double> mReferenceEnergyDensity;
    std::vector<Vector> mCurrentStresses;
    bool mIsInitialized;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void CalculateAll( MatrixType& rLeftHandSideMatrix,
                       VectorType& rRightHandSideVector,
                       ProcessInfo& rCurrentProcessInfo,
                       bool CalculateStiffnessMatrixFlag = true,
                       bool CalculateResidualVectorFlag = true,
                       bool MaterialUpdateFlag = true );

    void CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX );
    void CalculateStrain( Vector& StrainVector, const Matrix& B, const Matrix& Displacements );
    void CalculateElasticMatrix( Matrix& C, double E, double NU );
    void CalculateSecondDerivatives( std::vector<Vector>& rD2N_DX2, const CoordinatesArrayType& rPoint );

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
    //PhaseFieldFractureHybrid& operator=(const PhaseFieldFractureHybrid& rOther);

    /// Copy constructor.
    //PhaseFieldFractureHybrid(const PhaseFieldFractureHybrid& rOther);


    ///@}

}; // Class PhaseFieldFractureHybrid

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_PHASE_FIELD_APPLICATION_PHASE_FIELD_FRACTURE_H_INCLUDED defined 

