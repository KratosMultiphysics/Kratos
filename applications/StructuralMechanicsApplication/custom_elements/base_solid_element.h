// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_BASE_SOLID_ELEMENT_H_INCLUDED )
#define  KRATOS_BASE_SOLID_ELEMENT_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "utilities/integration_utilities.h"
#include "structural_mechanics_application_variables.h"

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
    
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  BaseSolidElement
    : public Element
{
public:

    ///@name Type Definitions
    ///@{
    
    // Counted pointer of BaseSolidElement
    KRATOS_CLASS_POINTER_DEFINITION( BaseSolidElement );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    BaseSolidElement()
    {};

    // Constructor using an array of nodes
    BaseSolidElement( IndexType NewId, GeometryType::Pointer pGeometry ):Element(NewId,pGeometry)
    {};

    // Constructor using an array of nodes with properties
    BaseSolidElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ):Element(NewId,pGeometry,pProperties)
    {};

    // Destructor
    virtual ~BaseSolidElement()
    {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     * @param rResult: The vector containing the equation id
     * @param rCurrentProcessInfo: The current process info instance
     */
    virtual void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo 
        );

    /**
     * Sets on rElementalDofList the degrees of freedom of the considered element geometry
     * @param rElementalDofList: The vector containing the dof of the element
     * @param rCurrentProcessInfo: The current process info instance
     */
    virtual void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo 
        );

    /**
     * Sets on rValues the nodal displacements
     * @param rValues: The values of displacements
     * @param Step: The step to be computed
     */
    virtual void GetValuesVector(
        Vector& rValues,
        int Step = 0 
        );
    
    /**
     * Sets on rValues the nodal velocities
     * @param rValues: The values of velocities
     * @param Step: The step to be computed
     */
    virtual void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0 
        );

    /**
     * Sets on rValues the nodal accelerations
     * @param rValues: The values of accelerations
     * @param Step: The step to be computed
     */
    virtual void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0 
        );

    /**
      * This is called during the assembling process in order to calculate the elemental mass matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    virtual void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo 
        );

   /**
      * This is called during the assembling process in order
      * to calculate the elemental damping matrix
      * @param rDampingMatrix: the elemental damping matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    virtual void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo 
        );

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo );

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{
    
protected:
    
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; // The vector containing the constitutive laws

    ///@}
    ///@name Protected Operators
    ///@{
    
    void  CalculateDerivativesOnReference(
        Matrix& J0, 
        Matrix& InvJ0, 
        Matrix& DN_DX, 
        double& detJ0, 
        const Matrix& DN_De
        )
    {
        J0.clear();
        
        for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
        {
            const array_1d<double, 3>& coords = GetGeometry()[i].GetInitialPosition(); //NOTE: here we refer to the original, undeformed position!!
            for(unsigned int k = 0; k < GetGeometry().WorkingSpaceDimension(); k++)
            {
                for(unsigned int m = 0; m < GetGeometry().LocalSpaceDimension(); m++)
                {
                    J0(k,m) += coords[k]*DN_De(i,m);
                }
            }
        }
        
        MathUtils<double>::InvertMatrix( J0, InvJ0, detJ0 );
        
        noalias( DN_DX ) = prod( DN_De, InvJ0);
    }


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

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{



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
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    }

}; // class BaseSolidElement.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_BASE_SOLID_ELEMENT_H_INCLUDED  defined 
