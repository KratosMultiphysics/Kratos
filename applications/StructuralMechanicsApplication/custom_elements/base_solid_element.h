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

    virtual void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        if (rResult.size() != dim*NumberOfNodes)
            rResult.resize(dim*NumberOfNodes,false);

        const unsigned int pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        if(dim == 2)
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                const unsigned int index = i * 2;
                rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
            }
        }
        else
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                const unsigned int index = i * 3;
                rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X,pos).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
                rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z,pos+2).EquationId();
            }
        }
        KRATOS_CATCH("")
    };


    virtual void GetDofList(
        DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim =  GetGeometry().WorkingSpaceDimension();
        ElementalDofList.resize(0);
        ElementalDofList.reserve(dim*NumberOfNodes);

        if(dim == 2)
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            }
        }
        else
        {
            for (unsigned int i = 0; i < NumberOfNodes; ++i)
            {
                ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
                ElementalDofList.push_back( GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            }
        }
        KRATOS_CATCH("")
    };

    virtual void GetValuesVector(
        Vector& values,
        int Step = 0 )
    {
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = NumberOfNodes * dim;
        if (values.size() != MatSize)
            values.resize(MatSize, false);
        for (unsigned int i = 0; i < NumberOfNodes; i++)
        {
            const array_1d<double, 3 > & Displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            unsigned int index = i * dim;
            for(unsigned int k=0; k<dim; ++k)
            {
                values[index+k] = Displacement[k];
            }
        }
    }


    virtual void GetFirstDerivativesVector(
        Vector& values,
        int Step = 0 )
    {
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = NumberOfNodes * dim;
        if (values.size() != MatSize)
        {
            values.resize(MatSize, false);
        }
        
        for (unsigned int i = 0; i < NumberOfNodes; i++)
        {
            const array_1d<double, 3 > & Velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            unsigned int index = i * dim;
            for(unsigned int k=0; k<dim; ++k)
            {
                values[index+k] = Velocity[k];
            }
        }
    }

    virtual void GetSecondDerivativesVector(
        Vector& values,
        int Step = 0 )
    {
        const unsigned int NumberOfNodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = NumberOfNodes * dim;
        if (values.size() != MatSize)
            values.resize(MatSize, false);
        for (unsigned int i = 0; i < NumberOfNodes; i++)
        {
            const array_1d<double, 3 > & Acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            unsigned int index = i * dim;
            for(unsigned int k=0; k<dim; ++k)
            {
                values[index+k] = Acceleration[k];
            }
        }
    }

    virtual void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        //lumped
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int NumberOfNodes = GetGeometry().size();
        unsigned int MatSize = dimension * NumberOfNodes;

        if ( rMassMatrix.size1() != MatSize )
            rMassMatrix.resize( MatSize, MatSize, false );

        rMassMatrix = ZeroMatrix( MatSize, MatSize );
        
        Matrix DN_DX( NumberOfNodes, dimension );
        Matrix J0(dimension,dimension), InvJ0(dimension,dimension);
        
        //reading integration points and local gradients
        IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(GetGeometry());
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( integration_method );
        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(integration_method);
        
        const double density = GetProperties()[DENSITY];
        double thickness = 1.0;
        if ( dimension == 2 && GetProperties().Has( THICKNESS )) 
            thickness = GetProperties()[THICKNESS];

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {

            const Matrix& DN_De = GetGeometry().ShapeFunctionsLocalGradients(integration_method)[PointNumber];

            double detJ0;
            CalculateDerivativesOnReference(J0, InvJ0, DN_DX, detJ0, DN_De);
            const double weight = integration_points[PointNumber].Weight()*detJ0*thickness;
            const Vector& N = row(Ncontainer,PointNumber);

            for ( unsigned int i = 0; i < NumberOfNodes; i++ )
            {
                unsigned int index_i = i * dimension;

                for ( unsigned int j = 0; j < NumberOfNodes; ++j )
                {
                    unsigned int index_j = j*dimension;
                    double NiNj_weight = N[i]*N[j]*weight*density;

                    for ( unsigned int k = 0; k < dimension; k++ )
                    {
                        
                        rMassMatrix( index_i+k, index_j+k ) += NiNj_weight;
                    }
                }
            }
        }



        KRATOS_CATCH( "" )
    }


    virtual void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        unsigned int NumberOfNodes = GetGeometry().size();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        //resizing as needed the LHS
        unsigned int MatSize = NumberOfNodes * dim;

        if ( rDampingMatrix.size1() != MatSize )
            rDampingMatrix.resize( MatSize, MatSize, false );

        noalias( rDampingMatrix ) = ZeroMatrix( MatSize, MatSize );

        KRATOS_CATCH( "" )
    }

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo )
    {
        if ( DISPLACEMENT.Key() == 0 )
        {
            KRATOS_ERROR <<  "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
        }

        //verify that the dofs exist
        for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
        {
            if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            {
                KRATOS_ERROR << "missing variable DISPLACEMENT on node " << this->GetGeometry()[i].Id() << std::endl;
            }

            if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false ||
                 this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false ||
                 this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
            {
                KRATOS_ERROR << "missing one of the dofs for the variable DISPLACEMENT on node " << GetGeometry()[i].Id() << " of condition " << Id() << std::endl;
            }
        }
        return 0;
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

    ///@}
    ///@name Friends
    ///@{
    
protected:
    
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    ///@}
    ///@name Protected Operators
    ///@{
    void  CalculateDerivativesOnReference(Matrix& J0, 
                             Matrix& InvJ0, 
                             Matrix& DN_DX, 
                             double& detJ0, 
                             const Matrix& DN_De
                                         )
    {
        J0.clear();
        for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
        {
            const auto& coords = GetGeometry()[i].GetInitialPosition(); //NOTE: here we refer to the original, undeformed position!!
            for(unsigned int k=0; k<GetGeometry().WorkingSpaceDimension(); k++)
            {
                for(unsigned int m=0; m<GetGeometry().LocalSpaceDimension(); m++)
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
