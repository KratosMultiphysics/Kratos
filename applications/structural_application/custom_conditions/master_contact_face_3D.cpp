/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
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
/* *********************************************************
*
*   Last Modified by:    $Author: hurga $
*   Date:                $Date: 2009-03-17 14:35:29 $
*   Revision:            $Revision: 1.4 $
*
* ***********************************************************/

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/master_contact_face_3D_newmark.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    MasterContactFace3D::MasterContactFace3D( IndexType NewId,
            GeometryType::Pointer pGeometry ) :
            Condition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
	GetValue( IS_CONTACT_SLAVE  )  = 0;
	GetValue( IS_CONTACT_MASTER )  = 1;  
    }

    //************************************************************************************
    //**** life cycle ********************************************************************
    //************************************************************************************
    MasterContactFace3D::MasterContactFace3D( IndexType NewId, GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pProperties ) :
            Condition( NewId, pGeometry, pProperties )
    {
//         for(unsigned int i = 0 ; i != GetGeometry().size() ; ++i)
//         {
//             (GetGeometry()[i].pAddDof(DISPLACEMENT_X));
//             (GetGeometry()[i].pAddDof(DISPLACEMENT_Y));
//             (GetGeometry()[i].pAddDof(DISPLACEMENT_Z));
//         }
        GetValue( LAMBDAS ).resize( GetGeometry().IntegrationPoints().size(), false );
        noalias( GetValue( LAMBDAS ) ) = ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( LAMBDAS_T ).resize( GetGeometry().IntegrationPoints().size(), 2, false );
        noalias( GetValue( LAMBDAS_T ) ) = ZeroMatrix( GetGeometry().IntegrationPoints().size(), 2 );
        GetValue( GAPS ).resize( GetGeometry().IntegrationPoints().size(), false );
        noalias( GetValue( GAPS ) ) = ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( DELTA_LAMBDAS ).resize( GetGeometry().IntegrationPoints().size() , false );
        noalias( GetValue( DELTA_LAMBDAS ) ) = ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( DELTA_LAMBDAS_T ).resize( GetGeometry().IntegrationPoints().size(), 2, false );
        noalias( GetValue( DELTA_LAMBDAS_T ) ) = ZeroMatrix( GetGeometry().IntegrationPoints().size(), 2 );
        GetValue( PENALTY ).resize( GetGeometry().IntegrationPoints().size(), false );
        noalias( GetValue( PENALTY ) ) =  ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( PENALTY_T ).resize( GetGeometry().IntegrationPoints().size(), false );
        noalias( GetValue( PENALTY_T ) ) = ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( IS_CONTACT_MASTER ) = 1;
        GetValue( IS_CONTACT_SLAVE ) = 0;
    }

    Condition::Pointer MasterContactFace3D::Create( IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties ) const
    {
        return Condition::Pointer( new MasterContactFace3D( NewId, GetGeometry().Create( ThisNodes ),
                                   pProperties ) );
    }

    /**
     * Destructor. Never to be called manually
     */
    MasterContactFace3D::~MasterContactFace3D()
    {
    }

    //************************************************************************************
    //************************************************************************************
    /**
     * returns condition type info
     */
//     int MasterContactFace3D::IsContactType()
//     {/*
//         return( 2 );*/
//     }
    //************************************************************************************
    //************************************************************************************
    /**
     * returns the tangential vectors of the current surface in an arbitrary point
     * TODO: TO BE REIMPLEMENTED!!!
     */
//     Matrix MasterContactFace3D::TangentialVectors( GeometryType::CoordinatesArrayType& rPoint )
//     {
//     }
    //************************************************************************************
    //************************************************************************************

    /**
     * returns normal vector in arbitrary point
     * calculates the normalized vector orthogonal to the current surface in given point
     * @param rPoint the given point in local coordinates
     * @return the normal vector
     */
//     Vector MasterContactFace3D::NormalVector( GeometryType::CoordinatesArrayType& rPoint )
//     {
//     }

    /**
     * returns closest point on current condition element with regard to given point in global coordinates
     * @param rResult a Point in global coordinates being overwritten by the desired information
     * @param rLocalResult a Point in local coordinates being overwritten by the desired information
     * @param rCandidate a Node of the current surface in global coordinates which lies closest to rPoint
     * @param rPoint the point in global coordinates the closest point on the current condition element is to
     * be calculated for
     * @return true if an orthogonal projection of the given point lies within the boundaries of the current
     * condition element
     */
    /**
    * returns closest point on current condition element with regard to given point in global coordinates
    * @param rResultGlobal a Point in global coordinates being overwritten by the desired information
    * @param rResultLocal a Point in global coordinates being overwritten by the desired information
    * @param rSlaveContactGlobalPoint the point in global coordinates the closest point on the current condition element is to
    * @param rCandidateGlobal the closest node to rSlaveContactGlobalPoint on current
    * surface
    * be calculated for
    * @return true if an orthogonal projection of the given point lies within the boundaries of the current
    * condition element
     */
    bool MasterContactFace3D::ClosestPoint( GeometryType::CoordinatesArrayType& rResultGlobal,
                                            GeometryType::CoordinatesArrayType& rResultLocal,
                                            const GeometryType::CoordinatesArrayType& rSlaveContactGlobalPoint,
                                            const GeometryType::CoordinatesArrayType& rCandidateGlobal
                                          )
    {
        return false;
    }

    //************************************************************************************
    //************************************************************************************

    /**
     * calculates only the RHS vector (certainly to be removed due to contact algorithm)
     */
    void MasterContactFace3D::CalculateRightHandSide( VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo )
    {
        unsigned int ndof = GetGeometry().size() * 3;

        if ( rRightHandSideVector.size() != ndof )
            rRightHandSideVector.resize( ndof, false );

        rRightHandSideVector = ZeroVector( ndof );

    }

    //************************************************************************************
    //************************************************************************************

    /**
     * calculates this contact element's local contributions
     */
    void MasterContactFace3D::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo )
    {
        unsigned int ndof = GetGeometry().size() * 3;

        if ( rRightHandSideVector.size() != ndof )
            rRightHandSideVector.resize( ndof, false );

        rRightHandSideVector = ZeroVector( ndof );

        if ( rLeftHandSideMatrix.size1() != ndof )
            rLeftHandSideMatrix( ndof, ndof );

        rLeftHandSideMatrix = ZeroMatrix( ndof, ndof );
    }

    //************************************************************************************
    //************************************************************************************
    /**
     * calculates the contact related contributions to the system
     * Does nothing as assembling is to be switched to linking objects
     */
    void MasterContactFace3D::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo,
                                            bool CalculateStiffnessMatrixFlag,
                                            bool CalculateResidualVectorFlag )
    {
    }

    //***********************************************************************
    //***********************************************************************
    /**
     * System matrix contribution due to contact energy
     * TO BE IMPLEMENTED
     */
    void MasterContactFace3D::CalculateAndAddKc( Matrix& K,
            const Vector& N,
            double weight,
            double dA,
            double penalty,
            Vector v
                                               )
    {
    }

    //***********************************************************************
    //***********************************************************************
    /**
     *
     */
    void MasterContactFace3D::CalculateAndAdd_PressureForce( Vector& residualvector,
            const Vector& N,
            Vector& v3,
            double pressure,
            double weight, double dA
                                                           )
    {
    }

    //************************************************************************************
    //************************************************************************************
    /**
     * REMOVED: the DOFs are managed by the linking conditions
     */
    void MasterContactFace3D::EquationIdVector( EquationIdVectorType& rResult,
            ProcessInfo& CurrentProcessInfo
                                              )
    {
        KRATOS_TRY
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = number_of_nodes * 3;

        if ( rResult.size() != dim )
            rResult.resize( dim );

        for ( unsigned int i = 0;i < number_of_nodes;i++ )
        {
            int index = i * 3;
            rResult[index]   = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[index+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
            rResult[index+2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************
    /**
     * REMOVED: the DOFs are managed by the linking conditions
     */
    void MasterContactFace3D::GetDofList( DofsVectorType& ConditionalDofList,
                                          ProcessInfo& CurrentProcessInfo )
    {
        ConditionalDofList.resize( 0 );

        for ( unsigned int i = 0;i < GetGeometry().size();i++ )
        {
            ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
            ConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }

    void MasterContactFace3D::GetValueOnIntegrationPoints( const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
        const GeometryType::ShapeFunctionsGradientsType& sf_gradients = GetGeometry().ShapeFunctionsLocalGradients();

        if ( rVariable == NORMAL )
        {
            if ( rValues.size() != integration_points.size() )
                rValues.resize( integration_points.size() );

            for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
            {
                //setting up result matrix
                Matrix T( 2, 3 );
                noalias( T ) = ZeroMatrix( 2, 3 );
                //shape function gradients
                Matrix DN = sf_gradients[PointNumber];
                //calculating tangential vectors

                for ( unsigned int n = 0; n < GetGeometry().PointsNumber(); n++ )
                {
                    T( 0, 0 ) += ( GetGeometry().GetPoint( n ).X0()
                                   + GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) ) * DN( n, 0 );
                    T( 0, 1 ) += ( GetGeometry().GetPoint( n ).Y0()
                                   + GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) ) * DN( n, 0 );
                    T( 0, 2 ) += ( GetGeometry().GetPoint( n ).Z0()
                                   + GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) ) * DN( n, 0 );
                    T( 1, 0 ) += ( GetGeometry().GetPoint( n ).X0()
                                   + GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_X ) ) * DN( n, 1 );
                    T( 1, 1 ) += ( GetGeometry().GetPoint( n ).Y0()
                                   + GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Y ) ) * DN( n, 1 );
                    T( 1, 2 ) += ( GetGeometry().GetPoint( n ).Z()
                                   + GetGeometry().GetPoint( n ).GetSolutionStepValue( DISPLACEMENT_Z ) ) * DN( n, 1 );
                }

                Vector Result( 3 );

                //calculating normal vector
                Result[0] = T( 0, 1 ) * T( 1, 2 ) - T( 0, 2 ) * T( 1, 1 );
                Result[1] = T( 0, 2 ) * T( 1, 0 ) - T( 0, 0 ) * T( 1, 2 );
                Result[2] = T( 0, 0 ) * T( 1, 1 ) - T( 0, 1 ) * T( 1, 0 );
                SD_MathUtils<double>::Normalize( Result );

                rValues[PointNumber][0] = Result[0];
                rValues[PointNumber][1] = Result[1];
                rValues[PointNumber][2] = Result[2];
            }
        }
    }

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int MasterContactFace3D::Check( const Kratos::ProcessInfo& rCurrentProcessInfo )
    {
        return 0;
    }


    /// Print information about this object.

    void MasterContactFace3D::PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << "Condition #" << Id();
    }

    /// Print object's data.

    void MasterContactFace3D::PrintData( std::ostream& rOStream ) const
    {
        rOStream << "MasterContactFace3D" << std::endl;
    }


} // Namespace Kratos
