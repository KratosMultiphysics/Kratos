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
*   Date:                $Date: 2009-02-26 13:49:16 $
*   Revision:            $Revision: 1.5 $
*
* ***********************************************************/

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_conditions/slave_contact_face_3D_newmark.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************
    SlaveContactFace3D::SlaveContactFace3D( IndexType NewId, 
                                  GeometryType::Pointer pGeometry) : 
            Condition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }
    
    //************************************************************************************
    //**** life cycle ********************************************************************
    //************************************************************************************
    SlaveContactFace3D::SlaveContactFace3D( IndexType NewId, GeometryType::Pointer pGeometry,
                                  PropertiesType::Pointer pProperties) : 
            Condition( NewId, pGeometry, pProperties )
    {
//         for(unsigned int i = 0 ; i != GetGeometry().size() ; ++i)
//         {
//             (GetGeometry()[i].pAddDof(DISPLACEMENT_X));
//             (GetGeometry()[i].pAddDof(DISPLACEMENT_Y));
//             (GetGeometry()[i].pAddDof(DISPLACEMENT_Z));
//         }
        mpMasterElements = ContactMasterContainerType::Pointer( new ContactMasterContainerType() );
        GetValue( LAMBDAS ).resize(GetGeometry().IntegrationPoints().size(),false);
        noalias(GetValue( LAMBDAS )) = ZeroVector( GetGeometry().IntegrationPoints().size());
        GetValue( LAMBDAS_T ).resize( GetGeometry().IntegrationPoints().size(), 2 ,false);
        noalias(GetValue( LAMBDAS_T )) = ZeroMatrix( GetGeometry().IntegrationPoints().size(), 2 );
        GetValue( GAPS ).resize( GetGeometry().IntegrationPoints().size(),false);
        noalias(GetValue( GAPS )) = ZeroVector( GetGeometry().IntegrationPoints().size());
        GetValue( DELTA_LAMBDAS ).resize( GetGeometry().IntegrationPoints().size(),false );
        noalias(GetValue( DELTA_LAMBDAS )) = ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( DELTA_LAMBDAS_T ).resize( GetGeometry().IntegrationPoints().size(), 2,false );
        noalias(GetValue( DELTA_LAMBDAS_T )) = ZeroMatrix( GetGeometry().IntegrationPoints().size(), 2 );
        GetValue( PENALTY ).resize( GetGeometry().IntegrationPoints().size() ,false);
        noalias(GetValue( PENALTY )) =  ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( PENALTY_T ).resize( GetGeometry().IntegrationPoints().size() ,false);
		noalias(GetValue( PENALTY_T )) = ZeroVector( GetGeometry().IntegrationPoints().size() );
        GetValue( IS_CONTACT_SLAVE ) = 1;
		GetValue( IS_CONTACT_MASTER ) = 0;    

        GetValue(NORMAL_STRESS).resize(GetGeometry().IntegrationPoints().size(),false);
        GetValue(NORMAL_STRESS)=ZeroVector( GetGeometry().IntegrationPoints().size());
        GetValue(TANGENTIAL_STRESS).resize(GetGeometry().IntegrationPoints().size(),false);
        GetValue(TANGENTIAL_STRESS)=ZeroVector( GetGeometry().IntegrationPoints().size());
        GetValue(STICK).resize(GetGeometry().IntegrationPoints().size(),false);
        GetValue(STICK)=ZeroVector( GetGeometry().IntegrationPoints().size());

		GetValue(NORMAL_CONTACT_STRESS) = 0.0; 
		GetValue(TANGENTIAL_CONTACT_STRESS) = 0.0; 
		GetValue(CONTACT_STICK) = 0.0; 
    }
    
    Condition::Pointer SlaveContactFace3D::Create( IndexType NewId, 
                                              NodesArrayType const& ThisNodes,  
                                              PropertiesType::Pointer pProperties) const
    {
        return Condition::Pointer( new SlaveContactFace3D(NewId, GetGeometry().Create(ThisNodes), 
                                   pProperties));
    }
    /**
     * Destructor. Never to be called manually
     */
    SlaveContactFace3D::~SlaveContactFace3D()
    {
    }
    
    //************************************************************************************
    //************************************************************************************
    
    /**
     * returns the tangential vectors of the current surface in an arbitrary point
     */
//     Matrix SlaveContactFace3D::TangentialVectors( GeometryType::CoordinatesArrayType& rPoint )
//     {
//     }
    //************************************************************************************
    //************************************************************************************
    void SlaveContactFace3D::CalculateOnIntegrationPoints(const Variable<double>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
	{
         //reading integration points and local gradients
            const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

            if(Output.size() != integration_points.size())
                Output.resize(integration_points.size(),false);
	    double result= 0.0;
	    double result_friction= 0.0;
	    double reference= 0.0;
   		if(rVariable==NORMAL_CONTACT_STRESS || rVariable==TANGENTIAL_CONTACT_STRESS)
        {
            for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
            {  
			Matrix TSlave = TangentialVectors_inOrigin(integration_points[PointNumber]);

			Vector vSlaveNonNormalized = ZeroVector(3);

        		vSlaveNonNormalized[0] = TSlave(0,1)*TSlave(1,2)-TSlave(0,2)*TSlave(1,1);
        		vSlaveNonNormalized[1] = TSlave(0,2)*TSlave(1,0)-TSlave(0,0)*TSlave(1,2);
        		vSlaveNonNormalized[2] = TSlave(0,0)*TSlave(1,1)-TSlave(0,1)*TSlave(1,0);

        		double dASlave = MathUtils<double>::Norm3( vSlaveNonNormalized );

                    	result+= (this->GetValue(NORMAL_STRESS)[PointNumber])*(this->GetGeometry().IntegrationPoints()[PointNumber].Weight())*dASlave;
			result_friction+= (this->GetValue(TANGENTIAL_STRESS)[PointNumber])*(this->GetGeometry().IntegrationPoints()[PointNumber].Weight())*dASlave;
			reference+= this->GetGeometry().IntegrationPoints()[PointNumber].Weight()*dASlave;
	    }
	 }
	
            for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
            {  
                if(rVariable==NORMAL_CONTACT_STRESS)
                {

                    	Output[PointNumber]= result/reference;

                }
                if(rVariable==TANGENTIAL_CONTACT_STRESS)
                {
                    Output[PointNumber]=  result_friction/reference;
                }
                if(rVariable==CONTACT_PENETRATION)
                {
                    //Output[PointNumber]=  this->GetValue(STICK)[PointNumber] ;
                    Output[PointNumber]=  this->GetValue(GAPS)[PointNumber] ;
                }
	    }
	}
    
    void SlaveContactFace3D::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
        Vector result = ZeroVector(GetGeometry().IntegrationPoints().size());
        CalculateOnIntegrationPoints( rVariable, result, rCurrentProcessInfo );
        if( rValues.size() != GetGeometry().IntegrationPoints().size() )
            rValues.resize(GetGeometry().IntegrationPoints().size());
        for( unsigned int i=0; i< result.size(); i++ )
            rValues[i] = result[i];
    }
    //************************************************************************************
    //************************************************************************************
    /**
     * returns normal vector in arbitrary point 
     * calculates the normalized vector orthogonal to the current surface in given point
     * @param rPoint the given point in local coordinates
     * @return the normal vector 
     */
//     Vector SlaveContactFace3D::NormalVector( GeometryType::CoordinatesArrayType& rPoint )
//     {
//     }
  
    //************************************************************************************
    //************************************************************************************
   
   
    /**
     * calculates only the RHS vector (certainly to be removed due to contact algorithm)
     */
    void SlaveContactFace3D::CalculateRightHandSide( VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo)
    {
        rRightHandSideVector.resize(0,false);
    }
    
    //************************************************************************************
    //************************************************************************************
    /**
     * calculates this contact element's local contributions
     */
    void SlaveContactFace3D::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, 
                                              VectorType& rRightHandSideVector, 
                                              ProcessInfo& rCurrentProcessInfo)
    {
        rRightHandSideVector.resize(0,false);
        rLeftHandSideMatrix(0,0);
    }

    //************************************************************************************
    //************************************************************************************
    /**
     * calculates the contact related contributions to the system
     * does nothing as assembling is switched to link condition
     */
    void SlaveContactFace3D::CalculateAll( MatrixType& rLeftHandSideMatrix, 
                                      VectorType& rRightHandSideVector, 
                                      ProcessInfo& rCurrentProcessInfo,
                                      bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag)
    {
    } // CalculateAll
    
    //***********************************************************************
    //***********************************************************************
    /**
     * System matrix contribution due to contact energy
     * TODO: implement mixed elementary contribution
     */
    void SlaveContactFace3D::CalculateAndAddKc( Matrix& K,
            const Vector& N,
            double weight,
            double dA,
            Vector v )
    {
    }
    
    //***********************************************************************
    //***********************************************************************
    /**
     * TO BE TESTED!!!
     */
    void SlaveContactFace3D::CalculateAndAdd_PressureForce( Vector& residualvector,
            const Vector& N,
            Vector& v3,
            double pressure,
            double weight,
            double DetJ )
    {
    }
    
    //************************************************************************************
    //************************************************************************************
    /**
     * REMOVED: the DOFs are managed by the linking conditions
     */
    void SlaveContactFace3D::EquationIdVector( EquationIdVectorType& rResult, 
                                          ProcessInfo& CurrentProcessInfo )
    {
        KRATOS_TRY
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = number_of_nodes*3;
        if(rResult.size() != dim)
            rResult.resize(dim);

        for (unsigned int i=0;i<number_of_nodes;i++)
        {
            int index = i*3;
            rResult[index]   = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        }
        KRATOS_CATCH("")
    }
    //************************************************************************************
    //************************************************************************************
    void SlaveContactFace3D::MasterElementsEquationIdVectors(EquationIdVectorContainerType& rResult,
            ProcessInfo& rCurrentProcessInfo )
    {
    } 
    
    //************************************************************************************
    //************************************************************************************
    /**
     * REMOVED: the DOFs are managed by the linking conditions
     */
    void SlaveContactFace3D::GetDofList( DofsVectorType& ConditionalDofList,
                                    ProcessInfo& CurrentProcessInfo)
    {
        ConditionalDofList.resize(0);

        for (unsigned int i=0;i<GetGeometry().size();i++)
        {
            ConditionalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            ConditionalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            ConditionalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }
    }

    Matrix SlaveContactFace3D::TangentialVectors_inOrigin( const GeometryType::CoordinatesArrayType& rPoint )
    {
        //setting up result matrix
        Matrix T = ZeroMatrix( 2, 3 );
        //shape function gradients
        Matrix DN = ZeroMatrix( this->GetGeometry().PointsNumber(),2);
        this->GetGeometry().ShapeFunctionsLocalGradients( DN, rPoint );
        //calculating tangential vectors
        for( unsigned int n=0; n<this->GetGeometry().PointsNumber(); n++ )
        {
            T(0,0) += this->GetGeometry().GetPoint(n).X0()*DN(n,0);
            T(0,1) += this->GetGeometry().GetPoint(n).Y0()*DN(n,0);
            T(0,2) += this->GetGeometry().GetPoint(n).Z0()*DN(n,0);
            T(1,0) += this->GetGeometry().GetPoint(n).X0()*DN(n,1);
            T(1,1) += this->GetGeometry().GetPoint(n).Y0()*DN(n,1);
            T(1,2) += this->GetGeometry().GetPoint(n).Z0()*DN(n,1);
        }
        return( T );
    }
} // Namespace Kratos
