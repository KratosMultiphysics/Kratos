///    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Zhiming Guo
//                   Riccardo Rossi
//




// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_mpm_element.h"
#include "includes/element.h"
#include "lagrangian_mpm_application_variables.h"

#include "utilities/math_utils.h"

#include "geometries/geometry.h"
//#include "custom_geometries/meshless_geometry.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
UpdatedLagrangianMPMElement::UpdatedLagrangianMPMElement(
	IndexType NewId,
	GeometryType::Pointer pGeometry)
    : MeshlessBaseElement(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//************************************************************************************
//************************************************************************************
UpdatedLagrangianMPMElement::UpdatedLagrangianMPMElement(
	IndexType NewId,
	GeometryType::Pointer pGeometry,  
	PropertiesType::Pointer pProperties)
    : MeshlessBaseElement(NewId, pGeometry, pProperties)
{
}

Element::Pointer UpdatedLagrangianMPMElement::Create(
	IndexType NewId, 
	NodesArrayType const& ThisNodes,  
	PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UpdatedLagrangianMPMElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Element::Pointer UpdatedLagrangianMPMElement::Create(Element::IndexType NewId,
                           Element::GeometryType::Pointer pGeom,
                           PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return Element::Pointer(new UpdatedLagrangianMPMElement(NewId, pGeom, pProperties));
        KRATOS_CATCH("");
    }

UpdatedLagrangianMPMElement::~UpdatedLagrangianMPMElement()
{
}

//************************************************************************************
//************************************************************************************
void UpdatedLagrangianMPMElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo)
    {
         const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetDomainSize();
        unsigned int matrix_size = number_of_nodes * dim;
       if (rLeftHandSideMatrix.size1() != matrix_size)
	  rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
        if (rRightHandSideVector.size() != matrix_size)
	  rRightHandSideVector.resize(matrix_size, false);
        
        //obtain shape functions
        
        //obtaine delta_disp and store it in a Matrix
        
        //Compute deltaF = I + D delta_disp / D xn
        
        //compute B
        
        //Compute Ftot
        
        //Compute material response sigma_cauchy, C
        
        //compute RHS
        
        //compute LHS
                
    }

    void UpdatedLagrangianMPMElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetDomainSize();
        unsigned int matrix_size = number_of_nodes * dim;
        
       if (rMassMatrix.size1() != matrix_size)
	  rMassMatrix.resize(matrix_size, matrix_size, false);
       
       //set it to zero
       rMassMatrix.clear();
       
       double integration_weight;
       Vector N;
       Matrix DN_Dx;
       this->GetGeometryData(integration_weight,N,DN_Dx);
       
       //fill the matrix
       for(unsigned int i=0; i<number_of_nodes; i++)
       {
           for(unsigned int j=0; j<number_of_nodes; j++)
            {
                for(unsigned int k=0; k<dim; k++)
                {
                    rMassMatrix(i*dim+k, j*dim+k) += integration_weight*N[i]*N[j];
                }
            }
       }
        
    }
    
    void UpdatedLagrangianMPMElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetDomainSize();
        unsigned int matrix_size = number_of_nodes * dim;

        if ( rResult.size() != matrix_size )
            rResult.resize( matrix_size, false );

        for ( int i = 0; i < number_of_nodes; i++ )
        {
            int index = i * dim;
            rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

            if ( dim == 3 )
                rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
        }

    }

//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianMPMElement::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo )
    {
        ElementalDofList.resize( 0 );

        for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

            if ( GetDomainSize() == 3 )
            {
                ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
            }
        }
    }
    
    void UpdatedLagrangianMPMElement::GetValuesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetDomainSize();
        unsigned int matrix_size = number_of_nodes * dim;

        if ( values.size() != matrix_size ) values.resize( matrix_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = i * dim;
            const array_1d<double,3>& disp = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT, Step );
            values[index] = disp[0];
            values[index + 1] = disp[1];;

            if ( dim == 3 )
                values[index + 2] = disp[2];;
        }
    }


//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianMPMElement::GetFirstDerivativesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetDomainSize();
        unsigned int matrix_size = number_of_nodes * dim;

        if ( values.size() != matrix_size ) values.resize( matrix_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = i * dim;
            const array_1d<double,3>& vel = GetGeometry()[i].GetSolutionStepValue( VELOCITY, Step );
            values[index] = vel[0];
            values[index + 1] = vel[1];

            if ( dim == 3 )
                values[index + 2] = vel[2];
        }
    }

//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianMPMElement::GetSecondDerivativesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetDomainSize();
        unsigned int matrix_size = number_of_nodes * dim;

        if ( values.size() != matrix_size ) values.resize( matrix_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = i * dim;
            const array_1d<double,3>& acc = GetGeometry()[i].GetSolutionStepValue( ACCELERATION, Step );
            values[index] = acc[0];
            values[index + 1] = acc[1];

            if ( dim == 3 )
                values[index + 2] = acc[2];
        }
    }
    
//************************************************************************************
//************************************************************************************
} // Namespace Kratos


