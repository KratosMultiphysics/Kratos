//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2013-05-23 00:32:00 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_BATHE_RECOVER_STRESS_UTILITY_H_INCLUDED)
#define  KRATOS_BATHE_RECOVER_STRESS_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <cmath>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "containers/weak_pointer_vector.h"
#include "includes/serializer.h"
#include "sd_math_utils.h"
#include "structural_application.h"


//#define UTILITY_DEBUG_LEVEL1


namespace Kratos
{

/*
 * This class implements recovery stress procedure as described in the paper "A stress improvement procedure" by Payen & Bathe
 */
class BatheRecoverStressUtility
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(BatheRecoverStressUtility);

    typedef double CoordinateType;
    typedef std::size_t IndexType;

    BatheRecoverStressUtility(unsigned int ExpansionLevel) : mExpansionLevel(ExpansionLevel) {}
    virtual ~BatheRecoverStressUtility() {}
    
    
    //recovery stress routine for 2d
    void CalculateImprovedStressOnIntegrationPoints( Element& rCurrentElement, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        //check dimension
        if( rCurrentElement.GetGeometry().WorkingSpaceDimension() != 2 )
            KRATOS_THROW_ERROR(std::logic_error, "This recovery routine currently works in 2d only", "");
    
        // build list of neighbours of the current element
        WeakPointerVector<Element>& neighb_elems = GetNeighborElements( rCurrentElement, mExpansionLevel );
        
        // from list of neighbour elements iterate through all the integration points and calculate the left hand side and right hand side contribution
        
        #ifdef UTILITY_DEBUG_LEVEL1
        std::cout << rCurrentElement.Id() << ":";
        for( WeakPointerVector<Element>::iterator ie = neighb_elems.begin(); ie != neighb_elems.end(); ++ie )
            std::cout << " " << ie->Id();
        std::cout << std::endl;
        #endif
        
        Matrix LeftHandSideMatrix( 18, 18 );
        Vector RightHandSideVector( 18 );
        
        noalias( LeftHandSideMatrix ) = ZeroMatrix( 18, 18 );
        noalias( RightHandSideVector ) = ZeroVector( 18 );
        
        Matrix Et_Operator( 3, 18 );
        Matrix DtEt_Operator( 2, 18 );
        Matrix Ebart_Operator( 3, 12 );
        Matrix Exi_Operator( 2, 6 );
        Matrix LHSPart1( 12, 18 );
        Matrix LHSPart2( 6, 18 );
        Vector RHSPart1( 12 );
        Vector RHSPart2( 6 );
        
        
        //calculate left hand side and right hand side
        for( WeakPointerVector<Element>::iterator ie = neighb_elems.begin(); ie != neighb_elems.end(); ++ie )
        {
            
            //contribution of left hand side
            noalias( LHSPart1 ) = ZeroMatrix( 12, 18 );
            noalias( LHSPart2 ) = ZeroMatrix( 6, 18 );
            CalculateLHS( ie, Et_Operator, DtEt_Operator, Ebart_Operator, Exi_Operator, LHSPart1, LHSPart2, rCurrentProcessInfo );
            
            //contribution of right hand side
            noalias( RHSPart1 ) = ZeroVector( 12 );
            noalias( RHSPart2 ) = ZeroVector( 6 );
            CalculateRHS( ie, Et_Operator, DtEt_Operator, Ebart_Operator, Exi_Operator, RHSPart1, RHSPart2, rCurrentProcessInfo );
            
            //assemble
            for(unsigned int i = 0; i < 12; i++)
            {
                for(unsigned int j = 0; j < 18; j++)
	    	    {
	    	        LeftHandSideMatrix( i, j ) += LHSPart1( i, j );
                }
                RightHandSideVector( i ) += RHSPart1( i );
            }
            
            for(unsigned int i = 0; i < 6; i++)
            {
                for(unsigned int j = 0; j < 18; j++)
	    	    {
	    	        LeftHandSideMatrix( 12 + i, j ) += LHSPart2( i, j );
                }
                RightHandSideVector( 12 + i ) += RHSPart2( i );
            }
        
        }
        
        
        #ifdef UTILITY_DEBUG_LEVEL1
        KRATOS_WATCH( LeftHandSideMatrix )
        KRATOS_WATCH( RightHandSideVector )
        #endif
        
        //calculate stress coefficients vector
        Matrix Inverse_LeftHandSideMatrix( 18, 18 );
        Vector StressCoefficients( 18 );
        
        int singular = SD_MathUtils<double>::InvertMatrix( LeftHandSideMatrix, Inverse_LeftHandSideMatrix );
        //lu_factorize() [SD_MathUtils<double>::InvertMatrix] returns 0 if it was successful. It returns (k+1) if it detects singularity after processing row k. So one should always check its return value.
        if( singular != 0 )
            KRATOS_THROW_ERROR(std::logic_error, "Singular matrix detected when recover stress for element", rCurrentElement.Id());
        noalias( StressCoefficients ) = prod( Inverse_LeftHandSideMatrix, RightHandSideVector );
        
        
        //recover stress at each integration points
        const Element::GeometryType::IntegrationPointsArrayType& integration_points = rCurrentElement.GetGeometry().IntegrationPoints( rCurrentElement.GetIntegrationMethod() );
        
        const Matrix& Ncontainer = rCurrentElement.GetGeometry().ShapeFunctionsValues( rCurrentElement.GetIntegrationMethod() );
        
        if ( rValues.size() != integration_points.size() )
			rValues.resize( integration_points.size() );
	    
	    unsigned int StrainSize = 3;

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
		{
		    if ( rValues[PointNumber].size() != StrainSize )
                rValues[PointNumber].resize( StrainSize );
		    
		    Vector Nvector = row( Ncontainer, PointNumber );

		    double RealX = 0.0;
		    double RealY = 0.0;

            for(unsigned int i = 0; i < rCurrentElement.GetGeometry().size(); i++)
		    {
                RealX += Nvector( i ) * rCurrentElement.GetGeometry()[i].X0();
                RealY += Nvector( i ) * rCurrentElement.GetGeometry()[i].Y0();
		    }
		    
		    CalculateEtOperator( Et_Operator, RealX, RealY );
		    noalias( rValues[PointNumber] ) = prod( Et_Operator, StressCoefficients );
		    
		}

    }
    
    
    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "A utility to re-compute the stress field in finite element" << std::endl;
        buffer << "Reference: Bathe, Payen, A stress improvement procedure, 2012" << std::endl;
        return buffer.str();
    }


    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }


    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }
    
    
private:
    
    unsigned int mExpansionLevel;
    
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        //TODO
        rSerializer.save("BatheRecoverStressUtility", *this);
    }

    virtual void load(Serializer& rSerializer)
    {
        //TODO
        rSerializer.load("BatheRecoverStressUtility", *this);
    }
    
    
    WeakPointerVector<Element>& GetNeighborElements(Element& rCurrentElement, unsigned int ExpansionLevel)
    {
        
        if( ExpansionLevel < 1 )
            KRATOS_THROW_ERROR(std::logic_error, "The expansion level of an element must be >= 1, detected", ExpansionLevel);
    
        if( ExpansionLevel == 1 )
        {
            WeakPointerVector<Element>& neighb_elems = rCurrentElement.GetValue(NEIGHBOUR_ELEMENTS);
            
            bool CurrentElementIsIncluded = false;
            
            for( unsigned int i = 0; i < rCurrentElement.GetGeometry().size(); i++ )
            {
                WeakPointerVector<Element>& tmp_elems = rCurrentElement.GetGeometry()[i].GetValue(NEIGHBOUR_ELEMENTS);
                    
                for( WeakPointerVector<Element>::iterator i_tmp = tmp_elems.begin(); i_tmp != tmp_elems.end(); ++i_tmp )
                {
                    if( i_tmp->Id() != rCurrentElement.Id() )
                    {
                        bool isCounted = false;
                        
                        for( WeakPointerVector<Element>::iterator i_tmp1 = neighb_elems.begin(); i_tmp1 != neighb_elems.end(); ++i_tmp1 )
                        {
                            if( i_tmp->Id() == i_tmp1->Id() )
                            {
                                isCounted = true;
                                break;
                            }
                        }
                        
                        if( isCounted == false )
                        {
                            neighb_elems.push_back( Element::WeakPointer( *(i_tmp.base()) ) );
                        }
                    }
                    else
                    {
                        if( !CurrentElementIsIncluded )
                        {
                            CurrentElementIsIncluded = true;
                            neighb_elems.push_back( Element::WeakPointer( *(i_tmp.base()) ) );
                        }
                    }
                }
            }
            
            return neighb_elems;
        }
        else
        {
            WeakPointerVector<Element>& neighb_elems = GetNeighborElements( rCurrentElement, ExpansionLevel - 1 );
            
            WeakPointerVector<Element> more_elems;
            
            for( WeakPointerVector<Element>::iterator ie = neighb_elems.begin(); ie != neighb_elems.end(); ++ie )
            {
                for( unsigned int i = 0; i < ie->GetGeometry().size(); i++ )
                {
                    WeakPointerVector<Element>& tmp_elems = ie->GetGeometry()[i].GetValue(NEIGHBOUR_ELEMENTS);
                    
                    for( WeakPointerVector<Element>::iterator i_tmp = tmp_elems.begin(); i_tmp != tmp_elems.end(); ++i_tmp )
                    {
                    
                        if( i_tmp->Id() != rCurrentElement.Id() )
                        {
                            bool isCounted = false;
                            
                            for( WeakPointerVector<Element>::iterator i_tmp1 = neighb_elems.begin(); i_tmp1 != neighb_elems.end(); ++i_tmp1 )
                            {
                                if( i_tmp->Id() == i_tmp1->Id() )
                                {
                                    isCounted = true;
                                    break;
                                }
                            }
                            
                            if( isCounted == false )
                            {
                                for( WeakPointerVector<Element>::iterator i_tmp2 = more_elems.begin(); i_tmp2 != more_elems.end(); ++i_tmp2 )
                                {
                                    if( i_tmp->Id() == i_tmp2->Id() )
                                    {
                                        isCounted = true;
                                        break;
                                    }
                                }
                            
                                if( isCounted == false )
                                {
                                    more_elems.push_back( Element::WeakPointer( *(i_tmp.base()) ) );
                                }
                            }
                        }
                    }
                }
            }
            
            for( WeakPointerVector<Element>::iterator i_tmp = more_elems.begin(); i_tmp != more_elems.end(); ++i_tmp )
            {
                neighb_elems.push_back( Element::WeakPointer( *(i_tmp.base()) ) );
            }
            
            return neighb_elems;
        }
    }
    
    
    void CalculateLHS( const WeakPointerVector<Element>::iterator& ie,
                        Matrix& Et_Operator,
                        Matrix& DtEt_Operator,
                        Matrix& Ebart_Operator,
                        Matrix& Exi_Operator,
                        Matrix& LHSPart1,
                        Matrix& LHSPart2,
                        const ProcessInfo& rCurrentProcessInfo )
    {
        GeometryData::IntegrationMethod ThisIntegrationMethod = ie->GetIntegrationMethod();
            
        //modify integration rule in case of 3-node triangle
        if( ie->GetGeometry().size() == 3 )
            ThisIntegrationMethod = GeometryData::GI_GAUSS_2;
            
        const Element::GeometryType::IntegrationPointsArrayType& integration_points = ie->GetGeometry().IntegrationPoints( ThisIntegrationMethod );
            
        const Matrix& Ncontainer = ie->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

        //calculate the Jacobian first
        Element::GeometryType::JacobiansType J0( integration_points.size() );
		J0 = ie->GetGeometry().Jacobian( J0, ThisIntegrationMethod );

        Vector DetJ0;
		DetJ0.resize( integration_points.size(), false );
		noalias( DetJ0 ) = ZeroVector( integration_points.size() );
            
        int dim = ie->GetGeometry().WorkingSpaceDimension();
        Matrix dummy(dim, dim);
		for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
		    MathUtils<double>::InvertMatrix( J0[PointNumber], dummy, DetJ0[PointNumber] );

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
		{
            //double xi = integration_points[PointNumber].X();
            //double eta = integration_points[PointNumber].Y();
            double weight = integration_points[PointNumber].Weight();
            
            Vector Nvector = row( Ncontainer, PointNumber );
            
            double RealX = 0.0;
            double RealY = 0.0;
            
            for(unsigned int i = 0; i < ie->GetGeometry().size(); i++)
            {
                RealX += Nvector( i ) * ie->GetGeometry()[i].X0();
                RealY += Nvector( i ) * ie->GetGeometry()[i].Y0();
            }
            
            CalculateEtOperator( Et_Operator, RealX, RealY );
            CalculateDtEtOperator( DtEt_Operator, RealX, RealY );
            CalculateEbartOperator( Ebart_Operator, RealX, RealY );
            CalculateExiOperator( Exi_Operator, RealX, RealY );
            
            noalias( LHSPart1 ) += weight * DetJ0[PointNumber] * prod( trans( Ebart_Operator ), Et_Operator );
            noalias( LHSPart2 ) += weight * DetJ0[PointNumber] * prod( trans( Exi_Operator ), DtEt_Operator );
        }
    }
    
    void CalculateRHS( const WeakPointerVector<Element>::iterator& ie,
                        Matrix& Et_Operator,
                        Matrix& DtEt_Operator,
                        Matrix& Ebart_Operator,
                        Matrix& Exi_Operator,
                        Vector& RHSPart1,
                        Vector& RHSPart2,
                        const ProcessInfo& rCurrentProcessInfo )
    {
        GeometryData::IntegrationMethod ThisIntegrationMethod = ie->GetIntegrationMethod();
            
        const Element::GeometryType::IntegrationPointsArrayType& integration_points = ie->GetGeometry().IntegrationPoints( ThisIntegrationMethod );
            
        const Matrix& Ncontainer = ie->GetGeometry().ShapeFunctionsValues( ThisIntegrationMethod );

        //calculate the Jacobian first
        Element::GeometryType::JacobiansType J0( integration_points.size() );
		J0 = ie->GetGeometry().Jacobian( J0, ThisIntegrationMethod );

        Vector DetJ0;
		DetJ0.resize( integration_points.size(), false );
		noalias( DetJ0 ) = ZeroVector( integration_points.size() );
            
        int dim = ie->GetGeometry().WorkingSpaceDimension();
        Matrix dummy(dim, dim);
		for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
		    MathUtils<double>::InvertMatrix( J0[PointNumber], dummy, DetJ0[PointNumber] );

        //get the stresses
        std::vector<Vector> Stresses;
        Stresses.resize( integration_points.size() );
        ie->GetValueOnIntegrationPoints( STRESSES, Stresses, rCurrentProcessInfo );
        
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
		{
            double weight = integration_points[PointNumber].Weight();
            
            Vector Nvector = row( Ncontainer, PointNumber );
            
            double RealX = 0.0;
            double RealY = 0.0;
            
            for(unsigned int i = 0; i < ie->GetGeometry().size(); i++)
            {
                RealX += Nvector( i ) * ie->GetGeometry()[i].X0();
                RealY += Nvector( i ) * ie->GetGeometry()[i].Y0();
            }
            
            CalculateEbartOperator( Ebart_Operator, RealX, RealY );
            CalculateExiOperator( Exi_Operator, RealX, RealY );
		    
            noalias( RHSPart1 ) += weight * DetJ0[PointNumber] * prod( trans( Ebart_Operator ), Stresses[PointNumber] );
            
            Vector BodyForce = ie->GetProperties()[BODY_FORCE];
            BodyForce.resize( dim, true );
            
            //check for gravity
            Vector gravity( dim );

        		double density = 0.0;
            if( ie->GetValue( USE_DISTRIBUTED_PROPERTIES ) )
            {
            	noalias( gravity ) = ie->GetValue(GRAVITY);
            	density = ie->GetValue(DENSITY);
            }
            else
            {
            	noalias( gravity ) = ie->GetProperties()[GRAVITY];
            	density = ie->GetProperties()[DENSITY];
            }
            
            noalias( BodyForce ) += density * gravity;
            
            noalias( RHSPart2 ) += - weight * DetJ0[PointNumber] * prod( trans( Exi_Operator ), BodyForce );
        }
    }
    
    void CalculateEtOperator( Matrix& Et_Operator, double x, double y )
    {
        noalias( Et_Operator ) = ZeroMatrix( 3, 18 );
        
        Et_Operator( 0, 0 ) = 1.0;
        Et_Operator( 0, 1 ) = x;
        Et_Operator( 0, 2 ) = y;
        Et_Operator( 0, 3 ) = x * y;
        Et_Operator( 0, 4 ) = pow( x, 2 );
        Et_Operator( 0, 5 ) = pow( y, 2 );
        
        Et_Operator( 1, 6 ) = 1.0;
        Et_Operator( 1, 7 ) = x;
        Et_Operator( 1, 8 ) = y;
        Et_Operator( 1, 9 ) = x * y;
        Et_Operator( 1, 10 ) = pow( x, 2 );
        Et_Operator( 1, 11 ) = pow( y, 2 );
        
        Et_Operator( 2, 12 ) = 1.0;
        Et_Operator( 2, 13 ) = x;
        Et_Operator( 2, 14 ) = y;
        Et_Operator( 2, 15 ) = x * y;
        Et_Operator( 2, 16 ) = pow( x, 2 );
        Et_Operator( 2, 17 ) = pow( y, 2 );
    }
    
    void CalculateDtEtOperator( Matrix& DtEt_Operator, double x, double y )
    {
        noalias( DtEt_Operator ) = ZeroMatrix( 2, 18 );
        
        DtEt_Operator( 0, 1 ) = 1.0;
        DtEt_Operator( 0, 3 ) = y;
        DtEt_Operator( 0, 4 ) = 2 * x;
        
        DtEt_Operator( 1, 8 ) = 1.0;
        DtEt_Operator( 1, 9 ) = x;
        DtEt_Operator( 1, 11 ) = 2 * y;
        
        DtEt_Operator( 0, 14 ) = 1.0;
        DtEt_Operator( 0, 15 ) = x;
        DtEt_Operator( 0, 17 ) = 2 * y;
        DtEt_Operator( 1, 13 ) = 1.0;
        DtEt_Operator( 1, 15 ) = y;
        DtEt_Operator( 1, 16 ) = 2 * x;
    }
    
    void CalculateEbartOperator( Matrix& Ebart_Operator, double x, double y )
    {
        noalias( Ebart_Operator ) = ZeroMatrix( 3, 12 );
        
        Ebart_Operator( 0, 0 ) = 1.0;
        Ebart_Operator( 1, 1 ) = 1.0;
        Ebart_Operator( 2, 2 ) = 1.0;
        
        Ebart_Operator( 0, 3 ) = x;
        Ebart_Operator( 0, 4 ) = y;
        Ebart_Operator( 0, 5 ) = 2 * x * y;
        Ebart_Operator( 2, 3 ) = -y;
        Ebart_Operator( 2, 5 ) = -pow( y, 2 );
        
        Ebart_Operator( 1, 6 ) = x;
        Ebart_Operator( 1, 7 ) = y;
        Ebart_Operator( 1, 8 ) = 2 * x * y;
        Ebart_Operator( 2, 7 ) = -x;
        Ebart_Operator( 2, 8 ) = -pow( x, 2 );
        
        Ebart_Operator( 0, 9 ) = pow( y, 2 );
        Ebart_Operator( 0, 11 ) = pow( x, 2 );
        Ebart_Operator( 1, 10 ) = pow( x, 2 );
        Ebart_Operator( 1, 11 ) = pow( y, 2 );
        Ebart_Operator( 2, 11 ) = -2 * x * y;
    }
    
    void CalculateExiOperator( Matrix& Exi_Operator, double x, double y )
    {
        noalias( Exi_Operator ) = ZeroMatrix( 2, 6 );
        
        Exi_Operator( 0, 0 ) = 1;
        Exi_Operator( 0, 1 ) = x;
        Exi_Operator( 0, 2 ) = y;
        
        Exi_Operator( 1, 3 ) = 1;
        Exi_Operator( 1, 4 ) = x;
        Exi_Operator( 1, 5 ) = y;
    }
    
};

}  // namespace Kratos.

#ifdef UTILITY_DEBUG_LEVEL1
#undef UTILITY_DEBUG_LEVEL1
#endif

#ifdef UTILITY_DEBUG_LEVEL2
#undef UTILITY_DEBUG_LEVEL2
#endif

#ifdef UTILITY_DEBUG_LEVEL3
#undef UTILITY_DEBUG_LEVEL3
#endif

#endif // KRATOS_EMBEDDED_DISCONTINUITIES_INFO_H_INCLUDED  defined 

