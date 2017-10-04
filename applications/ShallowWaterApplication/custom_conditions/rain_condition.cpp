//
//   Project Name:        Kratos
//   Last modified by:    Miguel Masó Sotomayor
//   Date:                July 31st 2017
//   Revision:            1.1
//
//

// Project includes
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "shallow_water_application.h"
#include "custom_conditions/rain_condition.hpp"

namespace Kratos
{

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void RainCondition<TNumNodes>::CalculateConsistentMassMatrix(boost::numeric::ublas::bounded_matrix<double,TNumNodes,TNumNodes>& rMassMatrix) 
    {
        KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
    }

    template< unsigned int TNumNodes >
    void RainCondition<TNumNodes>::CalculateLumpedMassMatrix(boost::numeric::ublas::bounded_matrix<double,TNumNodes,TNumNodes>& rMassMatrix) 
    {
        GeometryType& rGeom = GetGeometry();
        const unsigned int number_of_nodes = rGeom.PointsNumber();
        const unsigned int condition_size = number_of_nodes;
        rMassMatrix  = IdentityMatrix(condition_size, condition_size);
        rMassMatrix /= number_of_nodes;
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void RainCondition<TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& rGeom = GetGeometry();
        const unsigned int number_of_nodes = rGeom.PointsNumber();
        const unsigned int condition_size = number_of_nodes;
    
        //Resetting the LHS
        if (rLeftHandSideMatrix.size1() != condition_size)
            rLeftHandSideMatrix.resize(condition_size, condition_size, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(condition_size, condition_size);
    
        //Resetting the RHS
        if (rRightHandSideVector.size() != condition_size)
            rRightHandSideVector.resize(condition_size, false);
        noalias(rRightHandSideVector) = ZeroVector(condition_size);

        // Getting water height unit converter
        mHeightUnitConvert = rCurrentProcessInfo[WATER_HEIGHT_UNIT_CONVERTER];

        CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
        
        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void RainCondition<TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& rGeom = GetGeometry();
        const unsigned int number_of_nodes = rGeom.PointsNumber();
        const unsigned int condition_size = number_of_nodes;
        if(rRightHandSideVector.size() != condition_size)
            rRightHandSideVector.resize(condition_size,false);

        // Initialize variables
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,TNumNodes> mass_matrix = ZeroMatrix(condition_size,condition_size);
        boost::numeric::ublas::bounded_matrix<double,TNumNodes,2> DN_DX = ZeroMatrix(condition_size,2);
        array_1d<double,TNumNodes> N;                                   // Dimension = number of nodes. Position of the gauss point
        array_1d<double,TNumNodes> v_rain;                              // Nodal rain vector

        // Getting data for the given geometry
        double area;
        area = rGeom.Area();

        // Reading properties and conditions
        int counter = 0;
        for(unsigned int iii = 0; iii<TNumNodes; iii++){
            v_rain[counter++] = rGeom[iii].FastGetSolutionStepValue(RAIN);
        }
        
        // Compute parameters and derivatives matrices
        //~ CalculateConsistentMassMatrix(mass_matrix);
        CalculateLumpedMassMatrix(mass_matrix);
        mass_matrix *= mHeightUnitConvert * mHeightUnitConvert;
        // LHS = M*rain
        noalias(rRightHandSideVector) = prod(mass_matrix, v_rain);          // Add <q,rain>         to RHS (Mass Eq.)

        rRightHandSideVector *= area;

        KRATOS_CATCH("")
	}

//----------------------------------------------------------------------

    // This subroutine calculates the nodal contributions for the explicit steps of the
    // Fractional step procedure
    template< unsigned int TNumNodes >
    void RainCondition<TNumNodes>::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void RainCondition<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        GeometryType& rGeom = GetGeometry();
        const unsigned int number_of_nodes = rGeom.PointsNumber();
        const unsigned int condition_size = number_of_nodes;
        if(rResult.size() != condition_size)
            rResult.resize(condition_size,false);

        int counter=0;
        for (unsigned int i = 0; i<condition_size; i++){
            rResult[counter++] = rGeom[i].GetDof(HEIGHT).EquationId();
        }
        
        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------

    template< unsigned int TNumNodes >
    void RainCondition<TNumNodes>::GetDofList(DofsVectorType& rConditionDofList,ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& rGeom = GetGeometry();
        const unsigned int number_of_nodes = rGeom.PointsNumber();
        const unsigned int condition_size = number_of_nodes;
        if(rConditionDofList.size() != condition_size)
            rConditionDofList.resize(condition_size);
        
        int counter=0;
        for (unsigned int i = 0; i<condition_size; i++){
            rConditionDofList[counter++] = rGeom[i].pGetDof(HEIGHT);
        }
        
        KRATOS_CATCH("")
    }


template class RainCondition<3>;
template class RainCondition<4>;

} // namespace Kratos
