//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Veronika Singer
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_base_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    MPMParticleBaseLoadCondition::MPMParticleBaseLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
        : MPMBaseLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    MPMParticleBaseLoadCondition::MPMParticleBaseLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : MPMBaseLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer MPMParticleBaseLoadCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
    {
        return Kratos::make_shared<MPMParticleBaseLoadCondition>(NewId, pGeom, pProperties);
    }

    //************************************************************************************
    //************************************************************************************

    Condition::Pointer MPMParticleBaseLoadCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Kratos::make_shared<MPMParticleBaseLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    MPMParticleBaseLoadCondition::~MPMParticleBaseLoadCondition()
    {
    }

    //************************************************************************************
    //************************************************************************************


    double MPMParticleBaseLoadCondition::GetPointLoadIntegrationWeight()
    {
        return 1.0;
    }

    //*************************COMPUTE CURRENT DISPLACEMENT*******************************
    //************************************************************************************
    /*
    This function convert the computed nodal displacement into matrix of (number_of_nodes, dimension)
    */
    Matrix& MPMParticleBaseLoadCondition::CalculateCurrentDisp(Matrix & rCurrentDisp, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        GeometryType& rGeom = GetGeometry();
        const unsigned int number_of_nodes = rGeom.PointsNumber();
        const unsigned int dimension = rGeom.WorkingSpaceDimension();

        rCurrentDisp = ZeroMatrix(number_of_nodes, dimension);

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            const array_1d<double, 3 > & current_displacement  = rGeom[i].FastGetSolutionStepValue(DISPLACEMENT);

            for ( unsigned int j = 0; j < dimension; j++ )
            {
                rCurrentDisp(i,j) = current_displacement[j];
            }
        }

        return rCurrentDisp;

        KRATOS_CATCH( "" )
    }

    /**
     * Shape function values in given point. This method calculate the shape function
     * vector in given point.
     *
     * @param rPoint point which shape function values have to
     * be calculated in it.
     *
     * @return Vector of double which is shape function vector \f$ N \f$ in given point.
     *
    */
    Vector& MPMParticleBaseLoadCondition::MPMShapeFunctionPointValues( Vector& rResult, const array_1d<double,3>& rPoint )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        array_1d<double,3> rPointLocal = ZeroVector(3);
        rPointLocal = GetGeometry().PointLocalCoordinates(rPointLocal, rPoint);

        if (dimension == 2)
        {
            if (number_of_nodes == 3)
            {
                rResult.resize(3, false);

                rResult[0] = 1 - rPointLocal[0] - rPointLocal[1] ;
                rResult[1] = rPointLocal[0] ;
                rResult[2] = rPointLocal[1];
            }
            else if (number_of_nodes == 4)
            {
                rResult.resize(4, false);

                rResult[0] = 0.25 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) ;
                rResult[1] = 0.25 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) ;
                rResult[2] = 0.25 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) ;
                rResult[3] = 0.25 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) ;
            }
        }
        else if (dimension == 3)
        {
            if (number_of_nodes == 4)
            {
                rResult.resize(4, false);

                rResult[0] =  1.0-(rPointLocal[0]+rPointLocal[1]+rPointLocal[2]) ;
                rResult[1] = rPointLocal[0] ;
                rResult[2] = rPointLocal[1];
                rResult[3] = rPointLocal[2];
            }
            else if (number_of_nodes == 8)
            {
                rResult.resize(8, false);

                // Shape Functions (if the first node of the connettivity is the node at (-1,-1,-1))
                // NOTE: Implemented based on Carlos Felippa's Lecture on AFEM Chapter 11
                rResult[0] = 0.125 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) * (1 - rPointLocal[2]) ;
                rResult[1] = 0.125 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) * (1 - rPointLocal[2]) ;
                rResult[2] = 0.125 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) * (1 - rPointLocal[2]) ;
                rResult[3] = 0.125 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) * (1 - rPointLocal[2]) ;
                rResult[4] = 0.125 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) * (1 + rPointLocal[2]) ;
                rResult[5] = 0.125 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) * (1 + rPointLocal[2]) ;
                rResult[6] = 0.125 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) * (1 + rPointLocal[2]) ;
                rResult[7] = 0.125 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) * (1 + rPointLocal[2]) ;
            }
        }

        return rResult;

        KRATOS_CATCH( "" )
    }


} // Namespace Kratos
