
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


#if !defined(KRATOS_NONLOCAL_DAMAGE_UTILITIES )
#define  KRATOS_NONLOCAL_DAMAGE_UTILITIES

// System includes
#include <cmath>

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

class NonlocalDamageUtilities
{

protected:

    struct GaussPoint
    {
        // Default Constructor
        GaussPoint() {}

        // Constructor 1
        GaussPoint (ConstitutiveLaw::Pointer pConstitutiveLawPointer,
                    const array_1d<double,3>& Coords,
                    const double& IntegrationWeight)
        {
            pConstitutiveLaw = pConstitutiveLawPointer;
            noalias(Coordinates) = Coords;
            Weight = IntegrationWeight;
        }

        // Member variables
        ConstitutiveLaw::Pointer pConstitutiveLaw;
        array_1d<double,3> Coordinates;
        double Weight;
        std::vector<GaussPoint*> NeighbourPoints;
    };

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

public:

    KRATOS_CLASS_POINTER_DEFINITION( NonlocalDamageUtilities );

    /// Default Constructor
    NonlocalDamageUtilities() {}
    
    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~NonlocalDamageUtilities()
    {
        for (unsigned int i = 0; i < mGaussPointList.size(); i++)
            delete mGaussPointList[i];
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    virtual void SearchGaussPointsNeighbours (Parameters* pParameters, ModelPart& rModelPart)
    {
        KRATOS_ERROR << "Calling the default SearchGaussPointsNeighbours method" << std::endl;
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateNonlocalEquivalentStrain (Parameters* pParameters, const ProcessInfo& CurrentProcessInfo)
    {
        int NGPoints = static_cast<int>(mGaussPointList.size());
        double CharacteristicLength = (*pParameters)["characteristic_length"].GetDouble();
        
        // Loop through all Gauss Points
        #pragma omp parallel for
        for(int i = 0; i < NGPoints; i++)
        {
            const GaussPoint& ReceiverPoint = *(mGaussPointList[i]);
            double LocalEquivalentStrain;
            LocalEquivalentStrain = ReceiverPoint.pConstitutiveLaw->GetValue(LOCAL_EQUIVALENT_STRAIN,LocalEquivalentStrain);;
            double Numerator = ReceiverPoint.Weight*LocalEquivalentStrain;
            double WeightingFunctionDenominator = ReceiverPoint.Weight;
            
            //Loop through neighbours
            for(unsigned int j = 0; j < ReceiverPoint.NeighbourPoints.size(); j++)
            {
                const GaussPoint& SourcePoint = *(ReceiverPoint.NeighbourPoints[j]);
                double Distance;
                this->ComputeNeighbourDistance(Distance,ReceiverPoint,SourcePoint);
                LocalEquivalentStrain = SourcePoint.pConstitutiveLaw->GetValue(LOCAL_EQUIVALENT_STRAIN,LocalEquivalentStrain);

                Numerator += SourcePoint.Weight*exp(-4.0*Distance*Distance/(CharacteristicLength*CharacteristicLength))*LocalEquivalentStrain;
                WeightingFunctionDenominator += SourcePoint.Weight*exp(-4.0*Distance*Distance/(CharacteristicLength*CharacteristicLength));
            }
            double NonlocalEquivalentStrain = Numerator/WeightingFunctionDenominator;
            ReceiverPoint.pConstitutiveLaw->SetValue(NONLOCAL_EQUIVALENT_STRAIN,NonlocalEquivalentStrain,CurrentProcessInfo);
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
    std::vector<GaussPoint*> mGaussPointList;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual void ComputeNeighbourDistance(
        double& rDistance,
        const GaussPoint& ReceiverPoint,
        const GaussPoint& SourcePoint)
    {
        KRATOS_ERROR << "Calling the default ComputeNeighbourDistance method" << std::endl;
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
}; // Class NonlocalDamageUtilities

} // namespace Kratos.

#endif /* KRATOS_NONLOCAL_DAMAGE_UTILITIES defined */
