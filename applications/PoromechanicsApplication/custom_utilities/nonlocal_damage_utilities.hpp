//
//   Project Name:        KratosPoromechanicsApplication $
//   Last modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:               July 2016 $
//   Revision:            $Revision:                 1.0 $
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
        ConstitutiveLaw::Pointer pConstitutiveLaw;
        array_1d<double,3> Coordinates;
        double Weight;
    };

    ///------------------------------------------------------------------------------------

    struct NeighbourPoint
    {
        ConstitutiveLaw::Pointer pConstitutiveLaw;
        double Weight, Distance;
    };

    ///------------------------------------------------------------------------------------

    struct NonlocalPoint
    {
        ConstitutiveLaw::Pointer pConstitutiveLaw;
        std::vector<NeighbourPoint> NeighbourPoints;
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
        //mNonlocalPointList.clear();
        //mNonlocalPointList.shrink_to_fit();
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    virtual void SearchGaussPointsNeighbours (Parameters* pParameters, ModelPart& rModelPart)
    {
        KRATOS_THROW_ERROR( std::logic_error, "Calling the default SearchGaussPointsNeighbours method", "" )
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateNonlocalEquivalentStrain (Parameters* pParameters, const ProcessInfo& CurrentProcessInfo)
    {
        int NGPoints = static_cast<int>(mNonlocalPointList.size());
        double CharacteristicLength = (*pParameters)["characteristic_length"].GetDouble();
        
        // Loop through all Gauss Points
        #pragma omp parallel for
        for(int i = 0; i < NGPoints; i++)
        {
            double LocalEquivalentStrain;
            double NonlocalEquivalentStrain;
            double Numerator = 0.0;
            double WeightingFunctionDenominator = 0.0;
            
            //Loop through neighbours
            for(unsigned int j = 0; j < mNonlocalPointList[i].NeighbourPoints.size(); j++)
            {
                const NeighbourPoint& MyNeighbourPoint = mNonlocalPointList[i].NeighbourPoints[j];
                const double& Distance = MyNeighbourPoint.Distance;
                LocalEquivalentStrain = MyNeighbourPoint.pConstitutiveLaw->GetValue(LOCAL_EQUIVALENT_STRAIN,LocalEquivalentStrain);
                
                Numerator += MyNeighbourPoint.Weight*exp(-4.0*Distance*Distance/(CharacteristicLength*CharacteristicLength))*LocalEquivalentStrain;
                WeightingFunctionDenominator += MyNeighbourPoint.Weight*exp(-4.0*Distance*Distance/(CharacteristicLength*CharacteristicLength));
            }
            NonlocalEquivalentStrain = Numerator/WeightingFunctionDenominator;
            mNonlocalPointList[i].pConstitutiveLaw->SetValue(NONLOCAL_EQUIVALENT_STRAIN,NonlocalEquivalentStrain,CurrentProcessInfo);
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
    std::vector<NonlocalPoint> mNonlocalPointList;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
}; // Class NonlocalDamageUtilities

} // namespace Kratos.

#endif /* KRATOS_NONLOCAL_DAMAGE_UTILITIES defined */
