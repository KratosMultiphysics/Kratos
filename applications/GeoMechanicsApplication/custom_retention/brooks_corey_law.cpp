// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Samira Fazli
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_retention/brooks_corey_law.h"

namespace Kratos
{
//-------------------------------------------------------------------------------------------------
BrooksCoreyLaw::BrooksCoreyLaw()
    : RetentionLaw()
{
}

//-------------------------------------------------------------------------------------------------
BrooksCoreyLaw::BrooksCoreyLaw(const BrooksCoreyLaw& rOther)
    : RetentionLaw(rOther)
{
}

//-------------------------------------------------------------------------------------------------
RetentionLaw::Pointer BrooksCoreyLaw::Clone() const
{
    return Kratos::make_shared<BrooksCoreyLaw>(*this);
}

//-------------------------------------------------------------------------------------------------
BrooksCoreyLaw::~BrooksCoreyLaw()
{
}

//-------------------------------------------------------------------------------------------------
double BrooksCoreyLaw::
    CalculateSaturation(Parameters &rParameters)
{
    KRATOS_TRY;

    const double &p = rParameters.GetFluidPressure();
    const Properties &rMaterialProperties = rParameters.GetMaterialProperties();
    const double &pb     = rMaterialProperties[AIR_ENTRY_PRESSURE];
    const double &pe     = rMaterialProperties[AIR_EXPULSION_PRESSURE];
    if (p > 0.0 && p > pb )
      {
        if(p>=Lastp)
       {
        const double &satMax = rMaterialProperties[SATURATED_SATURATION];
        const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
        const double &pb     = rMaterialProperties[AIR_ENTRY_PRESSURE];
        const double &Lambda     = rMaterialProperties[BROOKS_COREY_PORE_SIZE_INDEX];
        

        double sat = satMin + (satMax - satMin) *  pow(pb/p, Lambda);
        return sat;
        }
        else if (p>pe)

        {
        const double &satMax = rMaterialProperties[SATURATED_SATURATION];
        const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
        const double &pe     = rMaterialProperties[AIR_EXPULSION_PRESSURE];
        const double &Lambdawet     = rMaterialProperties[BROOKS_COREY_PORE_SIZE_INDEX_WET];
        

        double sat = satMin + (satMax - satMin) *  pow(pe/p, Lambdawet);
        return sat;
        
        }
        else
        {
            double sat = rMaterialProperties[SATURATED_SATURATION];
        return sat;
        
        }

        } 
    else if(p > 0.0 && p > pe && p>=Lastp) 
    {   
        
        double sat = rMaterialProperties[SATURATED_SATURATION];
        return sat;
        
    }
    else if(p > 0.0 && p > pe && p<Lastp)
    {
        const double &satMax = rMaterialProperties[SATURATED_SATURATION];
        const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
        const double &pe     = rMaterialProperties[AIR_EXPULSION_PRESSURE];
        const double &Lambdawet     = rMaterialProperties[BROOKS_COREY_PORE_SIZE_INDEX_WET];
        
        double sat = satMin + (satMax - satMin) *  pow(pe/p, Lambdawet);
        return sat;
    }
    
    

    else 
    {
        double sat = rMaterialProperties[SATURATED_SATURATION];
        return sat;
        
    }

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double BrooksCoreyLaw::
    CalculateEffectiveSaturation(Parameters &rParameters)
{
    KRATOS_TRY;

    const double sat = CalculateSaturation(rParameters);

    const auto &rMaterialProperties = rParameters.GetMaterialProperties();
    const double &satMax = rMaterialProperties[SATURATED_SATURATION];
    const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];

    double effectiveSat = (sat - satMin) / (satMax - satMin);

    return effectiveSat;

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double BrooksCoreyLaw::
    CalculateDerivativeOfSaturation(Parameters &rParameters)
{
    KRATOS_TRY;
    const double &p = rParameters.GetFluidPressure();
    const Properties &rMaterialProperties = rParameters.GetMaterialProperties();
    const double &pb     = rMaterialProperties[AIR_ENTRY_PRESSURE];
    const double &pe     = rMaterialProperties[AIR_EXPULSION_PRESSURE];
    

    if (p > 0.0 && p > pb )
      {
        if(p>=Lastp)
       {
       	const auto &rMaterialProperties = rParameters.GetMaterialProperties();
        const double &satMax = rMaterialProperties[SATURATED_SATURATION];
        const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
        const double &pb     = rMaterialProperties[AIR_ENTRY_PRESSURE];
        const double &Lambda     = rMaterialProperties[BROOKS_COREY_PORE_SIZE_INDEX];
        

        double dSdp = (satMax - satMin) * (-Lambda) * pow(pb,Lambda)*pow(p, (-Lambda-1.0));
    
        return dSdp;
        
        }
        else if (p>=pe)

        {
        const auto &rMaterialProperties = rParameters.GetMaterialProperties();
        const double &satMax = rMaterialProperties[SATURATED_SATURATION];
        const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
        const double &pe     = rMaterialProperties[AIR_EXPULSION_PRESSURE];
        const double &Lambdawet     = rMaterialProperties[BROOKS_COREY_PORE_SIZE_INDEX_WET];
        

         double dSdp = (satMax - satMin) * (-Lambdawet) * pow(pe,Lambdawet)*pow(p, (-Lambdawet-1.0));
    
        return dSdp;
       
        }
        else
        {
            return 0.0;
       
        }

        } 
    else if(p > 0.0 && p > pe && p>=Lastp ) 
         {

       	  return 0.0;
          
          }
    else if(p > 0.0 && p > pe && p<Lastp)
         {
    	 const auto &rMaterialProperties = rParameters.GetMaterialProperties();
         const double &satMax = rMaterialProperties[SATURATED_SATURATION];
         const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
         const double &pe     = rMaterialProperties[AIR_EXPULSION_PRESSURE];
         const double &Lambdawet     = rMaterialProperties[BROOKS_COREY_PORE_SIZE_INDEX_WET];

         double dSdp = (satMax - satMin) * (-Lambdawet) * pow(pe,Lambdawet)*pow(p, (-Lambdawet-1.0));
    
         return dSdp;
         
         }
    

    else 
    {
         return 0.0;
        
    }

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double BrooksCoreyLaw::
    CalculateRelativePermeability(Parameters &rParameters)
{
    KRATOS_TRY;

    const double effSat = CalculateEffectiveSaturation(rParameters);
     const double sat = CalculateSaturation(rParameters);

    const auto &rMaterialProperties = rParameters.GetMaterialProperties();
    const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
    const double &Lambda  = rMaterialProperties[BROOKS_COREY_PORE_SIZE_INDEX];
    const double &n1  = rMaterialProperties[RELPERM_POWER];
    //double relPerm = pow(effSat, ((2+3*Lambda)/Lambda)); 
    //double relPerm = pow(effSat, n1);
    double relPerm = pow(sat, n1);
    

    const double &minRelPerm = rMaterialProperties[MINIMUM_RELATIVE_PERMEABILITY];

    relPerm = std::max(relPerm, minRelPerm);

    return relPerm;

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
double BrooksCoreyLaw::
    CalculateBishopCoefficient(Parameters &rParameters)
{
    KRATOS_TRY;
     const double &p = rParameters.GetFluidPressure();
     const auto &rMaterialProperties = rParameters.GetMaterialProperties();
     const double &pb     = rMaterialProperties[AIR_ENTRY_PRESSURE];
     const double &pe     = rMaterialProperties[AIR_EXPULSION_PRESSURE];
     
     
if (p > 0.0 && p > pb )
      {
        if(p>=Lastp)
       {
       	const auto &rMaterialProperties = rParameters.GetMaterialProperties();
        const double &Porosity = rMaterialProperties[POROSITY];
        const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
        const double &pb     = rMaterialProperties[AIR_ENTRY_PRESSURE];
        const double &Lambda  = rMaterialProperties[BROOKS_COREY_PORE_SIZE_INDEX];
        const double &Beta  = rMaterialProperties[BROOKS_COREY_FITTING_PARAMETER];// which is considered between 0.4 to 0.7

        double BishopCo = pow(pb/p, Beta)+pow(pb/p, 1+Beta)*Porosity*(Lambda/(Lambda-1))*(1-satMin)*(1-pow(pb/p, Lambda-1)); 
   
        return BishopCo;
        }
        else if (p>=pe)

        {
        const auto &rMaterialProperties = rParameters.GetMaterialProperties();
        const double &Porosity = rMaterialProperties[POROSITY];
        const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
        const double &pe     = rMaterialProperties[AIR_EXPULSION_PRESSURE];
        const double &Lambdawet  = rMaterialProperties[BROOKS_COREY_PORE_SIZE_INDEX_WET];
        const double &Betawet  = rMaterialProperties[BROOKS_COREY_FITTING_PARAMETER_WET];// which is considered between 0.4 to 0.7

        double BishopCo = pow(pe/p, Betawet)+pow(pe/p, 1+Betawet)*Porosity*(Lambdawet/(Lambdawet-1))*(1-satMin)*(1-pow(pe/p, Lambdawet-1)); 
   
        return BishopCo;
        
        }
        else
        {
           double BishopCo =1;
        return 1.0;
        }

        } 
    else if(p > 0.0 && p > pe && p>=Lastp) 
        {
         
         
          double BishopCo =1;
          return 1.0;
          
         }
    else if(p > 0.0 && p > pe && p<Lastp)
         {
    	 const auto &rMaterialProperties = rParameters.GetMaterialProperties();
         const double &Porosity = rMaterialProperties[POROSITY];
         const double &satMax = rMaterialProperties[SATURATED_SATURATION];
         const double &satMin = rMaterialProperties[RESIDUAL_SATURATION];
         const double &pe     = rMaterialProperties[AIR_EXPULSION_PRESSURE];
         const double &Lambdawet     = rMaterialProperties[BROOKS_COREY_PORE_SIZE_INDEX_WET];
         const double &Betawet  = rMaterialProperties[BROOKS_COREY_FITTING_PARAMETER_WET];// which is considered between 0.4 to 0.7

         double BishopCo = pow(pe/p, Betawet)+pow(pe/p, 1+Betawet)*Porosity*(Lambdawet/(Lambdawet-1))*(1-satMin)*(1-pow(pe/p, Lambdawet-1)); 
   
         return BishopCo;
         
         }
    

     else 
    {
        double BishopCo =1;
        return 1.0;
    }


    KRATOS_CATCH("")
}
//-------------------------------------------------------------------------------------------------
double BrooksCoreyLaw::
    CalculateIncrementOfSuction(Parameters &rParameters)
{
  KRATOS_TRY;
     const double &p = rParameters.GetFluidPressure();
     const auto &rMaterialProperties = rParameters.GetMaterialProperties();
     const double &pb     = rMaterialProperties[AIR_ENTRY_PRESSURE];

    if (p >= 0.0 )
    
    {
        double Suction = std::max(p, 0.0);
        double LastSuction= std::max(Lastp, 0.0);
        double IncSuction=Suction-LastSuction;
        return IncSuction;
       
    }
    else
    {
        return 0.0;
    }
    KRATOS_CATCH("")
}
//-------------------------------------------------------------------------------------------------
double& BrooksCoreyLaw::CalculateValue(RetentionLaw::Parameters& rParameterValues,
                                        const Variable<double>& rThisVariable,
                                        double& rValue)
{
    if (rThisVariable == DEGREE_OF_SATURATION)
    {
        rValue = this->CalculateSaturation(rParameterValues);
        return rValue;
    }
    else if (rThisVariable == EFFECTIVE_SATURATION)
    {
        rValue = this->CalculateEffectiveSaturation(rParameterValues);
        return rValue;
    }
    else if (rThisVariable == BISHOP_COEFFICIENT)
    {
        rValue = this->CalculateBishopCoefficient(rParameterValues);
        return rValue;
    }
    else if (rThisVariable == DERIVATIVE_OF_SATURATION)
    {
        rValue = this->CalculateDerivativeOfSaturation(rParameterValues);
        return rValue;
    }
    else if (rThisVariable == RELATIVE_PERMEABILITY)
    {
        rValue = this->CalculateRelativePermeability(rParameterValues);
        return rValue;
    }
     else if (rThisVariable == INCREMENT_OF_SUCTION)
    {
        rValue = this->CalculateIncrementOfSuction(rParameterValues);
        return rValue;
    }
    return( rValue );
}

//------------------------- RETENSION LAW GENERAL FEATURES ----------------------------------------
//-------------------------------------------------------------------------------------------------
void BrooksCoreyLaw::
    InitializeMaterial(const Properties& rMaterialProperties,
                       const GeometryType& rElementGeometry,
                       const Vector& rShapeFunctionsValues)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void BrooksCoreyLaw::
    Initialize( Parameters &rParameters)
{
   Lastp = rParameters.GetFluidPressure();
    
}

//-------------------------------------------------------------------------------------------------
void BrooksCoreyLaw::
    InitializeSolutionStep(Parameters &rParameters)
{
     // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void BrooksCoreyLaw::
    Finalize(Parameters &rParameters)
{
    // nothing is needed
}

//-------------------------------------------------------------------------------------------------
void BrooksCoreyLaw::
    FinalizeSolutionStep(Parameters &rParameters)
{
    
     Lastp = rParameters.GetFluidPressure();
     //KRATOS_INFO("Lastpfinalize") << Lastp << std::endl;
     

}

//-------------------------------------------------------------------------------------------------
int BrooksCoreyLaw::Check(const Properties& rMaterialProperties,
                           const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(!rMaterialProperties.Has(SATURATED_SATURATION))
                    << "SATURATED_SATURATION is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[SATURATED_SATURATION] < 0.0)
                    << "SATURATED_SATURATION cannot be less than 0 " << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(RESIDUAL_SATURATION))
                    << "RESIDUAL_SATURATION is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[RESIDUAL_SATURATION] < 0.0)
                    << "RESIDUAL_SATURATION cannot be less than 0 " << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(SATURATED_SATURATION))
                    << "SATURATED_SATURATION is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[SATURATED_SATURATION] > 1.0)
                    << "SATURATED_SATURATION cannot be greater than 1.0 " << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(RESIDUAL_SATURATION))
                    << "RESIDUAL_SATURATION is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[RESIDUAL_SATURATION] > 1.0)
                    << "RESIDUAL_SATURATION cannot be greater than 1.0 " << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(SATURATED_SATURATION))
                    << "SATURATED_SATURATION is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[SATURATED_SATURATION] < rMaterialProperties[RESIDUAL_SATURATION])
                    << "RESIDUAL_SATURATION cannot be greater than SATURATED_SATURATION " << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(AIR_ENTRY_PRESSURE))
                    << "AIR_ENTRY_PRESSURE is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(!(rMaterialProperties[AIR_ENTRY_PRESSURE] > 0.0))
                    << "AIR_ENTRY_PRESSURE must be greater than 0 " << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(BROOKS_COREY_PORE_SIZE_INDEX))
                    << "BROOKS_COREY_PORE_SIZE_INDEX is not availabe in material parameters" << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(BROOKS_COREY_FITTING_PARAMETER))
                    << "BROOKS_COREY_FITTING_PARAMETER is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(!(rMaterialProperties[BROOKS_COREY_FITTING_PARAMETER] > 0.4))
                    << "BROOKS_COREY_FITTING_PARAMETER must be greater than 0.4 " << std::endl;
    KRATOS_ERROR_IF(!(rMaterialProperties[BROOKS_COREY_FITTING_PARAMETER] < 0.7))
                    << "BROOKS_COREY_FITTING_PARAMETER must be smaller than 0.7 " << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(AIR_EXPULSION_PRESSURE))
                    << "AIR_EXPULSION_PRESSURE is not availabe in material parameters" << std::endl;

    KRATOS_ERROR_IF(!(rMaterialProperties[AIR_EXPULSION_PRESSURE] > 0.0))
                    << "AIR_EXPULSION_PRESSURE must be greater than 0 " << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(BROOKS_COREY_PORE_SIZE_INDEX_WET))
                    << "BROOKS_COREY_PORE_SIZE_INDEX_WET is not availabe in material parameters" << std::endl;

    KRATOS_ERROR_IF(!rMaterialProperties.Has(BROOKS_COREY_FITTING_PARAMETER_WET))
                    << "BROOKS_COREY_FITTING_PARAMETER_WET is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(!(rMaterialProperties[BROOKS_COREY_FITTING_PARAMETER_WET] > 0.4))
                    << "BROOKS_COREY_FITTING_PARAMETER_WET must be greater than 0.4 " << std::endl;
    KRATOS_ERROR_IF(!(rMaterialProperties[BROOKS_COREY_FITTING_PARAMETER_WET] < 0.7))
                    << "BROOKS_COREY_FITTING_PARAMETER_WET must be smaller than 0.7 " << std::endl;


    KRATOS_ERROR_IF(!rMaterialProperties.Has(MINIMUM_RELATIVE_PERMEABILITY))
                    << "MINIMUM_RELATIVE_PERMEABILITY is not availabe in material parameters" << std::endl;
    KRATOS_ERROR_IF(!(rMaterialProperties[MINIMUM_RELATIVE_PERMEABILITY] > 0.0))
                    << "MINIMUM_RELATIVE_PERMEABILITY must be greater than 0 " << std::endl;

    return 0;
}

} // Namespace Kratos