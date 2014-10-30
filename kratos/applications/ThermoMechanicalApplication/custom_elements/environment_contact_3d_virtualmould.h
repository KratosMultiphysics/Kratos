/*
==============================================================================
KratosPFEMApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


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

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-08-28 08:42:05 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_ENVIRONMENT_CONTACT_3D_VIRTUALMOULD_CONDITION_H_INCLUDED )
#define  KRATOS_ENVIRONMENT_CONTACT_3D_VIRTUALMOULD_CONDITION_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "environment_contact_3d.h"
#include "includes/convection_diffusion_settings.h"


#include "includes/serializer.h"

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

/// Short class definition.
/** Detail class definition.
*/
class EnvironmentContact3DVirtualMould
    : public EnvironmentContact3D
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of Monolithic2DNeumann
    KRATOS_CLASS_POINTER_DEFINITION(EnvironmentContact3DVirtualMould);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    EnvironmentContact3DVirtualMould():EnvironmentContact3D() {};
    EnvironmentContact3DVirtualMould(IndexType NewId, GeometryType::Pointer pGeometry)
        :EnvironmentContact3D(NewId, pGeometry) {} ;
    EnvironmentContact3DVirtualMould(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
        :EnvironmentContact3D(NewId, pGeometry, pProperties) {};

        
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
    {
    return Condition::Pointer(new EnvironmentContact3DVirtualMould(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }
       
//     virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
//     {
//         this->CalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
//     }

    /// Destructor.
    virtual ~EnvironmentContact3DVirtualMould(){};

    virtual void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
    {
        if(mgauss_temperatures.size1() == 0) //this means it is the first step
        {
            //initialize the database of local variables
            mgauss_temperatures.resize(3,3,false);

        
            //set each of the gauss temperatures to the value of the corresponding node
            for(unsigned int i=0; i<3; i++)
            {
                const double T = rCurrentProcessInfo[AMBIENT_TEMPERATURE];
                for(unsigned int j=0; j<mgauss_temperatures.size2(); j++)
                mgauss_temperatures(i,j) = T;
            }
        }
    }

    virtual void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
    {
        //save virtual 
        for(unsigned int i=0; i<3; i++)
        {
            double HTC_Alpha, heat_flux;
            ComputeGaussHeatFluxAndHTC(i, HTC_Alpha, heat_flux, rCurrentProcessInfo,true);
        }
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
//      virtual String Info() const;

    /// Print information about this object.
//      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{
    Matrix mgauss_temperatures;


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{
    virtual void ComputeGaussHeatFluxAndHTC( unsigned int igauss, 
                                             double& HTC_Alpha, 
                                             double& heat_flux, 
                                             ProcessInfo& rCurrentProcessInfo,
                                             bool save_internal_variables = false)
    {
        const double dist = GetGeometry()[igauss].FastGetSolutionStepValue(DISTANCE);
		double H_BAR=0;
		double B_BAR=0;

        if(dist > 0) //air
        {
            heat_flux = 0.0;
            HTC_Alpha = 0.0;
        }
        else
        {
            ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
            const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
            const Variable<double>& rTransferCoefVar = my_settings->GetTransferCoefficientVariable();
            
            
            //here we decide if the mould is made of sand or of steel
            const double htc = GetGeometry()[igauss].FastGetSolutionStepValue(rTransferCoefVar);
	
			
		    double DTIME= rCurrentProcessInfo[DELTA_TIME];
            
            
			// Now depending on HTC we select the material to be sand,steel or copper

			unsigned int  MOULD_MATERIAL; // 1=SAND, 2=STEEL, 3=COPPER
			MOULD_MATERIAL=2; // Default to Steel
			if(htc<750.0){MOULD_MATERIAL=1;} // We guess the material depending on the htc
			if(htc>2700.0){MOULD_MATERIAL=3;}// To be improved
			
			const double VDENS= rCurrentProcessInfo[MOULD_DENSITY]; 
            const double VSHEA= rCurrentProcessInfo[MOULD_SPECIFIC_HEAT]; 
            const double VTHIK= rCurrentProcessInfo[MOULD_THICKNESS]; 
            const double VFACT= rCurrentProcessInfo[MOULD_VFACT];
            const double SFACT= rCurrentProcessInfo[MOULD_SFACT];
            const double HTCCD= htc; 
            const double HTCCV= rCurrentProcessInfo[MOULD_HTC_ENVIRONMENT]; 
            const double VCOND = rCurrentProcessInfo[MOULD_CONDUCTIVITY]; 

           
			//switch(MOULD_MATERIAL)
			//{
			//	case 1: // Sand
			//		VDENS=1500; 
			//		VSHEA= 1000; 
			//		VTHIK= 0.05; 
			//		VFACT= 1.0; //*
			//		SFACT= 1.0; //*
			//		HTCCD= 300; //*
			//		HTCCV= 30; //*
			//		VCOND = 1.0; //*
			//		break;
			//	case 2: // Steel
			//		VDENS=7800; 
			//		VSHEA= 500; 
			//		VTHIK= 0.05; 
			//		VFACT= 1.0; //*
			//		SFACT= 1.0; //*
			//		HTCCD= 1000.0; //*
			//		HTCCV= 100.0; //*
			//		VCOND = 30.0; //*
			//		break;
			//	case 3: // Copper
			//		VDENS=8960; 
			//		VSHEA= 400; 
			//		VTHIK= 0.05; 
			//		VFACT= 1.0; //*
			//		SFACT= 1.0; //*
			//		HTCCD= 1500; //*
			//		HTCCV= 30.0; //*
			//		VCOND = 400.0; //*
			//		break;
			//}

			//
			//////default to sand mould
   ////         double VDENS=1500; 
   ////         double VSHEA= 1000; 
   ////         double VTHIK= 0.1; 
   ////         double VFACT= 1.0;
   ////         double SFACT= 1.0;
   ////         double HTCCD= 300; 
   ////         double HTCCV= 30; 
   ////         double VCOND = 1.0; 
   ////        
   ////         if(htc > 1100.0) //steel mould
   ////         {
   ////             VDENS=7800; 
   ////             VSHEA= 500; 
   ////             VTHIK= 0.1; 
   ////             HTCCD= 3000.0; 
   ////             HTCCV= 30.0; 
   ////             VCOND = 30.0;
   ////         }
			const unsigned int virtual_mould_type=2;

			// Linear virtual mould
			if(virtual_mould_type==1){
				double C  = (0.5*VDENS*VSHEA*VTHIK/DTIME)*VFACT;
				double K  = (VCOND/VTHIK)*VFACT;
				double H  = HTCCD;
				double HE = HTCCV*SFACT;
 
				//obtain the temperatures
				const double TEM2N = mgauss_temperatures(igauss,0);
				const double TEM3N = mgauss_temperatures(igauss,1);
				const double TGAUS = GetGeometry()[igauss].FastGetSolutionStepValue(rUnknownVar);
				const double ENVTE = rCurrentProcessInfo[AMBIENT_TEMPERATURE];

				const double DENON = (HE+H+2.0*C)*K+(H+C)*HE+C*H+C*C;

				H_BAR = (H*HE+2.0*C*H)*K+C*H*HE+C*C*H;
				H_BAR = H_BAR/DENON;

				B_BAR = H*HE*K*ENVTE+C*H*K*TEM3N+(C*H*K+C*H*HE+C*C*H)*TEM2N;
				B_BAR = B_BAR/DENON;

		//Compute HEAT FLUX
				heat_flux = -(H_BAR*TGAUS-B_BAR); //*DVOLU;    // 'VIRTUAL-MOULD convection heat flux'

				if(save_internal_variables == true)
				{
					double TEMP2 = (H*K+H*HE+C*H)*TGAUS+HE*K*ENVTE+C*K*TEM3N+(C*K+C*HE+C*C)*TEM2N;
					TEMP2 = TEMP2/DENON;
            
					double TEMP3 = H*K*TGAUS+(HE*K+(H+C)*HE)*ENVTE+(C*K+C*H+C*C)*TEM3N+C*K*TEM2N;
					TEMP3 = TEMP3/DENON;
            

					mgauss_temperatures(igauss,0) = TEMP2;
					mgauss_temperatures(igauss,1) = TEMP3;
				}
			}
			// Quadratic virtual mould
			if(virtual_mould_type==2){
				double C  = (VDENS*VSHEA*VTHIK/DTIME)*VFACT;
				double K  = (VCOND/VTHIK)*VFACT;
				double H  = HTCCD;
				double HE = HTCCV*SFACT;
 
				//obtain the temperatures
				const double TEM2N = mgauss_temperatures(igauss,0);
				const double TEM3N = mgauss_temperatures(igauss,1);
				const double TEM4N = mgauss_temperatures(igauss,2); //Ojo revisar que sea esto
				const double TGAUS = GetGeometry()[igauss].FastGetSolutionStepValue(rUnknownVar);
				const double ENVTE = rCurrentProcessInfo[AMBIENT_TEMPERATURE];

				const double DENON = 288.0*(HE+H+C)*K*K+((288.0*H+132.0*C)*HE+132.0*C*H+36.0*C*C)*K+(36.0*C*H+6.0*C*C)*HE+6.0*C*C*H+C*C*C;
				
				H_BAR = 288.0*(H*HE+C*H)*K*K+(132.0*C*H*HE+36.0*C*C*H)*K+6.0*C*C*H*HE+C*C*C*H;
				H_BAR = H_BAR/DENON;

				B_BAR = (288.0*H*HE*K*K-12.0*C*H*HE*K)*ENVTE+(48.0*C*H*K*K-2.0*C*C*H*K)*TEM4N+(192.0*C*H*K*K+(96.0*C*H*HE+16.0*C*C*H)*K)*TEM3N+(48.0*C*H*K*K+(48.0*C*H*HE+22.0*C*C*H)*K+6.0*C*C*H*HE+C*C*C*H)*TEM2N;
				B_BAR = B_BAR/DENON;

			//Compute HEAT FLUX
				heat_flux = -(H_BAR*TGAUS-B_BAR); //*DVOLU;    // 'VIRTUAL-MOULD convection heat flux'

				if(save_internal_variables == true)
				{
					double TEMP2 = (288.0*H*K*K+(288.0*H*HE+132.0*C*H)*K+36.0*C*H*HE+6.0*C*C*H)*TGAUS+(288.0*HE*K*K-12.0*C*HE*K)*ENVTE+(48.0*C*K*K-2.0*C*C*K)*TEM4N+(192.0*C*K*K+(96.0*C*HE+16.0*C*C)*K)*TEM3N+(48.0*C*K*K+(48*C*HE+22.0*C*C)*K+6.0*C*C*HE+C*C*C)*TEM2N;
					TEMP2 = TEMP2/DENON;
            
					double TEMP3 =(288.0*H*K*K+(144.0*H*HE+24.0*C*H)*K)*TGAUS+(288.0*HE*K*K+(144.0*H+24.0*C)*HE*K)*ENVTE+(48.0*C*K*K+(24.0*C*H+4.0*C*C)*K)*TEM4N+(192.0*C*K*K+(84.0*C*HE+84.0*C*H+28.0*C*C)*K+(36.0*C*H+6.0*C*C)*HE+6.0*C*C*H+C*C*C)*TEM3N+(48.0*C*K*K+(24.0*C*HE+4.0*C*C)*K)*TEM2N;
					TEMP3 = TEMP3/DENON;
            
					double TEMP4= (288.0*H*K*K-12.0*C*H*K)*TGAUS+(288.0*HE*K*K+(288.0*H+132.0*C)*HE*K+(36.0*C*H+6.0*C*C)*HE)*ENVTE+(48.0*C*K*K+(48.0*C*H+22.0*C*C)*K+6.0*C*C*H+C*C*C)*TEM4N+(192.0*C*K*K+(96.0*C*H+16.0*C*C)*K)*TEM3N+(48.0*C*K*K-2.0*C*C*K)*TEM2N;

					TEMP4 = TEMP4/DENON;

					mgauss_temperatures(igauss,0) = TEMP2;
					mgauss_temperatures(igauss,1) = TEMP3;
					mgauss_temperatures(igauss,2) = TEMP4;

				}
			}
        
    //... Compute linearization of the HEAT FLUX of the HEAT FLUX contribution

            HTC_Alpha = H_BAR; //*DVOLU;          // 'VIRTUAL-MOULD convection heat flux'
        }
        
  
// KRATOS_WATCH(heat_flux)
// KRATOS_WATCH(HTC_Alpha)

    }


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, EnvironmentContact3D);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, EnvironmentContact3D);
    }


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Privatmaterialse Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

// 	double ERFC( const double etta);
// 
// 	double CalcTempSemiInfiniteWall(ProcessInfo& rCurrentProcessInfo, const double T_0);
    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //Monolithic2DNeumann& operator=(const Monolithic2DNeumann& rOther);

    /// Copy constructor.
    //Monolithic2DNeumann(const Monolithic2DNeumann& rOther);


    ///@}

}; // Class Monolithic2DNeumann

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    Monolithic2DNeumann& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const Monolithic2DNeumann& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_ENVIRONMENT_CONTACT_3D_VIRTUALMOULD_CONDITION_H_INCLUDED   defined 


