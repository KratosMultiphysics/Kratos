//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_WAVEGENERATOR_H_INCLUDED )
#define  KRATOS_WAVEGENERATOR_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"


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
  class WaveGenerator
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of WaveGenerator
      KRATOS_CLASS_POINTER_DEFINITION(WaveGenerator);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      WaveGenerator(){}

      /// Destructor.
      virtual ~WaveGenerator(){}
      
      void GenerateWaveXYPlane(ModelPart::NodesContainerType& rNodes, const double d, const double  H, const double  T, const double  Z0, const double t, const double g)
      {
	double C0=g*T/( 6.2831853071795862 );
	double L0=C0*T;
	double teta = (-6.2831853071795862*t / T) - (3.1415926535897931 / 2.0);
	double aux = cosh(6.2831853071795862*d/L0);
	double ux= 0.5 * H * (g * T / L0) * cos(teta) / aux;
	double uy= 0.5 * H * (g * T / L0) * sin(teta) / aux;
	double h=0.5 * H * cos(teta);


        const PointerVector< Node<3> >::iterator it_begin = rNodes.begin();
	array_1d<double,3> temp;
	
	#pragma omp parallel for private(temp)
	for(int i=0; i<static_cast<int>(rNodes.size()); i++)
	{
	    PointerVector< Node<3> >::iterator it = it_begin + i;
	    const double Z = it->Z();
	    temp[0]  = ux * cosh(6.2831853071795862*(Z - Z0 + d)/L0);
	    temp[1]  = 0.0;
	    temp[2]  = uy * sinh(6.2831853071795862*(Z - Z0 + d)/L0);
	    
            const double distance = -h+Z - Z0;

	    it->FastGetSolutionStepValue(DISTANCE,1)= distance;

            if(distance <= 0) // applying velocity only to dense fluid part
                noalias(it->FastGetSolutionStepValue(VELOCITY))=temp;
	}
      }

      void GenerateVolumeWaveXYPlane(ModelPart::NodesContainerType& rNodes, const double d, const double  H, const double  T, const double  Z0, const double X0, const double t, const double g)
      {

        double C0=g*T/( 6.2831853071795862 );
        double L0=C0*T;

        const PointerVector< Node<3> >::iterator it_begin = rNodes.begin();
	array_1d<double,3> temp;
	array_1d<double,3> velocity;

	#pragma omp parallel for private(temp, velocity)
	for(int i=0; i<static_cast<int>(rNodes.size()); i++)
	{
 	    PointerVector< Node<3> >::iterator it = it_begin + i;
            const double x = it->X();
	    const double Z = it->Z();
            const double alpha = (x - X0) / L0;
            double beta = 0.8 + 0.199 * alpha * 2.;
            if(alpha < 0.01)
                beta = 0.00;
            double teta = (-6.2831853071795862*t / T) + (3.1415926535897931 * 2.0 * alpha) + (3.1415926535897931 / 2.0);
            double aux = cosh(6.2831853071795862*d/L0);
            double ux= 0.5 * H * (g * T / L0) * cos(teta) / aux;
            double uy= 0.5 * H * (g * T / L0) * sin(teta) / aux;
            double h=0.5 * H * cos(teta);

	    temp[0]  = ux * cosh(6.2831853071795862*(Z - Z0 + d)/L0);
	    temp[1]  = 0.0;
	    temp[2]  = uy * sinh(6.2831853071795862*(Z - Z0 + d)/L0);

            const double node_distance = it->FastGetSolutionStepValue(DISTANCE,1);

            const double distance = -h+Z - Z0;

            if((t / T) > alpha )
                it->FastGetSolutionStepValue(DISTANCE,1)= node_distance * beta + (1-beta) * distance;

//            if(distance <= 0) // correcting the pressure only to dense fluid part
//                it->FastGetSolutionStepValue(PRESSURE)= -distance * 1000;

            noalias(velocity) = it->FastGetSolutionStepValue(VELOCITY) * beta;

            if((t / T) > alpha )
            if(distance <= 0) // applying velocity only to dense fluid part
                noalias(it->FastGetSolutionStepValue(VELOCITY))=velocity + (1 - beta) * temp;
	}
        }

        void GenerateComposedVolumeWaveXYPlane(ModelPart::NodesContainerType& rNodes, const double d, const Vector& HVector, const Vector& TVector, const Vector& PhaseVector, const double Z0, const double X0, const double t, const double Length) {
            unsigned int number_of_waves;
            if (HVector.size() <= TVector.size())
                number_of_waves = HVector.size();
            else
                number_of_waves = TVector.size();

            const double g = 9.81;

            const PointerVector< Node < 3 > >::iterator it_begin = rNodes.begin();
            array_1d<double, 3 > wave_velocity;
            array_1d<double, 3 > velocity;

#pragma omp parallel for private(wave_velocity,velocity)
            for (int i = 0; i<static_cast<int> (rNodes.size()); i++) {
                PointerVector< Node < 3 > >::iterator it = it_begin + i;
                const double x = it->X();
                const double Z = it->Z();
                
                wave_velocity = ZeroVector(3);

                double wave_h = 0.00;
                double beta =  (x - X0) / Length;
                double gamma = 0.8 + 0.199 * beta * 2.0;
                if (gamma > 1.00)
                    gamma = 1.00;

                if (beta < 0.1)
                    gamma = 8.00 * beta;


                bool wave_arrived = false;

                for (unsigned int i_wave = 0; i_wave < number_of_waves; i_wave++) {
                    const double T = TVector[i_wave];
                    const double H = HVector[i_wave];

                    double C0 = g * T / (6.2831853071795862);
                    double L0 = C0*T;

                    const double alpha = (x - X0) / L0;
                    
                    if((t / T) > alpha)
                        wave_arrived = true;

                    double teta = (-6.2831853071795862 * t / T) + (3.1415926535897931 * 2.0 * alpha) + (3.1415926535897931 / 2.0) + PhaseVector[i_wave];
                    double aux = cosh(6.2831853071795862 * d / L0);
                    double ux = 0.5 * H * (g * T / L0) * cos(teta) / aux;
                    double uy = 0.5 * H * (g * T / L0) * sin(teta) / aux;
                    double h = 0.5 * H * cos(teta);

                    if ((t / T) > alpha)
                    {
                        wave_velocity[0] +=  ux * cosh(6.2831853071795862 * (Z - Z0 + d) / L0);
                        wave_velocity[1]  = 0.0;
                        wave_velocity[2] +=  uy * sinh(6.2831853071795862 * (Z - Z0 + d) / L0);
                    }
                     if ((t / T) > alpha)
                         wave_h += h;
                }

                const double node_distance = it->FastGetSolutionStepValue(DISTANCE, 1);

                double distance = -wave_h + Z - Z0;

                if(wave_arrived)
                    it->FastGetSolutionStepValue(DISTANCE, 1) = node_distance * gamma + (1 - gamma) * distance;

                //            if(distance <= 0) // correcting the pressure only to dense fluid part
                //                it->FastGetSolutionStepValue(PRESSURE)= -distance * 1000;

                  noalias(velocity) = it->FastGetSolutionStepValue(VELOCITY) * gamma;

                if(wave_arrived)
                       if (distance <= 0) // applying velocity only to dense fluid part
                            noalias(it->FastGetSolutionStepValue(VELOCITY)) = velocity + (1 - gamma) * wave_velocity;
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
      virtual std::string Info() const
      {
	std::stringstream buffer;
        buffer << "WaveGenerator" ;
        return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "WaveGenerator";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const {}
      
            
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
        
        
      ///@} 
      ///@name Protected Operators
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operations
      ///@{ 
        
        
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
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
        
      ///@} 
      ///@name Private  Access 
      ///@{ 
        
        
      ///@}    
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      WaveGenerator& operator=(WaveGenerator const& rOther){return *this;}

      /// Copy constructor.
      WaveGenerator(WaveGenerator const& rOther){};

        
      ///@}    
        
    }; // Class WaveGenerator 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    WaveGenerator& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const WaveGenerator& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_WAVEGENERATOR_H_INCLUDED  defined 


