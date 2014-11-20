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

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hurga $
//   Date:                $Date: 2007-10-18 16:23:41 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_PARTICLE_FENG_H_INCLUDED )
#define      KRATOS_PARTICLE_FENG_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "spheric_particle.h"

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

    /// Total Lagrangian element for 2D and 3D geometries.

    /**
     * Implements a total Lagrangian definition for structural analysis.
     * This works for arbitrary geometries in 2D and 3D
     */

    class DEM_FEM_Particle
    : public SphericParticle
    {
    public:
        ///@name Type Definitions
        ///@{
        ///Reference type definition for constitutive laws


        /// Counted pointer of Particle
        KRATOS_CLASS_POINTER_DEFINITION(DEM_FEM_Particle);


        typedef WeakPointerVector<Element > ParticleWeakVector;
        typedef WeakPointerVector<Element >::iterator ParticleWeakIterator;


        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        DEM_FEM_Particle(IndexType NewId, GeometryType::Pointer pGeometry);
        DEM_FEM_Particle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        /// Destructor.
        virtual ~DEM_FEM_Particle();


        ///@}
        ///@name Operators
        ///@{
        ///@}
        ///@name Operations
        ///@{
        /**
         * Returns the currently selected integration method
         * @return current integration method selected
         */

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

        void Initialize();

        void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

        void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);

        void ComputeBallToBallContactForce(   array_1d<double, 3>& rContactForce, array_1d<double, 3>& rContactMoment, array_1d<double, 3>& rElasticForce, 
               array_1d<double, 3>& InitialRotaMoment, double dt, ProcessInfo& rCurrentProcessInfo); 
			   
		void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, double dt, const array_1d<double,3>& gravity);
		
		
		void ComputeBallToRigidFaceContactForce(   array_1d<double, 3>& rContactForce, array_1d<double, 3>& rContactMoment, array_1d<double, 3>& rElasticForce, array_1d<double, 3>& InitialRotaMoment, ProcessInfo& rCurrentProcessInfo, double dt);

        ///@name Protected static Member Variables
 
	  void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo);
      void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo);
	  
	  
	  void SetInitialBallNeighbor();
	  void SetInitialRigidFaceNeighbor();

   protected:


    private:
     
		double mfcohesion;
		double mftension;
		
		vector<int> maInitialBallNeighborID;
		vector<int> maInitialRigidFaceNeighborID;
		
		vector<int> maInitialBallNeighborFailureType;
		vector<int> maInitialRigidFaceNeighborFailureType;
		
		//vector< array_1d<double, 3> > maPPContactForce;
		//vector< array_1d<double, 3> > maPBContactForce;
		
        ///@}
        ///@name Private Operations
        ///@{

	    
	    ///@}
            ///@name Private  Access
            ///@{
            ///@}
	    
	    ///@}
	    ///@name Serialization
	    ///@{	
	    friend class Serializer;
	    
	    // A private default constructor necessary for serialization 

        DEM_FEM_Particle() : SphericParticle()
        {

        }
   

        virtual void save(Serializer& rSerializer) const;
//        {
//            rSerializer.save("Name", "Particle");
//            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
//        }

        virtual void load(Serializer& rSerializer);
//        {
//            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
//        }



        ///@name Private Inquiry
        ///@{
        ///@}
        ///@name Un accessible methods
        ///@{
        /// Assignment operator.
        //Particle& operator=(const Particle& rOther);
        /// Copy constructor.
        //Particle(const Particle& rOther);
        ///@}

    }; // Class Particle

    ///@}
    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Input and output
    ///@{
    /// input stream function
    /*  inline std::istream& operator >> (std::istream& rIStream,
    Particle& rThis);
     */
    /// output stream function
    /*  inline std::ostream& operator << (std::ostream& rOStream,
    const Particle& rThis)
    {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
    }*/
    ///@}

} // namespace Kratos.
#endif // KRATOS_TOTAL_LAGRANGIAN_ELEMENT_H_INCLUDED  defined 
