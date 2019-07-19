//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                 LHauser $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_MOHR_COULOMB_PLASTIC_POTENTIAL_H_INCLUDED )
#define      KRATOS_MOHR_COULOMB_PLASTIC_POTENTIAL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_surfaces/yield_surface.hpp"
#include "custom_utilities/stress_invariants_utilities.hpp"
#include "custom_utilities/shape_deviatoric_plane_utilities.hpp"
//#include "custom_utilities/shape_deviatoric_plane_matsuoka_utilities.hpp"


namespace Kratos
{
   ///@addtogroup ConstitutiveModelsApplication
   ///@{

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
   template<class THardeningRule>
      class MohrCoulombPlasticPotential : public YieldSurface<THardeningRule>
   {
      public:

         ///@name Type Definitions
         ///@{

         typedef ConstitutiveModelData::MatrixType                          MatrixType;
         typedef ConstitutiveModelData::VectorType                          VectorType;
         typedef ConstitutiveModelData::ModelData                        ModelDataType;
         typedef ConstitutiveModelData::MaterialData                  MaterialDataType;

         typedef YieldSurface<THardeningRule>                                 BaseType;
         typedef typename YieldSurface<THardeningRule>::Pointer        BaseTypePointer;
         typedef typename BaseType::PlasticDataType                    PlasticDataType;

         /// Pointer definition of MohrCoulombPlasticPotential
         KRATOS_CLASS_POINTER_DEFINITION( MohrCoulombPlasticPotential );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         MohrCoulombPlasticPotential() : BaseType() {}

         /// Copy constructor.
         MohrCoulombPlasticPotential(MohrCoulombPlasticPotential const& rOther) : BaseType(rOther) {}


         /// Assignment operator.
         MohrCoulombPlasticPotential& operator=(MohrCoulombPlasticPotential const& rOther)
         {
            BaseType::operator=(rOther);
            return *this;
         }

         /// Clone.
         BaseTypePointer Clone() const override
         {
            return Kratos::make_shared<MohrCoulombPlasticPotential>(*this);
         }

         /// Destructor.
         ~MohrCoulombPlasticPotential() override {}


         ///@}
         ///@name Operators
         ///@{


         ///@}
         ///@name Operations
         ///@{


         //*************************************************************************************
         //*************************************************************************************
         // evaluation of the derivative of the yield surface respect the stresses
         VectorType& CalculateDeltaStressYieldCondition(const PlasticDataType& rVariables, VectorType& rDeltaStressYieldCondition) override
         {
            KRATOS_TRY

            const ModelDataType & rModelData = rVariables.GetModelData();
            const Properties& rMaterialProperties = rModelData.GetProperties();

            const MatrixType    & rStressMatrix = rModelData.GetStressMatrix();

            // Material Parameters
            const double & rDilatancy = rMaterialProperties[DILATANCY_ANGLE];

            double MeanStress, LodeAngle, J2;
            StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, J2, LodeAngle);


            VectorType V1, V2, V3;
            StressInvariantsUtilities::CalculateDerivativeVectors( rStressMatrix, V1, V2, V3);


            double K, dKdAngle;
            ShapeAtDeviatoricPlaneUtility::CalculateKLodeCoefficients( K, dKdAngle, LodeAngle, rDilatancy);

            double A1, A2, A3;
            A1 = sin( rDilatancy* Globals::Pi/180.0);
            A2 = K - dKdAngle * std::tan( 3.0 * LodeAngle);
            A3 = 2.0 * pow(J2, 2) * std::cos( 3.0 * LodeAngle);
            if ( fabs(A3) > 1e-12) {
               A3 = -sqrt(3)*dKdAngle/A3;
            } else {
               A3 = 0;
            }

            double a = 1.0;
            if ( rDilatancy > 0) {
               double alpha = sqrt( pow(J2*K, 2) + pow( a * std::sin(rDilatancy * Globals::Pi/180.0), 2) );
               alpha = J2*K/alpha;

               A2 *= alpha;
               A3 *= alpha;

               rDeltaStressYieldCondition  = A1 * V1 + A2 * V2 + A3 *V3 ;
            } else {
               const double & rCohesion = rMaterialProperties[COHESION];
               const double & rFriction = rMaterialProperties[INTERNAL_FRICTION_ANGLE];
               double Friction = rFriction* Globals::Pi / 180.0;
               double Point = MeanStress / std::tan(Friction) - a;
               double factor = std::exp(-Point + MeanStress);
               double Term = 2.0* fabs( std::sin(rDilatancy*Globals::Pi/180));

               
               rDeltaStressYieldCondition  = (A1 * V1 + A2 * V2 + A3 *V3)*(1.0-factor) + Term * factor * V1; 
            }

            /*std::cout << " p "  << MeanStress << " J2 " << J2 << " LodeAngle " << LodeAngle * 180.0/Globals::Pi << std::endl;
              std::cout << " A1 " << A1 << " V1 " << V1 << std::endl;
              std::cout << " A2 " << A2 << " V2 " << V2 << std::endl;
              std::cout << " A3 " << A3 << " V3 " << V3 << std::endl;
              std::cout << "                                                   SO: " << rDeltaStressYieldCondition << std::endl;
             */

            return rDeltaStressYieldCondition;

            KRATOS_CATCH(" ")
         }


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
         std::string Info() const override
         {
            std::stringstream buffer;
            buffer << "MohrCoulombPlasticPotential" ;
            return buffer.str();
         }

         /// Print information about this object.
         void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "MohrCoulombPlasticPotential";
         }

         /// Print object's data.
         void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "MohrCoulombPlasticPotential Data";
         }


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
         ///@name Serialization
         ///@{
         friend class Serializer;


         void save(Serializer& rSerializer) const override
         {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
         }

         void load(Serializer& rSerializer) override
         {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
         }

         ///@}
         ///@name Private Inquiry
         ///@{


         ///@}
         ///@name Un accessible methods
         ///@{

         ///@}

   }; // Class MohrCoulombPlasticPotential

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{


   ///@}

   ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MOHR_COULOMB_PLASTIC_POTENTIAL_H_INCLUDED  defined
