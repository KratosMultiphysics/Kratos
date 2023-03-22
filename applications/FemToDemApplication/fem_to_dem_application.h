//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//


#if !defined(KRATOS_FEM_TO_DEM_APPLICATION_H_INCLUDED )
#define  KRATOS_FEM_TO_DEM_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

//#include "custom_elements/zaratipito_element.hpp"
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "custom_constitutive/linear_plane_stress.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_plane_strain_2d.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_elements/generic_small_strain_femdem_element.hpp"
#include "custom_elements/generic_total_lagrangian_femdem_element.h"
#include "custom_elements/generic_total_lagrangian_mixtures_femdem_element.hpp"

#include "fem_to_dem_application_variables.h"

namespace Kratos {

class KRATOS_API(FEM_TO_DEM_APPLICATION) KratosFemToDemApplication : public KratosApplication 

{

public:
    
    KRATOS_CLASS_POINTER_DEFINITION(KratosFemToDemApplication);

    /// Default constructor.
    KratosFemToDemApplication();

    /// Destructor.
    virtual ~KratosFemToDemApplication(){}

    virtual void Register() override; 



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
    virtual std::string Info() const override {
        return "KratosFemToDemApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {
          KRATOS_WATCH("in my application");
          KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );

        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }


protected:

private:
    
    // Elements
    const GenericSmallStrainFemDemElement<2,0> mSmallStrainModifiedMohrCoulombFemDemElement2D;
    const GenericSmallStrainFemDemElement<3,0> mSmallStrainModifiedMohrCoulombFemDemElement3D;
    const GenericSmallStrainFemDemElement<2,1> mSmallStrainRankineFemDemElement2D;
    const GenericSmallStrainFemDemElement<3,1> mSmallStrainRankineFemDemElement3D;
    const GenericSmallStrainFemDemElement<2,2> mSmallStrainSimoJuFemDemElement2D;
    const GenericSmallStrainFemDemElement<3,2> mSmallStrainSimoJuFemDemElement3D;
    const GenericSmallStrainFemDemElement<2,3> mSmallStrainDruckerPragerFemDemElement2D;
    const GenericSmallStrainFemDemElement<3,3> mSmallStrainDruckerPragerFemDemElement3D;
    const GenericSmallStrainFemDemElement<2,4> mSmallStrainVonMisesFemDemElement2D;
    const GenericSmallStrainFemDemElement<3,4> mSmallStrainVonMisesFemDemElement3D;
    const GenericSmallStrainFemDemElement<2,5> mSmallStrainTrescaFemDemElement2D;
    const GenericSmallStrainFemDemElement<3,5> mSmallStrainTrescaFemDemElement3D;
    const GenericSmallStrainFemDemElement<2,6> mSmallStrainMohrCoulombFemDemElement2D;
    const GenericSmallStrainFemDemElement<3,6> mSmallStrainMohrCoulombFemDemElement3D;

    const GenericTotalLagrangianFemDemElement<2,0> mTotalLagrangianModifiedMohrCoulombFemDemElement2D;
    const GenericTotalLagrangianFemDemElement<3,0> mTotalLagrangianModifiedMohrCoulombFemDemElement3D;
    const GenericTotalLagrangianFemDemElement<2,1> mTotalLagrangianRankineFemDemElement2D;
    const GenericTotalLagrangianFemDemElement<3,1> mTotalLagrangianRankineFemDemElement3D;
    const GenericTotalLagrangianFemDemElement<2,2> mTotalLagrangianSimoJuFemDemElement2D;
    const GenericTotalLagrangianFemDemElement<3,2> mTotalLagrangianSimoJuFemDemElement3D;
    const GenericTotalLagrangianFemDemElement<2,3> mTotalLagrangianDruckerPragerFemDemElement2D;
    const GenericTotalLagrangianFemDemElement<3,3> mTotalLagrangianDruckerPragerFemDemElement3D;
    const GenericTotalLagrangianFemDemElement<2,4> mTotalLagrangianVonMisesFemDemElement2D;
    const GenericTotalLagrangianFemDemElement<3,4> mTotalLagrangianVonMisesFemDemElement3D;
    const GenericTotalLagrangianFemDemElement<2,5> mTotalLagrangianTrescaFemDemElement2D;
    const GenericTotalLagrangianFemDemElement<3,5> mTotalLagrangianTrescaFemDemElement3D;
    const GenericTotalLagrangianFemDemElement<2,6> mTotalLagrangianMohrCoulombFemDemElement2D;
    const GenericTotalLagrangianFemDemElement<3,6> mTotalLagrangianMohrCoulombFemDemElement3D;

    const GenericTotalLagrangianMixturesFemDemElement<2,0> mTotalLagrangianMixturesModifiedMohrCoulombFemDemElement2D;
    const GenericTotalLagrangianMixturesFemDemElement<3,0> mTotalLagrangianMixturesModifiedMohrCoulombFemDemElement3D;
    const GenericTotalLagrangianMixturesFemDemElement<2,1> mTotalLagrangianMixturesRankineFemDemElement2D;
    const GenericTotalLagrangianMixturesFemDemElement<3,1> mTotalLagrangianMixturesRankineFemDemElement3D;
    const GenericTotalLagrangianMixturesFemDemElement<2,2> mTotalLagrangianMixturesSimoJuFemDemElement2D;
    const GenericTotalLagrangianMixturesFemDemElement<3,2> mTotalLagrangianMixturesSimoJuFemDemElement3D;
    const GenericTotalLagrangianMixturesFemDemElement<2,3> mTotalLagrangianMixturesDruckerPragerFemDemElement2D;
    const GenericTotalLagrangianMixturesFemDemElement<3,3> mTotalLagrangianMixturesDruckerPragerFemDemElement3D;
    const GenericTotalLagrangianMixturesFemDemElement<2,4> mTotalLagrangianMixturesVonMisesFemDemElement2D;
    const GenericTotalLagrangianMixturesFemDemElement<3,4> mTotalLagrangianMixturesVonMisesFemDemElement3D;
    const GenericTotalLagrangianMixturesFemDemElement<2,5> mTotalLagrangianMixturesTrescaFemDemElement2D;
    const GenericTotalLagrangianMixturesFemDemElement<3,5> mTotalLagrangianMixturesTrescaFemDemElement3D;
    const GenericTotalLagrangianMixturesFemDemElement<2,6> mTotalLagrangianMixturesMohrCoulombFemDemElement2D;
    const GenericTotalLagrangianMixturesFemDemElement<3,6> mTotalLagrangianMixturesMohrCoulombFemDemElement3D;

    //Hiperelastic and elastic laws
   const LinearPlaneStrainFEMDEM mLinearPlaneStrainFEMDEM;
   const LinearPlaneStressFEMDEM mLinearPlaneStressFEMDEM;
   const ElasticIsotropic3DFEMDEM mElasticIsotropic3DFEMDEM;
   const HyperElasticIsotropicNeoHookeanPlaneStrain2DFEMDEM mHyperElasticIsotropicNeoHookeanPlaneStrain2DFEMDEM;
   const HyperElasticIsotropicNeoHookean3DFEMDEM mHyperElasticIsotropicNeoHookean3DFEMDEM;
    

    /// Assignment operator.
    KratosFemToDemApplication& operator=(KratosFemToDemApplication const& rOther);

    /// Copy constructor.
    KratosFemToDemApplication(KratosFemToDemApplication const& rOther);

}; // Class KratosFemToDemApplication

}  // namespace Kratos.

#endif // KRATOS_FEM_TO_DEM_APPLICATION_H_INCLUDED  defined
