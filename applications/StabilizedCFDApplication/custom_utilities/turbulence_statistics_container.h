//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#ifndef KRATOS_TURBULENT_STATISTICS_CONTAINER_H
#define KRATOS_TURBULENT_STATISTICS_CONTAINER_H

#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"


namespace Kratos
{

///@addtogroup StabilizedCFDApplication
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

class TurbulenceStatisticsContainer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FractionalStep
    KRATOS_CLASS_POINTER_DEFINITION(TurbulenceStatisticsContainer);

    enum StatisticalData
    {
        U,      // Means
        V,
        W,
        P,
        DU_DX,  // Gradients
        DU_DY,
        DU_DZ,
        DV_DX,
        DV_DY,
        DV_DZ,
        DW_DX,
        DW_DY,
        DW_DZ,
        Uss,    // Mean residuals
        Vss,
        Wss,
        Pss,
        UU,     // Large scale correlations
        UV,
        UW,
        VV,
        VW,
        WW,
        UP,     // Pressure correlations
        VP,
        WP,
        UUU,    // Triple correlations
        UVV,
        UWW,
        VUU,
        VVV,
        VWW,
        WUU,
        WVV,
        WWW,
        PDU_DX, // Pressure-gradient correlations
        PDU_DY,
        PDU_DZ,
        PDV_DX,
        PDV_DY,
        PDV_DZ,
        PDW_DX,
        PDW_DY,
        PDW_DZ,
        DU_DXDU_DX, // Dissipations
        DU_DYDU_DY,
        DU_DZDU_DZ,
        DV_DXDV_DX,
        DV_DYDV_DY,
        DV_DZDV_DZ,
        DW_DXDW_DX,
        DW_DYDW_DY,
        DW_DZDW_DZ,
        DU_DXDV_DX,
        DU_DYDV_DY,
        DU_DZDV_DZ,
        UssUss,     // Residual stresses
        UssVss,
        UssWss,
        VssVss,
        VssWss,
        WssWss,
        PssPss,
        UUss,       // Cross stresses
        UVss,
        UWss,
        VUss,
        VVss,
        VWss,
        WUss,
        WVss,
        WWss,
        NUM_DATA
    };

    enum InstantData
    {
        Time_Label,
        Tau_Momentum,
        Tau_Mass,
        Tau_FIC,
        Tau_VelNorm,
        H_effective,
        H_vel,
        H_min,
        H_max,
        Div_error,
//        DavgU_DX,
//        DavgU_DY,
//        DavgU_DZ,
//        DavgV_DX,
//        DavgV_DY,
//        DavgV_DZ,
//        DavgW_DX,
//        DavgW_DY,
//        DavgW_DZ,
//        MEAN_KINETIC_ENERGY,
//        TURBULENCE_KINETIC_ENERGY,
        NUM_INST_DATA
    };

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

//    TurbulenceStatisticsContainer();

    TurbulenceStatisticsContainer(unsigned int NumGauss);

    ~TurbulenceStatisticsContainer();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Initialize(unsigned int NumGauss);

    void AddStep(std::vector<double> &Values, unsigned int g, unsigned int NumSteps);

    void StoreInstantData(std::vector<double> &Values, unsigned int g);

    static void DumpData(ModelPart &rModelPart);

    static void DumpInstantData(ModelPart &rModelPart);

    ///@}
    ///@name Access
    ///@{

    double GetValue(StatisticalData Variable, unsigned int g, unsigned int NumSteps) const;

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

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

    static int mVersion;

    ///@}
    ///@name Member Variables
    ///@{

    std::vector< std::vector<double> > mData;

    std::vector< std::vector<double> > mInstantData;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    TurbulenceStatisticsContainer(){}

    void save(Serializer& rSerializer) const;

    void load(Serializer& rSerializer);


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


    TurbulenceStatisticsContainer(TurbulenceStatisticsContainer const& rOther);

    TurbulenceStatisticsContainer& operator=(TurbulenceStatisticsContainer const& rOther);

    ///@}


}; // Class TurbulentStatisticsContainer

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >>(std::istream& rIStream,
                                 TurbulenceStatisticsContainer& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator <<(std::ostream& rOStream,
                                 TurbulenceStatisticsContainer& rThis)
{
//    rThis.PrintInfo(rOStream);
//    rOStream << std::endl;
//    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_TURBULENT_STATISTICS_CONTAINER_H


