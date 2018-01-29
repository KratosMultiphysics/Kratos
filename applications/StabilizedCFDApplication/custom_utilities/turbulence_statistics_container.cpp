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

#include "turbulence_statistics_container.h"
#include "stabilized_cfd_application_variables.h"
#include "includes/model_part.h"

#include <string>
#include <iostream>
#include <fstream>

namespace Kratos
{

// Version 0: Initial set of variables
// Version 1: Adding P-DUi/DXj DUi/DXj-DUi/DXj and DU/DXj-DV/DXj correlations to analyze variance budgets
int TurbulenceStatisticsContainer::mVersion = 1;

TurbulenceStatisticsContainer::TurbulenceStatisticsContainer(unsigned int NumGauss):
    mData(NumGauss, std::vector<double>(TurbulenceStatisticsContainer::NUM_DATA,0.0)),
    mInstantData(NumGauss, std::vector<double>(TurbulenceStatisticsContainer::NUM_INST_DATA,0.0))
{
}


TurbulenceStatisticsContainer::~TurbulenceStatisticsContainer()
{
}

void TurbulenceStatisticsContainer::AddStep(std::vector<double> &Values, unsigned int g, unsigned int NumSteps)
{
    std::vector<double>& rData = mData[g];
    // covariances are computed using an online formula depending on averages, the UiUj values are not read from input
    // Sum{n+1}(U^2) = Sum{n}(U^2) + (n*U - Sum{n}(U))^2/(n*(n+1))
    // A similar formula is used for triple correlations
    if (NumSteps <= 1)
    {
        Values[UU] = 0.0;
        Values[VV] = 0.0;
        Values[WW] = 0.0;
        Values[UV] = 0.0;
        Values[UW] = 0.0;
        Values[VW] = 0.0;

        Values[UP] = 0.0;
        Values[VP] = 0.0;
        Values[WP] = 0.0;

        Values[UUU] = 0.0;
        Values[UVV] = 0.0;
        Values[UWW] = 0.0;
        Values[VUU] = 0.0;
        Values[VVV] = 0.0;
        Values[VWW] = 0.0;
        Values[WUU] = 0.0;
        Values[WVV] = 0.0;
        Values[WWW] = 0.0;
        
        Values[PDU_DX] = 0.0;
        Values[PDU_DY] = 0.0;
        Values[PDU_DZ] = 0.0;
        Values[PDV_DX] = 0.0;
        Values[PDV_DY] = 0.0;
        Values[PDV_DZ] = 0.0;
        Values[PDW_DX] = 0.0;
        Values[PDW_DY] = 0.0;
        Values[PDW_DZ] = 0.0;

        Values[DU_DXDU_DX] = 0.0;
        Values[DU_DYDU_DY] = 0.0;
        Values[DU_DZDU_DZ] = 0.0;

        Values[DV_DXDV_DX] = 0.0;
        Values[DV_DYDV_DY] = 0.0;
        Values[DV_DZDV_DZ] = 0.0;

        Values[DW_DXDW_DX] = 0.0;
        Values[DW_DYDW_DY] = 0.0;
        Values[DW_DZDW_DZ] = 0.0;

        Values[DU_DXDV_DX] = 0.0;
        Values[DU_DYDV_DY] = 0.0;
        Values[DU_DZDV_DZ] = 0.0;

        Values[UssUss] = 0.0;
        Values[UssVss] = 0.0;
        Values[UssWss] = 0.0;
        Values[VssVss] = 0.0;
        Values[VssWss] = 0.0;
        Values[WssWss] = 0.0;
        Values[PssPss] = 0.0;

        Values[UUss] = 0.0;
        Values[UVss] = 0.0;
        Values[UWss] = 0.0;
        Values[VUss] = 0.0;
        Values[VVss] = 0.0;
        Values[VWss] = 0.0;
        Values[WUss] = 0.0;
        Values[WVss] = 0.0;
        Values[WWss] = 0.0;
    }
    else
    {
        const double UpdateFactor = 1.0 / ((NumSteps-1)*NumSteps);
        const double Du = (NumSteps-1)*Values[U] - rData[U];
        const double Dv = (NumSteps-1)*Values[V] - rData[V];
        const double Dw = (NumSteps-1)*Values[W] - rData[W];

        Values[UU] = UpdateFactor * Du *Du;
        Values[VV] = UpdateFactor * Dv *Dv;
        Values[WW] = UpdateFactor * Dw *Dw;
        Values[UV] = UpdateFactor * Du *Dv;
        Values[UW] = UpdateFactor * Du *Dw;
        Values[VW] = UpdateFactor * Dv *Dw;

        const double Dp = (NumSteps-1)*Values[P] - rData[P];
        Values[UP] = UpdateFactor * Du *Dp;
        Values[VP] = UpdateFactor * Dv *Dp;
        Values[WP] = UpdateFactor * Dw *Dp;

        const double M2uu = rData[UU];
        const double M2vv = rData[VV];
        const double M2ww = rData[WW];
        const double M2uv = rData[UV];
        const double M2uw = rData[UW];
        const double M2vw = rData[VW];

        const double F1 = (NumSteps-2)*UpdateFactor*UpdateFactor;
        Values[UUU] = F1*Du*Du*Du - UpdateFactor*3.0*M2uu*Du;
        Values[UVV] = F1*Du*Dv*Dv - UpdateFactor*(M2vv*Du + 2.0*M2uv*Dv);
        Values[UWW] = F1*Du*Dw*Dw - UpdateFactor*(M2ww*Du + 2.0*M2uw*Dw);

        Values[VUU] = F1*Dv*Du*Du - UpdateFactor*(M2uu*Dv + 2.0*M2uv*Dv);
        Values[VVV] = F1*Dv*Dv*Dv - UpdateFactor*3.0*M2vv*Dv;
        Values[VWW] = F1*Dv*Dw*Dw - UpdateFactor*(M2ww*Dv + 2.0*M2vw*Dw);

        Values[WUU] = F1*Dw*Du*Du - UpdateFactor*(M2uu*Dw + 2.0*M2uw*Dw);
        Values[WVV] = F1*Dw*Dv*Dv - UpdateFactor*(M2vv*Dw + 2.0*M2vw*Dw);
        Values[WWW] = F1*Dw*Dw*Dw - UpdateFactor*3.0*M2ww*Dw;

        const double Ddudx = (NumSteps-1)*Values[DU_DX] - rData[DU_DX];
        const double Ddudy = (NumSteps-1)*Values[DU_DY] - rData[DU_DY];
        const double Ddudz = (NumSteps-1)*Values[DU_DZ] - rData[DU_DZ];
        const double Ddvdx = (NumSteps-1)*Values[DV_DX] - rData[DV_DX];
        const double Ddvdy = (NumSteps-1)*Values[DV_DY] - rData[DV_DY];
        const double Ddvdz = (NumSteps-1)*Values[DV_DZ] - rData[DV_DZ];
        const double Ddwdx = (NumSteps-1)*Values[DW_DX] - rData[DW_DX];
        const double Ddwdy = (NumSteps-1)*Values[DW_DY] - rData[DW_DY];
        const double Ddwdz = (NumSteps-1)*Values[DW_DZ] - rData[DW_DZ];

        Values[PDU_DX] = UpdateFactor * Dp *Ddudx;
        Values[PDU_DY] = UpdateFactor * Dp *Ddudy;
        Values[PDU_DZ] = UpdateFactor * Dp *Ddudz;
        Values[PDV_DX] = UpdateFactor * Dp *Ddvdx;
        Values[PDV_DY] = UpdateFactor * Dp *Ddvdy;
        Values[PDV_DZ] = UpdateFactor * Dp *Ddvdz;
        Values[PDW_DX] = UpdateFactor * Dp *Ddwdx;
        Values[PDW_DY] = UpdateFactor * Dp *Ddwdy;
        Values[PDW_DZ] = UpdateFactor * Dp *Ddwdz;

        Values[DU_DXDU_DX] = UpdateFactor * Ddudx*Ddudx;
        Values[DU_DYDU_DY] = UpdateFactor * Ddudy*Ddudy;
        Values[DU_DZDU_DZ] = UpdateFactor * Ddudz*Ddudz;

        Values[DV_DXDV_DX] = UpdateFactor * Ddvdx*Ddvdx;
        Values[DV_DYDV_DY] = UpdateFactor * Ddvdy*Ddvdy;
        Values[DV_DZDV_DZ] = UpdateFactor * Ddvdz*Ddvdz;

        Values[DW_DXDW_DX] = UpdateFactor * Ddwdx*Ddwdx;
        Values[DW_DYDW_DY] = UpdateFactor * Ddwdy*Ddwdy;
        Values[DW_DZDW_DZ] = UpdateFactor * Ddwdz*Ddwdz;

        Values[DU_DXDV_DX] = UpdateFactor * Ddudx*Ddvdx;
        Values[DU_DYDV_DY] = UpdateFactor * Ddudy*Ddvdy;
        Values[DU_DZDV_DZ] = UpdateFactor * Ddudz*Ddvdz;

        const double Duss = (NumSteps-1)*Values[Uss] - rData[Uss];
        const double Dvss = (NumSteps-1)*Values[Vss] - rData[Vss];
        const double Dwss = (NumSteps-1)*Values[Wss] - rData[Wss];
        const double Dpss = (NumSteps-1)*Values[Pss] - rData[Pss];

        Values[UssUss] = UpdateFactor * Duss *Duss;
        Values[UssVss] = UpdateFactor * Duss *Dvss;
        Values[UssWss] = UpdateFactor * Duss *Dwss;
        Values[VssVss] = UpdateFactor * Dvss *Dvss;
        Values[VssWss] = UpdateFactor * Dvss *Dwss;
        Values[WssWss] = UpdateFactor * Dwss *Dwss;

        Values[PssPss] = UpdateFactor * Dpss *Dpss;

        Values[UUss] = UpdateFactor * Du *Duss;
        Values[UVss] = UpdateFactor * Du *Dvss;
        Values[UWss] = UpdateFactor * Du *Dwss;
        Values[VUss] = UpdateFactor * Dv *Duss;
        Values[VVss] = UpdateFactor * Dv *Dvss;
        Values[VWss] = UpdateFactor * Dv *Dwss;
        Values[WUss] = UpdateFactor * Dw *Duss;
        Values[WVss] = UpdateFactor * Dw *Dvss;
        Values[WWss] = UpdateFactor * Dw *Dwss;
    }

    for (unsigned int i = 0; i < TurbulenceStatisticsContainer::NUM_DATA; i++)
        rData[i] += Values[i];
}

void TurbulenceStatisticsContainer::StoreInstantData(std::vector<double> &Values, unsigned int g)
{
    std::vector<double>& rData = mInstantData[g];

    for (unsigned int i = 0; i < TurbulenceStatisticsContainer::NUM_INST_DATA; i++)
        rData[i] = Values[i];
}

double TurbulenceStatisticsContainer::GetValue(StatisticalData Variable, unsigned int g, unsigned int NumSteps) const
{
    return (mData[g][Variable])/NumSteps;
}

void TurbulenceStatisticsContainer::DumpData(ModelPart &rModelPart)
{
    // Open output file
    std::stringstream FileName;
    FileName << "gp_statistics_" << rModelPart.GetCommunicator().MyPID() << ".csv";
    std::ofstream StatsFile;
    StatsFile.open(FileName.str().c_str(), std::ios::out | std::ios::trunc);

    unsigned int NumSteps = rModelPart.GetProcessInfo().GetValue(RECORDED_STEPS);

    for (ModelPart::ElementIterator it = rModelPart.GetCommunicator().LocalMesh().ElementsBegin();
         it != rModelPart.GetCommunicator().LocalMesh().ElementsEnd(); it++)
    {
        const Matrix& NContainer = it->GetGeometry().ShapeFunctionsValues(it->GetIntegrationMethod());
        const unsigned int NumNodes = it->GetGeometry().PointsNumber();
        const unsigned int NumGauss = it->GetGeometry().IntegrationPointsNumber(it->GetIntegrationMethod());

        const TurbulenceStatisticsContainer& rContainer = *( it->GetValue(TURBULENCE_STATISTICS) );

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            StatsFile << it->Id() << " " << g << " ";

            // Print integration point coordinates
            array_1d<double,3> Coords(3,0.0);
            for (unsigned int n = 0; n < NumNodes; n++)
                Coords += NContainer(g,n) * it->GetGeometry()[n].Coordinates();
            StatsFile << Coords[0] << " " << Coords[1] << " " << Coords[2] << " ";

            //  Print data
            for (unsigned int i = 0; i < NUM_DATA; i++)
                StatsFile << rContainer.mData[g][i]/NumSteps << " ";
            StatsFile << std::endl;
        }
    }

    StatsFile.close();
}

void TurbulenceStatisticsContainer::DumpInstantData(ModelPart &rModelPart)
{
    // Open output file
    std::stringstream FileName;
    const TurbulenceStatisticsContainer& rAuxContainer = *( (rModelPart.GetCommunicator().LocalMesh().ElementsBegin())->GetValue(TURBULENCE_STATISTICS) );
    FileName << "gp_instant_data_" << rAuxContainer.mInstantData[0][Time_Label] << "_" << rModelPart.GetCommunicator().MyPID() << ".csv";
    std::ofstream DataFile;
    DataFile.open(FileName.str().c_str(), std::ios::out | std::ios::trunc);

    for (ModelPart::ElementIterator it = rModelPart.GetCommunicator().LocalMesh().ElementsBegin();
         it != rModelPart.GetCommunicator().LocalMesh().ElementsEnd(); it++)
    {
        const Matrix& NContainer = it->GetGeometry().ShapeFunctionsValues(it->GetIntegrationMethod());
        const unsigned int NumNodes = it->GetGeometry().PointsNumber();
        const unsigned int NumGauss = it->GetGeometry().IntegrationPointsNumber(it->GetIntegrationMethod());

        const TurbulenceStatisticsContainer& rContainer = *( it->GetValue(TURBULENCE_STATISTICS) );

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            DataFile << it->Id() << " " << g << " ";

            // Print integration point coordinates
            array_1d<double,3> Coords(3,0.0);
            for (unsigned int n = 0; n < NumNodes; n++)
                Coords += NContainer(g,n) * it->GetGeometry()[n].Coordinates();
            DataFile << Coords[0] << " " << Coords[1] << " " << Coords[2] << " ";

            //  Print data
            for (unsigned int i = 1; i < NUM_INST_DATA; i++)
                DataFile << rContainer.mInstantData[g][i] << " ";
            DataFile << std::endl;
        }
    }

    DataFile.close();
}

void TurbulenceStatisticsContainer::save(Serializer &rSerializer) const
{
    rSerializer.save("mVersion",mVersion);
    rSerializer.save("NumGauss",mData.size());
    rSerializer.save("mData",mData);
}

void TurbulenceStatisticsContainer::load(Serializer &rSerializer)
{
    int LoadedDataVersion;
    int NumGauss;
    rSerializer.load("mVersion",LoadedDataVersion);
    rSerializer.load("NumGauss",NumGauss);

    if (LoadedDataVersion == mVersion)
    {
        rSerializer.load("mData",mData);
    }
    else
    {
        std::cout << "WARNING TurbulentStatisticsContainer:" << std::endl << "Ignoring read statististics: they were created using a different version" << std::endl;
        std::vector< std::vector<double> > LoadedData;
        rSerializer.load("mData",LoadedData);
        mData = std::vector< std::vector<double> >(NumGauss,std::vector<double>(NUM_DATA,0.0));
    }

    mInstantData = std::vector< std::vector<double> >(NumGauss,std::vector<double>(NUM_INST_DATA,0.0));
}

}
