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

#if !defined(KRATOS_INITIAL_STRESS_2D_UTILITIES )
#define  KRATOS_INITIAL_STRESS_2D_UTILITIES

// System includes
#include <fstream>
#include <iostream>
#include <cmath>

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/math_utils.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

class InitialStress2DUtilities
{

protected:

/// Basic Structs for the utility ---------------------------------------------------------------------------------------------------------------------------------------------

    struct UtilityVariables
    {
        double X_max, X_min, Y_max, Y_min;
        int NRows, NColumns;
        double RowSize, ColumnSize;
    };

public:

    KRATOS_CLASS_POINTER_DEFINITION( InitialStress2DUtilities );

    /// Constructor
    InitialStress2DUtilities() {}

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~InitialStress2DUtilities() {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void TransferInitialStresses (ModelPart& rInitialModelPart, ModelPart& rCurrentModelPart, const bool& constant_discretization)
    {
        if (constant_discretization == true) {
            this->AssignNodalVariables(rInitialModelPart,rCurrentModelPart);
        } else {
            // Define necessary variables
            UtilityVariables AuxVariables;

            this->ComputeCellMatrixDimensions(AuxVariables,rCurrentModelPart);

            this->NodalVariablesMapping(AuxVariables,rInitialModelPart,rCurrentModelPart);
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SaveInitialStresses (Parameters& rParameters, ModelPart& rCurrentModelPart)
    {
        std::fstream initial_stresses_mdpa;
        std::string initial_stress_mdpa_path = rParameters["initial_input_filename"].GetString() + ".mdpa";

        initial_stresses_mdpa.open(initial_stress_mdpa_path.c_str(), std::fstream::out | std::fstream::trunc);
        initial_stresses_mdpa.precision(15);

        initial_stresses_mdpa << "Begin Properties 0" << std::endl;
        initial_stresses_mdpa << "End Properties" << std::endl << std::endl;

        const int NNodes = static_cast<int>(rCurrentModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = rCurrentModelPart.NodesBegin();
        initial_stresses_mdpa << "Begin Nodes" << std::endl;
        // #pragma omp parallel for
        for(int i = 0; i < NNodes; i++) {
            ModelPart::NodesContainerType::iterator it_node = node_begin + i;
            initial_stresses_mdpa << it_node->Id() << " " <<
                it_node->X0() << " " <<
                it_node->Y0() << " " <<
                it_node->Z0() << std::endl;
        }
        initial_stresses_mdpa << "End Nodes" << std::endl << std::endl;

        this->WriteLinearElements(initial_stresses_mdpa,rCurrentModelPart);

        this->CalculateNodalStresses(rCurrentModelPart);
        Matrix InitialStressTensor(3,3);
        initial_stresses_mdpa << "Begin NodalData INITIAL_STRESS_TENSOR" << std::endl;
        // #pragma omp parallel for
        for(int i = 0; i < NNodes; i++) {
            ModelPart::NodesContainerType::iterator it_node = node_begin + i;
            noalias(InitialStressTensor) = it_node->FastGetSolutionStepValue(NODAL_EFFECTIVE_STRESS_TENSOR);
            initial_stresses_mdpa << it_node->Id() << " 0 [2,2] ((" <<
                InitialStressTensor(0,0) << "," << InitialStressTensor(0,1) << "),(" <<
                InitialStressTensor(1,0) << "," << InitialStressTensor(1,1) << "))" << std::endl;
        }
        initial_stresses_mdpa << "End NodalData" << std::endl;

        initial_stresses_mdpa.close();
    }

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void NodalVariablesMapping(
        const UtilityVariables& AuxVariables,
        ModelPart& rInitialModelPart,
        ModelPart& rCurrentModelPart)
    {
        // Define ElementOld Cell matrix
        std::vector< std::vector< std::vector<Element::Pointer> > > ElementOldCellMatrix;
        ElementOldCellMatrix.resize(AuxVariables.NRows);
        for(int i = 0; i < AuxVariables.NRows; i++) ElementOldCellMatrix[i].resize(AuxVariables.NColumns);

        // Locate Old Elments in cells
        double X_me;
        double Y_me;
        int PointsNumber;

        int NElems = static_cast<int>(rInitialModelPart.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = rInitialModelPart.ElementsBegin();

        #pragma omp parallel for private(X_me,Y_me,PointsNumber)
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElemOld = el_begin + i;

            double X_left = itElemOld->GetGeometry().GetPoint(0).X0();
            double X_right = X_left;
            double Y_top = itElemOld->GetGeometry().GetPoint(0).Y0();
            double Y_bot = Y_top;
            PointsNumber = itElemOld->GetGeometry().PointsNumber();

            for(int j = 1; j < PointsNumber; j++)
            {
                X_me = itElemOld->GetGeometry().GetPoint(j).X0();
                Y_me = itElemOld->GetGeometry().GetPoint(j).Y0();

                if(X_me > X_right) X_right = X_me;
                else if(X_me < X_left) X_left = X_me;
                if(Y_me > Y_top) Y_top = Y_me;
                else if(Y_me < Y_bot) Y_bot = Y_me;
            }

            int Column_left = int((X_left-AuxVariables.X_min)/AuxVariables.ColumnSize);
            int Column_right = int((X_right-AuxVariables.X_min)/AuxVariables.ColumnSize);
            int Row_top = int((AuxVariables.Y_max-Y_top)/AuxVariables.RowSize);
            int Row_bot = int((AuxVariables.Y_max-Y_bot)/AuxVariables.RowSize);

            if(Column_left < 0) Column_left = 0;
            else if(Column_left >= AuxVariables.NColumns) Column_left = AuxVariables.NColumns-1;
            if(Column_right >= AuxVariables.NColumns) Column_right = AuxVariables.NColumns-1;
            else if(Column_right < 0) Column_right = 0;

            if(Row_top < 0) Row_top = 0;
            else if(Row_top >= AuxVariables.NRows) Row_top = AuxVariables.NRows-1;
            if(Row_bot >= AuxVariables.NRows) Row_bot = AuxVariables.NRows-1;
            else if(Row_bot < 0) Row_bot = 0;

            for(int k = Row_top; k <= Row_bot; k++)
            {
                for(int l = Column_left; l <= Column_right; l++)
                {
                    #pragma omp critical
                    {
                        ElementOldCellMatrix[k][l].push_back((*(itElemOld.base())));
                    }
                }
            }
        }

        // Locate current nodes inside initial elements and interpolate nodal variables
        const int NNodes = static_cast<int>(rCurrentModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = rCurrentModelPart.NodesBegin();

        array_1d<double,3> GlobalCoordinates;
        array_1d<double,3> LocalCoordinates;

        #pragma omp parallel for private(X_me,Y_me,PointsNumber,GlobalCoordinates,LocalCoordinates)
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNodeNew = node_begin + i;

            X_me = itNodeNew->X0();
            Y_me = itNodeNew->Y0();

            int Column = int((X_me-AuxVariables.X_min)/AuxVariables.ColumnSize);
            int Row = int((AuxVariables.Y_max-Y_me)/AuxVariables.RowSize);

            if(Column >= AuxVariables.NColumns) Column = AuxVariables.NColumns-1;
            else if(Column < 0) Column = 0;
            if(Row >= AuxVariables.NRows) Row = AuxVariables.NRows-1;
            else if(Row < 0) Row = 0;

            noalias(GlobalCoordinates) = itNodeNew->Coordinates(); //Coordinates of new nodes are still in the original position
            bool IsInside = false;
            Element::Pointer pElementOld;

            for(unsigned int m = 0; m < (ElementOldCellMatrix[Row][Column]).size(); m++)
            {
                pElementOld = ElementOldCellMatrix[Row][Column][m];
                IsInside = pElementOld->GetGeometry().IsInside(GlobalCoordinates,LocalCoordinates); //Checks whether the global coordinates fall inside the original old element
                if(IsInside) break;
            }
            if(IsInside==false)
            {
                for(unsigned int m = 0; m < (ElementOldCellMatrix[Row][Column]).size(); m++)
                {
                    pElementOld = ElementOldCellMatrix[Row][Column][m];
                    IsInside = pElementOld->GetGeometry().IsInside(GlobalCoordinates,LocalCoordinates,1.0e-5);
                    if(IsInside) break;
                }
            }
            if(IsInside == false)
                std::cout << "ERROR!!, NONE OF THE OLD ELEMENTS CONTAINS NODE: " << itNodeNew->Id() << std::endl;

            PointsNumber = pElementOld->GetGeometry().PointsNumber();

            Vector ShapeFunctionsValuesVector(PointsNumber);
            pElementOld->GetGeometry().ShapeFunctionsValues(ShapeFunctionsValuesVector,LocalCoordinates);

            // Interpolation of nodal variables
            std::vector<Vector> ComponentsNodalVariableVector(3); // voigt_size
            for(int j = 0; j < 3; j++) { // voigt_size
                ComponentsNodalVariableVector[j].resize(PointsNumber);
            }
            Vector nodal_initial_stress_vector(3); // voigt_size
            Matrix nodal_initial_stress_tensor(2,2); // dimension
            for(int j = 0; j < PointsNumber; j++) {
                noalias(nodal_initial_stress_tensor) = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(INITIAL_STRESS_TENSOR);
                noalias(nodal_initial_stress_vector) = MathUtils<double>::StressTensorToVector(nodal_initial_stress_tensor);
                for (int k= 0; k < 3; k++) { // voigt_size
                    ComponentsNodalVariableVector[k][j] = nodal_initial_stress_vector[k];
                }
            }
            for (int k= 0; k < 3; k++) { // voigt_size
                nodal_initial_stress_vector[k] = inner_prod(ShapeFunctionsValuesVector,ComponentsNodalVariableVector[k]);
            }
            noalias(nodal_initial_stress_tensor) = MathUtils<double>::StressVectorToTensor(nodal_initial_stress_vector);
            Matrix& rNodalStress = itNodeNew->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR);
            if(rNodalStress.size1() != 2) // Dimension
                rNodalStress.resize(2,2,false);
            noalias(rNodalStress) = nodal_initial_stress_tensor;
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AssignNodalVariables(
        ModelPart& rInitialModelPart,
        ModelPart& rCurrentModelPart) {

        // Here rInitialModelPart and rCurrentModelPart have the same exact discretization 
        const int NNodes = static_cast<int>(rInitialModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin_old = rInitialModelPart.NodesBegin();
        ModelPart::NodesContainerType::iterator node_begin = rCurrentModelPart.NodesBegin();

        #pragma omp parallel for
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNodeOld = node_begin_old + i;
            ModelPart::NodesContainerType::iterator itNodeNew = node_begin + i;

            Matrix& rNodalStress = itNodeNew->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR);
            if(rNodalStress.size1() != 2) // Dimension
                rNodalStress.resize(2,2,false);
            noalias(rNodalStress) = itNodeOld->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR);
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

/// Common --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ComputeCellMatrixDimensions(
        UtilityVariables& rAuxVariables,
        ModelPart& rModelPart)
    {
        // Compute X and Y limits of the current geometry
        unsigned int NumThreads = ParallelUtilities::GetNumThreads();
        std::vector<double> X_max_partition(NumThreads);
        std::vector<double> X_min_partition(NumThreads);
        std::vector<double> Y_max_partition(NumThreads);
        std::vector<double> Y_min_partition(NumThreads);

        const int NNodes = static_cast<int>(rModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = rModelPart.NodesBegin();

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            X_max_partition[k] = node_begin->X();
            X_min_partition[k] = X_max_partition[k];
            Y_max_partition[k] = node_begin->Y();
            Y_min_partition[k] = Y_max_partition[k];

            double X_me, Y_me;

            #pragma omp for
            for(int i = 0; i < NNodes; i++)
            {
                ModelPart::NodesContainerType::iterator itNode = node_begin + i;

                X_me = itNode->X();
                Y_me = itNode->Y();

                if( X_me > X_max_partition[k] ) X_max_partition[k] = X_me;
                else if( X_me < X_min_partition[k] ) X_min_partition[k] = X_me;

                if( Y_me > Y_max_partition[k] ) Y_max_partition[k] = Y_me;
                else if( Y_me < Y_min_partition[k] ) Y_min_partition[k] = Y_me;
            }
        }

        rAuxVariables.X_max = X_max_partition[0];
        rAuxVariables.X_min = X_min_partition[0];
        rAuxVariables.Y_max = Y_max_partition[0];
        rAuxVariables.Y_min = Y_min_partition[0];

        for(unsigned int i=1; i < NumThreads; i++)
        {
            if(X_max_partition[i] > rAuxVariables.X_max) rAuxVariables.X_max = X_max_partition[i];
            if(X_min_partition[i] < rAuxVariables.X_min) rAuxVariables.X_min = X_min_partition[i];
            if(Y_max_partition[i] > rAuxVariables.Y_max) rAuxVariables.Y_max = Y_max_partition[i];
            if(Y_min_partition[i] < rAuxVariables.Y_min) rAuxVariables.Y_min = Y_min_partition[i];
        }

        // Calculate Average Element Length
        double AverageElementLength = 0.0;

        int NElems = static_cast<int>(rModelPart.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();

        #pragma omp parallel for reduction(+:AverageElementLength)
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = el_begin + i;

            AverageElementLength += itElem->GetGeometry().Length();
        }

        AverageElementLength = AverageElementLength/NElems;

        // Compute FracturePoints CellMatrix dimensions
        rAuxVariables.NRows = int((rAuxVariables.Y_max-rAuxVariables.Y_min)/(AverageElementLength));
        rAuxVariables.NColumns = int((rAuxVariables.X_max-rAuxVariables.X_min)/(AverageElementLength));
        if(rAuxVariables.NRows <= 0) rAuxVariables.NRows = 1;
        if(rAuxVariables.NColumns <= 0) rAuxVariables.NColumns = 1;
        rAuxVariables.RowSize = (rAuxVariables.Y_max-rAuxVariables.Y_min)/rAuxVariables.NRows;
        rAuxVariables.ColumnSize = (rAuxVariables.X_max-rAuxVariables.X_min)/rAuxVariables.NColumns;
    }

///------------------------------------------------------------------------------------

    void CalculateNodalStresses(ModelPart& rCurrentModelPart) {

        // Clear nodal Stresses
        const int NNodes = static_cast<int>(rCurrentModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = rCurrentModelPart.NodesBegin();

        #pragma omp parallel for
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator it_node = node_begin + i;
            it_node->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
            Matrix& rNodalStress = it_node->FastGetSolutionStepValue(NODAL_EFFECTIVE_STRESS_TENSOR);
            if(rNodalStress.size1() != 3)
                rNodalStress.resize(3,3,false);
            noalias(rNodalStress) = ZeroMatrix(3,3);
        }

        // Calculate and Extrapolate Stresses
        const ProcessInfo& r_current_process_info = rCurrentModelPart.GetProcessInfo();
        const auto it_elem_begin = rCurrentModelPart.ElementsBegin();
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rCurrentModelPart.Elements().size()); ++i) {
            auto it_elem = it_elem_begin + i;
            it_elem->FinalizeSolutionStep(r_current_process_info);
        }

        // Compute smoothed Stresses
        #pragma omp parallel for
        for(int n = 0; n < NNodes; n++)
        {
            ModelPart::NodesContainerType::iterator it_node = node_begin + n;
            const double& NodalArea = it_node->FastGetSolutionStepValue(NODAL_AREA);
            if (NodalArea>1.0e-20)
            {
                const double InvNodalArea = 1.0/NodalArea;
                Matrix& rNodalStress = it_node->FastGetSolutionStepValue(NODAL_EFFECTIVE_STRESS_TENSOR);
                for(unsigned int i = 0; i<3; i++) // Dimension
                {
                    for(unsigned int j = 0; j<3; j++)
                    {
                        rNodalStress(i,j) *= InvNodalArea;
                    }
                }
            }
        }
    }

///------------------------------------------------------------------------------------

    void WriteLinearElements(std::fstream& rInitialStressesMdpa, ModelPart& rCurrentModelPart) {

        int NElems = static_cast<int>(rCurrentModelPart.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = rCurrentModelPart.ElementsBegin();
        // NOTE: we are assuming that only one type of element is used (not quadratic elements)
        int PointsNumber = el_begin->GetGeometry().PointsNumber();
        if(PointsNumber == 3){
            rInitialStressesMdpa << "Begin Elements Element2D3N" << std::endl;
            // #pragma omp parallel for
            for(int i = 0; i < NElems; i++)
            {
                ModelPart::ElementsContainerType::iterator it_elem = el_begin + i;
                rInitialStressesMdpa << it_elem->Id() << " 0 " <<
                    it_elem->GetGeometry().GetPoint(0).Id() << " " <<
                    it_elem->GetGeometry().GetPoint(1).Id() << " " <<
                    it_elem->GetGeometry().GetPoint(2).Id() << std::endl;
            }
        } else {
            rInitialStressesMdpa << "Begin Elements Element2D4N" << std::endl;
            // #pragma omp parallel for
            for(int i = 0; i < NElems; i++)
            {
                ModelPart::ElementsContainerType::iterator it_elem = el_begin + i;
                rInitialStressesMdpa << it_elem->Id() << " 0 " <<
                    it_elem->GetGeometry().GetPoint(0).Id() << " " <<
                    it_elem->GetGeometry().GetPoint(1).Id() << " " <<
                    it_elem->GetGeometry().GetPoint(2).Id() << " " <<
                    it_elem->GetGeometry().GetPoint(3).Id() << std::endl;
            }
        }
        rInitialStressesMdpa << "End Elements" << std::endl << std::endl;
    }

}; // Class InitialStress2DUtilities

} // namespace Kratos.

#endif /* KRATOS_INITIAL_STRESS_2D_UTILITIES defined */
