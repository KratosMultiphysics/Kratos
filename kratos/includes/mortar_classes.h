//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_MORTAR_CLASSES )
#define  KRATOS_MORTAR_CLASSES

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "utilities/mortar_utilities.h"
#include "includes/mapping_variables.h"

namespace Kratos 
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
    typedef Point                                                PointType;
    typedef Node<3>                                               NodeType;
    typedef Geometry<NodeType>                                GeometryType;
    
    // Type definition for integration methods
    typedef GeometryType::IntegrationPointsArrayType IntegrationPointsType;
    
///@}
///@name  Enum's
///@{
/**
 * This enum defines a "hash" used to identify in which combination of cuts the point is found when the mortar segmentation is performed
 */
#if !defined(POINT_BELONGS)
#define POINT_BELONGS
    
    enum PointBelongs 
    {
        Master       = 0, 
        Slave        = 1, 
        Intersection = 2
    };
    
    enum PointBelongsLine2D2N 
    {
        MasterLine2D2N0      = 0, 
        MasterLine2D2N1      = 1, 
        SlaveLine2D2N0       = 2, 
        SlaveLine2D2N1       = 3, 
        IntersectionLine2D2N = 4
    };
    
    enum PointBelongsTriangle3D3N 
    {
        MasterTriangle3D3N0          = 0, 
        MasterTriangle3D3N1          = 1, 
        MasterTriangle3D3N2          = 2, 
        SlaveTriangle3D3N0           = 3, 
        SlaveTriangle3D3N1           = 4, 
        SlaveTriangle3D3N2           = 5, 
        IntersectionTriangle3D3N     = 6,
        IntersectionTriangle3D3N0101 = 10106,
        IntersectionTriangle3D3N1001 = 10016,
        IntersectionTriangle3D3N1201 = 10216,
        IntersectionTriangle3D3N2101 = 10126,
        IntersectionTriangle3D3N0201 = 10206,
        IntersectionTriangle3D3N2001 = 10026,
        IntersectionTriangle3D3N0110 = 1106,
        IntersectionTriangle3D3N1010 = 1016,
        IntersectionTriangle3D3N1210 = 1216,
        IntersectionTriangle3D3N2110 = 1126,
        IntersectionTriangle3D3N0210 = 1206,
        IntersectionTriangle3D3N2010 = 1026,
        IntersectionTriangle3D3N0112 = 21106,
        IntersectionTriangle3D3N1012 = 21016,
        IntersectionTriangle3D3N1212 = 21216,
        IntersectionTriangle3D3N2112 = 21126,
        IntersectionTriangle3D3N0212 = 21206,
        IntersectionTriangle3D3N2012 = 21026,
        IntersectionTriangle3D3N0121 = 12106,
        IntersectionTriangle3D3N1021 = 12016,
        IntersectionTriangle3D3N1221 = 12216,
        IntersectionTriangle3D3N2121 = 12126,
        IntersectionTriangle3D3N0221 = 12206,
        IntersectionTriangle3D3N2021 = 12026,
        IntersectionTriangle3D3N0102 = 20106,
        IntersectionTriangle3D3N1002 = 20016,
        IntersectionTriangle3D3N1202 = 20216,
        IntersectionTriangle3D3N2102 = 20126,
        IntersectionTriangle3D3N0202 = 20206,
        IntersectionTriangle3D3N2002 = 20026,
        IntersectionTriangle3D3N0120 = 2106,
        IntersectionTriangle3D3N1020 = 2016,
        IntersectionTriangle3D3N1220 = 2216,
        IntersectionTriangle3D3N2120 = 2126,
        IntersectionTriangle3D3N0220 = 2206,
        IntersectionTriangle3D3N2020 = 2026
    };
    
    enum PointBelongsQuadrilateral3D4N 
    {
        MasterQuadrilateral3D4N0          = 0, 
        MasterQuadrilateral3D4N1          = 1, 
        MasterQuadrilateral3D4N2          = 2, 
        MasterQuadrilateral3D4N3          = 3, 
        SlaveQuadrilateral3D4N0           = 4, 
        SlaveQuadrilateral3D4N1           = 5, 
        SlaveQuadrilateral3D4N2           = 6, 
        SlaveQuadrilateral3D4N3           = 7, 
        IntersectionQuadrilateral3D4N     = 8,
        IntersectionQuadrilateral3D4N0101 = 10108,
        IntersectionQuadrilateral3D4N1001 = 10018,
        IntersectionQuadrilateral3D4N1201 = 10218,
        IntersectionQuadrilateral3D4N2101 = 10128,
        IntersectionQuadrilateral3D4N2301 = 10328,
        IntersectionQuadrilateral3D4N3201 = 10238,
        IntersectionQuadrilateral3D4N3001 = 10038,
        IntersectionQuadrilateral3D4N0301 = 10308,
        IntersectionQuadrilateral3D4N0110 = 1108,
        IntersectionQuadrilateral3D4N1010 = 1018,
        IntersectionQuadrilateral3D4N1210 = 1218,
        IntersectionQuadrilateral3D4N2110 = 1128,
        IntersectionQuadrilateral3D4N2310 = 1328,
        IntersectionQuadrilateral3D4N3210 = 1238,
        IntersectionQuadrilateral3D4N3010 = 1038,
        IntersectionQuadrilateral3D4N0310 = 1308,
        IntersectionQuadrilateral3D4N0112 = 21108,
        IntersectionQuadrilateral3D4N1012 = 21018,
        IntersectionQuadrilateral3D4N1212 = 21218,
        IntersectionQuadrilateral3D4N2112 = 21128,
        IntersectionQuadrilateral3D4N2312 = 21328,
        IntersectionQuadrilateral3D4N3212 = 21238,
        IntersectionQuadrilateral3D4N3012 = 21038,
        IntersectionQuadrilateral3D4N0312 = 21308,
        IntersectionQuadrilateral3D4N0121 = 12108,
        IntersectionQuadrilateral3D4N1021 = 12018,
        IntersectionQuadrilateral3D4N1221 = 12218,
        IntersectionQuadrilateral3D4N2121 = 12128,
        IntersectionQuadrilateral3D4N2321 = 12328,
        IntersectionQuadrilateral3D4N3221 = 12238,
        IntersectionQuadrilateral3D4N3021 = 12038,
        IntersectionQuadrilateral3D4N0321 = 12308,
        IntersectionQuadrilateral3D4N0123 = 32108,
        IntersectionQuadrilateral3D4N1023 = 32018,
        IntersectionQuadrilateral3D4N1223 = 32218,
        IntersectionQuadrilateral3D4N2123 = 32128,
        IntersectionQuadrilateral3D4N2323 = 32328,
        IntersectionQuadrilateral3D4N3223 = 32238,
        IntersectionQuadrilateral3D4N3023 = 32038,
        IntersectionQuadrilateral3D4N0323 = 32308,
        IntersectionQuadrilateral3D4N0132 = 23108,
        IntersectionQuadrilateral3D4N1032 = 23018,
        IntersectionQuadrilateral3D4N1232 = 23218,
        IntersectionQuadrilateral3D4N2132 = 23128,
        IntersectionQuadrilateral3D4N2332 = 23328,
        IntersectionQuadrilateral3D4N3232 = 23238,
        IntersectionQuadrilateral3D4N3032 = 23038,
        IntersectionQuadrilateral3D4N0332 = 23308,
        IntersectionQuadrilateral3D4N0130 = 3108,
        IntersectionQuadrilateral3D4N1030 = 3018,
        IntersectionQuadrilateral3D4N1230 = 3218,
        IntersectionQuadrilateral3D4N2130 = 3128,
        IntersectionQuadrilateral3D4N2330 = 3328,
        IntersectionQuadrilateral3D4N3230 = 3238,
        IntersectionQuadrilateral3D4N3030 = 3038,
        IntersectionQuadrilateral3D4N0330 = 3308,
        IntersectionQuadrilateral3D4N0103 = 30108,
        IntersectionQuadrilateral3D4N1003 = 30018,
        IntersectionQuadrilateral3D4N1203 = 30218,
        IntersectionQuadrilateral3D4N2103 = 30128,
        IntersectionQuadrilateral3D4N2303 = 30328,
        IntersectionQuadrilateral3D4N3203 = 30238,
        IntersectionQuadrilateral3D4N3003 = 30038,
        IntersectionQuadrilateral3D4N0303 = 30308
    };
    
#endif
    
///@}
///@name  Functions
///@{
    
///@}
///@name Kratos Classes
///@{
    
/** \brief MortarKinematicVariables
 * This is the definition of the kinematic variables used on the mortar operators assemble, which means three shape functions (one for the slave , one for the master and the third for the Lagrange Multipliers), and the jacobian in the corresponding Gauss point
 */

template< const unsigned int TNumNodes>
class MortarKinematicVariables
{
public:
    ///@name Type Definitions
    ///@{
        
    /// Counted pointer of MortarKinematicVariables
    KRATOS_CLASS_POINTER_DEFINITION( MortarKinematicVariables );
         
    ///@}
    ///@name Life Cycle
    ///@{

    MortarKinematicVariables(){}
    
    virtual ~MortarKinematicVariables(){}
    
    // Shape functions for contact pair
    Vector NMaster, NSlave, PhiLagrangeMultipliers;

    // Determinant of slave cell's jacobian
    double DetjSlave;
        
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /**
     * This method initialized the operators
     */
    void Initialize()
    {
        // Shape functions
        NMaster                = ZeroVector(TNumNodes);
        NSlave                 = ZeroVector(TNumNodes);
        PhiLagrangeMultipliers = ZeroVector(TNumNodes);
        
        // Jacobian of slave
        DetjSlave = 0.0;
    }
    
    /**
     * Print information about this object
     */ 
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NSlave:" << NSlave << std::endl;
        rOStream << "NMaster: " <<NMaster << std::endl;
        rOStream << "PhiLagrangeMultipliers: "<< PhiLagrangeMultipliers << std::endl;
        rOStream << "DetjSlave: " << DetjSlave << std::endl;
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

    ///@}

}; // Class MortarKinematicVariables

/** \brief MortarKinematicVariablesWithDerivatives
 * This class derives from MortarKinematicVariables and includes additionally to the variables of the previous class, the variables needed to define the directional derivatives of the mortar operators, like the gradients of the shape functions and the jacobians
 */
template< const unsigned int TDim, const unsigned int TNumNodes>
class MortarKinematicVariablesWithDerivatives : public MortarKinematicVariables<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{
    
    typedef MortarKinematicVariables<TNumNodes> BaseClassType;    
    
    /// Counted pointer of MortarKinematicVariables
    KRATOS_CLASS_POINTER_DEFINITION( MortarKinematicVariablesWithDerivatives );
         
    ///@}
    ///@name Life Cycle
    ///@{

    MortarKinematicVariablesWithDerivatives()= default;
    
    ~MortarKinematicVariablesWithDerivatives() override= default;
  
    // Shape functions local derivatives for contact pair
    Matrix DNDeMaster, DNDeSlave;
    
    /*
    * Jacobians in current configuration on all integration points of slave segment
    * Only those two variables contain info on all GP
    * other variables contain info only on the currently-calculated GP
    */
    Matrix jSlave;
    Matrix jMaster;
        
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /**
     * This method initialized the operators
     */
    void Initialize()
    {
        BaseClassType::Initialize();
        
        // Shape functions local derivatives
        DNDeMaster = ZeroMatrix(TNumNodes, TDim - 1);
        DNDeSlave  = ZeroMatrix(TNumNodes, TDim - 1);
        
        // Jacobians on all integration points
        jSlave  = ZeroMatrix(TDim, TDim - 1);
        jMaster = ZeroMatrix(TDim, TDim - 1);
    }
    
    /**
     * Print information about this object
     */ 
    void PrintInfo(std::ostream& rOStream) const
    {
        BaseClassType::PrintInfo(rOStream);
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

    ///@}

}; // Class MortarKinematicVariablesWithDerivatives

/** \brief DerivativeData
 * This data will be used to compute the derivatives.
 * This class includes different information that is used in order to compute the directional derivatives in the mortar contact conditions
 */
template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation>
class DerivativeData
{
public:
    ///@name Type Definitions
    ///@{
    
    // Auxiliar types
    typedef bounded_matrix<int, 1, 1> DummyBoundedMatrixType;
    
    typedef array_1d<double, TNumNodes> GeometryArrayType;
    
    typedef bounded_matrix<double, TNumNodes, TDim> GeometryDoFMatrixType;
    
    typedef bounded_matrix<double, TNumNodes, TNumNodes> GeometryMatrixType;
    
    typedef typename std::conditional<TNumNodes == 2, DummyBoundedMatrixType, bounded_matrix<double, 3, 3>>::type VertexDerivativesMatrixType;
    
    // Auxiliar sizes
    static const unsigned int DummySize = 1;
    
    static const unsigned int DoFSizeGeometry = (TNumNodes * TDim);
    
    static const unsigned int DoFSizePairedGeometry = 2 * (TNumNodes * TDim);
    
    static const unsigned int DoFSizeDerivativesDependence = (TDim == 2) ? DoFSizeGeometry : DoFSizePairedGeometry;
    
    static const unsigned int DoFSizeDerivativesVertex = (TDim == 2) ? DummySize : DoFSizePairedGeometry;
    
    ///@}
    ///@name Life Cycle
    ///@{

    DerivativeData()= default;
    
    virtual ~DerivativeData()= default;
    
    // The ALM parameters
    array_1d<double, TNumNodes> PenaltyParameter;
    double ScaleFactor;
    
    // The normals of the nodes
    GeometryDoFMatrixType NormalMaster, NormalSlave;
    
    // Displacements and original coordinates
    GeometryDoFMatrixType X1, X2, u1, u2;
    
    // Derivatives    
    array_1d<double, DoFSizeDerivativesDependence> DeltaDetjSlave;
    array_1d<GeometryArrayType, DoFSizeDerivativesDependence> DeltaPhi;
    array_1d<GeometryArrayType, DoFSizePairedGeometry> DeltaN1, DeltaN2;
    array_1d<GeometryDoFMatrixType, DoFSizeGeometry> DeltaNormalSlave, DeltaNormalMaster;
    array_1d<VertexDerivativesMatrixType, DoFSizeDerivativesVertex> DeltaCellVertex;
    
    /**
     * Ae matrix. 
     * This matrix corresponds with a coefficient matrix that relates the standard shape functions with the dual shape function 
     * It is used to stabilize and for an assembling reduction
     * The name is taken from Popp's work
     */
    GeometryMatrixType Ae;
    
    // Derivatives Ae
    array_1d<GeometryMatrixType, DoFSizeDerivativesDependence> DeltaAe;
    
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /**
     * Initializer method 
     * @param SlaveGeometry The geometry of the slave 
     * @param rCurrentProcessInfo The process info from the system
     */
    
    virtual void Initialize(
        const GeometryType& SlaveGeometry,
        const ProcessInfo& rCurrentProcessInfo
        )
    {        
        // The normals of the nodes
        NormalSlave = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry,  NORMAL, 0);
        
        // Displacements and velocities of the slave       
        u1 = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 0) - MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 1);
        X1 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(SlaveGeometry, false, 1);
        
        // We get the ALM variables
        for (unsigned int i = 0; i < TNumNodes; ++i)
            PenaltyParameter[i] = SlaveGeometry[i].GetValue(INITIAL_PENALTY);
        ScaleFactor = rCurrentProcessInfo[SCALE_FACTOR];
        
        // We initialize the derivatives
        for (unsigned int i = 0; i < TNumNodes * TDim; ++i) {
            DeltaDetjSlave[i] = 0.0;
            if (TDim == 3) DeltaDetjSlave[i + TNumNodes * TDim] = 0.0;
            DeltaPhi[i] = ZeroVector(TNumNodes);
            if (TDim == 3) DeltaPhi[i + TNumNodes * TDim] = ZeroVector(TNumNodes);
            DeltaN1[i] = ZeroVector(TNumNodes);
            DeltaN1[i + TNumNodes * TDim] = ZeroVector(TNumNodes);
            DeltaN2[i] = ZeroVector(TNumNodes);
            DeltaN2[i + TNumNodes * TDim] = ZeroVector(TNumNodes);
        }
    }
    
    /**
     * This method reset to zero the cell vertex derivatives
     */
    
    virtual void ResetDerivatives()
    {                
        // Derivatives 
        if (TDim == 3) { // Derivative of the cell vertex
            for (unsigned int i = 0; i < TNumNodes * TDim; ++i) {
                DeltaCellVertex[i] = ZeroMatrix(3, 3);
                DeltaCellVertex[i + TNumNodes * TDim] = ZeroMatrix(3, 3);
            }
        }
    }
    
    /**
     * Initialize the DeltaAe components
     */
    
    void InitializeDeltaAeComponents()
    {
        // Ae
        Ae = ZeroMatrix(TNumNodes, TNumNodes);
        
        // Derivatives Ae
        for (unsigned int i = 0; i < DoFSizeDerivativesDependence; ++i)
            DeltaAe[i] = ZeroMatrix(TNumNodes, TNumNodes);
    }

    /**
     * Updating the Master pair
     * @param MasterGeometry The master geometry
     */
    
    virtual void UpdateMasterPair(const GeometryType& MasterGeometry)
    {        
        NormalMaster = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry,  NORMAL, 0);
        
        // Displacements, coordinates and normals of the master
        u2 = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry, DISPLACEMENT, 0)
           - MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry, DISPLACEMENT, 1);
        X2 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(MasterGeometry, false, 1);
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

    ///@}
    
};  // Class DerivativeData

/** \brief DerivativeDataFrictional
 * This class is a derived class of DerivativeData. Includes additionally the derivatives necessary to compute the directional derivatives for the frictional conditions
 */
template< unsigned int TDim, unsigned int TNumNodes, bool TNormalVariation>
class DerivativeDataFrictional : public DerivativeData<TDim, TNumNodes, TNormalVariation>
{
public:
    ///@name Type Definitions
    ///@{
    
    // Auxiliar types
    typedef DerivativeData<TDim, TNumNodes, TNormalVariation> BaseClassType;
    
    typedef bounded_matrix<double, TNumNodes, TDim> GeometryDoFMatrixType;
    
    // Auxiliar sizes
    static const unsigned int DoFSizeGeometry = (TNumNodes * TDim);
    
    static const unsigned int DoFSizePairedGeometry = 2 * (TNumNodes * TDim);
    
    ///@}
    ///@name Life Cycle
    ///@{

    DerivativeDataFrictional()= default;
    
    virtual ~DerivativeDataFrictional()= default;
    
    // The ALM parameters
    double TangentFactor;
    
    // Displacements and velocities
    GeometryDoFMatrixType u1old, u2old;
    
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /**
     * Initializer method 
     * @param SlaveGeometry The geometry of the slave 
     * @param rCurrentProcessInfo The process info from the system
     */
    
    void Initialize(
        const GeometryType& SlaveGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {        
        BaseClassType::Initialize(SlaveGeometry, rCurrentProcessInfo);
        
        TangentFactor = rCurrentProcessInfo[TANGENT_FACTOR];
        
        u1old = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 1) - MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 2);
    }
    
    /**
     * Updating the Master pair
     * @param MasterGeometry The geometry of the master
     */
    
    void UpdateMasterPair(const GeometryType& MasterGeometry) override
    {
        BaseClassType::UpdateMasterPair(MasterGeometry);
        
        u2old = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry, DISPLACEMENT, 1) - MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry, DISPLACEMENT, 2);
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

    ///@}
    
};  // Class DerivativeDataFrictional

/** \brief MortarOperator
 * This is the definition of the mortar operator according to the work of Alexander Popp: https://www.lnm.mw.tum.de/staff/alexander-popp/
 * In particular the thesis of contact mechanics based in mortar method available at: https://mediatum.ub.tum.de/?id=1109994
 * These mortar operator are assembled as mass matrices in order to transfer information from the slave side to the master side. 
 * The operators are DOperator for the slave side and MOperator for master respectively.
 * In order to compute these operators, the shape functions of both domains are necessary (using the slave side as reference), as well as the integration weight and jacobian in the integration point. 
 * Popp thesis pge 50 and following
 */
template< const unsigned int TNumNodes>
class MortarOperator
{
public:
    ///@name Type Definitions
    ///@{
        
    typedef MortarKinematicVariables<TNumNodes> KinematicVariables;
    
    // Auxiliar types
    typedef bounded_matrix<double, TNumNodes, TNumNodes> GeometryMatrixType;
    
    /// Counted pointer of MortarOperator
    KRATOS_CLASS_POINTER_DEFINITION( MortarOperator );
         
    ///@}
    ///@name Life Cycle
    ///@{

    MortarOperator(){}
    
    virtual ~MortarOperator(){}
    
    // Mortar condition matrices - DOperator and MOperator
    GeometryMatrixType DOperator, MOperator;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /**
     * This method initialized the operators
     */
    void Initialize()
    {
        // We initialize the D and M operators
        DOperator = ZeroMatrix(TNumNodes, TNumNodes);
        MOperator = ZeroMatrix(TNumNodes, TNumNodes);
    }
    
    /**
     * It calculates the mortar operators. Popp thesis page 56, equation 3.31 and 3.32
     * @param rKinematicVariables Corresponds with the kinematic variables
     * @param rIntegrationWeight The corresponding integration weight
     */
    void CalculateMortarOperators(
        KinematicVariables& rKinematicVariables,
        const double rIntegrationWeight
        )
    {
        /* DEFINITIONS */
        const double det_j_slave = rKinematicVariables.DetjSlave; 
        const Vector& phi_vector  = rKinematicVariables.PhiLagrangeMultipliers;
        const Vector& n1_vector   = rKinematicVariables.NSlave;
        const Vector& n2_vector   = rKinematicVariables.NMaster;
        
        for (unsigned int i_slave = 0; i_slave < TNumNodes; ++i_slave) {
            for (unsigned int j_slave = 0; j_slave < TNumNodes; ++j_slave) {
                const double phi = phi_vector[i_slave];
                
                DOperator(i_slave, j_slave) += det_j_slave * rIntegrationWeight * phi * n1_vector[j_slave];
                MOperator(i_slave, j_slave) += det_j_slave * rIntegrationWeight * phi * n2_vector[j_slave];
            }
        }
    }
    
    /**
     * It calculates the POperator (Inverse(D x M))
     * Popp thesis page 83 equation 3.88
     */
    GeometryMatrixType ComputePOperator()
    {
        // We calculate the inverse of D operator
        double auxdet;
        const GeometryMatrixType& inv_D_operator = MathUtils<double>::InvertMatrix<TNumNodes>(DOperator, auxdet);
        
        // We calculate the P operator
        const GeometryMatrixType POperator = prod(inv_D_operator, MOperator);
        
        return POperator;
    }
    
    /**
     * Print information about this object
     */ 
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DOperator: " << DOperator << std::endl;
        rOStream << "MOperator: " << MOperator << std::endl;
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

    ///@}

}; // Class MortarOperatorWithDerivatives

/** \brief MortarOperatorWithDerivatives
 * This class derives from the MortarOperator class and it includes the derived operators. 
 * The derived operators are defined in each DoF of each domain, which means TNumNodes x TDim x 2 derivatives definitions in order to compute all the necessary derivatives
 * Popp thesis page 102 and following 
 */
template< const unsigned int TDim, const unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
class MortarOperatorWithDerivatives : public MortarOperator<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{
        
    typedef MortarOperator<TNumNodes>                                                  BaseClassType;
    
    typedef MortarKinematicVariables<TNumNodes>                                   KinematicVariables;
    
    typedef DerivativeDataFrictional<TDim, TNumNodes, TNormalVariation> DerivativeDataFrictionalType;
    
    typedef DerivativeData<TDim, TNumNodes, TNormalVariation>        DerivativeFrictionalessDataType;
    
    typedef typename std::conditional<TFrictional, DerivativeDataFrictionalType, DerivativeFrictionalessDataType>::type DerivativeDataType;
    
    // Auxiliar types
    typedef bounded_matrix<double, TNumNodes, TNumNodes> GeometryMatrixType;
        
    // Auxiliar sizes
    static const unsigned int DoFSizeGeometry = (TNumNodes * TDim);
    
    static const unsigned int DoFSizePairedGeometry = 2 * (TNumNodes * TDim);
    
    static const unsigned int DoFSizeDerivativesDependence = (TDim == 2) ? DoFSizeGeometry : DoFSizePairedGeometry;
    
    /// Counted pointer of MortarOperatorWithDerivatives
    KRATOS_CLASS_POINTER_DEFINITION( MortarOperatorWithDerivatives );
    
    ///@}
    ///@name Life Cycle
    ///@{

    MortarOperatorWithDerivatives(){}
    
    ~MortarOperatorWithDerivatives() override{}
    
    // D and M directional derivatives
    array_1d<GeometryMatrixType, DoFSizePairedGeometry> DeltaDOperator, DeltaMOperator;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /**
     * This method initialized the operators
     */
    void Initialize()
    {
        BaseClassType::Initialize();
        
        // We initialize the D and M derivatives operators 
        for (unsigned int i = 0; i < TNumNodes * TDim; ++i)
        {
            DeltaDOperator[i] = ZeroMatrix(TNumNodes, TNumNodes);
            DeltaDOperator[i + TNumNodes * TDim] = ZeroMatrix(TNumNodes, TNumNodes);
            DeltaMOperator[i] = ZeroMatrix(TNumNodes, TNumNodes);
            DeltaMOperator[i + TNumNodes * TDim] = ZeroMatrix(TNumNodes, TNumNodes);
        }
    }
    
    /**
     * It calculates the mortar operators
     * Popp thesis page 102 equation equation 4.32 and 4.33 / 4.37 and 4.38
     * @param rKinematicVariables Corresponds with the kinematic variables
     * @param rIntegrationWeight The corresponding integration weight
     */
    void CalculateDeltaMortarOperators(
        KinematicVariables& rKinematicVariables,
        DerivativeDataType& rDerivativeData,
        const double rIntegrationWeight
        )
    {
        /* DEFINITIONS */
        const double det_j_slave = rKinematicVariables.DetjSlave; 
        const Vector& vector_phi = rKinematicVariables.PhiLagrangeMultipliers;
        const Vector& vector_n1  = rKinematicVariables.NSlave;
        const Vector& vector_n2  = rKinematicVariables.NMaster;
        
        // Derivatives
        const array_1d<double, DoFSizeDerivativesDependence>& delta_det_j_slave = rDerivativeData.DeltaDetjSlave;
        const array_1d<array_1d<double, TNumNodes >, DoFSizeDerivativesDependence>& delta_phi = rDerivativeData.DeltaPhi;
        const array_1d<array_1d<double, TNumNodes >, DoFSizePairedGeometry>& delta_n1  = rDerivativeData.DeltaN1;
        const array_1d<array_1d<double, TNumNodes >, DoFSizePairedGeometry>& delta_n2  = rDerivativeData.DeltaN2;
        
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
            const double phi = vector_phi[i_node];
            
            for (unsigned int j_node = 0; j_node < TNumNodes; ++j_node) {
                const double n1 = vector_n1[j_node];
                const double n2 = vector_n2[j_node];
                
                BaseClassType::DOperator(i_node, j_node) += det_j_slave * rIntegrationWeight * phi * n1;
                BaseClassType::MOperator(i_node, j_node) += det_j_slave * rIntegrationWeight * phi * n2;
                
                for (unsigned int i = 0; i < TDim * TNumNodes; ++i) {
                    DeltaDOperator[i](i_node, j_node) += delta_det_j_slave[i] * rIntegrationWeight * phi* n1        
                                                       + det_j_slave * rIntegrationWeight * delta_phi[i][i_node] * n1
                                                       + det_j_slave * rIntegrationWeight * phi* delta_n1[i][j_node];
                                                                                
                    DeltaMOperator[i](i_node, j_node) += delta_det_j_slave[i] * rIntegrationWeight * phi* n2        
                                                       + det_j_slave * rIntegrationWeight * delta_phi[i][i_node] * n2
                                                       + det_j_slave * rIntegrationWeight * phi* delta_n2[i][j_node];
                }
                for (unsigned int i = TDim * TNumNodes; i < 2 * TDim * TNumNodes; ++i) {
                    DeltaDOperator[i](i_node, j_node) += det_j_slave * rIntegrationWeight * phi * delta_n1[i][j_node];
                                                                                
                    DeltaMOperator[i](i_node, j_node) += det_j_slave * rIntegrationWeight * phi * delta_n2[i][j_node];
                }
                if (TDim == 3) {
                    for (unsigned int i = TDim * TNumNodes; i < 2 * TDim * TNumNodes; ++i) {
                        DeltaDOperator[i](i_node, j_node) += delta_det_j_slave[i] * rIntegrationWeight * phi * n1;
                        DeltaDOperator[i](i_node, j_node) += det_j_slave * rIntegrationWeight * delta_phi[i][i_node] * n1;
                                                                                    
                        DeltaMOperator[i](i_node, j_node) += delta_det_j_slave[i] * rIntegrationWeight * phi * n2;
                        DeltaMOperator[i](i_node, j_node) += det_j_slave * rIntegrationWeight * delta_phi[i][i_node] * n2;
                    }
                }
            }
        }
    }
    
    /**
     * Print information about this object
     */ 
    void PrintInfo(std::ostream& rOStream) const
    {
        BaseClassType::PrintInfo(rOStream);
        
        for (unsigned int i = 0; i < TNumNodes * TDim; ++i)
        {
            rOStream << "DeltaDOperator_" << i << ": " << DeltaDOperator[i] << std::endl;
            rOStream << "DeltaMOperator_" << i << ": " << DeltaMOperator[i] << std::endl;
        }
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

    ///@}

}; // Class MortarOperatorWithDerivatives

/** \brief DualLagrangeMultiplierOperators
 * This is the definition dual lagrange multiplier operators according to the work of Alexander Popp: https://www.lnm.mw.tum.de/staff/alexander-popp/
 * In particular the thesis of contact mechanics based in mortar method available at: https://mediatum.ub.tum.de/?id=1109994
 * In order to compute the dual LM shape function the Ae operator must be computed, which depends of the Me and De operators. Phi = Ae * NSlave.  In a similar way to the mortar operators, the De corresponds with a diagonal operator and Me with a sparse operator respectively. Ae = De * inv(Me) 
 * Popp thesis page 69 and following
 */

template< const unsigned int TNumNodes>
class DualLagrangeMultiplierOperators
{
public:
    ///@name Type Definitions
    ///@{
        
    typedef MortarKinematicVariables<TNumNodes> KinematicVariables;
    
    // Auxiliar types
    typedef bounded_matrix<double, TNumNodes, TNumNodes> GeometryMatrixType;
    
    /// Counted pointer of DualLagrangeMultiplierOperators
    KRATOS_CLASS_POINTER_DEFINITION( DualLagrangeMultiplierOperators );
         
    ///@}
    ///@name Life Cycle
    ///@{

    DualLagrangeMultiplierOperators(){}
    
    virtual ~DualLagrangeMultiplierOperators(){}
    
    GeometryMatrixType Me, De;
        
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /**
     * This method initialized the operators
     */
    void Initialize()
    {
        // We initialize the De and Me operators
        Me = ZeroMatrix(TNumNodes, TNumNodes);
        De = ZeroMatrix(TNumNodes, TNumNodes);
    }
    
    /**
     * Calculates the Ae components necessary to compute the Phi_LagrangeMultipliers shape functions. For that it integrates De and Me. Popp thesis page 70 eq. 3.65
     * @param rKinematicVariables The kinematic variables
     * @param rIntegrationWeight The integration weight considered
     */
    void CalculateAeComponents(
        KinematicVariables& rKinematicVariables,
        const double rIntegrationWeight
        )
    {
        /* DEFINITIONS */
        const Vector& n1 = rKinematicVariables.NSlave;
        const double det_j = rKinematicVariables.DetjSlave; 
        
        De += rIntegrationWeight * (ComputeDe(n1, det_j));
        Me += rIntegrationWeight * det_j * outer_prod(n1, n1);
    }
    
    /**
     * Calculates the matrix Ae. To avoid problems in the inversion the matrix is normalized
     * Popp thesis page 70. Equation 3.65
     */
    bool CalculateAe(GeometryMatrixType& Ae)
    {        
        const double tolerance = std::numeric_limits<double>::epsilon(); 
        
        // We compute the norm
        const double norm_me = norm_frobenius(Me);
        
        // Now we normalize the matrix
        const GeometryMatrixType normalized_Me = Me/norm_me;
        
        // We compute the normalized inverse 
        double aux_det = MathUtils<double>::DetMat<GeometryMatrixType>(normalized_Me); 
        if (std::abs(aux_det) >= tolerance) { 
            const GeometryMatrixType& normalized_inv_Me = MathUtils<double>::InvertMatrix<TNumNodes>(normalized_Me, aux_det, tolerance);  
             
            noalias(Ae) = (1.0/norm_me) * prod(De, normalized_inv_Me); 
            return true;
        } 
    #ifdef KRATOS_DEBUG
        else {
            std::cout << "WARNING:: Me matrix can not bee inverted. Determinant: " << aux_det << std::endl;
            KRATOS_WATCH(normalized_Me);
        }
    #endif
        
        noalias(Ae) = IdentityMatrix(TNumNodes);  
        return false;
    }   
    
    /**
     * Calculates the matrix De
     * @param N1 The shape function 
     * @param detJ The jacobian of the geometry 
     */
    GeometryMatrixType ComputeDe(        
        const Vector& N1, 
        const double detJ 
        ) const 
    {
        GeometryMatrixType De;
    
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            for (unsigned int j = 0; j < TNumNodes; ++j) {
                if (i == j) De(i,i) = detJ * N1[i];
                else De(i,j) = 0.0;
            }
        }
        
        return De;
    }
    
    /**
     * Print information about this object
     */ 
    void PrintInfo(std::ostream& rOStream) const {
        rOStream << "Me: " << Me << std::endl;
        rOStream << "De: " << De << std::endl;
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

    ///@}

}; // Class DualLagrangeMultiplierOperators

/** \brief DualLagrangeMultiplierOperatorsWithDerivatives
 * This is the definition dual lagrange multiplier operators including the derivatives. Is based in the same work as the previous class. In this case it computes the derivatives in order to compute the directionald erivative of the dual shape functions
 * Popp thesis page 111 and following
 */

template< const unsigned int TDim, const unsigned int TNumNodes, bool TFrictional, bool TNormalVariation>
class DualLagrangeMultiplierOperatorsWithDerivatives : public DualLagrangeMultiplierOperators<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{
        
    typedef DualLagrangeMultiplierOperators<TNumNodes>                                 BaseClassType;  
    
    typedef MortarKinematicVariablesWithDerivatives<TDim, TNumNodes>          KinematicVariablesType;
    
    typedef DerivativeDataFrictional<TDim, TNumNodes, TNormalVariation> DerivativeDataFrictionalType;
    
    typedef DerivativeData<TDim, TNumNodes, TNormalVariation>        DerivativeFrictionalessDataType;
    
    typedef typename std::conditional<TFrictional, DerivativeDataFrictionalType, DerivativeFrictionalessDataType>::type DerivativeDataType;
    
    // Auxiliar types
    typedef bounded_matrix<double, TNumNodes, TNumNodes> GeometryMatrixType;
    
    // Auxiliar sizes
    static const unsigned int DoFSizeGeometry = (TNumNodes * TDim);
    
    static const unsigned int DoFSizePairedGeometry = 2 * (TNumNodes * TDim);
    
    static const unsigned int DoFSizeDerivativesDependence = (TDim == 2) ? DoFSizeGeometry : DoFSizePairedGeometry;
    
    /// Counted pointer of DualLagrangeMultiplierOperatorsWithDerivatives
    KRATOS_CLASS_POINTER_DEFINITION( DualLagrangeMultiplierOperatorsWithDerivatives );
         
    ///@}
    ///@name Life Cycle
    ///@{

    DualLagrangeMultiplierOperatorsWithDerivatives(){}
    
    ~DualLagrangeMultiplierOperatorsWithDerivatives() override{}
    
    // Derivatives matrices
    array_1d<GeometryMatrixType, DoFSizeDerivativesDependence> DeltaMe;
    array_1d<GeometryMatrixType, DoFSizeDerivativesDependence> DeltaDe;
        
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /**
     * This method initialized the operators
     */
    void Initialize()
    {
        BaseClassType::Initialize();
        
        // Derivatives matrices
        for (unsigned int i = 0; i < DoFSizeDerivativesDependence; ++i) {
            DeltaMe[i] = ZeroMatrix(TNumNodes, TNumNodes);
            DeltaDe[i] = ZeroMatrix(TNumNodes, TNumNodes);
        }
    }
    
    /**
     * Calculates the Ae components and its derivatives necessary to compute the Phi_LagrangeMultipliers shape functions. Popp thesis page 112 eq. 4.59
     * @param rKinematicVariables The kinematic variables
     * @param rIntegrationWeight The integration weight considered
     */
    void CalculateDeltaAeComponents(
        KinematicVariablesType& rKinematicVariables,
        DerivativeDataType& rDerivativeData,
        const double rIntegrationWeight
        )
    {
        /* DEFINITIONS */
        const double det_j_slave = rKinematicVariables.DetjSlave;
        const Vector& n1 = rKinematicVariables.NSlave;
        
        BaseClassType::CalculateAeComponents(rKinematicVariables, rIntegrationWeight);
        
        for (unsigned int i = 0; i < DoFSizeDerivativesDependence; ++i) {
            const double delta_det_j = rDerivativeData.DeltaDetjSlave[i];
            const Vector& delta_n1 = rDerivativeData.DeltaN1[i];
            
            DeltaDe[i] += rIntegrationWeight * this->ComputeDe( n1, delta_det_j )
                       +  rIntegrationWeight * this->ComputeDe( delta_n1, det_j_slave );
            
            DeltaMe[i] += rIntegrationWeight * delta_det_j * outer_prod(n1, n1)
                       +  rIntegrationWeight * det_j_slave * (outer_prod(delta_n1, n1) + outer_prod(n1, delta_n1));
        }
    }
 
    /**
     * Calculates the matrix DeltaAe. Popp thesis page 112 equation 4.58
     */
    bool CalculateDeltaAe(DerivativeDataType& rDerivativeData)
    {        
        double aux_det;
        const double tolerance = std::numeric_limits<double>::epsilon();
        
        // We compute the norm
        const double norm_Me = norm_frobenius(BaseClassType::Me);
        
        // Now we normalize the matrix
        const GeometryMatrixType normalized_Me = BaseClassType::Me/norm_Me;
        
        // We compute the normalized inverse
        aux_det = MathUtils<double>::DetMat<GeometryMatrixType>(normalized_Me);
        if (std::abs(aux_det) < tolerance) return false;
        
        const GeometryMatrixType normalized_inv_Me = MathUtils<double>::InvertMatrix<TNumNodes>(normalized_Me, aux_det, tolerance); 
        
        // Now we compute the inverse
        const GeometryMatrixType inv_Me = normalized_inv_Me/norm_Me;
        
        noalias(rDerivativeData.Ae) = prod(BaseClassType::De, inv_Me);
        
        array_1d<GeometryMatrixType , DoFSizeDerivativesDependence>& delta_Ae = rDerivativeData.DeltaAe;
        
        for (unsigned int i = 0; i < DoFSizeDerivativesDependence; ++i) {            
            const GeometryMatrixType aux_matrix = DeltaDe[i] - prod(rDerivativeData.Ae, DeltaMe[i]);
            noalias(delta_Ae[i]) = prod(aux_matrix, inv_Me);
        }
        
        return true;
    }
    
    /**
     * This method prints the current operators
     */
    void PrintInfo(std::ostream& rOStream) const
    {
        BaseClassType::PrintInfo(rOStream);
        
        // Derivatives matrices
        for (unsigned int i = 0; i < DoFSizeDerivativesDependence; ++i) {
            rOStream << "DeltaMe_" << i << ": " << DeltaMe[i] << std::endl;
            rOStream << "DeltaDe_" << i << ": " << DeltaDe[i] << std::endl;
        }
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

    ///@}

}; // Class DualLagrangeMultiplierOperatorsWithDerivatives

/** @brief Custom Point container to be used by the mapper
 * This point which is a derived class of the standard point, contains the variable mBelongs. This variable is a "hash" that can be used to determine where in which intersections the point belongs
 */

template<unsigned int TNumNodes>
class PointBelong : public Point
{
public:
    ///@name Type Definitions
    ///@{
    
    typedef typename std::conditional<TNumNodes == 2, PointBelongsLine2D2N, typename std::conditional<TNumNodes == 3, PointBelongsTriangle3D3N, PointBelongsQuadrilateral3D4N>::type>::type BelongType;
    
    /// Counted pointer of PointBelong
    KRATOS_CLASS_POINTER_DEFINITION( PointBelong );

    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Default constructors
    PointBelong():
        Point()
    {}

    PointBelong(const array_1d<double, 3> Coords):
        Point(Coords)
    {}
    
    PointBelong(const array_1d<double, 3> Coords, const BelongType& ThisBelongs):
        Point(Coords),
        mBelongs(ThisBelongs)
    {}
    
    /// Destructor.
    ~PointBelong() override= default;
    
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This method allows to set where the point belongs
     */
    void SetBelong(BelongType ThisBelongs)
    {
        mBelongs = ThisBelongs;
    }
    
    /**
     * This method recovers where the point belongs
     */
    BelongType GetBelong() const
    {
        return mBelongs;
    }
    
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

    BelongType mBelongs; // To know if the point belongs to the master/slave/intersection (just 3D) side          

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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class PointBelong 

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_MORTAR_CLASSES  defined 
