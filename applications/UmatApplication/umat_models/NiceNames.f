c------------------------------------------------------ 
c------------------------------------------------------ 
C     William Fuentes, August 2011
C     w-fuente@uniandes.edu.co
C     KIT institut-IBFA
C     Karlsruhe, Germany
c------------------------------------------------------ 
c------------------------------------------------------ 
C     Library of structure variables.
c------------------------------------------------------ 
c------------------------------------------------------  
      MODULE nice_names
      implicit none
      save                                                              
      type MATERIALCONSTANTS                                        
      double precision ::
     1         K0,                                                    
     2         G0,                                           
     3         rK,                                            
     4         ec0,                                         
     5         hs,                                        
     6         n,                                                   
     7         phic,                                                   
     8         nd,                                                      
     9         Ad,                                                  
     1         Ydegree0,      
     2         nY, 
     3         ne, 
     4         hk, 
     5         mfactor, 
     6         nexp , 
     7         Kwater=0.0d0
      INTEGER :: surf=0   !=0, eight-curve equation, =1 distorted lemniscate
     8          
      end type MATERIALCONSTANTS
c------------------------------------------------------ 
c------------------------------------------------------       
      type STATEVARIABLES                                              
      double precision ::
     1        void, 
     2        pc,
     3        obs1,
     4        obs2,
     5        obs3,
     6        obs4,
     7        obs5,
     8        pcmin=0.0d0,        !Minimum preconsolidation pressure allowed 
     9        devEps
      double precision,dimension (1:3,1:3) ::                      
     1 alpha33
	  INTEGER :: maxnintS
	  
      end type STATEVARIABLES   
c------------------------------------------------------ 
c------------------------------------------------------       
      type ELEMENTVARIABLES                                          
      INTEGER ::
     1        NTENS, 
     2        NDI,
     3        NSHR,
     4        KSTEP,
     5        KINC
       double precision ::
     1        dtime
      double precision,dimension (2) :: 
     1        time
      end type ELEMENTVARIABLES   
c------------------------------------------------------ 
c------------------------------------------------------       
      type INTEGRATIONVARIABLES                                          
      DOUBLE PRECISION ::
     1        tolsubt=5.0d-5,    !Substepping tolerance
     2        maxdef=1.0d-2     !Maximum deformation allowed, if not, increment rejected
       INTEGER ::
     1        maxnint=10000,       !Maximum number of substepping
     2        scheme=1,          !1==substepping, 2==Sloan, 3==Roddeman
     3        maxiter=1000,      !Maximum number iterations
     4        errormsg=0,        !error message
     5        kstep,             !kstep
     6        jstep,             !jstep
     7        dtsub,              !dtsub
     8        dtsub2              !dtsub     
      end type INTEGRATIONVARIABLES 
c------------------------------------------------------ 
c------------------------------------------------------       
      type LoadDef  
      character*80 ::
     1 Elset,              !element set where the load is applying
     2 Loadname            !Name that identifies the load type                           
      DOUBLE PRECISION ::
     1        magnitude    !magnitude from the Load
      end type LoadDef   
c------------------------------------------------------       
      type UELDef         !user element definition
      character*80 ::
     1 USER, ELEMENT,     !Command name to save *user element        
     2 NODES,             !number of nodes per element
     3 TYPEU,             !name of the element for library
     4 PROPERTIES,        !number of parameters
     5 COORDINATES,       !number of the dimension space 2D=2
     6 VARIABLES          !number of state variables              
      end type UELDef        
c------------------------------------------------------
c------------------------------------------------------
c     Element nodes     
      type Element        !
      character*80 ::
     1 TYPEU            !ElementType       
     2 elset            !Name of the set of elements elset          
      end type Element        
c------------------------------------------------------ 
c------------------------------------------------------
c     NSET     
      type Nset        !
c     Structure that saves the name and the number of the Nset      
      character*10 ::
     1 NsetName            !Name of the Nset       
      integer, allocatable :: Nodes(:)       
      end type Nset        
c------------------------------------------------------ 
c------------------------------------------------------  
c     ELSET    
      type Elset        !
c     Structure that saves the name and the number of the Elset      
      character*10 ::
     1 ElsetName            !Name of the Elset       
      integer, allocatable :: Elements(:)       
      end type Elset        
c------------------------------------------------------ 
c------------------------------------------------------ 
c------------------------------------------------------  
c     Boundary
c     Saves the boundaries conditions   
      type Boundary       !
c     Structure that saves the boundaries conditions   
      character*10 ::
     1 NSetName            !Name of the Nset      
      integer  DOFs(2)     !Degree of freedoms 
      double precision value !value of the boundary condition   
      end type Boundary        
c------------------------------------------------------  
c     Modelchange
      type Modelchange       !
c     Structure that saves the boundaries conditions   
      character*10 ::
     1 ElSetName            !Name of the Elset      
      integer  Remove       !=0 for removing, =1 for Adding
      integer  Step         !step number
      integer, allocatable ::  NodesOut(:) , AllNodes(:) 
      end type Modelchange         
c------------------------------------------------------                
      end MODULE nice_names
