 -----------------------------------------------------
 BOUNDARY MAPPING
 -----------------------------------------------------

            Markers: 3
               Type: symmetry plane
               Name: SYM1
 block end
 ---------------------------

            Markers: 1
               Type: farfield
        Angle alpha (degree): 6
         Angle beta: 0
               Name: INLET
 block end
 ---------------------------

            Markers: 2
               Type: exit-pressure outflow
               Name: OUTLET
 block end
 ---------------------------

            Markers: 13
               Type: viscous wall
               Name: MEMBRANE
            Subtype: turbulent
      Monitor forces (0/1): 1
   Write surface data (0/1): 1
 block end
 ---------------------------

            Markers: 4
               Type: symmetry plane
               Name: SYM2
 block end
###########################################################################################
                                                                  : -
 Geometry --------------------------------------------------------: -
                                                                  : -
                                                        Grid scale: 1
                                           Reference relation area: 0.00045
                              Reference length (pitching momentum): 0.002
                        Reference length (rolling/yawing momentum): 0.225
                                               Origin coordinate x: 0
                                               Origin coordinate y: 0
                                               Origin coordinate z: 0
                                                                  : -
 References ------------------------------------------------------: -
                                                                  : -
                                             Reference temperature: 293.15
                                                   Reynolds length: 0.225
                                             Reference Mach number: 0.11658
                                                                  : -
 Perfect gas thermodynamic ---------------------------------------: -
                                                                  : -
                                                    Gas constant R: 287.06
                                                Gas constant gamma: 1.4
                                                                  : -
 Transport coefficients ------------------------------------------: -
                                                                  : -
                                                   Reynolds number: 613699.8287
                                                    Prandtl number: 0.72
                                               Sutherland constant: 110.4
                                    Sutherland reference viscosity: 1.716e-05
                                  Sutherland reference temperature: 273
                                                                  : -
 Flowfield initialization ----------------------------------------: -
                                                                  : -
                                       Init total conditions (0/1): 0
                                                                  : -
###########################################################################################
                                                                  : -
 Files/IO --------------------------------------------------------: -
                                                                  : -
                                         Boundary mapping filename: (thisfile)
                                             Primary grid filename: /home/inigo/software/kratosMerge/Kratos/applications/CoSimulationApplication/tests/TAU/Mesh/airfoil_Structured_scaliert.grid
                                                                  : -
                                               Output files prefix: /home/inigo/software/kratosMerge/Kratos/applications/CoSimulationApplication/tests/TAU/Outputs/airfoilSol
                                               Restart-data prefix: (none)
                                                                  : -
                                   Surface output description file: (thisfile)
                                     Field output description file: (thisfile)
                                   Profile output description file: (thisfile)
                                     Plane output description file: (thisfile)
                                Transition module description file: (thisfile)
                                                                  : -
 Universal -------------------------------------------------------: -
                                                                  : -
                                                       Solver type: Flow
                                                                  : -
 Controls --------------------------------------------------------: -
                                                                  : -
                                                      Output level: 25
                Reference system of forces and moments(tau/ln9300): tau
                               Write pointdata dimensionless (0/1): 0
                                                                  : -
                                  Automatic parameter update (0/1): 1
                             Automatic parameter update mode (0/1): 0
                                       Accumulate queue time (0/1): 0
                                                                  : -
 Memory management -----------------------------------------------: -
                                                                  : -
                                             Increase memory (0/1): 1
                                                                  : -
                                                                  : -
###########################################################################################
                                                                  : -
 Partitioning ----------------------------------------------------: -
                                                                  : -
                                       Type of partitioning (name): zoltan
                            Use parallel initial partitioner (0/1): 0
                                                                  : -
 #MPI--------------------------------------------------------------: -
 #                                      MPI-rank optimization (0/1): 1
 #                                     Number of CPU-cores per node: 28
                                                                  : -
###########################################################################################
                                                                  : -
 Preprocessing ---------------------------------------------------: -
                                                                  : -
                                        Number of multigrid levels: 3
                                                       Grid metric: Cell_Vertex
                                                                  : -
 Runtime optimisation --------------------------------------------: -
                                                                  : -
                             Cache-coloring (0/max_faces in color): 100000
                                      Bandwidth optimisation (0/1): 1
                                2D offset vector (0 / x=1,y=2,z=3): 0
                                       Compute lusgs mapping (0/1): 1
                                                                  : -
 Agglomeration ---------------------------------------------------: -
                                                                  : -
                                             Type of agglomeration: private
                                             Agglomeration version: 0
                                               Point fusing reward: 1.2
                                        Structured grid coarsening: 0
                                                                  : -
###########################################################################################
                                                                  : -
 Multigrid -------------------------------------------------------: -
                                                                  : -
				  	   MG description filename: 3w
                                         Multigrid indicator (0/1): 0
                                     SG start up steps (fine grid): 0
                          Turbulence equations use multigrid (0/1): 1
 		Skip dissipation update for forcing function (0/1): 0
								  : -
#     					 Number of MG-cycle levels: 3						# Settings for Improved Multigrid Scheme
#		   Number of relaxation steps between MG-transfers: 3						# Settings for Improved Multigrid Scheme
#	Use source terms in coarse grid turbulence equations (0/1): 1						# Settings for Improved Multigrid Scheme
#                Use face gradient correction on coarse grid (0/1): 1						# Settings for Improved Multigrid Scheme
#                   Set boundary fluxes for forcing function (0/1): 1       					# Settings for Improved Multigrid Scheme
                                                                  : -
 Full multigrid --------------------------------------------------: -
                                                                  : -
                                             Multigrid start level: 1
                           Maximal time step number (coarse grids): 100
                                   Minimum residual (coarse grids): 0.001
                   Full multigrid central scheme first-order (0/1): 1
                                                                  : -
###########################################################################################
                                                                  : -
 Timestepping Start/Stop -----------------------------------------: -
                                                                  : -
                                          Maximal time step number: 200
					          Minimum residual: 1e-10
							          : -
 Output options --------------------------------------------------: -
                                                                  : -
                                                     Output period: 200
                                                     Field output values: v_p_rho_cp_Ptot_mach_vort_Nk_heli_Hn_l2_Q_velgrad
                                                                  : -
                                             Surface output period: 200
                                             Surface output values: p_rho_cp_cfxyz_yplus_fnormal
							          : -
 Monitoring ------------------------------------------------------: -
                                                                  : -
                                             Monitor history (0/1): 1
                                    Residual monitoring type (0/1): 0
                          Residual evaluation for monitoring (0/1): 1
                                                                  : -
                                                 Monitoring values: Rrho_Max-Rrho_X-max-Rrho_Y-max-Rrho_Z-max-Rrho_Rrhovx_Rrhovy_Rrhovz_RrhoE_Rrhonue_C-lift_C-drag_C-sidef_C-mx_C-my_C-mz_Angle-a_Max-y+
                                    Monitoring significant figures: 3_3_3_3_3_3_3_3_3_3_3_3_3_3_3_3_2_3
                                                                  : -
                             Extended coefficient monitoring (0/1): 0
                                                                  : -
                                         Global monitoring (0/1/2): 1
                        Enable logfile output on all domains (0/1): 0
                                                                  : -
###########################################################################################
                                                                  : -
 Relaxation ------------------------------------------------------: -
                                                                  : -
                                                 Relaxation solver: Backward_Euler
                                  Hold static velocity field (0/1): 0
                                                                  : -
 Backward Euler --------------------------------------------------: -
                                                                  : -
                                                     Linear solver: Lusgs
                                     Implicit overrelaxation omega: 1
                                      Implicit overrelaxation beta: 1
                  Turbulence model equations eigenvalue correction: 0.15
                                                                  : -
 LUSGS -----------------------------------------------------------: -
                                                                  : -
                                                Sgs stages maximum: 3
Lusgs treat rotating frame of reference source terms implicitly (0/1): 0
                                                                  : -
 Timestepsize ----------------------------------------------------: -
                                                                  : -
                                                        CFL number: 4
                                         CFL number (coarse grids): 2
                                         CFL number (large grad p): 2bashrc.sh
                                                                  : -
                                        Time step smoothing factor: 0
                                                                  : -
 Dual time -------------------------------------------------------: -
							          : -
				            Unsteady time stepping: dual
                    Unsteady activate inner iteration output (0/1): 1
                             Unsteady show pseudo time steps (0/1): 1
                                  Unsteady physical time step size: 0.005
                                     Unsteady physical time offset: 0
                             Unsteady computational time step size: -1
                                      Unsteady physical time steps: 1
                           Unsteady inner iterations per time step: 10
                  Minimum number of inner iterations per time step: 10
                                    Unsteady implicit scheme order: 2
                                      Unsteady extrapolation order: 0
                       Compute harmonics of global forces (0/1..n): 0
                                Compute harmonics on surface (0/1): 1
                                                                  : -
 Flow time averaging ---------------------------------------------: -
							          : -
				           Compute flow statistics: (none)
					     Reset flow statistics: mean_variance_meanvelgrad
					            	          : -
 Smoother --------------------------------------------------------: -
                                                                  : -
                                                 Residual smoother: Point_explicit
                                               Correction smoother: Point_explicit
                                                                  : -
 Preconditioning -------------------------------------------------: -
                                                                  : -
                                                   Preconditioning: (none)
                                                     Cut-off value: 1
                                                                  : -
###########################################################################################

 Gradients  ------------------------------------------------------: -
                                                                  : -
                                       Reconstruction of gradients: Green_Gauss
                                                                  : -
 Flux main -------------------------------------------------------: -
                                                                  : -
                                 Inviscid flux discretization type: Central
                                  Viscous flux type TSL/Full (0/1): 1
                                       Mixed inviscid fluxes (0/1): 0
                                                                  : -
 Central flux ----------------------------------------------------: -
                                                                  : -
                                        Central dissipation scheme: Matrix_dissipation
                                  Central convective meanflow flux: Average_of_flux
                             Scaling of central scheme dissipation: Standard
                                Central convective turbulence flux: Average_of_flux
                                 2nd order dissipation coefficient: 0.5
                         Inverse 4th order dissipation coefficient: 64
                            Version of cell stretching coefficient: TAU
                             Use modified dissipation for 2D (0/1): 0
                                                                  : -
 Matrix Dissipation ----------------------------------------------: -
                      	          Pressure switch weighting factor: 1
                 Minimum artificial dissipation for acoustic waves: 0.2
                       Minimum artificial dissipation for velocity: 0.2
                   Couple mean flow and turbulence equations (0/1): 0
                              	         Hartens entropy fix (0/1): 0
							          : -
###########################################################################################

 Turbulenc ------------------------------------------------------: -
                                                                  : -
                                                   Turbulence mode: RANS
                                          Turbulence model version: SA
                                                  SA model version: SA-neg
                                                                  : -
                                 #SA attractor for zero value (0/1): 1
                                         General ratio mue-t/mue-l: 100
                                            Ratio Prandtl lam/turb: 0.8
                                         Maximum limit mue-t/mue-l: 20000
                                       General turbulent intensity: 0.1
                                            Reference bl-thickness: 0
                                                 Positivity scheme: 1
                                             #EARSM expansion order: 1
                                   # Vortical flow correction (0/1): 0
                                    Turbulence initialization type: walldist
                                   #Vortex generator modeling (0/1): 0
                                                                  : -

###########################################################################################

  -----------------------------------------------------
 Parameters for RBF surface deflection generation in netcdf format
 -----------------------------------------------------
                             #Modal amplitude filename: genericpitch
                                     Deformation type: RBF
                               Modal boundary markers: 13
                    Modal maximum x rotation (degree): 0.0
                    Modal maximum y rotation (degree): 0.0
                    Modal maximum z rotation (degree): 0.0
                    Modal maximum x translation: 0
                    Modal maximum y translation: 0
                    Modal maximum z translation: 0
 -----------------------------------------------------
 DEFORMATION - put general deformation params here...
 -----------------------------------------------------
                                  # Amplitude filename: genericpitch.nc
                               Modal amplitude factor: 1.0
      Parameters for RBF deformation -----------------:
                                  RBF number of groups: 1
     Parameters for RBF group start ------------------:
                          RBF markers specifying group: 13
        RBF basis coordinates and deflections filename: interface_deformfile.nc
                               RBF surface format name: netcdf
                            RBF walldistances filename: walldistances_matrix
                                       RBF matrix name: RBF_matrix
                                              RBF name: volume-spline
                                RBF radius full weight: 10
                                RBF radius zero weight: 100
                RBF maximum number of base points used: 100
  group end

###########################################################################################

 -----------------------------------------------------
  Surface output
 -----------------------------------------------------
            #Surface output period: 1
  #Surface output description file: (thisfile)
           # Surface output values: cp_mach_p_rho_cfxyz_yplus_fnormal

 -----------------------------------------------------
     Output Control ----------------------------------: -
					#Output format: tecplot
				       Output format: ensight_gold
		      Tecplot ascii output precision : 6
					  Ascii (0/1): 1
					    Precision: 6
					Tecplot title: interfaceData.dat
Chimera component output (Tecplot, Ensight Gold) (0/1): 0
	     Tecio separate grid and solution (0/1/2): 0
		   	     Tecio format szplt (0/1): 0
  -----------------------------------------------------
 Extra field pointdata output
 -----------------------------------------------------
  Field output description file: (thisfile)
  Field output values: cp_mach_visc
    solver at Fri Apr 10 15:29:27 2020
                                              Restart-data prefix: /home/inigo/software/kratosMerge/Kratos/applications/CoSimulationApplication/tests/TAU/Outputs/airfoilSol.pval.unsteady_i=1_t=5.000e-03
                                            Reset flow statistics: (none)

    solver at Fri Apr 10 15:29:27 2020
                                          Surface output filename: /home/inigo/software/kratosMerge/Kratos/applications/CoSimulationApplication/tests/TAU/Outputs/airfoilSol.surface.pval.unsteady_i=1_t=5.000e-03
