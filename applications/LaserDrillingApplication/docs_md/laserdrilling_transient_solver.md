Module laserdrilling_transient_solver
=====================================

Functions
---------

`CreateSolver(model, custom_settings)`
:   

Classes
-------

`LaserDrillingTransientSolver(model, custom_settings)`
:   The transient class for convection-diffusion solvers.
    
    Public member variables:
    transient_settings -- settings for the implicit dynamic solvers.
    
    See convection_diffusion_solver.py for more information.
    
    The constructor of the PythonSolver-Object.
    
    It is intended to be called from the constructor
    of deriving classes:
    super().__init__(settings)
    
    Keyword arguments:
    self -- It signifies an instance of a class.
    model -- The Model to be used
    settings -- The solver settings used

    ### Ancestors (in MRO)

    * KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver
    * KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_solver.ConvectionDiffusionSolver
    * KratosMultiphysics.python_solver.PythonSolver
    * abc.ABC

    ### Static methods

    `GetDefaultParameters()`
    :   This function returns the default-settings used by this class

    ### Methods

    `AddDecomposedNodesToSurfaceList(self)`
    :

    `AdjustTemperatureFieldAfterAblation(self)`
    :

    `AllocateKratosMemory(self)`
    :

    `ComputeIonizationEnergyPerUnitVolumeThreshold(self)`
    :   Must be implemented by child classes.

    `ComputeMaximumAblationRadius(self)`
    :

    `ComputeMaximumDepth(self)`
    :

    `ComputeOpticalPenetrationDepth(self)`
    :

    `ComputePeakFluence(self)`
    :   Computes the peak fluence of a gaussian pulse from its energy and waist radius
        Source: Woodfield 2024, eq (5)
        
        Parameters
        ----------
        None
        
        Returns
        -------
        The peak fluence

    `ComputePulseHoleAndAddToTotalHole(self)`
    :

    `ComputeSpotDiameter(self)`
    :

    `CountActiveElements(self)`
    :

    `CreateResultsFile(self, filename)`
    :

    `EnergyPerUnitArea1D(self, radius)`
    :

    `EvaporationDepth(self, r)`
    :   Calculates the depth of the ablated cavity as a function of radial position for a single gaussian pulse according to eq. (6) in Woodfield (2024).
        
        Parameters
        ----------
        r: float
            Radial coordinate with respect to the axis of the beam
        
        Returns
        -------
        float
            The depth of the cavity at the specified radial coordinate r

    `Finalize(self)`
    :   This function finalizes the PythonSolver
        Usage: It is designed to be called ONCE, AFTER the execution of the solution-loop

    `FinalizeSolutionStep(self)`
    :   This function performs all the required operations that should be executed
        (for each step) AFTER solving the solution step.

    `IdentifyInitialSurfaceNodes(self)`
    :

    `ImposeLaserDeltaTemperature(self)`
    :   Calls a function that applies the temperature increase due to the laser pulse.
        Depending on the parameters of the simulation, it decides which function to call.
        For example, it chooses between applying the pulse with or without refraction.
        Parameters
        ----------
        None
        
        Returns
        -------
        None

    `ImposeTemperatureIncreaseDueTo1DConduction(self)`
    :

    `ImposeTemperatureIncreaseDueToLaser(self)`
    :   Increases the temperature as an effect of the energy deposition by the laser pulse.
        Does not take into account the refraction of the ray at the boundary air-solid.
        Must be implemented by child classes.

    `ImposeTemperatureIncreaseDueToLaserWithRefraction(self)`
    :   Increases the temperature as an effect of the energy deposition by the laser pulse.
        Takes into account the refraction of the ray at the boundary air-solid.
        Must be implemented by child classes.

    `InitialThermalConductionTime(self, radius)`
    :

    `Initialize(self)`
    :   Perform initialization after adding nodal variables and dofs to the main model part.

    `InitializeSolutionStep(self)`
    :   This function performs all the required operations that should be executed
        (for each step) BEFORE solving the solution step.

    `MonitorDecomposedVolume(self)`
    :   Tallies up the volume of the decomposed elements
        
        Parameters
        ----------
        None
        
        Returns
        -------
        The sum of DECOMPOSED_ELEMENTAL_VOLUME

    `MonitorEnergy(self)`
    :   Tallies up the energy of the (not decomposed) elements
        
        Parameters
        ----------
        None
        
        Returns
        -------
        The sum of THERMAL_ENERGY over all ACTIVE elements

    `PenetrationDepthEstimation(self)`
    :

    `PrintDecomposedVolumeEvolution(self)`
    :

    `RemoveElementsByAblation(self)`
    :

    `RemoveElementsByEvaporation(self)`
    :

    `ResetTemperatureField(self)`
    :

    `ResidualHeatStage(self)`
    :   TODO: Currently, unused. It is overriden by LaserDrillingTransienSolverAblationPlusThermal.ResidualHeatStage

    `RetrievePreEvaporationTemperatureState(self)`
    :

    `SetParameters(self)`
    :

    `SetUpGNUPlotFiles(self)`
    :

    `SetUpHDF5Files(self)`
    :

    `SetUpResultsFiles(self)`
    :

    `SolveSolutionStep(self)`
    :   This function solves the current step.
        It can be called multiple times within one solution step
        Returns whether the problem is converged

    `StorePreEvaporationTemperature(self)`
    :

    `TemperatureVariationInZDueToLaser1D(self, radius, z)`
    :

    `UpdateLaserRelatedParameters(self)`
    :

    `WriteResults(self, filename, process_info)`
    :

    `sortSecond(self, val)`
    :   Helper function used, for instance, as a key for sorting a list. Returns the second element of the list.
        
        Parameters
        ----------
        val: list
             The list.
        
        Returns
        -------
        The second element of the list