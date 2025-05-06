Module laserdrilling_transient_solver_ablation_plus_thermal
===========================================================

Functions
---------

`CreateSolver(model, custom_settings)`
:   

Classes
-------

`LaserDrillingTransientSolverAblationPlusThermal(model, custom_settings)`
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

    * KratosMultiphysics.LaserDrillingApplication.laserdrilling_transient_solver.LaserDrillingTransientSolver
    * KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver
    * KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_solver.ConvectionDiffusionSolver
    * KratosMultiphysics.python_solver.PythonSolver
    * abc.ABC

    ### Methods

    `ComputeIonizationEnergyPerUnitVolumeThreshold(self)`
    :   Must be implemented by child classes.

    `ComputePulseVolume(self)`
    :

    `Finalize(self)`
    :   Overwrites LaserDrillingTransientSolver.Finalize

    `ImposeTemperatureIncreaseDueToLaser(self)`
    :   Increases the temperature as an effect of the energy deposition by the laser pulse.
        Does not take into account the refraction of the ray at the boundary air-solid.
        Must be implemented by child classes.

    `ImposeTemperatureIncreaseDueToLaserWithRefraction(self)`
    :   Increases the temperature as an effect of the energy deposition by the laser pulse.
        Takes into account the refraction of the ray at the boundary air-solid.
        Must be implemented by child classes.

    `RemoveElementsByAblation(self)`
    :   Removes elements by ablation. (this comment is a WIP)
        
        Overrides LaserDrillingTransientSolver.RemoveElementsByAblation
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None

    `RemoveElementsUsingEnergyPerVolumeThreshold(self)`
    :

    `ResidualHeatStage(self)`
    :   Overrides LaserDrillingTransientSolver.ResidualHeatStage

    `SolveSolutionStep(self)`
    :   TODO: Overrides LaserDrillingTransientSolver.SolveSolutionStep

    `TemperatureVariationDueToLaser(self, radius, z)`
    :   Computes the temperature increase caused by the laser in a specified position.
        
        Parameters
        ----------
        radius: float
            The radial coordinate of the point
        z: float
            The axial coordinate of the point
        
        Returns
        -------
        The temperature increase in kelvins