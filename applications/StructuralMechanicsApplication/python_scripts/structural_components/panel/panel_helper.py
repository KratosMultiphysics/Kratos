import KratosMultiphysics.StructuralMechanicsApplication as SMA



def NumberOfPlies(properties): #je Element

    number_plies = SMA.SHELL_ORTHOTROPIC_LAYERS # each row in the OTHOTROPIC_LAYERS matrix corresponds to a layer --> Number of rows = Number of layers

    if properties.Has(number_plies):
        return properties.GetValue(number_plies).Size1() 
        
    return 0


