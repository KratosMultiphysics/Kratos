import KratosMultiphysics.StructuralMechanicsApplication as SMA



def NumberOfPlies(properties): #je Element

    number_plies = SMA.SHELL_ORTHOTROPIC_LAYERS # Je Zeile in der Matrix aus OTHOTROPIC_LAYERS entspricht einer Lage --> Anzahl Zeilen = Anzahl Lagen

    if properties.Has(number_plies):
        return properties.GetValue(number_plies).Size1() #Checken ob len(hier den richtigen Output liefert) 
        
    return 0


