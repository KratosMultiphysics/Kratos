from KratosMultiphysics import *
from KratosMultiphysics.ThermoMechanicalApplication import *
def DetectSolidusLiquidusTemperature(model_part):
    """ 
    Finds both the solidus and liquidus temperature that are included in the
    Solid Fraction table (Table 3) of the ModelPart. It scans from 100 to 5,000 
    in 0.001 intervals.
    @param model_part
    """
    lowlimit=100.0
    uplimit=5000.0
    dtemp=0.001
    solidustemp=lowlimit
    liquidustemp=uplimit
    temp=lowlimit
    table=model_part.GetTable(3)
    while(temp<liquidustemp):
        SF=table.GetValue(temp)
        if SF<=0.0:
            liquidustemp=temp
        if SF>=1.0:
            solidustemp=temp
        temp+=dtemp
    print("")
    print("Found Solidus Temperature: ",solidustemp)
    print("Found Liquidus Temperature: ",liquidustemp)    
    
    model_part.ProcessInfo.SetValue(FLUID_TEMPERATURE, liquidustemp)#ProjectParameters.FLUID_TEMPERATURE)
    model_part.ProcessInfo.SetValue(SOLID_TEMPERATURE, solidustemp)#ProjectParameters.SOLID_TEMPERATURE)
    return liquidustemp, solidustemp


