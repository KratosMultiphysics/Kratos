from KratosMultiphysics import *
from KratosMultiphysics.ThermoMechanicalApplication import *
def DetectSolidusLiquidusTemperature(model_part):
	lowlimit=0.0
	uplimit=10000.0
	dtemp=0.01
	solidustemp=lowlimit
	liquidustemp=uplimit
	temp=lowlimit
	while(temp<liquidustemp):
		SF=model_part.GetTable(3).GetValue(temp)
		if SF<=0.0:
			liquidustemp=temp
		if SF>=1.0:
			solidustemp=temp
		temp+=dtemp
		print(temp)
		print(SF)
	print("SOLIDUS TEMP FOUND :=" +str(solidustemp))
	print("LIQUIDUS TEMP FOUND :=" +str(liquidustemp))
	model_part.ProcessInfo.SetValue(FLUID_TEMPERATURE, liquidustemp)#ProjectParameters.FLUID_TEMPERATURE)
	model_part.ProcessInfo.SetValue(SOLID_TEMPERATURE, solidustemp)#ProjectParameters.SOLID_TEMPERATURE)


