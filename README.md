# Heat-Pipe-Reactor
# Instructions (version 0.0.1)
# Before use: Make sure that there is no .vtu or .fld existing in folder !!!!!
# This might cause some unexpected bugs!!!!!

This is simulation code for heat pipe reactor (KRUSTY 1/8 model).

Completed:

    (a)Neutronics calculation 
    (b)Temperature calculation in solid area

To do:

    (a)Fluid calculation in heat pipe
    (b)Thermal expansion calculation
    (c)Burn-up calculation
    (d)3D-modeling


1.Unit: length(cm),temperature(K), others are SI units

2.Define_OpenMC.py includes functions for operating OpenMC:

    (a)define_Geo_Mat_Set: This function is used to generate  geometry.xml, material.xml and settings.xml
    
    (b)postProcess: This function is used to post-proccess heating source distribution and generate forcing function for Nektar++
    
    (c)editCellTemperature: This function is used to edit temperature of cells in geometry.xml

3.Define_Nektar.py includes functions for operating Nektar++:

    (a)runNektar: Just as its name. Of cause, It can generate .fld and convert it to .vtu that we need
    
    (b)postProcess: It can post-process results of temperature and generate temperature function in fuel area and coefficients of this function
    
    (c)fuelTemperature: It can use temperature function to generate avarage temperature in cells

For more details, please check codes respectively
