import xml.etree.ElementTree as ET
from subprocess import call
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
from scipy.optimize import curve_fit
import os

def editForceFunction(root,force_func):
    for functions in root.iter('FUNCTION'):
        if functions.attrib['NAME'] == 'Forcing':
            for function in functions:
                function.attrib['VALUE']=force_func
        else:
            continue
    return root

def runSolver(solver_name,file_name,iteration):
    call('$NEK/'+solver_name+' '+file_name+'.xml', shell=True)
    if iteration > 0:
        file_name_new = file_name + '.bak'+str(iteration-1)
    else:
        file_name_new = file_name
    call('$NEK/FieldConvert '+file_name+'.xml '+file_name_new+'.fld '+file_name_new+'.vtu',shell=True)
    return file_name_new

def readElements_num(root):
    counts = 0
    for ELEMENT in root.iter('ELEMENT'):
        for Q in ELEMENT:
            counts = counts +1
    return counts  

def editPipeBoundary(root,coolant_temp):
    for regions in root.iter('REGION'):
        if regions.attrib['REF'] == '3':
            for region in regions:
                region.attrib['VALUE']=str(coolant_temp)
        else:
            continue
    return root

def postProcess(file_name,controlRod_r=2.2,fuel_r=5.5,heat_pipe_R=0.635):
    # Read the result
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_name+'.vtu')
    reader.Update()
    output = reader.GetOutput()

    nodes_vtk_array = reader.GetOutput().GetPoints().GetData()
    temperature_vtk_array = reader.GetOutput().GetPointData().GetArray(0)
    
    # Get position 
    nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
    x,y= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1]

    # Get temperature
    temperature_numpy_array = vtk_to_numpy(temperature_vtk_array)
    T = temperature_numpy_array
    
    # Find data in fuel area
    index = np.where((x*x+y*y>=(controlRod_r**2)) & (x*x+y*y<=(fuel_r**2)) & ((x-fuel_r)*(x-fuel_r)+y*y>=(heat_pipe_R**2)))
    x_fuel = x[index]
    y_fuel = y[index]
    T_fuel = T[index]
    
    r_vec = np.sqrt(x_fuel*x_fuel + y_fuel*y_fuel)
    
    # Fit the radial distribution of temperature
    T_fuel_nor = T_fuel/1
    r_vec_nor = r_vec/1

    popt,pcov = curve_fit(func,r_vec_nor,T_fuel_nor)
    
    temp_func = '('+str(popt[0])+'*(sqrt(x*x+y*y))^3+'+str(popt[1])+')'
    
    return T,x,y,popt,temp_func

def runNektar(file_name,force_func,solver_name,temp_pipe,iteration):
    # input_file_name (string)
    # file name (string)
    # force_func (string)
    # solver_name (string)
    # temp_pipe (float)

    # Read mesh xml
    tree = ET.parse(file_name+'.xml')
    root = tree.getroot()

    # Edit forcing function and boundary condition in xml
    root = editForceFunction(root,force_func)
    root = editPipeBoundary(root,temp_pipe)
    tree.write(file_name+'.xml')

    #To be edited
    # Choose solver and run nektar++

    dir = os.getcwd()
    call('cd '+dir,shell=True)
    file_name_new = runSolver(solver_name,file_name,iteration)

    return file_name_new

def func(r,a,b):
    return a*r**3+b

def cellTemperature(r_cells,popt):
    n_cells = len(r_cells)
    fuel_temp = np.array([func(i,popt[0],popt[1]) for i in r_cells])
    return fuel_temp
