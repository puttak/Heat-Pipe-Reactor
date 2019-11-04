import xml.etree.ElementTree as ET
from subprocess import call
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import os

def runSolver(solver_name,file_name,iteration):
    call('$NEK/'+solver_name+' '+file_name+'.xml', shell=True)
    if iteration > 0:
        file_name_new = file_name + '.bak'+str(iteration-1)
    else:
        file_name_new = file_name
    call('$NEK/FieldConvert '+file_name+'.xml '+file_name_new+'.fld '+file_name_new+'.vtu',shell=True)
    return file_name_new

def editPipeBoundary(root,coolant_temp):
    for regions in root.iter('REGION'):
        if regions.attrib['REF'] == '5':
            for region in regions:
                region.attrib['VALUE']=str(coolant_temp)
        else:
            continue
    return root

def runNektar(file_name,solver_name,temp_pipe,iteration):
    # Read mesh xml
    tree = ET.parse(file_name+'.xml')
    root = tree.getroot()

    # Edit forcing function and boundary condition in xml
    root = editPipeBoundary(root,temp_pipe)
    tree.write(file_name+'.xml')

    #To be edited
    # Choose solver and run nektar++

    dir = os.getcwd()
    call('cd '+dir,shell=True)
    file_name_new = runSolver(solver_name,file_name,iteration)

    return file_name_new

def postProcess(file_name):
    # Read the result
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_name+'.vtu')
    reader.Update()
    output = reader.GetOutput()

    # Get temperature
    temperature_vtk_array = reader.GetOutput().GetPointData().GetArray(0)
    temp_nodes_vec = vtk_to_numpy(temperature_vtk_array)

    return temp_nodes_vec

def readNodesFromVtu(file_name):
    
    file_name = file_name+'_node'

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_name+'.vtu')
    reader.Update()
    output = reader.GetOutput()

    nodes_vtk_array = reader.GetOutput().GetPoints().GetData()

    # Get position 
    nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
    x,y,z= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1] , nodes_nummpy_array[:,2]
    
    
    return x,y,z






