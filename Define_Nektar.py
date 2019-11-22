import xml.etree.ElementTree as ET
from subprocess import call
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import os

# For simualting temperature field

def runSolver_Temp(solver_name,file_name,iteration):
    call('$NEK/'+solver_name+' '+file_name+'.xml', shell=True)
    if iteration > 0:
        file_name_new = file_name + '.bak'+str(iteration-1)
    else:
        file_name_new = file_name
    call('$NEK/FieldConvert '+file_name+'.xml '+file_name_new+'.fld '+file_name_new+'.vtu',shell=True)
    return file_name_new

def editPipeBoundary(root,temp_pipe):
    for regions in root.iter('REGION'):
        if regions.attrib['REF'] == '5':
            for region in regions:
                region.attrib['VALUE']=str(temp_pipe)
        else:
            continue
    return root

def runNektar_Temp(file_name,solver_name,temp_pipe,iteration):
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
    file_name_new = runSolver_Temp(solver_name,file_name,iteration)

    return file_name_new

def postProcess_Temp(file_name):
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

# For simulating structure displacement

def getTempDiff(temp_xml_name,fld_name_1,fld_name_2):

    call('$NEK/FieldConvert -m addfld:fromfld='+fld_name_1+'.fld:scale=-1 '+temp_xml_name+'.xml '+fld_name_2+'.fld \\temp_diff.fld', shell=True)
#    call('$NEK/FieldConvert -m gradient '+temp_xml_name+'.xml '+'temp_diff.fld '+'temp_grad.fld')

# def editForceFile_Structure(structure_xml_name,thermal_property_dic):

#     E = thermal_property_dic['E']
#     alpha = thermal_property_dic['alpha']
#     mu = thermal_property_dic['mu']
#     tree = ET.parse(file_name+'-conditions.xml')
#     root = tree.getroot()
#     for functions in root.findall('FUNCTION'):
#         if functions.attrib['NAME']=='Forcing':
#             for function in functions:
#                 function.attrib['VAR']='u,v,w'
#                 function.attrib['VALUE']=str(-E*alpha/(1-mu))
#                 function.attrib['FILE']='temp_grad.fld'
#     tree.write(file_name+'-conditions.xml',encoding="utf-8", xml_declaration=True)

# def editBoundary_Structure(structure_xml_name,thermal_property_dic,temp_surface):

#     E = thermal_property_dic['E']
#     alpha = thermal_property_dic['alpha']
#     mu = thermal_property_dic['mu']
#     tree = ET.parse(file_name+'-conditions.xml')
#     root = tree.getroot()
#     for regions in root.findall('REGION'):
#         if regions.attrib['REF']=='0':
#             for N in regions:
#                 if N.attrib['VAR']='u':
#                     N.attrib['VALUE']=str(E*alpha*temp_surface/(1-mu))+'*x/sqrt(x*x+y*y)'
#                 elif N.attrib['VAR']='v':
#                     N.attrib['VALUE']=str(E*alpha*temp_surface/(1-mu))+'*y/sqrt(x*x+y*y)'
#                 elif N.attrib['VAR']='w':
#                     N.attrib['VAR']=str(0)
#     tree.write(file_name+'-conditions.xml',encoding="utf-8", xml_declaration=True)
# def editForceFile_Structure(structure_xml_name,thermal_property_dic):
#     tree = ET.parse(file_name+'-conditions.xml')
#     root = tree.getroot()
#     for functions in root.findall('FUNCTION'):
#         if functions.attrib['NAME']=='Temperature':
#             for function in functions:
#                 function.atrrib['FILE']='temp_diff.fld'
#     tree.write(file_name+'-conditions.xml',encoding="utf-8", xml_declaration=True)

# def runNektar_Temp(structure_xml_name,solver_name,thermal_property_dic,iteration):
#     # Read mesh xml
#     tree = ET.parse(structure_xml_name+'.xml')
#     root = tree.getroot()

#     # Edit forcing function and boundary condition in xml
#     root = editThermalProperty(root,thermal_property_dic)
#     tree.write(file_name+'.xml')

#     #To be edited
#     # Choose solver and run nektar++

#     dir = os.getcwd()
#     call('cd '+dir,shell=True)
#     file_name_new = runSolver_Temp(solver_name,file_name,iteration)

#     return file_name_new

# def editThermalProperty(root,thermal_property_dic):
#     E = thermal_property_dic['E']
#     nu = thermal_property_dic['nu']
#     beta = thermal_property_dic['beta']



