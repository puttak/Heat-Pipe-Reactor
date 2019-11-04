import numpy as np
import openmc
import xml.etree.ElementTree as ET
import os

def define_Geo_Mat_Set(cells_num_dic,parameters_dic,settings_dic,temp_phy_mat,controlRod_deep):
    # Fuel U-10Mo
    fuel = openmc.Material(name='U-10Mo')
    fuel.set_density('g/cm3',16.8)
    fuel.add_element('Mo',0.1,'wo')
    # fuel.add_nuclide('U235',0.1773,'wo') # for LEU (19.7%)
    # fuel.add_nuclide('U238',0.7227,'wo')
    fuel.add_nuclide('U235',0.837,'wo') # for HEU (93.0%)
    fuel.add_nuclide('U238',0.063,'wo')

    # Structural Material HAYNES230
    structure_HAY = openmc.Material(name='HAYNES230')
    structure_HAY.set_density('g/cm3',8.97)
    structure_HAY.add_element('Ni',0.57,'wo')
    structure_HAY.add_element('Cr',0.22,'wo')
    structure_HAY.add_element('W',0.14,'wo')
    structure_HAY.add_element('Mo',0.02,'wo')
    structure_HAY.add_element('Fe',0.01875,'wo')
    structure_HAY.add_element('Co',0.03125,'wo')

    # Structural Material SS316
    structure_SS = openmc.Material(name='SS316')
    structure_SS.set_density('g/cm3',7.99)
    structure_SS.add_element('Ni',0.12,'wo')
    structure_SS.add_element('Cr',0.17,'wo')
    structure_SS.add_element('Mo',0.025,'wo')
    structure_SS.add_element('Mn',0.02,'wo')
    structure_SS.add_element('Fe',0.665,'wo')

    #Control Rod Material B4C
    ControlRod_B4C = openmc.Material(name='B4C')
    ControlRod_B4C.set_density('g/cm3',2.52)
    ControlRod_B4C.add_nuclide('B10',4,'ao')
    ControlRod_B4C.add_element('C',1,'ao')

    #Reflector Material BeO
    Reflector_BeO = openmc.Material(name='BeO')
    Reflector_BeO.set_density('g/cm3',3.025)
    Reflector_BeO.add_element('Be',1,'ao')
    Reflector_BeO.add_element('O',1,'ao')

    #Coolant Na
    Coolant_Na = openmc.Material(name='Na')
    Coolant_Na.set_density('g/cm3',0.76)
    Coolant_Na.add_element('Na',1,'ao')

    # Instantiate a Materials collection
    materials_file = openmc.Materials([fuel, structure_HAY, structure_SS, ControlRod_B4C, Reflector_BeO, Coolant_Na])

    # Export to "materials.xml"
    materials_file.export_to_xml()

    num_heat_pipe = 8

    # Parameters of reactor
    # Unit:cm

    fuel_r = parameters_dic['fuel_r']
    fuel_h = parameters_dic['fuel_h']

    controlRod_r = parameters_dic['controlRod_r']
    controlRod_h_max = parameters_dic['controlRod_h_max']
    controlRod_l = parameters_dic['controlRod_l']

    reflector_r = parameters_dic['reflector_r']
    reflector_h = parameters_dic['reflector_h']

    heat_pipe_R = parameters_dic['heat_pipe_R']
    heat_pipe_r = parameters_dic['heat_pipe_r']

    top_distance = parameters_dic['top_distance']
    bottom_distance = parameters_dic['bottom_distance']

    # Create cylinders for the fuel control rod and reflector
    fuel_OD = openmc.ZCylinder(x0=0.0, y0=0.0, r=fuel_r)
    controlRod_OD = openmc.ZCylinder(x0=0.0, y0=0.0, r=controlRod_r)
    reflector_OD = openmc.ZCylinder(x0=0.0, y0=0.0, r=reflector_r, boundary_type='vacuum')

    # Create planes for fuel control rod and reflector
    reflector_TOP = openmc.ZPlane(z0 = (top_distance+fuel_h/2),boundary_type='vacuum')
    reflector_BOTTOM = openmc.ZPlane(z0 = -(bottom_distance+fuel_h/2),boundary_type='vacuum')

    fuel_TOP = openmc.ZPlane(z0 = fuel_h/2)
    fuel_BOTTOM = openmc.ZPlane(z0 = -fuel_h/2)

    controlRodSpace_TOP = openmc.ZPlane(z0 = (controlRod_h_max-bottom_distance-fuel_h/2))

    # Create cylinders for heat pipes
    heat_pipe_OD = openmc.ZCylinder(x0=fuel_r, y0=0, r=heat_pipe_R)

    n_ang = num_heat_pipe
    ang_mesh = np.pi/n_ang

    # Number of cells
    n_r = cells_num_dic['n_r']
    n_r_outer = cells_num_dic['n_r_outer']
    n_h = cells_num_dic['n_h']

    r_mesh = np.linspace(controlRod_r,(fuel_r-heat_pipe_R),n_r+1)
    r_outer_mesh = np.linspace(fuel_r-heat_pipe_R,fuel_r,n_r_outer+1)
    h_mesh = np.linspace(-fuel_h/2,fuel_h/2,n_h+1)

    line_1 = openmc.Plane(a=np.tan(-ang_mesh),b=-1.0,c=0.0,d=0.0,boundary_type='reflective')
    line_2 = openmc.Plane(a=np.tan(ang_mesh),b=-1.0,c=0.0,d=0.0,boundary_type='reflective')

    # Get the vector of distance
    r_axe = np.zeros(n_r+n_r_outer)
    r_axe[0:n_r] = (r_mesh[0:n_r]+r_mesh[1:(n_r+1)])/2
    r_axe[n_r:(n_r+n_r_outer)]= (r_outer_mesh[0:n_r_outer]+r_outer_mesh[1:(n_r_outer+1)])/2

    h_axe = np.zeros(n_h)
    h_axe[0:n_h] = (h_mesh[0:n_h]+h_mesh[1:(n_h+1)])/2

    # Create volume vector and matrix
    volume_vec = np.zeros(n_r+n_r_outer)
    d_h = fuel_h/n_h

    for i in range(n_r+n_r_outer):
        if i >= n_r:
            d = heat_pipe_R*(i-n_r)/n_r_outer
            x_i = np.sqrt(2*heat_pipe_R*d-d*d)
            d = heat_pipe_R*(i-n_r+1)/n_r_outer
            x_i1 = np.sqrt(2*heat_pipe_R*d-d*d)
            s = (x_i+x_i1)*heat_pipe_R/n_r_outer
            volume_vec[i] = d_h*np.pi*(r_outer_mesh[i+1-n_r]*r_outer_mesh[i+1-n_r]-r_outer_mesh[i-n_r]*r_outer_mesh[i-n_r])/8-s
        else:
            volume_vec[i] = d_h*np.pi*(r_mesh[i+1]*r_mesh[i+1]-r_mesh[i]*r_mesh[i])/8

    volume_mat = np.zeros((n_h,(n_r+n_r_outer)))
    for i in range(n_h):
        volume_mat[i,:] = volume_vec

    # Create heat_pipe universe
    heat_pipe_Inner = openmc.ZCylinder(r=heat_pipe_r)

    coolant_cell = openmc.Cell(fill=Coolant_Na, region=(-heat_pipe_Inner & -reflector_TOP & +reflector_BOTTOM))
    pipe_cell = openmc.Cell(fill=structure_HAY, region=(+heat_pipe_Inner & -reflector_TOP & +reflector_BOTTOM))

    heat_pipe_universe = openmc.Universe(cells=(coolant_cell, pipe_cell))

    # Create a Universe to encapsulate a fuel pin
    pin_cell_universe = openmc.Universe(name='U-10Mo Pin')

    # Create fine-fuel-cell (num of cells: n_r + n_r_outer)
    fuel_cell_list = []
    fuel_cell_ID_list = []

    for j in range(n_h):
        cir_top = openmc.ZPlane(z0 = h_mesh[j+1])
        cir_bottom = openmc.ZPlane(z0 = h_mesh[j])
        for i in range(n_r):
            cir_in = openmc.ZCylinder(r=r_mesh[i])
            cir_out = openmc.ZCylinder(r=r_mesh[i+1])
            fuel_cell = openmc.Cell()
            fuel_cell.fill = fuel
            fuel_cell.region = +cir_in & -cir_out & +cir_bottom &-cir_top
            fuel_cell.temperature = temp_phy_mat[j,i]
            fuel_cell.id = (j+1)*10000+(i + 1)*100
            fuel_cell_ID_list.append((j+1)*10000+(i + 1)*100)
            fuel_cell_list.append(fuel_cell)
            pin_cell_universe.add_cell(fuel_cell)

        for i in range(n_r_outer):
            cir_in = openmc.ZCylinder(r=r_outer_mesh[i])
            cir_out = openmc.ZCylinder(r=r_outer_mesh[i+1])
            fuel_cell = openmc.Cell()
            fuel_cell.fill = fuel
            fuel_cell.region = +cir_in & -cir_out & +heat_pipe_OD & +cir_bottom &-cir_top
            fuel_cell.temperature = temp_phy_mat[j,i+n_r]
            fuel_cell.id = (j+1)*10000+(n_r + i + 1)*100
            fuel_cell_ID_list.append((j+1)*10000+(n_r + i + 1)*100)
            fuel_cell_list.append(fuel_cell)
            pin_cell_universe.add_cell(fuel_cell)

    # Create control rod Cell
    controlRod_TOP = openmc.ZPlane(z0 = (controlRod_deep-fuel_h/2-bottom_distance))
    controlRod_TOP.name = 'controlRod_TOP'
    if controlRod_deep < controlRod_h_max:
        controlRod_empty_cell = openmc.Cell(name='Control Rod Empty')
        controlRod_empty_cell.region = -controlRod_OD & +controlRod_TOP & -controlRodSpace_TOP
        pin_cell_universe.add_cell(controlRod_empty_cell)


    if controlRod_deep > 0:
        controlRod_cell = openmc.Cell(name='Control Rod')
        controlRod_cell.fill = ControlRod_B4C
        controlRod_cell.region = -controlRod_OD & +reflector_BOTTOM & -controlRod_TOP
        pin_cell_universe.add_cell(controlRod_cell)

    # Create heat pipe Cell
    heat_pipe_cell = openmc.Cell(name='Heat Pipe')
    heat_pipe_cell.fill = heat_pipe_universe
    heat_pipe_cell.region = -heat_pipe_OD & +reflector_BOTTOM & -reflector_TOP
    heat_pipe_cell.translation = (fuel_r,0,0)
    pin_cell_universe.add_cell(heat_pipe_cell)

    # Create reflector Cell
    reflector_radial_cell = openmc.Cell(name='Radial Reflector')
    reflector_radial_cell.fill = Reflector_BeO
    reflector_radial_cell.region = +fuel_OD & +heat_pipe_OD & +reflector_BOTTOM & -reflector_TOP
    pin_cell_universe.add_cell(reflector_radial_cell)

    reflector_bottom_cell = openmc.Cell(name='Bottom Reflector')
    reflector_bottom_cell.fill = Reflector_BeO
    reflector_bottom_cell.region = +controlRod_OD & +heat_pipe_OD & -fuel_OD  & -fuel_BOTTOM & +reflector_BOTTOM
    pin_cell_universe.add_cell(reflector_bottom_cell)

    reflector_top_cell = openmc.Cell(name='Top Reflector')
    reflector_top_cell.fill = Reflector_BeO
    reflector_top_cell.region = +heat_pipe_OD & -fuel_OD  & +controlRodSpace_TOP & -reflector_TOP
    reflector_top_cell.region = reflector_top_cell.region | (+controlRod_OD & +heat_pipe_OD & -fuel_OD & +fuel_TOP & -controlRodSpace_TOP)
    pin_cell_universe.add_cell(reflector_top_cell)

    # Create root Cell
    root_cell = openmc.Cell(name='root cell')
    root_cell.fill = pin_cell_universe

    # Add boundary planes
    root_cell.region = -reflector_OD & +line_2 & -line_1 & +reflector_BOTTOM & -reflector_TOP

    # Create root Universe
    root_universe = openmc.Universe(universe_id=0, name='root universe')
    root_universe.add_cell(root_cell)

    # Create Geometry and set root Universe
    geometry = openmc.Geometry(root_universe)

    # Export to "geometry.xml"
    geometry.export_to_xml()

    # Instantiate a Settings object
    settings_file = openmc.Settings()
    settings_file.batches = settings_dic['batches']
    settings_file.inactive = settings_dic['inactive']
    settings_file.particles = settings_dic['particles']
    settings_file.temperature['method']= 'interpolation'

    settings_file.source = openmc.Source(space=openmc.stats.Point((15,0,0)))
    # Export to "settings.xml"
    settings_file.export_to_xml()

    # Instantiate an empty Tallies object
    tallies_file = openmc.Tallies()

    for i in range(len(fuel_cell_ID_list)):
        tally = openmc.Tally(name='cell tally '+str(fuel_cell_ID_list[i]))
        tally.filters = [openmc.DistribcellFilter(fuel_cell_ID_list[i])]
        tally.scores = ['heating']
        tallies_file.append(tally)

    # Export to "tallies.xml"
    tallies_file.export_to_xml()
    
    return r_axe,h_axe,volume_mat,fuel_cell_ID_list

def postProcess(nodes_dic,volume_mat,temp_nodes_vec,parameters_dic,cells_num_dic,settings_dic,fuel_cell_ID_list):

    #To be edited


    # Parameters of heat pipes
    fuel_r = parameters_dic['fuel_r']
    controlRod_r = parameters_dic['controlRod_r']
    heat_pipe_R = parameters_dic['heat_pipe_R']
    fuel_h = parameters_dic['fuel_h']

    batches = settings_dic['batches']

    row = np.size(volume_mat,0)
    col = np.size(volume_mat,1)
    # Get tally data
    heat_tot_vec = np.zeros(row*col)

    sp = openmc.StatePoint('statepoint.{}.h5'.format(batches))
    k_vec = sp.k_generation[settings_dic['inactive']:]
    k_ave = k_vec.sum()/len(k_vec)

    for i in range(len(fuel_cell_ID_list)):
        t = sp.get_tally(name='cell tally '+str(fuel_cell_ID_list[i]))
        heat_tot_vec[i] = t.get_values(scores=['heating'],value='mean').item()

    heat_tot_mat = heat_tot_vec.reshape((row,col),order='C')
    # Define the power factor
    heat_power = parameters_dic['heat_power']
    heat_power = heat_power/8 #Power: 4kW 
    heat_power_origin = heat_tot_mat.sum()
    k_power = heat_power/(heat_power_origin*1.6022e-19) # Energy unit: eV---> J
    
    heat_ave_mat = heat_tot_mat*1.6022e-19*k_power/volume_mat

    #Get position of nodes
    x = nodes_dic['x']
    y = nodes_dic['y']
    z = nodes_dic['z']

    #Get cells mesh 
    n_r = cells_num_dic['n_r']
    n_r_outer = cells_num_dic['n_r_outer']
    n_h = cells_num_dic['n_h']

    r_inner_mesh = np.linspace(controlRod_r,(fuel_r-heat_pipe_R),n_r+1)
    r_outer_mesh = np.linspace(fuel_r-heat_pipe_R,fuel_r,n_r_outer+1)
    r_mesh = np.hstack((r_inner_mesh[0:len(r_inner_mesh)],r_outer_mesh[1:len(r_outer_mesh)]))

    h_mesh = np.linspace(-fuel_h/2,fuel_h/2,n_h+1)

    heat_nodes_vec = np.zeros(len(x))

    for i in range(row):
        for j in range(col):
            index = np.where(((x**2+y**2)>=(r_mesh[j]**2)) & ((x**2+y**2)<(r_mesh[j+1]**2)) & (z>=h_mesh[i]) & (z<h_mesh[i+1]) & ((x-fuel_r)*(x-fuel_r)+y*y>=(heat_pipe_R**2)))
            heat_nodes_vec[index] = heat_ave_mat[i,j]

    # Temp = 1173.5 # Temperature(approximation), unit: K
    # lamb = (0.606+0.0351*Temp)*0.01
    # Thermal conductivity. unit: W/(K.cm)

    lamb = (0.606+0.0351*temp_nodes_vec)*0.01
    points_force = -heat_nodes_vec/lamb

    editForceFile(x,y,z,points_force,'Force')

    return k_ave


def editCellTemperature(fuel_temp,fuel_cell_ID_list):
    #fuel_temp(matrix)
    #fuel_cell_ID_list(list)
    row = np.size(fuel_temp,0)
    col = np.size(fuel_temp,1)
    Y = fuel_temp.reshape(((row*col)),order='C')

    tree = ET.parse('geometry.xml')
    root = tree.getroot()

    k = 0
    for cell in root.iter('cell'):
        if cell.attrib['id']==str(fuel_cell_ID_list[k]):
            cell.attrib['temperature'] = str(Y[k])
            k = k+1
        else:
            continue
    tree.write('geometry.xml')

def editForceFile(x,y,z,points_force,file_name):
    mat = np.zeros((len(x),4))
    mat[:,0] = x
    mat[:,1] = y
    mat[:,2] = z
    mat[:,3] = points_force

    str_tot = " "
    for i in range(np.size(mat,0)):
        list_float = mat[i,:].tolist()
        list_str = map(str,list_float)
        str_1 = " ".join(list_str)
        str_tot = str_tot + '\n' +str_1
    str_tot = str_tot + '\n'+ " "

    if os.path.exists(file_name +'.pts'):
        tree = ET.parse(file_name+'.pts')
        root = tree.getroot()
        for points in root.findall('POINTS'):
            points.text = str_tot
        tree.write(file_name+'.pts')
    else:
        NEKTAR = ET.Element('NEKTAR')
        POINTS = ET.SubElement(NEKTAR,'POINTS',{'DIM':'3','FIELDS':'u'})
        POINTS.text = str_tot
        tree = ET.ElementTree(NEKTAR) 
        tree.write(file_name+'.pts',encoding="utf-8", xml_declaration=True)

def getCellTemperature(nodes_dic,temp_nodes_vec,parameters_dic,cells_num_dic):
    #Get position of nodes
    x = nodes_dic['x']
    y = nodes_dic['y']
    z = nodes_dic['z']

    #Get cells mesh 
    n_r = cells_num_dic['n_r']
    n_r_outer = cells_num_dic['n_r_outer']
    n_h = cells_num_dic['n_h']

    controlRod_r = parameters_dic['controlRod_r']
    fuel_r = parameters_dic['fuel_r']
    heat_pipe_R = parameters_dic['heat_pipe_R']
    fuel_h = parameters_dic['fuel_h']

    r_inner_mesh = np.linspace(controlRod_r,(fuel_r-heat_pipe_R),n_r+1)
    r_outer_mesh = np.linspace(fuel_r-heat_pipe_R,fuel_r,n_r_outer+1)
    r_mesh = np.hstack((r_inner_mesh[0:len(r_inner_mesh)],r_outer_mesh[1:len(r_outer_mesh)]))

    h_mesh = np.linspace(-fuel_h/2,fuel_h/2,n_h+1)

    temp_cells_mat = np.zeros((n_h,(n_r+n_r_outer)))

    for i in range(n_h):
        for j in range((n_r+n_r_outer)):
            index = np.where(((x**2+y**2)>=(r_mesh[j]**2)) & ((x**2+y**2)<(r_mesh[j+1]**2)) & (z>=h_mesh[i]) & (z<h_mesh[i+1]) & ((x-fuel_r)*(x-fuel_r)+y*y>=(heat_pipe_R**2)))
            if len(temp_nodes_vec[index])==0:
                if j>0:
                    temp_cells_mat[i,j] =  temp_cells_mat[i,j-1]
                else:
                    temp_cells_mat[i,j] =  temp_cells_mat[i-1,j]

            else:
                temp_cells_mat[i,j] = temp_nodes_vec[index].sum()/len(temp_nodes_vec[index])

    return temp_cells_mat


