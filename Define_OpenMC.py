import numpy as np
import pandas as pd
import openmc
from scipy.optimize import curve_fit
import xml.etree.ElementTree as ET

def define_Geo_Mat_Set(n_r,n_r_outer,temp_phy_vec,batches=30,inactive=10,particles=1000):
    # n_r(int)
    # n_r_outer(int)
    # temp_phy_vec(narray)
    # batches(int)
    # inactive(int)
    # particles(int)
    
    # Fuel U-10Mo
    fuel = openmc.Material(name='U-10Mo')
    fuel.set_density('g/cm3',16.8)
    fuel.add_element('Mo',0.1,'wo')
    # fuel.add_nuclide('U235',0.1773,'wo') # for LEU (19.7%)
    # fuel.add_nuclide('U238',0.7227,'wo')
    fuel.add_nuclide('U235',0.72,'wo') # for HEU (80.0%)
    fuel.add_nuclide('U238',0.18,'wo')

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
    
    # Parameters of reactor
    # Unit:cm

    fuel_r = 11/2
    controlRod_r = 4.4/2
    reflector_r = 33/2

    # Parameters of heat pipes
    num_heat_pipe = 8
    heat_pipe_R = 1.27/2
    heat_pipe_r = 1.27/2-0.1

    # Create cylinders for the fuel, control rod and reflector
    fuel_OD = openmc.ZCylinder(x0=0.0, y0=0.0, r=fuel_r)
    controlRod_OD = openmc.ZCylinder(x0=0.0, y0=0.0, r=controlRod_r)
    reflector_OD = openmc.ZCylinder(x0=0.0, y0=0.0, r=reflector_r, boundary_type='vacuum')

    # Create cylinders for heat pipes
    heat_pipe_OD = openmc.ZCylinder(x0=fuel_r, y0=0, r=heat_pipe_R)

    # Create mesh and boudary
    n_ang = num_heat_pipe

    #n_r = 20
    #n_r_outer = 10

    ang_mesh = np.pi/n_ang
    r_mesh = np.linspace(controlRod_r,(fuel_r-heat_pipe_R),n_r+1)
    r_outer_mesh = np.linspace(fuel_r-heat_pipe_R,fuel_r,n_r_outer+1)

    line_1 = openmc.Plane(a=np.tan(-ang_mesh),b=-1.0,c=0.0,d=0.0,boundary_type='reflective')
    line_2 = openmc.Plane(a=np.tan(ang_mesh),b=-1.0,c=0.0,d=0.0,boundary_type='reflective')

    # Get the vector of distance
    r_axe = np.zeros(n_r+n_r_outer)
    r_axe[0:n_r] = (r_mesh[0:n_r]+r_mesh[1:(n_r+1)])/2
    r_axe[n_r:(n_r+n_r_outer)]= (r_outer_mesh[0:n_r_outer]+r_outer_mesh[1:(n_r_outer+1)])/2

    
    # Create volume vector
    volume_vec = np.zeros(n_r+n_r_outer)
    for i in range(n_r+n_r_outer):
        if i >= n_r:
            d = heat_pipe_R*(i-n_r)/n_r_outer
            x_i = np.sqrt(2*heat_pipe_R*d-d*d)
            d = heat_pipe_R*(i-n_r+1)/n_r_outer
            x_i1 = np.sqrt(2*heat_pipe_R*d-d*d)
            s = (x_i+x_i1)*heat_pipe_R/n_r_outer
            volume_vec[i] = np.pi*(r_outer_mesh[i+1-n_r]*r_outer_mesh[i+1-n_r]-r_outer_mesh[i-n_r]*r_outer_mesh[i-n_r])/8-s
        else:
            volume_vec[i] = np.pi*(r_mesh[i+1]*r_mesh[i+1]-r_mesh[i]*r_mesh[i])/8

    # # Create fuel universe
    # fuel_cell0 = openmc.Cell()
    # fuel_cell0.fill = fuel
    # fuel_surf = openmc.ZCylinder(r=0)
    # fuel_cell0.region = +fuel_surf

    # fuel_universe = openmc.Universe()
    # fuel_universe.add_cell(fuel_cell0)

    # Create heat_pipe universe
    heat_pipe_Inner = openmc.ZCylinder(r=heat_pipe_r)
    coolant_cell = openmc.Cell(fill=Coolant_Na, region=-heat_pipe_Inner)
    pipe_cell = openmc.Cell(fill=structure_HAY, region=+heat_pipe_Inner)

    heat_pipe_universe = openmc.Universe(cells=(coolant_cell, pipe_cell))

    # Create a Universe to encapsulate a fuel pin
    pin_cell_universe = openmc.Universe(name='U-10Mo Pin')

    # Create fine-fuel-cell (num of cells: n_r + n_r_outer)
    fuel_cell_list = []
    fuel_cell_ID_list = []

    for i in range(n_r):
        cir_in = openmc.ZCylinder(r=r_mesh[i])
        cir_out = openmc.ZCylinder(r=r_mesh[i+1])
        fuel_cell = openmc.Cell()
        fuel_cell.fill = fuel
        fuel_cell.region = +cir_in & -cir_out 
        fuel_cell.temperature = temp_phy_vec[i]
        fuel_cell.id = (i + 1)*100
        fuel_cell_ID_list.append((i + 1)*100)
        fuel_cell_list.append(fuel_cell)
        pin_cell_universe.add_cell(fuel_cell)

    for i in range(n_r_outer):
        cir_in = openmc.ZCylinder(r=r_outer_mesh[i])
        cir_out = openmc.ZCylinder(r=r_outer_mesh[i+1])
        fuel_cell = openmc.Cell()
        fuel_cell.fill = fuel
        fuel_cell.region = +cir_in & -cir_out & +heat_pipe_OD
        fuel_cell.temperature = temp_phy_vec[i+n_r]
        fuel_cell.id = (n_r + i + 1)*100
        fuel_cell_ID_list.append((n_r + i + 1)*100)
        fuel_cell_list.append(fuel_cell)
        pin_cell_universe.add_cell(fuel_cell)

    # Create control rod Cell
    controlRod_cell = openmc.Cell(name='Control Rod')
    # controlRod_cell.fill = ControlRod_B4C
    controlRod_cell.region = -controlRod_OD 
    pin_cell_universe.add_cell(controlRod_cell)

    # Create heat pipe Cell
    heat_pipe_cell = openmc.Cell(name='Heat Pipe')
    heat_pipe_cell.fill = heat_pipe_universe
    heat_pipe_cell.region = -heat_pipe_OD
    heat_pipe_cell.translation = (fuel_r,0,0)
    pin_cell_universe.add_cell(heat_pipe_cell)

    # Create reflector Cell
    reflector_cell = openmc.Cell(name='Reflector')
    reflector_cell.fill = Reflector_BeO
    reflector_cell.region = +fuel_OD & +heat_pipe_OD 
    pin_cell_universe.add_cell(reflector_cell)


    # Create root Cell
    root_cell = openmc.Cell(name='root cell')
    root_cell.fill = pin_cell_universe

    # Add boundary planes
    root_cell.region = -reflector_OD & +line_2 & -line_1

    # Create root Universe
    root_universe = openmc.Universe(universe_id=0, name='root universe')
    root_universe.add_cell(root_cell)
    
    # Create Geometry and set root Universe
    geometry = openmc.Geometry(root_universe)
    
    # Export to "geometry.xml"
    geometry.export_to_xml()
    
    # OpenMC simulation parameters
#     batches = 30
#     inactive = 10
#     particles = 10000

    # Instantiate a Settings object
    settings_file = openmc.Settings()
    settings_file.batches = batches
    settings_file.inactive = inactive
    settings_file.particles = particles
    settings_file.temperature['method']= 'interpolation'

    settings_file.source = openmc.Source(space=openmc.stats.Point((15,0,0)))
    # Export to "settings.xml"
    settings_file.export_to_xml()
    
    # Instantiate an empty Tallies object
    tallies_file = openmc.Tallies()
    
    # Create cell tally to score flux and fission rate
    for i in range(n_r+n_r_outer):
        tally = openmc.Tally(name='cell tally '+str(fuel_cell_ID_list[i]))
        tally.filters = [openmc.DistribcellFilter(fuel_cell_ID_list[i])]
        tally.scores = ['heating']
        tallies_file.append(tally)

    # Create energy tally to score flux
    # energy_bins = np.logspace(np.log10(1e-2), np.log10(20.0e6), 20)
    # fine_energy_filter = openmc.EnergyFilter(energy_bins)
    # tally = openmc.Tally(name='energy tally')
    # tally.filters.append(fine_energy_filter)
    # tally.scores = ['flux']
    # tallies_file.append(tally)
    
    # Export to "tallies.xml"
    tallies_file.export_to_xml()
    
    return r_axe,volume_vec,fuel_cell_ID_list

def postProcess(heat_power,temp_func,r_axe,volume_vec,fuel_cell_ID_list,batches=30):
    # heat_power(float)
    # Temp_func(string)
    # r_axe(narray)
    # volume_vec(narray)
    # fuel_cell_ID_list(list)
    # batches(int)

    #To be edited
    fuel_r = 11/2
    controlRod_r = 4.4/2


    # Parameters of heat pipes
    heat_pipe_R = 1.27/2


    n_cells = len(r_axe)
    # Get tally data
    heat_tot = np.zeros(n_cells)
    sp = openmc.StatePoint('statepoint.{}.h5'.format(batches))
    for i in range(n_cells):
        t = sp.get_tally(name='cell tally '+str(fuel_cell_ID_list[i]))
        df = t.get_pandas_dataframe()

        heat_tot[i] = df.loc[0,["mean"]].values
    
    # Define the power factor
    heat_power = heat_power/8 #Power: 4kW 
    heat_power_origin = heat_tot.sum()
    k_power = heat_power/(heat_power_origin*1.6022e-19) # Energy unit: eV---> J
    
    heat_ave = heat_tot*1.6022e-19*k_power/volume_vec
    
    # Fitting for power radial distribution
    def func(r,a,b):
        return a/(-r**2+5.6**2)+b
    
    heat_ave_nor = heat_ave/1
    r_axe_nor = r_axe/1
    
    popt,pcov = curve_fit(func,r_axe_nor,heat_ave_nor)

    # The precision need to be edited in the future
    heat_func = '('+str(popt[0])+'/(-(x*x+y*y)+5.6*5.6)+'+str(popt[1])+')'

    domain_1 = '((x*x+y*y)>='+str(controlRod_r**2)+')'
    domain_2 = '(((x-5.5)*(x-5.5)+y*y)>='+str(heat_pipe_R**2)+')'
    domain_3 = '((x*x+y*y)<='+str(fuel_r**2)+')'

    # Temp = 1173.5 # Temperature(approximation), unit: K
    # lamb = (0.606+0.0351*Temp)*0.01
    # Thermal conductivity. unit: W/(K.cm)
    lamb = '(0.606+0.0351*'+ temp_func +')*0.01'

    force_func = domain_1+'*'+domain_2+'*'+domain_3+'*'+heat_func+'/(-'+lamb+')'
    
    return force_func


def editCellTemperature(fuel_temp,fuel_cell_ID_list):
    tree = ET.parse('geometry.xml')
    root = tree.getroot()

    k = 0
    for cell in root.iter('cell'):
        if cell.attrib['id']==str(fuel_cell_ID_list[k]):
            cell.attrib['temperature'] = str(fuel_temp[k])
            k = k+1
        else:
            continue
    tree.write('geometry.xml')
