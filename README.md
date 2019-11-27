# Heat-Pipe-Reactor
# Instructions (version 0.0.3)

This is simulation code for heat pipe reactor (KRUSTY 1/8 model).

![photo](https://github.com/bearsan/Heat-Pipe-Reactor/cross_section.jpg)
![photo](https://github.com/bearsan/Heat-Pipe-Reactor/coupling_method.jpg)

HeatPipeReactor.xml and HeatPipeReactor_node.vtu are files for the mesh with about 2000 elements.
HeatPipeReactor_3.xml is the file for the mesh with about 8000 elements. Because the size of HeatPipeReactor_3_node.vtu is too big to upload, please contact with me if you want this file for test.


Completed:

    (a)Neutronics simulation 
    (b)Temperature simulation in solid area
    (c)3D-Modeling

To do:

    (a)Fluid calculation in heat pipe
    (b)Thermal expansion simulation
    (c)Burn-up calculation



1.Unit: length(cm),temperature(K). Others are SI units.

2.Define_OpenMC.py includes functions for operating OpenMC:

    (a)define_Geo_Mat_Set: This function is used to generate  geometry.xml, material.xml and settings.xml.
    
    (b)postProcess: This function is used to post-proccess heating source distribution and generate forcing function for Nektar++.
    
    (c)editCellTemperature: This function is used to edit temperature of cells in geometry.xml.
    
    (d)editForcefile: This function is used to edit 'Force.pts', which can deliver forcing data to Nektar++
    
    (e)getCellTemperature: Get temperature of cells from tempearature of nodes

3.Define_Nektar.py includes functions for operating Nektar++:

    (a)runNektar: Just as its name. Of cause, It can generate .fld and convert it to .vtu that we need.
    
    (b)postProcess: It can post-process results of temperature and generate temperature function in fuel area and coefficients of this function.
    
    (c)readNodesFromVtu: Read positions(x,y,z) of nodes from HeatPipeReactor_nodes.vtu.

4.HeatPipeReactor_nodes.vtu is only used to store the information of nodes. HeatPipeReactor.xml is the setting for Nektar++.

For more details, please check codes respectively.

You can get the parameters of KRUSTY and thermal properties of fuel from the following reference papers. 

# Reference:
[1] D. E. Burkes, C. A. Papesch, A. P. Maddison, T. Hartmann, and F. J. Rice. Thermo-physical properties of DU-10 wt. Journal of Nuclear Materials, 403:160 – 166, 2010.

[2] Leonardo de Holanda Mencarini, Jeffrey C. King. Fuel geometry options for a moderated low-enriched uranium kilowatt-class space nuclear reactor. Nuclear Engineering and Design, 340:122 - 132, 2018.

[3] D. I. POSTON, et. al. Design of the KRUSTY Reactor. Proceedings NETS-2018, ANS (2018).

[4] M. S. EL-GENK and J. P. TOURNIER. Uses  of Liquid-Metal and Water Heat Pipes in Space Reactor Power Systems. Front. Heat Pipes,2, 013002 (2011).

[5] Paul K. Romano, Nicholas E. Horelik, Bryan R. Herman, Adam G. Nelson, Benoit Forget, and Kord Smith, “OpenMC: A State-of-the-Art Monte Carlo Code for Research and Development,” Ann. Nucl. Energy, 82, 90–97 (2015).

[6] Paul K. Romano and Benoit Forget, “The OpenMC Monte Carlo Particle Transport Code,” Ann. Nucl. Energy, 51, 274–281 (2013).

[7] C.Cantwell, D.Moxey, A.Comerford, A.Bolis, G.Rocco, G.Mengaldo, D.DeGrazia, S.Yakovlev, J.E.Lombard, D.Ekelschot, B.Jordi, H.Xu, Y.Mohamied, C.Eskilsson, B.Nelson,P.Vos, C.Biotto, R.Kirby, S.Sherwin, Nektar++:an open source spectral/hp element framework, Comput.Phys.Commun, 192, 205–219 (2015).

