# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 10:49:32 2024

@author: sergi
"""

import meep as mp
import matplotlib.pyplot as plt
import pickle
import argparse
import numpy as np

# this is just to ensure that I can create a tapered waveguide. Let's begin with the one from
# the previous example and take it from there. It should prove somewhat trivial. 
# Remember that the tapering length is variable, but I supposed it was 14.288 because...

# remember that first shoudl be straight waveguide, then we compare to the other structure;   

# function for defining the source as it appears we only change the center of the source in any case
def define_source(src_point, center_freq, waveguide_width):
    sources = [mp.EigenModeSource(src=mp.GaussianSource(center_freq, 0.1),
                              center=src_point,
                              size=mp.Vector3(y=2*waveguide_width),
                              eig_parity=mp.EVEN_Z, # TE mode : EVEN Z (?)
                              eig_band=1, # band number of eigenmode, 1 should be fundamental mode
                              eig_match_freq=True)] # launched spatial mode specified by frequency True
    return sources

# function to obtain the E-Field after x timesteps.
efield_data = list()
def store_Efield(sim):
    efield_data.append(sim.get_array(component=mp.Ex, center=mp.Vector3(), size=sim.cell_size))

def main():
    parser = argparse.ArgumentParser(description='A script to perform MEEP simulation of 1x2 MMI Coupler')
    
    # define the arguments here
    parser.add_argument('structure', choices=['straight_waveguide', 'tapered_waveguide', 'min_waveguide', 'min_straight',
                                              'min_waveguide_steady'],
                        help='Structure choice: straight waveguide for normalization, or waveguide MMI coupler.')
    parser.add_argument('--plot', action='store_true', help='Generate plot of simulation cell.', default=False)
    parser.add_argument('--gif', action='store_true', help='Generate animation of simulation.', default=False)
    args = parser.parse_args()  

    a = 1.0 # micron
    resolution = 80 # pixels per micron
    wavelength0 = 1.55 # 1550 nm
    f_meep = 1.0 / wavelength0
    
    Lmmi = 11.5 # width of the MMI coupler
    Wmmi = 3.6 # length of the MMI coupler
    
    Wt = 1.2 # tapered waveguide major 
    Wco = 0.5 # waveguide width
    Hco = 0.22 # waveguide height
    Wg = 0.6 # arm separation for the output waveguides
    
    wvg_x = 0.5 # try to keep the waveguide short
    
    theta = 1.0 # taper angle
    Li = 14.22 # calculated taper length; begin counting where waveguide broadens
    Sli = 2.0 # straight waveguide length; the shorter it is, the better? In theory, for computational resources at least
    
    # materials
    SiO2 = mp.Medium(index=1.44)
    Si = mp.Medium(index=2.85) # this is waveguide material
    
    # boundary layers ; remember should capture at least one cycle of your wave
    dpml_x = 4.0 # 2.0
    dpml_y =  4.0 # 2.0
    
    # deciding the thickness of the material; this is in the y direction I believe 
    dBurriedOxide = 2.0 # 2.0 mm of burried oxide, or SiO2 ...
    
    # defining the cell size dimensions
    sx = dpml_x + Li + wvg_x + Lmmi + Li + wvg_x + dpml_x + 2*Sli # this should be it; can make this more elegant
    sy = dpml_y + dBurriedOxide + Wmmi + dpml_y
    
    cellSize = mp.Vector3(sx,sy,0)
    
    boundary_layers = [mp.PML(dpml_x,direction=mp.X),
                   mp.PML(dpml_y,direction=mp.Y)]
    
    # add source you tube
    src_point = mp.Vector3(-0.5*sx+dpml_x+Sli)
    sources = define_source(src_point, f_meep, Wco)
    
    gifMonComp = mp.Ex 
    fluxMonComp = mp.Ex

    if args.structure == 'straight_waveguide':
        
        # this is base, for normalization purposes and whatnot
        straight_waveguide = [mp.Block(center=mp.Vector3(), size=mp.Vector3(sx,2*Wco), material=Si)]
        
        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cellSize,
            geometry=straight_waveguide,
            boundary_layers=boundary_layers,
            sources=sources,
            default_material=SiO2)
        
        # adding monitors to the simulation
        nfreq = 200
        incidentFluxRegion = mp.FluxRegion(center=mp.Vector3(-0.5*sx+dpml_x+Sli+1.0,),
                                           size=mp.Vector3(y=2*Wt))
        incidentFluxMonitor = sim.add_flux(f_meep, 0, 1, incidentFluxRegion) # 0,1, fluxregion
        
        # incident run
        sim.run(until_after_sources=mp.stop_when_fields_decayed(50,fluxMonComp,mp.Vector3(0.5*sx-dpml_x-Sli),1e-5))
        
        # obtain incident fluxes
        res = sim.get_eigenmode_coefficients(incidentFluxMonitor,[1],eig_parity=mp.EVEN_Z)
        incident_coeffs = res.alpha
        incident_flux = mp.get_fluxes(incidentFluxMonitor)[0]
        incident_flux_data = sim.get_flux_data(incidentFluxMonitor)
        
        print(f"Obtained flux : {incident_flux}")
        print(f"res: {res}")
        print(f"incident coeffs: {incident_coeffs}")
        
        # maybe try to obtain transmission from this straight version as well
        sim.restart_fields()
        
        transmissionFluxRegion = mp.FluxRegion(center=mp.Vector3(0.5*sx-dpml_x-Sli),
                                               size=mp.Vector3(y=2*Wt))
        transmissionFluxMonitor = sim.add_flux(f_meep, 0, 1, transmissionFluxRegion)
        
        reflectionFluxRegion = incidentFluxRegion
        
        reflectionMonitor = sim.add_flux(f_meep, 0, 1, reflectionFluxRegion)
        
        sim.load_minus_flux_data(reflectionMonitor, incident_flux_data)
        
        sim.run(until_after_sources=mp.stop_when_fields_decayed(50,fluxMonComp,mp.Vector3(0.5*sx-dpml_x),1e-5))
        
        transmission_flux = mp.get_fluxes(transmissionFluxMonitor)[0]
        reflection_flux = mp.get_fluxes(reflectionMonitor)[0]
        
        print(f"Flux at transmission monitor: {transmission_flux}")
        print(f"T:{transmission_flux/incident_flux}")
        print(f"R:{-reflection_flux/incident_flux}")
        
        with open('neg_flux.pickle', 'wb') as f:
            pickle.dump(incident_flux_data, f)
            print("Saved flux data for normalization.")
        
        with open('incident_flux.pickle', 'wb') as f:
            pickle.dump(incident_flux, f)
        
        with open('res.pickle', 'wb') as f:
            pickle.dump(res, f)
            
        with open('incident_coeffs.pickle', 'wb') as f:
            pickle.dump(incident_coeffs, f)
        
        if args.plot:
            sim.plot2D()
            plt.savefig('straight_waveguide.png')
            
        if args.gif:
            sim.restart_fields()
            file_name = 'straight_waveguide_2D.gif'
            f = plt.figure()
            Animate = mp.Animate2D(fields=gifMonComp, f=f)
            sim.run(mp.at_every(0.5, Animate), 
                until_after_sources=mp.stop_when_fields_decayed(50,gifMonComp,mp.Vector3(0.5*sx-dpml_x),1e-5))
            fps = 50
            Animate.to_gif(fps, file_name)
            plt.close()
            
    elif args.structure == 'tapered_waveguide':
        # let's start with the vertices of the waveguide and the tapers
        vertices = [mp.Vector3(-sx,0.5*Wco),
                    mp.Vector3(-0.5*sx+dpml_x+Sli,0.5*Wco),
                    mp.Vector3(-0.5*Li,0.5*Wt),
                    mp.Vector3(-0.5*Li,-0.5*Wt),
                     mp.Vector3(-0.5*sx+dpml_x+Sli,-0.5*Wco),
                    mp.Vector3(-sx,-0.5*Wco)]
        
        vertices_out_1 = [mp.Vector3(-0.5*Li+Lmmi*0.5, 0.5*Wg+Wt),
                          mp.Vector3(-0.5*Li+Lmmi*0.5+Li, 0.5*Wg+Wt-Wco),
                          mp.Vector3(sx, 0.5*Wg+Wt-Wco),
                          mp.Vector3(sx, 0.5*Wg+Wt-2*Wco),
                          mp.Vector3(-0.5*Li+Lmmi*0.5+Li, 0.5*Wg+Wt-2*Wco),
                          mp.Vector3(-0.5*Li+Lmmi*0.5, 0.5*Wg)]
        
        vertices_out_2 = [mp.Vector3(-0.5*Li+Lmmi*0.5, -0.5*Wg-Wt),
                          mp.Vector3(-0.5*Li+Lmmi*0.5+Li, -0.5*Wg-Wt+Wco),
                          mp.Vector3(sx, -0.5*Wg-Wt+Wco),
                          mp.Vector3(sx, -0.5*Wg-Wt+2*Wco),
                          mp.Vector3(-0.5*Li+Lmmi*0.5+Li, -0.5*Wg-Wt+2*Wco),
                          mp.Vector3(-0.5*Li+Lmmi*0.5, -0.5*Wg)]
        
        geometry = [mp.Prism(vertices,height=mp.inf,material=Si),
                    mp.Block(center=mp.Vector3(x=-0.5*Li+Lmmi*0.5), material=Si, size=mp.Vector3(Lmmi,Wmmi)),
                    mp.Prism(vertices_out_1, height=mp.inf, material=Si),
                    mp.Prism(vertices_out_2, height=mp.inf, material=Si)]
        
        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cellSize,
            geometry=geometry,
            boundary_layers=boundary_layers,
            sources=sources,
            default_material=SiO2)
        
        # we need incident again
        incidentFluxRegion = mp.FluxRegion(center=mp.Vector3(-0.5*sx+dpml_x+Sli+1.0,),
                                           size=mp.Vector3(y=2*Wt))
        incidentFluxMonitor = sim.add_flux(f_meep, 0, 1, incidentFluxRegion) # 0,1, fluxregion
        
        # # add transmission monitors on either arm
        nfreq=200
        transmissionFluxRegion1 = mp.FluxRegion(center=mp.Vector3(-0.5*Li+Lmmi*0.5+Li, (0.5*Wg+Wt-Wco)*0.5),
                                                size=mp.Vector3(y=2*Wt))
        transmissionFluxRegion2 = mp.FluxRegion(center=mp.Vector3(-0.5*Li+Lmmi*0.5+Li, (-0.5*Wg-Wt+Wco)*0.5),
                                                size=mp.Vector3(y=2*Wt))
        
        transmissionFluxMonitor1 = sim.add_flux(f_meep, 0, 1, transmissionFluxRegion1)
        transmissionFluxMonitor2 = sim.add_flux(f_meep, 0, 1, transmissionFluxRegion2)
        
        # add reflection monitor (same place as incident monitor)
        reflectionMonitor = sim.add_flux(f_meep, 0, 1, incidentFluxRegion)
        
        # I already have flux data saved to binary file from straight waveguide simulation. Load that
        with open('neg_flux.pickle', 'rb') as f:
            incident_flux_data = pickle.load(f)
        with open('incident_flux.pickle', 'rb') as f:
            incident_flux = pickle.load(f)
        with open('res.pickle', 'rb') as f:
            res = pickle.load(f)
        with open('incident_coeffs.pickle', 'rb') as f:
            incident_coeffs = pickle.load(f)
            
        # load negative flux for reflection
        sim.load_minus_flux_data(reflectionMonitor,incident_flux_data)
        
        # run this jazz and get your cheese
        sim.run(until_after_sources=mp.stop_when_fields_decayed(50, fluxMonComp, 
                                                                mp.Vector3(-0.5*Li+Lmmi*0.5+Li, (0.5*Wg+Wt-Wco)*0.5), 1e-5))
        
        # the stuff
        res = sim.get_eigenmode_coefficients(transmissionFluxMonitor1,[1],eig_parity=mp.EVEN_Z)
        taper_coeffs = res.alpha
        transmissionFlux1 = mp.get_fluxes(transmissionFluxMonitor1)[0]
        transmissionFlux2 = mp.get_fluxes(transmissionFluxMonitor2)[0]
        reflectionFlux = mp.get_fluxes(reflectionMonitor)[0]
        
        R_coeffs = (abs(taper_coeffs[0,0,1])**2/abs(incident_coeffs[0,0,0])**2)
        R_flux = (-transmissionFlux1/incident_flux)
        print("refl:, {:.8f}, {:.8f}".format(R_coeffs,R_flux))
        
        Transmission = (transmissionFlux1/incident_flux) + (transmissionFlux2/incident_flux)
        Reflection = -reflectionFlux/incident_flux
        
        print(f"Transmission1: {transmissionFlux1}, Transmission2: {transmissionFlux2}, Reflection: {reflectionFlux}")
        print(f"Obtained reflection {Reflection} and transmission {Transmission}.")
        print(f"These values should equal 1.0, and they equal {Transmission+Reflection}")
        
        # save the jazz
        with open('Transmission_tapered.pickle', 'wb') as f:
            pickle.dump(Transmission, f)
            print("Saved transmission values.")
            
        with open('Reflection_tapered.pickle', 'wb') as f:
            pickle.dump(Reflection, f)
            print("Saved reflection data")
            
        if args.plot:
            sim.plot2D()
            plt.savefig('tapered_waveguide.png')
        
        if args.gif:
            sim.restart_fields()
            file_name = 'tapered_waveguide_2D.gif'
            f = plt.figure(dpi=200)
            Animate = mp.Animate2D(fields=gifMonComp, f=f)
            sim.run(mp.at_every(0.5, Animate), 
                    until_after_sources=mp.stop_when_fields_decayed(50,
                                                                    gifMonComp,mp.Vector3(-0.5*Li+Lmmi*0.5+Li, (0.5*Wg+Wt-Wco)*0.5),1e-5))
            fps = 50
            Animate.to_gif(fps, file_name)
            plt.close()

    elif args.structure == 'min_waveguide':
        # let's redo waveguide length to min size, basically just MMI
        minWgL = 1.0
        # since we redefine waveguide length, we must redefine cell size and our lengths in x direction
        sx = dpml_x + minWgL + Lmmi + minWgL + dpml_x # x is pretty much MMI. No Li needed
        sy = dpml_y + dBurriedOxide + Wmmi + dpml_y
        cellSize = mp.Vector3(sx,sy,0)
        
        # add source you tube
        src_point = mp.Vector3(-0.5*sx+dpml_x+0.25*minWgL)
        sources = define_source(src_point, f_meep, Wco)
        
        # now we use Geometry, which ought to be tapered_waveguide - the tapering...
        geometry = [mp.Block(center=mp.Vector3(x=-0.5*sx+0.5*dpml_x), material=Si, size=mp.Vector3(dpml_x,Wco)),
                    mp.Block(center=mp.Vector3(x=sx*-0.5+dpml_x+0.5*minWgL), material=Si, size=mp.Vector3(minWgL,Wco)),
                    mp.Block(center=mp.Vector3(), material=Si, size=mp.Vector3(Lmmi,Wmmi)),
                    mp.Block(center=mp.Vector3(x=0+Lmmi*0.5+0.5*minWgL,y=0.5*Wg+Wco*0.5), material=Si, size=mp.Vector3(minWgL,Wco)),
                    mp.Block(center=mp.Vector3(x=0+Lmmi*0.5+0.5*minWgL,y=-0.5*Wg-Wco*0.5), material=Si, size=(minWgL,Wco)),
                    mp.Block(center=mp.Vector3(x=0.5*sx-0.5*dpml_x,y=0.5*Wg+Wco*0.5), material=Si, size=mp.Vector3(dpml_x,Wco)),
                    mp.Block(center=mp.Vector3(x=0.5*sx-0.5*dpml_x,y=-0.5*Wg-Wco*0.5), material=Si, size=mp.Vector3(dpml_x,Wco))]
        
        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cellSize,
            geometry=geometry,
            boundary_layers=boundary_layers,
            sources=sources,
            default_material=SiO2)
        
        reflectionFluxRegion = mp.FluxRegion(center=mp.Vector3(x=sx*-0.5+dpml_x+0.5*minWgL),
                                             size=mp.Vector3(y=2*Wco))
        transmissionFluxRegion1 = mp.FluxRegion(center=mp.Vector3(Lmmi*0.5+0.5*minWgL,y=0.5*Wg+Wco*0.5),
                                                size=mp.Vector3(y=2*Wco))
        transmissionFluxRegion2 = mp.FluxRegion(center=mp.Vector3(Lmmi*0.5+0.5*minWgL,y=-0.5*Wg-Wco*0.5),
                                                size=mp.Vector3(y=2*Wco))
        
        reflectionMonitor = sim.add_flux(f_meep, 0, 1, reflectionFluxRegion)
        transmissionMonitor1 = sim.add_flux(f_meep, 0, 1, transmissionFluxRegion1)
        transmissionMonitor2 = sim.add_flux(f_meep, 0, 1, transmissionFluxRegion2)
        
        with open('neg_flux_min.pickle', 'rb') as f:
            incident_flux_data = pickle.load(f)
        with open('incident_flux_min.pickle', 'rb') as f:
            incident_flux = pickle.load(f)
        with open('incident_coeffs_min.pickle', 'rb') as f:
            incident_coeffs = pickle.load(f)
        
        sim.load_minus_flux_data(reflectionMonitor, incident_flux_data)
        
        # TODO: we have to obtain the Efield and plot it. Unsure which point in time has what we need
        # therefore get a bunch of timesteps and plot. FMLA for that perhaps. 
        sim.run(mp.at_every(1, store_Efield), until_after_sources=mp.stop_when_fields_decayed(50, fluxMonComp, 
                                                                mp.Vector3(Lmmi*0.5+0.5*minWgL,y=0.5*Wg+Wco*0.5), 1e-2))
        
        res1 = sim.get_eigenmode_coefficients(transmissionMonitor1,[1],eig_parity=mp.EVEN_Z)
        res2 = sim.get_eigenmode_coefficients(transmissionMonitor2, [1], eig_parity=mp.EVEN_Z)
        coeffs1 = res1.alpha
        coeffs2 = res2.alpha
        
        # chatgpt says to time-average the field here. We should have waited enought to have equilibrium...
        # this also means using a continuous source... hm...
        
        # perhaps check our Efields here
        efield = np.array(efield_data)
        
        steady_state = np.mean(efield, axis=0)
        plt.imshow(steady_state.transpose())
        plt.savefig(f"Efields.png")
        plt.close()
        transmissionFlux1 = mp.get_fluxes(transmissionMonitor1)[0]
        transmissionFlux2 = mp.get_fluxes(transmissionMonitor2)[0]
        reflectionFlux = mp.get_fluxes(reflectionMonitor)[0]
        
        R1_coeffs = (abs(coeffs1[0,0,1])**2/abs(incident_coeffs[0,0,0])**2)
        R1_flux = (-transmissionFlux1/incident_flux)
        
        totalTransmission = (transmissionFlux1 + transmissionFlux2) / incident_flux
        reflection = -reflectionFlux / incident_flux
        
        print(f"Transmissions: {transmissionFlux1},{transmissionFlux2}. Total: {totalTransmission}.")
        print(f"Reflection: {reflectionFlux}. Total reflection: {reflection}.")
        
        print(f"Our R1 coeffs are: {R1_coeffs}")
        print(f"Our R1_flux are: {R1_flux}")
        
        if args.plot:
            sim.plot2D()
            plt.savefig('min_waveguide.png')
            
        if args.gif:
            sim.restart_fields()
            file_name = 'min_waveguide_2D.gif'
            f = plt.figure(dpi=200)
            Animate = mp.Animate2D(fields=gifMonComp, f=f)
            sim.run(mp.at_every(0.5,Animate),
                    until_after_sources=mp.stop_when_fields_decayed(50, gifMonComp,mp.Vector3(x=0+Lmmi*0.5+0.5*minWgL,y=0.5*Wg+Wco*0.5),1e-2))
            fps = 50
            Animate.to_gif(fps, file_name)
            plt.close()
    
    elif args.structure == 'min_straight':
        # let's redo waveguide length to min size, basically just MMI
        minWgL = 1.0
        # since we redefine waveguide length, we must redefine cell size and our lengths in x direction
        sx = dpml_x + minWgL + Lmmi + minWgL + dpml_x # x is pretty much MMI. No Li needed
        sy = dpml_y + dBurriedOxide + Wmmi + dpml_y
        cellSize = mp.Vector3(sx,sy,0)
        
        geometry = [mp.Block(center=mp.Vector3(), material=Si, size=mp.Vector3(sx,Wco))]
        
        src_point = mp.Vector3(-0.5*sx+dpml_x+0.25*minWgL)
        sources = define_source(src_point, f_meep, Wco)
        
        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cellSize,
            geometry=geometry,
            boundary_layers=boundary_layers,
            sources=sources,
            default_material=SiO2)
        
        incidentFluxRegion = mp.FluxRegion(center=mp.Vector3(-0.5*sx+dpml_x+0.5*minWgL),
                                           size=mp.Vector3(y=2*Wco))
        
        incidentFluxMonitor = sim.add_flux(f_meep, 0, 1, incidentFluxRegion)
        
        sim.run(until_after_sources=mp.stop_when_fields_decayed(50, fluxMonComp,
                                                                mp.Vector3(-0.5*sx+dpml_x+0.5*minWgL), 1e-2))
        
        res = sim.get_eigenmode_coefficients(incidentFluxMonitor,[1],eig_parity=mp.EVEN_Z)
        incident_coeffs = res.alpha
        
        incidentFlux = mp.get_fluxes(incidentFluxMonitor)[0]
        incidentFluxData = sim.get_flux_data(incidentFluxMonitor)
        
        print(f"Obtained incident flux: {incidentFlux}.")
        
        with open('neg_flux_min.pickle', 'wb') as f:
            pickle.dump(incidentFluxData, f)
        with open('incident_flux_min.pickle', 'wb') as f:
            pickle.dump(incidentFlux, f)
        with open('incident_coeffs_min.pickle', 'wb') as f:
            pickle.dump(incident_coeffs, f)
        
        sim.restart_fields()
        
        transmissionFluxRegion = mp.FluxRegion(center=mp.Vector3(0.5*sx-dpml_x-0.5*minWgL),
                                               size=mp.Vector3(y=2*Wco))
        reflectionFluxRegion = incidentFluxRegion
        
        transmissionMonitor = sim.add_flux(f_meep, 0, 1, transmissionFluxRegion)
        reflectionMonitor = sim.add_flux(f_meep, 0, 1, reflectionFluxRegion)
        
        sim.load_minus_flux_data(reflectionMonitor, incidentFluxData)
        
        sim.run(until_after_sources=mp.stop_when_fields_decayed(50, gifMonComp,
                                                                mp.Vector3(0.5*sx-dpml_x-0.5*minWgL), 1e-2))
        
        transmission = mp.get_fluxes(transmissionMonitor)[0]
        reflection = mp.get_fluxes(reflectionMonitor)[0]
        
        T = transmission/incidentFlux
        R = -reflection/incidentFlux
        
        print(f'Obtained: transmission: {transmission}, reflection: {reflection}.')
        print(f'Normalized: {T},{R}')
        
        if args.plot:
            sim.plot2D()
            plt.savefig('min_straight_waveguide.png')
            
        if args.gif:
            sim.restart_fields()
            file_name = 'min_straight_waveguide_2D.gif'
            f = plt.figure()
            Animate = mp.Animate2D(fields=gifMonComp, f=f)
            sim.run(mp.at_every(0.5, Animate), 
                until_after_sources=mp.stop_when_fields_decayed(50,gifMonComp,
                                                                mp.Vector3(0.5*sx-dpml_x-0.5*minWgL),1e-2))
            fps = 50
            Animate.to_gif(fps, file_name)
            plt.close()
            
    elif args.structure == 'min_waveguide_steady':
        minWgL = 1.0
        # since we redefine waveguide length, we must redefine cell size and our lengths in x direction
        sx = dpml_x + minWgL + Lmmi + minWgL + dpml_x # x is pretty much MMI. No Li needed
        sy = dpml_y + dBurriedOxide + Wmmi + dpml_y
        cellSize = mp.Vector3(sx,sy,0)
        
        geometry = [mp.Block(center=mp.Vector3(x=-0.5*sx+0.5*dpml_x), material=Si, size=mp.Vector3(dpml_x,Wco)),
                    mp.Block(center=mp.Vector3(x=sx*-0.5+dpml_x+0.5*minWgL), material=Si, size=mp.Vector3(minWgL,Wco)),
                    mp.Block(center=mp.Vector3(), material=Si, size=mp.Vector3(Lmmi,Wmmi)),
                    mp.Block(center=mp.Vector3(x=0+Lmmi*0.5+0.5*minWgL,y=0.5*Wg+Wco*0.5), material=Si, size=mp.Vector3(minWgL,Wco)),
                    mp.Block(center=mp.Vector3(x=0+Lmmi*0.5+0.5*minWgL,y=-0.5*Wg-Wco*0.5), material=Si, size=(minWgL,Wco)),
                    mp.Block(center=mp.Vector3(x=0.5*sx-0.5*dpml_x,y=0.5*Wg+Wco*0.5), material=Si, size=mp.Vector3(dpml_x,Wco)),
                    mp.Block(center=mp.Vector3(x=0.5*sx-0.5*dpml_x,y=-0.5*Wg-Wco*0.5), material=Si, size=mp.Vector3(dpml_x,Wco))]
        
        src_point = mp.Vector3(-0.5*sx+dpml_x+0.25*minWgL)
        # try continuous source, probably still need a eigenmodesource
        # sources = [mp.EigenModeSource(src=mp.ContinuousSource(frequency=f_meep, width=2),
        #             center=src_point,
        #             size=mp.Vector3(y=2*Wco),
        #             eig_band=1,
        #             eig_parity=mp.EVEN_Z)]
        sources = [mp.Source(mp.ContinuousSource(frequency=f_meep, width=0.08),
                              component=mp.Ex,
                              center=src_point,
                              size=mp.Vector3(y=2*Wco))]
        
        sim = mp.Simulation(
            resolution=resolution,
            cell_size=cellSize,
            geometry=geometry,
            boundary_layers=boundary_layers,
            sources=sources,
            default_material=SiO2,
            force_complex_fields=True)
        
        # we only need to keep tabs of the efield with our function and then plot that jazz
        # now we can sample to obtain electric field data
        # sim.run(until=2000)
        
        # last_field = sim.get_array(center=mp.Vector3(), size=cellSize, component=mp.Ex)
        
        # efield = np.array(efield_data)
        # last_field = np.array(last_field)
        # meanfield = np.mean(efield, axis=0)
        
        # plt.figure()
        # plt.imshow((np.abs(last_field.T)), interpolation='spline36', cmap='RdBu')
        # plt.colorbar()
        # plt.savefig('LastEField.png')
        # plt.close()
        
        # we try the frequency domain solver which uses iterative methods (no run)
        num_tols = 5
        tols = np.power(10, np.arange(-8.0, -8.0-num_tols, -1.0))
        ex_dat = np.zeros((int(sx*80),int(sy*80),num_tols), dtype=np.complex_)
        
        for i in range(num_tols):
            sim.init_sim()
            sim.solve_cw(tols[i], 1000, 10)
            ex_dat[:,:,i] = sim.get_array(center=mp.Vector3(), size=mp.Vector3(int(sx*80),int(sy*80)), component=mp.Ex)
            
        err_dat = np.zeros(num_tols-1)
        for i in range(num_tols-1):
            err_dat[i] = np.linalg.norm(ex_dat[:,:,i]-ex_dat[:,:,num_tols-1])
        
        print("Done.")
        
        plt.figure()
        plt.loglog(tols[:num_tols-1], err_dat, 'bo-')
        plt.xlabel("Frequency Domain Solver Tolerance")
        plt.ylabel("L2 norm of error in fields")
        plt.savefig("FDSolver_tolerances.png")
        plt.close()
        
        eps_data = sim.get_array(center=mp.Vector3(), size=cellSize, component=mp.Dielectric)
        ex_data = np.absolute(ex_dat[:,:,num_tols-1])
        
        plt.figure()
        plt.imshow(eps_data.T, interpolation='spline36', cmap='binary')
        plt.imshow(ex_data.T, interpolation='spline36', cmap='Reds')
        plt.axis('off')
        plt.savefig('FDfields.png')
        plt.close()
        
        if args.gif:
            sim.restart_fields()
            file_name = 'min_waveguide_steady_2D.gif'
            f = plt.figure()
            Animate = mp.Animate2D(fields=gifMonComp, f=f)
            sim.run(mp.at_every(0.5, Animate), 
                until=100)
            fps = 50
            Animate.to_gif(fps, file_name)
            plt.close()
        
if __name__ == "__main__":
    main()