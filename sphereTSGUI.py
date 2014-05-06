# -*- coding: utf-8 -*-
"""
Created on Tue May 06 11:09:21 2014

@author: gavinj
"""

from __future__ import division

from sphereTS import calcWaterProperties, materialProperties, sphereTS, sphereTSFreqResponse

import matplotlib.pyplot as plt

from traits.api import HasTraits, Str, Float, List, Bool
from traitsui.api import View, Item, Group, Handler, CheckListEditor, CSVListEditor, EnumEditor
from traitsui.menu import Action, CancelButton, HelpButton

class uiHandler(Handler):
    """
    """
    def calculate(self, info):
        material = info.object.sphere_material
        
        params = {'fstart': info.object.freq_start*1e3,
                  'fstop': info.object.freq_end*1e3,
                  'fstep': info.object.freq_step*1e3,
                  'a': info.object.sphere_diameter/2000.0,
                  'rho1': info.object.sphere_density,
                  'c1': info.object.sphere_c1,
                  'c2': info.object.sphere_c2,
                  'c': info.object.fluid_c,
                  'rho': info.object.fluid_density}

        f, TS = sphereTSFreqResponse(**params)
                                     
        plt.figure()
        plt.plot(f/1e3, TS)
        plt.xlabel('Frequency (kHz)')
        plt.ylabel('TS (dB re 1 m$^2$)')
        plt.grid()
        plt.title('{} sphere, {} mm'.format(material, params['a']*2*1000))
        
        # Then do the spot frequencies
        del params['fstart']
        del params['fstop']
        del params['fstep']
        # Merge the two spot frequency lists, then convert to a set to keep
        # just the unique values, then back to a list for sorting
        freqs = list(set(info.object.spot_freqs+info.object.extra_spot_freqs))
        freqs.sort()

        spot_TS_text = '  f (kHz)     TS (dB)'
        for f in freqs:
            TS = sphereTS(f*1000, **params)
            plt.plot(f, TS, 'o')
            spot_TS_text = spot_TS_text + '\n  {:>8g}     {:>8.1f}'.format(f, TS)
        
        ax = plt.gca()
        ylim = ax.get_ylim()
        xlim = ax.get_xlim()
        text_x_pos = (xlim[1] - xlim[0])*0.05 + xlim[0]
        text_y_pos = (ylim[1] - ylim[0])*0.05 + ylim[0]
        if len(freqs) > 0:
            plt.text(text_x_pos, text_y_pos, spot_TS_text, verticalalignment='bottom',
                     horizontalalignment='left',
                     bbox=dict(facecolor=(.7, .7, .7), alpha=0.5, edgecolor='None'))
        plt.show()
    
    def object_sphere_material_changed(self, info):
        """
        """
        m = materialProperties()
        s = m[info.object.sphere_material]
        info.object.sphere_density = s['rho1']
        info.object.sphere_c1 = s['c1']
        info.object.sphere_c2 = s['c2']
        
    def object_fluid_temperature_changed(self, info):
        self.updateFluidProperties(info)
        
    def object_fluid_salinity_changed(self, info):
        self.updateFluidProperties(info)
        pass

    def object_fluid_depth_changed(self, info):
        self.updateFluidProperties(info)
        pass
    
    def updateFluidProperties(self, info):
        c, rho = calcWaterProperties(info.object.fluid_salinity, 
                                     info.object.fluid_temperature,
                                     info.object.fluid_depth)
        info.object.fluid_c = c
        info.object.fluid_density = rho
        
    def object_use_ctd_changed(self, info):
        if info.object.use_ctd:
            self.updateFluidProperties(info)
            
    def close(self, info, is_ok):
        plt.close('all')
        return True
        
# TODO
# when the sphere properties are changed, deselect the selected material properties and change the plot title
#
# Limits on the input water temp, salinity, and depth numbers to sensible values
#
# Have multiple graph windows
#
# Print out the TS at spot frequencies
#
# Allow the specification of bandwidth for TS averaging
#
# Implement TS averaging over a bandwidth
#
# Implement the help text
#
# Tidy up the class/function structure
#
# Package as a windows executable
#
# Profile the code to speed up the TS calculations
#
# Add a waiting cursor when doing the calculations

class sphereTSGUI(HasTraits):
    """
    Calculate and show the sphere TS using a GUI
    """

       
    default_material = 'Tungsten carbide'
    m = materialProperties()
    params = m[default_material]
    params['material'] = default_material
    params['a'] = 0.0381/2.0
    params['rho'] = 1026.2
    params['c'] = 1490.0
    
    # Find a more clever way to do this...
    spot_freqs = [(12.0, '12'), (18.0, '18'), (38.0, '38'), (50.0, '50'), 
                   (70.0, '70'), (120.0, '120'), (200.0, '200'),
                   (333.0, '333'), (420.0, '420')]
                   
    extra_spot_freqs = List(Float,desc='comma separated frequencies [kHz]',
                            label='Additional spot freqs [kHz]')
    
    sphere_material = Str(params['material'], label='Material')
    sphere_diameter = Float(params['a']*2*1000, 
                            desc='sphere diameter [mm]',
                            label='Diameter [mm]')
    sphere_density = Float(params['rho1'], 
                           desc='sphere material density [kg/m^3]',
                           label='Density [kg/m^3]')
    sphere_c1 = Float(params['c1'], 
                      desc='longitudinal sound speed of the sphere material [m/s]',
                      label='Longitudal sound speed [m/s]')
    sphere_c2 = Float(params['c2'], 
                      desc='transverse sound speed of the sphere material [m/s]',
                      label='Transverse sound speed [m/s]')

    use_ctd = Bool(True)
    
    fluid_c = Float(params['c'],
                    desc='Sound speed in water [m/s]',
                    label='Sound speed in water [m/s]')
    fluid_density = Float(params['rho'],
                           desc='Density of water [kg/m^3]',
                           label='Density of water [kg/m^3]')
    
    fluid_temperature = Float(10.0, label='Temperature [degC]')
    fluid_salinity = Float(35.0, label='Salinity [PSU]')
    fluid_depth = Float(30.0, label='Depth [m]')

    freq_start = Float(12., label='Start frequency [kHz]')
    freq_end = Float(200., label='End frequency [kHz]')
    freq_step = Float(1, label='Step frequency [kHz]')

    CalculateButton = Action(name = 'Calculate', action = 'calculate')    

    view = View(Group(Group(
                        Item('sphere_material', style='custom', 
                             editor=EnumEditor(values=m.keys(), cols=2)), 
                        Item('sphere_diameter'),
                        Item('sphere_density'),
                        Item('sphere_c1'),
                        Item('sphere_c2'), 
                        label='Sphere properties', show_border=True),
                        '10', # some extra space
                      Group(
                        Item('use_ctd', label='Calculate from T, S, and D'),
                        Item('fluid_temperature', enabled_when = 'use_ctd'),
                        Item('fluid_salinity', enabled_when = 'use_ctd'),
                        Item('fluid_depth', enabled_when = 'use_ctd'),
                        Item('fluid_c', enabled_when = 'not use_ctd', format_str='%.2f'),
                        Item('fluid_density', enabled_when = 'not use_ctd',format_str='%.2f'),
                        label='Environment properties', show_border=True),
                        '10', # some extra space
                      Group(
                        Item('spot_freqs', style='custom',
                             label='Spot frequencies [kHz]', 
                             editor=CheckListEditor(values=spot_freqs, cols=3)),
                        Item('extra_spot_freqs', editor=CSVListEditor()),
                        Item('freq_start'),
                        Item('freq_end'),
                        Item('freq_step'),
                        label='Frequencies', show_border=True)),
            resizable=True, 
            title='Sphere TS calculator',
            buttons = [CalculateButton, CancelButton, HelpButton],
            handler = uiHandler())

ts = sphereTSGUI()
ts.configure_traits()


def testCall():
    """
    
    """
    
    #params = [0.0381/2, 1490, 6853, 4171, 1026.2, 14900]
    #f, TS = sphereTSFreqResponse(12e3,200e3, *params, fstep=200)

    m = ts.materialProperties()
    params = m['Tungsten carbide']
    
    params['a'] = 0.0381/2.0
    params['rho'] = 1026.2
    params['c'] = 1490.0
    
    f, TS = sphereTSFreqResponse(12e3, 200e3, fstep=200, **params)    
    
    plt.clf()
    plt.plot(f/1e3, TS)
    plt.xlabel('Frequency (kHz)')
    plt.ylabel('TS (dB re 1 m$^2$)')
    plt.show()