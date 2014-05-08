# -*- coding: utf-8 -*-
"""
Created on Tue May 06 11:09:21 2014

@author: gavinj
"""

from __future__ import division

from sphereTS import calcWaterProperties, materialProperties, sphereTSFreqResponse

import matplotlib.pyplot as plt
import math

import numpy as np

from traits.api import HasTraits, Str, Float, List, Bool
from traitsui.api import View, Item, Group, Handler, CheckListEditor, CSVListEditor, EnumEditor, HTMLEditor
from traitsui.menu import Action, CancelButton, OKButton

class AboutDialog(HasTraits):
    about_text = Str("""
    <h2>Sphere TS calculator</h2>
    <p>This program calculates the TS of acoustic calibration spheres.</p>
    
    <p>Developed by Gavin Macaulay.</p>
    
    <h3>Instructions</h3>
    <li>Set the appropriate sphere and environmental properties
    <li>Choose the desired requencies and frequency sweep limits
    <li>Click on the <b>Calculate</b> button to produce a plot of the TS values
    
    <h3>Methods</h3>
    <p>Equations 6-9 of MacLennan (1981) are used to calculate the TS between
    the given start and stop frequencies. An averaged 
    TS is also calculated via a linear mean over the frequencies within
    the given bandwidth.</p>
    
    <p>In addition the bandwidth averaged TS at spot frequencies is calculated. 
    This is done independently of the frequency sweep.
    
    <h3>Acknowledgements</h3>
    The style and functioning of this program closely follows the Matlab code 
    developed by Dezchang Chu, NOAA.
    
    <h3>References</h3>
    <p>MacLennan, D.N., 1981. The Theory of Solid Spheres 
        as Sonar Calibration Targets. Scottish Fisheries Research 
        Report No. 22, Department of Agriculture and Fisheries for Scotland.</p> 
    """)
    view = View(Item('about_text', editor=HTMLEditor(), show_label=False),
                resizable=True, title='About',
                buttons=[OKButton])

class uiHandler(Handler):
    """
    """
    def showAbout(self, info):
        info.object.aboutDialog.edit_traits()
    
    def calculate(self, info):

        fstart = info.object.freq_start*1e3 # [Hz]
        fstop = info.object.freq_end*1e3 # [Hz]
        bw = info.object.averaging_bandwidth*1e3 # [Hz]
        
        # Sort out the frequency step size for the TS calculations
        Ns = 10 # minimum number of sampes per bandwidth
        max_evaluations = 10000.0 # don't do more than this many sphere TS calculations (takes too long)

        fstep = min(bw/Ns, 1e3)

        if (fstop-fstart)/fstep > max_evaluations:
            fstep = (fstop-fstart)/max_evaluations
        
        params = {'fstart': (fstart-bw),
                  'fstop': (fstop+bw),
                  'fstep': fstep,
                  'a': info.object.sphere_diameter/2000.0,
                  'rho1': info.object.sphere_density,
                  'c1': info.object.sphere_c1,
                  'c2': info.object.sphere_c2,
                  'c': info.object.fluid_c,
                  'rho': info.object.fluid_density}

        f, TS = sphereTSFreqResponse(**params)

        # Do a running mean of length N.
        N = round(bw/fstep)
        #### NEED TO DO this in the linear domain
        TS_avg = 10*np.log10(np.convolve(np.power(10.0, TS/10.0), np.ones((N,))/N, mode='same'))
        
        # Since we added a little to the frequency range above, to give valid
        # averaging out to the supplied frequency limits, we now trim the data
        # back to the requested limits
        N = math.floor(N/2.0)
        f = f[N:-N]
        TS = TS[N:-N]
        TS_avg = TS_avg[N:-N]
                  
        plt.ion()
        plt.figure()
        plt.plot(f/1e3, TS, linewidth=1.5)
        plt.plot(f/1e3, TS_avg, linewidth=1.5)
        plt.xlabel('Frequency (kHz)')
        plt.ylabel('TS (dB re 1 m$^2$)')
        plt.xlim(fstart/1e3, fstop/1e3)
        plt.grid()
        plt.legend(('Unaveraged TS', 'Averaged TS'), loc='lower right', 
                   fancybox=True, fontsize=12)

        if info.object.use_another_material:
            title_text = 'User defined sphere, {} mm'.format(params['a']*2000)
        else:
            title_text = '{} sphere, {} mm'.format(info.object.sphere_material, 
                                                   params['a']*2000)

        plt.title(title_text)
        
        # Then do the spot frequencies

        # Merge the two spot frequency lists, then convert to a set to keep
        # just the unique values, then back to a list for sorting
        freqs = list(set(info.object.spot_freqs+info.object.extra_spot_freqs))
        freqs.sort()

        # Calculate a set of TS around each spot frequency and then average
        spot_TS_text = 'f (kHz)     Avg. TS (dB)'
        for f in freqs:
            params['fstart'] = (f*1e3-bw/2)
            params['fstop'] = (f*1e3+bw/2)
            params['fstep'] = bw/20
            
            ff, TS = sphereTSFreqResponse(**params)
            avgTS = 10.0*np.log10(np.average(np.power(10.0, TS/10.0)))
            
            plt.plot(f, avgTS, 'o')
            spot_TS_text = spot_TS_text + '\n{:>8g}     {:>8.1f}'.format(f, avgTS)
        
        # Put a table of freq and TS values on the plot
        if len(freqs) > 0:
            ax = plt.gca()
            ylim = ax.get_ylim()
            xlim = ax.get_xlim()
            text_x_pos = (xlim[1] - xlim[0])*0.05 + xlim[0]
            text_y_pos = (ylim[1] - ylim[0])*0.05 + ylim[0]

            plt.text(text_x_pos, text_y_pos, spot_TS_text, 
                     verticalalignment='bottom',
                     horizontalalignment='left',
                     bbox=dict(boxstyle='round,pad=0.5', 
                               facecolor='w', 
                               edgecolor='k'))
        plt.draw()
    
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
            
    def object_use_another_material_changed(self, info):
        # when we switch back to specific material, make sure
        # to reset the displayed material properties
        if not info.object.use_another_material:
            self.object_sphere_material_changed(info)
                
    def close(self, info, is_ok):
        plt.close('all')
        return True
        
# TODO
#
# Limits on the input water temp, salinity, and depth numbers to sensible values
#
# Implement the help text
#
# Implement the averaging as defined in MacLennans's report, rather than a simple mean
#
# Tidy up the class/function structure
#
# Package as a windows executable
#
# Profile the code to speed up the TS calculations

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
    use_another_material = Bool(False)
    
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
    
    averaging_bandwidth = Float(2.5, 
                                label='Bandwidth for averaged TS [kHz]')

    CalculateButton = Action(name = 'Calculate', action = 'calculate')   
    AboutButton = Action(name = 'About', action = 'showAbout')

    aboutDialog = AboutDialog()

    view = View(Group(Group(
                        Item('sphere_diameter'),
                        Item('sphere_material', style='custom', 
                             enabled_when = 'not use_another_material',
                             editor=EnumEditor(values=m.keys(), cols=2)), 
                        Item('use_another_material', label='Another material'),
                        Item('sphere_density', enabled_when = 'use_another_material'),
                        Item('sphere_c1', enabled_when = 'use_another_material'),
                        Item('sphere_c2', enabled_when = 'use_another_material'), 
                        label='Sphere properties', show_border=True),
                        '10', # some extra space
                      Group(
                        Item('use_ctd', label='Calculate from T, S, and D'),
                        Item('fluid_temperature', enabled_when = 'use_ctd'),
                        Item('fluid_salinity', enabled_when = 'use_ctd'),
                        Item('fluid_depth', enabled_when = 'use_ctd'),
                        Item('fluid_c', enabled_when = 'not use_ctd', format_str='%.2f'),
                        Item('fluid_density', enabled_when = 'not use_ctd',format_str='%.2f'),
                        label='Environmental properties', show_border=True),
                        '10', # some extra space
                      Group(
                        Item('spot_freqs', style='custom',
                             label='Spot frequencies [kHz]', 
                             editor=CheckListEditor(values=spot_freqs, cols=3)),
                        Item('extra_spot_freqs', editor=CSVListEditor()),
                        Item('freq_start'),
                        Item('freq_end'),
                        Item('averaging_bandwidth'),
                        label='Frequencies', show_border=True)),
            resizable=True, 
            title='Sphere TS calculator',
            buttons = [AboutButton, CalculateButton, CancelButton],
            handler = uiHandler())


if __name__ == "__main__":
    ts = sphereTSGUI()
    ts.configure_traits()
