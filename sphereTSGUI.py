# -*- coding: utf-8 -*-
"""
    Copyright 2014 Gavin Macaulay

    This file is part of SphereTS.

    SphereTS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SphereTS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SphereTS.  If not, see <http://www.gnu.org/licenses/>.
"""

# TODO:
#
# Package as a windows executable
#
# Compare results to Chu's program
#
# Provide a bulk TS calculation function
#
# Put in bitbucket and as a package on PyPI

from __future__ import division
print('Starting')
from sphereTS import calcWaterProperties, materialProperties, sphereTSFreqResponse

import matplotlib.pyplot as plt
import math

import numpy as np

from traits.api import HasTraits, Str, List, Bool, Range
from traitsui.api import View, Item, Group, Handler, CheckListEditor, CSVListEditor, EnumEditor, HTMLEditor
from traitsui.menu import Action, CancelButton, OKButton

class AboutDialog(HasTraits):
    """
    This class implements an About dialog box using the TraitsUI user interface
    framework. The text displayed in the dialog is read from a file containing 
    html-formatted text. 
    """
    
    def loadHelpText(self):
        """
        Loads the help text from the file.
        """
        file = open(self.helpFile, 'r')
        self.about_text = file.read()    
    
    helpFile = 'sphereTSGUIhelp.html'

    about_text = Str()
    
    view = View(Item('about_text', editor=HTMLEditor(), show_label=False),
                resizable=True, title='About',
                buttons=[OKButton])

class uiHandler(Handler):
    """
    This class handles user interface events.
    """
    def showAbout(self, info):
        """
        Reloads the help text and makes the About dialog box visible.
        """
        info.object.aboutDialog.loadHelpText()
        info.object.aboutDialog.edit_traits()
    
    def calculate(self, info):
        """
        Calculates the sphere target strength over a frequency range and at
        spot frequencies. 

        Displays the results using a Matplotlib figure window.        
        """
        fstart = info.object.freq_start*1e3 # [Hz]
        fstop = info.object.freq_end*1e3 # [Hz]
        bw = info.object.averaging_bandwidth*1e3 # [Hz]
        
        # Sort out the frequency step size for the TS calculations
        Ns = 10 # minimum number of sampes per bandwidth
        max_evaluations = 10000.0 # don't do more than this many sphere TS calculations (takes too long)

        fstep = min(bw/Ns, 1e3)

        if (fstop-fstart)/fstep > max_evaluations:
            fstep = (fstop-fstart)/max_evaluations
        
        params = {'fstart': max(1.0, (fstart-bw)),
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
        plt.gca().set_position((.1, .15, .8, .75))
        plt.plot(f/1e3, TS, linewidth=1.5, color='#1b9e77')
        plt.plot(f/1e3, TS_avg, linewidth=1.5, color='#d95f02')
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
            
            average_truncated = False
            if params['fstart'] < 0.0:
                params['fstart'] = 1
                average_truncated = True
                
            ff, TS = sphereTSFreqResponse(**params)
            avgTS = 10.0*np.log10(np.average(np.power(10.0, TS/10.0)))
            
            plt.plot(f, avgTS, 'o')
            spot_TS_text = spot_TS_text + '\n{:>8g}     {:>8.1f}'.format(f, avgTS)
            if average_truncated:
                spot_TS_text = spot_TS_text + '*'
        
        # Put a table of freq and TS values on the plot
        if len(freqs) > 0:
            ax = plt.gca()
            ylim = ax.get_ylim()
            xlim = ax.get_xlim()
            text_x_pos = (xlim[1] - xlim[0])*0.05 + xlim[0]
            text_y_pos = (ylim[1] - ylim[0])*0.05 + ylim[0]

            plt.text(text_x_pos, text_y_pos, spot_TS_text, 
                     verticalalignment='bottom', horizontalalignment='left',
                     bbox=dict(boxstyle='round,pad=0.5', 
                               facecolor='w', edgecolor='k'))
                    
        # Put the material properties on the plot too.
        material_text = '$\\rho = {:.1f} \/ kg/m^3$, $c = {:.1f} \/ m/s$, $\\rho_1 = {:.1f}$, $c_1 = {:.1f}$, $c_2 = {:.1f}$, $bw = {:.2f} \/ kHz$'\
        .format(params['rho'], params['c'], params['rho1'], \
                        params['c1'], params['c2'], bw/1e3)
        plt.figtext(0.02, 0.02, material_text)
        plt.draw()
    
    def object_sphere_material_changed(self, info):
        """
        Updates the sphere material variabls if the type of material is changed.
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
        
class sphereTSGUI(HasTraits):
    """
    Calculate and show the sphere TS using the TraitsUI framework.
    """
     
    default_material = 'Tungsten carbide'
    m = materialProperties()
    params = m[default_material]
    params['material'] = default_material
    params['a'] = 0.0381/2.0
    params['rho'] = 1026.2
    params['c'] = 1490.0
    
    spot_freqs = [12, 18, 38, 50, 70, 120, 200, 333, 420]
    spot_freqs = zip(spot_freqs, map(str, spot_freqs))
                   
    extra_spot_freqs = List(Range(low=0., exclude_low=True),desc='comma separated frequencies [kHz]',
                            label='Additional spot freqs [kHz]')
    
    sphere_material = Str(params['material'], label='Material')
    sphere_diameter = Range(low=0., value=params['a']*2*1000.0, exclude_low=True, label='Diameter [mm]')
    sphere_density = Range(low=0., value=params['rho1'], exclude_low=True, label='Density [kg/m^3]')
    sphere_c1 = Range(low=0., value=params['c1'], exclude_low=True, label='Longitudal sound speed [m/s]')
    sphere_c2 = Range(low=0., value=params['c2'], exclude_low=True, label='Transverse sound speed [m/s]')

    use_ctd = Bool(True)
    use_another_material = Bool(False, label='Another material')
    
    fluid_c = Range(low=0., value=params['c'], exclude_low=True, label='Sound speed in water [m/s]')
    fluid_density = Range(low=0., value=params['rho'], exclude_low=True, label='Density of water [kg/m^3]')
    
    fluid_temperature = Range(-2.0, 60, 10.0, label='Temperature [degC]')
    fluid_salinity = Range(0.0, 60.0, 35.0, label='Salinity [PSU]')
    fluid_depth = Range(0.0, 15000.0, 30.0, label='Depth [m]')

    freq_start = Range(low=0.0, value=12., label='Start frequency [kHz]')
    freq_end = Range(low=0.0, value=200., label='End frequency [kHz]')
    
    averaging_bandwidth = Range(low=0.1, value=2.5, label='Bandwidth for averaged TS [kHz]')

    CalculateButton = Action(name='Calculate', action='calculate')   
    AboutButton = Action(name='About', action='showAbout')

    aboutDialog = AboutDialog()

    view = View(Group(Group(
                        Item('sphere_diameter'),
                        Item('sphere_material', style='custom', 
                             enabled_when='not use_another_material',
                             editor=EnumEditor(values=m.keys(), cols=2)), 
                        Item('use_another_material'),
                        Item('sphere_density', enabled_when='use_another_material'),
                        Item('sphere_c1', enabled_when='use_another_material'),
                        Item('sphere_c2', enabled_when='use_another_material'), 
                        label='Sphere properties', show_border=True),
                        '10', # some extra space
                      Group(
                        Item('use_ctd', label='Calculate from T, S, and D'),
                        Item('fluid_temperature', enabled_when='use_ctd', style='text'),
                        Item('fluid_salinity', enabled_when='use_ctd', style='text'),
                        Item('fluid_depth', enabled_when='use_ctd', style='text'),
                        Item('fluid_c', enabled_when='not use_ctd', format_str='%.1f'),
                        Item('fluid_density', enabled_when='not use_ctd',format_str='%.1f'),
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
