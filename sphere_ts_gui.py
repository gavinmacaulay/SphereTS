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
# Add unittests

from __future__ import division
from __future__ import print_function

import sphere_ts

import matplotlib.pyplot as plt
import math

import numpy as np

from traits.api import HasTraits, Str, List, Bool, Range
from traitsui.api import View, Item, Group, Handler
from traitsui.api import CheckListEditor, CSVListEditor, EnumEditor, HTMLEditor
from traitsui.menu import Action, CancelButton, OKButton

import TableFactory as tf

class AboutDialog(HasTraits):
    """
    This class implements an About dialog box using the TraitsUI user interface
    framework. The text displayed in the dialog is read from a file containing
    html-formatted text.
    """

    helpFile = 'about.html'
    about_text = Str()
    view = View(Item('about_text', editor=HTMLEditor(), show_label=False),
                resizable=True, title='About',
                buttons=[OKButton])

    def load_help_text(self):
        """
        Loads the help text from the file.
        """
        f = open(self.helpFile, 'r')
        self.about_text = f.read()

class EK60Dialog(HasTraits):
    """
    This class implements a dialog box that displays bandwidth averaged sphere
    TS estimates.
    """
    html_text = Str()

    view = View(Item('html_text',
                     editor=HTMLEditor(format_text=False), show_label=False), 
                     title='EK60',
                     buttons=[OKButton],
                     resizable=True)

class UIHandler(Handler):
    """
    This class handles user interface events.
    """
    def show_about(self, info):
        """
        Reloads the help text and makes the About dialog box visible.
        """
        info.object.aboutDialog.load_help_text()
        info.object.aboutDialog.edit_traits()
        
    def show_ek60_ts(self, info):
        """
        Calculates sphere TS estimates at specified frequencies over a range of
        sound speeds averaged over the receive bandwidths used by the Simrad
        EK60 echosounders. 
        
        The results are presented in table form in a separate dialog box.
        """

        # The dialog box displays HTML. We use the TableFactory import to do 
        # that, but the appearance is controlled by some simple CSS...
        htmlheader = """\
        <!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN">
        <html>
        <head>
        <title>Sample table</title>
        <style type="text/css">
         body { font-family: Helvetica,Arial,FreeSans; }
         table.reporttable { border-collapse: collapse;}
         table.reporttable td { color: #00000; border: 1px solid #98bf21; padding: 3px 7px 2px 7px;}
         table.reporttable .odd td { background-color: white; }
         table.reporttable .even td { background-color: #EAF2D3; }
         table.reporttable th { background-color: #A7C942; color: #ffffff; }
        </style>
        </head>
        <body>
        """
        
        htmlfooter = "</body></html>"
        
        tables, params, ek60_params = self.calculate_ek60_ts_table(info)

        if len(tables) == 0:
            html = '<p>You need to select some EK60 spot frequencies '\
                   'to get results to show here.</p>'
        else:
            html = ''

        for freq, table in sorted(tables.iteritems()):
            # A header row with the global column labels
            sup_header = tf.RowSpec(
                         tf.ColumnSpec('', '', span=1),
                         tf.ColumnSpec('', 'Pulse length (&mu;s)/bandwidth (kHz)', 
                                       span=len(table[0])-1))
            
            # A header row with the bandwidths
            bw_cols = [tf.ColumnSpec(x[1]) for x in ek60_params[freq]]
            suf_header = tf.RowSpec(tf.ColumnSpec('',''), *bw_cols)

            # The columns definiton for the actual data
            cols = [tf.ColumnSpec(x) for x in table[0]]
            ts_row = tf.RowSpec(*cols)

            # Convert TS numbers into formatted text 
            for row in table:
                itercells = row.iteritems()
                next(itercells) # Don't want to format the sound speed column
                for key, item in itercells:
                    row[key] = '{:.2f}'.format(item)

            lines = ts_row.makeall(table)
            a = params['a']*2000.0 # convert from radius in m to diameter in mm
            title_text = 'Sphere target strength at {} kHz'.format(freq)
            details_text = '<p>'\
                           '&empty; = {0} mm, '\
                           '&rho;<sub>w</sub> = {1[rho]:.1f} kg/m<sup>3</sup>'\
                           '&rho;<sub>s</sub> = {1[rho1]} kg/m<sup>3</sup>, '\
                           'c<sub>c</sub> = {1[c1]} m/s, '\
                           'c<sub>s</sub> = {1[c2]} m/s, '\
                           '</p>'.format(a, params)
            html = html + tf.HTMLTable(title_text,
                            details_text,
                            headers=[sup_header, ts_row, suf_header]).render(lines)
                            
                            
        html = htmlheader + html + htmlfooter

        info.object.ek60Dialog.html_text = html
        info.object.ek60Dialog.edit_traits()
        

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
        max_evaluations = 10000.0

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

        f, ts = sphere_ts.freq_response(**params)

        # Do a running mean of length N.
        N = round(bw/fstep)
        # Note that this convolution somehow gives results that are fstep/2 too 
        # high in frequency, so when they are plotted, we adjust for that.
        ts_avg = 10*np.log10(np.convolve(np.power(10.0, ts/10.0),
                                         np.ones((N,))/N, mode='same'))

        # Since we added a little to the frequency range above, to give valid
        # averaging out to the supplied frequency limits, we now trim the data
        # back to the requested limits
        N = int(math.floor(N/2.0))
        f = f[N:-N]
        ts = ts[N:-N]
        ts_avg = ts_avg[N:-N]

        plt.ion()
        plt.figure()
        plt.gca().set_position((.1, .15, .8, .75))
        plt.plot(f/1e3, ts, linewidth=1.5, color='#1b9e77')
        plt.plot((f-fstep/2.)/1e3, ts_avg, linewidth=1.5, color='#d95f02')
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
        spot_ts_text = 'f (kHz)     Avg. TS (dB)'
        for f in freqs:
            params['fstart'] = (f*1e3-bw/2)
            params['fstop'] = (f*1e3+bw/2)
            params['fstep'] = bw/20.0

            average_truncated = False
            if params['fstart'] < 0.0:
                params['fstart'] = 1
                average_truncated = True

            ff, ts = sphere_ts.freq_response(**params)
            avg_ts = 10.0*np.log10(np.average(np.power(10.0, ts/10.0)))

            plt.plot(f, avg_ts, 'o')
            spot_ts_text = spot_ts_text + \
                '\n{:>8g}     {:>8.1f}'.format(f, avg_ts)
            if average_truncated:
                spot_ts_text = spot_ts_text + '*'

        # Put a table of freq and TS values on the plot
        if len(freqs) > 0:
            ax = plt.gca()
            ylim = ax.get_ylim()
            xlim = ax.get_xlim()
            text_x_pos = (xlim[1] - xlim[0])*0.05 + xlim[0]
            text_y_pos = (ylim[1] - ylim[0])*0.05 + ylim[0]

            plt.text(text_x_pos, text_y_pos, spot_ts_text,
                     verticalalignment='bottom', horizontalalignment='left',
                     bbox=dict(boxstyle='round,pad=0.5',
                               facecolor='w', edgecolor='k'))

        # Put the material properties on the plot too.
        material_text = '$\\rho_w = {:.1f} \/ kg/m^3$, $c_w = {:.1f} \/ m/s$, '\
                        '$\\rho_s = {:.1f}$, $c_c = {:.1f}$, $c_s = {:.1f}$, '\
                        '$b_f = {:.2f} \/ kHz$'\
                        .format(params['rho'], params['c'], params['rho1'], \
                            params['c1'], params['c2'], bw/1e3)
        plt.figtext(0.02, 0.02, material_text)
        plt.draw()

    def calculate_ek60_ts_table(self, info):
        """
        Calculates bandwidth averaged sphere TS at spot frequencies. The
        bandwidths are the same as those used by the Simrad EK60 echosounder.
        """
        ek60_params = {}
        ek60_params[18] = [(512, 1.73), (1024, 1.56), (2048, 1.17), (4096, 0.71), (8192, 0.38)]
        ek60_params[38] = [(256, 3.675), (512, 3.275), (1024, 2.425), (2048, 1.448), (4096, 0.766)]
        ek60_params[70] = [(128, 6.74), (256, 6.09), (512, 4.63), (1024, 2.83), (2048, 1.51)]
        ek60_params[120] = [(64, 11.66), (128, 10.79), (256, 8.61), (512, 5.49), (1024, 2.99)]
        ek60_params[200] = [(64, 18.54), (128, 15.55), (256, 10.51), (512, 5.90), (1024, 3.05)]
        ek60_params[333] = [(64, 27.944), (128, 20.078), (256, 11.720), (512, 6.145), (1024, 3.112)]
            
        ss = range(1450, 1525, 5)

        params = {'a': info.object.sphere_diameter/2000.0,
                  'rho1': info.object.sphere_density,
                  'c1': info.object.sphere_c1,
                  'c2': info.object.sphere_c2,
                  'rho': info.object.fluid_density}

        # Only pass on the spot frequencies that we have bandwidths for
        spot_freqs = set(info.object.spot_freqs).intersection(ek60_params.keys())

        t = sphere_ts.calculate_ts_table(spot_freqs, ss, ek60_params, params)
        
        return t, params, ek60_params
        
    def object_sphere_material_changed(self, info):
        """
        Updates the sphere material variables if the type of material is
        changed.
        """
        m = sphere_ts.material_properties()

        s = m[info.object.sphere_material]
        info.object.sphere_density = s['rho1']
        info.object.sphere_c1 = s['c1']
        info.object.sphere_c2 = s['c2']

    def object_fluid_temperature_changed(self, info):
        """ Recalculates the fluid properties if the temperature changes."""
        self.update_fluid_properties(info)

    def object_fluid_salinity_changed(self, info):
        """ Recalculates the fluid properties if the salinity changes."""
        self.update_fluid_properties(info)

    def object_fluid_depth_changed(self, info):
        """ Recalculates the fluid properties if the depth changes."""
        self.update_fluid_properties(info)

    def update_fluid_properties(self, info):
        """ Recalculates the fluid properties."""
        c, rho = sphere_ts.water_properties(info.object.fluid_salinity,
                                            info.object.fluid_temperature,
                                            info.object.fluid_depth)
        info.object.fluid_c = c
        info.object.fluid_density = rho

    def object_use_ctd_changed(self, info):
        """ Recalculate the fluid properties when switching between the two
            options for getting the fluid density and sound speed."""
        if info.object.use_ctd:
            self.update_fluid_properties(info)

    def object_use_another_material_changed(self, info):
        """ Reset the material properties when switching back to using
            a specified material."""
        if not info.object.use_another_material:
            self.object_sphere_material_changed(info)

    def close(self, info, is_ok):
        plt.close('all')
        return True
   
class SphereTSGUI(HasTraits):
    """
    Calculate and show the sphere TS using the TraitsUI framework.
    """

    default_material = 'Tungsten carbide'
    m = sphere_ts.material_properties()
    params = m[default_material]
    params['material'] = default_material
    params['a'] = 0.0381/2.0
    params['rho'] = 1026.2
    params['c'] = 1490.0

    spot_freqs = [12, 18, 38, 50, 70, 120, 200, 333, 420]
    spot_freqs = list(zip(spot_freqs, list(map(str, spot_freqs))))

    extra_spot_freqs = List(Range(low=0., exclude_low=True),
                            desc='comma separated frequencies [kHz]',
                            label='Additional spot freqs [kHz]')

    sphere_material = Str(params['material'], label='Material')
    sphere_diameter = Range(low=0., value=params['a']*2*1000.0,
                            exclude_low=True, label='Diameter [mm]')
    sphere_density = Range(low=0., value=params['rho1'], exclude_low=True,
                           label='Density [kg/m^3]')
    sphere_c1 = Range(low=0., value=params['c1'], exclude_low=True,
                      label='Longitudal sound speed [m/s]')
    sphere_c2 = Range(low=0., value=params['c2'], exclude_low=True,
                      label='Transverse sound speed [m/s]')

    use_ctd = Bool(True)
    use_another_material = Bool(False, label='Another material')

    fluid_c = Range(low=0., value=params['c'], exclude_low=True,
                    label='Sound speed in water [m/s]')
    fluid_density = Range(low=0., value=params['rho'], exclude_low=True,
                          label='Density of water [kg/m^3]')

    fluid_temperature = Range(-2.0, 60, 10.0, label='Temperature [degC]')
    fluid_salinity = Range(0.0, 60.0, 35.0, label='Salinity [PSU]')
    fluid_depth = Range(0.0, 15000.0, 30.0, label='Depth [m]')

    freq_start = Range(low=0.0, value=12., label='Start frequency [kHz]')
    freq_end = Range(low=0.0, value=200., label='End frequency [kHz]')

    averaging_bandwidth = Range(low=0.1, value=2.5,
                                label='Bandwidth for averaged TS [kHz]')


    CalculateButton = Action(name='Calculate', action='calculate')
    AboutButton = Action(name='About', action='show_about')
    EK60Button = Action(name='EK60 tables', action='show_ek60_ts')

    aboutDialog = AboutDialog()
    ek60Dialog = EK60Dialog()

    view = View(
        Group(
            Group(
                Item('sphere_diameter'),
                Item('sphere_material', style='custom',
                     enabled_when='not use_another_material',
                     editor=EnumEditor(values=list(m.keys()), cols=2)),
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
                Item('fluid_density', enabled_when='not use_ctd',
                     format_str='%.1f'),
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
                label='Frequencies', show_border=True)
            ),
        resizable=True,
        title='Sphere TS calculator',
        buttons=[EK60Button, AboutButton, CalculateButton, CancelButton],
        handler=UIHandler())

if __name__ == "__main__":
    gui = SphereTSGUI()
    gui.configure_traits()
