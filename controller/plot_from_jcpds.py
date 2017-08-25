import pymatgen as mg
from pymatgen import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ds_jcpds

class plot_vlines(object):
	'''plot a vlines figure from jcpds and optionally cif file.'''

	def __init__(self):
		self.jcpds_filename = ''
		self.cif_filename = ''
		self.pressure = 0.
		self.wavelength = 0.
		self.diff_lines = 0.

	def set_parameters(self):
		self.pressure = float(input('>pressure:'))
		self.wavelength = float(input('>wavelength:'))

	def plot_vlines_jcpds(self):
		'''plot a vlines figure from jcpds file'''
		#choose which jcpds file to import
		self.jcpds_filename = input('>jcpds filename:')
		#set parameters
		self.set_parameters()
		#determine whether to campare with a cif file
		plot_cif = input('>compare with a cif file[y/n]:')
		#import jcpds file
		pattern = ds_jcpds.JCPDS(filename = self.jcpds_filename)
		#calculate the pattern
		pattern.cal_dsp(pressure = self.pressure)
		pattern.get_DiffractionLines()
		tth, inten = pattern.get_tthVSint(wavelength = self.wavelength)
		#plot and compare with a cif file
		if plot_cif == 'n':
			plt.vlines(tth, 0., inten, color = 'r')
			plt.xlabel('two theta')
			plt.ylabel('intensity')
		elif plot_cif == 'y':
			plt.vlines(tth, 0., inten, color = 'r', label = 'jcpds')
			plt.xlabel('two theta')
			plt.ylabel('intensity')
			self.plot_cif_data()
			plt.vlines(self.diff_lines[:,0], 0., self.diff_lines[:,2], \
				       color = 'b', label = 'cif')
		else:
			print('"compare with a cif file" input error.')
		plt.legend()
		plt.show()

	def plot_cif_data(self):
		'''plot a vlines figure from cif file'''
		xrange = (0,40)
        #choose which cif file to import
		self.cif_filename = input('>cif filename:')
		#import cif file:
		material = mg.Structure.from_file(self.cif_filename)
		#get diffraction pattern
		cal = XRDCalculator(wavelength = self.wavelength)
		pattern_cif = cal.get_xrd_data(material, two_theta_range = xrange)
		#extract two theta, d-sp, int, hkl
		d_lines = []
		for values in pattern_cif:
			hkl_key = values[2].keys()
			hkl_txt = str(hkl_key)[12:-3].split(",")
			d_lines.append([values[0], values[3], values[1], \
				            int(hkl_txt[0]), int(hkl_txt[1]), int(hkl_txt[-1])])
		self.diff_lines = np.asarray(d_lines)


