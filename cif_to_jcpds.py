import pymatgen as mg
from pymatgen import Lattice, Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ds_jcpds

def import_cif():
    """specify the cif file you want to convert"""
    global cif_filename
    cif_filename = input('>cif filename:')

def set_parameters():
    global k0, k0p, alphat, wl_xray, xrange
    k0 = float(input('>bulk modulus K0(in GPa):'))
    k0p = float(input('>change in K0 with pressure K0P:'))
    alphat = float(input('>thermal expansion coefficient ALPHAT:'))
    wl_xray = 0.3344
    xrange = (0,40)

def convert_to_jcpds():
    """convert parameters and diffraction pattern to jcpds file"""
    #read cif file
    material = mg.Structure.from_file(cif_filename)
    #setup an jcpds object from a cif file
    material_jcpds = ds_jcpds.JCPDS()
    material_jcpds.set_from_cif(cif_filename, k0, k0p, \
                      thermal_expansion=alphat, two_theta_range=xrange)
    #save to a JCPDS file
    jcpds_filename = input('>jcpds filename:')
    comments = input('>comments:')
    material_jcpds.write_to_file(jcpds_filename, comments=comments)





