import sys
import os

import pymatgen as mg
from pymatgen import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from ds_jcpds import *

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

class Main(QMainWindow):
	"""plot a vlines figure from jcpds and optionally cif file."""
	def __init__(self):
		super().__init__()
		self.pressurej = 0.
		self.wavelengthj = 0.
		self.wavelengthc = 0.
		self.twothetamin = 0
		self.twothetamin = 0
		self.k0 = 0.
		self.k0p = 0.
		self.alphat = 0.
		self.diff_lines = 0.
		self.jcpds_filename = ''
		self.cif_filename = ''
		self.t_filename = ''
		self.tf_filename = ''
		self.comments = ''
		self.initUI()

	def initUI(self):
		self.resize(800, 600)
		self.center()

		'''
		self.figure = plt.figure()
		self.drawing = self.figure.add_subplot(111)
		self.canvas = FigureCanvas(self.figure)

		self.setCentralWidget(self.canvas)
		'''
		self.tabWidget = QTabWidget()
		self.tabWidget.setTabPosition(QTabWidget.West)

		self.infoj = QWidget()
		self.infoc = QWidget()
		self.infot = QWidget()
		self.plot = QWidget()
		self.tabWidget.addTab(self.infoj, 'info jcpds')
		self.tabWidget.addTab(self.infoc, 'info cif')
		self.tabWidget.addTab(self.infot, 'info jcpds_from_cif')
		self.tabWidget.addTab(self.plot, 'plot')

		self.jlabel = QLabel('original information of jcpds file')
		self.clabel = QLabel('original information of cif file')
		self.tlabel = QLabel('information of jcpds transferred from cif')
		self.jedit = QTextEdit()
		self.cedit = QTextEdit()
		self.tedit = QTextEdit()

		self.tbutton = QPushButton(QIcon('./icons/save.png'), 'Save as jcpds', self)
		self.tbutton.clicked.connect(self.savejcpds)

		self.vboxj = QVBoxLayout()
		self.vboxj.addWidget(self.jlabel)
		self.vboxj.addWidget(self.jedit)
		self.infoj.setLayout(self.vboxj)

		self.vboxc = QVBoxLayout()
		self.vboxc.addWidget(self.clabel)
		self.vboxc.addWidget(self.cedit)
		self.infoc.setLayout(self.vboxc)

		self.gridt = QGridLayout()
		self.gridt.addWidget(self.tlabel, 0, 0, 1, 1)
		self.gridt.addWidget(self.tbutton, 0, 3, 1, 1)
		self.gridt.addWidget(self.tedit, 1, 0, 1, 4)
		self.infot.setLayout(self.gridt)

		self.setCentralWidget(self.tabWidget)

		impjAct = QAction(QIcon('./icons/importj.png'), 'Import jcpds', self)
		impjAct.setShortcut('Ctrl+J')
		impjAct.triggered.connect(self.loadandreadj)

		impcAct = QAction(QIcon('./icons/importc.png'), 'Import cif', self)
		impcAct.setShortcut('Ctrl+F')
		impcAct.triggered.connect(self.loadandreadc)

		trsAct = QAction(QIcon('./icons/transfer.png'), 'cif to jcpds', self)
		trsAct.setShortcut('Ctrl+T')
		trsAct.triggered.connect(self.ciftojcpds)

		plotjAct = QAction(QIcon('./icons/plotj.png'), 'Plot jcpds', self)
		plotjAct.setShortcut('Ctrl+P')
		plotjAct.triggered.connect(self.plot_jcpds)

		plotcAct = QAction(QIcon('./icons/plotc.png'), 'Plot cif', self)
		plotcAct.setShortcut('Ctrl+O')
		plotcAct.triggered.connect(self.plot_cif)

		exitAct = QAction(QIcon('./icons/exit.png'), 'Exit', self)
		exitAct.setShortcut('Ctrl+Q')
		exitAct.triggered.connect(qApp.quit)

		tbimpj = self.addToolBar('Import jcpds')
		tbimpj.addAction(impjAct)

		tbimpc = self.addToolBar('Import cif')
		tbimpc.addAction(impcAct)

		tbtrs = self.addToolBar('cif to jcpds')
		tbtrs.addAction(trsAct)

		tbplotj = self.addToolBar('Plot jcpds')
		tbplotj.addAction(plotjAct)

		tbplotc = self.addToolBar('Plot cif')
		tbplotc.addAction(plotcAct)

		tbexit = self.addToolBar('Exit')
		tbexit.addAction(exitAct)

		menubar = self.menuBar()

		mnfile = menubar.addMenu('File')
		mnfile.addAction(impjAct)
		mnfile.addAction(impcAct)

		mnfc = menubar.addMenu('Functions')
		mnfc.addAction(trsAct)
		mnfc.addAction(plotjAct)
		mnfc.addAction(plotcAct)


		"""
		widget = QWidget()
		self.setCentralWidget(widget)

		self.vbox = QHBoxLayout()

		'''just for test
		self.jcpds = QLabel('jcpds')
		self.jcpdsEdit = QTextEdit()
		self.vbox.addWidget(self.jcpds)
		self.vbox.addWidget(self.jcpdsEdit)
		'''

		self.layout = QGridLayout()
		self.layout.addLayout(self.vbox, 0, 0)

		widget.setLayout(self.layout)
		"""

		self.setWindowTitle('plot jcpds demo')
		self.show()

	def center(self):
		qr = self.frameGeometry()
		cp = QDesktopWidget().availableGeometry().center()
		qr.moveCenter(cp)
		self.move(qr.topLeft())

	def loadandreadj(self):
		#choose which jcpds file to import and display information
		self.tabWidget.setCurrentWidget(self.infoj)
		fname = QFileDialog.getOpenFileName(self, 'Import jcpds', '/home', '*.jcpds')
		self.jcpds_filename = fname[0]

		self.jedit.setText('From: ' + self.jcpds_filename + '\n')

		crystal = JCPDS_extend(filename = self.jcpds_filename)
		crystal.read_file(self.jcpds_filename)

		if crystal.version_status == 'old':
			self.jedit.append('**old version jcpds file.**' + '\n')

		elif crystal.version_status == 'new':
			self.jedit.append(5*'-' + 'crystal information' + 5*'-')
			self.jedit.append('K0: {:.3f}'.format(crystal.k0))
			self.jedit.append('K0P: {:.5f}'.format(crystal.k0p))
			self.jedit.append('SYMMETRY: ' + crystal.symmetry.upper())
			self.jedit.append('A: {:.4f}'.format(crystal.a))
			self.jedit.append('B: {:.4f}'.format(crystal.b))
			self.jedit.append('C: {:.4f}'.format(crystal.c))
			self.jedit.append('ALPHA: {:.1f}'.format(crystal.alpha))
			self.jedit.append('BETA: {:.1f}'.format(crystal.beta))
			self.jedit.append('GAMMA: {:.1f}'.format(crystal.gamma))
			self.jedit.append('VOLUME: {:.4f}'.format(crystal.v))
			self.jedit.append('ALPHAT: ' + str(crystal.thermal_expansion))
			self.jedit.append('DIHKL: see below' + '\n')

		self.jedit.append(5*'-' + 'original information' + 5*'-')
		jcpdsfile = open(self.jcpds_filename, 'r')
		jline = jcpdsfile.read()
		self.jedit.append(jline)

	def loadandreadc(self):
		#choose which cif file to import and display information
		self.tabWidget.setCurrentWidget(self.infoc)
		fname = QFileDialog.getOpenFileName(self, 'Import cif', '/home', '*.cif')
		self.cif_filename = fname[0]

		self.cedit.setText('From: ' + self.cif_filename)
		self.cedit.append('\n')
		self.cedit.append(5*'-' + 'original information' + 5*'-')

		ciffile = open(self.cif_filename, 'r')
		cline = ciffile.read()
		self.cedit.append(cline)

	def setparametersc(self):
		setdialog = Paradialogc()
		if setdialog.exec_():
			self.twothetamin = float(setdialog.ttmin())
			self.twothetamax = float(setdialog.ttmax())
			self.k0 = float(setdialog.k0())
			self.k0p = float(setdialog.k0p())
			self.alphat = float(setdialog.alphat())
			self.wavelengthc = float(setdialog.wl())
			self.comments = setdialog.comment()

		setdialog.destroy()

	def ciftojcpds(self):
		'''
		import a cif file, transfer it to jcpds, display information and 
		optionally save jcpds file
		'''
		self.tabWidget.setCurrentWidget(self.infot)
		fname = QFileDialog.getOpenFileName(self, 'Import cif', '/home', '*.cif')
		self.t_filename = fname[0]

		self.setparametersc()
		ttrange = (self.twothetamin, self.twothetamax)

		self.tedit.setText('From: ' + self.t_filename + '\n')
		self.tedit.append(5*'-' + 'transferred jcpds information' + 5*'-')

		tfile = JCPDS_extend()
		tfile.set_from_cif(self.t_filename, self.k0, self.k0p,
				thermal_expansion = self.alphat, two_theta_range = ttrange)

		self.tedit.append('VERSION: {:d}'.format(tfile.version))
		self.tedit.append('COMMENT: ' + self.comments)
		self.tedit.append('K0: {:.3f}'.format(self.k0))
		self.tedit.append('K0P: {:.5f}'.format(self.k0p))
		self.tedit.append('SYMMETRY: ' + tfile.symmetry.upper())

		# 1 cubic, 2 hexagonal, 3 tetragonal, 4 orthorhombic
		# 5 monoclinic, 6 triclinic, 7 manual P, d-sp input

		if tfile.symmetry == 'cubic':  # cubic
			self.tedit.append('A: {:.4f}'.format(tfile.a0))
		elif tfile.symmetry == 'manual':  # P, d-sp input
			self.tedit.append('A: {:.4f}'.format(tfile.a0))
		elif tfile.symmetry == 'hexagonal' or tfile.symmetry == 'trigonal':
			self.tedit.append('A: {:.4f}'.format(tfile.a0))
			self.tedit.append('C: {:.4f}'.format(tfile.c0))
		elif tfile.symmetry == 'tetragonal':  # tetragonal
			self.tedit.append('A: {:.4f}'.format(tfile.a0))
			self.tedit.append('C: {:.4f}'.format(tfile.c0))
		elif tfile.symmetry == 'orthorhombic':  # orthorhombic
			self.tedit.append('A: {:.4f}'.format(tfile.a0))
			self.tedit.append('B: {:.4f}'.format(tfile.b0))
			self.tedit.append('C: {:.4f}'.format(tfile.c0))
		elif tfile.symmetry == 'monoclinic':  # monoclinic
			self.tedit.append('A: {:.4f}'.format(tfile.a0))
			self.tedit.append('B: {:.4f}'.format(tfile.b0))
			self.tedit.append('C: {:.4f}'.format(tfile.c0))
			self.tedit.append('BETA: {:.1f}'.format(tfile.beta0))
		elif tfile.symmetry == 'triclinic':  # triclinic
			self.tedit.append('A: {:.4f}'.format(tfile.a0))
			self.tedit.append('B: {:.4f}'.format(tfile.b0))
			self.tedit.append('C: {:.4f}'.format(tfile.c0))
			self.tedit.append('ALPHA: {:.1f}'.format(tfile.alpha0))
			self.tedit.append('BETA: {:.1f}'.format(tfile.beta0))
			self.tedit.append('GAMMA: {:.1f}'.format(tfile.gamma0))

		self.tedit.append('ALPHAT: ' + str(tfile.thermal_expansion))

		for line in tfile.DiffLines:
			self.tedit.append("DIHKL: {0:.4f} {1:.2f} {2:.2f} {3:.2f} {4:.2f}".format(
				line.dsp0, line.intensity, line.h, line.k, line.l))

	def savejcpds(self):
		#save transferred jcpds file
		fname = QFileDialog.getSaveFileName(self, 'Save as jcpds', '/home', '*.jcpds')
		self.tf_filename = fname[0]

		tffile = open(self.tf_filename, 'w')
		tftextlines = self.tedit.toPlainText().splitlines()
		for line in tftextlines[3:]:
			tffile.write(line + '\n')

		tffile.close()

	def plot_jcpds(self):
		#import jcpds file
		pattern = JCPDS_extend(filename = self.jcpds_filename)

		#calculate the pattern
		pattern.cal_dsp(pressure = self.pressure)
		pattern.get_DiffractionLines()
		tth, inten = pattern.get_tthVSint(wavelength = self.wavelength)

		#plot
		self.drawing.hold (False)
		self.drawing.plot(tth, inten, color = 'r', label = 'jcpds')
		self.canvas.draw()

	def plot_cif(self):
		pass

class JCPDS_extend(JCPDS):
	def __init__(self, filename=None):
		if filename is None:
			self.file = ' '
			self.name = ' '
			self.version = 0
			self.comments = ''
			self.symmetry = ''
			self.k0 = 0.
			self.k0p = 0.  # k0p at 298K
			self.dk0dt = 0.
			self.dk0pdt = 0.
			self.thermal_expansion = 0.  # alphat at 298K
			self.thermal_expansion_dt = 0.
			self.a0 = 0.
			self.b0 = 0.
			self.c0 = 0.
			self.alpha0 = 0.
			self.beta0 = 0.
			self.gamma0 = 0.
			self.v0 = 0.
			self.DiffLines = []
		else:
			self.file = filename
			self.dk0dt = 0.
			self.dk0pdt = 0.
			self.thermal_expansion_dt = 0.
			self.b0 = 0.
			self.c0 = 0.
			self.alpha0 = 0.
			self.beta0 = 0.
			self.gamma0 = 0.
			self.v0 = 0.
			self.read_file(self.file)
		self.a = 0.
		self.b = 0.
		self.c = 0.
		self.alpha = 0.
		self.beta = 0.
		self.gamma = 0.
		self.v = 0.

		self.version_status = ''

	def read_file(self, file):
		"""
		read a jcpds file
		"""
		self.file = file
		# Construct base name = file without path and without extension
		name = os.path.splitext(os.path.basename(self.file))[0]
		self.name = name
#		line = '', nd=0
		version = 0.
		self.comments = []
		self.DiffLines = []

		inp = open(file, 'r').readlines()
#		my_list = [] # get all the text first and throw into my_list

		if inp[0][0] in ('2', '3', '4'):
			version = int(inp[0])  # JCPDS version number
			self.version = version
			header = inp[1]  # header
			self.comments = header

			item = str.split(inp[2])
			crystal_system = int(item[0])
			if crystal_system == 1:
				self.symmetry = 'cubic'
			elif crystal_system == 2:
				self.symmetry = 'hexagonal'
			elif crystal_system == 3:
				self.symmetry = 'tetragonal'
			elif crystal_system == 4:
				self.symmetry = 'orthorhombic'
			elif crystal_system == 5:
				self.symmetry = 'monoclinic'
			elif crystal_system == 6:
				self.symmetry = 'triclinic'
			elif crystal_system == 7:
				self.symmetry = 'manual'
			# 1 cubic, 2 hexagonal, 3 tetragonal, 4 orthorhombic
			# 5 monoclinic, 6 triclinic, 7 manual P, d-sp input

			k0 = float(item[1])
			k0p = float(item[2])
			self.k0 = k0
			self.k0p = k0p

			item = str.split(inp[3])  # line for unit-cell parameters

			if crystal_system == 1:  # cubic
				a = float(item[0])
				b = a
				c = a
				alpha = 90.
				beta = 90.
				gamma = 90.
			elif crystal_system == 7:  # P, d-sp input
				a = float(item[0])
				b = a
				c = a
				alpha = 90.
				beta = 90.
				gamma = 90.
			elif crystal_system == 2:  # hexagonal
				a = float(item[0])
				c = float(item[1])
				b = a
				alpha = 90.
				beta = 90.
				gamma = 120.
			elif crystal_system == 3:  # tetragonal
				a = float(item[0])
				c = float(item[1])
				b = a
				alpha = 90.
				beta = 90.
				gamma = 90.
			elif crystal_system == 4:  # orthorhombic
				a = float(item[0])
				b = float(item[1])
				c = float(item[2])
				alpha = 90.
				beta = 90.
				gamma = 90.
			elif crystal_system == 5:  # monoclinic
				a = float(item[0])
				b = float(item[1])
				c = float(item[2])
				beta = float(item[3])
				alpha = 90.
				gamma = 90.
			elif crystal_system == 6:  # triclinic
				a = float(item[0])
				b = float(item[1])
				c = float(item[2])
				alpha = float(item[3])
				beta = float(item[4])
				gamma = float(item[5])

			self.a0 = a
			self.b0 = b
			self.c0 = c
			self.alpha0 = alpha
			self.beta0 = beta
			self.gamma0 = gamma

			item = str.split(inp[4])

			if self.version == 3:
				thermal_expansion = 0.
			else:
				thermal_expansion = float(item[0])
			self.thermal_expansion = thermal_expansion

			for line in inp[6:]:
				item = str.split(line)
				DiffLine = DiffractionLine()
				DiffLine.dsp0 = float(item[0])
				DiffLine.intensity = float(item[1])
				DiffLine.h = float(item[2])
				DiffLine.k = float(item[3])
				DiffLine.l = float(item[4])
				self.DiffLines.append(DiffLine)

			self._cal_v0()

			self.version_status = 'new'

		elif 'VERSION' in inp[0]:
			jcpdsfile = open(file, 'r')
			while True:
				jcpdsline = jcpdsfile.readline()
				if jcpdsline == '':
					break

				jlinespl = jcpdsline.split()

				if jlinespl[0] == 'VERSION:':
					version = int(jlinespl[1])
					self.version = version

				if jlinespl[0] == 'COMMENT:':
					header = ' '.join(jlinespl[1:])
					self.comments = header

				if jlinespl[0] == 'K0:':
					k0 = float(jlinespl[1])
					self.k0 = k0

				if jlinespl[0] == 'K0P:':
					k0p = float(jlinespl[1])
					self.k0p = k0p

				if jlinespl[0] == 'DK0DT:':
					dk0dt = float(jlinespl[1])
					self.dk0dt = dk0dt

				if jlinespl[0] == 'DK0PDT:':
					dk0pdt = float(jlinespl[1])
					self.dk0pdt = dk0pdt

				if jlinespl[0] == 'SYMMETRY:':
					self.symmetry = jlinespl[1].lower()

				if jlinespl[0] == 'A:':
					a = float(jlinespl[1])
					self.a0 = a

				if jlinespl[0] == 'B:':
					b = float(jlinespl[1])
					self.b0 = b

				if jlinespl[0] == 'C:':
					c = float(jlinespl[1])
					self.c0 = c

				if jlinespl[0] == 'ALPHA:':
					alpha = float(jlinespl[1])
					self.alpha0 = alpha

				if jlinespl[0] == 'BETA:':
					beta = float(jlinespl[1])
					self.beta0 = beta

				if jlinespl[0] == 'GAMMA:':
					gamma = float(jlinespl[1])
					self.gamma0 = gamma

				if jlinespl[0] == 'VOLUME:':
					v = float(jlinespl[1])
					self.v0 = v

				if jlinespl[0] == 'ALPHAT:':
					alphat = float(jlinespl[1])
					self.thermal_expansion = alphat

				if jlinespl[0] == 'DALPHATDT:':
					dalphatdt = float(jlinespl[1])
					self.thermal_expansion_dt = dalphatat

				if jlinespl[0] == 'DIHKL:':
					DiffLine = DiffractionLine()
					DiffLine.dsp0 = float(jlinespl[1])
					DiffLine.intensity = float(jlinespl[2])
					DiffLine.h = float(jlinespl[3])
					DiffLine.k = float(jlinespl[4])
					DiffLine.l = float(jlinespl[5])
					self.DiffLines.append(DiffLine)
					jcpdslines = jcpdsfile.readlines()
					for line in jcpdslines[0:]:
						item = line.split()
						DiffLine = DiffractionLine()
						DiffLine.dsp0 = float(item[1])
						DiffLine.intensity = float(item[2])
						DiffLine.h = float(item[3])
						DiffLine.k = float(item[4])
						DiffLine.l = float(item[5])
						self.DiffLines.append(DiffLine)

			if self.symmetry == 'cubic':
				self.b0 = self.a0
				self.c0 = self.a0
				self.alpha0 = 90.
				self.beta0 = 90.
				self.gamma0 = 90.
			elif self.symmetry == 'manual':
				self.b0 = self.a0
				self.c0 = self.a0
				self.alpha0 = 90.
				self.beta0 = 90.
				self.gamma0 = 90.
			elif self.symmetry == 'hexagonal' or self.symmetry == 'trigonal':
				self.b0 = self.a0
				self.alpha0 = 90.
				self.beta0 = 90.
				self.gamma0 = 120.
			elif self.symmetry == 'tetragonal':
				self.b0 = self.a0
				self.alpha0 = 90.
				self.beta0 = 90.
				self.gamma0 = 90.
			elif self.symmetry == 'orthorhombic':
				self.alpha0 = 90.
				self.beta0 = 90.
				self.gamma0 = 90.
			elif self.symmetry == 'monoclinic':
				self.alpha0 = 90.
				self.gamma0 = 90.
			#elif self.symmetry == 'triclinic':

			jcpdsfile.close()

			if self.v0 == 0.:
				self._cal_v0()

			self.version_status = 'new'

		else:
			self.version_status = 'old'

		if self.version_status == 'new':
			self.a = self.a0
			self.b = self.b0
			self.c = self.c0
			self.alpha = self.alpha0
			self.beta = self.beta0
			self.gamma = self.gamma0
			self.v = self.v0

class Paradialogj(QDialog):
	'''define the dialog for setting parameters'''
	def __init__(self, parent=None):
		QDialog.__init__(self, parent)
		self.resize(240, 200)

		grid = QGridLayout()

		grid.addWidget(QLabel('wavelength(nm)'), 0, 0, 1, 1)
		self.wwavelengthj = QLineEdit()
		grid.addWidget(self.wwavelengthj, 0, 1, 1, 1)

		grid.addWidget(QLabel('pressure(Gpa)'), 1, 0, 1, 1)
		self.wpressurej = QLineEdit()
		grid.addWidget(self.wpressurej, 1, 1, 1, 1)

		btbox = QDialogButtonBox()
		btbox.setOrientation(Qt.Horizontal)
		btbox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Save)
		btbox.accepted.connect(self.accept)
		btbox.rejected.connect(self.reject)

		layout = QVBoxLayout()
		layout.addLayout(grid)
		layout.addWidget(btbox)

		self.setWindowIcon(QIcon('./icons/set.png'))
		self.setWindowTitle('set parameters')
		self.setLayout(layout)

	def wl(self):
		return self.wwavelengthj.text()

	def pre(self):
		return self.wpressurej.text()

class Paradialogc(QDialog):
	'''define the dialog for setting parameters'''
	def __init__(self, parent=None):
		QDialog.__init__(self, parent)
		self.resize(240, 240)

		grid = QGridLayout()

		grid.addWidget(QLabel('2θ min value(degree)'), 0, 0, 1, 1)
		self.wtwothetamin = QLineEdit()
		grid.addWidget(self.wtwothetamin, 0, 1, 1, 1)

		grid.addWidget(QLabel('2θ max value(degree)'), 1, 0, 1, 1)
		self.wtwothetamax = QLineEdit()
		grid.addWidget(self.wtwothetamax, 1, 1, 1, 1)

		grid.addWidget(QLabel('k0(Gpa)'), 2, 0, 1, 1)
		self.wk0 = QLineEdit()
		grid.addWidget(self.wk0, 2, 1, 1, 1)

		grid.addWidget(QLabel('k0p'), 3, 0, 1, 1)
		self.wk0p = QLineEdit()
		grid.addWidget(self.wk0p, 3, 1, 1, 1)

		grid.addWidget(QLabel('alphat'), 4, 0, 1, 1)
		self.walphat = QLineEdit()
		grid.addWidget(self.walphat, 4, 1, 1, 1)

		grid.addWidget(QLabel('wavelength(nm)'), 5, 0, 1, 1)
		self.wwavelengthc = QLineEdit()
		grid.addWidget(self.wwavelengthc, 5, 1, 1, 1)

		grid.addWidget(QLabel('comments'), 6, 0, 1, 1)
		self.wcomments = QLineEdit()
		grid.addWidget(self.wcomments, 6, 1, 1, 1)

		btbox = QDialogButtonBox()
		btbox.setOrientation(Qt.Horizontal)
		btbox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Save)
		btbox.accepted.connect(self.accept)
		btbox.rejected.connect(self.reject)

		layout = QVBoxLayout()
		layout.addLayout(grid)
		layout.addWidget(btbox)

		self.setWindowIcon(QIcon('./icons/set.png'))
		self.setWindowTitle('set parameters')
		self.setLayout(layout)

	def ttmin(self):
		return self.wtwothetamin.text()

	def ttmax(self):
		return self.wtwothetamax.text()

	def k0(self):
		return self.wk0.text()

	def k0p(self):
		return self.wk0p.text()

	def alphat(self):
		return self.walphat.text()

	def wl(self):
		return self.wwavelengthc.text()

	def comment(self):
		return self.wcomments.text()

if __name__ == '__main__':

	app = QApplication(sys.argv)
	main = Main()
	sys.exit(app.exec_())
