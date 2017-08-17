import sys

import pymatgen as mg
from pymatgen import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from ds_jcpds import jcpds_modified

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

	#def loadandreadj(self):


	def loadandreadj(self):
		#choose which jcpds file to import and display information
		self.tabWidget.setCurrentWidget(self.infoj)
		fname = QFileDialog.getOpenFileName(self, 'Import jcpds', '/home', '*.jcpds')
		self.jcpds_filename = fname[0]
		self.jedit.setText('From: ' + self.jcpds_filename + '\n')
		self.jedit.append(5*'-' + 'crystal information' + 5*'-')
		crystal = jcpds_modified.JCPDS(filename = self.jcpds_filename)
		crystal.read_file(self.jcpds_filename)
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

		tfile = jcpds_modified.JCPDS()
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
		pattern = jcpds_modified.JCPDS(filename = self.jcpds_filename)

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
