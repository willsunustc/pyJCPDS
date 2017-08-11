import sys

import pymatgen as mg
from pymatgen import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

import ds_jcpds

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

class Plot_vlines(QMainWindow):
	"""plot a vlines figure from jcpds and optionally cif file."""
	def __init__(self):
		super().__init__()
		self.pressure = 0.
		self.wavelength = 0.
		self.twothetamin = 0
		self.twothetamin = 0
		self.diff_lines = 0.
		self.jcpds_filename = ''
		self.cif_filename = ''
		self.initUI()

	def initUI(self):
		self.resize(800, 600)
		self.center()

		self.figure = plt.figure()
		self.drawing = self.figure.add_subplot(111)
		self.canvas = FigureCanvas(self.figure)

		self.setCentralWidget (self.canvas)

		impAct = QAction(QIcon('import.png'), 'Import jcpds', self)
		impAct.setShortcut('Ctrl+O')
		impAct.triggered.connect(self.loadandread)

		setAct = QAction(QIcon('set.png'), 'Set parameters', self)
		setAct.setShortcut('Ctrl+L')
		setAct.triggered.connect(self.setparameters)

		plotAct = QAction(QIcon('plot.png'), 'Plot', self)
		plotAct.setShortcut('Ctrl+P')
		plotAct.triggered.connect(self.plot_vlines_jcpds)

		exitAct = QAction(QIcon('exit.png'), 'Exit', self)
		exitAct.setShortcut('Ctrl+Q')
		exitAct.triggered.connect(qApp.quit)

		tbimp = self.addToolBar('Import jcpds')
		tbimp.addAction(impAct)

		tbset = self.addToolBar('Set parameters')
		tbset.addAction(setAct)

		tbplot = self.addToolBar('Plot')
		tbplot.addAction(plotAct)

		tbexit = self.addToolBar('Exit')
		tbexit.addAction(exitAct)

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

	def loadandread(self):
		#choose which jcpds file to import
		fname = QFileDialog.getOpenFileName(self, 'Import jcpds', '/home')
		self.jcpds_filename = fname[0]

		'''just for test
		if fname[0]:
			f = open(fname[0], 'r')

			with f:
				jcpdscontent = f.read()
				self.jcpdsEdit.setText(jcpdscontent)
		'''

	def setparameters(self):
		setdialog = Paradialog()
		if setdialog.exec_():
			self.twothetamin = int(setdialog.ttmin())
			self.twothetamax = int(setdialog.ttmax())
			self.pressure = float(setdialog.pre())
			self.wavelength = float(setdialog.wl())

		setdialog.destroy()

		'''just for test
		self.jcpdsEdit.setText(str(self.twothetamin))
		self.jcpdsEdit.append(str(self.twothetamax))
		self.jcpdsEdit.append(str(self.pressure))
		self.jcpdsEdit.append(str(self.wavelength))
		'''

	def plot_vlines_jcpds(self):
		#import jcpds file
		pattern = ds_jcpds.JCPDS(filename = self.jcpds_filename)

		#calculate the pattern
		pattern.cal_dsp(pressure = self.pressure)
		pattern.get_DiffractionLines()
		tth, inten = pattern.get_tthVSint(wavelength = self.wavelength)

		#plot
		self.drawing.hold (False)
		self.drawing.plot(tth, inten, color = 'r', label = 'jcpds')
		self.canvas.draw()

class Paradialog(QDialog):
	'''define the dialog for setting parameters'''
	def __init__(self, parent=None):
		QDialog.__init__(self, parent)
		self.resize(240, 200)

		grid = QGridLayout()

		grid.addWidget(QLabel('2\Theta min value(degree)'), 0, 0, 1, 1)
		self.twothetamin = QLineEdit()
		grid.addWidget(self.twothetamin, 0, 1, 1, 1)

		grid.addWidget(QLabel('2\Theta max value(degree)'), 1, 0, 1, 1)
		self.twothetamax = QLineEdit()
		grid.addWidget(self.twothetamax, 1, 1, 1, 1)

		grid.addWidget(QLabel('pressure(Gpa)'), 2, 0, 1, 1)
		self.pressure = QLineEdit()
		grid.addWidget(self.pressure, 2, 1, 1, 1)

		grid.addWidget(QLabel('wavelength(nm)'), 3, 0, 1, 1)
		self.wavelength = QLineEdit()
		grid.addWidget(self.wavelength, 3, 1, 1, 1)

		btbox = QDialogButtonBox()
		btbox.setOrientation(Qt.Horizontal)
		btbox.setStandardButtons(QDialogButtonBox.Cancel|QDialogButtonBox.Save)
		btbox.accepted.connect(self.accept)
		btbox.rejected.connect(self.reject)

		layout = QVBoxLayout()
		layout.addLayout(grid)
		layout.addWidget(btbox)

		self.setWindowIcon(QIcon('set.png'))
		self.setWindowTitle('set parameters')
		self.setLayout(layout)

	def ttmin(self):
		return self.twothetamin.text()

	def ttmax(self):
		return self.twothetamax.text()

	def pre(self):
		return self.pressure.text()

	def wl(self):
		return self.wavelength.text()

if __name__ == '__main__':

	app = QApplication(sys.argv)
	plot = Plot_vlines()
	sys.exit(app.exec_())
