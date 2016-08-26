import Trajectory
import numpy as np
from io import StringIO
import time, OSC, socket
import threading, sys
from PyQt4 import QtGui, QtCore
import pyo64 as pyo
import matplotlib.pyplot as plt
from threading import Thread

from lib.Exptable import Exptable
from lib.DSP import DSP
from lib.DataGen import DataGen
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

"""
Set up Audio server
"""


FS = 44100/4
BLOCK = 1024
HOMEIP = '192.168.178.47'
OFFICEIP = "129.70.149.23"


# TODO put this in the GUI so that there is a button to stop.
s = pyo.Server(sr=FS, nchnls=2, buffersize=BLOCK, duplex=0).boot()
s.start()
# Need to adjust the maxsize to test the speed.

def linlin(x, smi, sma, dmi, dma): return (x-smi)/float(sma-smi)*(dma-dmi)+dmi
def linlinInvert(x, smi, sma, dmi, dma): return (x-dmi)*(sma-smi)/(dma-dmi) + smi


class Ptsgui(QtGui.QMainWindow):

    def __init__(self):
        super(Ptsgui, self).__init__()
        self.fifo = pyo.FIFOPlayer(maxsize=20, mul=[.5, .5])
        self.mixer = pyo.Mixer()
        self.mixer.addInput(0, self.fifo);
        self.mixer.setAmp(0, 0, 1)
        self.mixer.addInput(1, self.fifo);
        self.mixer.setAmp(1, 1, 1)
        self.mixer.out()
        self.fifo.out()

        self.data, self.pos = np.zeros(2), np.zeros(2)
        self.vel , self.t = 0, 1.
        self.exp_table = np.zeros(1) # Initialise the exp table
        self.norm_max = 0
        self.max_sigma = 0.4
        self.exp_resolution = 1000000  # Resolution for lookup table

        # TODO These parameter should become sliders
        self.m_comp, self.dt, self.resistant, self.sigma = 0.02, 0.005, 0.999, 0.15  # init
        self.t = 1.0  # Time in second per piece, TODO a para
        self.blockSize = 5000  # Buffer size for trajectory TODO a para
        self.audioVecSize = int(self.t * self.blockSize)  # I define that 5000 steps will return 1 second of audio
        self.windowing = DSP().makeWindow(self.audioVecSize, rampUp=0.05, rampDown=0.05)  # Windowing for audio
        self.initUI()

    def proc(self, stopevent, soundOutput):
        self.fifo.put(soundOutput)

    def datafigon_pick(self, event):
        ind = np.array(event.ind)  # Get the index of the clicked data
        self.pos = np.array(self.data[ind[0], :])
        self.vel = np.random.rand(self.dim) - 0.5

        trj, junk, forceSound = Trajectory.PTSM(self.pos, self.data, self.vel, self.exp_table, self.exp_resolution, \
                                                self.norm_max, sigma=self.sigma, dt=self.dt, r=self.resistant, \
                                                Nsamp=self.audioVecSize, compensation=self.m_comp)
        velSound = trj[:, 0] / np.max(np.absolute(trj[:, 0])) * self.windowing
        velSound = DSP().butter_lowpass_filter(velSound, 2000.0, FS, 6)  # 6th order
        # forceSound = forceSound / np.max(np.absolute(forceSound))
        stopevent = threading.Event()
        producer = threading.Thread(name="Compute audio signal", target=self.proc, args=[ stopevent, velSound])
        producer.start()
        self.drawTrj(trj)


    def drawTrj(self, trj): # Plot trajectory
        self.trjFigAx.clear()
        potmap = np.zeros((40, 40))
        for i in range(40):
            for j in range(40):
                x = float(i) / 40 - 0.5
                y = float(j) / 40 - 0.5
                potmap[j, i] = Trajectory.potential_ds(self.data[:, 0:2], np.array([x, y]), self.sigma)
        self.trjFigAx.matshow(potmap, cmap=plt.cm.gray ,  extent=(-0.5, 0.5, 0.5, -0.5))
        self.trjFigAx.plot(self.data[:, 0], self.data[:, 1], ".")
        self.trjFigAx.plot(trj[:, 1], trj[:, 2], "-", lw=0.7, c="green")
        self.trjFigAx.plot(trj[0, 1], trj[0, 2], "o", c="yellow")
        self.trjFigAx.plot(trj[self.audioVecSize - 1, 1], trj[self.audioVecSize - 1, 2], "x", c="red")
        self.trjFigAx.set_xlim([-0.5, 0.5])
        self.trjFigAx.set_ylim([-0.5, 0.5])
        self.trjFigCanvas.draw()

    def sendData(self):
        print "Send data"

    def genData(self):
        self.dataFigAx.clear()
        self.data = DataGen().datagen(4, 3, sigma=0.4, minnr=100, maxnr=300)
        self.N, self.dim = self.data.shape[0], self.data.shape[1]
        self.norm_max = self.dim
        self.max_sigma = np.log(self.dim) / 2.5 * 0.5 + 0.1
        delta_dis = np.linspace(0.0, self.dim, num=self.exp_resolution)  # Range of distance difference.
        self.exp_table = Exptable().updateExpTable(self.sigma, delta_dis)
        self.dataFigAx.scatter(self.data[:, 0], self.data[:, 1], c="blue", picker=2, s=[50] *self.N)
        self.dataFigCanvas.draw()

    def initUI(self):

        self.statusBar().showMessage('Ready, Move on each item to see user tip.')  # Tell user to wait while sending data
        self.setGeometry(50, 50, 1050, 650)
        self.setWindowTitle('PTS')

        mainlayout = QtGui.QVBoxLayout()
        self.widget = QtGui.QWidget()
        self.widget.setLayout(mainlayout)
        self.setCentralWidget(self.widget)  # Set main layout to layout of widget!!!



        # plotBox contents (sublay 1)
        self.dataFig = plt.figure()
        self.dataFigCanvas = FigureCanvas(self.dataFig)
        self.dataFigCanvas.mpl_connect('pick_event', self.datafigon_pick)
        self.dataFigAx = self.dataFig.add_subplot(111)
        self.dataFigAx.set_xlim([-0.5, 0.5])
        self.dataFigAx.set_ylim([-0.5, 0.5])
        self.dataFigToolbar = NavigationToolbar(self.dataFigCanvas,self)
        plotBoxL = QtGui.QVBoxLayout()
        plotBoxL.setSpacing(5)
        plotBoxL.addWidget(self.dataFigCanvas)
        plotBoxL.addWidget(self.dataFigToolbar)

        self.trjFig = plt.figure()
        self.trjFigCanvas = FigureCanvas(self.trjFig)
        self.trjFigAx = self.trjFig.add_subplot(111)
        self.trjFigAx.set_xlim([-0.5, 0.5])
        self.trjFigAx.set_ylim([-0.5, 0.5])
        self.trjFigToolbar = NavigationToolbar(self.trjFigCanvas,self)
        plotBoxR = QtGui.QVBoxLayout()
        plotBoxR.setSpacing(5)
        plotBoxR.addWidget(self.trjFigCanvas)
        plotBoxR.addWidget(self.trjFigToolbar)

        #-----------------
        # cltBox contents (sublay 1)
        genDataButton = QtGui.QPushButton('Generate Data', self)
        genDataButton.setToolTip('This button randomise a '
                                 'set of data based on the input dimension and cluster amount.')
        genDataButton.resize(genDataButton.sizeHint())
        genDataButton.clicked.connect(self.genData)
        sendDataButton = QtGui.QPushButton("Send Data", self)
        sendDataButton.resize(genDataButton.sizeHint())
        sendDataButton.clicked.connect(self.sendData)
        sendDataButton.setToolTip('Send the first 2 columns of data to Android device for plotting.')
        #----------------------

        #sigma
        sigmaTitle = QtGui.QLabel('Sigma')
        sigmaTitle.resize(sigmaTitle.sizeHint())
        self.sigmaSlider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.sigmaSlider.setToolTip('Large sigma results in 1 global minimum while small sigma results in local minimum '
                               'in each data point')

        self.sigmaSlider.setRange(0, 100)
        self.sigmaSlider.setValue(50)
        self.sigmaSlider.valueChanged[int].connect(self.changeValue)
        self.sigmaDisplay = QtGui.QLineEdit()
        self.sigmaDisplay.setFixedSize(50, 20)

        #dt
        dt_init = 0.005
        dtTitle = QtGui.QLabel('dt')
        self.dtSlider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.dtSlider.setToolTip('A larger dt will lead to a bigger distance change at each iteration.')
        self.dtSlider.setRange(0, 100)
        temp = linlinInvert(dt_init, 0, 100, 0.001, 0.01)
        self.dtSlider.setValue(temp)
        self.dtDisplay = QtGui.QLineEdit()
        self.dtDisplay.setFixedSize(50, 20)
        self.dtDisplay.setText(str(dt_init))
        self.dtSlider.valueChanged[int].connect(self.changeValue)

        #r
        rTitle = QtGui.QLabel('r')
        self.rSlider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.rSlider.setRange(0, 100)
        r_init = 0.999
        temp = linlinInvert(r_init, 0, 100, 0.99, 1.0)
        # print "rInit " + str(linlin(temp,0, 100, 0.99, 1.0))
        self.rSlider.setValue(temp)
        self.rSlider.setToolTip('r is the recipricate of friction. at r = 1, the particle will have no velocity lost.')
        self.rDisplay = QtGui.QLineEdit()
        self.rDisplay.setFixedSize(50, 20)
        self.rDisplay.setText(str(r_init))
        self.rSlider.valueChanged[int].connect(self.changeValue)
        #----------------------#




        # left pannel.
        dimTitle = QtGui.QLabel('Dimension')
        dimTitle.setFixedSize(70, 20)
        self.dimDisplay = QtGui.QSpinBox()

        self.dimDisplay.setRange(2, 1000)
        self.dimDisplay.setValue(3)
        self.dimDisplay.setFixedSize(50, 20)
        self.dimDisplay.resize(self.dimDisplay.sizeHint())
        self.dimDisplay.setToolTip('Choose the dimensionality of the data. Type = int, Min = 2, Max = 1000')
        self.dimDisplay.valueChanged[int].connect(self.changeValue)
        # self.rSlider.valueChanged[int].connect(self.changeValue)


        ncTitle = QtGui.QLabel('Clusters')

        ncTitle.setFixedSize(70, 20)
        self.ncDisplay = QtGui.QSpinBox()
        self.ncDisplay.setValue(4)
        self.ncDisplay.setRange(1, 1000)
        self.ncDisplay.setFixedSize(50, 20)
        self.ncDisplay.resize(self.ncDisplay.sizeHint())
        self.ncDisplay.setToolTip('Choose the number of clusters for the simulated data. Type = int, Min = 1, Max = 1000')
        self.ncDisplay.valueChanged[int].connect(self.changeValue)

        nrminTitle = QtGui.QLabel('Nr_min')
        nrminTitle.setFixedSize(70, 20)
        self.nrminDisplay = QtGui.QSpinBox()
        self.nrminDisplay.setRange(10, 300)
        self.nrminDisplay.setValue(50)
        self.nrminDisplay.setFixedSize(50, 20)
        self.nrminDisplay.resize(self.nrminDisplay.sizeHint())
        self.nrminDisplay.setToolTip('Choose the minimum number of row for each clusters. 10~300, type: int')
        # Use the changed value box to detect is nr_max < mr_min
        self.nrminDisplay.valueChanged[int].connect(self.changeValue)

        nrmaxTitle = QtGui.QLabel('Nr_max')
        nrmaxTitle.setFixedSize(50, 20)
        self.nrmaxDisplay = QtGui.QSpinBox()
        self.nrmaxDisplay.setRange(100, 600)
        self.nrmaxDisplay.setValue(200)
        self.nrmaxDisplay.setFixedSize(50, 20)
        self.nrmaxDisplay.resize(self.nrmaxDisplay.sizeHint())
        self.nrmaxDisplay.setToolTip('Choose the maximum number of row for each clusters. 100~600, type: int')
        # Use the changed value box to detect is nr_max < mr_min
        self.nrmaxDisplay.valueChanged[int].connect(self.changeValue)


        ipTitle = QtGui.QLabel('Tablet IP')
        ipTitle.setFixedSize(80, 20)
        # ipTitle.resize(dimTitle.sizeHint())
        self.ipDisplay = QtGui.QLineEdit()
        self.ipDisplay.setText(HOMEIP)
        # self.ipDisplay.resize(self.dimDisplay.sizeHint())
        self.ipDisplay.setFixedSize(106, 20)
        self.ipDisplay.returnPressed.connect(self.ipChangeValue)
        # self.ipDisplay.valueChanged[str].connect(self.changeValue)


        cltLeftBox = QtGui.QGridLayout()
        cltLeftBox.addWidget(ncTitle,1 , 0), cltLeftBox.addWidget(self.ncDisplay, 1, 1)
        cltLeftBox.addWidget(dimTitle, 1, 2), cltLeftBox.addWidget(self.dimDisplay, 1, 3)
        cltLeftBox.addWidget(nrminTitle, 2, 0), cltLeftBox.addWidget(self.nrminDisplay, 2, 1)
        cltLeftBox.addWidget(nrmaxTitle, 2, 2), cltLeftBox.addWidget(self.nrmaxDisplay, 2, 3)
        cltLeftBox.addWidget(genDataButton, 3, 0 )
        cltLeftBox.addWidget(ipTitle, 4,0), cltLeftBox.addWidget(self.ipDisplay, 4, 1)
        cltLeftBox.addWidget(sendDataButton, 6, 0, 6, 4)



        # Add right box
        cltRightBox = QtGui.QGridLayout()
        cltRightBox.addWidget(sigmaTitle, 1, 0)
        cltRightBox.addWidget(self.sigmaSlider, 1, 1)
        cltRightBox.addWidget(self.sigmaDisplay, 1, 2)
        cltRightBox.addWidget(dtTitle, 2, 0)
        cltRightBox.addWidget(self.dtSlider, 2, 1)
        cltRightBox.addWidget(self.dtDisplay, 2, 2)
        cltRightBox.addWidget(rTitle, 3, 0)
        cltRightBox.addWidget(self.rSlider, 3, 1)
        cltRightBox.addWidget(self.rDisplay, 3, 2)


        # Sub layout 1
        plotBox = QtGui.QHBoxLayout()
        plotBox.addLayout(plotBoxL)
        plotBox.addLayout(plotBoxR)

        cltBox = QtGui.QHBoxLayout()

        cltBox.addLayout(cltLeftBox)
        cltBox.addLayout(cltRightBox)

        mainlayout.addLayout(plotBox)
        mainlayout.addLayout(cltBox)
        self.show()

    def ipChangeValue(self):
        print self.ipDisplay.text()


    def changeValue(self,value):
        if self.sender() == self.sigmaSlider:
            smi, sma = 0, 100
            dmi, dma = 0.001, self.max_sigma
            self.sigma = linlin(value, smi, sma, dmi, dma)
            self.sigmaDisplay.setText(str(self.sigma))
        elif self.sender() == self.dtSlider:
            smi, sma = 0, 100
            dmi, dma = 0.001, 0.01
            self.dt = linlin(value, smi, sma, dmi, dma)
            self.dtDisplay.setText(str(self.dt))
        elif self.sender() == self.rSlider:
            smi, sma = 0, 100
            dmi, dma = 0.99, 1.0
            self.r = linlin(value, smi, sma, dmi, dma)
            self.rDisplay.setText(str(self.r))


        elif self.sender() == self.dimDisplay:
            print "dimension = " + str(value)
        elif self.sender() == self.ncDisplay:
            print "nc = " + str(value)
        elif self.sender() == self.nrminDisplay:
            print "nr min = " + str(value)
        elif self.sender() == self.nrmaxDisplay:
            print "nr max = " + str(value)
        else: pass



    def closeEvent(self, event):
		reply = QtGui.QMessageBox.question(self, 'Message',
            "Are you sure to quit?", QtGui.QMessageBox.Yes |
            QtGui.QMessageBox.No, QtGui.QMessageBox.Yes)

		if reply == QtGui.QMessageBox.Yes:
			event.accept()
		else:
			event.ignore()


def main():
	app = QtGui.QApplication(sys.argv)
	ptsgui = Ptsgui()
	sys.exit(app.exec_())

if __name__ == '__main__':
	main()