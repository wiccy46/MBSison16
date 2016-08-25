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
        self.statusBar().showMessage('Ready')  # Tell user to wait while sending data
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
        genDataButton.setToolTip('This button randomise a \
            set of data based on the input dimension and cluster amount.')
        genDataButton.resize(genDataButton.sizeHint())
        # genDataButton.move(50, 550)
        genDataButton.clicked.connect(self.genData)
        sendDataButton = QtGui.QPushButton("Send Data", self)
        sendDataButton.resize(genDataButton.sizeHint())
        # sendDataButton.move(200, 550)
        sendDataButton.clicked.connect(self.sendData)
        #----------------------

        #TODO need a display box.
        sigmaTitle = QtGui.QLabel('Sigma')
        sigmaSlider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        sigmaSlider.setRange(0, 100)
        sigmaSlider.setValue(50)
        sigmaSlider.valueChanged[int].connect(self.sigmaChangeValue)
        sigmaDisplay = QtGui.QLineEdit()
        sigmaDisplay.setFixedSize(50, 25)


        dt_init = 0.005
        dtTitle = QtGui.QLabel('dt')
        dtSlider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        dtSlider.setRange(0, 100)
        temp = linlinInvert(dt_init, 0, 100, 0.001, 0.01)
        # print "dtInit " +  str(linlin(temp,0, 100, 0.001, 0.01))
        dtSlider.setValue(temp)
        dtSlider.valueChanged[int].connect(self.dtChangeValue)
        dtDisplay = QtGui.QLineEdit()
        dtDisplay.setFixedSize(50, 25)
        dtDisplay.setText(str(dt_init))

        rTitle = QtGui.QLabel('r')
        rSlider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        rSlider.setRange(0, 100)
        r_init = 0.999
        temp = linlinInvert(r_init, 0, 100, 0.99, 1.0)
        # print "rInit " + str(linlin(temp,0, 100, 0.99, 1.0))
        dtSlider.setValue(temp)
        rSlider.valueChanged[int].connect(self.rChangeValue)
        rDisplay = QtGui.QLineEdit()
        rDisplay.setFixedSize(50, 25)
        rDisplay.setText(str(r_init))


        cltLeftBox = QtGui.QGridLayout()
        cltLeftBox.addWidget(genDataButton, 1, 0 )
        cltLeftBox.addWidget(sendDataButton, 1, 1)


        cltRightBox = QtGui.QGridLayout()
        cltRightBox.addWidget(sigmaTitle, 1, 0)
        cltRightBox.addWidget(sigmaSlider, 1, 1)
        cltRightBox.addWidget(sigmaDisplay, 1, 2)
        cltRightBox.addWidget(dtTitle, 2, 0)
        cltRightBox.addWidget(dtSlider, 2, 1)
        cltRightBox.addWidget(dtDisplay, 2, 2)
        cltRightBox.addWidget(rTitle, 3, 0)
        cltRightBox.addWidget(rSlider, 3, 1)
        cltRightBox.addWidget(rDisplay, 3, 2)


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

    def sigmaChangeValue(self, value):
        smi, sma = 0, 100
        dmi, dma = 0.001, self.max_sigma
        self.sigma = linlin(value, smi, sma, dmi, dma)
        print self.sigma

    def dtChangeValue(self, value):
        smi, sma = 0, 100
        dmi, dma = 0.001, 0.01
        self.dt = linlin(value, smi, sma, dmi, dma)
        print self.dt

    def rChangeValue(self, value):
        smi, sma = 0, 100
        dmi, dma = 0.99, 1.0
        self.r = linlin(value, smi, sma, dmi, dma)
        print self.r

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