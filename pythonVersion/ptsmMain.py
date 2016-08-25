import Trajectory
import numpy as np
from io import StringIO
import time, OSC, socket
import threading, sys
from PyQt4 import QtGui
import pyo64 as pyo
import matplotlib.pyplot as plt
from threading import Thread

from lib.Exptable import Exptable
from lib.DSP import DSP
from lib.DataGen import DataGen
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

def linlin(x, smi, sma, dmi, dma): return (x-smi)/(sma-smi)*(dma-dmi)+dmi
# Turn list into integer, not used.
def list2int(numbers):
    return int(''.join(["%d" % x for x in numbers]))

"""
Set up Audio server
"""
# Is 1024 too large for the bufferSize in this case.
fs = 44100/5
bs = 1024
s = pyo.Server(sr=fs, nchnls=2, buffersize=bs, duplex=0).boot()
s.start()
# Need to adjust the maxsize to test the speed.
fifo = pyo.FIFOPlayer(maxsize=20, mul=[.5,.5])
mixer = pyo.Mixer()
mixer.addInput(0, fifo);mixer.setAmp(0,0,1)
mixer.addInput(1, fifo);mixer.setAmp(1,1,1)
mixer.out()
fifo.out()



class soundOption(object):
    def velocity(self, event):
        global which2play
        which2play = "velocity"

    def force(self, event):
        global which2play
        which2play = "force"

def proc(fifo, stopevent, soundOutput):
    fifo.put(soundOutput)



# Well I am going to have new sliders anyway
# TODO slider update in a class. or not...
# Slider update function, for the interaction
def sliderUpdate(val):
    global dt, resistant
    sigma = ssigma.val
    dt = sdt.val
    resistant = sresistant.val

def sigmaUpdate(val):
    global sigma, exp_table, delta_dis
    sigma = ssigma.val
    exp_table = Exptable().updateExpTable(sigma, delta_dis)


# TODO all plotting in a class
def plotPotential(data2d, fig2, sigma=.1, Nx=40, Ny=40):
    # ensure the data sending for plotPotential is always in 2d.
    potmap = np.zeros((Nx, Ny))
    for i in range(Nx):
        for j in range(Ny):
            x = float(i) / Nx - 0.5
            y = float(j) / Ny - 0.5
            potmap[j, i] = Trajectory.potential_ds(data2d, np.array([x, y]), sigma)
    fig2.clf()
    plt.matshow(potmap, cmap=plt.cm.gray, extent=(-0.5, 0.5, 0.5, -0.5), fignum=2)
    fig2.gca().plot(data[:, 0], data[:, 1], ".")
    plt.show()

# Currently not used.
def plotTrajectory(data, sigma, trj, resultWindow, audioVecSize):
    plotPotential(data[:, 0:2], fig2=resultWindow, sigma=sigma)
    resultWindow.gca().plot(trj[:, 1], trj[:, 2], "-", lw=0.7, c="green")
    # Mark the beginning and the end of trajectory .
    resultWindow.gca().plot(trj[0, 1], trj[0, 2], "o", c="yellow")
    resultWindow.gca().plot(trj[audioVecSize - 1, 1], trj[audioVecSize - 1, 2], "x", c="red")
    resultWindow.gca().axis([-0.6, 0.6, -0.6, 0.6])
    resultWindow.canvas.draw()



def initializePot(data, N):  # It takes data and number of rows.
    sctPlot = ax.scatter(data[:, 0], data[:, 1], c="blue", picker=2, s=[50] * N)
    fig.subplots_adjust(bottom=0.25, left=0.1)
    plt.grid(False)
    plt.axis([-0.6, 0.6, -0.6, 0.6])
    return sctPlot

def spectrum(av, fs):
    NFFT = 1024
    plt.figure()
    plt.specgram(av, NFFT=NFFT, Fs=fs, noverlap=900, cmap=plt.cm.gist_heat)
    plt.show()

# TODO: a button to regenerate new data set. And send through
data = DataGen().datagen(4, 3,sigma = 0.4, minnr = 100, maxnr = 300)
N, dim = data.shape[0], data.shape[1]
# Create exp lookup

exp_resolution = 1000000 # Resolution for lookup table
delta_dis = np.linspace(0.0, dim, num = exp_resolution) # Range of distance difference.
# m_comp = 0.001 # mass compensation
norm_max = dim
sigma = 0.1

exp_table = Exptable().updateExpTable(sigma, delta_dis)

fig, ax = plt.subplots()
m_comp, dt, resistant, sigma = 0.02, 0.005 , 0.999, 0.15 # init
max_sigma = np.log(dim)/2.5 * 0.5 + 0.1
t = 1.0 # Time in second per piece
blockSize = 5000   # Buffer size for trajectory
audioVecSize  = int(t *  blockSize)  # I define that 5000 steps will return 1 second of audio
exp_table, norm_max = Exptable().createExpTable(dim, sigma, exp_resolution = exp_resolution)
windowing = DSP().makeWindow(audioVecSize, rampUp = 0.05, rampDown = 0.05) # Windowing for audio
velSound = np.zeros(audioVecSize)
forceSound = np.zeros(dim*audioVecSize) # This is only meant for visualisation.
drawResult = True
which2play = "velocity" # Initialise the button press
pos2d = np.zeros(2)
resultWindow = plt.figure(2, figsize=(8, 8))
# plotPotential(data[:, 0:2], fig2 = resultWindow,  sigma = sigma )

sctPlot = initializePot(data, N)


def on_pick(event):
    # In the continuous mode, window shouldn't be used.
    global data, sigma, resistant, dt, resultWindow, windowing, audioVecSize, \
        norm_max, velSound, drawResult, dim, forceSound, fs
    vel = (np.random.rand(dim) - 0.5)


    ind = np.array(event.ind)  # Get the index of the clicked data
    pos = np.array(data[ind[0], :])
    # Get PTSM trajectory information
    trj, junk, forceSound = Trajectory.PTSM(pos, data, vel, exp_table, exp_resolution, \
                                 norm_max, sigma=sigma, dt=dt, r=resistant, \
                                 Nsamp=audioVecSize, compensation=m_comp)
    # -------------------------------------------#
    # Process sound #----------
    velSound = trj[:, 0] / np.max(np.absolute(trj[:, 0])) * windowing
    velSound = DSP().butter_lowpass_filter(velSound, 2000.0, fs, 6)  # 6th order
    forceSound = forceSound / np.max(np.absolute(forceSound))
    stopevent = threading.Event()
    producer = threading.Thread(name="Compute audio signal", target=proc, args=[fifo, stopevent, velSound])
    producer.start()
    if drawResult == True:
        #         plotTrajectory(data = data, sigma = sigma, trj = trj,
        #                        resultWindow = resultWindow, audioVecSize = resultWindow)
        plotPotential(data[:, 0:2], fig2=resultWindow, sigma=sigma)
        resultWindow.gca().plot(trj[:, 1], trj[:, 2], "-", lw=0.7, c="green")
        # Mark the beginning and the end of trajectory .
        resultWindow.gca().plot(trj[0, 1], trj[0, 2], "o", c="yellow")
        resultWindow.gca().plot(trj[audioVecSize - 1, 1], trj[audioVecSize - 1, 2], "x", c="red")
        resultWindow.gca().axis([-0.6, 0.6, -0.6, 0.6])
        resultWindow.canvas.draw()


# ---------------#
# Create a slider for setting up the velocity
axcolor = 'lightgoldenrodyellow'
# Create sliders for sigma and dt.
axSigma = plt.axes([0.1, 0.1, 0.8, 0.02], axisbg=axcolor)
axDt = plt.axes([0.1, 0.06, 0.8, 0.02], axisbg=axcolor)
axR = plt.axes([0.1, 0.02, 0.8, 0.02], axisbg=axcolor)
axB1 = plt.axes([0.1, 0.13, 0.1, 0.075], axisbg=axcolor)
axB2 = plt.axes([0.21, 0.13, 0.1, 0.075], axisbg=axcolor)

ssigma = plt.Slider(axSigma, "Sigma", 0.001, max_sigma, valinit=sigma, color='blue')
sdt = plt.Slider(axDt, "dt", 0.001, 0.01, valinit=dt, color='blue')
sresistant = plt.Slider(axR, "r", 0.99, 1.0, valinit=resistant, color='blue')

choice = soundOption()

b1 = plt.Button(axB1, "Velocity")
b1.on_clicked(choice.velocity)
b2 = plt.Button(axB2, "Force")
b2.on_clicked(choice.force)

ssigma.on_changed(sigmaUpdate)
sdt.on_changed(sliderUpdate)
sresistant.on_changed(sliderUpdate)
resultWindow = plt.figure(2, figsize=(8, 8))
fig.canvas.mpl_connect('pick_event', on_pick)
# plt.show()


class Ptsgui(QtGui.QMainWindow):



    def __init__(self):
        super(Ptsgui, self).__init__()

        self.initUI()
        self.data = np.zeros(2)
        self.N = 0
        self.dim = 0


    def datafigon_pick(event):
        ind = np.array(event.ind)  # Get the index of the clicked data
        print ind


    def sendData(self):
        print "Send data"

    def genData(self):
        self.data = DataGen().datagen(4, 3, sigma=0.4, minnr=100, maxnr=300)
        self.N, self.dim = self.data.shape[0], self.data.shape[1]
        self.dataFigAx.scatter(data[:, 0], data[:, 1], c="blue", picker=2, s=[50] * N)
        self.dataFigCanvas.draw()


        print self.data

    def initUI(self):
        self.statusBar().showMessage('Ready')  # Tell user to wait while sending data
        self.setGeometry(50, 50, 1050, 650)
        self.setWindowTitle('PTS')

        mainlayout = QtGui.QVBoxLayout()

        self.widget = QtGui.QWidget()
        self.widget.setLayout(mainlayout)
        self.setCentralWidget(self.widget)  # Set main layout to layout of widget!!!

        # Sub layout 1
        plotBox = QtGui.QHBoxLayout()
        cltBox = QtGui.QHBoxLayout()
        cltBox.setSpacing(30)
        mainlayout.addLayout(plotBox)
        mainlayout.addLayout(cltBox)


        genDataButton = QtGui.QPushButton('Generate Data', self)
        genDataButton.setToolTip('This button randomise a \
            set of data based on the input dimension and cluster amount.')
        genDataButton.resize(genDataButton.sizeHint())
        # genDataButton.move(50, 550)
        genDataButton.clicked.connect(self.genData)
        cltBox.addWidget(genDataButton)
        sendDataButton = QtGui.QPushButton("Send Data", self)
        sendDataButton.resize(genDataButton.sizeHint())
        # sendDataButton.move(200, 550)
        sendDataButton.clicked.connect(self.sendData)
        cltBox.addWidget(sendDataButton)

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

        plotBox.addLayout(plotBoxL)
        plotBox.addLayout(plotBoxR)

        self.show()


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