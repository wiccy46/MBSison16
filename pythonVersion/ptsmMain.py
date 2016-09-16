import Trajectory
import numpy as np
import time, socket
import threading, sys
import pyo64 as pyo
import matplotlib.pyplot as plt
from lib.Exptable import Exptable
from lib.DSP import DSP
from lib.DataGen import DataGen
from lib.OSCsend import OSCsend
from lib.Listening import Listening
from lib.SerialRead import SerialListener
from PyQt4 import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import serial.tools.list_ports

"""
# Todo : 2. Decouple Send Data and listener. Dedicated listening button.
# TODO : 4. Closing the app needs to kill the thread as well.
# TODo: 1. the serial port name is different some time. Because I have other usb device connected before.
# TODO: 5. The potential plot only reflects the potentials of the 2 columns rather than the whole data set
# Todo: 7: Record sound button .
"""

FS = 44100/4
BLOCK = 1024
HOMEIP = '192.168.178.47'
SELFIP = socket.gethostbyname(socket.getfqdn())
OFFICEIP = "129.70.149.78"
LISTENPORT = 5678
SENDPORT = 7012
socketError = False
MEGA2560 = 'Arduino Leonardo'

"""
Set up Audio server
"""
s = pyo.Server(sr=FS, nchnls=2, buffersize=BLOCK, duplex=0).boot()
s.start()


def linlin(x, smi, sma, dmi, dma): return (x-smi)/float(sma-smi)*(dma-dmi)+dmi
def linlinInvert(x, smi, sma, dmi, dma): return (x-dmi)*(sma-smi)/(dma-dmi) + smi
def findMega(myList):
	for i in range(len(myList)):
		try:
			idx = myList[i].index(MEGA2560)
			return idx
		except ValueError:
			pass


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
        self.firstSend = True # First time sending doesn't need to clear the listner
        self.data, self.pos, self.N, self.dim = np.zeros(2), np.zeros(2), 0, 0
        self.ip = OFFICEIP
        self.genDim, self.genNrmin, self.genNrmax, self.genNc = 3, 50, 200,  4 # Init gen data para
        self.showPotential, self.connectArduino = False, False # Toggle buttons


        self.vel , self.t = 0, 1.
        self.exp_table = np.zeros(1) # Initialise the exp table
        self.norm_max = 0
        self.max_sigma = 0.4 # Use to limit sigma slider based on dim
        self.exp_resolution = 1000000  # Resolution for lookup table


        self.m_comp, self.dt, self.r, self.sigma = 0.02, 0.005, 0.999, 0.15  # init
        self.t = 1.0  # Time in second per piece, TODO a para
        self.blockSize = 5000  # Buffer size for trajectory TODO a para
        self.audioVecSize = int(self.t * self.blockSize)  # I define that 5000 steps will return 1 second of audio
        self.windowing = DSP().makeWindow(self.audioVecSize, rampUp=0.05, rampDown=0.05)  # Windowing for audio
        self.velSound = np.zeros(1)

        self.sigmaSlider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
        self.initUI()

    def androidUpdateSlider(self, value):
        self.dtSlider.setValue(int(value[0]))
        self.rSlider.setValue(int(value[1]))


    def androidUpdateSigma(self,value):
        self.sigmaSlider.setValue(int(value))

    def androidUpdateVelSound(self, value):
        # This function was meant for recording the soundfile. But it is not used.
        print "Velsound callbacked"
        print type(value)
        # self.velSound = value[0]

    def androidUpdateR(self, value):
        self.rSlider.setValue(int(value))


    def serialUpdate(self, value):
        temp =  linlin(value, 460, 490, 0, 100) # 540 is fulls queezed, 590 is relax state
        self.sigmaSlider.setValue(temp)

    # TODO, test whether the listern will also update the self.velSound
    def proc(self, stopevent):
        self.fifo.put(self.velSound)

    def recordSound(self, soundOutput):
        s.recstart(filename="rec1")
        self.fifo.put(soundOutput)
        time.sleep(self.t + 0.1)
        s.recstop()


    def datafigon_pick(self, event):
        ind = np.array(event.ind)  # Get the index of the clicked data
        self.pos = np.array(self.data[ind[0], :])
        self.vel = np.random.rand(self.dim) - 0.5
        trj, junk, forceSound = Trajectory.PTSM(self.pos, self.data, self.vel, self.exp_table, self.exp_resolution, \
                                                self.norm_max, sigma=self.sigma, dt=self.dt, r=self.r, \
                                                Nsamp=self.audioVecSize, compensation=self.m_comp)
        self.velSound = trj[:, 0] / np.max(np.absolute(trj[:, 0])) * self.windowing
        self.velSound = DSP().butter_lowpass_filter(self.velSound, 2000.0, FS, 6)  # 6th order
        # forceSound = forceSound / np.max(np.absolute(forceSound))
        stopevent = threading.Event()
        producer = threading.Thread(name="Compute audio signal", target=self.proc, args=[stopevent])
        producer.start()
        self.drawTrj(trj)

    def drawTrj(self, trj): # Plot trajectory
        self.trjFigAx.clear()
        potmap = np.zeros((40, 40))
        if self.showPotential:
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

        self.specFigAx.specgram(self.velSound, NFFT = 1024, Fs = FS, noverlap = 900, cmap= plt.cm.gist_heat)
        # self.specFigAx.set_yticks(np.arange(0 , max(self.specFigAx.get_yticks()/1000), 1))
        # print self.specFigAx.get_yticks()
        # print type(self.specFigAx.get_yticks())
        self.trjFigCanvas.draw()

    def sendData(self):
        global  socketError
        # If there is a listener already stop it.]
        print "Send data button clicked."
        if (self.firstSend):
            print "listener not created yet. "
            self.firstSend = False
        else:

            print "stopping listener."
            self.androidListener.stop()
            time.sleep(0.5)

        print "sendData called"
        self.statusBar().showMessage('Sending data, please wait ...')
        time.sleep(0.3)
        if self.N != 0:
            androidClient = OSCsend(self.ip, SENDPORT)
            nr = 100
            sortedData = self.data[self.data[:, 0].argsort()]
            for i in range(self.N / nr):
                androidClient.osc_msg(nr=nr, msg=sortedData[i * nr: i * nr + nr, 0:2])
                time.sleep(0.5)
            androidClient.osc_msg(nr=self.N % nr, msg=sortedData[self.N - self.N % nr: self.N, 0:2])
            # Start listening
            self.androidListener = Listening(gui = self, ip = SELFIP, sliderCallback=self.androidUpdateSlider,
                                             sigmaSliderCallback=self.androidUpdateSigma, rSliderCallback= self.androidUpdateR,
                                             velSoundCallback= self.androidUpdateVelSound, port = LISTENPORT)
            self.androidListener.spawn()
            if (socketError):
                print "socketError, server already exist. "
                self.statusBar().showMessage('Error, Click "Clear Listener" and try again.')
            else:
                self.androidListener.add_handler()
                self.androidListener.start()
                self.statusBar().showMessage('Finished sending, created listener.')

        else:
            self.statusBar().showMessage('Data not generated yet, please generate data first.')

    def genData(self):
        self.dataFigAx.clear()
        if (self.genDim > 2):
            self.data = DataGen().datagen(self.genDim, self.genNc, sigma=0.4, minnr=self.genNrmin, maxnr=self.genNrmax)
        elif (self.genDim == 2):
            self.data  = np.random.rand((self.genNrmax + self.genNrmin)/2 * self.genNc ,2)-[0.5, 0.5]
        else:
            self.statusBar().showMessage("ERROR! Data dimension needs to be greater than 1")
        self.N, self.dim = self.data.shape[0], self.data.shape[1]
        self.norm_max = self.dim
        self.max_sigma = np.log(self.dim) / 2.5 * 0.5 + 0.1
        delta_dis = np.linspace(0.0, self.dim, num=self.exp_resolution)  # Range of distance difference.
        self.exp_table = Exptable().updateExpTable(self.sigma, delta_dis)
        self.dataFigAx.scatter(self.data[:, 0], self.data[:, 1], c="blue", picker=2, s=[50] *self.N)
        self.dataFigAx.set_xlim([-0.5, 0.5])
        self.dataFigAx.set_ylim([-0.5, 0.5])
        self.dataFigCanvas.draw()
        self.statusBar().showMessage('Data generated. Click on the data on the left graph to probe a particle trajectory.')

    def clearListener(self):
        try:
            self.androidListener.stop()
            time.sleep(0.3)
            self.androidListener.stop()
            self.statusBar().showMessage("Clear Listener")
        except AttributeError:
            pass



    def getfiles(self):
        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.AnyFile)
        dlg.setFilter("Text files (*.txt *csv)") # And then it depends whether txt or csv do thing differently
        filenames = QtCore.QStringList()
        if dlg.exec_():
            filenames = dlg.selectedFiles()
            f = open(filenames[0], 'r')
            with f:
                temp = f.read()
                self.contents.setText(temp) # send data to data here. but might need processing.


    def printDataInfo(self):
        if(self.data):
            self.statusBar().showMessage("Data: Row: " + str(self.N) + " , Dim: " + str(self.dim) + " , Clusters: "
                                         + str(self.genNc) + ".")
        else:
            self.statusBar().showMessage("No Data!")

    def toggleButtons(self, pressed):
        source = self.sender()  # could be shared with other similar toggleButton
        if pressed: temp = True
        else: temp = False
        if source.text() == "Show Potential":
            self.showPotential = temp
        elif source.text() == "Squeeze Ball":
            self.connectArduino = temp
            if (self.connectArduino == True):
                print "Create Serial Connection"
                self.portList = list(serial.tools.list_ports.comports())
                result = findMega(self.portList)
                if result:
                    self.serialThread = SerialListener("Serial1", self.portList[result][0], self.serialUpdate)
                    self.serialThread.start()
                    self.statusBar().showMessage("Serial Connection Success.")
                else:
                    self.statusBar().showMessage("Connection Failed. Squeeze ball is not connected to the computer.")
            else:
                self.serialThread.stop()
        else: pass



    def initUI(self):
        self.statusBar().showMessage('Move on each item to see user tip.')  # Tell user to wait while sending data
        self.setGeometry(50, 50, 1050, 650)
        self.setWindowTitle('Particle Trajectory Sonification')

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
        # self.dataFigAx.set_title("Pick a data to as a particle")
        self.dataFigToolbar = NavigationToolbar(self.dataFigCanvas,self)

        plotBoxL = QtGui.QVBoxLayout()
        plotBoxL.setSpacing(5)
        plotBoxL.addWidget(self.dataFigCanvas)
        plotBoxL.addWidget(self.dataFigToolbar)

        self.trjFig = plt.figure()
        self.trjFigCanvas = FigureCanvas(self.trjFig)
        self.trjFigAx = self.trjFig.add_subplot(121)
        self.trjFigAx.set_xlim([-0.5, 0.5])
        self.trjFigAx.set_ylim([-0.5, 0.5])
        self.trjFigAx.set_title('Trajectory')

        self.specFigAx = self.trjFig.add_subplot(122)
        self.specFigAx.set_title('Spectrogram')
        self.specFigAx.set_ylabel("Frequency (Hz)")

        self.specFigAx.set_xlabel("Time (s)")
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

        openFileButton = QtGui.QPushButton("Load File")
        openFileButton.clicked.connect(self.getfiles)

        printInfoButton = QtGui.QPushButton("Info")
        printInfoButton.setToolTip('Display data information.')
        printInfoButton.clicked.connect(self.printDataInfo)

        sendDataButton = QtGui.QPushButton("Send Data", self)
        sendDataButton.clicked.connect(self.sendData)
        sendDataButton.setToolTip('Send the first 2 columns of data to Android device for plotting.')

        clearListenerButton = QtGui.QPushButton('Clear Listener', self)
        clearListenerButton.setToolTip('Click to clear the osc listener.')
        clearListenerButton.clicked.connect(self.clearListener)

        recButton = QtGui.QPushButton("Record", self)
        recButton.setToolTip('Record the last particle trajectory sound, require a particle.')
        recButton.clicked.connect(self.recordSound)

        self.showPotentialButton = QtGui.QPushButton("Show Potential", self)
        self.showPotentialButton.setCheckable(True)
        self.showPotentialButton.clicked[bool].connect(self.toggleButtons)



        #----------------------
        #sigma
        sigmaTitle = QtGui.QLabel('Sigma')
        sigmaTitle.resize(sigmaTitle.sizeHint())

        self.sigmaSlider.setToolTip('Large sigma results in 1 global minimum while small sigma results in local minimum '
                               'in each data point')
        self.sigmaSlider.setRange(0, 100)
        self.sigmaSlider.setValue(50)
        self.sigmaSlider.valueChanged[int].connect(self.changeValue)
        # self.sigmaSlider.valueChanged[int].connect(self.serialUpdate)
        self.sigmaDisplay = QtGui.QLineEdit()
        self.sigmaDisplay.setFixedSize(50, 20)
        self.sigmaDisplay.setText('0.197') # Using hard input is not clever.

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
        self.ipDisplay.setText(self.ip)
        # self.ipDisplay.resize(self.dimDisplay.sizeHint())
        self.ipDisplay.setFixedSize(106, 20)
        self.ipDisplay.returnPressed.connect(self.ipChangeValue)
        # self.ipDisplay.valueChanged[str].connect(self.changeValue)


        self.connectArduinoButton = QtGui.QPushButton("Squeeze Ball", self)
        self.connectArduinoButton.setCheckable(True)
        self.connectArduinoButton.clicked[bool].connect(self.toggleButtons)


        cltLeftBox = QtGui.QGridLayout()
        cltLeftBox.addWidget(ncTitle,1 , 0), cltLeftBox.addWidget(self.ncDisplay, 1, 1)
        cltLeftBox.addWidget(dimTitle, 1, 2), cltLeftBox.addWidget(self.dimDisplay, 1, 3)
        cltLeftBox.addWidget(nrminTitle, 2, 0), cltLeftBox.addWidget(self.nrminDisplay, 2, 1)
        cltLeftBox.addWidget(nrmaxTitle, 2, 2), cltLeftBox.addWidget(self.nrmaxDisplay, 2, 3)
        cltLeftBox.addWidget(genDataButton, 3, 0 ), cltLeftBox.addWidget(openFileButton, 3, 1) , \
            cltLeftBox.addWidget(printInfoButton, 3, 2), cltLeftBox.addWidget(self.showPotentialButton, 3, 3)
        cltLeftBox.addWidget(ipTitle, 4,0), cltLeftBox.addWidget(self.ipDisplay, 4, 1), cltLeftBox.addWidget(self.connectArduinoButton, 4, 2)
        cltLeftBox.addWidget(sendDataButton, 5, 0, 5, 2), cltLeftBox.addWidget(clearListenerButton, 5, 2), cltLeftBox.addWidget(recButton, 5, 3)


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
        plotBox.addLayout(plotBoxL, 3)
        plotBox.addLayout(plotBoxR, 4)

        cltBox = QtGui.QHBoxLayout()

        cltBox.addLayout(cltLeftBox)
        cltBox.addLayout(cltRightBox)

        mainlayout.addLayout(plotBox, 8)
        mainlayout.addLayout(cltBox, 2)
        self.show()

    def ipChangeValue(self):
        self.ip = self.ipDisplay.text()
        self.statusBar().showMessage('Android IP: ' + self.ip)


    def changeValue(self,value):
        if self.sender() == self.sigmaSlider:
            smi, sma = 0, 100
            dmi, dma = 0.001, self.max_sigma
            self.sigma = linlin(value, smi, sma, dmi, dma)
            self.sigmaDisplay.setText(str(format(self.sigma, '.3f')))

        elif self.sender() == self.dtSlider:
            smi, sma = 0, 100
            dmi, dma = 0.001, 0.01
            self.dt = linlin(value, smi, sma, dmi, dma)
            self.dtDisplay.setText(str(format(self.dt, '.3f')))

        elif self.sender() == self.rSlider:
            smi, sma = 0, 100
            dmi, dma = 0.99, 1.0
            self.r = linlin(value, smi, sma, dmi, dma)
            self.rDisplay.setText(str(format(self.r, '.3f')))

        elif self.sender() == self.dimDisplay:
            print "dimension = " + str(value)
            self.genDim = value
        elif self.sender() == self.ncDisplay:
            print "nc = " + str(value)
            self.genNc = value
        elif self.sender() == self.nrminDisplay:
            if (value > self.genNrmax):
                self.statusBar().showMessage('Warning, input value higher than Nr Max.')
                self.genNrmin = self.genNrmax
                self.nrminDisplayDisplay.setValue(self.genNrmin)
            else:
                print "nr min = " + str(value)
                self.genNrmin = value
        elif self.sender() == self.nrmaxDisplay:
            if (value < self.genNrmin):
                self.statusBar().showMessage('Warning, input value lower than Nr Min.')
                self.genNrmax = self.genNrmin
                self.nrminDisplayDisplay.setValue(self.genNrmax)
            else:
                print "nr max = " + str(value)
                self.genNrmax = value
        else: print "Unrecognised input source."

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