import numpy as np
from lib.DSP import DSP
import threading, OSC, socket
import Trajectory
FS = 44100/4

# Hasn't tested yet.
class Listening(object):
    def __init__(self, gui,  ip , callback, port=5678):
        self.gui = gui
        self.receive_address = ip, port
        self.callback = callback

    def printpara(self):
        print self.gui.sigma


    def trigger_handler(self, addr, tags, stuff, source):
        print "received trigger."
        self.sortedPos = np.array([stuff[0], stuff[1]])
        # # Find nearest neighbour
        # # TODO, exclude playing sound if there is no neighbouring point < dis_thres
        # # if nearest d < dis_thres, play , else print "no neighbour"
        #
        deltas = self.gui.data[:, 0:2] - self.sortedPos
        dist_2 = np.einsum('ij,ij->i', deltas, deltas)
        idx = np.argmin(dist_2)
        vel = np.random.rand(self.gui.dim) - 0.5  # Need to be controllable by pressure
        pos = np.array(self.gui.data[idx, :])
        #

        trj, junk, forceSound = Trajectory.PTSM(pos, self.gui.data, vel, self.gui.exp_table,
                                                self.gui.exp_resolution, self.gui.norm_max, sigma=self.gui.sigma,
                                                dt=self.gui.dt, r=self.gui.r,
                                                Nsamp=self.gui.audioVecSize, compensation=self.gui.m_comp)
        velSound = trj[:, 0] / np.max(np.absolute(trj[:, 0])) * self.gui.windowing
        velSound = DSP().butter_lowpass_filter(velSound, 2000.0, FS, 6)  # 6th order
        # forceSound = forceSound / np.max(np.absolute(forceSound))
        stopevent = threading.Event()
        producer = threading.Thread(name="Compute audio signal", target=self.gui.proc,
                                    args=[stopevent, velSound])
        producer.start()
        # Try to draw trj here. might fail though
        self.gui.drawTrj(trj)

    def slider_handler(self, addr, tags, stuff, source):
        self.callback(np.array([stuff[0], stuff[1]]))



    def spawn(self):
        global  socketError
        try:
            self.receiveServer = OSC.OSCServer(self.receive_address)  # create a serve to receive OSC from the tablet
            self.receiveServer.addDefaultHandlers()
            print"Server Created."
        except socket.error:
            print "Socket Error."
            socketError = True
            print socketError

    def stop_handler(self, addr, tags, stuff, source):
        # Close the OSC server
        print "\nClosing OSCServer."
        self.receiveServer.close()
        print "Waiting for Server-thread to finish"
        try:
            self.emorating_oscServer.join()  ##!!!
            print "Done"
        except AttributeError:
            print "Done"

    def add_handler(self):
        try:
            self.receiveServer.addMsgHandler("/trigger", self.trigger_handler)
            self.receiveServer.addMsgHandler("/stop", self.stop_handler)
            self.receiveServer.addMsgHandler("/sliders", self.slider_handler)
        except AttributeError:
            pass

    def start(self):
        try:
            self.emorating_oscServer = threading.Thread(target=self.receiveServer.serve_forever)
            self.emorating_oscServer.start()
            print "\nOSCServer established."
        except AttributeError:
            pass

    def stop(self):
        # Close the OSC server
        print "\nClosing OSCServer."
        self.receiveServer.close()
        print "Waiting for Server-thread to finish"
        try:
            self.emorating_oscServer.join()  ##!!!
            print "Done"
        except AttributeError:
            print "Done"