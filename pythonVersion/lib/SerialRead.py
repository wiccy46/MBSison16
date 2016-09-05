import serial
import threading

class SerialListener(threading.Thread):
    def __init__(self, threadID, portName, callback):
        self.threadID = threadID
        self.ser = serial.Serial(portName)
        self.flag = True
        threading.Thread.__init__(self)
        self.value = 0
        self.callback = callback

    def run(self):
        while(self.flag == True):
            self.callback(float(self.ser.readline()))

    def stop(self):
        self.flag = False

#
# class GUI():
#     def __init__(self):
#         self.sigma = 0.0
#         self.serialThread = SerialListener("Serial1", '/dev/cu.usbmodem1411')
#
#
#     def run(self):
#         self.serialThread.start()
