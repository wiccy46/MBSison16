import serial
import threading

class SerialListener(threading.Thread):
    def __init__(self, threadID, arduinoPort, callback):
        self.threadID = threadID
        self.ser = serial.Serial(arduinoPort)
        self.flag = True
        threading.Thread.__init__(self)
        self.value = 0
        self.callback = callback

    def run(self):
        while(self.flag == True):
            self.callback(float(self.ser.readline()))

    def stop(self):
        self.flag = False


