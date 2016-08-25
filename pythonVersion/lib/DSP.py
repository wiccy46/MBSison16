import numpy as np
from scipy.signal import butter, lfilter

# TODO all audio stuff in a separate class.
class DSP():
    def butter_lowpass(self, cutoff, fs, order=5):
        nyq = 0.5 * fs
        normal_cutoff = cutoff / nyq
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        return b, a

    def butter_lowpass_filter(self, data, cutoff, fs, order=5):
        b, a = self.butter_lowpass(cutoff, fs, order=order)
        y = lfilter(b, a, data)
        return y

    # Create a rigid window to ramp on/off the audio vector to prevent clipping.
    # But I still hear clipping.
    def makeWindow(self, s, rampUp=0.02, rampDown=0.02):
        w = np.append(np.linspace(0, 1.0, num=int(s * rampUp)), \
                      np.ones(s - int(s * rampUp) - int(s * rampDown)))
        w = np.append(w, np.linspace(1.0, 0, num=int(s * rampDown)))
        return w
