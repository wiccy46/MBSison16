import OSC

class OSCsend():

    def __init__(self, ip, port):
        self.ip = ip
        self.port = port
        self.clientAndroid = OSC.OSCClient();
        self.clientAndroid.connect((self.ip, self.port))

    def osc_msg(self, nr, msg):
        bundle = OSC.OSCBundle()
        #     bundle.setAddress("/resetData")
        #     bundle.append("reset")
        bundle.setAddress("/getNr")
        bundle.append(nr)
        bundle.setAddress("/getData")
        bundle.append(msg)
        self.clientAndroid.send(bundle)
        print "sending OSC data"


# Do this outside of the class
    nr = 100

    sortedData = data[data[:, 0].argsort()]

    for i in range(N / nr):
        osc_msg(nr=nr, msg=sortedData[i * nr: i * nr + nr, 0:2])
        time.sleep(2)

    osc_msg(nr=N % nr, msg=sortedData[N - N % nr: N, 0:2])