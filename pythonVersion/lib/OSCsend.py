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
