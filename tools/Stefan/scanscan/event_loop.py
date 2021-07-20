import time

class Abstact_event_loop:

    # can probably be replaced by a library that implements
    # a generic event loop

    # This class should be inherited and extended by f()

    def f(self):
        """This function is called every iteration of the event loop.
        It should be overriden by a concrete function"""
        raise NotImplementedError("Please implement f()")

    def start(self):
        pass

    def __init__(self,filehandle=None):
        """ Class for running a function every delta seconds 
        delta -- the time between each iteration
        f -- function to be run each iteration
        start -- function to run before starting the loop
        signals -- dictionary of special events and responses, see below.

        signals:
        signals are a dictionary on the form
        {string : function}
        where the string is a name of a signal and the function is executed
        before the next loop if the named signal has been emitted.
        """
        
        self.filehandle = filehandle
        self.emitted_signals = []

    def print(self,s):
        if self.filehandle is None:
            print(s)
        else:
            self.filehandle.write(str(s) + "\n")
            self.filehandle.flush() # to immediately update file

    def event_loop(self,delta, signals):
        self.signals = signals
        self.delta = delta
        
        starttime = time.time()
        self.start()
        while True:
            self.consume_signals()
            self.f()
            time.sleep(self.delta - ((time.time() - starttime) % self.delta))

    def emit_signal(self,signal_name,*args):
        if signal_name in self.signals:
            self.emitted_signals.append(lambda : self.signals[signal_name](*args))
        else:
            raise ValueError("Unrecongized signal: '" + signal_name + "'")

    def consume_signals(self):
        # copy and clear signals so that new signals can be emitted in the
        # functions called as a result of the previous signals
        _signals = self.emitted_signals
        self.emitted_signals = []
        for s in _signals:
            s()
