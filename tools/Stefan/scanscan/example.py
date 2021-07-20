from event_loop import Abstact_event_loop

class Example_event_loop(Abstact_event_loop):
    def start(self):
        print("started")
        self.i = 0

    def f(self):
        self.i = self.i + 1
        print(self.i)
        if (self.i%3 == 0):
            self.emit_signal("glomp",self.i)
        if (self.i == 1):
            self.emit_signal("lol")

    def glomp(self,a):
        print("glomp" + str(a))
        if a == 6:
            self.emit_signal("lol")

    def l(self):
        print("lol")

if __name__ == "__main__":
    el = Example_event_loop()
    el.event_loop(2.0,{"glomp":el.glomp,"lol":el.l})
