import curses
import time
from datetime import datetime

from threading import Thread



class ViewEntry(object):
    """A displayable object placed in a list by the Viewer object. Viewer is responslbe for placement, but ViewEntry controls what is rendered.
    Wraps around Simulation and Simulgroup objects and gives a common rendering interface."""
    
    def __init__(self,simul,level = 0):
        self.simul = simul
        if self.simul.typestring == "group":
            self.expandable = True
        else:
            self.expandable = False
        self.expanded = False
        # level at which object is rendered
        self.level = level
        self.force = False

    def render(self,n):
        return "  "*self.level + self.simul.typestring +" [" + str(n).zfill(2) + "] " + self.simul.short_dirname.rjust(10) + " | "  +  self.simul.statestring.rjust(4) + " | " + self.simul.info + "|" + str(self.auto)

    @property
    def auto(self):
        return self.simul.auto

    @auto.setter
    def auto(self,a):
        self.simul.auto = a

    def run(self):
        ret = self.simul.run(force=self.force)
        self.force = False
        return ret


class Viewer(object):

    def __init__(self,simul_list,logfile = "log.txt"):
        """Requires a list of simulgroup or simul elements. Can be a mixed list."""
        if type(simul_list) is not list:
            raise ValueError("Expected list input")
        if len(simul_list) <= 0:
            raise ValueError("Expected list input of length>0")

        new_list = []
        if type(simul_list[0]) is list:
            # flatten list
            for session in simul_list:
                for simul in session:
                    new_list.append(simul)
        else:
            new_list = simul_list    

        self.entry_list = [ViewEntry(nl) for nl in new_list]
        self.Nsimuls = len(self.entry_list)


        # related to how often view is refreshed when no input is given
        self.delta = 5
        
        # print to this file
        self.logfile = logfile

        # properties related to the display area
        ## header area
        self.N_headers = 1
        self.datetime_y = 0
        self.datetime_x = 0    
        self.scoll_indicator_y = 0
        self.scoll_indicator_x = 20

        ## footer area
        self.N_footers = 2
        self.input_x = 0
        self.keys_x = 0



        ## main area    
        self.simul_list_first_y = self.N_headers
        self.simul_list_x = 0
        self.selector_y = self.simul_list_first_y - 1

        

        # properties related to the data being displayed
        ## header
        self.header_cmd = lambda : datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "  " + str(min([self.simul_list_last,self.Nsimuls])) + "/" + str(self.Nsimuls)


        ## footer
        self.keys_string = "up/down: scroll | (e)xpand | (r)elaunch | (a)uto | (q)uit "
        

        self.simul_list_first = 0 # index of first simul to display
        

        # curses magic
        curses.wrapper(self.main)

    def print(self,txt):
        txt = str(txt) + "\n"
        with open(self.logfile,'a') as f:
            f.write(txt)

    def screen_to_list(self,i):
        # simul_list_first: depends on current scroll
        return i -  self.simul_list_first + 1

    def list_to_screen(self,i):
        # simul_list_first: depends on current scroll
        return i + self.simul_list_first - 1

    def get_parent_index(self,n):
        parentlevel = self.entry_list[n].level - 1
        for i in range(n-1,-1,-1):
            this_level = self.entry_list[i].level
            if this_level == parentlevel:
                # found parent
                return i
        return None

    def get_children_indices(self,n):
        parentlevel = self.entry_list[n].level
        j = n
        for i in range(n+1,self.Nsimuls):
            if self.entry_list[i].level == parentlevel:
                j = i
                break
        else:
            j = self.Nsimuls
        return range(n+1,j)
        

    def remove_entry(self,i):
        self.entry_list.pop(i)
        if self.selector_y >=i:
            self.selector_y = self.selector_y - 1
        self.Nsimuls = self.Nsimuls - 1


    def insert_entry(self,i,e):
        self.entry_list.insert(i,e)
        if self.selector_y >=i:
            # this should probably never happen
            # and is thus untested
            self.selector_y = self.selector_y + 1
        self.Nsimuls = self.Nsimuls + 1

    def main(self,stdscr):
        # These properties can only be set in main
        # since they rely on the size of the curses window.
        # footer
        self.input_y = curses.LINES - 1 - 1
        self.keys_y = curses.LINES - 1
        # main area
        self.simul_list_last_y = curses.LINES - 1 - self.N_footers
        self.simul_list_last = self.simul_list_last_y # index of last simul to display
        self.N_simulations_to_display = 1 + self.simul_list_last_y - self.simul_list_first_y
        # footer input
        self.input_string = " " * (curses.COLS-1)

        stdscr.timeout(5000)

        # re-run simulations that need it.
        for entry in self.entry_list:
            if entry.auto:
                entry.run()
                self.print([s.dirname for s in entry.simul.simullist])
            
        while True:
            # Clear screen
            stdscr.clear()
            # print header rows
            stdscr.addstr(self.datetime_y, self.datetime_x, str(self.header_cmd()), curses.A_REVERSE)

            i = self.simul_list_first_y
            n = i + self.simul_list_first - 1
            for entry in self.entry_list[self.simul_list_first:(self.simul_list_last+1)]:
                # if simul.error:
                #     attr=curses.A_STANDOUT
                # else:
                #     attr=curses.A_NORMAL
                if n == self.selector_y:
                    s="*"
                else:
                    s=" "
                stdscr.addstr(i, self.simul_list_x, "[" + s + "] " + entry.render(n))
                i = i + 1
                n = n + 1

            # print footer rows
            stdscr.addstr(self.keys_y, self.keys_x, self.keys_string)
            stdscr.addstr(self.input_y, self.input_x, self.input_string, curses.A_REVERSE)
            stdscr.move(self.input_y,self.input_x)

            try:
                key = stdscr.getkey()
            except:
                curses.error
                key = -1
            if key == "q":
                break
            elif key == "KEY_DOWN":
                if self.selector_y < (self.Nsimuls - 1):
                    self.selector_y =  self.selector_y + 1
                if self.simul_list_last < (self.Nsimuls):
                    self.simul_list_first = self.simul_list_first + 1
                    self.simul_list_last = self.simul_list_last + 1

            elif key == "KEY_UP":
                if self.selector_y > 0:
                    self.selector_y =  self.selector_y - 1
                if self.simul_list_first > 0 :
                    self.simul_list_first = self.simul_list_first - 1
                    self.simul_list_last = self.simul_list_last - 1

            elif key == " ":
                if self.selector_y < (self.Nsimuls - self.N_simulations_to_display - 1):
                    self.selector_y =  self.selector_y  + self.N_simulations_to_display
                elif self.selector_y < (self.Nsimuls - 1):
                    self.selector_y = self.Nsimuls - 1

                if self.simul_list_last < (self.Nsimuls - self.N_simulations_to_display - 1):
                    # jump to next window of simulation
                    self.simul_list_first = self.simul_list_first + self.N_simulations_to_display
                    self.simul_list_last = self.simul_list_last + self.N_simulations_to_display
                elif self.simul_list_last < (self.Nsimuls - 1):
                    # jump to last
                    self.simul_list_first = self.Nsimuls - self.N_simulations_to_display 
                    self.simul_list_last = self.Nsimuls 

            elif key == "v":
                # TODO: visualize
                input_string = stdscr.getstr(self.input_y, self.input_x, curses.COLS-1)
            
            elif key == "e":
                sel = self.entry_list[self.selector_y]
                #self.entry_list[self.selector_y].expand()
                if sel.expandable:
                    if sel.expanded:
                        cs = self.get_children_indices(self.selector_y)
                        for i in cs[::-1]:
                            self.remove_entry(i)
                        sel.expanded = False
                    else:
                        i = self.selector_y
                        child_simuls = sel.simul.simullist
                        Nchildren = len(child_simuls)
                        if Nchildren > 0:
                            for e in [ViewEntry(s,level = sel.level + 1) for s in child_simuls]:
                                i = i + 1
                                self.insert_entry(i,e)
                            sel.expanded = True
                else:
                    if sel.level > 0:
                        p = self.get_parent_index(self.selector_y)
                        cs = self.get_children_indices(p)
                        for i in cs[::-1]:
                            self.remove_entry(i)
                        self.entry_list[p].expanded = False
                
            elif key == "a":
                self.entry_list[self.selector_y].auto = not self.entry_list[self.selector_y].auto
            elif key == "r":
                self.entry_list[self.selector_y].force = not self.entry_list[self.selector_y].force
                self.entry_list[self.selector_y].run()

            elif key == -1:
                # the case when no input is pressed
                # assume user is idle and re-run simulations that need it.
                for entry in self.entry_list:
                    if entry.auto:
                        entry.run()
            
                

if __name__ == "__main__":
    #fake simulation class to mock up output of the real one
    import random
    class simul(object):
        i = 0
        def __init__(self):
            simul.i = simul.i + 1
            self.auto=True
            self.state = 0
            self.force = False

        def run(self,force=False):
            self.state = self.state + 1
            if force:
                self.state = 0

        @property
        def statestring(self):
            if self.state < 2:
                return "OOM?"
            elif self.state <5:
                return "WHAT"
            else:
                return "DONE"

        @property
        def short_dirname(self):
            return "02/Nx10"

        @property
        def info(self):
            return "placeholder info"
        
        @property
        def typestring(self):
            return "group"



    Nsimuls = 100 # a large number
    simul_list = [None]*Nsimuls
    for i in range(Nsimuls):
        simul_list[i] = simul()

    v = Viewer(simul_list)

# def viewer(simul_list):
#     while True:
    
#         for session in simul_list:
#             for i,simul in enumerate(session):
#                 print "=================="
#                 print simul.destination
#                 simul.job_status()
#                 print simul.status
#                 print simul.time
#                 print simul.last_out
#                 simul.copy_to_destination()

#         time.sleep(600) #check every 10 min
        
