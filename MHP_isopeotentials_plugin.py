#! /usr/bin/env python

#AUTHOR
#  Jean-Marc Crowet (jmcrowet@ulg.ac.be); 23/11/2010

"""
DESCRIPTION

  This plugin computes the MHP (Mean Hydrophobicity Potential) values 
  for the specified pdb file and draw the isopotentials in PyMol.
  MHP values are computed for each positions on a user defined grid and 
  saved in a *.dx file. 
  MHP isopotentials are distant from the molecular surface so you have to 
  define a distance between the molecule and the box walls.
  
USAGE 

  Through the interface

     If you have not a *.dx file yet, use the MHP Compute procedure 
     of this plugin. Fill in the parameters (default values are 
     generally OK) and click on OK.
     MHP computaions take some time so if you have *.dx file, 
     you should use the MHP Display procedure.

  On the command line

     MHP_pymol *.pdb, *.dx, [d, delta, iso_phi, iso_pho]

        Terms in brackets are optional and their default values 
        are d = 5, delta = 1, iso_phi = 0.1, iso_pho = -0.1
        d correspond to the distance between the solute and the box (A)
        delat correspond to the mesh size (A)
        iso_phi and iso_pho correspond to the isopotential values

     MHP_display *.dx [iso_phi, iso_pho]

        MHP_display draw the isopotentials in PyMol from the *.dx file.

     MHP_Help

Created by Jean-Marc Crowet (jmcrowet@ulg.ac.be)
CBMN - University of Liege (www.fsagx.ac.be/bp/)

"""

import Pmw
import Tkinter
from Tkinter import *
import tkFileDialog
import tkMessageBox
from pymol import cmd
import string, os, re, math, sys

def __init__(self):
    self.menuBar.addcascademenu('Plugin', 'MHP', 'MHP', label = 'MHP isopotentials')
    self.menuBar.addmenuitem('MHP', 'command', 'MHP_Compute',
        label = 'MHP Compute', command = lambda s=self: MHP_Compute(s))
    self.menuBar.addmenuitem('MHP', 'command', 'MHP_Display',
        label = 'MHP Display', command = lambda s=self: MHP_Display(s))
    self.menuBar.addmenuitem('MHP', 'command', 'Help',
        label = 'Help', command = lambda s=self: MHP_doc(s))

Rvdw = { 'H' : 1.0, 'O' : 1.4, 'N' : 1.55, 'C' : 1.7, 'S' : 1.8, 'P' : 1.8}
Rbond = { 'H' : 0.5, 'O' : 0.7, 'N' : 0.775, 'C' : 0.85, 'S' : 1.1, 'P' : 0.9}
defaults = { 'd' : 5.0, 'delta' : 1.0, 'iso_phi' : 0.1, 'iso_pho' : -0.1}

class MHP_Compute:

    def __init__(self, app):
        self.parent = app.root
        self.dialog = Pmw.Dialog(self.parent, buttons = ('OK','Cancel'),
            title = 'PyMOL MHP Compute - Parameters', command = self.execute) 
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        frame = Tkinter.Frame(self.dialog.interior())
        frame.pack(fill='both',expand=1)
        self.entry_f = Pmw.EntryField(frame,labelpos='w',
              label_text = 'PDB filename : ', value = '')
        self.entry_f.pack(side='left',fill='x',expand=1, padx=10)
        self.buttonBox = Pmw.ButtonBox(frame, labelpos = None)
        self.buttonBox.pack(expand=1)
        self.buttonBox.add('...', command = self.getFileName)
        self.entry_d = Pmw.EntryField(self.dialog.interior(),labelpos='w',
              label_text = 'Distance between the solute and the box (A) : ',
              value = str(5.0), validate = {'validator' : 'real','min':0,})
        self.entry_d.pack(fill = 'x', padx = 10, pady = 10)
        self.entry_delta = Pmw.EntryField(self.dialog.interior(),labelpos='w',
              label_text = 'Mesh size (A) : ', value = str(1.0),
              validate = {'validator' : 'real','min':0,})
        self.entry_delta.pack(fill = 'x', padx = 10, pady = 10)
        self.entry_iso_phi = Pmw.EntryField(self.dialog.interior(),labelpos='w',
              label_text = 'Hydrophilic isopothential value : ',
              value = str(0.1), validate = {'validator' : 'real',})
        self.entry_iso_phi.pack(fill = 'x', padx = 10, pady = 10)
        self.entry_iso_pho = Pmw.EntryField(self.dialog.interior(),labelpos='w',
              label_text = 'Hydrophobic isopotential value : ',
              value = str(-0.1), validate = {'validator' : 'real',})
        self.entry_iso_pho.pack(fill = 'x', padx = 10, pady = 10)
        frame2 = Tkinter.Frame(self.dialog.interior())
        frame2.pack(fill='both',expand=1)
        self.entry_dx = Pmw.EntryField(frame2,labelpos='w',
              label_text = 'Output dx file as : ', value = '')
        self.entry_dx.pack(side='left',fill='x',expand=1, padx=10)
        self.buttonBox2 = Pmw.ButtonBox(frame2, labelpos = None)
        self.buttonBox2.pack(expand=1)
        self.buttonBox2.add('...', command = self.getSaveAsName)

        self.dialog.show()

    def execute(self, result):
        if result == 'OK':
            file_name = self.entry_f.getvalue()
            d = self.entry_d.getvalue()
            delta = self.entry_delta.getvalue()
            iso_phi = self.entry_iso_phi.getvalue()
            iso_pho = self.entry_iso_pho.getvalue()
            dx_file = self.entry_dx.getvalue()
            dir_name = os.path.dirname(dx_file)
            dx_name = os.path.basename(dx_file)
            dx_name_root = os.path.splitext(dx_name)[0]
            dx_name_ext = os.path.splitext(dx_name)[1]
            if os.path.exists(file_name) and os.path.isdir(dir_name) and dx_name_root != '' and dx_name_ext == '.dx' and d != '' and delta != '' and iso_phi != '' and iso_pho != '':
                MHP_pymol(file_name,dx_file,d,delta,iso_phi,iso_pho)
            else:
                tkMessageBox.showerror('Invalid path or empty field', 
                   'You must enter valid path names and fill in all fields')
        else:
            if __name__ == '__main__':
                self.parent.destroy()
            else:
                self.dialog.withdraw()

    def getFileName(self):
        file_types = [ ("PDB files", "*.pdb", "TEXT"),
                      ("All files", "*"),]
        dialog = tkFileDialog.Open(filetypes=file_types)
        fname = dialog.show()
        if fname != "":
            entry_f_value = fname
            self.entry_f.setvalue(fname)

    def getSaveAsName(self):
        file_types = [ ("dx files", "*.dx", "TEXT")]
        dialog = tkFileDialog.SaveAs(filetypes=file_types)
        fname = dialog.show()
        if fname != "":
            self.entry_dx.setvalue(fname)


class MHP_Display:

    def __init__(self, app):
        self.parent = app.root
        self.dialog = Pmw.Dialog(self.parent, buttons = ('OK','Cancel'),
            title = 'PyMOL MHP Display - Parameters', command = self.execute) 
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        frame = Tkinter.Frame(self.dialog.interior())
        frame.pack(fill='both',expand=1)
        self.entry_f = Pmw.EntryField(frame,labelpos='w',
              label_text = 'dx filename : ', value = '')
        self.entry_f.pack(side='left',fill='x',expand=1, padx=10)
        self.buttonBox = Pmw.ButtonBox(frame, labelpos = None)
        self.buttonBox.pack(expand=1)
        self.buttonBox.add('...', command = self.getFileName)
        self.entry_iso_phi = Pmw.EntryField(self.dialog.interior(),labelpos='w',
              label_text = 'Hydrophilic isopothential value : ',
              value = str(0.1), validate = {'validator' : 'real',})
        self.entry_iso_phi.pack(fill = 'x', padx = 10, pady = 10)
        self.entry_iso_pho = Pmw.EntryField(self.dialog.interior(),labelpos='w',
              label_text = 'Hydrophobic isopotential value : ',
              value = str(-0.1), validate = {'validator' : 'real',})
        self.entry_iso_pho.pack(fill = 'x', padx = 10, pady = 10)

        self.dialog.show()

    def execute(self, result):
        if result == 'OK':
            file_name = self.entry_f.getvalue()
            iso_phi = self.entry_iso_phi.getvalue()
            iso_pho = self.entry_iso_pho.getvalue()
            if os.path.exists(file_name) and iso_phi != '' and iso_pho != '':
                MHP_display(file_name,iso_phi,iso_pho)
            else:
                tkMessageBox.showerror('Invalid path or empty field', 
                   'You must enter valid path names and fill in all fields')
        else:
            if __name__ == '__main__':
                self.parent.destroy()
            else:
                self.dialog.withdraw()

    def getFileName(self):
        file_types = [ ("dx files", "*.dx", "TEXT"),
                      ("All files", "*"),]
        dialog = tkFileDialog.Open(filetypes=file_types)
        fname = dialog.show()
        if fname != "":
            self.entry_f.setvalue(fname)

class MHP_doc:

    def __init__(self, app):
        self.parent = app.root
        self.dialog = Pmw.Dialog(self.parent, buttons = ('OK',),
            title = 'PyMOL MHP Help', command = self.execute) 
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        bar = Scrollbar(self.dialog.interior())
        text = Text(self.dialog.interior(),yscrollcommand=bar.set)
        bar.config(command=text.yview)
        text.insert(END,__doc__)
        text.pack(side=LEFT,expand="no",fill="both")
        bar.pack(side=LEFT,expand="no",fill="y")

        self.dialog.show()

    def execute(self, result):
        if __name__ == '__main__':
            self.parent.destroy()
        else:
            self.dialog.withdraw()


def MHP_Help():

    print __doc__

cmd.extend('MHP_Help', MHP_Help)


def load_pdb(file_name,d,delta):

    h_match = re.compile('[0-9]?[HhDd].*')
    o_match = re.compile('^[Oo].*')
    n_match = re.compile('^[Nn].*')
    c_match = re.compile('^[Cc].*')
    s_match = re.compile('^[Ss].*')
    p_match = re.compile('^[Pp].*')

    x = []
    y = []
    z = []
    a = []
    Etr = []
    Ri = []
    x_min=y_min=z_min=99999
    x_max=y_max=z_max=-99999
    pdb_lines = open(file_name).readlines()
    for line in pdb_lines:
        if line[0:4] == 'ATOM' or line[0:4] == 'HETA':
            aname = line[12:17].strip()
            resname = line[17:21].strip()
            x_temp = float(line[30:39])
            x.append(x_temp)
            y_temp = float(line[38:47])
            y.append(y_temp)
            z_temp = float(line[46:55])
            z.append(z_temp)
            
            if o_match.match(aname):
                a.append('O')
            elif n_match.match(aname):
                a.append('N')
            elif c_match.match(aname):
                a.append('C')
            elif s_match.match(aname):
                a.append('S')
            elif p_match.match(aname):
                a.append('P')
            elif h_match.match(aname):
                a.append('H')
            else:
                a.append('Unknow_Atom')
                print 'Unknow_Atom:' +  str(len(x)+1) + ' ' + resname + ' ' + aname           

            if x_temp <= x_min:
                x_min = x_temp
            if x_temp >= x_max:
                x_max = x_temp
            if y_temp <= y_min:
                y_min = y_temp
            if y_temp >= y_max:
                y_max = y_temp
            if z_temp <= z_min:
                z_min = z_temp
            if z_temp >= z_max:
                z_max = z_temp

    x_min = x_min - d
    x_max = x_max + d
    y_min = y_min - d
    y_max = y_max + d
    z_min = z_min - d
    z_max = z_max + d

    dx = (int((x_max - x_min)/delta)+1) * delta - (x_max - x_min)
    dy = (int((y_max - y_min)/delta)+1) * delta - (y_max - y_min)
    dz = (int((z_max - z_min)/delta)+1) * delta - (z_max - z_min)

    x_max += dx
    y_max += dy
    z_max += dz

    return x, y, z, a, x_min, x_max, y_min, y_max, z_min, z_max


def pdb_connectivity(x, y, z, a):

    Connect = [' ']*len(x)
    NbC = [0]*len(x)
    for i in range (0,len(x),1):
        for j in range (i+1, len(x),1):
            Test = 'OK'
#            if a[i] == 'H' and NbC[i] == 1:
#                Test = 'KO'
#            elif a[i] == 'O' and NbC[i] == 2:
#                Test = 'KO'
#            elif a[i] == 'N' and NbC[i] == 4:
#                Test = 'KO'
#            elif a[i] == 'C' and NbC[i] == 4:
#                Test = 'KO'
#            elif a[i] == 'S' and NbC[i] == 2:
#                Test = 'KO'
#            elif a[i] == 'P' and NbC[i] == 4:
#                Test = 'KO'
            if Test == 'OK':
                Dist = ((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2 + ( z[i] - z[j]) ** 2) ** 0.5
                Dist_vdw = Rbond[a[i]] + Rbond [a[j]]
                if Dist < Dist_vdw:
                    Connect[i] = Connect[i] + a[j] + ' '
                    Connect[j] = Connect[j] + a[i] + ' '
                    NbC[i] += 1
                    NbC[j] += 1

    return Connect


def atom_properties(Connect, a):

    Atom_Etr = { 'HC' : -0.537, 'Hq' : 1.03, 'O' : 2.833, 'N' : 3.035, 'Csp2' : -1.513, 'Csp3' : -2.436, 'CH3' : -4.047, 'S' : -2.751, 'P' : 1.0 }

    Etr = []
    for i in range (0,len(a),1): 
        data = Connect[i].strip().split()   
        if a[i] == 'H':
            if data[0] == 'C' or data[0] == 'S':
                Etr.append(-0.537)
            else:
                Etr.append(1.03)
        elif a[i] == 'O':
            Etr.append(2.833)
        elif a[i] == 'N':
            Etr.append(3.035)
        elif a[i] == 'C':
            if len(data) <= 3:
                Etr.append(-1.513)
            else:
                k = 0
                for j in range(1,len(data),1):
                    if data[j] == 'H':
                        k += 1
                if k >= 3:
                    Etr.append(-4.047)
                else:
                    Etr.append(-2.436)
        elif a[i] == 'S':
            Etr.append(-2.751)
        elif a[i] == 'P':
            Etr.append(1.0)
        else:
            Etr.append(0.0)

    return Etr


def compute_mhp(file_name, x, y, z, a, Etr, x_min, x_max, y_min, y_max, z_min, z_max, delta):

    x_nb = int((x_max - x_min) / delta) + 1
    y_nb = int((y_max - y_min) / delta) + 1
    z_nb = int((z_max - z_min) / delta) + 1
    nb = x_nb * y_nb * z_nb

    dx_data = open(file_name, 'w')

    dx_data.write("object 1 class gridpositions counts %d %d %d\n" % (x_nb,y_nb,z_nb))
    dx_data.write("origin %f %f %f\n" % (x_min,y_min,z_min))
    dx_data.write("delta %f 0 0\n" % (delta))
    dx_data.write("delta 0 %f 0\n" % (delta))
    dx_data.write("delta 0 0 %f\n" % (delta))
    dx_data.write("object 2 class gridconnections counts %d %d %d\n" % (x_nb,y_nb,z_nb))
    dx_data.write("object 3 class array type double rank 0 items %d data follows\n" % (nb))

    j = 0
    k = 0
    MHP = [0]*3
    for xgi in range(0,x_nb,1):
        xg = x_min + xgi * delta
        for ygi in range(0,y_nb,1):
            yg = y_min + ygi * delta
            for zgi in range(0,z_nb,1):
                zg = z_min + zgi * delta
                MHP[j] = 0
                for i in range (0,len(x),1):
                    if (xg - x[i]) < 10 or (yg - y[i]) < 10 or (zg - z[i]) < 10:
                        Di = ( (xg - x[i]) ** 2 + (yg - y[i]) ** 2 + (zg - z[i]) ** 2 ) ** 0.5
                        MHP[j] += Etr[i] * math.exp( Rvdw[a[i]] - Di )
                if j == 2:
                    dx_data.write("%e %e %e\n" % (float(MHP[0]),float(MHP[1]),float(MHP[2])))
                j += 1
                if j == 3:
                    j = 0
                if xgi == x_nb-1 and ygi == y_nb-1 and zgi == z_nb-1:
                    if j == 1:
                        dx_data.write("%e\n" % (float(MHP[0])))
                    if j == 2:
                        dx_data.write("%e %e\n" % (float(MHP[0]),float(MHP[1])))
                k += 1
                temp = 100 * k / nb
                sys.stdout.write(str(int(temp)) + " %\r")

    dx_data.write("attribute \"dep\" string \"positions\"\n")
    dx_data.write("object \"regular positions regular connections\" class field\n")
    dx_data.write("component \"positions\" value 1\n")
    dx_data.write("component \"connections\" value 2\n")
    dx_data.write("component \"data\" value 3\n")

    dx_data.close()


def MHP_pymol(file_name,dx_file,d=5.0,delta=1.0,iso_phi=0.1,iso_pho=-0.1):

    d = float(d)                #Distance between the solute and the box
    delta = float(delta)        #Mesh size
    iso_phi = float(iso_phi)  #Isopotential for the hydrophilic surface
    iso_pho = float(iso_pho)  #Isopotential for the hydrophobic surface

    x, y, z, a, x_min, x_max, y_min, y_max, z_min, z_max = load_pdb(file_name,d,delta)
    Connect = pdb_connectivity(x, y, z, a)
    Etr = atom_properties(Connect, a)
    compute_mhp(dx_file, x, y, z, a, Etr, x_min, x_max, y_min, y_max, z_min, z_max, delta)
    MHP_display(dx_file,iso_phi,iso_pho)

cmd.extend('MHP_pymol', MHP_pymol)


def MHP_display(dx_file,iso_phi=0.1,iso_pho=-0.1):

    dx_file_name = os.path.basename(dx_file)
    Obj_name = os.path.splitext(dx_file_name)[0]
    cmd.load(dx_file,Obj_name)
    Obj_name_phi = Obj_name + "_phi"
    cmd.isosurface(Obj_name_phi, Obj_name, iso_phi)
    cmd.color('green', Obj_name_phi)
    Obj_name_pho = Obj_name + "_pho"
    cmd.isosurface(Obj_name_pho, Obj_name, iso_pho)
    cmd.color('orange', Obj_name_pho)

cmd.extend('MHP_display', MHP_display)


# Create demo in root window for testing.
if __name__ == '__main__':
    class App:
        def my_show(self,*args,**kwargs):
            pass
    app = App()
    app.root = Tkinter.Tk()
    Pmw.initialise(app.root)
    app.root.title('Some Title')

    widget = MHP_Compute(app)
#    widget = MHP_Display(app)
#    widget = MHP_doc(app)
    exitButton = Tkinter.Button(app.root, text = 'Exit', command = app.root.destroy)
    exitButton.pack()
    app.root.mainloop()

