import numpy as np
import sympy as sp
from tkinter import *
from tkinter import messagebox
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from matplotlib import *

# Define a Vector2D class

class Vec2:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def __str__(self):
        return f"({self.x}, {self.y})"

    def __add__(self, v):
        return Vec2(self.x + v.x, self.y + v.y)

    def __radd__(self, v):
        return Vec2(self.x + v.x, self.y + v.y)

    def __sub__(self, v):
        return Vec2(self.x - v.x, self.y - v.y)

    def __rsub__(self, v):
        return Vec2(v.x- self.x , v.y - self.y)

    def __mul__(self, n):
        return Vec2(self.x * n, self.y * n)

    def __rmul__(self, n):
        return Vec2(self.x * n, self.y * n)

    def dot(self, v):
        return self.x*v.x + self.y*v.y

    def get_length(self):
        return np.sqrt(self.dot(self) )

# Define a Particle class. The particles are the bodies attracting each other

class Particle():
    n = 0
    def __init__(self,initial_pos,initial_vel, mass):

        self.i = Particle.n
        Particle.n += 1

        self.m = mass
        self.G = 1

        self.pos = Vec2(sp.symbols("x_"+str(self.i)),sp.symbols("y_"+str(self.i)))
        self.vel = Vec2(sp.symbols("vx_"+str(self.i)),sp.symbols("vy_"+str(self.i)))
        self.acc = Vec2(0,0)

        self.lamb_vel = Vec2(None,None)
        self.lamd_acc = Vec2(None,None)

        self.initial_pos = initial_pos
        self.initial_vel = initial_vel

        self.vf_vel = Vec2(0,0)
        self.vf_acc = Vec2(0,0)

        self.sol_pos = Vec2(None,None)
        self.sol_vel = Vec2(None,None)

    def calculate_acc(self,particles):
        for j in range(len(particles)):
            if self.i !=j:
                self.acc += (particles[j].pos - self.pos)*particles[j].m*self.G*(1/(((self.pos.x-particles[j].pos.x)**2 + (self.pos.y-particles[j].pos.y)**2)**(3/2)))


    def lambdify_vel(self,particles):
        self.lamb_vel.x = sp.lambdify(self.vel.x, self.vel.x)
        self.lamb_vel.y = sp.lambdify(self.vel.y, self.vel.y)
   

    def lambdify_acc(self,particles):
        
        var = []
        for j in range(len(particles)):           
            var.append(particles[j].pos.x)
            var.append(particles[j].pos.y)
               
        self.lamd_acc.x = sp.lambdify([var], self.acc.x)
        self.lamd_acc.y = sp.lambdify([var], self.acc.y)



#The initial conditions of the particles and their masses through a GUI
################################################################################################################################

par = []

def checkNum(number):
    for i in number:
        if i.isalpha():
            return True
    return False


def addParticle():
    if checkNum(initialPositionEntryX.get()):
        messagebox.showerror('Error', 'Wrong initial x position')
        return
    if checkNum(initialVelocityEntryX.get()):
        messagebox.showerror('Error', 'Wrong initial x velocity')
        return
    if checkNum(initialPositionEntryY.get()):
        messagebox.showerror('Error', 'Wrong initial y position')
        return
    if checkNum(initialVelocityEntryY.get()):
        messagebox.showerror('Error', 'Wrong initial y position')
        return
    if checkNum(massEntry.get()):
        messagebox.showerror('Error', 'Wrong mass')
        return
    xPosition = float(initialPositionEntryX.get())
    xVelocity = float(initialVelocityEntryX.get())
    yPosition = float(initialPositionEntryY.get())
    yVelocity = float(initialVelocityEntryY.get())
    mass = float(massEntry.get())
    par.append(Particle(initial_pos=Vec2(xPosition, yPosition), initial_vel=Vec2(xVelocity, yVelocity), mass=mass))
    return

def run():
    body.destroy()


body = Tk()
body.iconbitmap(default="Simulation.ico")
body.title("N-body simulation 2D")
body.geometry("410x200")
initialPositionLabelX = Label(body, text="Initial X Position")
initialPositionLabelX.grid(row=0, column=0, padx=5, pady=5)
initialVelocityLabelX = Label(body, text="Initial X Velocity")
initialVelocityLabelX.grid(row=0, column=1, padx=5, pady=5)
massLabel = Label(body, text="Mass")
massLabel.grid(row=0, column=2, padx=5, pady=5)
initialPositionEntryX = Entry(body)
initialPositionEntryX.grid(row=1, column=0, padx=5, pady=5)
initialVelocityEntryX = Entry(body)
initialVelocityEntryX.grid(row=1, column=1, padx=5, pady=5)
initialPositionLabelY = Label(body, text="Initial Y Position")
initialPositionLabelY.grid(row=2, column=0, padx=5, pady=5)
initialVelocityLabelY = Label(body, text="Initial Y Velocity")
initialVelocityLabelY.grid(row=2, column=1, padx=5, pady=5)
initialPositionEntryY = Entry(body)
initialPositionEntryY.grid(row=3, column=0, padx=5, pady=5)
initialVelocityEntryY = Entry(body)
initialVelocityEntryY.grid(row=3, column=1, padx=5, pady=5)
massEntry = Entry(body)
massEntry.grid(row=1, column=2, padx=5, pady=5)
addButton = Button(body, text="Add body", width=17, command=addParticle)
addButton.grid(row=3, column=2, padx=5, pady=15)
runButton = Button(body, text="Run", width=17, command=run)
runButton.grid(row=4, column=2, padx=5, pady=15)
body.mainloop()

# Simulation time and number of steps
t_end = 60.0
steps = 800


################################################################################################################################



n = len(par)


#create the functions to integrate
for i in range(n):
    par[i].calculate_acc(par)

for i in range(n):
    par[i].lambdify_vel(par)
    par[i].lambdify_acc(par)




import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def vectorfield(var, t):
    
    pos = var[0:2*n] 
    vel = var[2*n:4*n] 
    f = []
    
    for i in range(0,n):        
        par[i].vf_vel.x = par[i].lamb_vel.x(vel[2*i])
        par[i].vf_vel.y = par[i].lamb_vel.y(vel[2*i + 1])
        f.append(par[i].vf_vel.x)
        f.append(par[i].vf_vel.y)
        
    for i in range(0,n):        
        par[i].vf_acc.x = par[i].lamd_acc.x(pos)
        par[i].vf_acc.y = par[i].lamd_acc.y(pos)
        f.append(par[i].vf_acc.x)
        f.append(par[i].vf_acc.y)

    return f

from scipy.integrate import odeint

var = []
for i in range(len(par)):
    var.append(par[i].initial_pos.x)
    var.append(par[i].initial_pos.y)
    
for i in range(len(par)):
    var.append(par[i].initial_vel.x)
    var.append(par[i].initial_vel.y)


t = np.linspace(0,t_end,steps+1)

sol = odeint(vectorfield, var, t)
sol = np.transpose(sol)


for i in range(n):
    par[i].sol_pos.x = sol[2*i]
    par[i].sol_pos.y = sol[2*i+1]
    
for i in range(n):
    par[i].sol_vel.x = sol[2*n + 2*i]
    par[i].sol_vel.y = sol[2*n + 2*i+1]

Energy = 0 
for i in range(0,n):
    for j in range(i+1,n):
        Energy += (-1/(((par[i].sol_pos.x-par[j].sol_pos.x)**2 + (par[i].sol_pos.y-par[j].sol_pos.y)**2)**(1/2)))

for i in range(0,n):
    Energy += 0.5*(par[i].sol_vel.x*par[i].sol_vel.x + par[i].sol_vel.y*par[i].sol_vel.y)


plt.style.use('dark_background')
fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(1,1,1)

plt.subplots_adjust(bottom=0.2,left=0.15)

ax.axis('equal')
ax.axis([-1, 30, -1, 30])
ax.set_title('Energy =' + str(Energy[0]))
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)


circle = [None]*n
line  = [None]*n
for i in range(n):
    circle[i] = plt.Circle((par[i].sol_pos.x[0], par[i].sol_pos.y[0]), 0.08, ec="w", lw=2.5, zorder=20)
    ax.add_patch(circle[i])
    line[i] = ax.plot(par[i].sol_pos.x[:0],par[i].sol_pos.y[:0])[0]
    
from matplotlib.widgets import Slider

slider_ax = plt.axes([0.1, 0.05, 0.8, 0.05])
slider = Slider(slider_ax,'t',0,t_end,valinit=0,color = '#5c05ff')

def update(time):
    i = int(np.rint(time*steps/t_end))
    
    ax.set_title('Energy =' + str(Energy[i]))
    for j in range(n):
        circle[j].center = par[j].sol_pos.x[i], par[j].sol_pos.y[i]
        line[j].set_xdata(par[j].sol_pos.x[:i+1])
        line[j].set_ydata(par[j].sol_pos.y[:i+1])
        
slider.on_changed(update)
plt.show()