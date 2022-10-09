import warnings
import tkinter as tk
import matplotlib.pyplot as pl
from numpy import zeros, concatenate, ravel, diff, array, ones
import numpy as np
from scipy.interpolate import fitpack
from scipy.interpolate import dfitpack
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,  
NavigationToolbar2Tk) 

# Windows
Window = tk.Tk() 

# Judul & Ukuran
Window.title("Spline Interpolation")
Window.geometry('1300x700')
  
# Label
Array1 = [] #List 1
Array2 = [] #List 2
def Boxes():
    Data2 = tk.Frame(Data1)
    Data2.pack()
    
    xx=size.get()
    Array1.clear()
    Array2.clear()
    for i in range(xx):      
        listA = tk.Entry(Data2)
        listA.grid(row=2+i,column=0,sticky="w")
        Array1.append(listA)
        listB = tk.Entry(Data2) 
        listB.grid(row=2+i,column=1,sticky="w")
        Array2.append(listB)
    
    def Clear():
        Data2.destroy()
        
    ClearArray=tk.Button(Array,text="Clear",command=Clear)
    ClearArray.grid(row=0,column=3,sticky="w")

# Frame1
Array = tk.Frame(Window)
Array.pack()

text1=tk.Label(Array,text="Enter the Size of Array:", font="Arial 10 bold")
text1.grid(row=0,column=0,sticky="w")

size=tk.IntVar()

ArraySize=tk.Entry(Array,textvariable=size)
ArraySize.grid(row=0,column=1,sticky="w")

SizeofArray=tk.Button(Array,text="Submit",command=Boxes)
SizeofArray.grid(row=0,column=2,sticky="w")

# Frame2
Data1 = tk.Frame(Window)
Data1.pack()

# Frame3
Data2 = tk.Frame(Data1)
Data2.pack()

# Frame4
Number = tk.Frame(Window)
Number.pack()

text2=tk.Label(Number,text="Enter New Number to Interpolate:",font="Arial 10 bold")
text2.grid(row=1,column=0,sticky="w")

xs=tk.DoubleVar()

NNumber=tk.Entry(Number,textvariable=xs)
NNumber.grid(row=1,column=1,sticky="w")

NNumberB=tk.Button(Number,text="Submit",command=lambda:[Plot1(), Plot2(), Plot3(), Plot4()])
NNumberB.grid(row=1,column=2,sticky="w")

# Frame5
Plot = tk.Frame(Window)
Plot.pack()

Window.mainloop()