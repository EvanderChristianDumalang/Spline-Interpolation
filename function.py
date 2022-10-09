# Plot1
def Plot1():
    xx=size.get()
    xnew = float(xs.get())
    x = []
    y = []
    try:
        for i in range (xx):
            x.append(float(Array1[i].get()))
            y.append(float(Array2[i].get()))
    except:
        tk.messagebox.showwarning("Warning",f"There's Something Wrong, Element of Array Empty.")

    start = 0
    i = 0
    while x[start] < xnew:
        start += 1
        i += 1

    spl1 = InterpolatedUnivariateSpline(x, y, k=1)
    Spline1=tk.Label(Number, text ="Spline Orde 1: " + str(spl1(xnew)))
    Spline1.grid(row=2,column=0,sticky="w")

    fig = Figure(figsize = (3, 3), 
                 dpi = 100) 

    pl1 = fig.add_subplot(111) 

    x.insert(i, xnew)
    y.insert(i, spl1(xnew))

    pl1.plot(x,y,'bo-')
    pl1.set_title('Spline Orde 1')

    canvas1 = FigureCanvasTkAgg(fig, master = Plot)   
    canvas1.get_tk_widget().grid(row=6,column=0,sticky="w")

# Plot2
def Plot2():
    xx=size.get()
    xnew = float(xs.get())
    x = []
    y = []
    try:
        for i in range (xx):
            x.append(float(Array1[i].get()))
            y.append(float(Array2[i].get()))
    except:
        tk.messagebox.showwarning("Warning",f"There's Something Wrong, Element of Array Empty.")

    start = 0
    i = 0
    while x[start] < xnew:
        start += 1
        i += 1

    spl2 = InterpolatedUnivariateSpline(x, y, k=2)
    Spline2=tk.Label(Number, text ="Spline Orde 2: " + str(spl2(xnew)))
    Spline2.grid(row=3,column=0,sticky="w")

    fig = Figure(figsize = (3, 3), 
                 dpi = 100) 

    pl2 = fig.add_subplot(111) 

    x.insert(i, xnew)
    y.insert(i, spl2(xnew))

    pl2.plot(x,y,'ro-')
    pl2.set_title('Spline Orde 2')

    canvas2 = FigureCanvasTkAgg(fig, master = Plot)   
    canvas2.get_tk_widget().grid(row=6,column=1,sticky="w")

# Plot3
def Plot3():
    xx=size.get()
    xnew = float(xs.get())
    x = []
    y = []
    try:
        for i in range (xx):
            x.append(float(Array1[i].get()))
            y.append(float(Array2[i].get()))
    except:
        tk.messagebox.showwarning("Warning",f"There's Something Wrong, Element of Array Empty.")

    start = 0
    i = 0
    while x[start] < xnew:
        start += 1
        i += 1

    spl3 = InterpolatedUnivariateSpline(x, y, k=3)
    Spline3=tk.Label(Number, text ="Spline Orde 3: " + str(spl3(xnew)))
    Spline3.grid(row=4,column=0,sticky="w")

    fig = Figure(figsize = (3, 3), 
                 dpi = 100) 

    pl3 = fig.add_subplot(111) 

    x.insert(i, xnew)
    y.insert(i, spl3(xnew))

    pl3.plot(x,y,'go--')
    pl3.set_title('Spline Orde 3')

    canvas3 = FigureCanvasTkAgg(fig, master = Plot)   
    canvas3.get_tk_widget().grid(row=6,column=2,sticky="w")

# Plot1&Plot2&Plot3
def Plot4():
    a = []
    b = []
    x = []
    y = []
    xx=size.get()
    xnew = float(xs.get())
    try:
        for i in range (xx):
            x.append(float(Array1[i].get()))
            y.append(float(Array2[i].get()))
            a.append(float(Array2[i].get()))
            b.append(float(Array2[i].get()))
    except:
        tk.messagebox.showwarning("Warning",f"There's Something Wrong, Element of Array Empty.")

    spl1 = InterpolatedUnivariateSpline(x, y, k=1)
    spl2 = InterpolatedUnivariateSpline(x, y, k=2)
    spl3 = InterpolatedUnivariateSpline(x, y, k=3)

    start = 0
    i = 0
    while x[start] < xnew:
        start += 1
        i += 1

    fig = Figure(figsize = (3, 3), 
                 dpi = 100) 

    pl4 = fig.add_subplot(111) 

    x.insert(i, xnew)
    y.insert(i, spl1(xnew))
    a.insert(i, spl2(xnew))
    b.insert(i, spl3(xnew))

    pl4.set_title('Differences')
    pl4.plot(x,y,'bo-')
    pl4.plot(x,a,'ro-')
    pl4.plot(x,b,'go--')

    canvas4 = FigureCanvasTkAgg(fig, master = Plot)   
    canvas4.get_tk_widget().grid(row=6,column=3,sticky="w")