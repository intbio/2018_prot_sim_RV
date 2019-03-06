#Проект. Построение графика зависимости концентрации белка от времени.
#Используемые библиотеки
from bokeh.plotting import Figure, output_file, show
import math
import scipy
import numpy as np
from bokeh.layouts import widgetbox, row, column
from bokeh.models import CustomJS, ColumnDataSource, Slider, Button, PreText, TextInput
output_file("Project_(bokeh_version).html")

#Переменные
CmrnaO = (2*math.pow(10,16)/(6.02*math.pow(10,23)))*math.pow(10,6) #В мкмоль/л - начальная концентрация мРНК в клетке E.Coli 
CprotO = 0 #В мкмоль/л - начальная концентрация белка в клетке E.Coli
CRNAP =  (22*math.pow(10,17)/(6.02*math.pow(10,23)))*math.pow(10,6) #В мкмоль/л - концентрация РНК-полимеразы II в клетке E.Coli
Cribosomes = (27*math.pow(10,18)/(6.02*math.pow(10,23)))*math.pow(10,6) #В мкмоль/л - концентрация рибосом в клетке E.Coli
koc = 0.03 #В 1/c - Open complex formation rate
Np0 = 3 #Число копий плазмиды в клетке E.Coli
Kp = 10 #В µM - константа диссоциации для связывания РНК-полимеразы II с ДНК
kdm = 0.00278 #В 1/c - фактор деградации мРНК
kdp0 = 0.00028 #В 1/c - фактор деградации белка
kt0 = 0.0667 #В 1/c - фактор инициации трансляции
dt1=[dt*0.04 for dt in range(0, 250)]
y = [0]*250

#Создание осей графика и ползунков
source = ColumnDataSource(data=dict(dt1=dt1, y=y))
fig = Figure(plot_width=500, plot_height=500, title="Зависимость концентрации белка в клетке от времени",x_axis_label='Время, с', y_axis_label='Концентрация белка,µМоль/л')
#Построение графика
for i in range(0,250):
    if (i == 0):
        Cmrna = CmrnaO
        Cprot = CprotO
    Cmrna = Cmrna*(1 - kdm*dt1[i]) + koc*Np0*(CRNAP/(CRNAP + Kp))*dt1[i]
    Cprot = Cprot*(1-kdp0*dt1[i]) + kt0*Cmrna*Cribosomes*dt1[i]
    y[i] = Cprot
    
fig.circle('dt1','y', source=source, size = 3, color = 'red')

#Перестройка графика при изменении параметров

#Общий callback
callback = CustomJS(args=dict(source=source), code="""
    var Cml = Cmrna.value.length
    var Cpl = Cprot.value.length
    var CRl = CRNAPin.value.length
    var Cril = Cribosomesin.value.length
    var kol = kocin.value.length
    var Kpl = Kpin.value.length
    var kdl = kdmin.value.length != 0 
    if(Cml != 0){CmrnaO = parseFloat(Cmrna.value)}else{CmrnaO = (2*Math.pow(10,16)/(6.02*Math.pow(10,23)))*Math.pow(10,6)}
    if(Cpl != 0){CprotO = parseFloat(Cprot.value)}else{CprotO = 0}
    if(CRl != 0){CRNAP = parseFloat(CRNAPin.value)}else{CRNAP =  (22*Math.pow(10,17)/(6.02*Math.pow(10,23)))*Math.pow(10,6)}
    if(Cril != 0){Cribosomes = parseFloat(Cribosomesin.value)}else{Cribosomes = (27*Math.pow(10,18)/(6.02*Math.pow(10,23)))*Math.pow(10,6)} 
    if(kol != 0){koc = parseFloat(kocin.value)}else{koc = 0.03}
    Np0 = Np.value
    if(Kpl != 0){Kp = parseFloat(Kpin.value)}else{Kp = 10}
    kt0 = kt.value
    if(kdl != 0){kdm = parseFloat(kdmin.value)}else{kdm = 0.00278} 
    var data = source.data
    var dt1 = data['dt1']
    var y = data['y']
    kdp0 = kdp.value
    for (var i = 0; i < 250; i++)
    {
        if (i == 0)
        {
            Cmrna = CmrnaO
            Cprot = CprotO
        }
        Cmrna = Cmrna*(1 - kdm*dt1[i]) + koc*Np0*(CRNAP/(CRNAP + Kp))*dt1[i]
        Cprot = Cprot*(1-kdp0*dt1[i]) + kt0*Cmrna*Cribosomes*dt1[i]
        y[i] = Cprot
    }
    Cmrna.value = ((2*Math.pow(10,16)/(6.02*Math.pow(10,23)))*Math.pow(10,6)).toString()
    Cprot.value = (0).toString()
    CRNAPin.value = ((22*Math.pow(10,17)/(6.02*Math.pow(10,23)))*Math.pow(10,6)).toString()
    Cribosomesin = ((27*Math.pow(10,18)/(6.02*Math.pow(10,23)))*Math.pow(10,6)).toString()
    kocin.value = (0.03).toString()
    Kpin.value = (10).toString()
    kdmin.value = (0.00278).toString()
    source.change.emit();
""")

#Слайдеры
Npsl = Slider(title = 'Np', start=1, end=10, value=Np0, step=1, callback=callback)
callback.args["Np"] = Npsl
ktsl = Slider(title = 'kt', start=0.04, end=0.08, value=kt0, step = 0.01, callback=callback)
callback.args["kt"] = ktsl
kdpsl = Slider(title = 'kdp', start=0, end=0.04, value=kdp0, step =.001, callback=callback)
callback.args["kdp"] = kdpsl

#Textbox'ы
ICmrnaO = TextInput(callback=callback, value=str(CmrnaO))
callback.args["Cmrna"] = ICmrnaO
TCmrnaO = PreText(text = "CmrnaO")
ICprotO = TextInput(callback=callback, value=str(CprotO))
callback.args["Cprot"] = ICprotO
TCprotO = PreText(text = "CprotO")
ICRNAP = TextInput(callback=callback, value=str(CRNAP))
callback.args["CRNAPin"] = ICRNAP
TCRNAP = PreText(text = "CRNAP")
ICribosomes = TextInput(callback=callback, value=str(Cribosomes))
callback.args["Cribosomesin"] = ICribosomes
TCribosomes = PreText(text = "Cribosomes")
IKp = TextInput(callback=callback, value=str(Kp))
callback.args["Kpin"] = IKp
TKp = PreText(text = "Kp")
Ikoc = TextInput(callback=callback, value=str(koc))
callback.args["kocin"] = Ikoc
Tkoc = PreText(text = "koc")
Ikdm = TextInput(callback=callback, value=str(kdm))
callback.args["kdmin"] = Ikdm
Tkdm = PreText(text = "kdm")

#Кнопка Reset
reset = CustomJS(args=dict(source=source), code="""
    CmrnaO = (2*Math.pow(10,16)/(6.02*Math.pow(10,23)))*Math.pow(10,6) 
    CprotO = 0
    CRNAP = (22*Math.pow(10,17)/(6.02*Math.pow(10,23)))*Math.pow(10,6)
    Cribosomes = (27*Math.pow(10,18)/(6.02*Math.pow(10,23)))*Math.pow(10,6) 
    koc = 0.03
    Np0 = 3
    Kp = 10
    kt0 = 0.0667
    kdm = 0.00278 
    var data = source.data
    var dt1 = data['dt1']
    var y = data['y']
    var kdp0 =  0.00028 
    for (var i = 0; i < 250; i++)
    {
        if (i == 0)
        {
            Cmrna = CmrnaO
            Cprot = CprotO
        }
        Cmrna = Cmrna*(1 - kdm*dt1[i]) + koc*Np0*(CRNAP/(CRNAP + Kp))*dt1[i]
        Cprot = Cprot*(1-kdp0*dt1[i]) + kt0*Cmrna*Cribosomes*dt1[i]
        y[i] = Cprot
    }
    source.change.emit();
""")
resetbt = Button(callback=reset,label="Reset")

show(column(row(fig,widgetbox(ICmrnaO,TCmrnaO,ICprotO,TCprotO,ICRNAP,TCRNAP,ICribosomes,TCribosomes,IKp,TKp,Ikoc,Tkoc,Ikdm,Tkdm)),widgetbox(Npsl,ktsl,kdpsl,resetbt)))


