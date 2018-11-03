#Проект. Построение графика зависимости концентрации белка от времени.

#Используемые библиотеки
from pylab import *
import matplotlib.pyplot as plt
import math
import scipy

#Переменные - ошибка в них, из-за этого в клетке слишком большая концентрация белка
CellVolume = 7*math.pow(10,-15) #В литрах - объём клетки E.Coli
GeneNum = 4288 #Число генов у E.Coli
CmrnaO = 5.48*math.pow(10,-6)/(CellVolume*GeneNum) #В мкмоль/л - начальная концентрация мРНК в клетке E.Coli 
TransNum = 60 #Число трансляций с одной мРНК в клетке E.Coli
#RNAPNum = 3000 #Число РНК-полимераз II в клетке E.Coli
CprotO = CmrnaO*TransNum #В мкмоль/л - начальная концентрация белка в клетке E.Coli
CRNAP = 5.48*math.pow(10,-6)/(CellVolume*GeneNum) #В мкмоль/л - концентрация РНК-полимеразы II в клетке E.Coli
Cribosomes = (420*math.pow(10,-15)/(27*math.pow(10,5)*1.66*math.pow(10,-24)))*math.pow(10,6)/CellVolume #В мкмоль/л - концентрация рибосом в клетке E.Coli - наибольшая ошибка, вероятно, здесь
koc = 0.03 #В 1/c - Open complex formation rate
Np = 1 #Число копий плазмиды в клетке E.Coli
Kp = 10 #В µM - константа диссоциации для связывания РНК-полимеразы II с ДНК
kdm = 0.00278 #В 1/c - фактор деградации мРНК
kdp = 0.00028 #В 1/c - фактор деградации белка
kt = 0.0667 #В 1/c - фактор инициации трансляции

#Построение графика
fig = plt.figure('MyProj',[8, 6])
upbd = 10 #Верхняя граница, нижняя - 0
k = 0
for dt in np.arange(0,upbd,0.02):
    if (k == 0):
        Cmrna = CmrnaO
        Cprot = CprotO
    Cmrna = Cmrna*(1 - kdm*dt) + koc*Np*(CRNAP/(CRNAP + Kp))*dt
    Cprot = Cprot*(1-kdp*dt) + kt*Cmrna*Cribosomes*dt
    plt.scatter(dt,Cprot, c = 'red')
    k += 1

grid1 = plt.grid(True)
plt.ylabel('Концентрация белка,µМоль/л')
plt.xlabel('Время, с')
plt.title('Зависимость концентрации белка в клетке от времени')
plt.show()
