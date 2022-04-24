import pandas as pd
import math
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap

import sys
from fractions import Fraction
import matplotlib.ticker as mticker
import tkinter as tk
from tkinter import Canvas



import urllib.request

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

values = pd.read_csv('dipole_moments.csv', sep = ';', encoding='latin-1')
data = pd.DataFrame(values)


def magnetic_moment(J,nucleon,sign,q):
    J = float(J)
    if nucleon == 'p':
        if sign == '+':
            return J - 0.5 + 2.793*q
        if sign == '-':
            return J*(J + 3/2 - 2.793*q)/(J+1)
    if nucleon == 'n':
        if sign == '+':
            return -1.913*q
        if sign == '-':
            return J*1.913*q/(J+1)


Mass_Number = data['A']
Atomic_Number = data['Z']

#Here you select the range of isotopes that you want to plot

Isotopes = data.loc[(data["A"] > 0) &  (data["Z"] >= 0)  & (data["Duplicate"] == 1) ] #& (data["Ex [keV]"] == 0)

Isotopes['J'] = Isotopes['J'].astype(float)
Isotopes['m'] = Isotopes['m'].astype(float)
Isotopes['Ex [keV]'] = Isotopes['Ex [keV]'].astype(float)
Isotopes['error'] = Isotopes['error'].astype(float)
Isotopes['Z'] = Isotopes['Z'].astype(int)
Isotopes['A'] = Isotopes['A'].astype(int)
Isotopes['N'] = Isotopes['N'].astype(int)


Isotopes_odd_proton = Isotopes.loc[(Isotopes['Z'] % 2 != 0) & (Isotopes['A'] % 2 != 0) ] # Selects values for the odd proton case
Isotopes_odd_neutron = Isotopes.loc[(Isotopes['A'] % 2 != 0) & (Isotopes['Z'] % 2 == 0)] # Selects values for the odd neutron case


Schmidt_proton_plus = [];
Schmidt_proton_minus = [];
Schmidt_proton_J = [];

Schmidt_proton_plus_shorter = [];
Schmidt_proton_minus_shorter = [];
Schmidt_proton_J_shorter = [];


Schmidt_neutron_plus = [];
Schmidt_neutron_minus = [];
Schmidt_neutron_J = [];

Schmidt_neutron_plus_shorter = [];
Schmidt_neutron_minus_shorter = [];
Schmidt_neutron_J_shorter = [];


for i in sorted(Isotopes_odd_neutron['J']):
    Schmidt_neutron_plus.append( magnetic_moment(i, 'n','+',1) )

for i in sorted(Isotopes_odd_neutron['J']):
    Schmidt_neutron_minus.append( magnetic_moment(i, 'n','-',1) )

for i in sorted(Isotopes_odd_neutron['J']):
    Schmidt_neutron_plus_shorter.append( magnetic_moment(i, 'n','+',0.7) )

for i in sorted(Isotopes_odd_neutron['J']):
    Schmidt_neutron_minus_shorter.append( magnetic_moment(i, 'n','-',0.7) )

for i in sorted(Isotopes_odd_neutron['J']):
    Schmidt_neutron_J.append(i)


for i in sorted(Isotopes_odd_proton['J']):
    Schmidt_proton_plus.append( magnetic_moment(i, 'p','+',1) )

for i in sorted(Isotopes_odd_proton['J']):
    Schmidt_proton_minus.append( magnetic_moment(i, 'p','-',1) )

for i in sorted(Isotopes_odd_proton['J']):
    Schmidt_proton_plus_shorter.append( magnetic_moment(i, 'p','+',0.7) )

for i in sorted(Isotopes_odd_proton['J']):
    Schmidt_proton_minus_shorter.append( magnetic_moment(i, 'p','-',0.7) )

for i in sorted(Isotopes_odd_proton['J']):
    Schmidt_proton_J.append(i)

X_tick_p = [str(Fraction(item).limit_denominator()) for item in Schmidt_proton_J]
X_tick_n = [str(Fraction(item).limit_denominator()) for item in Schmidt_neutron_J]


step = 1./2


def fractions(x,pos):
    if np.isclose((x/step)%(1./step),0.):
        # x is an integer, so just return that
        return '{:.0f}'.format(x)
    else:
        # this returns a latex formatted fraction
        return '$\\frac{{{:2.0f}}}{{{:2.0f}}}$'.format(x/step,1./step)
        # if you don't want to use latex, you could use this commented
        # line, which formats the fraction as "1/13"
        ### return '{:2.0f}/{:2.0f}'.format(x/step,1./step)


Isotopes_odd_proton.plot(x = 'J', y = 'm', yerr = 'error', s = 100, kind = 'scatter', color = 'mediumpurple', alpha = 0.3, figsize=(15, 12))

ax = plt.gca()

ax.xaxis.set_major_locator(mticker.MultipleLocator(step))
ax.xaxis.set_major_formatter(mticker.FuncFormatter(fractions))

plt.xlabel('J')
plt.ylabel('Magnetic Moment ($\mu_N$)')
plt.plot(Schmidt_proton_J, Schmidt_proton_plus, color = 'Firebrick', label = 'J = l + 1/2', linewidth=3.0)
plt.plot(Schmidt_proton_J, Schmidt_proton_minus, color = 'slategray', label = 'J = l - 1/2', linewidth=3.0)
plt.plot(Schmidt_proton_J, Schmidt_proton_plus_shorter, color = 'Firebrick', label = 'J = l + 1/2, 70 %', linewidth=2.0, linestyle='--')
plt.plot(Schmidt_proton_J, Schmidt_proton_minus_shorter, color = 'slategray', label = 'J = l - 1/2, 70 %', linewidth=2.0, linestyle='--')
plt.legend()
plt.xticks(np.arange(0.5, max(Schmidt_proton_J)+1, step=1), fontsize=22)
plt.title('Odd proton')
plt.show()



Isotopes_odd_neutron.plot(x = 'J', y = 'm', yerr = 'error', s = 100, kind = 'scatter', color = 'mediumpurple', alpha = 0.3, figsize=(15, 12)) #Changing J to A will plot them as a function of the mass number

ax = plt.gca()

ax.xaxis.set_major_locator(mticker.MultipleLocator(step))
ax.xaxis.set_major_formatter(mticker.FuncFormatter(fractions))

plt.xlabel('J')
plt.ylabel('Magnetic Moment ($\mu_N$)')
plt.plot(Schmidt_neutron_J, Schmidt_neutron_plus, color = 'Firebrick', label = 'J = l + 1/2', linewidth=3.0)
plt.plot(Schmidt_neutron_J, Schmidt_neutron_minus, color = 'slategray', label = 'J = l - 1/2', linewidth=3.0)
plt.plot(Schmidt_neutron_J, Schmidt_neutron_plus_shorter, color = 'Firebrick', label = 'J = l + 1/2, 70 %', linewidth=2.0, linestyle='--')
plt.plot(Schmidt_neutron_J, Schmidt_neutron_minus_shorter, color = 'slategray', label = 'J = l - 1/2, 70 %', linewidth=2.0, linestyle='--')
plt.legend()
plt.xticks(np.arange(0.5, max(Schmidt_neutron_J)+1, step=1), fontsize=22)
plt.title('Odd neutron')
plt.show()



#########################################################################################################

def lc_read_csv(url):
    req = urllib.request.Request(url)
    req.add_header('User-Agent', 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0')
    return pd.read_csv(urllib.request.urlopen(req))

livechart = "https://nds.iaea.org/relnsd/v0/data?"

df = lc_read_csv(livechart+"fields=decay_rads&nuclides=all&rad_types=g") # all returns only g.s.
df = df.replace(' ', '123456')

pd.set_option('display.max_rows', 100)
print(Isotopes_odd_proton)
print(Isotopes_odd_neutron)


Livechart = pd.read_csv('all_data.csv', sep = ',', encoding='latin-1')
Livechart = pd.DataFrame(Livechart)

Livechart['halflife(Seconds)'] = Livechart['halflife(Seconds)'].astype(str)
Livechart['n'] = Livechart['n'].astype(int)
Livechart['z'] = Livechart['z'].astype(int)

Livechart_v2 = pd.read_csv('all_data_v2.csv', sep = ',', encoding='latin-1')
Livechart_v2 = pd.DataFrame(Livechart_v2)

Livechart_v2['halflife(Seconds)'] = Livechart_v2['halflife(Seconds)'].astype(str)
Livechart_v2['n'] = Livechart_v2['n'].astype(int)
Livechart_v2['z'] = Livechart_v2['z'].astype(int)

app = tk.Tk()
app.title("Chart of Nuclides")

canvas = Canvas(app, bg = 'lavender blush')
canvas.pack()

canvas.config(width=2000, height=1800)

def rectangle(x1,y1,x2,y2,color):
    return canvas.create_rectangle(x1, y1, x2, y2, fill=color, outline = 'white')

x_min = 20
y_min = 20

count = 0

Atomic_number = 118
Neutron_number = 0

alpha = 0.8


count_v2 = 0

for i in np.arange(0,Atomic_number+1): # Don't touch 0
    y1 = y_min + (i-1)*10*alpha
    y2 = y_min + (i)*10*alpha

    for j in np.arange(1,182): # Don't touch 1
        x1 = x_min + (j-1)*10*alpha
        x2 = x_min + (j)*10*alpha

        for l,k,p in zip(df['z'],df['n'],df['magnetic_dipole']):

            if (( Atomic_number == int(l) and Neutron_number == int(k) ) and p == '123456' ):
                rectangle(x1,y1,x2,y2,'grey75')
                count += 1
                break


            if (( Atomic_number == int(l) and Neutron_number == int(k) ) and p != '123456' ):
                rectangle(x1,y1,x2,y2,'navy')
                count_v2 += 1
                break


        for l,k,p in zip(Isotopes['Z'],Isotopes['A'],Isotopes['Ex [keV]']):

            if (Atomic_number == int(l) and Neutron_number == int(k - l) and str(p) == '0.0'):
                rectangle(x1,y1,x2,y2,'Firebrick')
                break

        Neutron_number += 1

    Atomic_number -= 1
    Neutron_number = 0



canvas.create_rectangle(80*alpha + x_min,950 + y_min,90*alpha + x_min,800 + y_min, outline = 'black') # N = 8
canvas.create_rectangle(200*alpha + x_min,890 + y_min,210*alpha + x_min,690 + y_min, outline = 'black') # N = 20
canvas.create_rectangle(280*alpha + x_min,870 + y_min,290*alpha + x_min,650 + y_min, outline = 'black') # N = 28
canvas.create_rectangle(500*alpha + x_min,750 + y_min,510*alpha + x_min,500 + y_min, outline = 'black') # N = 50
canvas.create_rectangle(820*alpha + x_min,600 + y_min,830*alpha + x_min,320 + y_min, outline = 'black') # N = 82
canvas.create_rectangle(1280*alpha + x_min,350 + y_min,1290*alpha + x_min,150 + y_min, outline = 'black') # N = 126

canvas.create_rectangle(10 + x_min,1090*alpha + y_min,170 + x_min,1100*alpha + y_min,outline = 'black') # Z = 8
canvas.create_rectangle(100 + x_min,970*alpha + y_min,340 + x_min,980*alpha + y_min,outline = 'black') # Z = 20
canvas.create_rectangle(140 + x_min,890*alpha + y_min,450 + x_min,900*alpha + y_min,outline = 'black') # Z = 28
canvas.create_rectangle(370 + x_min,670*alpha + y_min,750 + x_min,680*alpha + y_min,outline = 'black') # Z = 50
canvas.create_rectangle(750 + x_min,350*alpha + y_min,1150 + x_min,360*alpha + y_min,outline = 'black') # Z = 82

app.mainloop()
