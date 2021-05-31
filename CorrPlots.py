# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
#from scipy.optimize import curve_fit
from scipy import stats 
sns.set(color_codes=True)

dfAll = pd.read_csv (r'~/Desktop/USCorrAnalysis.csv')
#df = df.drop(36)
dfAll["alpha"] = np.log(2)/dfAll["alpha"]
df = dfAll
df = df.drop(58)
df = df.drop(57)
df = df.drop(56)

##############
### plot 1 ###
##############
h = sns.jointplot(data=df, x="DaysAfter", y="I0")
#plt.ylim([0,1.05*max(df["xi"])])
h.set_axis_labels('x', 'y', fontsize=16)

# or set labels via the axes objects
h.ax_joint.set_xlabel(r'First Reported Case (Days After Jan. 21)')
h.ax_joint.set_ylabel(r'$I_0$')
a,p = stats.spearmanr(df['DaysAfter'],df['I0'])
plt.show()
print('FRC, I0')
print([a,p])

##############
### plot 2 ###
##############
h = sns.jointplot(data=df, x="psi", y="CumCase")
#plt.ylim([0,1.05*max(df["xi"])])
h.set_axis_labels('x', 'y', fontsize=16)

# or set labels via the axes objects
h.ax_joint.set_xlabel(r'$\psi$')
h.ax_joint.set_ylabel(r'Cum. Case Est.')
a,p = stats.spearmanr(df['psi'],df['CumCase'])
print('psi,cum case')
print([a,p])
plt.show()
 ##############
### plot 3 ###
##############   
h = sns.jointplot(data=df, x="psi", y="PeakCase")
#plt.ylim([0,1.05*max(df["xi"])])
h.set_axis_labels('x', 'y', fontsize=16)

# or set labels via the axes objects
h.ax_joint.set_xlabel(r'$\psi$')
h.ax_joint.set_ylabel(r'Peak Case Est.')
a,p = stats.spearmanr(df['psi'],df['PeakCase'])
print('psi,peak case')
print([a,p])
plt.show()
##############
### plot 4 ###
##############
h = sns.jointplot(data=df, x="xi", y="RepCase")
#plt.ylim([0,1.05*max(df["xi"])])
h.set_axis_labels('x', 'y', fontsize=16)

# or set labels via the axes objects
h.ax_joint.set_xlabel(r'IFR')
h.ax_joint.set_ylabel(r'Cum. Reported Cases')
a,p = stats.spearmanr(df['xi'],df['RepCase'])
print('ifr, cum rep case')
print([a,p])
plt.show()
##############
### plot 5 ###
##############

h = sns.jointplot(data=df, x="xi", y="Ratio")
#plt.ylim([0,1.05*max(df["xi"])])
h.set_axis_labels('x', 'y', fontsize=16)

# or set labels via the axes objects
h.ax_joint.set_xlabel(r'IFR')
h.ax_joint.set_ylabel(r'Cum. Cases/Repot. Cases')
a,p = stats.spearmanr(df['xi'],df['Ratio'])
print('ifr, ratio')
print([a,p])
plt.show()
##############
### plot 6 ###
##############

h = sns.jointplot(data=df, x="psi", y="Deaths")
#plt.ylim([0,1.05*max(df["xi"])])
h.set_axis_labels('x', 'y', fontsize=16)

# or set labels via the axes objects
h.ax_joint.set_xlabel(r'$\psi$')
h.ax_joint.set_ylabel(r'Deaths')
a,p = stats.spearmanr(df['psi'],df['Deaths'])
print('psi,death')
print([a,p])
plt.show()
##############
### plot 7 ###
##############
dfND = df
dfND.drop(35)
h = sns.jointplot(data=df, x="psi", y="alpha")
#plt.ylim([0,1.05*max(df["xi"])])
h.set_axis_labels('x', 'y', fontsize=16)
# or set labels via the axes objects
h.ax_joint.set_xlabel(r'$\psi$')
h.ax_joint.set_ylabel(r'fatigue half-life (days)')
a,p = stats.spearmanr(df['psi'],df['alpha'])
print('psi,alpha')
print([a,p])
plt.show()
##############
### plot 8 ###
##############

df = df.assign(RelPeak = df["DaysAfter"]+df["TimePeak"])

#df.assign(RelPeak = df["DaysAfter"]+df["TimePeak"])
h = sns.jointplot(data=df, x="R0", y="TimePeak")
#plt.ylim([0,1.05*max(df["xi"])])
h.set_axis_labels('x', 'y', fontsize=16)
# or set labels via the axes objects
h.ax_joint.set_xlabel(r'$\mathcal{R}_0$')
h.ax_joint.set_ylabel(r'Time of Peak Cases (Days after First Reported Case)')
a,p = stats.spearmanr(df['R0'],df['TimePeak'])
print('Time Peak,R0')
print([a,p])
plt.show()

#df.assign(RelPeak = df["DaysAfter"]+df["TimePeak"])
h = sns.jointplot(data=df, x="R0", y="RelPeak")
#plt.ylim([0,1.05*max(df["xi"])])
h.set_axis_labels('x', 'y', fontsize=16)
# or set labels via the axes objects
h.ax_joint.set_xlabel(r'$\mathcal{R}_0$')
h.ax_joint.set_ylabel(r'Time of Peak Cases (Days after Jan 21)')
a,p = stats.spearmanr(df['R0'],df['RelPeak'])
print('Rel Peak,R0')
print([a,p])
plt.show()
##############
### plot 9 ###
##############
df = df.assign(CumAsc = 1/df['Ratio'])

#df.assign(RelPeak = df["DaysAfter"]+df["TimePeak"])
h = sns.jointplot(data=df, x="CumAsc", y="RepCase")
#plt.ylim([0,1.05*max(df["xi"])])
h.set_axis_labels('x', 'y', fontsize=16)

# or set labels via the axes objects
h.ax_joint.set_xlabel(r'End Cum. Asc. Ratio')
h.ax_joint.set_ylabel(r'Cum. Reported Cases')
a,p = stats.spearmanr(df['CumAsc'],df['RepCase'])
print('ascer,reported')
print([a,p])


###############
### plot 10 ###
###############
#df.assign(RelPeak = df["DaysAfter"]+df["TimePeak"])
h = sns.jointplot(data=df, x="CumCase", y="xi")
#plt.ylim([0,1.05*max(df["xi"])])
h.set_axis_labels('x', 'y', fontsize=16)

# or set labels via the axes objects
h.ax_joint.set_xlabel(r'Cum. Case Est.')
h.ax_joint.set_ylabel(r'IFR')
a,p = stats.spearmanr(df['CumCase'],df['xi'])
print('cum case,ifr')
print([a,p])
plt.show()
###############
### plot 11 ###
###############
#df.assign(RelPeak = df["DaysAfter"]+df["TimePeak"])
h = sns.jointplot(data=df, x="R0", y="psi")
#plt.ylim([0,1.05*max(df["xi"])])
h.set_axis_labels('x', 'y', fontsize=16)
h.ax_joint.set_xlabel(r'$\mathcal{R}_0$')
h.ax_joint.set_ylabel(r'$\psi$')
a,p = stats.spearmanr(df['psi'],df['R0'])
print('R0,psi')
print([a,p])
plt.show()

#print(df)
dfPeak = df
dfPeak.drop(49)
dfPeak.drop(36)
dfPeak.drop(41)
dfPeak.drop(11)
h = sns.jointplot(data=dfPeak, x="R0", y="PeakCase")
#plt.ylim([0,1.05*max(df["xi"])])
h.set_axis_labels('x', 'y', fontsize=16)
h.ax_joint.set_xlabel(r'$\mathcal{R}_0$')
h.ax_joint.set_ylabel(r'Peak Daily Cases')
a,p = stats.pearsonr(dfPeak['PeakCase'],dfPeak['R0'])
print('R0,peak case')
print([a,p])
plt.show()


###############
### plot 12 ###
###############




df2 = dfAll
df2["R0"] = df2["R0"].round(2)
df2["psi"] = df2["psi"].round(1)
#df2["alpha"] = np.log(2)/df2["alpha"]
df2["alpha"] = df2["alpha"].round(1)
df2["xi"] = df2["xi"].round(3)
df2["I0"] = df2["I0"].round(1)
df2["CumCase"] = df2["CumCase"].round(2)
df2["Ratio"] = df2["Ratio"].round(2)
df2["EndCases"] = df2["EndCases"].round(2)
df3 = df2[["State","R0","psi","alpha","xi","I0","CumCase","Ratio","EndCases"]]
df3

df2 = df2.assign(bg = 21/df['ag '])
df2["k"] = df2["k"].round(3)
df2["A"] = df2["A"].round(1)
df2["ag "] = df2["ag "].round(3)
df2["bg"] = df2["bg"].round(2)
df2["Deaths"] = df2["Deaths"].round(3)
df4 = df2[["State","k","A","ag ","bg","Deaths"]]
df4

print(df3.to_latex(index=False))

dfARE = pd.read_csv (r'~/Desktop/USAREs.csv')
dfARE = dfARE.round(2)
dfARE["R0"] = dfARE["R0"].round(2)
dfARE2 = dfARE[["State","R0","psi","alpha","xi","I0","CumCase","Ratio","EndCases"]]
print(dfARE2.to_latex(index=False))

df2 = dfAll 
xdata = df2["psi"]
ydata = df2["CumCase"]
zdata = df2["alpha"]
xdata = xdata.to_numpy()
ydata = ydata.to_numpy()
#xdata = (1-zdata)*xdata
zdata = zdata.to_numpy()
def func(x):
    return ((max(ydata)*(1+min(xdata))*1)/(1+x))

#popt, pcov = curve_fit(func, xdata, ydata,bounds=(0, 1))
m = min(xdata)
M = max(xdata)
xs = np.linspace(m,M)
plt.scatter(xdata,ydata)
plt.plot(xs,func(xs),'r')
plt.xlabel('$\psi$',fontsize=16)
plt.ylabel('True Cum. Case Est. (% pop)',fontsize=16)
plt.title('Inverse Proportionality Plot Cum. Cases',fontsize=16)
plt.legend(['$C_{NY}(1+\psi_{NY})(1+\psi)^{-1}$','Case Estimates for states'],fontsize=16)
plt.show()

df2 = dfAll
xdata = df2["psi"]
ydata = df2["PeakCase"]
zdata = df2["alpha"]
xdata = xdata.to_numpy()
ydata = ydata.to_numpy()
#xdata = (1-zdata)*xdata
zdata = zdata.to_numpy()
def func(x):
    return ((max(ydata)*(1+min(xdata)))/(1+x))
#xdata

m = min(xdata)
M = max(xdata)
xs = np.linspace(m,M)
plt.scatter(xdata,ydata)
plt.plot(xs,func(xs),'r')
plt.xlabel('$\psi$',fontsize=16)
plt.ylabel('True Peak Case Est. (% pop)',fontsize=16)
plt.title('Inverse Proportionality Plot Peak Cases',fontsize=16)
plt.legend(['$P_{NY}(1+\psi_{NY})(1+\psi)^{-1}$','Case Estimates for states'],fontsize=16)
plt.show()

rho = df.corr(method='spearman') 
fig, ax = plt.subplots(figsize=(14,14)) 
matrix = np.triu(rho)
A = sns.heatmap(rho,annot=True,mask=matrix)
A.set_yticklabels(A.get_yticklabels(), rotation=35)
plt.show()

