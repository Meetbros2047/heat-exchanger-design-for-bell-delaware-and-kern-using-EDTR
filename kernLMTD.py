# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 00:30:15 2024

@author: mitansh waghmare
"""
import tkinter as tk
from tkinter import Label, Entry, Button
from math import exp
from math import pi
import pandas as pd
import math
import time
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
start_time = time.time()
# given quantities
def solve_equations():
    mc = float(entry_mc.get())
    mh = float(entry_mh.get())
    tc1 = float(entry_tc1.get())
    th1 = float(entry_th1.get())
    th2 = float(entry_th2.get())
    tc2 = float(entry_tc2.get())
    l = float(entry_l.get())
    kw = float(entry_kw.get())
    hot = entry_hot.get()
    cold = entry_cold.get()
    Rd = float(entry_Rd.get())
    p_allowed = float(entry_p_allowed.get())
    p1_allowed = float(entry_p1_allowed.get())
    sgh = float(entry_sgh.get())
    sgc = float(entry_sgc.get())
    

    tcavg = (tc1 + tc2)/2
    thavg = (th1+th2)/2
# defining useful function that can be used in both the type of heat exchanger
    global cpc,cph,uh,uc,uwh,uwc,kh,kc,b
    def input_func():
        global cpc
        cpc = float(entry_cpc.get())
        global cph
        cph = float(entry_cph.get())
        global uh
        uh = float(entry_uh.get())
        global uc
        uc = float(entry_uc.get())
        global uwh
        uwh = float(entry_uwh.get())
        global uwc
        uwc = float(entry_uwc.get())
        global kh
        kh = float(entry_kh.get())
        global kc
        kc = float(entry_kc.get())
    data_input = data1.get()
    heatexchanger = heat1.get()
    if data_input=='input':
        input_func()         
        
    elif data_input=='database':
        def get_specific_heat_capacity(T,liquid_name):
            excel_file = 'kerndataforcode.xlsx'
            specific_heat_sheet = 'Sheet3'
            specific_heat_df = pd.read_excel(excel_file, sheet_name=specific_heat_sheet)
            specific_heat_constants = specific_heat_df.loc[specific_heat_df['Name'] == liquid_name, 'C1':'C5'].values.flatten().astype(float)
            specific_heat_equation = specific_heat_df.loc[specific_heat_df['Name'] == liquid_name, 'Eqn'].values
            mw = specific_heat_df.loc[specific_heat_df['Name'] == liquid_name, 'Mol. wt.'].values
            mw = mw[0]
            c1s = (specific_heat_constants[0])
            c2s = (specific_heat_constants[1])
            c3s = (specific_heat_constants[2])
            c4s = (specific_heat_constants[3])
            c5s = (specific_heat_constants[4])
            print(specific_heat_equation)
            if specific_heat_equation == 100.0:
                # Use Equation1 for specific heat capacity
                cpo = c1s + c2s*T + c3s*(T**2) + c4s*(T**3) + c5s*(T**4)
            elif specific_heat_equation == 114.0:
                # Use Equation2 for specific heat capacity
                Tc = specific_heat_df.loc[specific_heat_df['Name'] == liquid_name, 'Tc'].values.flatten().astype(float)
                Tr  = T/Tc
                tou = 1-Tr
                cpo = (c1s**2)/tou + c2s -2*c1s*c3s*tou-c1s*c4s*(tou**2)-(c3s**2)*(tou**3)/3-c3s*c4s*(tou**4)/2-(c4s**2)*(tou**5)/2
            return cpo/mw
        def get_viscosity_value(T,liquid_name):
            excel_file = 'kerndataforcode.xlsx'
            viscosity_sheet = 'Sheet1'
            viscosity_df = pd.read_excel(excel_file, sheet_name=viscosity_sheet)
            viscosity_constants = viscosity_df.loc[viscosity_df['Name'] == liquid_name, 'C1':'C5'].values.flatten().astype(float)
            viscosity_equation = viscosity_df.loc[viscosity_df['Name'] == liquid_name, 'Eqn'].values
            c1v = (viscosity_constants[0])
            c2v = (viscosity_constants[1])
            c3v = (viscosity_constants[2])
            c4v = (viscosity_constants[3])
            c5v = (viscosity_constants[4])
            if viscosity_equation == 101:
                # Use Equation1 for viscosity
                mu = math.e**(c1v+(c2v/T)+c3v*math.log(T)+c4v*(T**c5v))
            else:
                # Use Equation2 for viscosity
                mu = c1v+c2v*T+c3v*(T**2)+c4v*(T**3)+c5v**(T**4) 
            return mu
        def get_conductivity_value(T,liquid_name):
            excel_file = 'kerndataforcode.xlsx'
            thermal_conductivity_sheet = 'Sheet2'
            thermal_conductivity_df = pd.read_excel(excel_file, sheet_name=thermal_conductivity_sheet)
            thermal_conductivity_constants = thermal_conductivity_df.loc[thermal_conductivity_df['Name'] == liquid_name, 'C1':'C5'].values.flatten().astype(float)
            c1c = (thermal_conductivity_constants[0])
            c2c = (thermal_conductivity_constants[1])
            c3c = (thermal_conductivity_constants[2])
            c4c = (thermal_conductivity_constants[3])
            c5c = (thermal_conductivity_constants[4])
            k = c1c+c2c*T+c3c*(T**2)+c4c*(T**3)+c5c*(T**4)
            return k

        # the specific heat capacity at temperature tcavg and thavg
        cpc = get_specific_heat_capacity(tcavg,cold)
        cph = get_specific_heat_capacity(thavg,hot)
    q = mc*cpc*(tc2-tc1)
    tam = ((th1+th2)-(tc1+tc2))/2

    def doublepipeheat_exchanger():
        global doubletype
        doubletype = entry_doubletype.get()
        global D 
        D = float(entry_D.get())
        global D1
        D1 = float(entry_D1.get())
        global D2
        D2 = float(entry_D2.get())
        global up, uwp, ua, uwa,kp,ka,tw,pipe,tpavg,annulus,ap,aa,mp,ma,cpp,cpa,tpavg,taavg,sgp,sga,pipe,annulus,gp,ga,De
        ap = (pi*(D**2))/4
        aa = (pi*((D2**2)-(D1**2)))/4
        if (ap>aa and mh>mc) or (ap<aa and mh<mc):
            mp = mh
            ma = mc
            cpp = cph
            cpa = cpc
            tpavg = thavg
            taavg = tcavg
            sgp = sgh
            sga = sgc
            pipe = hot
            annulus = cold
        else:
            mp = mc
            ma = mh
            cpp = cpc
            cpa = cph
            tpavg = tcavg
            taavg = tcavg
            pipe = cold
            annulus = hot
            sgp = sgc
            sga = sgh
        gp = mp/ap
        ga = ma/aa
        De = ((D2**2)-(D1**2))/D1
        tw = (th1+th2+tc1+tc2)/4
        if data_input=='input':
            if (ap>aa and mh>mc) or (ap<aa and mh<mc):
                up = uh
                uwp = uwh
                ua = uc
                uwa = uwc
            else:
                up = uc
                uwp = uwc
                ua = uh
                uwa = uwh
        elif data_input=='database':  
            # Call the function to get the viscosity value
            up = get_viscosity_value(tpavg,pipe)
            uwp = get_viscosity_value(tw,pipe)
            # Call the function to get the viscosity value
            ua = get_viscosity_value(taavg,annulus)
            uwa = get_viscosity_value(tw,annulus)
        if data_input=='input':
            if (ap>aa and mh>mc) or (ap<aa and mh<mc):
                kp = kh
                ka = kc
            else:
                kp = kc
                ka = kh
        elif data_input=='database':
            # Call the function to get the conductivity value
            kp = get_conductivity_value(tpavg,pipe)
            # Call the function to get the conductivity value
            ka = get_conductivity_value(taavg,annulus)
    def shellandtubeheat_exchanger():
        global M
        M = float(entry_M.get())
        global nt
        nt = float(entry_nt.get())
        global n1
        n1 = float(entry_n1.get())
        global d
        d = float(entry_d.get())
        global do
        do = float(entry_do.get())
        global ds
        ds= float(entry_ds.get())
        global b
        b = float(entry_b.get())
        global pt
        pt = float(entry_pt.get())
        global rdg
        rdg = float(entry_rdg.get())
        global tp
        tp = entry_tp.get()
        global ut, uwt, us, uws, kt, ks,ttavg,tube,tw,tsavg,shell,ms,mt,cps,cpt,tsavg,ttavg,tube,shell,sgs,sgt,delpa,delpp,q,n,A,ud,Rd,l
        ms = mh
        mt = mc
        cps = cph
        cpt = cpc
        tsavg = thavg
        ttavg = tcavg
        tube = cold
        shell = hot
        sgs = sgh
        sgt = sgc
        if data_input=='input':
            ut = uc
            uwt = uwc
            us = uh
            uws = uwh
            
        elif data_input=='database':
            # Call the function to get the viscosity value
            ut = get_viscosity_value(ttavg,tube)
            uwt = get_viscosity_value(tw,tube)
            us = get_viscosity_value(tsavg,shell)
            uws = get_viscosity_value(tw,shell)
        if data_input=='input':
            kt = kc
            ks = kh 
        elif data_input=='database':
            #conductivity value for tube and shell  
            kt = get_conductivity_value(ttavg,tube)
            ks = get_conductivity_value(tsavg,shell)
        
    
    # for double pipe heat exchange+
    if heatexchanger.lower()=="doublepipeheatexchanger":
        doublepipeheat_exchanger()
        if doubletype.lower() == "parallel":
              epsi = (1/(mh*cph))+(1/(mc*cpc))
        elif doubletype.lower() == "countercurrent":
              epsi = (1/(mh*cph))-(1/(mc*cpc))
        gp = mp/ap
        ga = ma/aa
        De = ((D2**2)-(D1**2))/D1
        tw = (th1+th2+tc1+tc2)/4
        while True:
            phip = (up/uwp)**0.14
            phia = (ua/uwa)**0.14
            rep = (D*gp)/up
            rea = (De*ga)/ua
            # calculation of hi
            prp = (cpp*up)/kp
            if rep<2100:
                nup = 1.86*((rep*prp*D)/l)*phip
            elif rep>4000:
                nup = 0.027*(rep**0.8)*(prp**(1/3))*phip
            hi = (nup*kp)/D
            #calculation of ho
            pra = (cpa*ua)/ka
            if rea<2100:
                nua = 1.86*((rea*pra*D)/l)*phia
            elif rea>4000:
                nua = 0.027*(rea**0.8)*(pra**(1/3))*phia
            ho = (nua*ka)/De
            tw1 =((hi*tpavg)+(ho*tpavg*(D1/D)))/((hi+ho*(D1/D)))
            if data_input=='input':
                tw=tw1
            elif data_input=='database':
                if tw1!=tw:
                    tw=tw1
                    continue
    
            #to calclate the uc value
            if kw==-1:
                Uc=(hi*ho)/(hi+ho)
            else:
                Uc = 1/((D1/hi*D)+((D1/2*kw)*(math.log((D1/D))))+(1/ho))
            ud = 1/((1/Uc)+Rd)
            A = (math.log(abs((((2*tam)/q)+epsi)/(((2*tam)/q)-epsi))))/(ud*epsi)
            l1=A/(3.14*D1)
            n = l1/(2*l)
            n = math.floor(n)+1
            l3 = n*(2*l)
            if rep<2100 or rea<2100:
                if l!=l3:
                    l=l3
                    continue
                elif l1==l3 and tw==tw1:
                    break
            else:
                l=l3
                break
        A=3.14*D1*l
        ud=math.log(abs((((2*tam)/q)+epsi)/(((2*tam)/q)-epsi)))/(A*epsi)
        Rd=(Uc-ud)/(Uc*ud)
        De1=D2-D1
        rea1=De1*ga/ua
        # calculation for pipe
        if rep<2000:
            fp = 16/rep
        elif rep>4000:
            fp = 0.0035+(0.264/(D*gp/(up**0.42)))
        rhop = (sgp)*1000
        Fp = (4*fp*(gp**2)*l)/(2*(rhop**2)*D)
        #calculation of annulus
        if rea1<2000:
            fa = 16/rea
        elif rea1>4000:
            fa = 0.0035+(0.264/(De1*ga/(ua**0.42)))
        rhoa = (sga)*1000
        Fa = (4*fa*(ga**2)*l)/(2*(rhoa**2)*De1)
        #calculation of entrance loss
        v = ga/rhoa
        Fi=n*(v**2)/(2)
        delpp = Fp*rhop
        delpa = (Fa+Fi)*rhoa
        if delpp<p1_allowed and delpa<p_allowed:
            global desc_9, desc_10, desc_11, desc_12, desc_13, desc_14, desc_15, desc_16, desc_17, desc_18
            desc_9 = f"the calculated pressure drop for pipe is {delpp} and for for annulus is {delpa}." 
            desc_10 = f"heat duty is {q}"
            desc_11 = f"pressure drop of fluid in pipe is {delpp}"
            desc_12 =f"pressure drop of fluid in pipe is {delpa}"
            desc_13 = f"length of the heat exchanger is {l}"
            desc_14 = f"number of hairpins are {n}"
            desc_15 = f"area of heatexchanger is {A}"
            desc_16 = f"ud is {ud}"
            desc_17 = f"Rd is {Rd}"
            desc_18 = "the values of this is less than the allowed pressure.hence given heat exchanger is allowed to work."
            result_label1.config(text=f"{desc_9}\n{desc_10}\n{desc_11}\n{desc_12}\n{desc_13}\n{desc_14}\n{desc_15}\n{desc_16}\n{desc_17}\n{desc_18}")
            result_label1.grid(row=36, column=0, columnspan=2, pady=10)
            result_label2.grid_forget()
            result_label3.grid_forget()
            result_label4.grid_forget()
            
        else:
            desc_9 = "given heat exchanger is not allowed to work"
            desc_10 = f"heat duty is {q}"
            desc_11 = f"pressure drop of fluid in pipe is {delpp}"
            desc_12 = f"pressure drop of fluid in pipe is {delpa}"
            desc_13 = f"length of the heat exchanger is {l}"
            desc_14 = f"number of hairpins are {n}"
            desc_15 = f"area of heatexchanger is {A}"
            desc_16 = f"ud is {ud}"
            desc_17 = f"Rd is {Rd}"
            result_label2.config(text=f"{desc_9}\n{desc_10}\n{desc_11}\n{desc_12}\n{desc_13}\n{desc_14}\n{desc_15}\n{desc_16}\n{desc_17}")
            result_label2.grid(row=37, column=0, columnspan=2, pady=10)
            result_label1.grid_forget()
            result_label3.grid_forget()
            result_label4.grid_forget()
        if doubletype.lower() == "parallel":
            fig = Figure()
            ax = fig.add_subplot(111)
            ax.set_title("Tkinter Embedded Matplotlib Plot")
            b = 1/(mc*cpc)
            B = 1/(mh*cph)
            x = np.arange(math.ceil(l))
            thx = (((b/B)*np.ones(math.ceil(l))+np.exp(-ud*(B+b)*3.14*D*x))*th1+(np.ones(math.ceil(l))-np.exp(-ud*(B+b)*3.14*D*x))*tc1)/(1+(b/B))
            tcx = (((B/b)-1)*th1+(((B/b)*np.exp(ud*(B-b)*3.14*D*l))-(B/b))*tc1)/((B/b)*np.exp(ud*(B-b)*3.14*D*l)-np.ones(math.ceil(l)))+th1+(((B/b)-1)*th1+(((B/b)*np.exp(ud*(B-b)*3.14*D*l))-(B/b))*tc1)/((B/b)*np.exp(ud*(B-b)*3.14*D*l)-np.ones(math.ceil(l)))
            ax.plot(x, thx,label='hot')
            ax.plot(x, tcx, label='cold')
            ax.set_xlabel("X-axis")
            ax.set_ylabel("Y-axis")

            # Create a canvas widget to embed the plot into the Tkinter window
            canvas1 = FigureCanvasTkAgg(fig, master=window)
            canvas1.draw()

# Place the canvas widget using the grid layout
            canvas1.get_tk_widget().grid(row=40, column=0)
            plt.show()
        elif doubletype.lower() == "countercurrent":
            b = 1/(mc*cpc)
            B = 1/(mh*cph)
            x = np.arange(math.ceil(l))
            th1 = th1*np.ones(math.ceil(l))
            th2 = th2*np.ones(math.ceil(l))
            tc1 = tc1*np.ones(math.ceil(l))
            thx = ((((B/b)-1)*th1)+(((B/b)*np.exp(ud*(B-b)*3.14*D*x)-(B/b))*tc1))/((B/b)*np.exp(ud*(B-b)*3.14*D*x)-np.ones(math.ceil(l)))
            tcx = (b/B)*(thx-th2)+tc1
            plt.plot(x, thx,'r', linewidth=2.0)
            plt.plot(x, tcx,'b', linewidth=2.0)
            plt.show()
    # for shell and tube heat exchanger     
    elif heatexchanger.lower()=="shellandtubeheatexchanger":
        shellandtubeheat_exchanger()
        tlm = ((th1-tc2)-(th2-tc1))/(math.log((th1-tc2)/(th2-tc1)))
        R = (th1-th2)/(tc2-tc1)
        P1 = (tc2-tc1)/(th1-tc1)
        pn = (1-((1-P1*R)/(1-P1))**M)/(R-((1-P1*R)/(1-P1))**M)
        W = (1-pn*R)/(1-pn)
        c = (((R**2)+1)**(1/2))/(R-1)
        Fn = c*math.log(W**(1/M))/math.log((1+(W**(1/M))-c+c*(W**(1/M)))/(1+(W**(1/M))+c-c*(W**(1/M))))
        # F need to be calculated
        C =  pt-do
        at = (nt*3.14*(d**2))/(4*n1)
        As = (ds*C*b)/pt
        gt = mt/at
        gs = ms/As
        tw = (th1+th2+tc1+tc2)/4
        while True:
            phis = (us/uws)**0.14
            phit = (ut/uwt)**0.14
            if tp =="triangular":
                de=(((pt**2)*(3**(1/2)))-(1.57*(do**2)))/(1.57*do)
            elif tp =="square":
                de=(4*(pt**2)-3.14*(do**2))/(3.14*do)
            ret = (d*gt)/ut
            res = (de*gs)/us
            prt = (cpt*ut)/kt
            if ret<2100:
                global nut
                nut = 1.86*phit*((ret*prt*(d/l))**(1/3))
            elif ret>4000:
                nut = 0.027*phit*(prt**(1/3))*(ret**0.8)
            hi=(kt*nut)/d
            prs = (us*cps)/ks
            if res<2100:
                global nus
                nus = 1.86*phis*((res*prs*(d/l))**(1/3))
            elif res>4000:
                nus = 0.027*phis*(prs**(1/3))*(res**0.8)
            ho=(ks*nus)/de
            tw1 = (hi*ttavg+ho*tsavg*(do/d))/(hi+ho*(do/d))
            if data_input=='input':
                tw=tw1
            elif data_input=='database':
                if tw!=tw1:
                    tw=tw1
                    continue
            if kw==-1:
                Uc=(hi*ho)/(hi+ho)
            else:
                Uc=1/((do/(hi*d))+((do/(2*kw))*(math.log(do/d)))+(1/ho))
            A=(3.14*do*nt*l)
            ud = q/(A*Fn*tlm)
            rdc = (Uc-ud)/(Uc*ud)
            if rdg>rdc:
                print("calculated dirt factor is less than the given dirt factor hence the heatexchanger with the given specification is not allowed to work")
                break
            elif rdc>=rdg and tw==tw1:
                rhot = sgt*1000
                if ret<2100:
                    global ft
                    ft=16/ret
                elif ret>4000:
                    ft=(0.0035+(0.264/((d*gt/ut)**0.42)))
                PT = (4*ft*(gt**2)*l*n1)/(2*rhot*d*phit)
                pr = (4*n1*(gt**2)*(0.000532))/sgt
                pT = pr+PT
                rhos = sgs*1000
                N=(l/b)-1
                if 400<res<1000000:
                    global fs
                    fs = exp(0.576-0.19*(math.log(res)))
                I=math.floor(N)+1
                ps = (fs*(gs**2)*ds*(I+1))/(2*rhos*de*rhos)
                if ps<p_allowed and pT<p1_allowed:
                    global desc_1, desc_2, desc_3, desc_4, desc_5, desc_6, desc_7
                    desc_1 = "pressure drop level then the heat exchanger is allowed to work"
                    desc_2 = f'the area of the heat exchanger is {A}'
                    desc_3 = f"dirty heat transfer coefficient is {ud}"
                    desc_4 = f"clean heat transfer coefficient is {Uc}"
                    desc_5 = f"the heat duty is {q}"
                    desc_6 = f"pressure drop at tube side is {pT}"
                    desc_7 = f"pressure drop at shell side is {ps}"
                    result_label3.config(text=f"{desc_1}\n{desc_2}\n{desc_3}\n{desc_4}\n{desc_5}\n{desc_6}\n{desc_7}")
                    result_label3.grid(row=38, column=0, columnspan=2, pady=10)
                    result_label2.grid_forget()
                    result_label1.grid_forget()
                    result_label4.grid_forget()
                else:
                   global desc_8
                   desc_8 = "the heat exchanger is not allowed to work"
                   result_label4.config(text=f"{desc_8}")
                   result_label4.grid(row=39, column=0, columnspan=2, pady=10)
                   result_label2.grid_forget()
                   result_label3.grid_forget()
                   result_label1.grid_forget()
                break
def show_inputfunc():
    
    label_cpc.grid(row=16, column=0, padx=10, pady=5, sticky="e")
    entry_cpc.grid(row=16, column=1, padx=10, pady=5, sticky="w")
    label_cph.grid(row=17, column=0, padx=10, pady=5, sticky="e")
    entry_cph.grid(row=17, column=1, padx=10, pady=5, sticky="w")
    label_uh.grid(row=18, column=0, padx=10, pady=5, sticky="e")
    entry_uh.grid(row=18, column=1, padx=10, pady=5, sticky="w")
    label_uc.grid(row=19, column=0, padx=10, pady=5, sticky="e")
    entry_uc.grid(row=19, column=1, padx=10, pady=5, sticky="w")
    label_uwh.grid(row=20, column=0, padx=10, pady=5, sticky="e")
    entry_uwh.grid(row=20, column=1, padx=10, pady=5, sticky="w")
    label_uwc.grid(row=21, column=0, padx=10, pady=5, sticky="e")
    entry_uwc.grid(row=21, column=1, padx=10, pady=5, sticky="w")
    label_kh.grid(row=22, column=0, padx=10, pady=5, sticky="e")
    entry_kh.grid(row=22, column=1, padx=10, pady=5, sticky="w")
    label_kc.grid(row=23, column=0, padx=10, pady=5, sticky="e")
    entry_kc.grid(row=23, column=1, padx=10, pady=5, sticky="w")
def show_database_func():
    label_cpc.grid_forget()
    entry_cpc.grid_forget()
    label_cph.grid_forget()
    entry_cph.grid_forget()
    label_uh.grid_forget()
    entry_uh.grid_forget()
    label_uc.grid_forget()
    entry_uc.grid_forget()
    label_uwh.grid_forget()
    entry_uwh.grid_forget()
    label_uwc.grid_forget()
    entry_uwc.grid_forget()
    label_kh.grid_forget()
    entry_kh.grid_forget()
    label_kc.grid_forget()
    entry_kc.grid_forget()
def show_doublepipeheat_exchanger():
    label_doubletype.grid(row=25, column=0, padx=10, pady=5, sticky="e")
    entry_doubletype.grid(row=25, column=1, padx=10, pady=5, sticky="w")
    label_D.grid(row=26, column=0, padx=10, pady=5, sticky="e")
    entry_D.grid(row=26, column=1, padx=10, pady=5, sticky="w")
    label_D1.grid(row=27, column=0, padx=10, pady=5, sticky="e")
    entry_D1.grid(row=27, column=1, padx=10, pady=5, sticky="w")
    label_D2.grid(row=28, column=0, padx=10, pady=5, sticky="e")
    entry_D2.grid(row=28, column=1, padx=10, pady=5, sticky="w")


    label_M.grid_forget()
    entry_M.grid_forget()
    label_nt.grid_forget()
    entry_nt.grid_forget()
    label_n1.grid_forget()
    entry_n1.grid_forget()
    label_d.grid_forget()
    entry_d.grid_forget()
    label_do.grid_forget()
    entry_do.grid_forget()
    label_ds.grid_forget()
    entry_ds.grid_forget()
    label_b.grid_forget()
    entry_b.grid_forget()
    label_pt.grid_forget()
    entry_pt.grid_forget()
    label_rdg.grid_forget()
    entry_rdg.grid_forget()
    label_tp.grid_forget()
    entry_tp.grid_forget()
    
    

def show_shellandtubeheatexchanger():
    label_doubletype.grid_forget()
    entry_D.grid_forget()
    label_D.grid_forget()
    entry_D.grid_forget()
    label_D1.grid_forget()
    entry_D1.grid_forget()
    label_D2.grid_forget()
    entry_D2.grid_forget()

    label_M.grid(row=25, column=0, padx=10, pady=5, sticky="e")
    entry_M.grid(row=25, column=1, padx=10, pady=5, sticky="w")
    label_nt.grid(row=26, column=0, padx=10, pady=5, sticky="e")
    entry_nt.grid(row=26, column=1, padx=10, pady=5, sticky="w")
    label_n1.grid(row=27, column=0, padx=10, pady=5, sticky="e")
    entry_n1.grid(row=27, column=1, padx=10, pady=5, sticky="w")
    label_d.grid(row=28, column=0, padx=10, pady=5, sticky="e")
    entry_d.grid(row=28, column=1, padx=10, pady=5, sticky="w")
    label_do.grid(row=29, column=0, padx=10, pady=5, sticky="e")
    entry_do.grid(row=29, column=1, padx=10, pady=5, sticky="w")
    label_ds.grid(row=30, column=0, padx=10, pady=5, sticky="e")
    entry_ds.grid(row=30, column=1, padx=10, pady=5, sticky="w")
    label_b.grid(row=31, column=0, padx=10, pady=5, sticky="e")
    entry_b.grid(row=31, column=1, padx=10, pady=5, sticky="w")
    label_pt.grid(row=32, column=0, padx=10, pady=5, sticky="e")
    entry_pt.grid(row=32, column=1, padx=10, pady=5, sticky="w")
    label_rdg.grid(row=33, column=0, padx=10, pady=5, sticky="e")
    entry_rdg.grid(row=33, column=1, padx=10, pady=5, sticky="w")
    label_tp.grid(row=34, column=0, padx=10, pady=5, sticky="e")
    entry_tp.grid(row=34, column=1, padx=10, pady=5, sticky="w")


# building Gui
window = tk.Tk()
window.title("Heat Exchanger Design")
# Configure window size and background color

screen_width = window.winfo_screenwidth()
screen_height = window.winfo_screenheight()

# Set the window size to match the screen size
window.geometry(f"{screen_width}x{screen_height}")


window.configure(bg="#F0F0F0")

canvas = tk.Canvas(window)
canvas.grid(row=0, column=0, sticky="nsew")

frame = tk.Frame(canvas)
canvas.create_window((0, 0), window=frame, anchor="nw")
scrollbar =tk.Scrollbar(window, orient=tk.VERTICAL, command=canvas.yview)
scrollbar.grid(row=0, column=20, sticky="ns")
canvas.config(yscrollcommand=scrollbar.set)


# Function to set a consistent style for labels and entries
def set_style(widget):
    widget.configure(bg="#F0F0F0", font=("Arial", 12))
choice = tk.IntVar()
label_mc = Label(frame, text="Mass flow rate of cold fluid(Kg/s):")
set_style(label_mc)
label_mc.grid(row=0, column=0, padx=10, pady=5, sticky="e")
entry_mc = Entry(frame)
set_style(entry_mc)
entry_mc.grid(row=0, column=1, padx=10, pady=5, sticky="w")

label_mh = Label(frame, text="Mass flow rate of hot fluid(kg/s):")
set_style(label_mh)
label_mh.grid(row=1, column=0, padx=10, pady=5, sticky="e")
entry_mh = Entry(frame)
set_style(entry_mh)
entry_mh.grid(row=1, column=1, padx=10, pady=5, sticky="w")
  
label_tc1 = Label(frame, text="inlet temperature of cold fluid(K):")
set_style(label_tc1)
label_tc1.grid(row=2, column=0, padx=10, pady=5, sticky="e")
entry_tc1 = Entry(frame)
set_style(entry_tc1)
entry_tc1.grid(row=2, column=1, padx=10, pady=5, sticky="w")

label_th1 = Label(frame, text="inlet temperature of hot fluid(K):")
set_style(label_th1)
label_th1.grid(row=3, column=0, padx=10, pady=5, sticky="e")
entry_th1 = Entry(frame)
set_style(entry_th1)
entry_th1.grid(row=3, column=1, padx=10, pady=5, sticky="w")

label_th2 = Label(frame, text="outlet temperature of hot fluid(K):")
set_style(label_th2)
label_th2.grid(row=4, column=0, padx=10, pady=5, sticky="e")
entry_th2 = Entry(frame)
set_style(entry_th2)
entry_th2.grid(row=4, column=1, padx=10, pady=5, sticky="w")

label_tc2 = Label(frame, text="outlet temperature of cold fluid(K):")
set_style(label_tc2)
label_tc2.grid(row=5, column=0, padx=10, pady=5, sticky="e")
entry_tc2 = Entry(frame)
set_style(entry_tc2)
entry_tc2.grid(row=5, column=1, padx=10, pady=5, sticky="w")

label_l = Label(frame, text="length of one hairpin(m):")
set_style(label_l)
label_l.grid(row=6, column=0, padx=10, pady=5, sticky="e")
entry_l = Entry(frame)
set_style(entry_l)
entry_l.grid(row=6, column=1, padx=10, pady=5, sticky="w")

label_kw = Label(frame, text="kw:")
set_style(label_kw)
label_kw.grid(row=7, column=0, padx=10, pady=5, sticky="e")
entry_kw = Entry(frame)
set_style(entry_kw)
entry_kw.grid(row=7, column=1, padx=10, pady=5, sticky="w")


label_hot = Label(frame, text="name of hot fluid:")
set_style(label_hot)
label_hot.grid(row=8, column=0, padx=10, pady=5, sticky="e")
entry_hot = Entry(frame)
set_style(entry_hot)
entry_hot.grid(row=8, column=1, padx=10, pady=5, sticky="w")

label_cold = Label(frame, text="name of cold fluid:")
set_style(label_cold)
label_cold.grid(row=9, column=0, padx=10, pady=5, sticky="e")
entry_cold = Entry(frame)
set_style(entry_cold)
entry_cold.grid(row=9, column=1, padx=10, pady=5, sticky="w")

label_Rd = Label(frame, text="allowable dirt factor(K(m^2)/W):")
set_style(label_Rd)
label_Rd.grid(row=10, column=0, padx=10, pady=5, sticky="e")
entry_Rd = Entry(frame)
set_style(entry_Rd)
entry_Rd.grid(row=10, column=1, padx=10, pady=5, sticky="w")

label_p_allowed = Label(frame, text="enter the allowed pressure shell or annulus(pascal):")
set_style(label_p_allowed)
label_p_allowed.grid(row=11, column=0, padx=10, pady=5, sticky="e")
entry_p_allowed = Entry(frame)
set_style(entry_p_allowed)
entry_p_allowed.grid(row=11, column=1, padx=10, pady=5, sticky="w")

label_p1_allowed = Label(frame, text="enter the allowed pressure tube or pipe(pascal):")
set_style(label_p1_allowed)
label_p1_allowed.grid(row=12, column=0, padx=10, pady=5, sticky="e")
entry_p1_allowed = Entry(frame)
set_style(entry_p1_allowed)
entry_p1_allowed.grid(row=12, column=1, padx=10, pady=5, sticky="w")

label_sgh = Label(frame, text="specific gravity of hot fluid:")
set_style(label_sgh)
label_sgh.grid(row=13, column=0, padx=10, pady=5, sticky="e")
entry_sgh = Entry(frame)
set_style(entry_sgh)
entry_sgh.grid(row=13, column=1, padx=10, pady=5, sticky="w")

label_sgc = Label(frame, text="specific gravity of cold fluid:")
set_style(label_sgc)
label_sgc.grid(row=14, column=0, padx=10, pady=5, sticky="e")
entry_sgc = Entry(frame)
set_style(entry_sgc)
entry_sgc.grid(row=14, column=1, padx=10, pady=5, sticky="w")

label_data_type = Label(frame, text="data type")
set_style(label_data_type)
label_data_type.grid(row=15, column=0, padx=10, pady=5, sticky="e")


data1 = tk.StringVar()

radio_input = tk.Radiobutton(frame, text="input", variable=data1, value="input", command = show_inputfunc)
radio_input.grid(row=15,column=1, padx=10, pady=5)

radio_database = tk.Radiobutton(frame, text="database", variable=data1, value="database", command = show_database_func)
radio_database.grid(row=15,column=2,padx=10,pady=5)

label_cpc = Label(frame, text="specific heat capacity of cold fluid:")
set_style(label_cpc)
entry_cpc = Entry(frame)
set_style(entry_cpc)


label_cph = Label(frame, text="specific heat capacity of hot fluid:")
set_style(label_cph)
entry_cph = Entry(frame)
set_style(entry_cph)


label_uh = Label(frame, text="viscosity of hot fluid at average temperature(kg/m*s):")
set_style(label_uh)
entry_uh = Entry(frame)
set_style(entry_uh)

label_uc = Label(frame, text="viscosity of cold fluid at average temperature(kg/m*s):")
set_style(label_uc)
entry_uc = Entry(frame)
set_style(entry_uc)


label_uwh = Label(frame, text="viscosity of hot fluid at wall temperature(kg/m*s):")
set_style(label_uwh)
entry_uwh = Entry(frame)
set_style(entry_uwh)


label_uwc = Label(frame, text="viscosity of cold fluid at wall temperature(kg/m*s):")
set_style(label_uwc)
entry_uwc = Entry(frame)
set_style(entry_uwc)


label_kh = Label(frame, text="thermal conductivity of hot fluid at average temperature(W/m*K):")
set_style(label_kh)
entry_kh = Entry(frame)
set_style(entry_kh)


label_kc = Label(frame, text="thermal conductivity of cold fluid at average temperature(W/m*K):")
set_style(label_kc)
entry_kc = Entry(frame)
set_style(entry_kc)

label_heatexchanger = Label(frame, text="heat exchanger")
set_style(label_heatexchanger)
label_heatexchanger.grid(row=24, column=0, padx=10, pady=5, sticky="e")


heat1 = tk.StringVar()

radio_doublepipe = tk.Radiobutton(frame, text="double pipe", variable=heat1, value="doublepipeheatexchanger", command = show_doublepipeheat_exchanger)
radio_doublepipe.grid(row=24, column=1, padx=10, pady=5)

radio_shellandtube = tk.Radiobutton(frame, text="shell and tube", variable=heat1, value="shellandtubeheatexchanger", command = show_shellandtubeheatexchanger)
radio_shellandtube.grid(row=24, column=2, padx=10, pady=5)

label_doubletype = Label(frame, text="parallel or countercurrent:")
set_style(label_doubletype)
entry_doubletype = Entry(frame)
set_style(entry_doubletype)


label_D = Label(frame, text="inner diameter of pipe(m):")
set_style(label_D)
entry_D = Entry(frame)
set_style(entry_D)


label_D1 = Label(frame, text="outer diameter of pipe(m):")
set_style(label_D1)
entry_D1 = Entry(frame)
set_style(entry_D1)


label_D2 = Label(frame, text="inner diameter of annulus(m):")
set_style(label_D2)
entry_D2 = Entry(frame)
set_style(entry_D2)


label_M = Label(frame, text="Number of shell passes:")
set_style(label_M)
entry_M = Entry(frame)
set_style(entry_M)


label_nt = Label(frame, text="number of tubes:")
set_style(label_nt)
entry_nt = Entry(frame)
set_style(entry_nt)

label_n1 = Label(frame, text="number of tube passes(m):")
set_style(label_n1)
entry_n1 = Entry(frame)
set_style(entry_n1)


label_d = Label(frame, text="inner diameter of tube(m):")
set_style(label_d)
entry_d = Entry(frame)
set_style(entry_d)


label_do = Label(frame, text="outer diameter of tube(m):")
set_style(label_do)
entry_do = Entry(frame)
set_style(entry_do)


label_ds = Label(frame, text="inner diameter of shell(m):")
set_style(label_ds)
entry_ds = Entry(frame)
set_style(entry_ds)


label_b = Label(frame, text="baffle spacing:")
set_style(label_b)
entry_b = Entry(frame)
set_style(entry_b)

label_pt = Label(frame, text="tube pitch:")
set_style(label_pt)
entry_pt = Entry(frame)
set_style(entry_pt)

label_rdg = Label(frame, text="allowable dirt factor:")
set_style(label_rdg)
entry_rdg = Entry(frame)
set_style(entry_rdg)

label_tp = Label(frame, text="type of pitch(square or triangular):")
set_style(label_tp)
entry_tp = Entry(frame)
set_style(entry_tp)

def on_frame_configure(event):
    canvas.configure(scrollregion=canvas.bbox("all"))

frame.bind("<Configure>", on_frame_configure)


calculate_button = Button(frame, text="Calculate", command=solve_equations, bg="#4CAF50", fg="white", font=("Arial", 14))
calculate_button.grid(row=35, column=0, columnspan=2, pady=20)
   

window.grid_rowconfigure(0, weight=1)
window.grid_columnconfigure(0, weight=1)

result_label1 = Label(frame, text="", bg="#F0F0F0", font=("Arial", 14))


result_label2 = Label(frame, text="", bg="#F0F0F0", font=("Arial", 14))


result_label3 = Label(frame, text="", bg="#F0F0F0", font=("Arial", 14))


result_label4 = Label(frame, text="", bg="#F0F0F0", font=("Arial", 14))




window.mainloop()

end_time = time.time()
runtime = end_time-start_time
print(f"runtime is {runtime}")
