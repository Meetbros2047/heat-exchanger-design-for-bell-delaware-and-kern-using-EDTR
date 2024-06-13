# -*- coding: utf-8 -*-
"""
Created on Mon May 27 00:55:02 2024

@author: mitansh waghmare
"""
import tkinter as tk
from tkinter import Label, Entry, Button
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
        global l
        l = float(entry_l.get())
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
        global question
        question = entry_question.get()
        global Np
        Np = float(entry_Np.get())
        global Wp
        Wp = float(entry_Wp.get())
        global lsb
        lsb = float(entry_lsb.get())
        global ltb
        ltb = float(entry_ltb.get())
        global lbb
        lbb = float(entry_lbb.get())
        global ntp
        ntp = float(entry_ntp.get())
        global n1
        n1 = float(entry_n1.get())
        global d
        d = float(entry_d.get())
        global do
        do = float(entry_do.get())
        global ds
        ds = float(entry_ds.get())
        global lb
        lb = float(entry_lb.get())
        global pt
        pt = float(entry_pt.get())
        global rds
        rds = float(entry_rds.get())
        global rdt
        rdt = float(entry_rdt.get())
        global bpercent
        bpercent = float(entry_bpercent.get())
        global layout_type
        layout_type = float(entry_layout_type.get())
        global lbout
        lbout = float(entry_lbout.get())
        global lbin
        lbin = float(entry_lbin.get())
        global nss
        nss = float(entry_nss.get())
        global pitch
        pitch = entry_pitch.get()
        global l5
        l5 = float(entry_l5.get())
        global ms,mt,cps,cpt,tsavg,ttavg,tube,shell,sgs,sgt,ut,uwt,us,uws,kt,ks
        # F need to be calculated
        if question=='shell':
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
        elif question=='tube':
            ms = mc
            mt = mh
            cps = cpc
            cpt = cph
            tsavg = tcavg
            ttavg = thavg
            tube = hot
            shell = cold
            sgs = sgc
            sgt = sgh
        if data_input=='input':
            if question=='shell':
                ut = float(entry_uc.get())
                uwt = float(entry_uwc.get())
                us = float(entry_uh.get())
                uws = float(entry_uwh.get())
            elif question=='tube':
                us = float(entry_uc.get())
                uws = float(entry_uwc.get())
                ut = float(entry_uh.get())
                uwt = float(entry_uwh.get())
        elif data_input=='database':
            # Call the function to get the viscosity value
            ut = get_viscosity_value(ttavg,tube)
            uwt = get_viscosity_value(tw,tube)
            us = get_viscosity_value(tsavg,shell)
            uws = get_viscosity_value(tw,shell)
            
        if data_input=='input':
            if question=='shell':
                kt = float(entry_kc.get())
                ks = float(entry_kh.get())    
            elif question=='tube':
                ks = float(entry_kc.get())
                kt= float(entry_kh.get())
        elif data_input=='database':
            #conductivity value for tube and shell   
            kt = get_conductivity_value(ttavg,tube)
            ks = get_conductivity_value(tsavg,shell)
    # for double pipe heat exchange
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
            desc_9 = f"the calculated pressure drop for pipe is {delpp} pascal and for for annulus is {delpa} pascal." 
            desc_10 = f"heat duty is {q} watt"
            desc_11 = f"pressure drop of fluid in pipe is {delpp} pascal"
            desc_12 =f"pressure drop of fluid in pipe is {delpa} pascal"
            desc_13 = f"length of the heat exchanger is {l} m"
            desc_14 = f"number of hairpins are {n}"
            desc_15 = f"area of heatexchanger is {A} m^2"
            desc_16 = f"ud is {ud} watt/(m^2)*K"
            desc_17 = f"Rd is {Rd} K*(m^2)/W"
            desc_18 = "the values of this is less than the allowed pressure.hence given heat exchanger is allowed to work."
            result_label1.config(text=f"{desc_9}\n{desc_10}\n{desc_11}\n{desc_12}\n{desc_13}\n{desc_14}\n{desc_15}\n{desc_16}\n{desc_17}\n{desc_18}")
            result_label1.grid(row=48, column=0, columnspan=2, pady=10)
            result_label2.grid_forget()
            result_label3.grid_forget()
            result_label4.grid_forget()
            result_label5.grid_forget()
            
        else:
            desc_9 = "given heat exchanger is not allowed to work"
            desc_10 = f"heat duty is {q} Watt"
            desc_11 = f"pressure drop of fluid in pipe is {delpp} Pascal"
            desc_12 = f"pressure drop of fluid in pipe is {delpa} pascal"
            desc_13 = f"length of the heat exchanger is {l} m"
            desc_14 = f"number of hairpins are {n}"
            desc_15 = f"area of heatexchanger is {A} m^2"
            desc_16 = f"ud is {ud} watt/(m^2)*K"
            desc_17 = f"Rd is {Rd} K*(m^2)/W"
            result_label2.config(text=f"{desc_9}\n{desc_10}\n{desc_11}\n{desc_12}\n{desc_13}\n{desc_14}\n{desc_15}\n{desc_16}\n{desc_17}")
            result_label2.grid(row=49, column=0, columnspan=2, pady=10)
            result_label1.grid_forget()
            result_label3.grid_forget()
            result_label4.grid_forget()
            result_label5.grid_forget()
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
        
    # for shell and tube heat exchanger
    elif heatexchanger.lower()=="shellandtubeheatexchanger":
        shellandtubeheat_exchanger()
        epsi = math.sqrt(((mc*cpc)**(-2))+((mh*cph)**(-2)))
        nt = ntp*n1
        #C =  pt-do
        at = (nt*3.14*(d**2))/(4*n1)
        gt = mt/at
        tw = (th1+th2+tc1+tc2)/4
        N = ((l5-0.2*ds)/lb)-1
        N = math.floor(N)
        while True:
            global Np,Wp,lsb,ltb,lbb,do,pt,rds,rdt,bpercent,layout_type,lbout,lbin,nss,pitch,ms,cps,cpt,tsavg,ttavg,tube,shell,sgs,sgt,ut,uwt,us,uws,kt,ks
            phis = (us/uws)**0.14
            phit = (ut/uwt)**0.14
            ret = (d*gt)/ut
            prt = (cpt*ut)/kt
            if ret<2100:
                nut = 1.86*phit*((ret*prt*(d/l5))**(1/3))
            elif ret>4000:
                nut = 0.027*phit*(prt**(1/3))*(ret**0.8)
            hi=(kt*nut)/d
            prs = (us*cps)/ks
            def re_constant(res):
                if 10000<res<=100000:
                    rek = '100000-10000'
                elif 1000<res<=10000:
                    rek = '10000-1000'
                elif 100<res<=1000:
                    rek = '1000-100'
                elif 10<res<=100:
                    rek = '100-10'
                elif res<10:
                    rek = '<10'
                return rek
                    
            def diff_constants(res2,layout_type):
                excel_file = 'value_of_constants.xlsx'
                constants = 'Sheet1'
                constants_df = pd.read_excel(excel_file, sheet_name=constants)
                mask = (constants_df['Reynolds Number'] == res2) & (constants_df['layout angle'] == layout_type)
                constants_value = constants_df.loc[mask, 'a1':'b4'].values.flatten().astype(float)
                return constants_value
            def k1n2(n1,pitch):
                excel_file = 'value_of_constants.xlsx'
                kconstant = 'Sheet4'
                kconstant_df = pd.read_excel(excel_file, sheet_name=kconstant)
                if pitch=='triangular':
                    kconstantk1 = kconstant_df.loc[kconstant_df['Number Passes'] == n1, 'TPK1'].values.flatten().astype(float)
                    kconstantn2 = kconstant_df.loc[kconstant_df['Number Passes'] == n1, 'TPN2'].values.flatten().astype(float)
                elif pitch=='square':
                    kconstantk1 = kconstant_df.loc[kconstant_df['Number Passes'] == n1, 'SPK1'].values.flatten().astype(float)
                    kconstantn2 = kconstant_df.loc[kconstant_df['Number Passes'] == n1, 'SPN2'].values.flatten().astype(float)
                return kconstantk1, kconstantn2
            new_constant = k1n2(n1,pitch)
            k1 = new_constant[0]
            n2 = new_constant[1]
            dotl = do*(nt/k1)**(1/n2)
            if layout_type==30:
                pn = (pt/2.0)
                sm = lb*(ds-dotl+((dotl-do)/pt)*(pt-do))
                pp = (math.sqrt(3)/2)*pt
            elif layout_type==45:
                pn = pt/math.sqrt(2)
                sm = lb*(ds-dotl+((dotl-do)/pn)*(pt-do))
                pp = pt/math.sqrt(2)
            elif layout_type==90:
                pn = pt
                sm = lb*(ds-dotl+((dotl-do)/pt)*(pt-do))
                pp = pt
            res = (do*ms)/(us*sm) 
            lc = (bpercent*ds)/100
            nc=math.floor((ds/pp)*(1-(2*(bpercent/100))))
            res2 = re_constant(res)
            gs = ms/sm
            constant = diff_constants(res2,layout_type) 
            a1 = (constant[0])
            a2 = (constant[1])
            a3 = (constant[2])
            a4 = (constant[3])
            b1 = (constant[4])
            b2 = (constant[5])
            b3 = (constant[6])
            b4 = (constant[7])
            a = a3/(1+0.14*(res**a4))
            jh = a1*(res**a2)*(((1.33*do)/pt)**a)
            hid = jh*gs*cps*((prs)**(-2.0/3))*phis
            ftc = (1/math.pi)*(math.pi+2*((ds-2*lc)/(dotl-do))*math.sin(math.acos((ds-2*lc)/(dotl-do)))-2*math.acos((ds-2*lc)/(dotl-do)))
            jc = 0.55+0.72*ftc
            stb = (pi*ltb*do*nt*(1+ftc))/(4.0)
            ssb = ((ds* lsb)/2)*(math.pi-math.acos(1-((2*lc)/ds)))
            alpha = 0.44 * (stb/(stb+ssb))
            jl = alpha + (1-alpha)*(math.exp(-2.2*((stb+ssb)/sm)))
            fbp = ((ds-dotl+0.5*Np*Wp)/sm)*lb 
            if res<=100:
                cbh = 1.35
            elif res>100:
                cbh = 1.25
            if (nss/nc)>=0.5:
                jb = 1
            elif (nss/nc)<0.5:
                jb = math.exp(-cbh*fbp*(1-((2*nss/nc)**(1/3))))
            ncw = (0.8*(lc-0.5*(ds+do-dotl)))/(pp)
            nr = (nc+ncw)*(N+1)
            jrr = 1.51/(nr**0.18)
            if res>100:
                jr = 1
            elif 20<=res<=100:
                jr = jrr-((1-jrr)*(0.25-0.0125*res))
            elif res<20:
                jr = jrr
            if res<100:
                n3 = 0.333
            elif res>100:
                n3 = 0.6
            lin = (lbin/lb)
            lout = (lbout/lb)
            if ((0.95*lb)<=lbin<=(1.05*lb)) and ((0.95*lb)<=lbout<=(1.05*lb)):
                js = 1
            else:
                js = ((N-1)+(lin**(1-n3))+(lout**(1-n3)))/(N-1+lin+lout)
            ho = hid*jc*jl*jb*jr*js
            tw1 = (hi*ttavg+ho*tsavg*(do/d))/(hi+ho*(do/d))
            if data_input=='input':
                tw=tw1
            elif data_input=='database':
                if tw!=tw1:
                    tw=tw1
                    continue
            if kw==-1:
                Uc=1/((do/(hi*d))+(1/ho))
            else:
                Uc=1/((do/(hi*d))+((do/(2*kw))*(math.log(do/d)))+(1/ho))
            if kw==-1:
                ud2=1/((do/(hi*d))+(1/ho)+((do/d)*rdt)+rds)
            else:
                ud2=1/((do/(hi*d))+((do/(2*kw))*(math.log(do/d)))+(1/ho)+((do/d)*rdt)+rds)
            A=(3.14*do*nt*l5)
            ud1 = math.log(abs((((2*tam)/q)+epsi)/(((2*tam)/q)-epsi)))/(A*epsi)
            rdo = (1/ud2)-(1/ud1)
            if rdo<0:
                desc_10 = "The given heat exchanger is underdesigned according to the specified fouling factor."
                result_label5.config(text=f"{desc_10}")
                result_label5.grid(row=50, column=0, columnspan=2, pady=10)
                result_label2.grid_forget()
                result_label1.grid_forget()
                result_label4.grid_forget()
                result_label3.grid_forget()
                break
            elif rdo>0:
                rhot = sgt*1000
                if ret<2100:
                    ft=16/ret
                elif ret>4000:
                    ft=(0.0035+(0.264/((d*gt/ut)**0.42)))
                delpt = (2*ft*(gt**2)*l5*n1)/(2*rhot*d*phit)
                delpr = (2*n1*(gt**2))/(rhot)
                delpT = delpt + delpr
                rhos = sgs*1000
                I=(l5/lb)-1
                I=math.floor(I)
                b = b3/(1+0.14*(res**b4))
                fs = b1*(res**b2)*(((1.33*do)/pt)**b)
                pid = (2*fs*nc*(gs**2))/(rhos*phis)
                sw = ((ds**2)/4)*((math.acos(1-2*(lc/ds)))-(1-2*(lc/ds))*(math.sqrt(1-(1-2*(lc/ds))**2)))-((nt/8)*(1-ftc)*math.pi*(do**2))
                dw = (4*sw)/((math.pi/2)*nt*(1-ftc)*do+(ds*math.acos(1-2*(lc/ds))))
                if res<=100:
                    pw = 0.026*(((ms)*us)/(rhos*math.sqrt(sm*sw)))*((ncw/(pt-do))+(lb/(dw**2)))+(ms**2)/(rhos*math.sqrt(sm*sw))
                elif res>100:
                    pw = (2+(0.6*ncw))*((ms**2)/(2*rhos*math.sqrt(sm*sw)))
                n4 = -0.15*(1+(ssb/(ssb+stb)))+0.8
                rl = math.exp(-1.33*(1+(ssb/(ssb+stb)))*(((stb+ssb)/sm)**n4))
                if res<=100:
                    cbp = 4.5
                elif res>100:
                    cbp = 3.7
                if (nss/nc)>=0.5:
                    rb = 1 
                elif (nss/nc)<0.5:
                    rb = math.exp(-cbp*fbp*(1-((nss/nc)**(1/3))))
                if res<=100:
                    n5 = 1
                elif res>100:
                    n5 = 0.2
                if ((0.95*lb)<=lbin<=(1.05*lb)) and ((0.95*lb)<=lbout<=(1.05*lb)):
                    rs = 2
                else:
                    rs = (lin**(n5-2))+(lout**(n5-2))
                delps = rb*(pid)*((N-1)*rl+2*(1+(ncw/nc))*rs)+N*rl*pw
                if delps<p_allowed and delpT<p1_allowed:
                    global desc_1, desc_2, desc_3, desc_4, desc_5, desc_6, desc_7
                    desc_1 = "calculated pressure drop is within the limit then the give heat exchanger is allowed to work"
                    desc_2 = f'the area of the heat exchanger is {A} m^2'
                    desc_3 = f"dirty heat transfer coefficient is {ud2} watt/m^2*K"
                    desc_4 = f"clean heat transfer coefficient is {Uc} watt/m^2*K"
                    desc_5 = f"the heat duty is {q} watt"
                    desc_6 = f"pressure drop at tube side is {delpT} pascal"
                    desc_7 = f"pressure drop at shell side is {delps} pascal"
                    result_label3.config(text=f"{desc_1}\n{desc_2}\n{desc_3}\n{desc_4}\n{desc_5}\n{desc_6}\n{desc_7}")
                    result_label3.grid(row=50, column=0, columnspan=2, pady=10)
                    result_label2.grid_forget()
                    result_label1.grid_forget()
                    result_label4.grid_forget()
                    result_label5.grid_forget()
                else:
                   global desc_8
                   desc_8 = "calculated pressure drop is not within the  permissibe limit then the give heat exchanger is not allowed to work"
                   result_label4.config(text=f"{desc_8}")
                   result_label4.grid(row=51, column=0, columnspan=2, pady=10)
                   result_label2.grid_forget()
                   result_label3.grid_forget()
                   result_label1.grid_forget()
                   result_label5.grid_forget()
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
    label_l.grid(row=29, column=0, padx=10, pady=5, sticky="e")
    entry_l.grid(row=29, column=1, padx=10, pady=5, sticky="w")

   
    label_M.grid_forget()
    entry_M.grid_forget()
    label_question.grid_forget()
    entry_question.grid_forget()
    label_Np.grid_forget()
    entry_Np.grid_forget()
    label_Wp.grid_forget()
    entry_Wp.grid_forget()
    label_lsb.grid_forget()
    entry_lsb.grid_forget()
    label_ltb.grid_forget()
    entry_ltb.grid_forget()
    label_lbb.grid_forget()
    entry_lbb.grid_forget()
    label_ntp.grid_forget()
    entry_ntp.grid_forget()
    label_n1.grid_forget()
    entry_n1.grid_forget()
    label_d.grid_forget()
    entry_d.grid_forget()
    label_do.grid_forget()
    entry_do.grid_forget()
    label_ds.grid_forget()
    entry_ds.grid_forget()
    label_lb.grid_forget()
    entry_lb.grid_forget()
    label_pt.grid_forget()
    entry_pt.grid_forget()
    label_rds.grid_forget()
    entry_rds.grid_forget()
    label_rdt.grid_forget()
    entry_rdt.grid_forget()
    label_bpercent.grid_forget()
    entry_bpercent.grid_forget()
    label_layout_type.grid_forget()
    entry_layout_type.grid_forget()
    label_lbout.grid_forget()
    entry_lbout.grid_forget()
    label_lbin.grid_forget()
    entry_lbin.grid_forget()
    label_nss.grid_forget()
    entry_nss.grid_forget()
    label_pitch.grid_forget()
    entry_pitch.grid_forget()
    label_l5.grid_forget()
    entry_l5.grid_forget()
def show_shellandtubeheatexchanger():
    label_doubletype.grid_forget()
    entry_doubletype.grid_forget()
    label_D.pack_forget()
    entry_D.pack_forget()
    label_D1.pack_forget()
    entry_D1.pack_forget()
    label_D2.pack_forget()
    entry_D2.pack_forget()
    label_l.grid_forget()
    entry_l.grid_forget()
    

    label_M.grid(row=25, column=0, padx=10, pady=5, sticky="e")
    entry_M.grid(row=25, column=1, padx=10, pady=5, sticky="w")
    label_question.grid(row=26, column=0, padx=10, pady=5, sticky="e")
    entry_question.grid(row=26, column=1, padx=10, pady=5, sticky="w")
    label_Np.grid(row=27, column=0, padx=10, pady=5, sticky="e")
    entry_Np.grid(row=27, column=1, padx=10, pady=5, sticky="w")
    label_Wp.grid(row=28, column=0, padx=10, pady=5, sticky="e")
    entry_Wp.grid(row=28, column=1, padx=10, pady=5, sticky="w")
    label_ntp.grid(row=29, column=0, padx=10, pady=5, sticky="e")
    entry_ntp.grid(row=29, column=1, padx=10, pady=5, sticky="w")
    label_lsb.grid(row=30, column=0, padx=10, pady=5, sticky="e")
    entry_lsb.grid(row=30, column=1, padx=10, pady=5, sticky="w")
    label_ltb.grid(row=31, column=0, padx=10, pady=5, sticky="e")
    entry_ltb.grid(row=31, column=1, padx=10, pady=5, sticky="w")
    label_lbb.grid(row=32, column=0, padx=10, pady=5, sticky="e")
    entry_lbb.grid(row=32, column=1, padx=10, pady=5, sticky="w")
    label_n1.grid(row=33, column=0, padx=10, pady=5, sticky="e")
    entry_n1.grid(row=33, column=1, padx=10, pady=5, sticky="w")
    label_d.grid(row=34, column=0, padx=10, pady=5, sticky="e")
    entry_d.grid(row=34, column=1, padx=10, pady=5, sticky="w")
    label_do.grid(row=35, column=0, padx=10, pady=5, sticky="e")
    entry_do.grid(row=35, column=1, padx=10, pady=5, sticky="w")
    label_ds.grid(row=36, column=0, padx=10, pady=5, sticky="e")
    entry_ds.grid(row=36, column=1, padx=10, pady=5, sticky="w")
    label_lb.grid(row=37, column=0, padx=10, pady=5, sticky="e")
    entry_lb.grid(row=37, column=1, padx=10, pady=5, sticky="w")
    label_pt.grid(row=38, column=0, padx=10, pady=5, sticky="e")
    entry_pt.grid(row=38, column=1, padx=10, pady=5, sticky="w")
    label_rds.grid(row=39, column=0, padx=10, pady=5, sticky="e")
    entry_rds.grid(row=39, column=1, padx=10, pady=5, sticky="w")
    label_rdt.grid(row=40, column=0, padx=10, pady=5, sticky="e")
    entry_rdt.grid(row=40, column=1, padx=10, pady=5, sticky="w")
    label_bpercent.grid(row=41, column=0, padx=10, pady=5, sticky="e")
    entry_bpercent.grid(row=41, column=1, padx=10, pady=5, sticky="w")
    label_layout_type.grid(row=42, column=0, padx=10, pady=5, sticky="e")
    entry_layout_type.grid(row=42, column=1, padx=10, pady=5, sticky="w")
    label_lbout.grid(row=43, column=0, padx=10, pady=5, sticky="e")
    entry_lbout.grid(row=43, column=1, padx=10, pady=5, sticky="w")
    label_lbin.grid(row=44, column=0, padx=10, pady=5, sticky="e")
    entry_lbin.grid(row=44, column=1, padx=10, pady=5, sticky="w")
    label_nss.grid(row=45, column=0, padx=10, pady=5, sticky="e")
    entry_nss.grid(row=45, column=1, padx=10, pady=5, sticky="w")
    label_pitch.grid(row=46, column=0, padx=10, pady=5, sticky="e")
    entry_pitch.grid(row=46, column=1, padx=10, pady=5, sticky="w")
    label_l5.grid(row=47, column=0, padx=10, pady=5, sticky="e")
    entry_l5.grid(row=47, column=1, padx=10, pady=5, sticky="w")
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

label_p_allowed = Label(frame, text="enter the allowable pressure drop shell side or annulus(pascal):")
set_style(label_p_allowed)
label_p_allowed.grid(row=11, column=0, padx=10, pady=5, sticky="e")
entry_p_allowed = Entry(frame)
set_style(entry_p_allowed)
entry_p_allowed.grid(row=11, column=1, padx=10, pady=5, sticky="w")

label_p1_allowed = Label(frame, text="enter the allowed pressure drop tube side or pipe(pascal):")
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

label_cpc = Label(frame, text="specific heat capacity of cold fluid(J/Kg*K):")
set_style(label_cpc)
entry_cpc = Entry(frame)
set_style(entry_cpc)


label_cph = Label(frame, text="specific heat capacity of hot fluid(J/Kg*K):")
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
# for double pipe heat exchanger
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

label_l = Label(frame, text="enter the lenth of one hairpin(m):")
set_style(label_l)
entry_l = Entry(frame)
set_style(entry_l)

# for shell and tube heat exchanger
label_M = Label(frame, text="Number of shell passes:")
set_style(label_M)
entry_M = Entry(frame)
set_style(entry_M)

label_question = Label(frame, text="enter if the hot should go in 'tube' or shell:")
set_style(label_question)
entry_question = Entry(frame)
set_style(entry_question)

label_Np = Label(frame, text="enter the number of pass divider lanes:")
set_style(label_Np)
entry_Np = Entry(frame)
set_style(entry_Np)

label_Wp = Label(frame, text="enter the width of bypass lane(m):")
set_style(label_Wp)
entry_Wp = Entry(frame)
set_style(entry_Wp)

label_ntp = Label(frame, text="Number of tube per pass:")
set_style(label_ntp)
entry_ntp = Entry(frame)
set_style(entry_ntp)

label_lsb = Label(frame, text="enter the shell to baffle clearance(m):")
set_style(label_lsb)
entry_lsb = Entry(frame)
set_style(entry_lsb)

label_ltb = Label(frame, text="enter the tube outside diameter to baffle hole clearance(m):")
set_style(label_ltb)
entry_ltb = Entry(frame)
set_style(entry_ltb)

label_lbb = Label(frame, text="enter the shell diameter to bundle bypass clearance(m):")
set_style(label_lbb)
entry_lbb = Entry(frame)
set_style(entry_lbb)

label_n1 = Label(frame, text="number of tube passes:")
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


label_lb = Label(frame, text="baffle spacing(m):")
set_style(label_lb)
entry_lb = Entry(frame)
set_style(entry_lb)

label_pt = Label(frame, text="tube pitch(m):")
set_style(label_pt)
entry_pt = Entry(frame)
set_style(entry_pt)

label_rds = Label(frame, text="shell side dirt factor((K*(m^2))/W):")
set_style(label_rds)
entry_rds = Entry(frame)
set_style(entry_rds)

label_rdt = Label(frame, text="tube side dirt factor((K*(m^2))/W):")
set_style(label_rdt)
entry_rdt = Entry(frame)
set_style(entry_rdt)

label_bpercent = Label(frame, text="enter the percentage of baffle cut:")
set_style(label_bpercent)
entry_bpercent = Entry(frame)
set_style(entry_bpercent)

label_layout_type = Label(frame, text="enter the layout type in degrees:")
set_style(label_layout_type)
entry_layout_type = Entry(frame)
set_style(entry_layout_type)

label_lbout = Label(frame, text="enter the baffle spacing in outlet(m):")
set_style(label_lbout)
entry_lbout = Entry(frame)
set_style(entry_lbout)

label_lbin = Label(frame, text="enter the baffle spacing in inlet(m):")
set_style(label_lbin)
entry_lbin = Entry(frame)
set_style(entry_lbin)

label_nss = Label(frame, text="number of sealing strip:")
set_style(label_nss)
entry_nss = Entry(frame)
set_style(entry_nss)

label_pitch = Label(frame, text="what is the pitch type:")
set_style(label_pitch)
entry_pitch = Entry(frame)
set_style(entry_pitch)

label_l5 = Label(frame, text="length of tube per pass:")
set_style(label_l5)
entry_l5 = Entry(frame)
set_style(entry_l5)


def on_frame_configure(event):
    canvas.configure(scrollregion=canvas.bbox("all"))

frame.bind("<Configure>", on_frame_configure)


calculate_button = Button(frame, text="Calculate", command=solve_equations, bg="#4CAF50", fg="white", font=("Arial", 14))
calculate_button.grid(row=48, column=0, columnspan=2, pady=20)
   

window.grid_rowconfigure(0, weight=1)
window.grid_columnconfigure(0, weight=1)

result_label1 = Label(frame, text="", bg="#F0F0F0", font=("Arial", 14))


result_label2 = Label(frame, text="", bg="#F0F0F0", font=("Arial", 14))


result_label3 = Label(frame, text="", bg="#F0F0F0", font=("Arial", 14))


result_label4 = Label(frame, text="", bg="#F0F0F0", font=("Arial", 14))


result_label5 = Label(frame, text="", bg="#F0F0F0", font=("Arial", 14))



window.mainloop()

end_time = time.time()
runtime = end_time-start_time
print(f"runtime is {runtime}")