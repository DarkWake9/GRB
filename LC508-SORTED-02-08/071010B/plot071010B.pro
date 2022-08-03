;****************************************************************
;ANALISI OTTICO - X: 
;                        STEP 6 (Dati Ottici)
;****************************************************************
;------------------------------------------------------------------------------
;OAB Merate, 20 GENNAIO 2012
;E.Z. 
;------------------------------------------------------------------------------
;PLOT dei dati ottici e X con i relativi fit, sia divisi per filtri
;che per la LC combinata
;Input file: GRB+'_LC/'+GRB+'_bands.dat' con le colonne con tmin,
;tmax, tmin_tmax
;Output file: grb+'_LC'+grb+'_OPX2.eps'
;------------------------------------------------------------------------------
;FIT FUNCTIONS
;------------------------------------------------------------------------------


function bro1,x,P
         common ciao, t90
         ;Single Broken Power Law
         ;P[0]=tb
         ;P[1]=a
         ;P[2]=b
         ;P[3]=s
         ;P[4]=f0
         Y=P[4]*(((x-t90)/P[0])^(-P[1]/P[3])+((x-t90)/P[0])^(-P[2]/P[3]))^(P[3])          
         return,Y
end
;..............................................................................
function bro2,x,P
         common ciao, t90
         ;Double Broken Power Law
         ;P[0],P[5]=tb
         ;P[1],P[6]=a
         ;P[2],P[7]=b
         ;P[3],P[8]=s
         ;P[4],P[8]=f0
         Y1=P[4]*(((x-t90)/P[0])^(-P[1]/P[3])+((x-t90)/P[0])^(-P[2]/P[3]))^(P[3])   
         Y2=P[9]*(((x-t90)/P[5])^(-P[6]/P[8])+((x-t90)/P[5])^(-P[7]/P[8]))^(P[8])    
         Y=Y1+Y2        
         return,Y
end
;..............................................................................
function bro3,x,P
         common ciao, t90
         ;Rise + Double Broken Power Law 
         ;P[0],P[5],P[10]=tb
         ;P[1],P[6]=a
         ;P[2],P[7],P[11]=b
         ;P[3],P[8],P[12]=s
         ;P[4],P[9]=f0 
         Y1=P[4]*(((x-t90)/P[0])^(-P[1]/P[3])+((x-t90)/P[0])^(-P[2]/P[3]))^(P[3])   
         Y2=P[9]*(((x-t90)/P[5])^(-P[6]/P[8])+((x-t90)/P[5])^(-P[7]/P[8]))^(P[8])  
         Y3=(((x-t90)/P[10])^(-P[11]/P[12])+1.d)^(P[12])  
         Y=Y1+(Y2*Y3)        
         return,Y
      end
;..............................................................................
function bro4,x,P
         common ciao, t90
         ;P[0]=tb
         ;P[1]=a
         ;P[2]=b
         ;P[3]=s
         ;P[4],P[5]=f0
         ;P[6]=alpha 
         Y1=P[4]*(((x-t90)/P[0])^(-P[1]/P[3])+((x-t90)/P[0])^(-P[2]/P[3]))^(P[3])   
         Y2=P[5]*(x-t90)^(-P[6])
         Y=Y1+Y2     
         return,Y
      end
;..............................................................................
function bro5,x,P
         common ciao, t90
         ;P[0],P[5]=tb
         ;P[1],P[6],P[11]=a
         ;P[2],P[7]=b
         ;P[3],P[8]=s
         ;P[4],P[9],P[10]=f0 
         Y1=P[4]*(((x-t90)/P[0])^(-P[1]/P[3])+((x-t90)/P[0])^(-P[2]/P[3]))^(P[3])  
         Y2=(((x-t90)/P[5])^(-P[6]/P[7])+1.d)^(P[7])  
         Y3=P[8]*(x-t90)^(-P[9])
         Y=Y1*Y2+Y3     
         return,Y
end
;..............................................................................
function plaw,x,P 
         common ciao, t90
         ;P[0]=alpha
         ;P[1]=f0
         return,P[0]*(x-t90)^(-P[1])
end
;..............................................................................
function bro6,x,P
         common ciao, t90
         Y1=P[4]*(((x-t90)/P[0])^(-P[1]/P[3])+((x-t90)/P[0])^(-P[2]/P[3]))^(P[3])  
         Y2=(((x-t90)/P[5])^(-P[6]/P[7])+1.d)^(P[7])
         return,Y1*Y2
end  
;..............................................................................
function bro7,x,P
         common ciao, t90
         return,(((x-t90)/P[0])^(-P[1]/P[2])+1.d)^(P[2])
end  
;..........................................................................................................
function bro8,x,P
         common ciao, t90
         ;Double Broken Power Law
         ;P[0],P[5]=tb
         ;P[1],P[6]=a
         ;P[2],P[7]=b
         ;P[3],P[8]=s
         ;P[4],P[8]=f0
         Y1=P[4]*(((x-t90)/P[0])^(-P[1]/P[3])+((x-t90)/P[0])^(-P[2]/P[3]))^(P[3])   
         Y2=P[9]*(((x-t90)/P[5])^(-P[6]/P[8])+((x-t90)/P[5])^(-P[7]/P[8]))^(P[8])
         Y3=P[14]*(((x-t90)/P[10])^(-P[11]/P[13])+((x-t90)/P[10])^(-P[12]/P[13]))^(P[13])     
         Y=Y1+Y2+Y3        
         return,Y
      end
;..........................................................................................................
function bro9,x,P
         common ciao, t90
         ;Double Broken Power Law
         ;P[0],P[5]=tb
         ;P[1],P[6]=a
         ;P[2],P[7]=b
         ;P[3],P[8]=s
         ;P[4],P[8]=f0
         Y1=P[4]*(((x-t90)/P[0])^(-P[1]/P[3])+((x-t90)/P[0])^(-P[2]/P[3]))^(P[3])   
         Y2=P[9]*(((x-t90)/P[5])^(-P[6]/P[8])+((x-t90)/P[5])^(-P[7]/P[8]))^(P[8])
         Y3=P[10]*(x-t90)^(-P[11])     
         Y=Y1+Y2+Y3        
         return,Y
      end
;..........................................................................................................
;..........................................................................................................
;..........................................................................................................
pro plot071010B
common ciao, t90



;---------------
;minx=30.
;maxx=5.e7
;miny=2.e-8
;maxy=3.e-2
;---------------

t90=0.
;Stupid thinks...
basic_colors,black,white,red,green,blue,yellow,cyan,magenta,orange,mint,purple,pink,olive,lightblue,gray
cm2=teXtoIDL("cm^{2}")
flusso=teXtoIDL("Flux (erg cm^{-2} s^{-1} A^{-1})")
flusso10=teXtoIDL("^{10} erg cm^{-2} s^{-1} A^{-1}")
plotsym,0,0.5,/fill
;The program asks the GRB name, on shell
nome=' '
read,nome,prompt='GRB name: GRB',format='(a)'
GRB=nome

readcol,GRB+'_LC/'+GRB+'_plotrange.dat',minx,maxx,miny,maxy,format='d,d,d,d'


pp=[0.1,0.14,0.99,0.98]
pp1=[0.14,0.50,0.97,0.98]
pp2=[0.14,0.14,0.97,0.50]


readcol,GRB+'_LC/comb/'+GRB+'_filterOK.dat',filter,format='a'
readcol,GRB+'_LC/comb/'+GRB+'_comb_flux.dat',time,err_time,flux,err_flux,format='d,d,d,d'
informations=GRB+'_LC/'+GRB+'_comb.par'
guessdata=GRB+'_LC/'+GRB+'_comb_guess.dat'
finalfit=GRB+'_LC/'+GRB+'_comb_result.dat'
openr,4,informations
nihil=' '
readf,4,nihil,nihil,nihil,z,nihil,dl,nihil,fit_func
close,4
multiplot,[1,2]      
set_plot,'ps'
device,filename=GRB+'_LC/'+GRB+'_OPX2.eps',/COLOR,/encapsulated
;plot,[time[0],time[1]],[flux[0]*10000.,flux[1]*10000.],/xlog,/ylog,psym=3,ytitle='Flux (Jy)',charsize=1.,/xst,/yst,xthick=4.,ythick=4.,xr=[1.d,1.d8],yr=[min(flux)*0.1,max(flux)*100000.],charthick=4.,pos=pp1

plot,[time[0],time[1]],[1.e-6,1.],/xlog,/ylog,psym=3,ytitle='Flux (Jy)',charsize=1.,/xst,/yst,xthick=4.,ythick=4.,xr=[minx,maxx],yr=[miny,maxy],charthick=4.,pos=pp1

;Bande colorate per indicare dove fare le SED

bande=' '
read,bande,prompt='Colored bands (SED) (no=0, yes=1=>GRB_bande.dat [x1,x2])?',format='(a)'
if bande eq '0' then begin
   bb=-1. 
   goto,dopoo
endif
collo=[yellow,pink,orange,green,lightblue,cyan,magenta]
readcol,GRB+'_LC/'+GRB+'_bands.dat',x1,x2,format='d,d'
print,x1
print,x2
bb=n_elements(x1)
for p=0,bb-1 do begin
   POLYFILL, [x1[p], x1[p], x2[p], x2[p]], [miny,maxy,maxy,miny], COL = collo[p],/LINE_FILL,orientation=45,spacing=0.1  
endfor
dopoo:

;; oploterror,time,flux*10000.,err_time,err_flux*10000.,psym=8,color=red,errcolor=red,thick=4.
;; if (fit_func eq 0) then begin
;;    readcol,guessdata,aa,norm,format='d,d,d,d,d'
;;    guess=[aa,norm] 
;;    FF=plaw(time,guess)*1.d-10
;;    oplot,time,FF*10000.,linestyle=0,thick=4.,color=black
;; endif
;; if (fit_func eq 1) then begin
;;    readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
;;    guess=[tb,a1,a2,s,norm] 
;;    FF=bro1(time,guess)*1.d-10
;;    oplot,time,FF*10000.,linestyle=0,thick=4.,color=black
;; endif
;; if (fit_func eq 2) then begin
;;    readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
;;    guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1],a2[1],s[1],norm[1]] 
;;    FF=bro2(time,guess)*1.d-10
;;    oplot,time,FF*10000.,linestyle=0,thick=4.,color=black
;; endif
;; if (fit_func eq 3) then begin
;;    readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
;;    guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1],a2[1],s[1],norm[1],tb[2],a1[2],a2[2]] 
;;    FF=bro3(time,guess)*1.d-10
;;    oplot,time,FF*10000.,linestyle=0,thick=4.,color=black
;; endif
;; if (fit_func eq 4 or fit_func eq 6) then begin
;;    readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
;;    guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1]] 
;;    FF=bro4(time,guess)*1.d-10
;;    oplot,time,FF*10000.,linestyle=0,thick=4.,color=black
;; endif
;; if (fit_func eq 5) then begin
;;    readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
;;    guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1],a2[1],tb[2],a1[2]] 
;;    FF=bro1(time,guess)*1.d-10
;;    oplot,time,FF*10000.,linestyle=0,thick=4.,color=black
;; endif
colore=[black,green,blue,yellow,cyan,magenta,orange,mint,purple,pink,olive,lightblue,gray]
g=0
kk=1
;filtro=strarr(n_elements(filter)+2)
filtro=strarr(n_elements(filter)+1)
;filtro[1]='Combined'
filtro[0]='GRB '+GRB
;col=strarr(n_elements(filter)+2)
col=strarr(n_elements(filter)+1)
col[0]=white
;col[1]=red
;piss=strarr(n_elements(filter)+2)
piss=strarr(n_elements(filter)+1)
piss[0]=0
;piss[1]=8
for k=0,n_elements(filter)-1 do begin
    readcol,GRB+'_LC/comb/'+GRB+'_'+filter[k]+'_flux.dat',time,err_time,flux,err_flux,format='d,d,d,d'
    informations=GRB+'_LC/'+GRB+'_'+filter[k]+'.par'
    guessdata=GRB+'_LC/'+GRB+'_'+filter[k]+'_guess.dat'
    finalfit=GRB+'_LC/'+GRB+'_'+filter[k]+'_result.dat'
    openr,4,informations
    nihil=' '
    readf,4,nihil,nihil,nihil,z,nihil,dl,nihil,fit_func
    close,4

    oploterror,time,flux,err_time,err_flux,psym=kk,color=colore[g],errcolor=colore[g],thick=4.

    if (fit_func eq 0) then begin
       readcol,guessdata,aa,norm,format='d,d,d,d,d'
       guess=[aa,norm] 
       FF=plaw(time,guess)
       oplot,time,FF*1.d-10,linestyle=0,thick=4.,color=black
    endif
    if (fit_func eq 1) then begin
       readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
       guess=[tb,a1,a2,s,norm*1.d-10] 
       FF=bro1(time,guess)
       oplot,time,FF,linestyle=0,thick=4.,color=black
    endif
    if (fit_func eq 2) then begin
       readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
       guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1],a2[1],s[1],norm[1]] 
       FF=bro2(time,guess)*1.d-10
       oplot,time,FF,linestyle=0,thick=4.,color=black
    endif
    if (fit_func eq 3) then begin
       readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
       guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1],a2[1],s[1],norm[1],tb[2],a1[2],a2[2]] 
       FF=bro3(time,guess)*1.d-10
       oplot,time,FF,linestyle=0,thick=4.,color=black
    endif
    if (fit_func eq 4 or fit_func eq 6) then begin
       readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
       guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1]] 
       FF=bro4(time,guess)*1.d-10
       oplot,time,FF,linestyle=0,thick=4.,color=black
    endif
    if (fit_func eq 5) then begin
       readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
       guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1],a2[1],tb[2],a1[2]] 
       FF=bro1(time,guess)*1.d-10
       oplot,time,FF,linestyle=0,thick=4.,color=black
    endif
    if (fit_func eq 8) then begin
       readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
       guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1],a2[1],s[1],norm[1],tb[2],a1[2],a2[2],s[2],norm[2]] 
       FF=bro8(time,guess)*1.d-10
       oplot,time,FF,linestyle=0,thick=4.,color=black
    endif
    if (fit_func eq 9) then begin
       readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
       guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1],a2[1],s[1],norm[1],tb[2],a1[2]] 
       FF=bro9(time,guess)*1.d-10
       oplot,time,FF,linestyle=0,thick=4.,color=black
    endif
    filtro[k+1]=filter[k]
    col[k+1]=colore[g]
    piss[k+1]=kk
    g=g+1.
    kk=kk+1.
    if kk eq 3 then kk=4
    if kk eq 9 then kk=0
    if g eq 13 then g=0
 endfor
legend,filtro,color=col,/right,box=0,charthick=4.,psym=piss
 
 

;####################################################################
;####################################################################
;####################################################################
; XRT DATA
plotsym,0,0.5,/fill
;readcol,'/Users/elena/Desktop/XLC_CATALOG/archivio/grb'+GRB+'/lightcurve/lc_SN/flux/flux_NH_fixed/SN'+GRB+'_4_nhtot_obs_flux.txt',timef,err_timef,flux,err_flux_l,err_flux_r,format='d,d,d,d,d',skipline=2
assex=textoidl("Obs Time since BAT trigger (s)")
assey=textoidl("Flux (erg s^{-1} cm^{-2})")

if GRB eq '080913A' then begin
   GRB='080913'
endif
readcol,'/Users/elena/Desktop/ANALISIOTTICA2011/ANALISI/dati_ottici_FASE1/GRB2011.dat',grb11,format='a'
kk=where(GRB eq grb11,kkkk)
if kkkk eq 1 then goto,phildata
readcol,'/Users/elena/Desktop/XLC_CATALOG/TABELLA/TABLE_FIT_120417.dat',name,type,zz,tmin,tmax,a1,ea1,a2,ea2,$
a3,ea3,a4,ea4,a5,ea5,tb1,etb1,tb2,etb2,$
tb3,etb3,tb4,etb4,norm1,enorm1,norm2,enorm2,norm3,enorm3,$
s1,es1,s2,es2,s3,es3,chi2,dof,pval,$
format='a,a,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d'

;readcol,'/Users/elena/Desktop/XLC_CATALOG/TABELLA/GRB'+GRB+'/'+GRB+'/GRB'+GRB+'_continuum.dat',timef,err_time,flux,err_flux,format='d,d,d,d,d'
;readcol,'/Users/elena/Desktop/XLC_CATALOG/TABELLA/GRB'+GRB+'/'+GRB+'/GRB'+GRB+'_excess.dat',te,ete,fe,efe,format='d,d,d,d,d'
readcol,GRB+'_LC/SN071010B_4_complete.txt',timef,err_time,flux,err_flux,dur,sn,format='d,d,d,d,d,d,d',skipline=2
flux=flux*5.8865e-11
err_flux=err_flux*5.8865e-11
;-----------
tn=fltarr(15000)
tn[0]=0.1
for j=1,14999 do begin
    tn[j]=tn[j-1]*1.01
endfor
yyy=where(tn ge min(timef) and tn le max(timef),yy)
ttt=dblarr(yy)
for j=0,yy-1 do begin
    ttt[j]=tn[yyy[j]]
endfor
;-----------

grbok=where(GRB eq name)
;POWER-LAW
if (type[grbok] eq '0UN' or type[grbok] eq '0UF' or type[grbok] eq '0CN' or type[grbok] eq '0CF') then begin
   guess=[norm1[grbok],a1[grbok]]
   fitfunc=plaw(ttt,guess)
   t1=0.
   t2=0.
   t3=0.
endif
;BRO1
if (type[grbok] eq '1UN' or type[grbok] eq '1UF' or type[grbok] eq '1CN' or type[grbok] eq '1CF') then begin
   guess=[tb1[grbok],a1[grbok],a2[grbok],s1[grbok],norm1[grbok]]
   fitfunc=bro1(ttt,guess)
   t1=tb1[grbok]
   t2=0.
   t3=0.
endif
;BRO2
if (type[grbok] eq '2UN' or type[grbok] eq '2UF' or type[grbok] eq '2CN' or type[grbok] eq '2CF') then begin
   guess=[tb1[grbok],a1[grbok],a2[grbok],s1[grbok],norm1[grbok],tb3[grbok],a3[grbok],a4[grbok],s2[grbok],norm2[grbok]]
   fitfunc=bro2(ttt,guess)
   t1=tb1[grbok]
   t2=tb2[grbok]
   t3=tb3[grbok]
endif
;BRO4
if (type[grbok] eq '4UN' or type[grbok] eq '4UF' or type[grbok] eq '4CN' or type[grbok] eq '4CF') then begin
   guess=[tb2[grbok],a1[grbok],a2[grbok],s1[grbok],norm1[grbok],norm2[grbok],a3[grbok]]
   fitfunc=bro4(ttt,guess)   
   t1=tb1[grbok]
   t2=tb2[grbok]
   t3=0.
endif
;BRO6
if (type[grbok] eq '6UN' or type[grbok] eq '6UF' or type[grbok] eq '6CN' or type[grbok] eq '6CF') then begin
   guess=[tb1[grbok],a1[grbok],a2[grbok],s1[grbok],norm1[grbok],norm2[grbok],a3[grbok]]
   fitfunc=bro4(ttt,guess)
   t1=tb1[grbok]
   t2=tb2[grbok]
   t3=0.  
endif

goto,dopo

phildata:
readcol,'XRAY_EVANS/'+GRB+'_evans.dat',timef,ETXP,ETXM,flux,err_flux,format='d,d,d,d,d'
err_time=dblarr(n_elements(timef))
err_time=(ETXP+ETXM)/2.
flux=flux*1.d10
err_flux=err_flux*1.d10
dopo:
multiplot
minf=min(flux*0.1*1.d-10)
maxf=max(flux*10.*1.d-10)
if GRB eq '050820A' then maxf=1.d-6
if GRB eq '060124' then maxf=1.d-6
if GRB eq '060526' then maxf=1.d-6
if GRB eq '060904B' then maxf=1.d-7
plot,timef,flux,/xlog,/ylog,psym=8,xtitle=assex,ytitle=assey,charsize=1,/xst,/yst,xthick=4.,ythick=4.,xr=[minx,maxx],charthick=4.,pos=pp2,yr=[minf,maxf]

if (bb gt 0.) then begin
   for p=0,bb-1 do begin
      POLYFILL, [x1[p], x1[p], x2[p], x2[p]], [minf,maxf, maxf,minf], COL = collo[p],/LINE_FILL,orientation=45,spacing=0.1  
   endfor
endif
oploterror,timef,flux*1.d-10,err_time,err_flux*1.d-10,psym=8,color=black,thick=3.,errcolor=black
if (kkkk ne 1) then begin 
; if (n_elements(te) le 2) then begin
;    oploterror,te,fe*1.d-10,ete,efe*1.d-10,psym=8,color=black,thick=3.,errcolor=black
; endif else begin
;    oploterror,te,fe*1.d-10,ete,efe*1.d-10,psym=8,color=red,thick=3.,errcolor=red
; endelse

;HANDREMOVED
; readcol,'/Users/elena/Desktop/XLC_CATALOG/TABELLA/GRB'+GRB+'/'+GRB+'/GRB'+GRB+'_hand_removed.dat',th,eth,fh,efh1,efh2,format='d,d,d,d,d'
; for k=0,n_elements(th)-1 do begin 
;     th[k]=th[k]*(zz[grbok]+1.d)
;     eth[k]=eth[k]*(zz[grbok]+1.d)
;  endfor
; oploterror,th,fh*1.d-10,eth,(efh1*1.d-10+efh2*1.d-10)/2.,psym=8,color=red,thick=3.,errcolor=red
;++++++++++++++++++++++++
oplot,ttt,fitfunc*1.d-10,linestyle=0,color=blue,thick=3.
oplot,[t1,t1],[1.d-30,1.],linestyle=2,color=gray,thick=3.
oplot,[t2,t2],[1.d-30,1.],linestyle=2,color=gray,thick=3.
oplot,[t3,t3],[1.d-30,1.],linestyle=2,color=gray,thick=3.
endif
multiplot,/reset    
device,/close
stop
end
