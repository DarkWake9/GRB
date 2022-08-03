;****************************************************************
;ANALISI OTTICO - X: 
;                        STEP 5 (Dati Ottici)
;****************************************************************
;------------------------------------------------------------------------------
;OAB Merate,17 GENNAIO 2012
;E.Z. 
;------------------------------------------------------------------------------
;OPTICAL LIGHT CURVES ANALYSIS
;Program that fits the OPTICAL light curves (the original program is 'xcurve.pro').
;/LC_FIT/: folder where there is this program and where to launch IDL
;/LC_FIT/#_LC/: folder where to put the input files and where the
;program puts the output files.
;
;In the file name _filter or _comb!!!
;
;Input files:
;1) GRB+'_LC/comb/'+GRB+'_'+filter+'_combLC.dat': the data (comb.pro)
;2) GRB+'_LC/'+GRB+'_'+filter+'.par': where I put some useful
;informations for the fit (vd README)
;3) GRB+'_LC/'+GRB+'_'+filter+'_guess.dat':Initial guess
;
;Output files:
;1) GRB+'_LC/comb/'+GRB+'_'+filter+'_flux.dat': the data in flux units
;2) GRB+'_LC/'+GRB+'_'+filter+'_result.dat': The results (parameters)
;3) GRB+'_LC/'+GRB+'_'+filter+'.eps': Plot with the data and the fit function
;4) GRB+'_LC/'+GRB+'_'+filter+'_cova.dat': covariance matrix
;5) GRB+'_LC/'+GRB+'_'+filter+'_cova_useful.dat': covariance matrix
;written in a useful way
;
; + It is possible to erase some points:
;   write # befor the data that you want eliminate in the data file
;   (_combLC.dat) and copy these line in the file 
;   GRB+'_LC/comb/'+GRB+'_'+filter+'_scartati.dat'
;   In this case will be another output file, with the eliminated
;   points calibrated in flux:
;   GRB+'_LC/comb/'+GRB+'_'+filter+'_scartatiFLUX.dat'
;
;Functions for the fit: 
;     0)'plaw': power law
;     1)'bro1': single broken power law 
;     2)'bro2': double broken power law
;     3)'bro3': rise + double broken power law
;     4)'bro4': plaw + broken power law
;     5)'bro5': plaw + double broken power law
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
;##########################################################################################################
pro lcuno

  
;##########################################################################################################
;##########################################################################################################
;##########################################################################################################
;THE BEGIN
;##########################################################################################################
;##########################################################################################################
;##########################################################################################################
common ciao, t90
t90=0.
;Stupid thinks...
basic_colors,black,white,red,green,blue,yellow,cyan,magenta,orange,mint,purple,pink,olive,lightblue,gray
cm2=teXtoIDL("cm^{2}")
flusso=teXtoIDL("erg cm^{-2} s^{-1} !A^{-1}")
flusso10=teXtoIDL("^{10} erg cm^{-2} s^{-1} !A^{-1}")
plotsym,0,0.8,/fill

;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
;FILES
;Complete or single filter LC?
;INPUTFILE,PARFILE...
TEMPOOK=700.
FILTROO='RC'
GRB='080607'
nome=GRB
inputfile='../comb/'+GRB+'_'+FILTROO+'_flux.dat'
informations=GRB+'_'+FILTROO+'.par' ;Initial guess
guessdata=GRB+'_'+FILTROO+'_guess.dat' ;Initial guess
finalfit=GRB+'_'+FILTROO+'_result.dat' ;The results
epsfile=GRB+'_'+FILTROO+'.eps' ;Plot with the data and the fit function
covafile=GRB+'_'+FILTROO+'b_cova.dat'
covafile_useful=GRB+'_'+FILTROO+'_cova_useful.dat'

readcol,inputfile,t,et,m,em,format='d,d,d'
lung=n_elements(t)

;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
for j=0,lung-2 do begin
    for l=j+1,lung-1 do begin
        if t[l] lt t[j] then begin 
           temp=t[j]
           temp1=et[j]
           temp2=m[j]
           temp3=em[j]
           t[j]=t[l]
           et[j]=et[l]
           m[j]=m[l]
           em[j]=em[l]
           t[l]=temp
           et[l]=temp1
           m[l]=temp2
           em[l]=temp3             
        endif
    endfor
endfor

;On screen plot
bene=where(t gt TEMPOOK)
time=t[bene]
err_time=et[bene]
flux=m[bene]*1.d10
err_flux=em[bene]*1.d10

set_plot,'x'
window,0,retain=1
plot,time,flux,/xlog,/ylog,psym=6,title='GRB'+GRB,xtitle='Obs Time since BAT trigger',ytitle='10'+flusso10,xrange=[min(time),max(time)],yrange=[min(flux)*0.1,max(flux)*10.]
oploterror,time,flux,err_time,err_flux,psym=6

print,' '
print,' '
print,'*****************************************************'
print,'UPDATE THE FOLLOWING FILES:'
print,'#_comb.par or #_filter_comb.par' 
print,'(see the example, GRB_***.par, it depends on the function type)'
print,'#_comb_guess.dat #_filter_guess.dat(read README.txt)'
print,' '
print,'            THEN...                '
print,' '
print,'TO CONTINUE WRITE: .cont'
print,' '
print,' '
print,'*****************************************************'
stop

;------------------------------------------------------------------------------------------------

;Read the GRB#.par
openr,9,informations
nnnn=' '
readf,9,nnnn,nnnn,nnnn,z,nnnn,dlo,nnnn,fit_func
close,9

if (fit_func eq 0) then goto, fit_function_0  ;Power Law
if (fit_func eq 1) then goto, fit_function_1  ;1BRO
if (fit_func eq 2) then goto, fit_function_2  ;2BRO
if (fit_func eq 3) then goto, fit_function_3  ;3BRO
if (fit_func eq 4 or fit_func eq 6) then goto, fit_function_4  ;4BRO
if (fit_func eq 5) then goto, fit_function_5  ;5BRO
if (fit_func eq 8) then goto, fit_function_8  ;8BRO
if (fit_func eq 9) then goto, fit_function_9  ;9BRO
;##########################################################################################################
;##########################################################################################################
;##########################################################################################################
;THE FIT
;##########################################################################################################
;##########################################################################################################
;##########################################################################################################

;##########################################################################################################
;###################################         0              ###############################################
;##########################################################################################################
fit_function_0:
readcol,guessdata,aa,norm,format='d,d,d,d,d'
guess=[aa,norm] 
FF=plaw(time,guess)
window,0,retain=1
plot,time-t90,flux,/xlog,/ylog,psym=8,title='GRB'+nome,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.2,$
               charsize=1.3,thick=2.,xrange=[min(time),max(time)],yrange=[min(flux),max(flux)]
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,FF,linestyle=0,color=red,thick=3.

print,'####################################'
print,' '
print,' '
print,'WRITE .cont TO CONTINUE'
print,' '
print,' '
print,'####################################'
stop

;openw,4,excess
;printf,4,'#Time',' ','Err_time',' ','Flux',' ','Err_Flux',' ','Residual',' ','Err_Residual',format='(a,a,a,a,a,a,a,a,a,a,a)'
;printf,4,'-1.',' ','-1.',' ','-1.',' ','-1.',' ','-1.',' ','-1.',format='(a,a,a,a,a,a,a,a,a,a,a)'
;-----------------------------------------------------------
;Read the GRB#.par
openr,2,informations
nihil=' '
readf,2,nihil,nihil,nihil,z,nihil,dlo,nihil,fit_func,nihil,nihil,zero,nihil,uno
close,2
;-----------------------------------------------------------
;If we want to fix or not the parameters...
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},2)
if (zero eq 1) then begin
   pi(0).fixed = 1
   guess(0) = guess[0]
endif   
if (uno eq 1) then begin
   pi(1).fixed = 1
   guess(1) = guess[1]
endif  
;-----------------------------------------------------------
for j=0,n_elements(flux) do begin
    para=dblarr(5)
    para=MPFITFUN('plaw',time,flux,err_flux,guess,best_fjac=best_fjac, pfree_index=pfree_index, /calc_fjac, $
                 covar=pcovar,DOF=dof,BESTNORM=chi,PERROR=errore,PARINFO=pi)
    
    cova0=pcovar[0,1]
    fit=fltarr(n_elements(time))
    fit=plaw(time,para)
    residual=flux-fit
    residual2=(flux-fit)/fit
    ;1 sigma error obtained with 'mpfit_properr.pro'
    ycovar = mpfit_properr(best_fjac, pcovar, pfree_index)
    sy_prop = sqrt(abs(ycovar))
    rapporto=residual/sy_prop
    err_res=fltarr(n_elements(fit))
    for oo=0,n_elements(fit)-1 do begin
        err_res[oo]=sqrt((err_flux[oo]/flux[oo])^(2.d)+(sy_prop[oo,oo]/fit[oo])^(2.d))*abs(residual2[oo])
    endfor
     guess=para
  plotsym,0,0.8,/fill
  fit=plaw(time,para)
  set_plot,'x'
  plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
  oploterror,time-t90,flux,err_time,err_flux,psym=8
  oplot,time-t90,fit,linestyle=0,color=red,thick=3.
  pvalue,n_elements(time),chi,dof,p_val
  if (p_val ge 0.3) then goto, fuori0
  if j eq 10 then goto, fuori0
endfor
fuori0:

close,4
openw,1,guessdata
printf,1,para[0],' ',para[1],format='(e,a,e)'
close,1
window,0,retain=1
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,fit,linestyle=0,color=red,thick=3.

print,'*******************************************'
print,'FIT PARAMETERS:'
print,''
print,'CHI^2',chi
print,'DOF',DOF
print,'p_value',p_val
print,'Alpha = ',para[1],'+/-',errore[1]
print,'Normalization = ',para[0],'+/-',errore[0],' 10^(-10) erg/s/cm^2/A'
print,'*******************************************'
tbmi0=-1.
tbmi1=-1.
tbmi2=-1.
tbmi3=-1.
openw,2,finalfit
printf,2,'Norm (10^-10 erg/s/cm2)',' ','Err_Norm',' ','Alpha',' ','Err_alpha',' ','CHI2',' ','DOF',' ','P_value',format='(a,a,a,a,a,a,a,a,a,a,a,a,a)'
printf,2,' '
printf,2,para[0],' ',errore[0],' ',para[1],' ',errore[1],' ',chi,' ',dof,' ',p_val,format='(e,a,e,a,e,a,e,a,e,a,e,a,e)'
close,2
goto,final_part
;##########################################################################################################
;###################################         1              ###############################################
;##########################################################################################################
fit_function_1:
readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
guess=[tb,a1,a2,s,norm] 
FF=bro1(time,guess)
window,0,retain=1
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,FF,linestyle=0,color=red,thick=3.
print,'####################################'
print,' '
print,' '
print,'WRITE .cont TO CONTINUE'
print,' '
print,' '
print,'####################################'
stop
;-----------------------------------------------------------
;Read the GRB#.par
openr,2,informations
nihil=' '
readf,2,nihil,nihil,nihil,z,nihil,dloo,nihil,fit_func,nihil,nihil,zero,nihil,uno,nihil,due,nihil,tre,nihil,quattro
close,2
;-----------------------------------------------------------
;If we want to fix or not the parameters...
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},5)
numero=[zero,uno,due,tre,quattro]
for k=0,n_elements(numero)-1 do begin
    if (numero[k] eq 1) then begin
       pi(k).fixed = 1
       guess(k) = guess[k]
    endif
endfor
;-----------------------------------------------------------
for j=0,1 do begin;n_elements(flux) do begin
    para=dblarr(5)
    para=MPFITFUN('bro1',time,flux,err_flux,guess,best_fjac=best_fjac, pfree_index=pfree_index, /calc_fjac, $
                 covar=pcovar,DOF=dof,BESTNORM=chi,PERROR=errore,PARINFO=pi)
    ;............................................
    kk=n_elements(para)
    nn=fltarr(kk)
    nn[0]=kk-1.
    for i=1,kk-2. do begin
        nn[i]=nn[i-1]-1.
    endfor
    gg=total(nn)
    cova1=dblarr(gg) 
    for w=0,gg-1 do begin
        for i=0,n_elements(para)-1 do begin
            for j=i+1,n_elements(para)-1 do begin
                cova1[w]=pcovar[i,j]
            endfor
        endfor
    endfor
    ;............................................
    Fit=fltarr(n_elements(time))
    fit=bro1(time,para)
    residual=flux-fit
    residual2=(flux-fit)/fit
    ;1 sigma error obtained with 'mpfit_properr.pro'
    ycovar = mpfit_properr(best_fjac, pcovar, pfree_index)
    sy_prop = sqrt(abs(ycovar))
    rapporto=residual/sy_prop
    err_res=fltarr(n_elements(fit))
    for oo=0,n_elements(fit)-1 do begin
        err_res[oo]=sqrt((err_flux[oo]/flux[oo])^(2.d)+(sy_prop[oo,oo]/fit[oo])^(2.d))*abs(residual2[oo])
    endfor
     guess=para
  plotsym,0,0.8,/fill
  fit=bro1(time,para)
  set_plot,'x'
  plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
  oploterror,time-t90,flux,err_time,err_flux,psym=8
  oplot,time-t90,fit,linestyle=0,color=red,thick=3.
  pvalue,n_elements(time-t90),chi,dof,p_val
  if (j eq 2) then goto, fuori1
endfor
fuori1:
window,0,retain=1
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,fit,linestyle=0,color=red

window,2,retain=1
norma=bro1(para[0],para)
par1=[max(flux)*2.,para[1]]
f1=plaw(time,par1)
par2=[max(flux)*2.d3,para[2]]
f2=plaw(time,par2)
par3=[para[4],-0.829]
f3=plaw(time,par3)
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,f1,linestyle=0,color=blue,thick=3.
oplot,time-t90,f2,linestyle=0,color=green,thick=3.
oplot,[para[0],para[0]],[1.d-10,1.d10],linestyle=2,color=red,thick=3.
oplot,time-t90,f3,linestyle=0,color=pink,thick=4.
openw,1,guessdata
printf,1,para[0],' ',para[1],' ',para[2],' ',para[3],' ',para[4],format='(e,a,e,a,e,a,e,a,e)'
close,1

print,'*******************************************'
print,'FIT PARAMETERS:'
print,''
print,'CHI^2',chi
print,'DOF',DOF
print,'p_value',p_val
print,'Break time = ',para[0],'+/-',errore[0] 
print,'Alpha 1 = ',para[1],'+/-',errore[1]
print,'Alpha 2 = ',para[2],'+/-',errore[2]
print,'Smoothness = ',para[3],'+/-',errore[3]
print,'Normalization = ',para[4],'+/-',errore[4],' 10^(-10) erg/s/cm^2/A'
print,'*******************************************'
tbmi0=para[0]
tbmi1=-1.
tbmi2=-1.
tbmi3=-1.
openw,2,finalfit
printf,2,'Break Time',' ','Err_BreakTime',' ','Alpha1',' ','Err_alpha1',' ','Alpha2',' ','Err_alpha2',' ','s',' ','Err_s',' ','Norm (10^48 erg/s/cm2)',' ','Err_Norm',' ','CHI2',' ','DOF',' ','P_value',format='(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
printf,2,' '
printf,2,para[0],' ',errore[0],' ',para[1],' ',errore[1],' ',para[2],' ',errore[2],' ',para[3],' ',errore[3],' ',para[4],' ',errore[4],' ',chi,' ',dof,' ',p_val,format='(e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e)'
close,2
goto,final_part
;##########################################################################################################
;###################################         2              ###############################################
;##########################################################################################################
fit_function_2:
readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1],a2[1],s[1],norm[1]] 
FF=bro2(time,guess)
set_plot,'x'
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,FF,linestyle=0,color=red,thick=3.
print,'####################################'
print,' '
print,' '
print,'WRITE .cont TO CONTINUE'
print,' '
print,' '
print,'####################################'
stop
;; openw,4,excess
;; printf,4,'Time',' ','Err_time',' ','Flux',' ','Err_Flux',' ','Residual',' ','Err_Residual',format='(a,a,a,a,a,a,a,a,a,a,a)'

;-----------------------------------------------------------
;Read the GRB#.par
openr,2,informations
nihil=' '
readf,2,nihil,nihil,nihil,z,nihil,dlo,nihil,fit_func,nihil,nihil,zero,nihil,uno,nihil,due,nihil,tre,nihil,quattro,nihil,cinque,nihil,sei,nihil,sette,nihil,otto,nihil,nove,nihil,tmin,nihil,tmax
close,2
;-----------------------------------------------------------
;If we want to fix or not the parameters...
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},10)
numero=[zero,uno,due,tre,quattro,cinque,sei,sette,otto,nove]
for k=0,n_elements(numero)-1 do begin
    if (numero[k] eq 1) then begin
       pi(k).fixed = 1
       guess(k) = guess[k]
    endif
endfor
;-----------------------------------------------------------
for j=0,n_elements(flux) do begin
    para=dblarr(10)
    para=MPFITFUN('bro2',time,flux,err_flux,guess,best_fjac=best_fjac, pfree_index=pfree_index, /calc_fjac, $
                 covar=pcovar,DOF=dof,BESTNORM=chi,PERROR=errore,PARINFO=pi)
    ;............................................
    kk=n_elements(para)
    nn=fltarr(kk)
    nn[0]=kk-1.
    for i=1,kk-2. do begin
        nn[i]=nn[i-1]-1.
    endfor
    gg=total(nn)
    cova2=dblarr(gg) 
    for w=0,gg-1 do begin
        for i=0,n_elements(para)-1 do begin
            for j=i+1,n_elements(para)-1 do begin
                cova2[w]=pcovar[i,j]
            endfor
        endfor
    endfor
    ;............................................
    fit=fltarr(n_elements(time))
    fit=bro2(time,para)
    residual=(flux-fit)
    residual2=(flux-fit)/fit
    ;1 sigma error obtained with 'mpfit_properr.pro'
    ycovar = mpfit_properr(best_fjac, pcovar, pfree_index)
    sy_prop = sqrt(abs(ycovar))
    rapporto=residual/sy_prop
    err_res=fltarr(n_elements(fit))
    for oo=0,n_elements(fit)-1 do begin
        err_res[oo]=sqrt((err_flux[oo]/flux[oo])^(2.d)+(sy_prop[oo,oo]/fit[oo])^(2.d))*abs(residual2[oo])
    endfor
     guess=para
  plotsym,0,0.8,/fill
  fit=bro2(time,para)
  set_plot,'x'
  plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
  oploterror,time-t90,flux,err_time,err_flux,psym=8
  oplot,time-t90,fit,linestyle=0,color=red,thick=3.
  pvalue,n_elements(time),chi,dof,p_val
  if (p_val ge 0.3) then goto, fuori2
  if (j ge 10.) then goto, fuori2
endfor
fuori2:
window,0,retain=1
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,fit,linestyle=0,color=red,thick=3.

tbmid=0.
goto,nuovo
risposta=' '
read,risposta,prompt='The fit is ok? NO=0, YES=1:       ',format='(a)'

if (risposta eq '0') then goto, stopfinale

;===================================================================================
pp0=[para[5],para[6],para[7],para[8],para[9]]
pp1=[para[0],para[1],para[2],para[3],para[4]]
fu0=bro1(time,pp0)
fu1=bro1(time,pp1)
wh0=where(fu0 le fu1)
wh1=where(fu0 ge fu1)
for j=0,n_elements(wh1)-1 do begin
if (wh1[j] gt [-1]) then tmax=time[min(wh1)]
endfor
for j=0,n_elements(wh0)-1 do begin
if (wh0[j] gt [-1]) then tmin=time[max(wh0)]
endfor
dt=10.
for i=0,1000 do begin
    num=floor((tmax-tmin)/dt)
    if (num le 1) then num=3.
    if (num gt 1000) then dt=dt*10000.
    tnew=dblarr(num)
    tnew[0]=tmin
    for j=1,num-1 do begin
       tnew[j]=tnew[j-1]+dt
    endfor
    fu0_bis=bro1(tnew,pp0)
    fu1_bis=bro1(tnew,pp1)
    
    rapporto=where(fu0_bis/fu1_bis ne 0,cc)
    rappo=dblarr(cc)
    for j=0,cc-1 do begin
        rappo[j]=fu0_bis[j]/fu1_bis[j]
    endfor
    rap=min(rappo)

    if (rap le rap/20.) then goto,vai2   ;5%
    dt=dt/2.   
endfor
vai2:
tbmid=tnew(where(fu0_bis/fu1_bis eq rap))
;=====================================================================================
nuovo:
window,2,retain=1
par1=[max(flux)*10.,para[1]]
f1=plaw(time,par1)
par2=[max(flux)*1.d2,para[2]]
f2=plaw(time,par2)
par3=[max(flux)*10.,para[6]]
f3=plaw(time,par1)
par4=[max(flux)*1.d2,para[6]]
f4=plaw(time,par2)
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,f1,linestyle=0,color=blue,thick=3.
oplot,time-t90,f2,linestyle=0,color=green,thick=3.
oplot,time-t90,f3,linestyle=0,color=orange,thick=3.
oplot,time-t90,f4,linestyle=0,color=pink,thick=3.
oplot,[para[0]-t90,para[0]-t90],[1.d-10,1.d10],linestyle=2,color=red,thick=3.
oplot,[para[5]-t90,para[5]-t90],[1.d-10,1.d10],linestyle=2,color=red,thick=3.
oplot,[tbmid-t90,tbmid-t90],[1.d-10,1.d10],linestyle=2,color=pink,thick=3.

openw,1,guessdata
printf,1,para[0],' ',para[1],' ',para[2],' ',para[3],' ',para[4],format='(e,a,e,a,e,a,e,a,e)'
printf,1,para[5],' ',para[6],' ',para[7],' ',para[8],' ',para[9],format='(e,a,e,a,e,a,e,a,e)'
close,1
print,'*******************************************'
print,'FIT PARAMETERS:'
print,''
print,'CHI^2',chi
print,'DOF',DOF
print,'p_value',p_val
print,'Break time A= ',para[0],'+/-',errore[0] 
print,'Alpha 1 A= ',para[1],'+/-',errore[1]
print,'Alpha 2 A= ',para[2],'+/-',errore[2]
print,'Smoothness A= ',para[3],'+/-',errore[3]
print,'Normalization A= ',para[4],'+/-',errore[4],' 10^(-10) erg/s/cm^2/A'
print,' '
print,'Break time B= ',para[5],'+/-',errore[5] 
print,'Alpha 1 B= ',para[6],'+/-',errore[6]
print,'Alpha 2 B= ',para[7],'+/-',errore[7]
print,'Smoothness B= ',para[8],'+/-',errore[8]
print,'Normalization B= ',para[9],'+/-',errore[9],' 10^(-10) erg/s/cm^2/A'
print,' '
print,'Break time in between the 2 curves= ',tbmid;,'+/-',err_tbmid 
print,'*******************************************'
tbmi0=tbmid
tbmi1=para[0]
tbmi2=para[5]
tbmi3=-1.
openw,2,finalfit
printf,2,'Break Time',' ','Err_BreakTime',' ','Alpha1',' ','Err_alpha1',' ','Alpha2',' ','Err_alpha2',' ','s',' ','Err_s',' ','Norm (10^48 erg/s/cm2)',' ','Err_Norm',' ','CHI2',' ','DOF',' ','P_value',format='(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
printf,2,' '
printf,2,para[0],' ',errore[0],' ',para[1],' ',errore[1],' ',para[2],' ',errore[2],' ',para[3],' ',errore[3],' ',para[4],' ',errore[4],' ',chi,' ',dof,' ',p_val,format='(e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e)'
printf,2,para[5],' ',errore[5],' ',para[6],' ',errore[6],' ',para[7],' ',errore[7],' ',para[8],' ',errore[8],' ',para[9],' ',errore[9],' ','0.',' ','0.',' ','0.',format='(e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,a,a,a,a,a)'
printf,2,tbmid,' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',format='(e,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
close,2
goto,final_part
;##########################################################################################################
;###################################         3              ###############################################
;##########################################################################################################
fit_function_3:
readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1],a2[1],s[1],norm[1],tb[2],a1[2],a2[2]] 
FF=bro3(time,guess)
set_plot,'x'
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,FF,linestyle=0,color=red,thick=3.
print,'####################################'
print,' '
print,' '
print,'WRITE .cont TO CONTINUE'
print,' '
print,' '
print,'####################################'
stop
;;openw,4,excess
;;printf,4,'Time',' ','Err_time',' ','Flux',' ','Err_Flux',' ','Residual',' ','Err_Residual',format='(a,a,a,a,a,a,a,a,a,a,a)'
;-----------------------------------------------------------
;Read the GRB#.par
openr,2,informations
nihil=' '
readf,2,nihil,nihil,nihil,z,nihil,dlo,nihil,fit_func,nihil,nihil,zero,nihil,uno,nihil,due,nihil,tre,nihil,quattro,nihil,cinque,nihil,sei,nihil,sette,nihil,otto,nihil,nove,nihil,dieci,nihil,undici,nihil,dodici,tmin,nihil,tmax,nihil
close,2
;-----------------------------------------------------------
;If we want to fix or not the parameters...
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},5)
numero=[zero,uno,due,tre,quattro,cinque,sei,sette,otto,nove,dieci,undici,dodici]
for k=0,n_elements(numero)-1 do begin
    if (numero[k] eq 1) then begin
       pi(k).fixed = 1
       guess(k) = guess[k]
    endif
endfor
;-----------------------------------------------------------
for j=0,n_elements(flux) do begin
    para=dblarr(13)
    para=MPFITFUN('bro3',time,flux,err_flux,guess,best_fjac=best_fjac, pfree_index=pfree_index, /calc_fjac, $
                 covar=pcovar,DOF=dof,BESTNORM=chi,PERROR=errore,PARINFO=pi)
    ;............................................
    kk=n_elements(para)
    nn=fltarr(kk)
    nn[0]=kk-1.
    for i=1,kk-2. do begin
        nn[i]=nn[i-1]-1.
    endfor
    gg=total(nn)
    cova3=dblarr(gg) 
    for w=0,gg-1 do begin
        for i=0,n_elements(para)-1 do begin
            for j=i+1,n_elements(para)-1 do begin
                cova3[w]=pcovar[i,j]
            endfor
        endfor
    endfor
    ;............................................
    fit=fltarr(n_elements(time))
    fit=bro3(time,para)
    residual=(flux-fit)
    residual2=(flux-fit)/fit
    ;1 sigma error obtained with 'mpfit_properr.pro'
    ycovar = mpfit_properr(best_fjac, pcovar, pfree_index)
    sy_prop = sqrt(abs(ycovar))
    rapporto=residual-fit
    err_res=fltarr(n_elements(fit))
    for oo=0,n_elements(fit)-1 do begin
        err_res[oo]=sqrt((err_flux[oo]/flux[oo])^(2.d)+(sy_prop[oo,oo]/fit[oo])^(2.d))*abs(residual2[oo])
    endfor
     guess=para
  plotsym,0,0.8,/fill
  fit=bro3(time,para)
  set_plot,'x'
  plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
  oploterror,time-t90,flux,err_time,err_flux,psym=8
  oplot,time-t90,fit,linestyle=0,color=red,thick=3.
  pvalue,n_elements(time),chi,dof,p_val
  if (p_val ge 0.3) then goto, fuori3
  if (j ge 10.) then goto, fuori3
endfor
fuori3:
;;close,4
para3=para
errore3=errore
;;readcol,excess,timesc,err_timesc,fluxsc,err_flux_sc,res_ex,err_res_exc,format='d,d,d,d,d,d'
set_plot,'x'
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
;;oploterror,timesc-t90,fluxsc,err_timesc,err_flux_sc,psym=8,color=green,errcolor=green
oplot,time-t90,fit,linestyle=0,color=red,thick=3.
;oploterror,t-t90,f,et,ef,psym=8,color=orange,errcolor=orange
;;legend,['Continuum','Removed','Hand-Removed'],color=[white,green,orange],psym=[8,8,8],box=0,/right,/top

tbmid=0.
goto,nuovo4

risposta=' '
read,risposta,prompt='The fit is ok? NO=0, YES=1:       ',format='(a)'
if (risposta eq '0') then goto, stopfinale

openw,1,guessdata
printf,1,para[0],' ',para[1],' ',para[2],' ',para[3],' ',para[4],format='(e,a,e,a,e,a,e,a,e)'
printf,1,para[5],' ',para[6],' ',para[7],' ',para[8],' ',para[9],format='(e,a,e,a,e,a,e,a,e)'
printf,1,para[10],' ',para[11],' ',para[12],' ','0.',' ','0.',format='(e,a,e,a,e,a,a,a,a)'
close,1

;===================================================================================
pp0=[para[5],para[6],para[7],para[8],para[9],para[10],para[11],para[12]]
pp1=[para[0],para[1],para[2],para[3],para[4]]
fu0=bro6(time,pp0)
fu1=bro1(time,pp1)
wh0=where(fu0 le fu1)
wh1=where(fu0 ge fu1)
tmin=time[max(wh0)]
tmax=time[min(wh1)]
dt=10.
for i=0,10000 do begin
    num=floor((tmax-tmin)/dt)
    if (num le 1) then num=3.
    if (num gt 1000) then dt=dt*2.
    tnew=dblarr(num)
    tnew[0]=tmin
    for j=1,num-1 do begin
       tnew[j]=tnew[j-1]+dt
    endfor
    fu0_bis=bro6(tnew,pp0)
    fu1_bis=bro1(tnew,pp1)
    rap=min(fu0_bis/fu1_bis)
    if (rap le rap/20.) then goto,vai3   ;5%
    dt=dt/2.   
endfor
vai3:
tbmid=tnew(where(fu0_bis/fu1_bis eq rap))
;=====================================================================================
nuovo4:
window,2,retain=1
par1=[max(flux)*10.,para[1]]
f1=plaw(time,par1)
par2=[max(flux)*1.d2,para[2]]
f2=plaw(time,par2)
par3=[max(flux)*10.,para[6]]
f3=plaw(time,par3)
par4=[max(flux)*1.d2,para[7]]
f4=plaw(time,par4)
par5=[max(flux)*1.d2,para[7]+para[11]]
f5=plaw(time,par5)
plot,time,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,f1,linestyle=0,color=blue,thick=3.
oplot,time-t90,f2,linestyle=0,color=green,thick=3.
oplot,time-t90,f3,linestyle=0,color=orange,thick=3.
oplot,time-t90,f4,linestyle=0,color=pink,thick=3.
oplot,time-t90,f5,linestyle=0,color=yellow,thick=3.
oplot,[para[0]-t90,para[0]-t90],[1.d-10,1.d10],linestyle=2,color=red,thick=3.
oplot,[para[5]-t90,para[5]-t90],[1.d-10,1.d10],linestyle=2,color=red,thick=3.
oplot,[tbmid-t90,tbmid-t90],[1.d-10,1.d10],linestyle=2,color=red,thick=3.
print,'*******************************************'
print,'FIT PARAMETERS:'
print,''
print,'CHI^2',chi
print,'DOF',DOF
print,'p_value',p_val
print,'Break time A= ',para[0],'+/-',errore[0] 
print,'Alpha 1 A= ',para[1],'+/-',errore[1]
print,'Alpha 2 A= ',para[2],'+/-',errore[2]
print,'Smoothness A= ',para[3],'+/-',errore[3]
print,'Normalization A= ',para[4],'+/-',errore[4],' 10^48 erg/s/cm^2'
print,' '
print,'Break time B= ',para[5],'+/-',errore[5] 
print,'Alpha 1 B= ',para[6],'+/-',errore[6]
print,'Alpha 2 B= ',para[7],'+/-',errore[7]
print,'Smoothness B= ',para[8],'+/-',errore[8]
print,'Normalization B= ',para[9],'+/-',errore[9],' 10^48 erg/s/cm^2'
print,' '
print,'Break time C= ',para[10],'+/-',errore[10] 
print,'Alpha 1 C= ',para[11],'+/-',errore[11]
print,'Smoothness C= ',para[12],'+/-',errore[12]
print,' '
print,'WARNING!!! The true third slope is: ',para[7]+para[6],'+/-',errore[7]+para[6]
print,' '
print,'Break time in between the 2 curves= ',tbmid;,'+',err_pos,',','-',err_neg 
print,'*******************************************'
tbmi0=tbmid
tbmi1=para[0]
tbmi2=para[5]
tbmi3=para[10]
openw,2,finalfit
printf,2,'Break Time',' ','Err_BreakTime',' ','Alpha1',' ','Err_alpha1',' ','Alpha2 (line 1,2)/s (line3)',' ','Err_alpha2 (line 1,2)/s (line 3)',' ','s',' ','Err_s',' ','Norm (10^48 erg/s/cm2)',' ','Err_Norm',' ','CHI2',' ','DOF',' ','P_value',format='(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
printf,2,' '
printf,2,para[0],' ',errore[0],' ',para[1],' ',errore[1],' ',para[2],' ',errore[2],' ',para[3],' ',errore[3],' ',para[4],' ',errore[4],' ',chi,' ',dof,' ',p_val,format='(e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e)'
printf,2,para[5],' ',errore[5],' ',para[6],' ',errore[6],' ',para[7],' ',errore[7],' ',para[8],' ',errore[8],' ',para[9],' ',errore[9],' ','0.',' ','0',' ','0.',format='(e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,a,a,a,a,a)'
printf,2,para[10],' ',errore[10],' ',para[11],' ',errore[11],' ',para[12],' ',errore[12],' ','0.',' ','0.','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',format='(e,a,e,a,e,a,e,a,e,a,e,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
printf,2,tbmid,' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',format='(e,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
close,2
goto,final_part
;##########################################################################################################
;###################################         4              ###############################################
;##########################################################################################################
fit_function_4:
readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1]] 
FF=bro4(time,guess)
set_plot,'x'
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,FF,linestyle=0,color=red,thick=3.
print,'####################################'
print,' '
print,' '
print,'WRITE .cont TO CONTINUE'
print,' '
print,' '
print,'####################################'
stop
;;openw,4,excess
;;printf,4,'Time',' ','Err_time',' ','Flux',' ','Err_Flux',' ','Residual',' ','Err_Residual',format='(a,a,a,a,a,a,a,a,a,a,a)'
;-----------------------------------------------------------
;Read the GRB#.par
openr,2,informations
nihil=' '
readf,2,nihil,nihil,nihil,z,nihil,dlo,nihil,fit_func,nihil,nihil,zero,nihil,uno,nihil,due,nihil,tre,nihil,quattro,nihil,cinque,nihil,sei,nihil,tmin,nihil,tmax
close,2
;-----------------------------------------------------------
;If we want to fix or not the parameters...
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},7)
numero=[zero,uno,due,tre,quattro,cinque,sei]
for k=0,n_elements(numero)-1 do begin
    if (numero[k] eq 1) then begin
       pi(k).fixed = 1
       guess(k) = guess[k]
    endif
endfor
;-----------------------------------------------------------
for j=0,5. do begin
    para=dblarr(7)
    para=MPFITFUN('bro4',time,flux,err_flux,guess,best_fjac=best_fjac, pfree_index=pfree_index, /calc_fjac, $
                 covar=pcovar,DOF=dof,BESTNORM=chi,PERROR=errore,PARINFO=pi)
    ;............................................
    kk=n_elements(para)
    nn=fltarr(kk)
    nn[0]=kk-1.
    for i=1,kk-2. do begin
        nn[i]=nn[i-1]-1.
    endfor
    gg=total(nn)
    cova4=fltarr(gg) 
    for w=0,gg-1 do begin
        for i=0,kk-1 do begin
            for j=i+1,kk-1 do begin
                cova4[w]=pcovar[i,j]
            endfor
        endfor
    endfor
    ;............................................
    fit=fltarr(n_elements(time))
    fit=bro4(time,para)
    residual=(flux-fit)
    residual2=(flux-fit)/fit
    ;1 sigma error obtained with 'mpfit_properr.pro'
    ycovar = mpfit_properr(best_fjac, pcovar, pfree_index)
    sy_prop = sqrt(abs(ycovar))
    rapporto=dblarr(n_elements(flux))
    for i=0,n_elements(flux)-1 do begin
    rapporto[i]=residual[i]/sy_prop[i,i]
    endfor
    err_res=fltarr(n_elements(fit))
    for oo=0,n_elements(fit)-1 do begin
        err_res[oo]=sqrt((err_flux[oo]/flux[oo])^(2.d)+(sy_prop[oo,oo]/fit[oo])^(2.d))*abs(residual2[oo])
     endfor
     guess=para
  plotsym,0,0.8,/fill
  fit=bro4(time,para)
  set_plot,'x'
  plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=2.,charsize=1.3,thick=2.
  oploterror,time-t90,flux,err_time,err_flux,psym=8
  oplot,time-t90,fit,linestyle=0,color=red,thick=3.
  pvalue,n_elements(time),chi,dof,p_val
  if (p_val ge 0.3) then goto, fuori4
  print,'j = ',j
  if (j ge 10.) then goto, fuori4
endfor
fuori4:
para4=para
errore4=errore
set_plot,'x'
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,fit,linestyle=0,color=red,thick=3.

tbmid=0.
goto,nuovo3

risposta=' '
read,risposta,prompt='The fit is ok? NO=0, YES=1:       ',format='(a)'
if (risposta eq '0') then goto, stopfinale

openw,1,guessdata
printf,1,para[0],' ',para[1],' ',para[2],' ',para[3],' ',para[4],format='(e,a,e,a,e,a,e,a,e)'
printf,1,para[5],' ',para[6],' ','0.',' ','0.',' ','0.',format='(e,a,e,a,a,a,a,a,a)'
close,1
;===================================================================================
pp0=[para[0],para[1],para[2],para[3],para[4]]
pp1=[para[5],para[6]]
if (fit_func eq 4) then begin
   fu0=bro1(time,pp0)
   fu1=plaw(time,pp1)
endif
if (fit_func eq 6) then begin
   fu0=plaw(time,pp1)
   fu1=bro1(time,pp0)
endif
wh0=where(fu0 le fu1)
wh1=where(fu0 ge fu1)
for j=0,n_elements(wh1)-1 do begin
if (wh1[j] gt -1) then tmax=time[min(wh1)]
endfor
for j=0,n_elements(wh0)-1 do begin
if (wh0[j] gt -1) then tmin=time[max(wh0)]
endfor
dt=10.
for i=0,1000 do begin
    num=floor((tmax-tmin)/dt)
    if (num le 1) then num=3.
    if (num gt 1000) then dt=dt*10000.
    tnew=dblarr(num)
    tnew[0]=tmin
    for j=1,num-1 do begin
       tnew[j]=tnew[j-1]+dt
    endfor
    fu0_bis=bro1(tnew,pp0)
    fu1_bis=plaw(tnew,pp1)
    rap=min(fu0_bis/fu1_bis)
    if (rap le rap/20.) then goto,vai4   ;5%
    dt=dt/2.   
endfor
vai4:
tbmid=tnew(where(fu0_bis/fu1_bis eq rap))
;=====================================================================================
nuovo3:
window,2,retain=1
par1=[max(flux)*10.,para[1]]
f1=plaw(time,par1)
par2=[max(flux)*1.d2,para[2]]
f2=plaw(time,par2)
par3=[max(flux)*10.,para[6]]
f3=plaw(time,par3)
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,f1,linestyle=0,color=blue,thick=3.
oplot,time-t90,f2,linestyle=0,color=green,thick=3.
oplot,time-t90,f3,linestyle=0,color=orange,thick=3.
oplot,[para[0]-t90,para[0]-t90],[1.d-10,1.d10],linestyle=2,color=red,thick=3.
oplot,[tbmid-t90,tbmid-t90],[1.d-10,1.d10],linestyle=2,color=pink,thick=3.
openw,1,guessdata
printf,1,para[0],' ',para[1],' ',para[2],' ',para[3],' ',para[4],format='(e,a,e,a,e,a,e,a,e)'
printf,1,para[5],' ',para[6],' ','0.',' ','0.',' ','0.',format='(e,a,e,a,a,a,a,a,a)'
close,1

print,'*******************************************'
print,'FIT PARAMETERS:'
print,''
print,'CHI^2',chi
print,'DOF',DOF
print,'p_value',p_val
print,'Break time A= ',para[0],'+/-',errore[0] 
print,'Alpha 1 A= ',para[1],'+/-',errore[1]
print,'Alpha 2 A= ',para[2],'+/-',errore[2]
print,'Smoothness A= ',para[3],'+/-',errore[3]
print,'Normalization A= ',para[4],'+/-',errore[4],' 10^(-10) erg/s/cm^2/A'
print,' '
print,'Normalization B= ',para[5],'+/-',errore[5],' 10^(-10) erg/s/cm^2/A'
print,'Alpha 1 B= ',para[6],'+/-',errore[6]
print,' '
print,'Break time in between the 2 curves= ',tbmid;,'+',err_pos,',','-',err_neg 
print,'*******************************************'
tbmi0=tbmid
tbmi1=-1.
tbmi2=para[0]
tbmi3=-1.
openw,2,finalfit
printf,2,'Break Time',' ','Err_BreakTime',' ','Alpha1',' ','Err_alpha1',' ','Alpha2 (line 1,2)/s (line3)',' ','Err_alpha2 (line 1,2)/s (line 3)',' ','s',' ','Err_s',' ','Norm (10^(-10) erg/s/cm2)',' ','Err_Norm',' ','CHI2',' ','DOF',' ','P_value',format='(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
printf,2,' '
printf,2,para[0],' ',errore[0],' ',para[1],' ',errore[1],' ',para[2],' ',errore[2],' ',para[3],' ',errore[3],' ',para[4],' ',errore[4],' ',chi,' ',dof,' ',p_val,format='(e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e)'
printf,2,para[5],' ',errore[5],' ',para[6],' ',errore[6],' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0',' ','0.',format='(e,a,e,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
printf,2,tbmid,' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',format='(e,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
close,2
goto,final_part
;##########################################################################################################
;##########################################################################################################
;###################################         5              ###############################################
;##########################################################################################################
;##########################################################################################################
fit_function_5:
readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1],a2[1],tb[2],a1[2]] 
FF=bro1(time,guess)
set_plot,'x'
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='Flux (erg/s/)'+cm2,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,FF,linestyle=0,color=red,thick=3.
print,'####################################'
print,' '
print,' '
print,'WRITE .cont TO CONTINUE'
print,' '
print,' '
print,'####################################'
stop
;;openw,4,excess
;;printf,4,'Time',' ','Err_time',' ','Flux',' ','Err_Flux',' ','Residual',' ','Err_Residual',format='(a,a,a,a,a,a,a,a,a,a,a)'
;-----------------------------------------------------------
;Read the GRB#.par
openr,2,informations
nihil=' '
readf,2,nihil,nihil,nihil,z,nihil,dlo,nihil,fit_func,nihil,nihil,zero,nihil,uno,nihil,due,nihil,tre,nihil,quattro,nihil,cinque,nihil,sei,nihil,sette,nihil,otto,nihil,nove,tmin,nihil,tmax,nihil
close,2
;-----------------------------------------------------------
;If we want to fix or not the parameters...
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},5)
numero=[zero,uno,due,tre,quattro,cinque,sei,sette,otto,nove]
for k=0,n_elements(numero)-1 do begin
    if (numero[k] eq 1) then begin
       pi(k).fixed = 1
       guess(k) = guess[k]
    endif
endfor
;-----------------------------------------------------------
for j=0,n_elements(flux) do begin
    para=dblarr(12)
    para=MPFITFUN('bro5',time,flux,err_flux,guess,best_fjac=best_fjac, pfree_index=pfree_index, /calc_fjac, $
                 covar=pcovar,DOF=dof,BESTNORM=chi,PERROR=errore,PARINFO=pi)
    ;............................................
    kk=n_elements(para)
    nn=fltarr(kk)
    nn[0]=kk-1.
    for i=1,kk-2. do begin
        nn[i]=nn[i-1]-1.
    endfor
    gg=total(nn)
    cova5=dblarr(gg) 
    for w=0,gg-1 do begin
        for i=0,n_elements(para)-1 do begin
            for j=i+1,n_elements(para)-1 do begin
                cova5[w]=pcovar[i,j]
            endfor
        endfor
    endfor
    ;............................................
    fit=fltarr(n_elements(time))
    fit=bro5(time,para)
    residual=flux-fit
    residual2=(flux-fit)/fit
    ;1 sigma error obtained with 'mpfit_properr.pro'
    ycovar = mpfit_properr(best_fjac, pcovar, pfree_index)
    sy_prop = sqrt(abs(ycovar))
    err_res=fltarr(n_elements(fit))
    for oo=0,n_elements(fit)-1 do begin
        err_res[oo]=sqrt((err_flux[oo]/flux[oo])^(2.d)+(sy_prop[oo,oo]/fit[oo])^(2.d))*abs(residual2[oo])
    endfor
    rapporto=residual/sy_prop
     guess=para
  plotsym,0,0.8,/fill
  fit=bro5(time,para)
  set_plot,'x'
  plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
  oploterror,time-t90,flux,err_time,err_flux,psym=8
  oplot,time-t90,fit,linestyle=0,color=red,thick=3.
  pvalue,n_elements(time),chi,dof,p_val
  if (p_val ge 0.3) then goto, fuori5
  if (j ge 10.) then goto, fuori5
endfor
fuori5:
para5=para
errore5=errore
set_plot,'x'
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oploterror,timesc-t90,fluxsc,err_timesc,err_flux_sc,psym=8,color=green,errcolor=green
oplot,time-t90,fit,linestyle=0,color=red,thick=3.
tbmid=0.
goto,nuovo5
risposta=' '
read,risposta,prompt='The fit is ok? NO=0, YES=1:       ',format='(a)'
if (risposta eq '0') then goto, stopfinale

openw,1,guessdata
printf,1,para[0],' ',para[1],' ',para[2],' ',para[3],' ',para[4],format='(e,a,e,a,e,a,e,a,e)'
printf,1,para[5],' ',para[6],' ',para[7],' ','0.',' ','0. ',format='(e,a,e,a,e,a,a,a,a)'
printf,1,para[8],' ',para[9],' ','0.',' ','0.',' ','0.',format='(e,a,e,a,a,a,a,a,a)'
close,1
;===================================================================================
pp0=[para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7]]
pp1=[para[5],para[6]]
fu0=bro6(time,pp0)
fu1=plaw(time,pp1)
wh0=where(fu0 le fu1)
wh1=where(fu0 ge fu1)
tmin=time[max(wh0)]
tmax=time[min(wh1)]
dt=10.
for i=0,10000 do begin
    num=floor((tmax-tmin)/dt)
    if (num le 1) then num=3.
    if (num gt 1000) then dt=dt*2.
    tnew=dblarr(num)
    tnew[0]=tmin
    for j=1,num-1 do begin
       tnew[j]=tnew[j-1]+dt
    endfor
    fu0_bis=bro6(tnew,pp0)
    fu1_bis=plaw(tnew,pp1)
    rap=min(fu0_bis/fu1_bis)
    if (rap le rap/20.) then goto,vai5   ;5%
    dt=dt/2.   
endfor
vai5:
tbmid1=tnew(where(fu0_bis/fu1_bis eq rap))
;=====================================================================================
nuovo5:

window,2,retain=1
par1=[max(flux)*10.,para[1]]
f1=plaw(time,par1)
par2=[max(flux)*1.d2,para[2]]
f2=plaw(time,par2)
par3=[max(flux)*10.,para[9]]
f3=plaw(time,par3)
par4=[max(flux)*10.,para[1]+para[2]]
f4=plaw(time,par4)
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,f1,linestyle=0,color=blue,thick=3.
oplot,time-t90,f2,linestyle=0,color=green,thick=3.
oplot,time-t90,f3,linestyle=0,color=orange,thick=3.
oplot,time-t90,f4,linestyle=0,color=pink,thick=3.
oplot,[para[0]-t90,para[0]-t90],[1.d-10,1.d10],linestyle=2,color=red,thick=3.
oplot,[para[5]-t90,para[5]-t90],[1.d-10,1.d10],linestyle=2,color=red,thick=3.
oplot,[tbmid-t90,tbmid-t90],[1.d-10,1.d10],linestyle=2,color=red,thick=3.
print,'*******************************************'
print,'FIT PARAMETERS:'
print,''
print,'CHI^2',chi
print,'DOF',DOF
print,'p_value',p_val
print,'Break time A= ',para[0],'+/-',errore[0] 
print,'Alpha 1 A= ',para[1],'+/-',errore[1]
print,'Alpha 2 A= ',para[2],'+/-',errore[2]
print,'Smoothness A= ',para[3],'+/-',errore[3]
print,'Normalization A= ',para[4],'+/-',errore[4],' 10^(-10) erg/s/cm^2/A'
print,' '
print,'Break time B= ',para[5],'+/-',errore[5] 
print,'Alpha 1 B= ',para[6],'+/-',errore[6]
print,'Normalization B= ',para[7],'+/-',errore[7],' 10^(-10) erg/s/cm^2/A'
print,' '
print,'Normalization C= ',para[8],'+/-',errore[8] 
print,'Alpha 1 C= ',para[9],'+/-',errore[9]
print,' '
print,'WARNING!!! The true third slope is: ',para[1]+para[2],'+/-',errore[1]+para[2]
print,' '
print,'Break time in between the 2 curves= ',tbmid1,'+',err_pos,',','-',err_neg
print,'*******************************************'
tbmi0=tbmid
tbmi1=tbmid1
tbmi2=para[0]
tbmi3=para[1]
openw,2,finalfit
printf,2,'Break Time',' ','Err_BreakTime',' ','Alpha1',' ','Err_alpha1',' ','Alpha2 (line 1,2)/s (line3)',' ','Err_alpha2 (line 1,2)/s (line 3)',' ','s',' ','Err_s',' ','Norm (10^48 erg/s/cm2)',' ','Err_Norm',' ','CHI2',' ','DOF',' ','P_value',format='(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
printf,2,' '
printf,2,para[0],' ',errore[0],' ',para[1],' ',errore[1],' ',para[2],' ',errore[2],' ',para[3],' ',errore[3],' ',para[4],' ',errore[4],' ',chi,' ',dof,' ',p_val,format='(e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e)'
printf,2,para[5],' ',errore[5],' ',para[6],' ',errore[6],' ',para[7],' ',errore[7],' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0',' ','0.',format='(e,a,e,a,e,a,e,a,e,a,e,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
printf,2,para[8],' ',errore[8],' ',para[9],' ',errore[9],' ','0.',' ','0.',' ','0.',' ','0.','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',format='(e,a,e,a,e,a,e,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
printf,2,tbmid1,' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',format='(e,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
close,2
goto,final_part
;##########################################################################################################
;###################################         8              ###############################################
;##########################################################################################################
fit_function_8:
readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1],a2[1],s[1],norm[1],tb[2],a1[2],a2[2],s[2],norm[2]] 
FF=bro8(time,guess)
set_plot,'x'
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,FF,linestyle=0,color=red,thick=3.
print,'####################################'
print,' '
print,' '
print,'WRITE .cont TO CONTINUE'
print,' '
print,' '
print,'####################################'
stop
;; openw,4,excess
;; printf,4,'Time',' ','Err_time',' ','Flux',' ','Err_Flux',' ','Residual',' ','Err_Residual',format='(a,a,a,a,a,a,a,a,a,a,a)'

;-----------------------------------------------------------
;Read the GRB#.par
openr,2,informations
nihil=' '
readf,2,nihil,nihil,nihil,z,nihil,dlo,nihil,fit_func,nihil,nihil,zero,nihil,uno,nihil,due,nihil,tre,nihil,quattro,nihil,cinque,nihil,sei,nihil,sette,nihil,otto,nihil,nove,nihil,dieci,nihil,undici,nihil,dodici,nihil,tredici,nihil,quattordici,nihil,tmin,nihil,tmax
close,2
;-----------------------------------------------------------
;If we want to fix or not the parameters...
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},15)
numero=[zero,uno,due,tre,quattro,cinque,sei,sette,otto,nove,dieci,undici,dodici,tredici,quattordici]
for k=0,n_elements(numero)-1 do begin
    if (numero[k] eq 1) then begin
       pi(k).fixed = 1
       guess(k) = guess[k]
    endif
endfor
;-----------------------------------------------------------
for j=0,n_elements(flux) do begin
    para=dblarr(15)
    para=MPFITFUN('bro8',time,flux,err_flux,guess,best_fjac=best_fjac, pfree_index=pfree_index, /calc_fjac, $
                 covar=pcovar,DOF=dof,BESTNORM=chi,PERROR=errore,PARINFO=pi)
    ;............................................
    kk=n_elements(para)
    nn=fltarr(kk)
    nn[0]=kk-1.
    for i=1,kk-2. do begin
        nn[i]=nn[i-1]-1.
    endfor
    gg=total(nn)
    cova2=dblarr(gg) 
    for w=0,gg-1 do begin
        for i=0,n_elements(para)-1 do begin
            for j=i+1,n_elements(para)-1 do begin
                cova2[w]=pcovar[i,j]
            endfor
        endfor
    endfor
    ;............................................
    fit=fltarr(n_elements(time))
    fit=bro8(time,para)
    residual=(flux-fit)
    residual2=(flux-fit)/fit
    ;1 sigma error obtained with 'mpfit_properr.pro'
    ycovar = mpfit_properr(best_fjac, pcovar, pfree_index)
    sy_prop = sqrt(abs(ycovar))
    rapporto=residual/sy_prop
    err_res=fltarr(n_elements(fit))
    for oo=0,n_elements(fit)-1 do begin
        err_res[oo]=sqrt((err_flux[oo]/flux[oo])^(2.d)+(sy_prop[oo,oo]/fit[oo])^(2.d))*abs(residual2[oo])
    endfor
  guess=para
  plotsym,0,0.8,/fill
  fit=bro8(time,para)
  set_plot,'x'
  plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
  oploterror,time-t90,flux,err_time,err_flux,psym=8
  oplot,time-t90,fit,linestyle=0,color=red,thick=3.
  pvalue,n_elements(time),chi,dof,p_val
  if (p_val ge 0.3) then goto, fuori8
  if (j ge 10.) then goto, fuori8
endfor
fuori8:
window,0,retain=1
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,fit,linestyle=0,color=red,thick=3.

risposta=' '
read,risposta,prompt='The fit is ok? NO=0, YES=1:       ',format='(a)'
if (risposta eq '0') then goto, stopfinale

;===================================================================================

openw,1,guessdata
printf,1,para[0],' ',para[1],' ',para[2],' ',para[3],' ',para[4],format='(e,a,e,a,e,a,e,a,e)'
printf,1,para[5],' ',para[6],' ',para[7],' ',para[8],' ',para[9],format='(e,a,e,a,e,a,e,a,e)'
printf,1,para[10],' ',para[11],' ',para[12],' ',para[13],' ',para[14],format='(e,a,e,a,e,a,e,a,e)'
close,1
print,'*******************************************'
print,'FIT PARAMETERS:'
print,''
print,'CHI^2',chi
print,'DOF',DOF
print,'p_value',p_val
print,'Break time A= ',para[0],'+/-',errore[0] 
print,'Alpha 1 A= ',para[1],'+/-',errore[1]
print,'Alpha 2 A= ',para[2],'+/-',errore[2]
print,'Smoothness A= ',para[3],'+/-',errore[3]
print,'Normalization A= ',para[4],'+/-',errore[4],' 10^(-10) erg/s/cm^2/A'
print,' '
print,'Break time B= ',para[5],'+/-',errore[5] 
print,'Alpha 1 B= ',para[6],'+/-',errore[6]
print,'Alpha 2 B= ',para[7],'+/-',errore[7]
print,'Smoothness B= ',para[8],'+/-',errore[8]
print,'Normalization B= ',para[9],'+/-',errore[9],' 10^(-10) erg/s/cm^2/A'
print,' '
print,'Break time A= ',para[10],'+/-',errore[10] 
print,'Alpha 1 A= ',para[11],'+/-',errore[11]
print,'Alpha 2 A= ',para[12],'+/-',errore[12]
print,'Smoothness A= ',para[13],'+/-',errore[13]
print,'Normalization A= ',para[14],'+/-',errore[14],' 10^(-10) erg/s/cm^2/A'
print,' '
print,' '
;print,'Break time in between the 2 curves= ',tbmid;,'+/-',err_tbmid 
print,'*******************************************'
tbmi0=para[0]
tbmi1=para[5]
tbmi2=para[10]
tbmi3=-1.
openw,2,finalfit
printf,2,'Break Time',' ','Err_BreakTime',' ','Alpha1',' ','Err_alpha1',' ','Alpha2',' ','Err_alpha2',' ','s',' ','Err_s',' ','Norm (10^48 erg/s/cm2)',' ','Err_Norm',' ','CHI2',' ','DOF',' ','P_value',format='(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
printf,2,' '
printf,2,para[0],' ',errore[0],' ',para[1],' ',errore[1],' ',para[2],' ',errore[2],' ',para[3],' ',errore[3],' ',para[4],' ',errore[4],' ',chi,' ',dof,' ',p_val,format='(e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e)'
printf,2,para[5],' ',errore[5],' ',para[6],' ',errore[6],' ',para[7],' ',errore[7],' ',para[8],' ',errore[8],' ',para[9],' ',errore[9],' ','0.',' ','0.',' ','0.',format='(e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,a,a,a,a,a)'
printf,2,para[10],' ',errore[10],' ',para[11],' ',errore[11],' ',para[12],' ',errore[12],' ',para[13],' ',errore[13],' ',para[14],' ',errore[14],' ','0.',' ','0.',' ','0.',format='(e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,a,a,a,a,a)'
close,2
goto,final_part
;##########################################################################################################
;##########################################################################################################
;##########################################################################################################
;##########################################################################################################
;###################################         9              ###############################################
;##########################################################################################################
fit_function_9:
readcol,guessdata,tb,a1,a2,s,norm,format='d,d,d,d,d'
guess=[tb[0],a1[0],a2[0],s[0],norm[0],tb[1],a1[1],a2[1],s[1],norm[1],tb[2],a1[2]] 
FF=bro9(time,guess)
set_plot,'x'
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,FF,linestyle=0,color=red,thick=3.
print,'####################################'
print,' '
print,' '
print,'WRITE .cont TO CONTINUE'
print,' '
print,' '
print,'####################################'
stop
;; openw,4,excess
;; printf,4,'Time',' ','Err_time',' ','Flux',' ','Err_Flux',' ','Residual',' ','Err_Residual',format='(a,a,a,a,a,a,a,a,a,a,a)'

;-----------------------------------------------------------
;Read the GRB#.par
openr,2,informations
nihil=' '
readf,2,nihil,nihil,nihil,z,nihil,dlo,nihil,fit_func,nihil,nihil,zero,nihil,uno,nihil,due,nihil,tre,nihil,quattro,nihil,cinque,nihil,sei,nihil,sette,nihil,otto,nihil,nove,nihil,dieci,nihil,undici,nihil,tmin,nihil,tmax
close,2
;-----------------------------------------------------------
;If we want to fix or not the parameters...
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},12)
numero=[zero,uno,due,tre,quattro,cinque,sei,sette,otto,nove,dieci,undici]
for k=0,n_elements(numero)-1 do begin
    if (numero[k] eq 1) then begin
       pi(k).fixed = 1
       guess(k) = guess[k]
    endif
endfor
;-----------------------------------------------------------
for j=0,n_elements(flux) do begin
    para=dblarr(12)
    para=MPFITFUN('bro9',time,flux,err_flux,guess,best_fjac=best_fjac, pfree_index=pfree_index, /calc_fjac, $
                 covar=pcovar,DOF=dof,BESTNORM=chi,PERROR=errore,PARINFO=pi)
    ;............................................
    kk=n_elements(para)
    nn=fltarr(kk)
    nn[0]=kk-1.
    for i=1,kk-2. do begin
        nn[i]=nn[i-1]-1.
    endfor
    gg=total(nn)
    cova2=dblarr(gg) 
    for w=0,gg-1 do begin
        for i=0,n_elements(para)-1 do begin
            for j=i+1,n_elements(para)-1 do begin
                cova2[w]=pcovar[i,j]
            endfor
        endfor
    endfor
    ;............................................
    fit=fltarr(n_elements(time))
    fit=bro9(time,para)
    residual=(flux-fit)
    residual2=(flux-fit)/fit
    ;1 sigma error obtained with 'mpfit_properr.pro'
    ycovar = mpfit_properr(best_fjac, pcovar, pfree_index)
    sy_prop = sqrt(abs(ycovar))
    rapporto=residual/sy_prop
    err_res=fltarr(n_elements(fit))
    for oo=0,n_elements(fit)-1 do begin
        err_res[oo]=sqrt((err_flux[oo]/flux[oo])^(2.d)+(sy_prop[oo,oo]/fit[oo])^(2.d))*abs(residual2[oo])
    endfor
     guess=para
  plotsym,0,0.8,/fill
  fit=bro9(time,para)
  set_plot,'x'
  plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
  oploterror,time-t90,flux,err_time,err_flux,psym=8
  oplot,time-t90,fit,linestyle=0,color=red,thick=3.
  pvalue,n_elements(time),chi,dof,p_val
  if (p_val ge 0.3) then goto, fuori9
  if (j ge 10.) then goto, fuori9
endfor
fuori9:
window,0,retain=1
plot,time-t90,flux,/xlog,/ylog,psym=8,xtitle='Time (s)',ytitle='10'+flusso10,charthick=1.5,charsize=1.3,thick=2.
oploterror,time-t90,flux,err_time,err_flux,psym=8
oplot,time-t90,fit,linestyle=0,color=red,thick=3.
risposta=' '
read,risposta,prompt='The fit is ok? NO=0, YES=1:       ',format='(a)'
if (risposta eq '0') then goto, stopfinale

;===================================================================================
openw,1,guessdata
printf,1,para[0],' ',para[1],' ',para[2],' ',para[3],' ',para[4],format='(e,a,e,a,e,a,e,a,e)'
printf,1,para[5],' ',para[6],' ',para[7],' ',para[8],' ',para[9],format='(e,a,e,a,e,a,e,a,e)'
printf,1,para[10],' ',para[11],' ','0.',' ','0.',' ','0.',format='(e,a,e,a,a,a,a,a,a)'
close,1
print,'*******************************************'
print,'FIT PARAMETERS:'
print,''
print,'CHI^2',chi
print,'DOF',DOF
print,'p_value',p_val
print,'Break time A= ',para[0],'+/-',errore[0] 
print,'Alpha 1 A= ',para[1],'+/-',errore[1]
print,'Alpha 2 A= ',para[2],'+/-',errore[2]
print,'Smoothness A= ',para[3],'+/-',errore[3]
print,'Normalization A= ',para[4],'+/-',errore[4],' 10^(-10) erg/s/cm^2/A'
print,' '
print,'Break time B= ',para[5],'+/-',errore[5] 
print,'Alpha 1 B= ',para[6],'+/-',errore[6]
print,'Alpha 2 B= ',para[7],'+/-',errore[7]
print,'Smoothness B= ',para[8],'+/-',errore[8]
print,'Normalization B= ',para[9],'+/-',errore[9],' 10^(-10) erg/s/cm^2/A'
print,' '
print,'Normalization A= ',para[10],'+/-',errore[10],' 10^(-10) erg/s/cm^2/A'
print,'Alpha 1 A= ',para[11],'+/-',errore[11]
print,' '
print,' '
print,'*******************************************'
tbmi0=para[0]
tbmi1=para[5]
tbmi2=para[10]
tbmi3=-1.
openw,2,finalfit
printf,2,'Break Time',' ','Err_BreakTime',' ','Alpha1',' ','Err_alpha1',' ','Alpha2',' ','Err_alpha2',' ','s',' ','Err_s',' ','Norm (10^48 erg/s/cm2)',' ','Err_Norm',' ','CHI2',' ','DOF',' ','P_value',format='(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
printf,2,' '
printf,2,para[0],' ',errore[0],' ',para[1],' ',errore[1],' ',para[2],' ',errore[2],' ',para[3],' ',errore[3],' ',para[4],' ',errore[4],' ',chi,' ',dof,' ',p_val,format='(e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e)'
printf,2,para[5],' ',errore[5],' ',para[6],' ',errore[6],' ',para[7],' ',errore[7],' ',para[8],' ',errore[8],' ',para[9],' ',errore[9],' ','0.',' ','0.',' ','0.',format='(e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,e,a,a,a,a,a,a)'
printf,2,para[10],' ',errore[10],' ',para[11],' ',errore[11],' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',' ','0.',format='(e,a,e,a,e,a,e,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)'
close,2
goto,final_part
;##########################################################################################################
;##########################################################################################################
;##########################################################################################################
;THE END
;##########################################################################################################
;##########################################################################################################
;##########################################################################################################
final_part:
set_plot,'ps'
device,filename=epsfile,/COLOR
plotsym,0,0.8,/fill

xmin=min(time)*0.1
xmax=max(time)*10.
ymin=min(flux)*10.^(-10.)*0.1
ymax=max(flux)*10.^(-10.)*10.

plot,time-t90,flux*10.^(-10.),/xlog,/ylog,psym=8,title='GRB'+nome,xtitle='Time - Observer Frame (s)',ytitle=flusso,charthick=2.,charsize=1.3,thick=2.,$
yrange=[ymin,ymax],xr=[xmin,xmax],/xst,/yst
oploterror,time-t90,flux*10.^(-10.),err_time,err_flux*10.^(-10.),psym=8
oplot,time-t90,fit*10.^(-10.),linestyle=0,color=red,thick=3.
oplot,[tbmi0-t90,tbmi0-t90],[1.d-40,1.d30],linestyle=2,color=gray,thick=3.
oplot,[tbmi1-t90,tbmi1-t90],[1.d-40,1.d30],linestyle=2,color=gray,thick=3.
oplot,[tbmi2-t90,tbmi2-t90],[1.d-40,1.d30],linestyle=2,color=gray,thick=3.
oplot,[tbmi3-t90,tbmi3-t90],[1.d-40,1.d30],linestyle=2,color=gray,thick=3.
oplot,time,para[4]*1.d-10*time^(0.829),linestyle=0,color=pink,thick=4.
device,/close
;...................................................................................
openw,6,covafile
printf,6,pcovar
close,6
;...................................................................................
openw,8,covafile_useful
for i=0,n_elements(para)-1 do begin
    for j=i+1,n_elements(para)-1 do begin
        printf,8,pcovar[i,j]
    endfor
endfor
close,8

print,'**************************************************************'
print,' '
print,'THE END :)'
print,' '
print,'If the fit is OK, you can continue with "ottico.pro".'
print,' '
print,'**************************************************************'

numero=1000.
indice=1000.
intervallo=1000.

;@@@@@@@@@@@@@@@@@@@@@
stopfinale:
stop
end
