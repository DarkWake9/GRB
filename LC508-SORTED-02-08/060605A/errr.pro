pro errr
readcol,'errori.dat',e1,e2,format='d,d'
openw,1,'erroriok.dat'
for j=0,n_elements(e1)-1 do begin
     printf,1,(e1[j]+e2[j])/2.,format='(d)'
endfor
close,1
end
