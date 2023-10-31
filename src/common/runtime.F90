SUBROUTINE runtime(ifl)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
integer imin,isec,itime,i1,i2,ifl
real etime,timediff,timearray(2)

timediff = etime(timearray)
itime = anint(timediff)
isec = mod(itime,60)
imin = (itime-isec)/60
i2 = mod(isec,10)
i1 = (isec - i2)/10
write(*,1200)'runtime (mango): ',imin,':',i1,i2,timearray(1),'u',timearray(2),'s'
write(ifl,1200)'%runtime (mango): ',imin,':',i1,i2,timearray(1),'u',timearray(2),'s'
1200  format(a,i6,a1,i1,i1,2(f10.3,a1))
 
return
end
