program STM_GENERATOR
implicit none
integer :: i,j,k,kk,iz,nmx,nmy,res
integer :: ntyp,nx,ny,nz
integer :: ibgn,iend,incr
real    :: a0,a1(3),a2(3),a3(3),area,volume
real    :: crt,zp1,zpn,zfine,zchg,aaa,bbb,splint

integer, allocatable :: image(:,:),img(:,:),natm(:)
real   , allocatable :: chg(:,:,:),cline(:),zgrid(:),zpp(:)

open(11,file='PARCHG',status='old')
open(12,file='STM.bmp',form='unformatted',status='replace',ACCESS='STREAM')

!image resolution
res=255

!print *, 'Supercell (NX, NY) :'
!read *, nmx,nmy
nmx=1; nmy=1

print *, 'How many kinds of atoms?'
read *, ntyp
allocate(natm(ntyp))

read (11,*)
read (11,*) a0
read (11,*) a1
read (11,*) a2
read (11,*) a3

area=abs(a1(1)*a2(2)-a1(2)*a2(1))
volume=area*abs(a3(3))*a0*a0*a0
print *, 'volume: ', volume

!read (11,*) ! for vasp5 elements
read (11,*) natm
read (11,*)
do i=1,sum(natm)
  read (11,*)
end do
read (11,*)
print *, 'number of atoms: ', natm

read (11,*) nx,ny,nz
print *, 'nx:',nx,'   ny:',ny,'   nz:',nz

allocate(chg(nx,ny,nz))
allocate(cline(nz))
allocate(zgrid(nz))
allocate(zpp(nz))
allocate(image(nx,ny))
allocate(img(nx*nmx,ny*nmy))

read (11,*) chg
chg=chg/volume
print *, 'total charge: ', sum(chg*volume/(nx*ny*nz))

!print *, 'Enter Z for the plane 6 Angstrom above surface:'
!read *, iz
!write(*,'(2E14.4)') maxval(chg(:,:,iz)),minval(chg(:,:,iz))

print *, 'Enter the isosurface value:'
read *, crt
print *, 'Starting and ending for Z-scanning :'
read *, ibgn, iend
if (ibgn.gt.iend) then
  incr=-1
else
  incr=1
end if

image=0
do i=1,nx
  do j=1,ny
    do k=1,nz
      zgrid(k)=k
      cline(k)=chg(i,j,k)
    end do
    zp1=(cline(2)-cline(nz))*0.5
    zpn=(cline(1)-cline(nz-1))*0.5
    call spline(nz,zgrid,cline,zp1,zpn,zpp)

    zchg=0.
    outer: do k=ibgn,iend,incr
      do kk=0,9
        zfine=k+kk*0.1*incr
        zchg=splint(nz,zgrid,cline,zpp,zfine)
        if (zchg.gt.crt.and.k.ne.ibgn) then
          image(i,j)=k*10+kk
          exit outer
        end if
        if (zchg.gt.crt.and.k.eq.ibgn) then
          write(*,*) 'Warning! ibgn not well-selected!', i, j
        end if
!!!! this part may not be always necessary.
        if (k.eq.iend) then
          image(i,j)=(iend+incr)*10
          exit outer
        end if
!!!! this part may not be always necessary.
      end do
    end do outer
  end do
end do

print *, 'maximum - minimum gray scale before scaling: ', maxval(image),minval(image)

aaa=FLOAT(res)/real(maxval(image)-minval(image))
bbb=-aaa*minval(image)
image=image*aaa+bbb
if (incr.eq.1) image=res-image
!image=image-minval(image)
!image=res-image*res/maxval(image)
print *, 'maximum - minimum gray scale after scaling : ', maxval(image),minval(image)

do i=1,nmx
  do j=1,nmy
    img((i-1)*nx+1:(i*nx),(j-1)*ny+1:(j*ny))=image
  end do
end do

call wrtbmp(12,nx*nmx,ny*nmy,img)

end program STM_GENERATOR

subroutine wrtbmp(nfile,nx,ny,img)
implicit none
integer :: nfile,nx,ny,img(nx,ny)

integer :: res
character :: identifier(2)
integer :: biSize
integer :: Reserved
integer :: Offset
integer :: HeadSize
integer :: biWidth
integer :: biHeight
character :: biPlanes
character :: biBitCount
integer :: biCompression
integer :: biSizeImage
integer :: biXPelsPerMeter
integer :: biYPelsPerMeter
integer :: biClrUsed
integer :: biClrImportant

character :: swap
integer :: i,j,nbit

if (mod(ny,4).eq.0) then
  nbit=0
else
  nbit=4-mod(ny,4)
end if

res=255
identifier(1)=CHAR(66)
identifier(2)=CHAR(77)
biSize=nx*(ny+nbit)+1080
Reserved=0
Offset=1078
HeadSize=40
biWidth=ny
biHeight=nx
biPlanes=CHAR(1)
biBitCount=CHAR(8)
biCompression=0
biSizeImage=0
biXPelsPerMeter=2834
biYPelsPerMeter=2834
biClrUsed=0
biClrImportant=0

swap=CHAR(0)
!write BITMAPINFOHEADER (54 bytes)
write(nfile) identifier(1)
write(nfile) identifier(2)
write(nfile) biSize
write(nfile) Reserved
write(nfile) Offset
write(nfile) HeadSize
write(nfile) biWidth
write(nfile) biHeight
write(nfile) biPlanes
write(nfile) swap
write(nfile) biBitCount
write(nfile) swap
write(nfile) biCompression
write(nfile) biSizeImage
write(nfile) biXPelsPerMeter
write(nfile) biYPelsPerMeter
write(nfile) biClrUsed
write(nfile) biClrImportant

!write RGBQUAD (1024 bytes)
do i=0,res
  do j=1,3
    swap=CHAR(i)
    write(nfile) swap
  end do
  swap=CHAR(0)
  write(nfile) swap
end do

!write BMP
do i=1,nx
  do j=1,ny
    swap=CHAR(img(i,j))
    write(nfile) swap
  end do
  swap=CHAR(0)
  do j=1,nbit
    write(nfile) swap
  end do
end do

! ending BMP file (2 bytes)
swap=CHAR(0)
write(nfile) swap
write(nfile) swap

end subroutine wrtbmp

!==============================================================================!
!              C u b i c   S p l i n e   I n t e r p o l a t i o n             !
!                                                                              !
! N      --- Number of tabulated data points                                   !
! X(N)   --- X(i)                                                              !
! Y(N)   --- Y(i)                                                              !
! YP1    --- First derivative at X(1)                                          !
! YPN    --- First derivative at X(N)                                          !
! YPP(N) --- The output second derivatives at all data points                  !
!==============================================================================!
SUBROUTINE SPLINE(N,X,Y,YP1,YPN,YPP)
INTEGER :: N
REAL    :: YP1,YPN
REAL    :: X(N),Y(N),YPP(N)

INTEGER :: I
REAL    :: C1,CN
REAL    :: SUB(3:N-1),DIG(2:N-1),SUP(2:N-2),RHS(2:N-1)

! Initialize the diagonal (DIG), subdiagonal (SUB) and superdiagonal (SUP)
! elements and the RHS vector(RHS). Natural cubic spline is used first.
DIG(2)=(X(3)-X(1))/3.0
SUP(2)=(X(3)-X(2))/6.0
RHS(2)=(Y(3)-Y(2))/(X(3)-X(2))-(Y(2)-Y(1))/(X(2)-X(1))
DO I=3,N-2
  SUB(I)=(X(I)-X(I-1))/6.0
  DIG(I)=(X(I+1)-X(I-1))/3.0
  SUP(I)=(X(I+1)-X(I))/6.0
  RHS(I)=(Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))/(X(I)-X(I-1))
END DO
SUB(N-1)=(X(N-1)-X(N-2))/6.0
DIG(N-1)=(X(N)-X(N-2))/3.0
RHS(N-1)=(Y(N)-Y(N-1))/(X(N)-X(N-1))-(Y(N-1)-Y(N-2))/(X(N-1)-X(N-2))

! Special cares are required at point 2 and point N-1 if the first
! derivatives are provided and less than 0.99E30.
IF(YP1.LT.0.99E30) THEN
  DIG(2)=DIG(2)-(X(2)-X(1))/12.0
  C1=-0.5*YP1+0.5*(Y(2)-Y(1))/(X(2)-X(1))
  RHS(2)=RHS(2)-C1
END IF
IF(YPN.LT.0.99E30) THEN
  DIG(N-1)=DIG(N-1)-(X(N)-X(N-1))/12.0
  CN=0.5*YPN-0.5*(Y(N)-Y(N-1))/(X(N)-X(N-1))
  RHS(N-1)=RHS(N-1)-CN
END IF

! Gaussian elimination for the tridiagonal system.
DO I=3,N-1
  DIG(I)=DIG(I)-SUB(I)*SUP(I-1)/DIG(I-1)
  RHS(I)=RHS(I)-SUB(I)*RHS(I-1)/DIG(I-1)
END DO

! Back-substitution to solve the values of YPP at points from N-1 to 2.
YPP(N-1)=RHS(N-1)/DIG(N-1)
DO I=N-2,2,-1
  YPP(I)=(RHS(I)-SUP(I)*YPP(I+1))/DIG(I)
END DO

! Finally, calculate YPP at points 1 and N.
YPP(1)=0.0
IF(YP1.LT.0.99E30) THEN
  C1=C1-YPP(2)*(X(2)-X(1))/12.0
  YPP(1)=C1*6.0/(X(2)-X(1))
END IF

YPP(N)=0.0
IF(YPN.LT.0.99E30) THEN
  CN=CN-YPP(N-1)*(X(N)-X(N-1))/12.0
  YPP(N)=CN*6.0/(X(N)-X(N-1))
END IF

END SUBROUTINE SPLINE

REAL FUNCTION SPLINT(N,X,Y,YPP,XI)
IMPLICIT NONE
INTEGER :: N
REAL    :: X(N),Y(N),YPP(N)
REAL    :: XI

INTEGER :: I,J,K
REAL    :: DX,A,B

I=1
J=N
DO WHILE (J-I.GT.1)
  K=(J+I)/2
  IF (X(K).GT.XI) THEN
    J=K
  ELSE
    I=K
  END IF
END DO
DX=X(J)-X(I)
A=(X(J)-XI)/DX
B=(XI-X(I))/DX
SPLINT=A*Y(I)+B*Y(J)+((A**3-A)*YPP(I)+(B**3-B)*YPP(J))*(DX*DX)/6.0

END FUNCTION SPLINT

