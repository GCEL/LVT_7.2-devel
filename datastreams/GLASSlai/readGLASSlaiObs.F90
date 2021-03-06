!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readGLASSlaiObs
! \label{readGLASSlaiObs}
!
! !INTERFACE: 
subroutine readGLASSlaiObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use GLASSlaiobsMod
          
  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  NOTES: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 Mar 2015: Sujay Kumar, Initial Specification
! 
!EOP

  integer                :: c,r, tindex
  integer                :: flag
  integer                :: ftn
  character*100          :: fname
  logical                :: file_exists
  integer                :: iret  
  integer                :: nc, nr
  real                   :: varfield(LVT_rc%lnc,LVT_rc%lnr)

  varfield = LVT_rc%udef

  nc = GLASSlaiobs(source)%nc
  nr = GLASSlaiobs(source)%nr

  LVT_rc%resetFlag(source) = .false. 
  GLASSlaiobs(source)%startFlag = .false. 
  
  call create_GLASSlai_filename(glasslaiobs(Source)%odir, &
       LVT_rc%dyr(source),LVT_rc%ddoy(source),&
       fname)
  
  inquire(file=trim(fname),exist=file_exists) 
  
  ! Check if both files exist:
  if( file_exists ) then 
     write(LVT_logunit,*) '[INFO] Reading GLASS LAI file ',trim(fname)
     
     call read_GLASS_LAI_data(source, fname, varfield)
     
  else
     write(LVT_logunit,*)'[WARN] GLASS LAI file missing: ',trim(fname)
     varfield  = LVT_rc%udef
  endif
  
  call LVT_logSingleDataStreamVar(LVT_MOC_LAI,source,varfield,&
       vlevel=1,units="-")
  
end subroutine readGLASSlaiObs

!BOP
! 
! !ROUTINE: read_GLASS_LAI_data
! \label{read_GLASS_LAI_data}
!
! !INTERFACE:
subroutine read_GLASS_LAI_data(source, fname, laiobs_ip)
! 
! !USES:   

  use LVT_coreMod
  use LVT_logMod
  use LVT_timeMgrMod
  use GLASSlaiobsMod

  implicit none
#if (defined USE_HDF4) 
#include "hdf.f90"
#endif
!
! !INPUT PARAMETERS: 
! 
  integer                       :: source
  character (len=*)             :: fname
  real                          :: laiobs_ip(LVT_rc%lnc*LVT_rc%lnr)

! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the GLASS LAI file and applies the data
!  quality flags to filter the data. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the RTGLASS AMSR-E file
!  \item[laiobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP

#if (defined USE_HDF4)
  integer,  parameter     :: nc=7200, nr=3600
  integer                 :: gdopen,gdattach,gdrdfld
  integer                 :: gddetach,gdclose
  integer                 :: file_id,grid_id,region_id,iret,c,r,c1,r1
  character*50            :: grid_name,lai_name
  integer                 :: start(2), edge(2), stride(2)
  integer*1, allocatable  :: lai_raw(:)
  integer                 :: lat_off, lon_off
  real                    :: lai_in(GLASSlaiobs(source)%nc*GLASSlaiobs(source)%nr)
  logical*1               :: lai_data_b(GLASSlaiobs(source)%nc*GLASSlaiobs(source)%nr)
  logical*1               :: laiobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)

  if(LVT_rc%dyr(source).ge.2001) then 
     grid_name ="GLASS01B01"
  else
     grid_name ="GLASS03LAI02"
  endif

  file_id = gdopen(trim(fname),DFACC_READ)
  if (file_id.eq.-1)then
     write(LVT_logunit,*) "[ERR] Failed to open hdf file",fname
     call LVT_endrun()
  end if
  
  lai_name = "LAI"

  grid_id = gdattach(file_id,grid_name)

  start(1)=0  !hdfeos lib uses 0-based count
  start(2)=0
  edge(1)=nc
  edge(2)=nr
  stride(1)=1
  stride(2)=1
  
  allocate(lai_raw(nc*nr))
  iret = gdrdfld(grid_id,lai_name,start,stride,edge,lai_raw)
  
  if(iret.ne.0) then 
     write(LVT_logunit,*) "[ERR] Failed to read file",fname
     call LVT_endrun()
  endif

  lai_data_b = .false. 

  do r=1,GLASSlaiobs(source)%nr
     do c=1,GLASSlaiobs(source)%nc
        c1 = c 
        r1 = nr - r + 1
        if(lai_raw(c1+(r1-1)*nc).ne.-1) then 
           lai_in(c+(r-1)*GLASSlaiobs(source)%nc) =&
                lai_raw(c1+(r1-1)*nc)*0.1
           lai_data_b(c+(r-1)*GLASSlaiobs(source)%nc) =  .true. 
        else
           lai_in(c+(r-1)*GLASSlaiobs(source)%nc) = -9999.0
           lai_data_b(c+(r-1)*GLASSlaiobs(source)%nc) = .false. 
        endif
     enddo
  enddo
  deallocate(lai_raw)
  
  iret=gddetach(grid_id)
  iret=gdclose(file_id)

  if(LVT_isAtAfinerResolution(GLASSlaiobs(source)%datares)) then
!--------------------------------------------------------------------------
! Interpolate to the LVT analysis domain
!-------------------------------------------------------------------------- 
     call bilinear_interp(LVT_rc%gridDesc,&
          lai_data_b, lai_in, laiobs_b_ip, laiobs_ip, &
          GLASSlaiobs(source)%nc*GLASSlaiobs(source)%nr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          GLASSlaiobs(source)%rlat,GLASSlaiobs(source)%rlon,&
          GLASSlaiobs(source)%w11,GLASSlaiobs(source)%w12,&
          GLASSlaiobs(source)%w21,GLASSlaiobs(source)%w22,&
          GLASSlaiobs(source)%n11,GLASSlaiobs(source)%n12,&
          GLASSlaiobs(source)%n21,GLASSlaiobs(source)%n22,LVT_rc%udef,iret)

  else
     call upscaleByAveraging(GLASSlaiobs(source)%nc*GLASSlaiobs(source)%nr,&
          LVT_rc%lnc*LVT_rc%lnr, &
          LVT_rc%udef, GLASSlaiobs(source)%n11,&
          lai_data_b,lai_in, laiobs_b_ip, laiobs_ip)
  endif
  
#endif

end subroutine read_GLASS_LAI_data


!BOP
! !ROUTINE: create_GLASSlai_filename
! \label{create_GLASSlai_filename}
! 
! !INTERFACE: 
subroutine create_GLASSlai_filename(ndir,  yr, doy, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, doy
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the GLASS LAI filename
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the GLASS LAI data directory
!  \item[yr]  current year
!  \item[mo]  current doy
!  \item[filename] Generated GLASS LAI filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=3) :: fdoy
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy

  filename = trim(ndir)//'/'//trim(fyr)//'/GLASS01B02.V03.A'//&
       trim(fyr)//trim(fdoy)//'.hdf'

end subroutine create_GLASSlai_filename


