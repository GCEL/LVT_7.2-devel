!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readALEXIObs
! \label{readALEXIObs}
!
! !INTERFACE: 
subroutine readALEXIObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_timeMgrMod
  use LVT_logMod
  use LVT_histDataMod
  use ALEXI_obsMod

  implicit none
  integer,   intent(in)   :: source
!
! !INPUT PARAMETERS: 

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  NOTES: 
!   The ALEXI output is available at monthly intervals. So 
!   the comparisons against model data should use at least a 
!   24 hour (1day) averaging interval. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  10 Dec 2010: Sujay Kumar, Initial Specification
! 
!EOP
  real                   :: currTime
  logical                :: alarmCheck 
  logical                :: file_exists
  character*100          :: lh_filename, sh_filename, gh_filename
  integer                :: c,r,ios,ftn
  logical*1              :: li(ALEXIobs(source)%nc*ALEXIobs(source)%nr)
  logical*1              :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: lht(ALEXIobs(source)%nc,ALEXIobs(source)%nr)
  real                   :: lh(ALEXIobs(source)%nc*ALEXIobs(source)%nr)
  real                   :: qle(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: sht(ALEXIobs(source)%nc,ALEXIobs(source)%nr)
  real                   :: sh(ALEXIobs(source)%nc*ALEXIobs(source)%nr)
  real                   :: qh(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: ght(ALEXIobs(source)%nc,ALEXIobs(source)%nr)
  real                   :: gh(ALEXIobs(source)%nc*ALEXIobs(source)%nr)
  real                   :: qg(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: varfield(LVT_rc%lnc, LVT_rc%lnr)

  currTime = float(LVT_rc%dhr(source))*3600+ &
       60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmCheck = (mod(currtime,86400.0).eq.0)

  if(ALEXIobs(source)%startFlag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false. 
     ALEXIobs(source)%startFlag = .false. 
     
     call create_ALEXI_lh_filename(ALEXIobs(source)%odir, ALEXIobs(source)%res,&
          LVT_rc%dyr(source), &
          LVT_rc%ddoy(source), lh_filename)

     inquire(file=lh_filename, exist=file_exists) 
     if(file_exists) then 
        ftn = LVT_getNextUnitNumber()
        write(LVT_logunit,*) '[INFO] Reading ',trim(lh_filename)
        open(ftn, file=lh_filename, form='unformatted',access='direct',&
             convert='little_endian',recl=ALEXIobs(source)%nc*ALEXIobs(source)%nr*4)
        read(ftn,rec=1) lht
        call LVT_releaseUnitNumber(ftn)
        
        li = .false. 
        do r=1,ALEXIobs(source)%nr
           do c=1,ALEXIobs(source)%nc
              lh(c+(r-1)*ALEXIobs(source)%nc) = lht(c,ALEXIobs(source)%nr-r+1)
              if(lht(c,ALEXIobs(source)%nr-r+1).ne.LVT_rc%udef.and.&
                   .not.isNaN(lht(c,ALEXIobs(source)%nr-r+1))) then 
                 li(c+(r-1)*ALEXIobs(source)%nc) = .true. 
              else
                 lh(c+(r-1)*ALEXIobs(source)%nc) = LVT_rc%udef
              endif
           enddo
        enddo
        call neighbor_interp(LVT_rc%gridDesc,li,lh,&
             lo, qle, ALEXIobs(source)%nc*ALEXIobs(source)%nr, &
             LVT_rc%lnc*LVT_rc%lnr,&
             ALEXIobs(source)%rlat, ALEXIobs(source)%rlon, &
             ALEXIobs(source)%n11,LVT_rc%udef, ios)
     else
        qle = LVT_rc%udef
     endif

     call create_ALEXI_sh_filename(ALEXIobs(source)%odir,  ALEXIobs(source)%res,&
          LVT_rc%dyr(source), &
          LVT_rc%ddoy(source), sh_filename)

     inquire(file=sh_filename, exist=file_exists) 
     if(file_exists) then 
        ftn = LVT_getNextUnitNumber()
        write(LVT_logunit,*) '[INFO] Reading ',trim(sh_filename)
        open(ftn, file=sh_filename, form='unformatted',access='direct',&
             convert='little_endian',recl=ALEXIobs(source)%nc*ALEXIobs(source)%nr*4)
        read(ftn,rec=1) sht
        call LVT_releaseUnitNumber(ftn)
        
        li = .false. 
        do r=1,ALEXIobs(source)%nr
           do c=1,ALEXIobs(source)%nc
              sh(c+(r-1)*ALEXIobs(source)%nc) = sht(c,ALEXIobs(source)%nr-r+1)
              if(sht(c,ALEXIobs(source)%nr-r+1).ne.LVT_rc%udef.or.&
                   .not.isNaN(sht(c,ALEXIobs(source)%nr-r+1))) then 
                 li(c+(r-1)*ALEXIobs(source)%nc) = .true. 
              endif
           enddo
        enddo

        call neighbor_interp(LVT_rc%gridDesc,li,sh,&
             lo, qh, ALEXIobs(source)%nc*ALEXIobs(source)%nr, &
             LVT_rc%lnc*LVT_rc%lnr,&
             ALEXIobs(source)%rlat, ALEXIobs(source)%rlon, &
             ALEXIobs(source)%n11,LVT_rc%udef, ios)
     else
        qh = LVT_rc%udef
     endif

     call create_ALEXI_gh_filename(ALEXIobs(source)%odir, ALEXIobs(source)%res,&
          LVT_rc%dyr(source), &
          LVT_rc%ddoy(source), gh_filename)

     inquire(file=gh_filename, exist=file_exists) 
     if(file_exists) then 
        ftn = LVT_getNextUnitNumber()
        write(LVT_logunit,*) '[INFO] Reading ',trim(gh_filename)
        open(ftn, file=gh_filename, form='unformatted',access='direct',&
             convert='little_endian',recl=ALEXIobs(source)%nc*ALEXIobs(source)%nr*4)
        read(ftn,rec=1) sht
        call LVT_releaseUnitNumber(ftn)
        
        li = .false. 
        do r=1,ALEXIobs(source)%nr
           do c=1,ALEXIobs(source)%nc
              gh(c+(r-1)*ALEXIobs(source)%nc) = ght(c,ALEXIobs(source)%nr-r+1)
              if(ght(c,ALEXIobs(source)%nr-r+1).ne.LVT_rc%udef.or.&
                   .not.isNaN(ght(c,ALEXIobs(source)%nr-r+1))) then 
                 li(c+(r-1)*ALEXIobs(source)%nc) = .true. 
              endif
           enddo
        enddo

        call neighbor_interp(LVT_rc%gridDesc,li,sh,&
             lo, qg, ALEXIobs(source)%nc*ALEXIobs(source)%nr, &
             LVT_rc%lnc*LVT_rc%lnr,&
             ALEXIobs(source)%rlat, ALEXIobs(source)%rlon, &
             ALEXIobs(source)%n11,LVT_rc%udef, ios)
     else
        qg = LVT_rc%udef
     endif

  else
     qle = LVT_rc%udef
     qh = LVT_rc%udef
     qg = LVT_rc%udef
  endif

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(qle(c+(r-1)*LVT_rc%lnc).ne.LVT_rc%udef) then
           varfield(c,r) = qle(c+(r-1)*LVT_rc%lnc)*1E6/86400.0
        else
           varfield(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_QLE,source,varfield,vlevel=1,units="W/m2")

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(qh(c+(r-1)*LVT_rc%lnc).ne.LVT_rc%udef) then
           varfield(c,r) = qh(c+(r-1)*LVT_rc%lnc)*1E6/86400.0
        else
           varfield(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_QH,source,varfield,vlevel=1,units="W/m2")

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(qg(c+(r-1)*LVT_rc%lnc).ne.LVT_rc%udef) then
           varfield(c,r) = qg(c+(r-1)*LVT_rc%lnc)*1E6/86400.0
        else
           varfield(c,r) = LVT_rc%udef
        endif
     enddo
  enddo

  call LVT_logSingleDataStreamVar(LVT_MOC_QG,source,varfield,vlevel=1,units="W/m2")

  
end subroutine readALEXIObs

!BOP
! 
! !ROUTINE: create_ALEXI_lh_filename
! \label{create_ALEXI_lh_filename}
!
! !INTERFACE: 
subroutine create_ALEXI_lh_filename(odir,res,yr,doy,filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for ALEXI_LH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the ALEXI_LH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr
  integer                      :: res
  integer                      :: doy
  character(len=*)             :: filename
!EOP

  character*4             :: fyr
  character*3             :: fdoy
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy,fmt='(i3.3)') doy

  if(res.eq.10) then 
     filename = trim(odir)//'/'//trim(fyr)//'/'//&
          'EDAY_EF'//trim(fyr)//trim(fdoy)//'.dat'
  elseif(res.eq.4) then 
     filename = trim(odir)//'/'//trim(fyr)//'/'//&
          'EDAY_FSUN_CFSR_'//trim(fyr)//trim(fdoy)//'.dat'
  endif

end subroutine create_ALEXI_lh_filename

!BOP
! 
! !ROUTINE: create_ALEXI_sh_filename
! \label{create_ALEXI_sh_filename}
!
! !INTERFACE: 
subroutine create_ALEXI_sh_filename(odir,res,yr,doy,filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for ALEXI_SH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the ALEXI_SH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr
  integer                      :: res
  integer                      :: doy
  character(len=*)             :: filename
!EOP

  character*4             :: fyr
  character*3             :: fdoy

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy,fmt='(i3.3)') doy

  if(res.eq.10) then 
     filename = trim(odir)//'/'//trim(fyr)//'/'//&
          'HDAY_EF'//trim(fyr)//trim(fdoy)//'.dat'
     
  elseif(res.eq.4) then 
     filename = 'NONE'
  endif
end subroutine create_ALEXI_sh_filename

!BOP
! 
! !ROUTINE: create_ALEXI_gh_filename
! \label{create_ALEXI_gh_filename}
!
! !INTERFACE: 
subroutine create_ALEXI_gh_filename(odir,res,yr,doy,filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for ALEXI_GH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the ALEXI_GH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr
  integer                      :: doy
  integer                      :: res
  character(len=*)             :: filename
!EOP

  character*4             :: fyr
  character*3             :: fdoy

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy,fmt='(i3.3)') doy

  if(res.eq.10) then 
     filename = trim(odir)//'/'//trim(fyr)//'/'//&
          'GDAY'//trim(fyr)//trim(fdoy)//'.dat'     
  elseif(res.eq.4) then 
     filename = 'NONE'
  endif
end subroutine create_ALEXI_gh_filename


