!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: SMOPSsm_obsMod
! \label(SMOPSsm_obsMod)
!
! !INTERFACE:
module SMOPSsm_obsMod
! 
! !USES: 
  use ESMF
  use map_utils

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the observation plugin for the Land Parameter
!  Retrieval Model (LPRM) AMSR-E soil moisture product
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  12 Dec 2014: Sujay Kumar, Initial Specification
! 
!EOP
! 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOPSsm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOPSsmobs
!EOP
  type, public :: rtsmopsl2smobsdec

     character*100          :: odir
     integer                :: useASCAT
     integer                :: useWindSat
     integer                :: useSMOS
     integer                :: useSMAP
     integer                :: mo
     logical                :: startmode 
     integer                :: rtsmopsnc, rtsmopsnr
     type(proj_info)        :: rtsmopsproj
     integer, allocatable       :: n11(:)
     real,  allocatable         :: rlat(:)
     real,  allocatable         :: rlon(:)
     real,  allocatable         :: maxv(:,:)
     real,  allocatable         :: minv(:,:)
  end type rtsmopsl2smobsdec

  type(rtsmopsl2smobsdec), allocatable:: SMOPSsmobs(:)

contains
  
!BOP
! 
! !ROUTINE: SMOPSsm_obsInit
! \label{SMOPSsm_obsInit}
!
! !INTERFACE: 
  subroutine SMOPSsm_obsinit(i)
! 
! !USES: 
    use ESMF
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading LPRM AMSRE soil moisture data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: npts, status
    real                  :: gridDesci(50)

    if(.not.allocated(SMOPSsmobs)) then 
       allocate(SMOPSsmobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, SMOPSsmobs(i)%odir, &
         label='SMOPS soil moisture observation directory:', rc=status)
    call LVT_verify(status, &
         'SMOPS soil moisture observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, SMOPSsmobs(i)%useASCAT, &
         label='SMOPS soil moisture use ASCAT data:', rc=status)
    call LVT_verify(status, &
         'SMOPS soil moisture use ASCAT data: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, SMOPSsmobs(i)%useWindSat, &
         label='SMOPS soil moisture use WindSat data:', rc=status)
    call LVT_verify(status, &
         'SMOPS soil moisture use WindSat data: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, SMOPSsmobs(i)%useSMOS, &
         label='SMOPS soil moisture use SMOS data:', rc=status)
    call LVT_verify(status, &
         'SMOPS soil moisture use SMOS data: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, SMOPSsmobs(i)%useSMAP, &
         label='SMOPS soil moisture use SMAP data:', rc=status)
    call LVT_verify(status, &
         'SMOPS soil moisture use SMAP data: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    if(SMOPSsmobs(i)%useASCAT + SMOPSsmobs(i)%useWindSat + &
         SMOPSsmobs(i)%useSMOS + SMOPSsmobs(i)%useSMAP.gt.1) then
       write(LVT_logunit,*) '[ERR] Please do not select multiple sensor sources'
       write(LVT_logunit,*) '[ERR] simultaneously for LVT preprocessing..'
       write(LVT_logunit,*) '[ERR] If concurrent use of these data sources are desired,'
       write(LVT_logunit,*) '[ERR] please generate the CDF for each source separately '
       write(LVT_logunit,*) '[ERR] (using LVT) and then supply them to LIS'
       call LVT_endrun()
    endif

    SMOPSsmobs(i)%startmode = .true. 

    SMOPSsmobs(i)%rtsmopsnc = 1440
    SMOPSsmobs(i)%rtsmopsnr = 720
    
    call map_set(PROJ_LATLON, -89.875,-179.875,&
         0.0, 0.25,0.25, 0.0,&
         SMOPSsmobs(i)%rtsmopsnc,SMOPSsmobs(i)%rtsmopsnr,&
         SMOPSsmobs(i)%rtsmopsproj)
    
    gridDesci = 0
    gridDesci(1) = 0
    gridDesci(2) = 1440
    gridDesci(3) = 720
    gridDesci(4) = -89.875
    gridDesci(5) = -179.875
    gridDesci(6) = 128
    gridDesci(7) = 89.875
    gridDesci(8) = 179.875
    gridDesci(9) = 0.25
    gridDesci(10) = 0.25
    gridDesci(20) = 64

    allocate(SMOPSsmobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(SMOPSsmobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(SMOPSsmobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(SMOPSsmobs(i)%maxv(LVT_rc%lnc,LVT_rc%lnr))
    allocate(SMOPSsmobs(i)%minv(LVT_rc%lnc,LVT_rc%lnr))

    SMOPSsmobs(i)%maxv  = -1000000.0
    SMOPSsmobs(i)%minv  = 1000000.0
    call neighbor_interp_input(gridDesci,LVT_rc%gridDesc(:),&
         LVT_rc%lnc*LVT_rc%lnr,SMOPSsmobs(i)%rlat, &
         SMOPSsmobs(i)%rlon,SMOPSsmobs(i)%n11)

  end subroutine SMOPSsm_obsinit


end module SMOPSsm_obsMod
