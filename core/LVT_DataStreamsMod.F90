!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
#include "LVT_NetCDF_inc.h"
module LVT_DataStreamsMod
!BOP
! 
! !MODULE: LVT_DataStreamMod
! \label(LVT_DataStreamMod)
!
! !INTERFACE:
! 
! !USES:   
  use LVT_histDataMod
  use LVT_coreMod
  use LVT_logMod
  use LVT_LISoutputHandlerMod
  use LVT_timeMgrMod
  use map_utils
  use grib_api
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  The code in this file contains the basic datastructures and 
!  control routines for handling the operations associated 
!  with the datastreams. The invocation to read, perform 
!  temporal averaging and resetting of the datastreams are 
!  performed from this module. The calculations of derived
!  variables are also performed in this module.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_DataStreamsInit
  public :: LVT_readDataStreams
  public :: LVT_writeDataStreams
  public :: LVT_tavgDataStreams
  public :: LVT_resetDataStreams
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  
!EOP
  
contains

!BOP
! 
! !ROUTINE: LVT_DataStreamsInit
! \label{LVT_DataStreamsInit}
!
! !INTERFACE:   
  subroutine LVT_DataStreamsInit
! 
! !USES:   
    use LVT_datastream_pluginMod, only : LVT_datastream_plugin
    implicit none
!
! !DESCRIPTION:
! 
!  This subroutine invokes the call to initialize each datastream.
!
!   The routines invoked are: 
!   \begin{description}
!    \item[LVT\_datastream\_plugin] (\ref{LVT_datastream_plugin}) \newline
!      routine to register all the supported datastream plugin implementations
!    \item[LVT\_LISoutputInit] (\ref{LVT_LISoutputInit}) \newline
!      routine to initialize the handling of LIS output data (if LIS output
!      data is one of the datastreams)
!   \end{description} 
!EOP
    integer                          :: kk
    type(LVT_metadataEntry), pointer :: ds1, ds2
    real                             :: gridDesci(50)

    call LVT_datastream_plugin
    
    call observationsetup(trim(LVT_rc%obssource(1))//char(0),1)
    call observationsetup(trim(LVT_rc%obssource(2))//char(0),2)

    if(LVT_rc%nDatastreams.gt.2) then 
       call observationsetup(trim(LVT_rc%obssource(3))//char(0),3)
    endif

    call LVT_LISoutputInit

! checking for duplicate entries in a given datastream
! Note that this check is not enabled for three datastrems. 
! The responsibility of ensuring non-duplicate entries is
! on the user. 

    LVT_rc%ds1_dup = .false. 
    ds1 => LVT_histData%head_ds1_list       
    do while(associated(ds1))
       ds2 => ds1%next
       do while(associated(ds2))
          if(ds2%index.ne.ds1%index.and.&
               ds1%short_name.eq.ds2%short_name) then 
             LVT_rc%ds1_dup = .true. 
          endif
          ds2 => ds2%next
       enddo
       ds1 => ds1%next
    enddo

    LVT_rc%ds2_dup = .false. 
    ds1 => LVT_histData%head_ds2_list       
    do while(associated(ds1))
       ds2 => ds1%next
       do while(associated(ds2))
          if(ds2%index.ne.ds1%index.and.&
               ds1%short_name.eq.ds2%short_name) then 
             LVT_rc%ds2_dup = .true. 
          endif
          ds2 => ds2%next
       enddo
       ds1 => ds1%next
    enddo

!-------------------------------------------------------------------
! for 557 post, the HYCOM data is processed to include the water
! temperature fields
!-------------------------------------------------------------------
    if(LVT_rc%runmode.eq."557 post") then 
       if(LVT_rc%processHYCOM.eq.1) then 

          allocate(LVT_rc%HYCOM_cind(LVT_rc%lnc,LVT_rc%lnr))
          allocate(LVT_rc%HYCOM_rind(LVT_rc%lnc,LVT_rc%lnr))
          LVT_rc%HYCOM_proc_start = .true. 

          LVT_rc%HYCOM_nc = 4500
          LVT_rc%HYCOM_nr = 2001

          gridDesci = 0 
          gridDesci(1) = 0 
          gridDesci(2) = LVT_rc%HYCOM_nc
          gridDesci(3) = LVT_rc%HYCOM_nr
          gridDesci(4) = -80.0
          gridDesci(5) = -180.0
          gridDesci(7) = 80.0
          gridDesci(8) = 180.0
          gridDesci(6) = 128
          gridDesci(9) = 0.08
          gridDesci(10) = 0.08
          gridDesci(20) = 64

          allocate(LVT_rc%HYCOM_n11(LVT_rc%HYCOM_nc*LVT_rc%HYCOM_nr))

          call upscaleByAveraging_input(gridDesci, LVT_rc%gridDesc,&
               LVT_rc%HYCOM_nc*LVT_rc%HYCOM_nr, &
               LVT_rc%lnc*LVT_rc%lnr, LVT_rc%HYCOM_n11)

          LVT_histData%watertemp%short_name = "water_temp"
          LVT_histData%watertemp%long_name = "water_temp"
          LVT_histData%watertemp%standard_name = "water temperature"
          LVT_histData%watertemp%units = "K"
          LVT_histData%watertemp%nunits = 1
          LVT_histData%watertemp%format = 'F'
          LVT_histData%watertemp%vlevels = 1
          LVT_histData%watertemp%timeAvgOpt = 0 
          LVT_histData%watertemp%startNlevs = 1
          LVT_histData%watertemp%endNlevs = 1
          allocate(LVT_histData%watertemp%value(LVT_rc%ngrid,&
               1,LVT_histData%watertemp%vlevels))
          allocate(LVT_histData%watertemp%unittypes(1))
          LVT_histData%watertemp%unittypes(1) = "K"
          
       endif
    endif
  end subroutine LVT_DataStreamsInit

!BOP
! 
! !ROUTINE: LVT_readDataStreams
! \label{LVT_readDataStreams}
!
! !INTERFACE: 
  subroutine LVT_readDataStreams
! 
! !USES:   
    implicit none
!
!
! !DESCRIPTION: 
!  This subroutine invokes the routines that read the datastreams
! 
!EOP

    call readObservationSource(trim(LVT_rc%obssource(1))//char(0),1)      
    call readObservationSource(trim(LVT_rc%obssource(2))//char(0),2)      

    if(LVT_rc%nDatastreams.gt.2) then 
       call readObservationSource(trim(LVT_rc%obssource(3))//char(0),3)      
    endif
  end subroutine LVT_readDataStreams

!BOP
! 
! !ROUTINE: LVT_writeDataStreams
! \label{LVT_writeDataStreams}
!
! !INTERFACE: 
  subroutine LVT_writeDataStreams
! 
! !USES:   
    use LVT_logMod

    implicit none
!
!
! !DESCRIPTION: 
!  This subroutine invokes the routines that writes the datastream values
!  to an external file. Currently this feature is only supported
!  in the '557 post' runmode, for the processing of LIS outputs to the
!  grib format. The datastream1 output must be set to 'LIS output'. 
! 
!EOP


    integer, parameter                   :: nsoillayers = 4
    character*200                        :: fname_mean,fname_ssdev
    character(len=8)                     :: cdate2
    character(len=4)                     :: cdate3
    integer                              :: ftn_mean,ftn_ssdev
    integer                              :: time_unit
    integer                              :: time_past
    integer                              :: time_curr
    integer                              :: iret
    integer                              :: timeRange
    integer                              :: yr,mo,da,hr,mn,ss
    character*10                         :: stepType
    type(LVT_metadataEntry), pointer     :: dataEntry
    type(LVT_LISmetadataEntry), pointer  :: lisdataEntry
    real, dimension(nsoillayers)         :: toplev, botlev, depscale
    real, dimension(1)                   :: toplev0, botlev0
    real                                 :: lyrthk(nsoillayers)
    integer                              :: i,k,t,m
    real*8                               :: time
    real                                 :: gmt
    integer                              :: doy
    integer                              :: pdTemplate
    integer                              :: c,r,gid
    real                                 :: gtmp1_1d(LVT_rc%lnc*LVT_rc%lnr)
    real                                 :: gtmp1_ss(LVT_rc%lnc*LVT_rc%lnr)
    integer                              :: ngtmp1_1d(LVT_rc%lnc*LVT_rc%lnr)
    real                                 :: variance
    real                                 :: lat(LVT_rc%lnc,LVT_rc%lnr)
    real                                 :: lon(LVT_rc%lnc,LVT_rc%lnr)
    character*20                         :: output_fmt
    integer                              :: shuffle, deflate, deflate_level
    integer                              :: xlatID,xlonID,xtimeID
    integer                              :: xlat_ss_ID,xlon_ss_ID,xtime_ss_ID
    integer                              :: dimID(4),tdimID
    character*8                          :: xtime_begin_date
    character*6                          :: xtime_begin_time
    character*50                         :: xtime_units
    character*50                         :: xtime_timeInc
    character(len=8)                     :: date
    character(len=10)                    :: stime
    character(len=5)                     :: zone
    integer, dimension(8)                :: values
    type(LVT_metadataEntry)              :: xlat,xlon

!    output_fmt = "grib2"

    lyrthk(1) = 0.1*100.0
    lyrthk(2) = 0.3*100.0
    lyrthk(3) = 0.6*100.0
    lyrthk(4) = 1.0*100.0


    if(LVT_rc%runmode.eq."557 post") then 
       if(LVT_rc%lvt_out_format.eq."grib1") then  

          write(unit=cdate2,fmt='(i4.4,i2.2,i2.2)') &
               LVT_rc%yr, LVT_rc%mo, LVT_rc%da
          write(unit=cdate3,fmt='(i2.2,i2.2)') &
               LVT_rc%hr,LVT_rc%mn

          fname_mean = trim(LVT_rc%statsodir)//&
               '/PS.AFWA_SC.'//trim(LVT_rc%security_class)//&
               '_DI.'//trim(LVT_rc%distribution_class)//&
               '_DC.'//trim(LVT_rc%data_category)//&
               '_GP.LIS_GR.C0P25DEG_AR.'//&
               trim(LVT_rc%area_of_data)//&
               '_PA.03-HR-SUM_DD.'//&
               trim(cdate2)//'_DT.'//trim(cdate3)//'_DF.GR1'

          fname_ssdev = trim(LVT_rc%statsodir)//&
               '/PS.AFWA_SC.'//trim(LVT_rc%security_class)//&
               '_DI.'//trim(LVT_rc%distribution_class)//&
               '_DC.'//trim(LVT_rc%data_category)//&
               '_GP.LIS_GR.C0P25DEG_AR.'//&
               trim(LVT_rc%area_of_data)//&
               '_PA.03-HR-SUM_DD.'//&
               trim(cdate2)//'_DT.'//trim(cdate3)//'_DF_SSDEV.GR1'


          ! Setup of GRIB-1 and GRIB-2 Metadata Section
          
          ! toplev is the depth of the top of each soil layer
          ! botlev is the depth of the bottom of each soil layer
          toplev(1) = 0.0
          botlev(1) = lyrthk(1)
          
          ! determine bounding levels for each soil moisture layer
          do i = 2, nsoillayers
             toplev(i) = toplev(i-1) + lyrthk(i-1)
             botlev(i) = botlev(i-1) + lyrthk(i)
          enddo
          !hardcoded to zero for now
          depscale = 0

          ! Set values for non layered fields (Fluxes, Sfc Fields, etc.)
          toplev0 = 0
          botlev0 = 0
          
          yr = LVT_rc%yr
          mo = LVT_rc%mo
          da = LVT_rc%da
          hr = LVT_rc%hr
          mn = LVT_rc%mn
          ss = LVT_rc%ss
          
          call LVT_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,-1*LVT_rc%statswriteint)
          
          if(LVT_rc%statswriteint .GT. 0) then
             time_unit = 254     ! seconds
             time_curr = 0
             time_past = LVT_rc%statswriteint
          endif
          if(LVT_rc%statswriteint .GE. 60) then
             time_unit = 0      ! minutes       
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 60)
          endif
          if(LVT_rc%statswriteint .GE. 3600) then
             time_unit = 1    ! hours
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 3600)
          endif
          if(LVT_rc%statswriteint .GE. 86400) then
             time_unit = 2   ! days
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 86400)
          endif
          
          !time_past: from LVT_grib1_finalize
          !time_P1 (Negative Time Unit for avg, or 0 for analysis) 
          !According to the in-line comments, time_past must be negative or 0.
          !Here we are setting it to a positive value.  This produces bad output.
          !Setting it to a negative value also produces bad output.
          !So I am resetting it to zero.  This produces output that matches
          !the binary output.
          !    time_past=0

          call grib_open_file(ftn_mean,fname_mean,'w',iret)
          call LVT_verify(iret, 'failed to open grib file '//trim(fname_mean))

          call grib_open_file(ftn_ssdev,fname_ssdev,'w',iret)
          call LVT_verify(iret, 'failed to open grib file '//trim(fname_ssdev))

       elseif(LVT_rc%lvt_out_format.eq."grib2") then  
          write(unit=cdate2,fmt='(i4.4,i2.2,i2.2)') &
               LVT_rc%yr, LVT_rc%mo, LVT_rc%da
          write(unit=cdate3,fmt='(i2.2,i2.2)') &
               LVT_rc%hr,LVT_rc%mn

          fname_mean = trim(LVT_rc%statsodir)//&
               '/PS.AFWA_SC.'//trim(LVT_rc%security_class)//&
               '_DI.'//trim(LVT_rc%distribution_class)//&
               '_DC.'//trim(LVT_rc%data_category)//&
               '_GP.LIS_GR.C0P25DEG_AR.'//&
               trim(LVT_rc%area_of_data)//&
               '_PA.03-HR-SUM_DD.'//&
               trim(cdate2)//'_DT.'//trim(cdate3)//'_DF.GR2'

          fname_ssdev = trim(LVT_rc%statsodir)//&
               '/PS.AFWA_SC.'//trim(LVT_rc%security_class)//&
               '_DI.'//trim(LVT_rc%distribution_class)//&
               '_DC.'//trim(LVT_rc%data_category)//&
               '_GP.LIS_GR.C0P25DEG_AR.'//&
               trim(LVT_rc%area_of_data)//&
               '_PA.03-HR-SUM_DD.'//&
               trim(cdate2)//'_DT.'//trim(cdate3)//'_DF_SSDEV.GR2'


          ! Setup of GRIB-1 and GRIB-2 Metadata Section
          
          ! toplev is the depth of the top of each soil layer
          ! botlev is the depth of the bottom of each soil layer
          toplev(1) = 0.0
          botlev(1) = lyrthk(1)
          
          ! determine bounding levels for each soil moisture layer
          do i = 2, nsoillayers
             toplev(i) = toplev(i-1) + lyrthk(i-1)
             botlev(i) = botlev(i-1) + lyrthk(i)
          enddo
          !hardcoded to zero for now
          depscale = 0

          ! Set values for non layered fields (Fluxes, Sfc Fields, etc.)
          toplev0 = 0
          botlev0 = 0
          
          yr = LVT_rc%yr
          mo = LVT_rc%mo
          da = LVT_rc%da
          hr = LVT_rc%hr
          mn = LVT_rc%mn
          ss = LVT_rc%ss
          
          call LVT_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,-1*LVT_rc%statswriteint)
          
          if(LVT_rc%statswriteint .GT. 0) then
             time_unit = 254     ! seconds
             time_curr = 0
             time_past = LVT_rc%statswriteint
          endif
          if(LVT_rc%statswriteint .GE. 60) then
             time_unit = 0      ! minutes       
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 60)
          endif
          if(LVT_rc%statswriteint .GE. 3600) then
             time_unit = 1    ! hours
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 3600)
          endif
          if(LVT_rc%statswriteint .GE. 86400) then
             time_unit = 2   ! days
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 86400)
          endif
          
          !time_past: from LVT_grib1_finalize
          !time_P1 (Negative Time Unit for avg, or 0 for analysis) 
          !According to the in-line comments, time_past must be negative or 0.
          !Here we are setting it to a positive value.  This produces bad output.
          !Setting it to a negative value also produces bad output.
          !So I am resetting it to zero.  This produces output that matches
          !the binary output.
          !    time_past=0

          call grib_open_file(ftn_mean,fname_mean,'w',iret)
          call LVT_verify(iret, 'failed to open grib file '//trim(fname_mean))

          call grib_open_file(ftn_ssdev,fname_ssdev,'w',iret)
          call LVT_verify(iret, 'failed to open grib file '//trim(fname_ssdev))

       elseif(LVT_rc%lvt_out_format.eq."netcdf") then  

          call date_and_time(date,stime,zone,values)

          write(unit=cdate2,fmt='(i4.4,i2.2,i2.2)') &
               LVT_rc%yr, LVT_rc%mo, LVT_rc%da
          write(unit=cdate3,fmt='(i2.2,i2.2)') &
               LVT_rc%hr,LVT_rc%mn

          fname_mean = trim(LVT_rc%statsodir)//&
               '/PS.AFWA_SC.'//trim(LVT_rc%security_class)//&
               '_DI.'//trim(LVT_rc%distribution_class)//&
               '_DC.'//trim(LVT_rc%data_category)//&
               '_GP.LIS_GR.C0P25DEG_AR.'//&
               trim(LVT_rc%area_of_data)//&
               '_PA.03-HR-SUM_DD.'//&
               trim(cdate2)//'_DT.'//trim(cdate3)//'_DF.nc'

          fname_ssdev = trim(LVT_rc%statsodir)//&
               '/PS.AFWA_SC.'//trim(LVT_rc%security_class)//&
               '_DI.'//trim(LVT_rc%distribution_class)//&
               '_DC.'//trim(LVT_rc%data_category)//&
               '_GP.LIS_GR.C0P25DEG_AR.'//&
               trim(LVT_rc%area_of_data)//&
               '_PA.03-HR-SUM_DD.'//&
               trim(cdate2)//'_DT.'//trim(cdate3)//'_DF_SSDEV.nc'


          ! Setup of GRIB-1 and GRIB-2 Metadata Section
          
          ! toplev is the depth of the top of each soil layer
          ! botlev is the depth of the bottom of each soil layer
          toplev(1) = 0.0
          botlev(1) = lyrthk(1)
          
          ! determine bounding levels for each soil moisture layer
          do i = 2, nsoillayers
             toplev(i) = toplev(i-1) + lyrthk(i-1)
             botlev(i) = botlev(i-1) + lyrthk(i)
          enddo
          !hardcoded to zero for now
          depscale = 0

          ! Set values for non layered fields (Fluxes, Sfc Fields, etc.)
          toplev0 = 0
          botlev0 = 0
          
          yr = LVT_rc%yr
          mo = LVT_rc%mo
          da = LVT_rc%da
          hr = LVT_rc%hr
          mn = LVT_rc%mn
          ss = LVT_rc%ss
          
          call LVT_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,-1*LVT_rc%statswriteint)
          
          if(LVT_rc%statswriteint .GT. 0) then
             time_unit = 254     ! seconds
             time_curr = 0
             time_past = LVT_rc%statswriteint
          endif
          if(LVT_rc%statswriteint .GE. 60) then
             time_unit = 0      ! minutes       
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 60)
          endif
          if(LVT_rc%statswriteint .GE. 3600) then
             time_unit = 1    ! hours
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 3600)
          endif
          if(LVT_rc%statswriteint .GE. 86400) then
             time_unit = 2   ! days
             time_curr = 0
             time_past = (LVT_rc%statswriteint / 86400)
          endif
          
          !time_past: from LVT_grib1_finalize
          !time_P1 (Negative Time Unit for avg, or 0 for analysis) 
          !According to the in-line comments, time_past must be negative or 0.
          !Here we are setting it to a positive value.  This produces bad output.
          !Setting it to a negative value also produces bad output.
          !So I am resetting it to zero.  This produces output that matches
          !the binary output.
          !    time_past=0



          shuffle = NETCDF_shuffle
          deflate = NETCDF_deflate
          deflate_level =NETCDF_deflate_level

          xlat%short_name = "latitude"
          xlat%long_name = "latitude"
          xlat%standard_name = "latitude"
          xlat%units = "degree_north"
          xlat%nunits = 1
          xlat%format = 'F'
          xlat%vlevels = 1
          xlat%timeAvgOpt = 0 
          xlat%startNlevs = 1
          xlat%endNlevs = 1
          allocate(xlat%value(LVT_rc%ngrid,&
               1,xlat%vlevels))
          allocate(xlat%unittypes(1))
          xlat%unittypes(1) = "degree_north"

          xlon%short_name = "longitude"
          xlon%long_name = "longitude"
          xlon%standard_name = "longitude"
          xlon%units = "degree_east"
          xlon%nunits = 1
          xlon%format = 'F'
          xlon%vlevels = 1
          xlon%timeAvgOpt = 0 
          xlon%startNlevs = 1
          xlon%endNlevs = 1
          allocate(xlon%value(LVT_rc%ngrid,&
               1,xlon%vlevels))
          allocate(xlon%unittypes(1))
          xlon%unittypes(1) = "degree_east"

#if (defined USE_NETCDF4)
          iret = nf90_create(path=trim(fname_mean), cmode =nf90_hdf5, &
               ncid = ftn_mean)
          call LVT_verify(iret, 'failed to open grib file '//trim(fname_mean))

          iret = nf90_create(path=trim(fname_ssdev), cmode =nf90_hdf5, &
               ncid = ftn_ssdev)
          call LVT_verify(iret, 'failed to open grib file '//trim(fname_ssdev))
#endif
#if (defined USE_NETCDF3)
          iret = nf90_create(path=trim(fname_mean), cmode =nf90_clobber, &
               ncid = ftn_mean)
          call LVT_verify(iret, 'failed to open grib file '//trim(fname_mean))

          iret = nf90_create(path=trim(fname_ssdev), cmode =nf90_clobber, &
               ncid = ftn_ssdev)
          call LVT_verify(iret, 'failed to open grib file '//trim(fname_ssdev))
#endif
          !Headers
          call LVT_verify(nf90_def_dim(ftn_mean,'east_west',LVT_rc%gnc,dimID(1)))
          call LVT_verify(nf90_def_dim(ftn_mean,'north_south',LVT_rc%gnr,dimID(2)))

          call LVT_verify(nf90_def_dim(ftn_mean,'time',1,tdimID))
          call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"missing_value",&
               LVT_rc%udef))
          
          call LVT_verify(nf90_def_var(ftn_mean,&
               trim(xlat%short_name),&
               nf90_float,&
               dimids = dimID(1:2), varID=xlatID))
#if(defined USE_NETCDF4) 
          call LVT_verify(nf90_def_var_deflate(ftn_mean,&
               xlatID,&
               shuffle,deflate,deflate_level))
#endif
          call LVT_verify(nf90_def_var(ftn_mean,&
               trim(xlon%short_name),&
               nf90_float,&
               dimids = dimID(1:2), varID=xlonID))
#if(defined USE_NETCDF4) 
          call LVT_verify(nf90_def_var_deflate(ftn_mean,&
               xlonID,&
               shuffle,deflate,deflate_level))
#endif
          call LVT_verify(nf90_put_att(ftn_mean,xlatID,&
               "units",trim(xlat%units)))
          call LVT_verify(nf90_put_att(ftn_mean,xlatID,&
               "standard_name",trim(xlat%standard_name)))
          call LVT_verify(nf90_put_att(ftn_mean,xlatID,&
               "long_name",trim(xlat%long_name)))
          call LVT_verify(nf90_put_att(ftn_mean,xlatID,&
               "scale_factor",1.0))
          call LVT_verify(nf90_put_att(ftn_mean,xlatID,&
               "add_offset",0.0))
          call LVT_verify(nf90_put_att(ftn_mean,xlatID,&
               "missing_value",LVT_rc%udef))
          call LVT_verify(nf90_put_att(ftn_mean,xlatID,&
               "_FillValue",LVT_rc%udef))

          call LVT_verify(nf90_put_att(ftn_mean,xlonID,&
               "units",trim(xlon%units)))
          call LVT_verify(nf90_put_att(ftn_mean,xlonID,&
               "standard_name",trim(xlon%standard_name)))
          call LVT_verify(nf90_put_att(ftn_mean,xlonID,&
               "long_name",trim(xlon%long_name)))
          call LVT_verify(nf90_put_att(ftn_mean,xlonID,&
               "scale_factor",1.0))
          call LVT_verify(nf90_put_att(ftn_mean,xlonID,&
               "add_offset",0.0))
          call LVT_verify(nf90_put_att(ftn_mean,xlonID,&
               "missing_value",LVT_rc%udef))
          call LVT_verify(nf90_put_att(ftn_mean,xlonID,&
               "_FillValue",LVT_rc%udef))

!define time field
          call LVT_verify(nf90_def_var(ftn_mean,'time',&
               nf90_float,dimids = tdimID,varID=xtimeID))
          write(xtime_units,200) LVT_rc%yr, LVT_rc%mo, LVT_rc%da, &
               LVT_rc%hr, LVT_rc%mn, LVT_rc%ss
200       format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', &
               I2.2,':',I2.2)
          write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') &
               LVT_rc%yr, LVT_rc%mo, LVT_rc%da
          write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') &
               LVT_rc%hr, LVT_rc%mn, LVT_rc%ss
          write(xtime_timeInc, fmt='(I20)') &
               LVT_rc%ts
          
          call LVT_verify(nf90_put_att(ftn_mean,xtimeID,&
               "units",trim(xtime_units)))
          call LVT_verify(nf90_put_att(ftn_mean,xtimeID,&
               "long_name","time"))
          call LVT_verify(nf90_put_att(ftn_mean,xtimeID,&
               "time_increment",trim(adjustl(xtime_timeInc))))
          call LVT_verify(nf90_put_att(ftn_mean,xtimeID,&
               "begin_date",xtime_begin_date))
          call LVT_verify(nf90_put_att(ftn_mean,xtimeID,&
               "begin_time",xtime_begin_time))

          call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"title", &
               "LVT land surface analysis output"))
          call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"institution", &
               trim(LVT_rc%institution)))
          call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"history", &
               "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
               date(7:8)//"T"//stime(1:2)//":"//stime(3:4)//":"//stime(5:10)))
          call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"references", &
               "Kumar_etal_GMD_2012"))
          call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"comment", &
               "website: http://lis.gsfc.nasa.gov/"))

          !grid information
          if(trim(LVT_rc%domain).eq."latlon") then !latlon
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"MAP_PROJECTION", &
                  "EQUIDISTANT CYLINDRICAL"))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(9)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(10)))       
          elseif(trim(LVT_rc%domain).eq."mercator") then 
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"MAP_PROJECTION", &
                  "MERCATOR"))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(9)))
          elseif(trim(LVT_rc%domain).eq."lambert") then !lambert conformal
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"MAP_PROJECTION", &
                  "LAMBERT CONFORMAL"))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"TRUELAT2", &
                  LVT_rc%gridDesc(7)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(9)))
             
          elseif(trim(LVT_rc%domain).eq."polar") then ! polar stereographic
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"MAP_PROJECTION", &
                  "POLAR STEREOGRAPHIC"))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"ORIENT", &
                  LVT_rc%gridDesc(7)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn_mean,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(9)))
          endif


          !Headers
          call LVT_verify(nf90_def_dim(ftn_ssdev,'east_west',LVT_rc%gnc,dimID(1)))
          call LVT_verify(nf90_def_dim(ftn_ssdev,'north_south',LVT_rc%gnr,dimID(2)))

          call LVT_verify(nf90_def_dim(ftn_ssdev,'time',1,tdimID))
          call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"missing_value",&
               LVT_rc%udef))
          
          call LVT_verify(nf90_def_var(ftn_ssdev,&
               trim(xlat%short_name),&
               nf90_float,&
               dimids = dimID(1:2), varID=xlat_ss_ID))
#if(defined USE_NETCDF4) 
          call LVT_verify(nf90_def_var_deflate(ftn_ssdev,&
               xlat_ss_ID,&
               shuffle,deflate,deflate_level))
#endif
          call LVT_verify(nf90_def_var(ftn_ssdev,&
               trim(xlon%short_name),&
               nf90_float,&
               dimids = dimID(1:2), varID=xlon_ss_ID))
#if(defined USE_NETCDF4) 
          call LVT_verify(nf90_def_var_deflate(ftn_ssdev,&
               xlon_ss_ID,&
               shuffle,deflate,deflate_level))
#endif
          call LVT_verify(nf90_put_att(ftn_ssdev,xlat_ss_ID,&
               "units",trim(xlat%units)))
          call LVT_verify(nf90_put_att(ftn_ssdev,xlat_ss_ID,&
               "standard_name",trim(xlat%standard_name)))
          call LVT_verify(nf90_put_att(ftn_ssdev,xlat_ss_ID,&
               "long_name",trim(xlat%long_name)))
          call LVT_verify(nf90_put_att(ftn_ssdev,xlat_ss_ID,&
               "scale_factor",1.0))
          call LVT_verify(nf90_put_att(ftn_ssdev,xlat_ss_ID,&
               "add_offset",0.0))
          call LVT_verify(nf90_put_att(ftn_ssdev,xlat_ss_ID,&
               "missing_value",LVT_rc%udef))
          call LVT_verify(nf90_put_att(ftn_ssdev,xlat_ss_ID,&
               "_FillValue",LVT_rc%udef))

          call LVT_verify(nf90_put_att(ftn_ssdev,xlon_ss_ID,&
               "units",trim(xlon%units)))
          call LVT_verify(nf90_put_att(ftn_ssdev,xlon_ss_ID,&
               "standard_name",trim(xlon%standard_name)))
          call LVT_verify(nf90_put_att(ftn_ssdev,xlon_ss_ID,&
               "long_name",trim(xlon%long_name)))
          call LVT_verify(nf90_put_att(ftn_ssdev,xlon_ss_ID,&
               "scale_factor",1.0))
          call LVT_verify(nf90_put_att(ftn_ssdev,xlon_ss_ID,&
               "add_offset",0.0))
          call LVT_verify(nf90_put_att(ftn_ssdev,xlon_ss_ID,&
               "missing_value",LVT_rc%udef))
          call LVT_verify(nf90_put_att(ftn_ssdev,xlon_ss_ID,&
               "_FillValue",LVT_rc%udef))

!define time field
          call LVT_verify(nf90_def_var(ftn_ssdev,'time',&
               nf90_float,dimids = tdimID,varID=xtime_ss_ID))
          write(xtime_units,201) LVT_rc%yr, LVT_rc%mo, LVT_rc%da, &
               LVT_rc%hr, LVT_rc%mn, LVT_rc%ss
201       format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', &
               I2.2,':',I2.2)
          write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') &
               LVT_rc%yr, LVT_rc%mo, LVT_rc%da
          write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') &
               LVT_rc%hr, LVT_rc%mn, LVT_rc%ss
          write(xtime_timeInc, fmt='(I20)') &
               LVT_rc%ts
          
          call LVT_verify(nf90_put_att(ftn_ssdev,xtime_ss_ID,&
               "units",trim(xtime_units)))
          call LVT_verify(nf90_put_att(ftn_ssdev,xtime_ss_ID,&
               "long_name","time"))
          call LVT_verify(nf90_put_att(ftn_ssdev,xtime_ss_ID,&
               "time_increment",trim(adjustl(xtime_timeInc))))
          call LVT_verify(nf90_put_att(ftn_ssdev,xtime_ss_ID,&
               "begin_date",xtime_begin_date))
          call LVT_verify(nf90_put_att(ftn_ssdev,xtime_ss_ID,&
               "begin_time",xtime_begin_time))

          call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"title", &
               "LVT land surface analysis output"))
          call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"institution", &
               trim(LVT_rc%institution)))
          call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"history", &
               "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
               date(7:8)//"T"//stime(1:2)//":"//stime(3:4)//":"//stime(5:10)))
          call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"references", &
               "Kumar_etal_GMD_2012"))
          call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"comment", &
               "website: http://lis.gsfc.nasa.gov/"))

          !grid information
          if(trim(LVT_rc%domain).eq."latlon") then !latlon
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"MAP_PROJECTION", &
                  "EQUIDISTANT CYLINDRICAL"))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(9)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(10)))       
          elseif(trim(LVT_rc%domain).eq."mercator") then 
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"MAP_PROJECTION", &
                  "MERCATOR"))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(9)))
          elseif(trim(LVT_rc%domain).eq."lambert") then !lambert conformal
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"MAP_PROJECTION", &
                  "LAMBERT CONFORMAL"))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"TRUELAT2", &
                  LVT_rc%gridDesc(7)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(9)))
             
          elseif(trim(LVT_rc%domain).eq."polar") then ! polar stereographic
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"MAP_PROJECTION", &
                  "POLAR STEREOGRAPHIC"))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"ORIENT", &
                  LVT_rc%gridDesc(7)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn_ssdev,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(9)))
          endif

          dataEntry => LVT_histData%head_ds1_list
          
          do while(associated(dataEntry))
             !reset the pointers to the head of the linked list
             if(LVT_LIS_rc(1)%anlys_data_class.eq."LSM") then 
                lisdataEntry => LVT_LISoutput(1)%head_lsm_list
             elseif(LVT_LIS_rc(1)%anlys_data_class.eq."Routing") then 
                lisdataEntry => LVT_LISoutput(1)%head_routing_list
             elseif(LVT_LIS_rc(1)%anlys_data_class.eq."RTM") then 
                lisdataEntry => LVT_LISoutput(1)%head_rtm_list
             elseif(LVT_LIS_rc(1)%anlys_data_class.eq."Irrigation") then 
                lisdataEntry => LVT_LISoutput(1)%head_irrig_list
             endif
             do while(associated(lisdataEntry)) 
                if(lisdataEntry%short_name.eq.dataEntry%short_name) then

                   call defineNETCDFheaderVar(ftn_mean,dimID, lisdataEntry)  

                   call defineNETCDFheaderVar_ss(ftn_ssdev,dimID, lisdataEntry)  

                endif
                lisdataEntry => lisdataEntry%next
             enddo
             dataEntry => dataEntry%next
          enddo

          if(LVT_rc%processHYCOM.eq.1) then 

             call LVT_verify(nf90_def_var(ftn_mean,&
                  trim(LVT_histData%watertemp%short_name),&
                  nf90_float,&
                  dimids = dimID(1:2), &
                  varID=LVT_histData%watertemp%varId_def),&
                  'nf90_def_var for '//&
                  trim(LVT_histData%watertemp%short_name)//&
                  'failed in defineNETCDFheadervar')                     
#if(defined USE_NETCDF4)
             call LVT_verify(nf90_def_var_deflate(ftn_mean,&
                  LVT_histData%watertemp%varId_def,&
                  shuffle, deflate, deflate_level),&
                  'nf90_def_var_deflate for '//&
                  trim(LVT_histData%watertemp%short_name)//&
                  'failed in defineNETCDFheadervar')    
#endif             
          endif
          
          call LVT_verify(nf90_enddef(ftn_mean))
          call LVT_verify(nf90_put_var(ftn_mean,xtimeID,0.0))

          call LVT_verify(nf90_enddef(ftn_ssdev))
          call LVT_verify(nf90_put_var(ftn_ssdev,xtime_ss_ID,0.0))
      
          lat = LVT_rc%udef
          lon = LVT_rc%udef

          do r=1,LVT_rc%gnr
             do c=1,LVT_rc%gnc
                gid = LVT_domain%gindex(c,r)
                if(gid.ne.-1) then 
                   lat(c,r) = LVT_domain%grid(gid)%lat
                   lon(c,r) = LVT_domain%grid(gid)%lon
                endif
             enddo
          enddo
          call LVT_verify(nf90_put_var(ftn_mean,xlatID, &
               lat, (/1,1/),&
               (/LVT_rc%gnc,LVT_rc%gnr/)),&
               'nf90_put_var failed for lat')
          call LVT_verify(nf90_put_var(ftn_mean,xlonID, &
               lon, (/1,1/),&
               (/LVT_rc%gnc,LVT_rc%gnr/)),&
               'nf90_put_var failed for lon')

          call LVT_verify(nf90_put_var(ftn_ssdev,xlat_ss_ID, &
               lat, (/1,1/),&
               (/LVT_rc%gnc,LVT_rc%gnr/)),&
               'nf90_put_var failed for lat')
          call LVT_verify(nf90_put_var(ftn_ssdev,xlon_ss_ID, &
               lon, (/1,1/),&
               (/LVT_rc%gnc,LVT_rc%gnr/)),&
               'nf90_put_var failed for lon')
       endif
       
       dataEntry => LVT_histData%head_ds1_list
          
       do while(associated(dataEntry))
!reset the pointers to the head of the linked list
          if(LVT_LIS_rc(1)%anlys_data_class.eq."LSM") then 
             lisdataEntry => LVT_LISoutput(1)%head_lsm_list
          elseif(LVT_LIS_rc(1)%anlys_data_class.eq."Routing") then 
             lisdataEntry => LVT_LISoutput(1)%head_routing_list
          elseif(LVT_LIS_rc(1)%anlys_data_class.eq."RTM") then 
             lisdataEntry => LVT_LISoutput(1)%head_rtm_list
          elseif(LVT_LIS_rc(1)%anlys_data_class.eq."Irrigation") then 
             lisdataEntry => LVT_LISoutput(1)%head_irrig_list
          endif
          do while(associated(lisdataEntry)) 
             
             if(lisdataEntry%short_name.eq.dataEntry%short_name) then
                
                ! Set timerange indicator equal to 133 for AFWA's specifications
                ! for surface runoff, baseflow, and total precipitation
                ! to make the LIS-7 output match the LIS-6 style. - dmm
                      
                if(lisdataEntry%timeAvgOpt.eq.0) then !instantaneous only 
                   stepType = "instant"
                   timeRange = 1
                   pdTemplate = 0
                elseif(lisdataEntry%timeAvgOpt.eq.1) then !time averaged only 
                   stepType = "avg"
                   
                   timeRange = 7
                   pdTemplate = 8
                elseif(lisdataEntry%timeAvgOpt.eq.3) then !accumulated
                   stepType = "accum"
                   timeRange = 4
                   pdTemplate = 8
                elseif(lisdataEntry%timeAvgOpt.eq.2) then !pick the time avg field
                   stepType = "avg"
                   timeRange = 7
                   pdTemplate = 8
                endif
                if ((lisdataEntry%index.eq.LVT_LIS_MOC_QS(1)).or.   &
                     (lisdataEntry%index.eq.LVT_LIS_MOC_QSB(1)).or.   &
                     (lisdataEntry%index.eq.LVT_LIS_MOC_TOTALPRECIP(1))) then
                   timeRange = 133
                endif
                
                do k=1,dataEntry%vlevels
                   gtmp1_1d = 0.0
                   ngtmp1_1d = 0
                   gtmp1_ss = 0.0

                   do m=1,LVT_rc%nensem
                      do r=1,LVT_rc%lnr
                         do c=1,LVT_rc%lnc
                            if(LVT_domain%gindex(c,r).ne.-1) then 
                               gid = LVT_domain%gindex(c,r)
                               gtmp1_1d(c+(r-1)*LVT_rc%lnc) = &
                                    gtmp1_1d(c+(r-1)*LVT_rc%lnc) + &
                                    dataEntry%value(gid,m,k) 
                               ngtmp1_1d(c+(r-1)*LVT_rc%lnc) = &
                                    ngtmp1_1d(c+(r-1)*LVT_rc%lnc) + 1

                               gtmp1_ss(c+(r-1)*LVT_rc%lnc) = &
                                    gtmp1_ss(c+(r-1)*LVT_rc%lnc) + & 
                                    dataEntry%value(gid,m,k)**2
                               
                            endif
                         enddo
                      enddo
                   enddo

                   do r=1,LVT_rc%lnr
                      do c=1,LVT_rc%lnc
                         if(ngtmp1_1d(c+(r-1)*LVT_rc%lnc).gt.0) then
                            gtmp1_1d(c+(r-1)*LVT_rc%lnc) = &
                                 gtmp1_1d(c+(r-1)*LVT_rc%lnc)/&
                                    ngtmp1_1d(c+(r-1)*LVT_rc%lnc)
                            variance = gtmp1_ss(c+(r-1)*LVT_rc%lnc)/&
                                 ngtmp1_1d(c+(r-1)*LVT_rc%lnc) - & 
                                 gtmp1_1d(c+(r-1)*LVT_rc%lnc)**2

                            if(variance.gt.0) then 
                               gtmp1_ss(c+(r-1)*LVT_rc%lnc) = & 
                                    sqrt(variance)
                            else
                               gtmp1_ss(c+(r-1)*LVT_rc%lnc) = LVT_rc%udef
                            endif                            
                         else
                            gtmp1_1d(c+(r-1)*LVT_rc%lnc) = LVT_rc%udef
                            gtmp1_ss(c+(r-1)*LVT_rc%lnc) = LVT_rc%udef
                         endif
                      enddo
                   enddo
!                   open(100,file='test.bin',form='unformatted')
!                   write(100) gtmp1_1d
                   if(.not.((dataEntry%short_name.eq."Landcover").or.&
                        (dataEntry%short_name.eq."Landmask").or.&
                        (dataEntry%short_name.eq."Soiltype").or.&
                        (dataEntry%short_name.eq."Greenness").or.&
                        (dataEntry%short_name.eq."Tair_f_min").or.&
                        (dataEntry%short_name.eq."Tair_f_max").or.&
                        (dataEntry%short_name.eq."TotalPrecip").or.&
                        (dataEntry%short_name.eq."Wind_f").or.&
                        (dataEntry%short_name.eq."Tair_f").or.&
                        (dataEntry%short_name.eq."SWdown_f").or.&
                        (dataEntry%short_name.eq."LWdown_f").or.&
                        (dataEntry%short_name.eq."RHMin"))) then 
                      if(LVT_rc%applyNoiseReductionFilter.eq.1) then 
                         call applyNoiseReductionFilter(gtmp1_1d)
                      endif
                   endif

!                   write(100) gtmp1_1d
!                   close(100)
!                   stop
                   
                   
                   if(LVT_rc%lvt_out_format.eq."grib2") then 
                      
                      call writeSingleGrib2Var(ftn_mean,&
                           gtmp1_1d,&
                           lisdataentry%varid_def,&
                           lisdataentry%gribSF,&
                           lisdataentry%gribSfc,&
                           lisdataentry%gribLvl,&
                           lisdataentry%gribDis,&
                           lisdataentry%gribCat,&
                           pdTemplate,&
                           stepType,&
                           time_unit,&
                           time_past,&
                           time_curr,&
                           timeRange,&
                           k,&
                           toplev(k:k),&
                           botlev(k:k),&
                           depscale(k:k))

                      call writeSingleGrib2Var(ftn_ssdev,&
                           gtmp1_ss,&
                           lisdataentry%varid_def,&
                           lisdataentry%gribSF,&
                           lisdataentry%gribSfc,&
                           lisdataentry%gribLvl,&
                           lisdataentry%gribDis,&
                           lisdataentry%gribCat,&
                           pdTemplate,&
                           stepType,&
                           time_unit,&
                           time_past,&
                           time_curr,&
                           timeRange,&
                           k,&
                           toplev(k:k),&
                           botlev(k:k),&
                           depscale(k:k))
                   elseif(LVT_rc%lvt_out_format.eq."grib1") then 
                      call writeSingleGrib1Var(ftn_mean,&
                           gtmp1_1d,&
                           lisdataentry%varid_def,&
                           lisdataentry%gribSF,&
                           lisdataentry%gribSfc,&
                           lisdataentry%gribLvl,&
                           stepType,&
                           time_unit,&
                           time_past,&
                           time_curr,&
                           timeRange,&
                           k,&
                           toplev(k:k),&
                           botlev(k:k))

                      call writeSingleGrib1Var(ftn_ssdev,&
                           gtmp1_ss,&
                           lisdataentry%varid_def,&
                           lisdataentry%gribSF,&
                           lisdataentry%gribSfc,&
                           lisdataentry%gribLvl,&
                           stepType,&
                           time_unit,&
                           time_past,&
                           time_curr,&
                           timeRange,&
                           k,&
                           toplev(k:k),&
                           botlev(k:k))

                   elseif(LVT_rc%lvt_out_format.eq."netcdf") then 
                      call writeSingleNetcdfVar(ftn_mean,&
                           gtmp1_1d,&
                           lisdataentry%varid_def,&
                           k)
                      call writeSingleNetcdfVar(ftn_ssdev,&
                           gtmp1_ss,&
                           lisdataentry%varid_ss,&
                           k)
                   endif
                   
                enddo
                exit
             endif
             lisdataEntry => lisdataEntry%next
          enddo
          dataEntry => dataEntry%next
       enddo

       call LVT_append_HYCOM_fields(ftn_mean,&
            time_unit,&
            time_past,&
            time_curr,&
            timeRange,&
            toplev(1),&
            botlev(1))
       
       if(LVT_rc%lvt_out_format.eq."grib1") then  
          call grib_close_file(ftn_mean,iret)
          call grib_close_file(ftn_ssdev,iret)
       elseif(LVT_rc%lvt_out_format.eq."grib2") then  
          call grib_close_file(ftn_mean,iret)
          call grib_close_file(ftn_ssdev,iret)
       elseif(LVT_rc%lvt_out_format.eq."netcdf") then  
          call LVT_verify(nf90_close(ftn_mean))
          call LVT_verify(nf90_close(ftn_ssdev))
       endif
    endif
       
  end subroutine LVT_writeDataStreams

!BOP
! 
! !ROUTINE: LVT_append_HYCOM_fields
! \label{LVT_append_HYCOM_fields}
!
! !INTERFACE: 
  subroutine LVT_append_HYCOM_fields(ftn_mean,time_unit, time_past, time_curr,&
       timeRange, toplev,botlev)

! 
! !DESCRIPTION: 
!  This subroutine read the water temperature fields from the HYCOM output, 
!  reprojects it to the LVT/LIS grid and appends to the grib1 formatted file. 
!
!EOP
    
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    integer                 :: ftn_mean
    integer                 :: time_unit
    integer                 :: time_past
    integer                 :: time_curr
    integer                 :: timeRange
    real                    :: toplev(1)
    real                    :: botlev(1)
    
    character*100           :: hycom_fname
    character*10            :: cdate
    logical                 :: file_exists
    integer                 :: nid,ios
    integer                 :: c,r,c1,r1,k,cindex,rindex
    integer                 :: watertid
    real                    :: watert_ip(LVT_rc%lnc*LVT_rc%lnr)
    logical*1               :: lo(LVT_rc%lnc*LVT_rc%lnr)
    logical*1               :: lb(LVT_rc%HYCOM_nc*LVT_rc%HYCOM_nr)
    real                    :: watert(LVT_rc%HYCOM_nc,LVT_rc%HYCOM_nr,1,1)
    real                    :: watert_1d(LVT_rc%HYCOM_nc*LVT_rc%HYCOM_nr)
    character*10            :: stepType
    integer                 :: varid_def
    integer                 :: gribSF, gribSfc,gribLvl,gribCat,gribDis
    real                    :: depscale(1)
    integer                 :: pdTemplate

    ! find the filename, open the file, read the field

    if(LVT_rc%processHYCOM.eq.1) then 
       
       write(unit=cdate,fmt='(i4.4,i2.2,i2.2,i2.2)') &
            LVT_rc%yr, LVT_rc%mo, LVT_rc%da, LVT_rc%hr

       hycom_fname = trim(LVT_rc%HYCOMdir)//'/'//&
           'hycom_glb_928_'//trim(cdate)//'_t012_ts3z.nc'
           
       watert  = LVT_rc%udef
       print*, hycom_fname
       inquire(file=hycom_fname,exist=file_exists)

       if(file_exists) then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
          write(LVT_logunit,*) '[INFO] Reading HYCOM data ',trim(hycom_fname)

          ios = nf90_open(path=trim(hycom_fname),mode=NF90_NOWRITE,ncid=nid)
          call LVT_verify(ios, 'Error opening file'//trim(hycom_fname))

!variable ids
          ios = nf90_inq_varid(nid, 'water_temp',watertid)
          call LVT_verify(ios, 'Error nf90_inq_varid: water_temp')
          
!values
          ios = nf90_get_var(nid,watertid, watert,&
               start=(/1,1,1,1/), count=(/LVT_rc%HYCOM_nc,LVT_rc%HYCOM_nr,1,1/))
          call LVT_verify(ios, 'Error nf90_get_var: water_temp')

          ios = nf90_close(nid)
          call LVT_verify(ios, 'Error in nf90_close')
#endif
          watert_1d = -9999.0
          lb = .false. 

          do r=1,LVT_rc%HYCOM_nr 
             do c=1,LVT_rc%HYCOM_nc
                if(watert(c,r,1,1).ne.-30000) then 
                   if(c.gt.2250) then 
                      c1 = c-2250
                      r1 = r
                   else
                      c1 = c+2250
                      r1 = r
                   endif
                   watert_1d(c1+(r1-1)*LVT_rc%HYCOM_nc) = watert(c,r,1,1)*0.001+20.0
                   lb(c1+(r1-1)*LVT_rc%HYCOM_nc) = .true.
                   
                endif
             enddo
          enddo
                   
          call upscaleByAveraging(&
               LVT_rc%HYCOM_nc*LVT_rc%HYCOM_nr, &
               LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
               LVT_rc%HYCOM_n11, lb, &
               watert_1d, lo, watert_ip)

          varid_def = 17
          gribSF    = 1
          gribSfc   = 1
          gribLvl   = 1
          gribCat   = 0
          gribDis   = 0
          pdTemplate = 0 
          stepType  = "avg"
          
          if(LVT_rc%lvt_out_format.eq."grib2") then 
     ! add to the grib file
             call writeSingleGrib2Var(ftn_mean,&
                  watert_ip,&
                  varid_def,&
                  gribSF,&
                  gribSfc,&
                  gribLvl,&
                  gribDis,&
                  gribCat,&
                  pdTemplate,&
                  stepType,&
                  time_unit,&
                  time_past,&
                  time_curr,&
                  timeRange,&
                  1,&
                  toplev(1),&
                  botlev(1),&
                  depscale(1))
             
          elseif(LVT_rc%lvt_out_format.eq."grib1") then 
             call writeSingleGrib1Var(ftn_mean,&
                  watert_ip,&
                  varid_def,&
                  gribSF,&
                  gribSfc,&
                  gribLvl,&
                  stepType,&
                  time_unit,&
                  time_past,&
                  time_curr,&
                  timeRange,&
                  1,&
                  toplev(1),&
                  botlev(1))
          elseif(LVT_rc%lvt_out_format.eq."netcdf") then 
             call writeSingleNetcdfVar(ftn_mean,&
                  watert_ip,&
                  LVT_histData%watertemp%varId_def,&
                  1)
             
          endif
       endif
    endif

  end subroutine LVT_append_HYCOM_fields

  subroutine applyNoiseReductionFilter(gvar)
    
    real   :: gvar(LVT_rc%lnc*LVT_rc%lnr)
    real   :: gtmp(LVT_rc%lnc*LVT_rc%lnr)
    
    integer :: c,r, c1,r1,c_s, c_e, r_s, r_e
    real    :: avg_val
    real    :: navg_val
    real    :: sigma,wt

    gtmp = LVT_rc%udef
    

    if(LVT_rc%smoothingFilterType.eq."box filter") then 
       do r=1,LVT_rc%lnr
          do c=1,LVT_rc%lnc
             
             c_s = max(1,c-2)
             c_e = min(LVT_rc%lnc,c+2)
             r_s = max(1,r-2)
             r_e = min(LVT_rc%lnr,r+2)
             
             avg_val = 0.0
             navg_val = 0
             do c1=c_s, c_e
                do r1=r_s,r_e
                   if(gvar(c1+(r1-1)*LVT_rc%lnc).ne.LVT_rc%udef) then 
                      avg_val = avg_val + gvar(c1+(r1-1)*LVT_rc%lnc)
                      navg_val = navg_val + 1
                   endif
                enddo
             enddo
             if(navg_val.gt.0) then 
                avg_val = avg_val/navg_val
             else
                avg_val = LVT_rc%udef
             endif
             
             gtmp(c+(r-1)*LVT_rc%lnc) = avg_val
             
          enddo
       enddo
       

    elseif(LVT_rc%smoothingFilterType.eq."gaussian filter") then 
       sigma = 2.0
       do r=1,LVT_rc%lnr
          do c=1,LVT_rc%lnc
             
             c_s = max(1,c-5)
             c_e = min(LVT_rc%lnc,c+5)
             r_s = max(1,r-5)
             r_e = min(LVT_rc%lnr,r+5)
             
             avg_val = 0.0
             navg_val = 0
             do c1=c_s, c_e
                do r1=r_s,r_e
                   if(gvar(c1+(r1-1)*LVT_rc%lnc).ne.LVT_rc%udef) then 
                      wt = exp(-((c1-c)**2+(r1-r)**2)/(2*sigma**2))/&
                           (2*LVT_CONST_PI*sigma**2)
                      avg_val = avg_val + wt*gvar(c1+(r1-1)*LVT_rc%lnc)
                      navg_val = navg_val + wt
                   endif
                enddo
             enddo
             if(navg_val.gt.0) then 
                if(gvar(c+(r-1)*LVT_rc%lnc).ne.LVT_rc%udef) then 
                   avg_val = avg_val/navg_val
                else
                   avg_val = LVT_rc%udef
                endif
             else
                avg_val = LVT_rc%udef
             endif
             
             gtmp(c+(r-1)*LVT_rc%lnc) = avg_val
             
          enddo
       enddo       

    endif

    gvar = gtmp

  end subroutine applyNoiseReductionFilter
!BOP
! 
! !ROUTINE: writeSingleGrib1Var
! \label{writeSingleGrib1Var}
!
! !INTERFACE: 
  subroutine writeSingleGrib1Var(ftn,gtmp,gribId,gribSF,gribSfc,gribLvl,&
       sType, time_unit, time_p1, time_p2, &
       timeRange,k,toplev,botlev)
!
! !DESCRIPTION: 
!  This subroutine writes a single variable to a grib file
! 
!EOP
    
    integer                       :: ftn
    real                          :: gtmp(LVT_rc%lnc*LVT_rc%lnr)
    integer, intent(in)           :: gribid
    integer, intent(in)           :: gribSF
    integer, intent(in)           :: gribSfc
    integer, intent(in)           :: gribLvl
    character(len=*), intent(in)  :: sType
    integer, intent(in)           :: timeRange
    integer, intent(in)           :: time_unit
    integer, intent(in)           :: time_p1
    integer, intent(in)           :: time_p2
    integer, intent(in)           :: k
    real                          :: toplev(1)
    real                          :: botlev(1)
 

    character*8                   :: date
    integer                       :: yr1, mo1,da1,hr1,mn1
    integer                       :: idate,idate1
    real                          :: lat_ur, lon_ur
    real                          :: lat_ll, lon_ll
    integer                       :: igrib,iret
    integer                       :: decimalPrecision,gribSFtemp
    character*100                 :: message(20)

    ! Note passing string of defined points only to output
    ! because bitmap in GRIB-1 file will fill in the rest
    
    call grib_new_from_template(igrib,"GRIB1",iret)
    call LVT_verify(iret, 'grib_new_from_template failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'table2Version',LVT_rc%grib_table,iret)
    call LVT_verify(iret,'grib_set:table2version failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'generatingProcessIdentifier',LVT_rc%grib_process_id,iret)
    call LVT_verify(iret,'grib_set:generatingProcessIdentifier failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'gridDefinition',LVT_rc%grib_grid_id,iret)
    call LVT_verify(iret,'grib_set:grid ID failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'indicatorOfParameter',gribid, iret)
    call LVT_verify(iret,'grib_set:indicatorOfParameter failed in LVT_DataStreamsMod')
    
    !    call grib_set(igrib,'paramId',gribid, iret)
    !    call LVT_verify(iret,'grib_set:paramId failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'indicatorOfTypeOfLevel',gribSfc, iret)
    call LVT_verify(iret,'grib_set:indicatorOfTypeOfLevel failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'level',gribLvl, iret)
    call LVT_verify(iret,'grib_set:level failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'topLevel',toplev(1), iret)
    call LVT_verify(iret,'grib_set:topLevel failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'bottomLevel',botlev(1), iret)
    call LVT_verify(iret,'grib_set:bottomLevel failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'stepType',sType, iret)
    call LVT_verify(iret,'grib_set:stepType failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'stepUnits',time_unit, iret)
    call LVT_verify(iret,'grib_set:stepUnits failed in LVT_DataStreamsMod')

    call grib_set(igrib,'startStep',time_p1, iret)
    call LVT_verify(iret,'grib_set:startStep failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'endStep',time_p2, iret)
    call LVT_verify(iret,'grib_set:endStep failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'timeRangeIndicator',timeRange, iret)
    call LVT_verify(iret,'grib_set:timeRangeIndicator failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'swapScanningLat',1, iret)
    call LVT_verify(iret,'grib_set:swapScanningLat failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'Ni',LVT_rc%gnc,iret)
    call LVT_verify(iret, 'grib_set:Ni failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'Nj',LVT_rc%gnr,iret)
    call LVT_verify(iret, 'grib_set:Ni failed in LVT_DataStreamsMod')
    
    call ij_to_latlon(LVT_domain%lvtproj,float(LVT_rc%gnc),&
         float(LVT_rc%gnr),lat_ur,lon_ur)      
    call ij_to_latlon(LVT_domain%lvtproj,1.0, 1.0, &
         lat_ll,lon_ll)      
    
    call grib_set(igrib, 'latitudeOfFirstGridPointInDegrees',lat_ll,iret)
    call LVT_verify(iret, 'grib_set:latitudeOfFirstGridPointInDegrees failed in LVT_DataStreamsMod')
    
    call grib_set(igrib, 'longitudeOfFirstGridPointInDegrees',lon_ll,iret)
    call LVT_verify(iret, 'grib_set:longitudeOfFirstGridPointInDegrees failed in LVT_DataStreamsMod')
    
    call grib_set(igrib, 'latitudeOfLastGridPointInDegrees',lat_ur,iret)
    call LVT_verify(iret, 'grib_set:latitudeOfLastGridPointInDegrees failed in LVT_DataStreamsMod')
    
    call grib_set(igrib, 'longitudeOfLastGridPointInDegrees',lon_ur,iret)
    call LVT_verify(iret, 'grib_set:longitudeOfLastGridPointInDegrees failed in LVT_DataStreamsMod')
    
    call grib_set(igrib, 'missingValue',LVT_rc%udef,iret)
    call LVT_verify(iret, 'grib_set:missingValue failed in LVT_DataStreamsMod')
    
! Should not need to fix the "num bits" value for each parameter
! if the "decimalPrecision" (aka, "DecScale") is set properly. - dmm
!     call grib_set(igrib, 'bitsPerValue',12,iret)
!     call LVT_verify(iret, 'grib_set:bitsPerValue failed in LVT_DataStreamsMod')

! Set the "decimalPrecision" (aka, "DecScale") based on the
! gribSF (grib scale factor) set in the MODEL OUTPUT TBL. - dmm
    gribSFtemp = gribSF
    decimalPrecision = 0
    do while (gribSFtemp.ge.10)
       decimalPrecision = decimalPrecision + 1
       gribSFtemp = gribSFtemp / 10
    enddo
    call grib_set(igrib, 'decimalPrecision',decimalPrecision,iret)
    call LVT_verify(iret, 'grib_set:decimalPrecision failed in LVT_DataStreamsMod')
    
    call grib_set(igrib, 'bitmapPresent',1,iret)
    call LVT_verify(iret, 'grib_set:bitmapPresent failed in LVT_DataStreamsMod')
 
    if (LVT_rc%domain.eq."latlon") then 
       call grib_set(igrib,'gridType','regular_ll',iret)
       call LVT_verify(iret,'grib_set: gridType failed in LVT_DataStreamsMod')
       
       call grib_set(igrib,'iDirectionIncrementInDegrees',LVT_rc%gridDesc(9),iret)
       call LVT_verify(iret,'grib_set:iDirectionIncrementInDegrees failed in LVT_DataStreamsMod')
       
       call grib_set(igrib,'jDirectionIncrementInDegrees',LVT_rc%gridDesc(10),iret)
       call LVT_verify(iret,'grib_set:jDirectionIncrementInDegrees failed in LVT_DataStreamsMod')
       
    else  !Unsupported Map Projection for GRIB output
       
       message(1)='program:  LVT_DataStreamsMod'
       message(2)=' subroutine:  writevar_grib1_withstats_real'
       message(3)='  Unsupported map projection for GRIB1 output!'
       call lvt_abort(message)
       stop
       
    endif
    
    da1=LVT_rc%da
    mo1=LVT_rc%mo
    yr1=LVT_rc%yr
    
    write(unit=date,fmt='(i4.4,i2.2,i2.2)') yr1,mo1,da1
    read(date,'(I8)') idate
    
    call grib_set(igrib,'dataDate',idate,iret)
    call LVT_verify(iret, 'grib_set:dataDate failed in LVT_DataStreamsMod')
    
    hr1=LVT_rc%hr
    mn1=LVT_rc%mn
    
    write(unit=date,fmt='(i2.2,i2.2)') hr1,mn1
    read(date,'(I4)') idate1
    
    call grib_set(igrib,'dataTime',idate1,iret)
    call LVT_verify(iret, 'grib_set:dataTime failed in LVT_DataStreamsMod')
    

    call grib_set(igrib,'values',gtmp,iret)
    call LVT_verify(iret, 'grib_set:values failed in LVT_DataStreamsMod')
    
    ! Move setting of centre and subCentre to the end of the settings.
    ! The order these are written is important and will affect output. - dmm
    call grib_set(igrib,'centre',LVT_rc%grib_center_id,iret)
    call LVT_verify(iret,'grib_set:centre failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'subCentre',LVT_rc%grib_subcenter_id,iret)
    call LVT_verify(iret,'grib_set:subCentre failed in LVT_DataStreamsMod')
    
    call grib_write(igrib,ftn,iret)
    call LVT_verify(iret, 'grib_write failed in LVT_DataStreamsMod')
    
    call grib_release(igrib,iret)
    call LVT_verify(iret,'grib_release failed in LVT_DataStreamsMod')


  end subroutine writeSingleGrib1Var

!BOP
! 
! !ROUTINE: writeSingleGrib2Var
! \label{writeSingleGrib2Var}
!
! !INTERFACE: 
  subroutine writeSingleGrib2Var(ftn,gtmp,gribId,gribSF,gribSfc,gribLvl,&
       gribDis, gribCat, pdTemplate, &
       sType, time_unit, time_p1, time_p2, &
       timeRange,k,toplev,botlev,depscale)
!
! !DESCRIPTION: 
!  This subroutine writes a single variable to a grib2 file based on 
!  the implementation by Hiroko Kato within LIS. 
! 
!
!EOP
    
    integer                       :: ftn
    real                          :: gtmp(LVT_rc%lnc*LVT_rc%lnr)
    integer, intent(in)           :: gribid
    integer, intent(in)           :: gribSF
    integer, intent(in)           :: gribSfc
    integer, intent(in)           :: gribLvl
    integer, intent(in)           :: gribDis
    integer, intent(in)           :: gribCat
    integer, intent(in)           :: pdTemplate
    character(len=*), intent(in)  :: sType
    integer, intent(in)           :: timeRange
    integer, intent(in)           :: time_unit
    integer, intent(in)           :: time_p1
    integer, intent(in)           :: time_p2
    integer, intent(in)           :: k
    real                          :: toplev(1)
    real                          :: botlev(1)
    real                          :: depscale(1)
 
    integer                       :: sType_int
    character*8                   :: date
    integer                       :: yr1, mo1,da1,hr1,mn1
    integer                       :: idate,idate1
    real                          :: lat_ur, lon_ur
    real                          :: lat_ll, lon_ll
    integer                       :: igrib,iret
    integer                       :: decimalPrecision,gribSFtemp
    character*100                 :: message(20)

    ! Note passing string of defined points only to output
    ! because bitmap in GRIB-1 file will fill in the rest
    
    call grib_new_from_template(igrib,"GRIB2",iret)
    call LVT_verify(iret, 'grib_new_from_template failed in LVT_DataStreamsMod')

! Section 0: Indicator     
    
    call grib_set(igrib, 'discipline', gribDis, iret)
    call LVT_verify(iret, 'grib_set: discipline failed in LVT_DataStreamsMod')

! Section 1: Identification
    call grib_set(igrib,'centre',LVT_rc%grib_center_id,iret)
    call LVT_verify(iret,'grib_set:centre failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'subCentre',LVT_rc%grib_subcenter_id,iret)
    call LVT_verify(iret,'grib_set:subCentre failed in LVT_DataStreamsMod')

    call grib_set(igrib,'tablesVersion',LVT_rc%grib_table,iret)
    call LVT_verify(iret,'grib_set:tablesversion failed in LVT_DataStreamsMod')
    
 ! Hard-coded: significance of reference time 0=Analysis, 1.2.table
    call grib_set(igrib,'significanceOfReferenceTime',0,iret)
    call LVT_verify(iret,'grib_set:significanceOfReferenceTime failed in LVT_DataStreamsMod')

    da1=LVT_rc%da
    mo1=LVT_rc%mo
    yr1=LVT_rc%yr
    
    write(unit=date,fmt='(i4.4,i2.2,i2.2)') yr1,mo1,da1
    read(date,'(I8)') idate
    
    call grib_set(igrib,'dataDate',idate,iret)
    call LVT_verify(iret, 'grib_set:dataDate failed in LVT_DataStreamsMod')
    
    hr1=LVT_rc%hr
    mn1=LVT_rc%mn
    
    write(unit=date,fmt='(i2.2,i2.2)') hr1,mn1
    read(date,'(I4)') idate1
    
    call grib_set(igrib,'dataTime',idate1,iret)
    call LVT_verify(iret, 'grib_set:dataTime failed in LVT_DataStreamsMod')

      ! Hard-coded: production status of data 2=Research products, 1.3.table
    call grib_set(igrib,'productionStatusOfProcessedData',2,iret)
    call LVT_verify(iret,'grib_set:productionStatusOfProcessedData failed in LVT_DataStreamsMod')
    ! Hard-coded: type of data 0=Analysis, 1.4.table
    call grib_set(igrib,'typeOfProcessedData',0,iret)
    call LVT_verify(iret,'grib_set:typeOfProcessedData failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'stepType',sType, iret)
    call LVT_verify(iret,'grib_set:stepType failed in LVT_DataStreamsMod')


! Section 2: Local Use Section (Optional) --none for now

! Section 3: Grid

    call grib_set(igrib,'gridDefinitionTemplateNumber',LVT_rc%grib_grid_id,iret)
    call LVT_verify(iret,'grib_set:gridDefinitionTemplateNumber failed in LVT_DataStreamsMod')

    ! Hard-coded: shape of the Earth 0=radius = 6,367,470.0 m; 3.2.table 
    call grib_set(igrib,'shapeOfTheEarth',0,iret)
    call LVT_verify(iret,'grib_set:shapeOfTheEarth failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'swapScanningLat',1, iret)
    call LVT_verify(iret,'grib_set:swapScanningLat failed in LVT_DataStreamsMod')

    call grib_set(igrib,'Ni',LVT_rc%gnc,iret)
    call LVT_verify(iret, 'grib_set:Ni failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'Nj',LVT_rc%gnr,iret)
    call LVT_verify(iret, 'grib_set:Ni failed in LVT_DataStreamsMod')

    call ij_to_latlon(LVT_domain%lvtproj,float(LVT_rc%gnc),&
         float(LVT_rc%gnr),lat_ur,lon_ur)      
    call ij_to_latlon(LVT_domain%lvtproj,1.0, 1.0, &
         lat_ll,lon_ll)      
    
    call grib_set(igrib, 'latitudeOfFirstGridPointInDegrees',lat_ll,iret)
    call LVT_verify(iret, 'grib_set:latitudeOfFirstGridPointInDegrees failed in LVT_DataStreamsMod')
    
    call grib_set(igrib, 'longitudeOfFirstGridPointInDegrees',lon_ll,iret)
    call LVT_verify(iret, 'grib_set:longitudeOfFirstGridPointInDegrees failed in LVT_DataStreamsMod')
    
    call grib_set(igrib, 'latitudeOfLastGridPointInDegrees',lat_ur,iret)
    call LVT_verify(iret, 'grib_set:latitudeOfLastGridPointInDegrees failed in LVT_DataStreamsMod')
    
    call grib_set(igrib, 'longitudeOfLastGridPointInDegrees',lon_ur,iret)
    call LVT_verify(iret, 'grib_set:longitudeOfLastGridPointInDegrees failed in LVT_DataStreamsMod')


    if (LVT_rc%domain.eq."latlon") then 
       call grib_set(igrib,'gridType','regular_ll',iret)
       call LVT_verify(iret,'grib_set: gridType failed in LVT_DataStreamsMod')
       
       call grib_set(igrib,'iDirectionIncrementInDegrees',LVT_rc%gridDesc(9),iret)
       call LVT_verify(iret,'grib_set:iDirectionIncrementInDegrees failed in LVT_DataStreamsMod')
       
       call grib_set(igrib,'jDirectionIncrementInDegrees',LVT_rc%gridDesc(10),iret)
       call LVT_verify(iret,'grib_set:jDirectionIncrementInDegrees failed in LVT_DataStreamsMod')
       
    else  !Unsupported Map Projection for GRIB output
       
       message(1)='program:  LVT_DataStreamsMod'
       message(2)=' subroutine:  writevar_grib1_withstats_real'
       message(3)='  Unsupported map projection for GRIB1 output!'
       call lvt_abort(message)
       stop
       
    endif

! Section 4: Product Definition Section

    call grib_set(igrib,'productDefinitionTemplateNumber',pdTemplate, iret)
    call LVT_verify(iret,'grib_set:productDefinitionTemplateNumber failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'parameterCategory',gribCat, iret)
    call LVT_verify(iret,'grib_set:parameterCategory failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'parameterNumber',gribid, iret)
    call LVT_verify(iret,'grib_set:parameterNumber failed in LVT_DataStreamsMod')

!Hard-coded: Type of generating process 0=Analysis, 4.3.table
     call grib_set(igrib,'typeOfGeneratingProcess',0, iret)
     call LVT_verify(iret,'grib_set:typeOfGeneratingProcess failed in LVT_DataStreamsMod')

    call grib_set(igrib,'generatingProcessIdentifier',LVT_rc%grib_process_id,iret)
    call LVT_verify(iret,'grib_set:generatingProcessIdentifier failed in LVT_DataStreamsMod')
    
    call grib_set(igrib,'indicatorOfUnitOfTimeRange',time_unit, iret)
    call LVT_verify(iret,'grib_set:indicatorOfUnitOfTimeRange failed in LVT_DataStreamsMod')    
    
    if(sType.eq."avg") then 
       sType_int = 0
    elseif(sType.eq."accum") then 
       sType_int = 1
    endif

    if(stype.eq."instant") then 
       ! no need to specify statistical information template.4.statistcal.def
    else
       call grib_set(igrib,'typeOfStatisticalProcessing',sType_int, iret)
       call LVT_verify(iret,'grib_set:typeOfStatisticalProcessing failed in LVT_DataStreamsMod')
      ! Hard-coded: Type of time increment between successive fields used in the
      ! statistical processing (4.11.table)
      ! 6=local use for the backward (past) avg/accum/max/min in LVT
      call grib_set(igrib,'typeOfTimeIncrement',6, iret)
      call LVT_verify(iret,'grib_set:typeOfTimeIncrement failed in LVT_DataStreamsMod')
      
      call grib_set(igrib,'indicatorOfUnitForTimeRange',time_unit, iret)
      call LVT_verify(iret,'grib_set:indicatorOfUnitForTimeRange failed in LVT_DataStreamsMod')
      call grib_set(igrib,'lengthOfTimeRange',timeRange, iret)
      call LVT_verify(iret,'grib_set:lengthOfTimeRange failed in LVT_DataStreamsMod')
      call grib_set(igrib,'indicatorOfUnitForTimeIncrement',time_unit, iret)
      call LVT_verify(iret,'grib_set:indicatorOfUnitForTimeIncrement failed in LVT_DataStreamsMod')
      ! Hard-coded: Time increment between successive fields, in units defined
      ! by the previous Octets.  0 means that the statistical processing is the 
      ! result of a continuous processes, not the processing of a number of
      ! descrete samples.
      call grib_set(igrib,'timeIncrement',0, iret)
      call LVT_verify(iret,'grib_set:timeIncrement failed in LVT_DataStreamsMod')
   endif
   ! First and second surface types need to be set prior to scale and depth
   ! topLevel and bottomLevel or level will be derived by the settings below.
   ! Levels other than Surface and Soil layers are not considered yet--hkb
   call grib_set(igrib,'typeOfFirstFixedSurface',gribSfc, iret)
   call LVT_verify(iret,'grib_set:typeOfFirstFixedSurface failed in LVT_DataStreamsMod')     

   if ( gribSfc .eq. 106 ) then   ! soil layers
      call grib_set(igrib,'typeOfSecondFixedSurface',gribSfc, iret)
      call LVT_verify(iret,'grib_set:typeOfFirstFixedSurface failed in LVT_DataStreamsMod')
      ! Removed Hard-coded scale: depscale is passed from writeSingleGrib2Var
      call grib_set(igrib,'scaleFactorOfFirstFixedSurface',depscale(1), iret)
      call LVT_verify(iret,'grib_set:scaledFactorOfFirstFixedSurface failed in LVT_DataStreamsMod')   

      call grib_set(igrib,'scaledValueOfFirstFixedSurface',toplev(1), iret)
      call LVT_verify(iret,'grib_set:scaledValueOfFirstFixedSurface failed in LVT_DataStreamsMod')
      call grib_set(igrib,'scaleFactorOfSecondFixedSurface',depscale(1), iret)
      call LVT_verify(iret,'grib_set:scaledFactorOfFirstFixedSurface failed in LVT_DataStreamsMod')
      call grib_set(igrib,'scaledValueOfSecondFixedSurface',botlev(1), iret)
      call LVT_verify(iret,'grib_set:scaledValueOfSecondFixedSurface failed in LVT_DataStreamsMod')
   elseif ( gribSfc .eq. 1 ) then    ! surface
      call grib_set(igrib,'typeOfSecondFixedSurface',255, iret)
      call LVT_verify(iret,'grib_set:typeOfFirstFixedSurface failed in LVT_DataStreamsMod')
      call grib_set(igrib,'scaleFactorOfFirstFixedSurface',0, iret)
      call LVT_verify(iret,'grib_set:scaledFactorOfFirstFixedSurface failed in LVT_DataStreamsMod')
      call grib_set(igrib,'scaledValueOfFirstFixedSurface',toplev(1), iret)
      call LVT_verify(iret,'grib_set:scaledValueOfFirstFixedSurface failed in LVT_DataStreamsMod')

      call grib_set(igrib,'scaleFactorOfSecondFixedSurface',255, iret)
      call LVT_verify(iret,'grib_set:scaledFactorOfFirstFixedSurface failed in LVT_DataStreamsMod')
      call grib_set(igrib,'scaledValueOfSecondFixedSurface',255, iret)
      call LVT_verify(iret,'grib_set:scaledValueOfSecondFixedSurface failed in LVT_DataStreamsMod')
     else   ! 114 (snow level) or old 112 ??
       write(LVT_logunit,*) 'Warning: special surface type !! '//&
                     'verify scale/depth for ',gribSfc
      call grib_set(igrib,'typeOfSecondFixedSurface',gribSfc, iret)
      call LVT_verify(iret,'grib_set:typeOfFirstFixedSurface failed in LVT_DataStreamsMod')
      call grib_set(igrib,'scaleFactorOfFirstFixedSurface',0, iret)
      call LVT_verify(iret,'grib_set:scaledFactorOfFirstFixedSurface failed in LVT_DataStreamsMod')
      call grib_set(igrib,'scaledValueOfFirstFixedSurface',toplev(1), iret)
      call LVT_verify(iret,'grib_set:scaledValueOfFirstFixedSurface failed in LVT_DataStreamsMod')

      call grib_set(igrib,'scaleFactorOfSecondFixedSurface',0, iret)
      call LVT_verify(iret,'grib_set:scaledFactorOfFirstFixedSurface failed in LVT_DataStreamsMod')
      call grib_set(igrib,'scaledValueOfSecondFixedSurface',botlev(1), iret)
      call LVT_verify(iret,'grib_set:scaledValueOfSecondFixedSurface failed in LVT_DataStreamsMod')
   endif

! Section 5: Data Representation

     call grib_set(igrib,'packingType',LVT_rc%grib_packing_type,iret)
     call LVT_verify(iret, 'grib_set:packingType failed in LVT_DataStreamsMod')
     call grib_set(igrib, 'missingValue',LVT_rc%udef,iret)
     call LVT_verify(iret, 'grib_set:missingValue failed in LVT_DataStreamsMod')
     
! Should not need to fix the "num bits" value for each parameter
! if the "decimalPrecision" (aka, "DecScale") is set properly. - dmm
!     call grib_set(igrib, 'bitsPerValue',12,iret)
!     call LVT_verify(iret, 'grib_set:bitsPerValue failed in LVT_DataStreamsMod')

! Set the "decimalPrecision" (aka, "DecScale") based on the
! gribSF (grib scale factor) set in the MODEL OUTPUT TBL. - dmm
     gribSFtemp = gribSF
     decimalPrecision = 0
     do while (gribSFtemp.ge.10)
        decimalPrecision = decimalPrecision + 1
        gribSFtemp = gribSFtemp / 10
     enddo
     call grib_set(igrib, 'decimalPrecision',decimalPrecision,iret)
     call LVT_verify(iret, 'grib_set:decimalPrecision failed in LVT_DataStreamsMod')

! Section 6: Bit-Map

    call grib_set(igrib, 'bitmapPresent',1,iret)
    call LVT_verify(iret, 'grib_set:bitmapPresent failed in LVT_DataStreamsMod')
 
    call grib_set(igrib,'values',gtmp,iret)
    call LVT_verify(iret, 'grib_set:values failed in LVT_DataStreamsMod')
    
    call grib_write(igrib,ftn,iret)
    call LVT_verify(iret, 'grib_write failed in LVT_DataStreamsMod')
    
    call grib_release(igrib,iret)
    call LVT_verify(iret,'grib_release failed in LVT_DataStreamsMod')


  end subroutine writeSingleGrib2Var

!BOP
! !ROUTINE: defineNETCDFheaderVar
! \label{defineNETCDFheaderVar}
! 
! !INTERFACE: 
  subroutine defineNETCDFheaderVar(ftn, dimID, dataEntry)  

! !USES: 

! !ARGUMENTS:     
    integer                                 :: ftn
    type(LVT_lismetadataEntry), pointer     :: dataEntry
    integer                                 :: dimID(4)
! 
! !DESCRIPTION: 
!    This routine writes the required NETCDF header for a single variable
! 
!   The arguments are: 
!   \begin{description}
!   \item[n]
!    index of the nest
!   \item[ftn]
!    NETCDF file unit handle
!   \item[dimID]
!    NETCDF dimension ID corresponding to the variable
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be 
    !    written
!   \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LIS\_endrun](\ref{LIS_endrun})
!     call to abort program when a fatal error is detected. 
!   \item[LIS\_verify](\ref{LVT_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!EOP
    logical       :: nmodel_status
    integer       :: data_index
    integer       :: shuffle, deflate, deflate_level
    integer       :: fill_value

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    data_index = dataEntry%index
    fill_value = LVT_rc%udef

    shuffle = NETCDF_shuffle
    deflate = NETCDF_deflate
    deflate_level =NETCDF_deflate_level

    if(dataEntry%selectOpt.eq.1)then 
       if(dataEntry%vlevels.gt.1) then 
          call LVT_verify(nf90_def_dim(ftn,&
               trim(dataEntry%short_name)//'_profiles',&
               dataEntry%vlevels, dimID(3)),&
               'nf90_def_dim failed (2d gridspace) in LVT_DataStreamsMod')
       endif
       if(dataEntry%vlevels.gt.1) then 
          call LVT_verify(nf90_def_var(ftn,trim(dataEntry%short_name),&
               nf90_float,&
               dimids = dimID(1:3), varID=dataEntry%varId_def),&
               'nf90_def_var for '//trim(dataEntry%short_name)//&
               'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
          call LVT_verify(nf90_def_var_fill(ftn,&
               dataEntry%varId_def, &
               1,fill_value), 'nf90_def_var_fill failed for '//&
               dataEntry%short_name)
          
          call LVT_verify(nf90_def_var_deflate(ftn,&
               dataEntry%varId_def,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
               'failed in defineNETCDFheadervar')                     
#endif
       else
          call LVT_verify(nf90_def_var(ftn,trim(dataEntry%short_name),&
               nf90_float,&
               dimids = dimID(1:2), varID=dataEntry%varId_def),&
               'nf90_def_var for '//trim(dataEntry%short_name)//&
               'failed in defineNETCDFheadervar')                     
#if(defined USE_NETCDF4)
          call LVT_verify(nf90_def_var_fill(ftn,&
               dataEntry%varId_def, &
               1,fill_value), 'nf90_def_var_fill failed for '//&
               dataEntry%short_name)                
          
          call LVT_verify(nf90_def_var_deflate(ftn,&
               dataEntry%varId_def,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
               'failed in defineNETCDFheadervar')                     
#endif                
       endif
       
       call LVT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "units",trim(dataEntry%units)),&
            'nf90_put_att for units failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "standard_name",trim(dataEntry%standard_name)),&
            'nf90_put_att for standard_name failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "long_name",trim(dataEntry%long_name)),&
            'nf90_put_att for long_name failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "scale_factor",1.0),&
            'nf90_put_att for scale_factor failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "add_offset",0.0),&
            'nf90_put_att for add_offset failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "missing_value",LVT_rc%udef),&
            'nf90_put_att for missing_value failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "_FillValue",LVT_rc%udef),&
            'nf90_put_att for _FillValue failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "vmin",dataEntry%valid_min),&
            'nf90_put_att for vmin failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varId_def,&
            "vmax",dataEntry%valid_max),&
            'nf90_put_att for vmax failed in defineNETCDFheaderVar')

    endif
#endif
  end subroutine defineNETCDFheaderVar


!BOP
! !ROUTINE: defineNETCDFheaderVar_SS
! \label{defineNETCDFheaderVar_SS}
! 
! !INTERFACE: 
  subroutine defineNETCDFheaderVar_SS(ftn, dimID, dataEntry)  

! !USES: 

! !ARGUMENTS:     
    integer                                 :: ftn
    type(LVT_lismetadataEntry), pointer     :: dataEntry
    integer                                 :: dimID(4)
! 
! !DESCRIPTION: 
!    This routine writes the required NETCDF header for a single variable
! 
!   The arguments are: 
!   \begin{description}
!   \item[n]
!    index of the nest
!   \item[ftn]
!    NETCDF file unit handle
!   \item[dimID]
!    NETCDF dimension ID corresponding to the variable
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be 
    !    written
!   \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LIS\_endrun](\ref{LIS_endrun})
!     call to abort program when a fatal error is detected. 
!   \item[LIS\_verify](\ref{LVT_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!EOP
    logical       :: nmodel_status
    integer       :: data_index
    integer       :: shuffle, deflate, deflate_level
    integer       :: fill_value

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    data_index = dataEntry%index
    fill_value = LVT_rc%udef

    shuffle = NETCDF_shuffle
    deflate = NETCDF_deflate
    deflate_level =NETCDF_deflate_level

    if(dataEntry%selectOpt.eq.1)then 
       if(dataEntry%vlevels.gt.1) then 
          call LVT_verify(nf90_def_dim(ftn,&
               trim(dataEntry%short_name)//'_profiles',&
               dataEntry%vlevels, dimID(3)),&
               'nf90_def_dim failed (2d gridspace) in LVT_DataStreamsMod')
       endif
       if(dataEntry%vlevels.gt.1) then 
          call LVT_verify(nf90_def_var(ftn,trim(dataEntry%short_name),&
               nf90_float,&
               dimids = dimID(1:3), varID=dataEntry%varid_ss),&
               'nf90_def_var for '//trim(dataEntry%short_name)//&
               'failed in defineNETCDFheadervar')
#if(defined USE_NETCDF4)
          call LVT_verify(nf90_def_var_fill(ftn,&
               dataEntry%varid_ss, &
               1,fill_value), 'nf90_def_var_fill failed for '//&
               dataEntry%short_name)
          
          call LVT_verify(nf90_def_var_deflate(ftn,&
               dataEntry%varid_ss,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
               'failed in defineNETCDFheadervar')                     
#endif
       else
          call LVT_verify(nf90_def_var(ftn,trim(dataEntry%short_name),&
               nf90_float,&
               dimids = dimID(1:2), varID=dataEntry%varid_ss),&
               'nf90_def_var for '//trim(dataEntry%short_name)//&
               'failed in defineNETCDFheadervar')                     
#if(defined USE_NETCDF4)
          call LVT_verify(nf90_def_var_fill(ftn,&
               dataEntry%varid_ss, &
               1,fill_value), 'nf90_def_var_fill failed for '//&
               dataEntry%short_name)                
          
          call LVT_verify(nf90_def_var_deflate(ftn,&
               dataEntry%varid_ss,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for '//trim(dataEntry%short_name)//&
               'failed in defineNETCDFheadervar')                     
#endif                
       endif
       
       call LVT_verify(nf90_put_att(ftn,dataEntry%varid_ss,&
            "units",trim(dataEntry%units)),&
            'nf90_put_att for units failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varid_ss,&
            "standard_name",trim(dataEntry%standard_name)),&
            'nf90_put_att for standard_name failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varid_ss,&
            "long_name",trim(dataEntry%long_name)),&
            'nf90_put_att for long_name failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varid_ss,&
            "scale_factor",1.0),&
            'nf90_put_att for scale_factor failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varid_ss,&
            "add_offset",0.0),&
            'nf90_put_att for add_offset failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varid_ss,&
            "missing_value",LVT_rc%udef),&
            'nf90_put_att for missing_value failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varid_ss,&
            "_FillValue",LVT_rc%udef),&
            'nf90_put_att for _FillValue failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varid_ss,&
            "vmin",dataEntry%valid_min),&
            'nf90_put_att for vmin failed in defineNETCDFheaderVar')
       call LVT_verify(nf90_put_att(ftn,dataEntry%varid_ss,&
            "vmax",dataEntry%valid_max),&
            'nf90_put_att for vmax failed in defineNETCDFheaderVar')

    endif
#endif
  end subroutine defineNETCDFheaderVar_SS

!BOP
! 
! !ROUTINE: writeSingleNetcdfVar
! \label{writeSingleNetcdfVar}
!
! !INTERFACE: 
  subroutine writeSingleNetcdfVar(ftn,gtmp,varID,k)
!
! !DESCRIPTION: 
!  This subroutine writes a single variable to a grib file
! 
!EOP

! !ARGUMENTS:     
    integer                       :: ftn
    real                          :: gtmp(LVT_rc%lnc*LVT_rc%lnr)
    integer, intent(in)           :: varid
    integer, intent(in)           :: k

    integer                       :: c,r
    real                          :: gtmp2d(LVT_rc%lnc,LVT_rc%lnr)
#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    do r=1,LVT_rc%lnr
       do c=1,LVT_rc%lnc
          gtmp2d(c,r) = gtmp(c+(r-1)*LVT_rc%lnc)
       enddo
    enddo

    call LVT_verify(nf90_put_var(ftn,varID, gtmp2d,(/1,1,k/),&
         (/LVT_rc%gnc,LVT_rc%gnr,1/)),&
         'nf90_put_var failed for in LVT_DataStreamsMod')
    
#endif

  end subroutine writeSingleNetcdfVar

!BOP
! 
! !ROUTINE: LVT_tavgDataStreams
! \label{LVT_tavgDataStreams}
!
! !INTERFACE: 
  subroutine LVT_tavgDataStreams
! 
! !USES:   
    use LVT_statsDataMod

    implicit none
!
!
! !DESCRIPTION: 
!   This routine invokes the calls to compute temporal averages of 
!   desired set of variables, based on the specified 
!   temporal averaging frequency. 
!  
!   The routines invoked are: 
!   \begin{description}
!    \item[tavgSingleDataStream](\ref{tavgSingleDataStream})
!     computes the temporal average for a single variable
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer      :: kk
    type(LVT_metadataEntry), pointer :: dataEntry
    type(LVT_metadataEntry), pointer :: ds1, ds2

    if(LVT_rc%computeFlag) then            

!data stream 1
       do kk=1,LVT_rc%nDataStreams
          if(kk.eq.1) then 
             dataEntry => LVT_histData%head_ds1_list
          elseif(kk.eq.2) then 
             dataEntry => LVT_histData%head_ds2_list
          elseif(kk.eq.3) then 
             dataEntry => LVT_histData%head_ds3_list
          endif
       
          do while(associated(dataEntry))
             call tavgSingleDataStream(dataEntry)
             dataEntry => dataEntry%next
          enddo

! copy duplicate entries
! Note that this check is not enabled for three datastrems. 
! The responsibility of ensuring non-duplicate entries is
! on the user. 
          if(LVT_rc%ds1_dup) then 
             ds1 => LVT_histData%head_ds1_list       
             do while(associated(ds1))
                ds2 => ds1%next
                do while(associated(ds2))
                   if(ds2%index.ne.ds1%index.and.&
                        ds1%short_name.eq.ds2%short_name) then 
                      ds2%value = ds1%value
                      ds2%count = ds1%count
                   endif
                   ds2 => ds2%next
                enddo
                ds1 => ds1%next
             enddo
             
          endif

          if(LVT_rc%ds2_dup) then 
             ds1 => LVT_histData%head_ds2_list       
             do while(associated(ds1))
                ds2 => ds1%next
                do while(associated(ds2))
                   if(ds2%index.ne.ds1%index.and.&
                        ds1%short_name.eq.ds2%short_name) then 
                      ds2%value = ds1%value
                      ds2%count = ds1%count
                   endif
                   ds2 => ds2%next
                enddo
                ds1 => ds1%next
             enddo
             
          endif

          if(LVT_rc%var_based_strat .gt. 0) then
             call tavgSingleDataStream(LVT_histData%strat_varEntry)
             LVT_stats%strat_var(:,:,:) = &
                  LVT_histData%strat_varEntry%value(:,:,:)
          endif
       enddo
    endif
  end subroutine LVT_tavgDataStreams

!BOP
! 
! !ROUTINE: tavgSingleDataStream
! \label{tavgSingleDataStream}
!
! !INTERFACE:
  subroutine tavgSingleDataStream( dataEntry)
! 
! !USES: 

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine temporally averages the accumulated data in a
!   given datastream
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    type(LVT_metadataEntry) :: dataEntry
!EOP
    integer :: k,t,c,r,m,gid

    ! EMK Do not average the data if the raw input are accumulations.
    ! This is indicated in the lvt.config file
    if (dataEntry%timeAvgOpt .eq. 3) return

    if(dataEntry%selectNlevs.ge.1) then 
       if(LVT_rc%computeEnsMetrics.eq.1) then 
          do t=1,LVT_LIS_rc(1)%ntiles
             do k=1,dataEntry%vlevels
                c = LVT_LIS_domain(1)%tile(t)%col
                r = LVT_LIS_domain(1)%tile(t)%row
                if(LVT_LIS_domain(1)%gindex(c,r).ne.-1) then 
                   gid = LVT_LIS_domain(1)%gindex(c,r)
                   do m=1,LVT_rc%nensem
                      if(dataEntry%count(t,m,k).ne.0) then 
                         dataEntry%value(t,m,k) = &
                              dataEntry%value(t,m,k)/dataEntry%count(t,m,k)
                         
                      endif
                   enddo
                endif
             enddo
          enddo         
       else
          do r=1,LVT_rc%lnr
             do c=1,LVT_rc%lnc
                do k=1,dataEntry%vlevels
                   if(LVT_domain%gindex(c,r).ne.-1) then 
                      gid = LVT_domain%gindex(c,r)
                      do m=1,LVT_rc%nensem
                         if(dataEntry%count(gid,m,k).ne.0) then 
                            dataEntry%value(gid,m,k) = &
                                 dataEntry%value(gid,m,k)/&
                                 dataEntry%count(gid,m,k)
                            
                         endif
                      enddo
                   endif
                enddo
             enddo
          enddo
       endif
    endif

  end subroutine tavgSingleDataStream


!BOP
! 
! !ROUTINE: LVT_resetDataStreams
! \label{LVT_resetDataStreams}
!
! !INTERFACE: 
  subroutine LVT_resetDataStreams
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine reinitializes the data structures that hold the observational
!   data
! 
!   The routines invoked are: 
!   \begin{description}
!    \item[resetSingleDataStream2](\ref{resetSingleDataStream2})
!     resets the datastructures for a single variable
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    type(LVT_metadataEntry), pointer :: ds1
    type(LVT_metadataEntry), pointer :: ds2
    type(LVT_metadataEntry), pointer :: ds3
    
    if(LVT_rc%computeFlag) then            

!data stream 1
       ds1 => LVT_histData%head_ds1_list
       
       do while(associated(ds1))
          call resetSingleDataStream(ds1)
          ds1 => ds1%next
       enddo
       
!data stream 2
       ds2 => LVT_histData%head_ds2_list
       
       do while(associated(ds2))
          call resetSingleDataStream(ds2)
          ds2 => ds2%next
       enddo
       
       if(LVT_rc%nDataStreams.gt.2) then 
          
!data stream 3
          ds3 => LVT_histData%head_ds3_list
          
          do while(associated(ds3))
             call resetSingleDataStream(ds3)
             ds3 => ds3%next
          enddo
       endif
!need special handler for LIS output
       if(LVT_rc%lis_output_obs) then 
          if(LVT_rc%obssource(1).eq."LIS output") then 
             call LVT_resetLISoutputContainers(1)
          endif
          if(LVT_rc%obssource(2).eq."LIS output") then 
             call LVT_resetLISoutputContainers(2)
          endif
          if(LVT_rc%nDataStreams.gt.2) then 
             if(LVT_rc%obssource(3).eq."LIS output") then 
                call LVT_resetLISoutputContainers(3)
             endif
          endif
       endif
    endif
  end subroutine LVT_resetDataStreams


!BOP
! 
! !ROUTINE: resetSingleDataStream
! \label{resetSingleDataStream}
!
! !INTERFACE: 
  subroutine resetSingleDataStream(dataEntry)
! 
! !USES:   
    implicit none 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine resets the data structures that hold the observational 
!  data and the temporal averaging counters
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 
    type(LVT_metadataEntry) :: dataEntry
! 
!EOP

    integer                 :: k 

    if(dataEntry%selectNlevs.ge.1) then 
       do k=1,dataEntry%vlevels
          dataEntry%value(:,:,k) = 0 
          dataEntry%count(:,:,k) = 0 
          if(dataEntry%stdev_flag) then 
             dataEntry%count_stdev(:,k)= 0 
             dataEntry%stdev(:,k) = 0
          endif
       enddo      
    endif
  end subroutine resetSingleDataStream
end module LVT_DataStreamsMod
