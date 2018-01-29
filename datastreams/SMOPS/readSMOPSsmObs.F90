!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readSMOPSsmObs
! \label{readSMOPSsmObs}
!
! !INTERFACE: 
subroutine readSMOPSsmObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use SMOPSsm_obsMod, only : SMOPSsmobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer,   intent(in) :: source
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!
! !REVISION HISTORY: 
!  11 Dec 2014: Sujay Kumar, Initial Specification
! 
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j, status
  integer           :: yr, mo, da, hr, mn, ss
  character*100     :: fname
  real              :: smobs(LVT_rc%lnc*LVT_rc%lnr)
  real              :: smobs_2d(LVT_rc%lnc,LVT_rc%lnr)
  real              :: relsmc(LVT_rc%lnc,LVT_rc%lnr)
  real              :: smobs_av(LVT_rc%lnc, LVT_rc%lnr)
  integer           :: count_av(LVT_rc%lnc, LVT_rc%lnr)
  type(ESMF_Time)   :: cTime, stTime, enTime

!-----------------------------------------------------------------------
! It is assumed that CDF is computed using daily observations. 
!-----------------------------------------------------------------------

  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  smobs= LVT_rc%udef
  smobs_2d = LVT_rc%udef

  if(SMOPSsmobs(source)%startmode.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 

     LVT_rc%resetFlag(source) = .false.      
     SMOPSsmobs(source)%startmode = .false. 

     if(LVT_rc%smoothObs.eq.1) then 
        smobs_av = 0.0
        count_av = 0

        call ESMF_TimeSet(cTime, yy = LVT_rc%dyr(source), &
             mm = LVT_rc%dmo(source), &
             dd = LVT_rc%dda(source), &
             h  = LVT_rc%dhr(source), &
             m  = LVT_rc%dmn(source), & 
             s  = LVT_rc%dss(source), &
             calendar = LVT_calendar, & 
             rc = status)
        call LVT_verify(status)
        
        stTime = cTime - LVT_obsSmTwL
        enTime = cTime + LVT_obsSmTwL

        cTime = stTime
        do while (cTime.le.enTime) 
           call ESMF_TimeGet(cTime, yy = yr, mm=mo, dd=da,&
                h = hr, m = mn, s =ss, calendar=LVT_calendar, &
                rc=status)
           call LVT_verify(status)
           
           call create_SMOPSsm_filename(SMOPSsmobs(source)%odir, &
                yr, mo, da, fname)
           
           inquire(file=trim(fname),exist=file_exists)
           if(file_exists) then
              smobs = LVT_rc%udef
              write(LVT_logunit,*) '[INFO] Reading ',trim(fname)
              call read_SMOPS_data(source, fname,smobs)
           endif
           
           cTime = cTime + LVT_obsSmTwI
           
           do r=1,LVT_rc%lnr
              do c=1,LVT_rc%lnc
                 if(smobs(c+(r-1)*LVT_rc%lnc).ne.LVT_rc%udef) then 
                    smobs_av(c,r) = smobs_av(c,r) + smobs(c+(r-1)*LVT_rc%lnc)
                    count_av(c,r) = count_av(c,r) + 1
                 endif
              enddo
           enddo
        enddo
        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(count_av(c,r).gt.0) then 
                 smobs_2d(c,r) = smobs_av(c,r)/count_av(c,r)
              else
                 smobs_2d(c,r) = LVT_rc%udef
              endif
           enddo
        enddo
        
     else
        call create_SMOPSsm_filename(SMOPSsmobs(source)%odir, &
             LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source), fname)
        
        inquire(file=trim(fname),exist=file_exists)
        if(file_exists) then
           
           write(LVT_logunit,*) '[INFO] Reading ',trim(fname)
           call read_SMOPS_data(source, fname,smobs)
        endif
        
        do r=1,LVT_rc%lnr
           do c=1,LVT_rc%lnc
              if(smobs(c+(r-1)*LVT_rc%lnc).ne.-9999.0) then 
                 smobs_2d(c,r) = smobs(c+(r-1)*LVT_rc%lnc)
              endif
           enddo
        enddo

     endif
  else
     smobs_2d = LVT_rc%udef     
  endif

  relsmc = LVT_rc%udef

  if(LVT_rc%curr_pass.eq.1) then 
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(smobs_2d(c,r).ne.LVT_rc%udef) then 
              if(smobs_2d(c,r).gt.SMOPSsmobs(source)%maxv(c,r)) then
                 SMOPSsmobs(source)%maxv(c,r) = smobs_2d(c,r)
              endif
              
              if(smobs_2d(c,r).lt.SMOPSsmobs(source)%minv(c,r)) then
                 SMOPSsmobs(source)%minv(c,r) = smobs_2d(c,r)
              endif
           endif
        enddo
     enddo
  elseif(LVT_rc%curr_pass.eq.2) then 

     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc

           if(smobs_2d(c,r).ne.LVT_rc%udef.and.&
                (SMOPSsmobs(source)%maxv(c,r).ne.SMOPSsmobs(source)%minv(c,r))) then 
              relsmc(c,r) = (smobs_2d(c,r) - SMOPSsmobs(source)%minv(c,r))/&
                   (SMOPSsmobs(source)%maxv(c,r) - SMOPSsmobs(source)%minv(c,r))
           else
              relsmc(c,r) = LVT_rc%udef
           endif
              
        enddo
     enddo
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_RELSMC, source, &
       relsmc,vlevel=1,units="%")

  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist, source, &
       smobs_2d,vlevel=1,units="m3/m3")

  open(100,file='test.bin',form='unformatted')
  write(100) smobs_2d
  close(100)
  stop

end subroutine readSMOPSsmObs


!BOP
! 
! !ROUTINE: read_SMOPS_data
! \label(read_SMOPS_data)
!
! !INTERFACE:
subroutine read_SMOPS_data(source, fname, smobs_ip)
! 
! !USES:   
  use grib_api
  use LVT_coreMod,  only : LVT_rc
  use LVT_logMod
  use map_utils,    only : latlon_to_ij
  use SMOPSsm_obsMod, only : SMOPSsmobs

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: source 
  character (len=*)             :: fname
  real                          :: smobs_ip(LVT_rc%lnc*LVT_rc%lnr)


! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the SMOPS grib2 file and applies the data
!  quality flags to filter the data. The retrievals are rejected when 
!  the estimated error is above a predefined threshold (the recommeded
!  value is 5%). 
! 
! The information from gribtab is: 
!{ 2, 0, 7, 1, 3, 210, "BLENDEDSM", "NOAA Blended Liquid Volumetric Soil Moisture (Non-Frozen)", "m3m-3"},
!{ 2, 0, 7, 1, 3, 211, "AMSRESM", "NOAA AMSR-E Liquid Volumetric Soil Moisture (Non-Frozen)", "m3m-3"},
!{ 2, 0, 7, 1, 3, 212, "SMOSSM", "SMOS Liquid Volumetric Soil Moisture (Non-Frozen)", "m3m-3"},
!{ 2, 0, 7, 1, 3, 213, "ASCATSM", "ASCAT Liquid Volumetric Soil Moisture (Non-Frozen)", "m3m-3"},
!{ 2, 0, 7, 1, 3, 214, "WINDSATSM", "NOAA WindSat Liquid Volumetric Soil Moisture (Non-Frozen)", "m3m-3"},
!{ 2, 0, 7, 1, 3, 215, "RESSM1", "Reserved Liquid Volumetric Soil Moisture 1 (Non-Frozen)", "m3m-3"},
!{ 2, 0, 7, 1, 3, 216, "RESSM2", "Reserved Liquid Volumetric Soil Moisture 2 (Non-Frozen)", "m3m-3"},
!{ 2, 0, 7, 1, 3, 217, "BLENDEDSMHR", "Observation Hour for BLENDEDSM", "-"},
!{ 2, 0, 7, 1, 3, 218, "BLENDEDSMMN", "Observation Minute for BLENDEDSM", "-"},
!{ 2, 0, 7, 1, 3, 219, "AMSRESMHR", "Observation Hour for AMSRESM", "-"},
!{ 2, 0, 7, 1, 3, 220, "AMSRESMMN", "Observation Minute for AMSRESM", "-"},
!{ 2, 0, 7, 1, 3, 221, "SMOSSMHR", "Observation Hour for SMOSSM", "-"},
!{ 2, 0, 7, 1, 3, 222, "SMOSSMMN", "Observation Minute for SMOSSM", "-"},
!{ 2, 0, 7, 1, 3, 223, "ASCATSMHR", "Observation Hour for ASCATSM", "-"},
!{ 2, 0, 7, 1, 3, 224, "ASCATSMMN", "Observation Minute for ASCATSM", "-"},
!{ 2, 0, 7, 1, 3, 225, "WINDSATSMHR", "Observation Hour for WINDSATSM", "-"},
!{ 2, 0, 7, 1, 3, 226, "WINDSATSMMN", "Observation Minute for WINDSATSM", "-"},
!{ 2, 0, 7, 1, 3, 227, "RESSM1HR", "Observation Hour for RESSM1", "-"},
!{ 2, 0, 7, 1, 3, 228, "RESSM1MN", "Observation Minute for RESSM1", "-"},
!{ 2, 0, 7, 1, 3, 229, "RESSM2HR", "Observation Hour for RESSM2", "-"},
!{ 2, 0, 7, 1, 3, 230, "RESSM2MN", "Observation Minute for RESSM2", "-"},
!{ 2, 0, 7, 1, 3, 231, "BLENDEDSMQA", "NOAA Blended Liquid Volumetric Soil Moisture QA", "-"},
!{ 2, 0, 7, 1, 3, 232, "AMSRESMQA", "NOAA AMSR-E Liquid Volumetric Soil Moisture QA", "-"},
!{ 2, 0, 7, 1, 3, 233, "SMOSSMQA", "SMOS Liquid Volumetric Soil Moisture QA", "-"},
!{ 2, 0, 7, 1, 3, 234, "ASCATSMQA", "ASCAT Liquid Volumetric Soil Moisture QA", "-"},
!{ 2, 0, 7, 1, 3, 235, "WINDSATSMQA", "NOAA WindSat Liquid Volumetric Soil Moisture QA", "-"},
!{ 2, 0, 7, 1, 3, 236, "RESSM1QA", "Reserved Liquid Volumetric Soil Moisture 1 QA", "-"},
!{ 2, 0, 7, 1, 3, 237, "RESSM2QA", "Reserved Liquid Volumetric Soil Moisture 2 QA", "-"},

!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the SMOPS AMSR-E file
!  \item[smobs\_ip]    soil moisture data processed to the LVT domain
! \end{description}
!
!
!EOP
  INTEGER*2, PARAMETER :: FF = 255
  real,    parameter  :: err_threshold = 5 ! in percent
!  real,    parameter  :: err_threshold = 100 ! in percent
! pre-3.0 versions
!  integer, parameter  :: param_ASCAT= 213, param_ASCAT_qa = 234
!  integer, parameter  :: param_windsat = 214, param_windsat_qa = 235
!  integer, parameter  :: param_SMOS = 212, param_SMOS_qa = 233

  integer, parameter  :: param_ASCAT= 213, param_ASCAT_qa = 234
  integer, parameter  :: param_windsat = 214, param_windsat_qa = 235
  integer, parameter  :: param_SMOS = 212, param_SMOS_qa = 242
!SMAPSM
  integer, parameter  :: param_SMAP = 218, param_SMAP_qa = 248
!NRTSMAP
!  integer, parameter  :: param_SMAP = 217, param_SMAP_qa = 247
  real                :: sm_ascat(SMOPSsmobs(source)%rtsmopsnc*&
       SMOPSsmobs(source)%rtsmopsnr)
  real                :: sm_ascat_t(SMOPSsmobs(source)%rtsmopsnc*&
  SMOPSsmobs(source)%rtsmopsnr)
  real                :: sm_ascat_qa(SMOPSsmobs(source)%rtsmopsnc*&
       SMOPSsmobs(source)%rtsmopsnr)
  integer*2           :: sm_ascat_qa_t(SMOPSsmobs(source)%rtsmopsnc*&
  SMOPSsmobs(source)%rtsmopsnr)
  real                :: sm_windsat(SMOPSsmobs(source)%rtsmopsnc*&
       SMOPSsmobs(source)%rtsmopsnr)
  real                :: sm_windsat_t(SMOPSsmobs(source)%rtsmopsnc*&
  SMOPSsmobs(source)%rtsmopsnr)
  real                :: sm_windsat_qa(SMOPSsmobs(source)%rtsmopsnc*&
       SMOPSsmobs(source)%rtsmopsnr)
  integer*2           :: sm_windsat_qa_t(SMOPSsmobs(source)%rtsmopsnc*&
  SMOPSsmobs(source)%rtsmopsnr)
  real                :: sm_smos(SMOPSsmobs(source)%rtsmopsnc*&
       SMOPSsmobs(source)%rtsmopsnr)
  real                :: sm_smos_t(SMOPSsmobs(source)%rtsmopsnc*&
  SMOPSsmobs(source)%rtsmopsnr)
  real                :: sm_smos_qa(SMOPSsmobs(source)%rtsmopsnc*&
       SMOPSsmobs(source)%rtsmopsnr)
  integer*2           :: sm_smos_qa_t(SMOPSsmobs(source)%rtsmopsnc*&
  SMOPSsmobs(source)%rtsmopsnr)
  real                :: sm_smap(SMOPSsmobs(source)%rtsmopsnc*&
       SMOPSsmobs(source)%rtsmopsnr)
  real                :: sm_smap_t(SMOPSsmobs(source)%rtsmopsnc*&
  SMOPSsmobs(source)%rtsmopsnr)
  real                :: sm_smap_qa(SMOPSsmobs(source)%rtsmopsnc*&
       SMOPSsmobs(source)%rtsmopsnr)
  integer*2           :: sm_smap_qa_t(SMOPSsmobs(source)%rtsmopsnc*&
  SMOPSsmobs(source)%rtsmopsnr)
  real                :: sm_data(SMOPSsmobs(source)%rtsmopsnc*&
       SMOPSsmobs(source)%rtsmopsnr)
  logical*1           :: sm_data_b(SMOPSsmobs(source)%rtsmopsnc*&
       SMOPSsmobs(source)%rtsmopsnr)
  logical*1           :: smobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)

  integer             :: c,r,i,j,kk
  integer             :: ftn,iret,igrib,nvars
  integer             :: param_num
  logical             :: var_found_ascat
  logical             :: var_found_windsat
  logical             :: var_found_smos
  logical             :: var_found_smap
  real                :: err, ql

  call grib_open_file(ftn,trim(fname), 'r',iret)
  if(iret.ne.0) then 
     write(LVT_logunit,*) '[ERR] Could not open file: ',trim(fname)
     call LVT_endrun()
  endif
  call grib_multi_support_on

  do
     call grib_new_from_file(ftn,igrib,iret)

     if ( iret == GRIB_END_OF_FILE ) then
        exit
     endif

     call grib_get(igrib, 'parameterNumber',param_num, iret)
     call LVT_verify(iret, &
          'grib_get: parameterNumber failed in readSMOPSsmobs')

     var_found_ascat = .false. 
     if(SMOPSsmobs(source)%useASCAT.eq.1) then 
        if(param_num.eq.param_ascat) then 
           var_found_ascat = .true.
        endif
     endif
     
     if(var_found_ascat) then
        call grib_get(igrib, 'values',sm_ascat,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%rtsmopsnr
           do c=1,SMOPSsmobs(source)%rtsmopsnc
              sm_ascat_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = &
                   sm_ascat(c+((SMOPSsmobs(source)%rtsmopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%rtsmopsnc)
           enddo
        enddo     
        
     endif

     var_found_ascat = .false. 
     if(SMOPSsmobs(source)%useASCAT.eq.1) then 
        if(param_num.eq.param_ascat_qa) then 
           var_found_ascat = .true.
        endif
     endif
     
     if(var_found_ascat) then
        call grib_get(igrib, 'values',sm_ascat_qa,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%rtsmopsnr
           do c=1,SMOPSsmobs(source)%rtsmopsnc
              sm_ascat_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = &
                   INT(sm_ascat_qa(c+((SMOPSsmobs(source)%rtsmopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%rtsmopsnc))
           enddo
        enddo       
     endif

!windsat
     var_found_windsat = .false. 
     if(SMOPSsmobs(source)%useWINDSAT.eq.1) then 
        if(param_num.eq.param_windsat) then 
           var_found_windsat = .true.
        endif
     endif
     
     if(var_found_windsat) then
        call grib_get(igrib, 'values',sm_windsat,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%rtsmopsnr
           do c=1,SMOPSsmobs(source)%rtsmopsnc
              sm_windsat_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = &
                   sm_windsat(c+((SMOPSsmobs(source)%rtsmopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%rtsmopsnc)
           enddo
        enddo     
        
     endif

     var_found_windsat = .false. 
     if(SMOPSsmobs(source)%useWINDSAT.eq.1) then 
        if(param_num.eq.param_windsat_qa) then 
           var_found_windsat = .true.
        endif
     endif
     
     if(var_found_windsat) then
        call grib_get(igrib, 'values',sm_windsat_qa,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%rtsmopsnr
           do c=1,SMOPSsmobs(source)%rtsmopsnc
              sm_windsat_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = &
                   INT(sm_windsat_qa(c+((SMOPSsmobs(source)%rtsmopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%rtsmopsnc))
           enddo
        enddo       
     endif
!smos
     var_found_smos = .false. 
     if(SMOPSsmobs(source)%useSMOS.eq.1) then 
        if(param_num.eq.param_smos) then 
           var_found_smos = .true.
        endif
     endif
     
     if(var_found_smos) then
        call grib_get(igrib, 'values',sm_smos,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%rtsmopsnr
           do c=1,SMOPSsmobs(source)%rtsmopsnc
              sm_smos_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = &
                   sm_smos(c+((SMOPSsmobs(source)%rtsmopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%rtsmopsnc)
           enddo
        enddo     
        
     endif

     var_found_smos = .false. 
     if(SMOPSsmobs(source)%useSMOS.eq.1) then 
        if(param_num.eq.param_smos_qa) then 
           var_found_smos = .true.
        endif
     endif
     
     if(var_found_smos) then
        call grib_get(igrib, 'values',sm_smos_qa,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%rtsmopsnr
           do c=1,SMOPSsmobs(source)%rtsmopsnc
              sm_smos_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = &
                   INT(sm_smos_qa(c+((SMOPSsmobs(source)%rtsmopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%rtsmopsnc))
           enddo
        enddo       
     endif

!smap
     var_found_smap = .false. 
     if(SMOPSsmobs(source)%useSMAP.eq.1) then 
        if(param_num.eq.param_smap) then 
           var_found_smap = .true.
        endif
     endif
     
     if(var_found_smap) then
        call grib_get(igrib, 'values',sm_smap,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%rtsmopsnr
           do c=1,SMOPSsmobs(source)%rtsmopsnc
              sm_smap_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = &
                   sm_smap(c+((SMOPSsmobs(source)%rtsmopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%rtsmopsnc)
           enddo
        enddo     
        
     endif

     var_found_smap = .false. 
     if(SMOPSsmobs(source)%useSMAP.eq.1) then 
        if(param_num.eq.param_smap_qa) then 
           var_found_smap = .true.
        endif
     endif
     
     if(var_found_smap) then
        call grib_get(igrib, 'values',sm_smap_qa,iret)
        call LVT_warning(iret,'error in grib_get:values in readSMOPSsmObs')
        
        do r=1,SMOPSsmobs(source)%rtsmopsnr
           do c=1,SMOPSsmobs(source)%rtsmopsnc
              sm_smap_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = &
                   INT(sm_smap_qa(c+((SMOPSsmobs(source)%rtsmopsnr-r+1)-1)*&
                   SMOPSsmobs(source)%rtsmopsnc))
           enddo
        enddo       
     endif
     
     call grib_release(igrib,iret)
     call LVT_verify(iret, 'error in grib_release in readSMOPSsmObs')
  enddo

  call grib_close_file(ftn)

  if(SMOPSsmobs(source)%useASCAT.eq.1) then 
     do r=1, SMOPSsmobs(source)%rtsmopsnr
        do c=1, SMOPSsmobs(source)%rtsmopsnc
           if(sm_ascat_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc).ne.9999) then 
              !estimated error
              err = ISHFT(sm_ascat_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc),-8)
              !quality flag - not used currently
              ql = IAND(sm_ascat_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc),FF)
              
              if(err.lt.err_threshold) then 
                 sm_data_b(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = .true. 
              else
                 sm_data_b(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = .false.
                 sm_ascat_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = LVT_rc%udef
              endif
           else
              sm_ascat_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = LVT_rc%udef
              sm_data_b(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = .false. 
           endif
        enddo
     enddo
!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
     call neighbor_interp(LVT_rc%gridDesc(:),&
          sm_data_b, sm_ascat_t, smobs_b_ip, smobs_ip, &
          SMOPSsmobs(source)%rtsmopsnc*SMOPSsmobs(source)%rtsmopsnr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          SMOPSsmobs(source)%rlat, SMOPSsmobs(source)%rlon, &
          SMOPSsmobs(source)%n11, LVT_rc%udef, iret)
  endif

  if(SMOPSsmobs(source)%useWINDSAT.eq.1) then 
     do r=1, SMOPSsmobs(source)%rtsmopsnr
        do c=1, SMOPSsmobs(source)%rtsmopsnc
           if(sm_windsat_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc).ne.9999) then 
              !estimated error
              err = ISHFT(sm_windsat_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc),-8)
              !quality flag - not used currently
              ql = IAND(sm_windsat_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc),FF)
              
              if(err.lt.err_threshold) then 
                 sm_data_b(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = .true. 
              else
                 sm_data_b(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = .false.
                 sm_windsat_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = LVT_rc%udef
              endif
           else
              sm_windsat_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = LVT_rc%udef
              sm_data_b(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = .false. 
           endif
        enddo
     enddo
!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
     call neighbor_interp(LVT_rc%gridDesc(:),&
          sm_data_b, sm_windsat_t, smobs_b_ip, smobs_ip, &
          SMOPSsmobs(source)%rtsmopsnc*SMOPSsmobs(source)%rtsmopsnr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          SMOPSsmobs(source)%rlat, SMOPSsmobs(source)%rlon, &
          SMOPSsmobs(source)%n11, LVT_rc%udef, iret)
  endif

  if(SMOPSsmobs(source)%useSMOS.eq.1) then 
     do r=1, SMOPSsmobs(source)%rtsmopsnr
        do c=1, SMOPSsmobs(source)%rtsmopsnc
           if(sm_smos_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc).ne.9999) then 
              !estimated error
              err = ISHFT(sm_smos_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc),-8)
              !quality flag - not used currently
              ql = IAND(sm_smos_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc),FF)
              
              if(err.lt.err_threshold) then 
                 sm_data_b(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = .true. 
              else
                 sm_data_b(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = .false.
                 sm_smos_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = LVT_rc%udef
              endif
           else
              sm_smos_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = LVT_rc%udef
              sm_data_b(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = .false. 
           endif
        enddo
     enddo
!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 
     call neighbor_interp(LVT_rc%gridDesc(:),&
          sm_data_b, sm_smos_t, smobs_b_ip, smobs_ip, &
          SMOPSsmobs(source)%rtsmopsnc*SMOPSsmobs(source)%rtsmopsnr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          SMOPSsmobs(source)%rlat, SMOPSsmobs(source)%rlon, &
          SMOPSsmobs(source)%n11, LVT_rc%udef, iret)
  endif

!  open(100,file='test_inp1.bin',form='unformatted')
!  write(100) sm_smap_t
!  close(100)

  if(SMOPSsmobs(source)%useSMAP.eq.1) then 
     do r=1, SMOPSsmobs(source)%rtsmopsnr
        do c=1, SMOPSsmobs(source)%rtsmopsnc
!           if((sm_smap_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc).eq.9999).and.&
!                (sm_smap_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc).ne.9999)) then
!              print*, sm_smap_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc),&
!                   sm_smap_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc)
!           endif
           if(sm_smap_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc).ne.9999) then 
!           if(sm_smap_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc).ne.9999) then 
              !estimated error
              err = ISHFT(sm_smap_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc),-8)
              !quality flag - not used currently
              ql = IAND(sm_smap_qa_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc),FF)
              
!              if(err.lt.err_threshold) then 
                 sm_data_b(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = .true. 
!              else
!                 sm_data_b(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = .false.
!                 sm_smap_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = LVT_rc%udef
!              endif
           else
              sm_smap_t(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = LVT_rc%udef
              sm_data_b(c+(r-1)*SMOPSsmobs(source)%rtsmopsnc) = .false. 
           endif
        enddo
     enddo

!  open(100,file='test_inp.bin',form='unformatted')
!  write(100) sm_smap_t
!  close(100)
!  stop
!--------------------------------------------------------------------------
! Interpolate to the LVT running domain
!-------------------------------------------------------------------------- 

     call neighbor_interp(LVT_rc%gridDesc(:),&
          sm_data_b, sm_smap_t, smobs_b_ip, smobs_ip, &
          SMOPSsmobs(source)%rtsmopsnc*SMOPSsmobs(source)%rtsmopsnr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          SMOPSsmobs(source)%rlat, SMOPSsmobs(source)%rlon, &
          SMOPSsmobs(source)%n11, LVT_rc%udef, iret)
  endif

!  print*, 'writing '
!  open(100,file='test_ip.bin',form='unformatted')
!  write(100) smobs_ip
!  close(100)
!  stop
end subroutine read_SMOPS_data

!BOP
! !ROUTINE: create_SMOPSsm_filename
! \label{create_SMOPSsm_filename}
! 
! !INTERFACE: 
subroutine create_SMOPSsm_filename(ndir, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the SMOPS filename based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the SMOPS soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated RT SMOPS filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  filename = trim(ndir)//'/'//trim(fyr)//'/NPR_SMOPS_CMAP_D' &
       //trim(fyr)//trim(fmo)//trim(fda)//'.gr2'
  
end subroutine create_SMOPSsm_filename




