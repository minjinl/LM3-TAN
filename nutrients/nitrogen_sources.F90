module nitrogen_sources_mod

use constants_mod, only : PI
use time_manager_mod, only : time_type, get_date, operator(/=), operator(-), &
     operator(<), valid_calendar_types, get_calendar_type, time_type_to_real
use fms_mod, only : &
    write_version_number, field_exist, file_exist, open_namelist_file, &
     check_nml_error, close_file, stdlog, mpp_pe, mpp_root_pe, error_mesg, &
     FATAL, NOTE, string
use time_interp_external_mod, only : init_external_field, time_interp_external, &
     time_interp_external_init
!use time_interp_external_mod, only : init_external_field, time_interp_external
use horiz_interp_mod, only : horiz_interp_type, horiz_interp_init, &
     horiz_interp_new, horiz_interp_del, horiz_interp
use time_interp_mod, only : time_interp
use get_cal_time_mod, only : get_cal_time

use nfu_mod, only : nfu_validtype, nfu_inq_var, nfu_get_dim_bounds, nfu_get_rec, &
     nfu_get_dim, nfu_get_valid_range, nfu_is_valid

use land_data_mod, only: land_state_type, lnd
use land_io_mod, only : read_field
use land_tile_io_mod, only : print_netcdf_error
use nitrogen_data_mod, only : do_nitrogen

implicit none
private

! ==== public interfaces =====================================================
public :: nitrogen_sources_init
public :: nitrogen_sources_end
public :: nitrogen_deposition
public :: nitrogen_fertilization
! ==== end of public interfaces ===============================================

! ==== module constants =======================================================
character(len=*), private, parameter :: &
   version = '$Id: nitrogen_sources.F90,v 1.1.2.5 2010/05/10 20:52:18 slm Exp $', &
   tagname = '$Name: nitro_20111003_slm $' ,&
   module_name = 'nitrogen_sources_mod', &
   diag_mod_name = 'vegn'

integer, parameter :: & ! types of deposition calculations
   NDEP_PRESCRIBED = 1, & ! single prescribed value
   NDEP_INTERP2    = 2, & ! interpolated between natural and anthropogenic
   NDEP_TIMESERIES = 3    ! regular time series

! ==== NetCDF declarations ===================================================
include 'netcdf.inc'
#define __NF_ASRT__(x) call print_netcdf_error((x),module_name,__LINE__)

! ==== module variables ======================================================

! ---- namelist --------------------------------------------------------------
real    :: wetdepfrac = 0.0 ! fraction of deposition associated with rainfall
real    :: ndeprate   = 0.0 ! prescribed deposition rate
character(len=32) :: n_deposition_type = 'interpolate-two-maps' ! or 'prescribed',
   ! or 'time-series'
integer :: year_nat   = 1860! year before which natural deposition rate is used
integer :: year_ant   = 1990! year after each anthropogenic deposition rate is used

logical :: do_crop_fertilization = .FALSE. ! if true, fertilization is applied to crops
logical :: do_past_fertilization = .FALSE. ! if true, fertilization is applied to pastures

character(512) :: &
   ndep_nit_file = 'INPUT/ndep_nit.nc', & ! input file name for nitrate deposition
   ndep_amm_file = 'INPUT/ndep_amm.nc', &              ! input file name for ammonium deposition
   ndep_don_file = 'INPUT/ndep_don.nc'                 ! input file name for DON deposition
character(512) :: &
   fm_nit_file = 'INPUT/fm_nit.nc', & ! input file name for nitrate deposition
   fm_amm_file = 'INPUT/fm_amm.nc', &              ! input file name for ammonium deposition
   fm_don_file = 'INPUT/fm_don.nc'                 ! input file name for DON deposition   
   
character(32) :: &
   ndep_nit_field = 'ndep_nit', &         ! input field name for nitrate deposition
   ndep_amm_field = 'ndep_amm', &             ! input field name for ammonium deposition
   ndep_don_field = 'ndep_don'                ! input field name for ammonium deposition
character(32) :: &
   fm_nit_field = 'fm_nit', &         ! input field name for nitrate deposition
   fm_amm_field = 'fm_amm', &             ! input field name for ammonium deposition
   fm_don_field = 'fm_don'                ! input field name for ammonium deposition
      
namelist/nitrogen_deposition_nml/n_deposition_type, &
   ndep_nit_file, ndep_nit_field, &
   ndep_amm_file, ndep_amm_field, &
   ndep_don_file, ndep_don_field, &
   fm_nit_file, fm_nit_field, &
   fm_amm_file, fm_amm_field, &
   fm_don_file, fm_don_field, &
   wetdepfrac, year_nat,year_ant,ndeprate, &
   do_crop_fertilization, do_past_fertilization
   
! ---- end of namelist -------------------------------------------------------

integer :: n_deposition_option = 0 !selector for deposition options
real, allocatable :: & ! buffers for input data
   ndep_nat(:,:),        & ! natural deposition rate
   ndep_ant(:,:),        & ! anthropogenic deposition rate
   ndep_nit_buffer(:,:), & ! nitrate deposition
   ndep_amm_buffer(:,:), & ! ammonium deposition
   ndep_don_buffer(:,:), &  ! DON source
   fm_nit_buffer(:,:), & ! nitrate deposition
   fm_amm_buffer(:,:), & ! ammonium deposition
   fm_don_buffer(:,:)    ! DON source
      
type(time_type) :: current_time ! current time for time series reading/interpolation
integer :: &
   id_ndep_nit = -1, & ! id of external field for time interpolation
   id_ndep_amm = -1, & ! id of external field for time interpolation
   id_ndep_don = -1, & ! id of external field for time interpolation
   id_fm_nit = -1, & ! id of external field for time interpolation
   id_fm_amm = -1, & ! id of external field for time interpolation
   id_fm_don = -1    ! id of external field for time interpolation
   
type :: fert_data_type ! type describing input fertilization data
   integer :: ncid = -1  ! netcdf ID of the input data file
   integer :: nlon_in = -1, nlat_in = -1 ! input grid size 
   type(time_type), pointer :: time_in(:) => NULL() ! input time axis
   real, pointer :: lonb_in(:,:) => NULL() ! input longitude boundaries
   real, pointer :: latb_in(:,:) => NULL() ! input latitude boundaries
   
   integer :: curr_rec = -1 ! number of the current record
   real, pointer :: nfert(:,:) => NULL() ! amount of fertilizer
   real, pointer :: tfert(:,:) => NULL() ! time of application
   real, pointer :: mask (:,:) => NULL() ! mask of valid data
   type(nfu_validtype) :: nfert_v ! valid range for amount of fertilizer
   type(nfu_validtype) :: tfert_v ! valid range for time of fertilizer
end type

type(fert_data_type) :: &
   fert_crop, & ! crop fertilization
   fert_past    ! pasture fertilization

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
 
! ============================================================================
subroutine nitrogen_sources_init(time)
  type(time_type), intent(in) :: time      ! current time

  ! ---- local vars
  integer :: unit, ierr, io
  
  call write_version_number(version, tagname)
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=nitrogen_deposition_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'nitrogen_deposition_nml')
     enddo
10   continue
     call close_file (unit)
  endif
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=nitrogen_deposition_nml)
     call close_file (unit)
  endif
  
  if (trim(n_deposition_type)=='interpolate-two-maps') then
    n_deposition_option = NDEP_INTERP2
  else if (trim(n_deposition_type)=='prescribed') then
    n_deposition_option = NDEP_PRESCRIBED
  else if (trim(n_deposition_type)=='time-series') then
    n_deposition_option = NDEP_TIMESERIES
  else
    call error_mesg('nitrogen_sources_init','n_deposition_type = "'//trim(n_deposition_type)&
          //'" is invalid, use "prescribed" or "interpolate-two-maps"',&
          FATAL)
  endif
  
  ! do nothing further if no nitrogen in the model
  if (.not.do_nitrogen) then
     call error_mesg('nitrogen_sources_init',&
       'do_nitrogen=.FALSE.: NOT initializing the nitrogen deposition sources', &
       NOTE)
     return
  endif
  
  ! read the data, if necessary
  if (n_deposition_option==NDEP_INTERP2) then
    allocate(&
        ndep_nat(lnd%is:lnd%ie,lnd%js:lnd%je), &
        ndep_ant(lnd%is:lnd%ie,lnd%js:lnd%je)  )
    call read_field( 'INPUT/nbiodata_nat.nc','ndep', &
          lnd%glon(lnd%is:lnd%ie,lnd%js:lnd%je), &
          lnd%glat(lnd%is:lnd%ie,lnd%js:lnd%je), &
          ndep_nat, interp='nearest')
    call read_field( 'INPUT/nbiodata_ant.nc','ndep', &
          lnd%glon(lnd%is:lnd%ie,lnd%js:lnd%je), &
          lnd%glat(lnd%is:lnd%ie,lnd%js:lnd%je), &
          ndep_ant, interp='nearest')
  endif
  
  if (n_deposition_option == NDEP_TIMESERIES) then
    ! allocate and initialize buffer for data
    allocate(ndep_nit_buffer(lnd%is:lnd%ie,lnd%js:lnd%je))
    allocate(ndep_amm_buffer(lnd%is:lnd%ie,lnd%js:lnd%je))
    allocate(ndep_don_buffer(lnd%is:lnd%ie,lnd%js:lnd%je))
    allocate(fm_nit_buffer(lnd%is:lnd%ie,lnd%js:lnd%je))
    allocate(fm_amm_buffer(lnd%is:lnd%ie,lnd%js:lnd%je))
    allocate(fm_don_buffer(lnd%is:lnd%ie,lnd%js:lnd%je))
    ndep_nit_buffer(:,:) = 0.0
    ndep_amm_buffer(:,:) = 0.0
    ndep_don_buffer(:,:) = 0.0
    fm_nit_buffer(:,:) = 0.0
    fm_amm_buffer(:,:) = 0.0
    fm_don_buffer(:,:) = 0.0
    
    ! initialize external fields 
    call time_interp_external_init()
    if (trim(ndep_nit_file) /= '') then
       ! initilaize external field
       id_ndep_nit = init_external_field(trim(ndep_nit_file),trim(ndep_nit_field),&
                                         domain=lnd%domain,use_comp_domain=.TRUE.)
       ! read data for current time
       call time_interp_external(id_ndep_nit,time,ndep_nit_buffer)
    endif
    if (trim(ndep_amm_file) /= '') then
       ! initilaize external field
              id_ndep_amm = init_external_field(trim(ndep_amm_file),trim(ndep_amm_field),&
                                         domain=lnd%domain,use_comp_domain=.TRUE.)
       ! read data for current time
       call time_interp_external(id_ndep_amm,time,ndep_amm_buffer)
    endif
    if (trim(ndep_don_file) /= '') then
       ! initilaize external field
       id_ndep_don = init_external_field(trim(ndep_don_file),trim(ndep_don_field),&
                                         domain=lnd%domain,use_comp_domain=.TRUE.)
       ! read data for current time
       call time_interp_external(id_ndep_don,time,ndep_don_buffer)
    endif
    if (trim(fm_nit_file) /= '') then
       ! initilaize external field
       id_fm_nit = init_external_field(trim(fm_nit_file),trim(fm_nit_field),&
                                         domain=lnd%domain,use_comp_domain=.TRUE.)
       ! read data for current time
       call time_interp_external(id_fm_nit,time,fm_nit_buffer)
    endif
    if (trim(fm_amm_file) /= '') then
       ! initilaize external field
              id_fm_amm = init_external_field(trim(fm_amm_file),trim(fm_amm_field),&
                                         domain=lnd%domain,use_comp_domain=.TRUE.)
       ! read data for current time
       call time_interp_external(id_fm_amm,time,fm_amm_buffer)
    endif
    if (trim(fm_don_file) /= '') then
       ! initilaize external field
       id_fm_don = init_external_field(trim(fm_don_file),trim(fm_don_field),&
                                         domain=lnd%domain,use_comp_domain=.TRUE.)
       ! read data for current time
       call time_interp_external(id_fm_don,time,fm_don_buffer)
    endif
    
    ! set the current time
    current_time = time
  endif

  if (do_crop_fertilization) &
     call fert_data_init("INPUT/nfert_crop.nc",fert_crop)
  if (do_past_fertilization) &
     call fert_data_init("INPUT/nfert_past.nc",fert_past)

end subroutine nitrogen_sources_init


! ============================================================================
subroutine nitrogen_sources_end()
  ! clean up deposition fields
  if (allocated(ndep_nat)) deallocate(ndep_nat)
  if (allocated(ndep_ant)) deallocate(ndep_ant)
  if (allocated(ndep_nit_buffer)) deallocate(ndep_nit_buffer)
  if (allocated(ndep_amm_buffer)) deallocate(ndep_amm_buffer)
  if (allocated(ndep_don_buffer)) deallocate(ndep_don_buffer)
  if (allocated(fm_nit_buffer)) deallocate(fm_nit_buffer)
  if (allocated(fm_amm_buffer)) deallocate(fm_amm_buffer)
  if (allocated(fm_don_buffer)) deallocate(fm_don_buffer)
  call fert_data_end(fert_crop)
  call fert_data_end(fert_past)
end subroutine


! ============================================================================
subroutine nitrogen_deposition(time, i, j, p_ann, precip, ndep_nit, ndep_amm, ndep_don, fm_nit, fm_amm, fm_don)
  type(time_type), intent(in)  :: time ! current time, for time interpolation
  integer, intent(in)  :: i,j    ! coordinates of current point
  real   , intent(in)  :: p_ann  ! annual-mean precipitation, kg/(m2 s)
  real   , intent(in)  :: precip ! current precipitation, kg/(m2 s)
  real   , intent(out) :: ndep_nit ! nitrate deposition rate, kg N/(m2 yr)
  real   , intent(out) :: ndep_amm ! ammonium deposition rate, kg N/(m2 yr)
  real   , intent(out) :: ndep_don ! DON deposition rate, kg N/(m2 yr)
  real   , intent(out) :: fm_nit ! nitrate deposition rate, kg N/(m2 yr)
  real   , intent(out) :: fm_amm ! ammonium deposition rate, kg N/(m2 yr)
  real   , intent(out) :: fm_don ! DON deposition rate, kg N/(m2 yr)

  real :: w ! averaging weight
  integer :: second, minute, hour, day, month, year
  
  if (n_deposition_option == NDEP_PRESCRIBED) then
    ndep_nit = ndeprate
  else if (n_deposition_option == NDEP_INTERP2) then
    call get_date(time, year,month,day,hour,minute,second)
    w = real(year - year_nat)/real(year_ant-year_nat)
    w = min(w,1.0); w = max(w,0.0)
    ndep_nit = (1-w)*ndep_nat(i,j)+w*ndep_ant(i,j)
  else if (n_deposition_option == NDEP_TIMESERIES) then
    if (time/=current_time) then
      ! get the proper data from external file. The majority of calls
      ! are done with current_time == time, so that time_interp_external is
      ! called only every fast time step, not for every i,j
      if (id_ndep_nit>0) call time_interp_external(id_ndep_nit,time,ndep_nit_buffer)
      if (id_ndep_amm>0) call time_interp_external(id_ndep_amm,time,ndep_amm_buffer)
      if (id_ndep_don>0) call time_interp_external(id_ndep_don,time,ndep_don_buffer)
      if (id_fm_nit>0) call time_interp_external(id_fm_nit,time,fm_nit_buffer)
      if (id_fm_amm>0) call time_interp_external(id_fm_amm,time,fm_amm_buffer)
      if (id_fm_don>0) call time_interp_external(id_fm_don,time,fm_don_buffer)
      current_time = time
    endif
    ndep_nit = ndep_nit_buffer(i,j)
    ndep_amm = ndep_amm_buffer(i,j)
    ndep_don = ndep_don_buffer(i,j)
    fm_nit = fm_nit_buffer(i,j)
    fm_amm = fm_amm_buffer(i,j)
    fm_don = fm_don_buffer(i,j)
  else
    call error_mesg('nitrogen_deposition','n_deposition_option '//&
          string(n_deposition_option)//' is invalid"', FATAL)
  endif
  ! scale the deposition by precipitation
  if (p_ann>0) then
     ndep_nit = ndep_nit*(1-wetdepfrac)+wetdepfrac*ndep_nit*precip/p_ann
     ndep_amm = ndep_amm*(1-wetdepfrac)+wetdepfrac*ndep_amm*precip/p_ann
  endif

end subroutine nitrogen_deposition


! ============================================================================
subroutine fert_data_init(filename, fert)
  character(*), intent(in) :: filename
  type(fert_data_type), intent(inout) :: fert

  ! ---- local vars
  integer :: dimids(NF_MAX_VAR_DIMS) ! IDs of dimensions
  integer :: dimlens(NF_MAX_VAR_DIMS) ! sizes of dimensions
  integer :: timedim ! id of the record (time) dimension
  integer :: timevar ! id of the time variable
  character(len=NF_MAX_NAME) :: timename  ! name of the time variable
  character(len=256)         :: timeunits ! units of time in the file
  character(len=24) :: calendar ! model calendar
  real, allocatable :: time(:)  ! time line from fertilization file
  integer :: i
  integer :: ierr ! error code
  integer :: nrec ! number of records in the fertilization data
  
  ierr=nf_open(filename,NF_NOWRITE,fert%ncid)

  if(ierr/=NF_NOERR) call error_mesg('nitrogen_sources_init', &
       'do_fertilization is requested, but fertilization "'// &
       trim(filename)//'" could not be opened because '//nf_strerror(ierr), &
       FATAL)

  __NF_ASRT__(nfu_inq_var(fert%ncid,'nfert',dimids=dimids,dimlens=dimlens,nrec=nrec))
  allocate(time(nrec), fert%time_in(nrec))
  
  ! read fertilization time line 
  ! get the time axis 
  __NF_ASRT__(nf_inq_unlimdim(fert%ncid, timedim))
  __NF_ASRT__(nfu_get_dim(fert%ncid, timedim, time))
  ! get units of time
  __NF_ASRT__(nf_inq_dimname(fert%ncid, timedim, timename))
  __NF_ASRT__(nf_inq_varid(fert%ncid,timename, timevar))
  timeunits = ' '
  __NF_ASRT__(nf_get_att_text(fert%ncid,timevar,'units',timeunits))
  ! get model calendar
  calendar=valid_calendar_types(get_calendar_type())

  ! loop through the time axis and store time_type values in time_in
  do i = 1,size(time)
     fert%time_in(i) = get_cal_time(time(i),timeunits,calendar)
  end do
  
  ! get the boundaries of the horizontal axes
  allocate(fert%lonb_in(dimlens(1)+1,1), &
           fert%latb_in(1,dimlens(2)+1)  )
  __NF_ASRT__(nfu_get_dim_bounds(fert%ncid, dimids(1), fert%lonb_in(:,1)))
  __NF_ASRT__(nfu_get_dim_bounds(fert%ncid, dimids(2), fert%latb_in(1,:)))
  ! convert them to radian
  fert%lonb_in = fert%lonb_in*PI/180.0
  fert%latb_in = fert%latb_in*PI/180.0
  ! set up the input buffer sizes
  fert%nlon_in = dimlens(1)
  fert%nlat_in = dimlens(2)
  ! set up the valid ranges
  __NF_ASRT__(nfu_get_valid_range(fert%ncid,'nfert',fert%nfert_v))
  ierr = nfu_get_valid_range(fert%ncid,'tfert',fert%tfert_v) ! ignore errors here
  
  ! allocate buffers for data
  allocate(fert%nfert(lnd%is:lnd%ie,lnd%js:lnd%je))
  allocate(fert%tfert(lnd%is:lnd%ie,lnd%js:lnd%je))
  ! and mask
  allocate(fert%mask (lnd%is:lnd%ie,lnd%js:lnd%je))
  ! set the current record indicator to non-existing record number
  fert%curr_rec = -1
  
  ! deallocate temporary variables
  deallocate(time)

end subroutine


! ============================================================================
#define __DEALLOC__(x) if (associated(x)) deallocate(x)
subroutine fert_data_end(fert)
  type(fert_data_type), intent(inout) :: fert
  
  integer :: ierr ! error code
  
  __DEALLOC__(fert%time_in)
  __DEALLOC__(fert%lonb_in)
  __DEALLOC__(fert%latb_in)
  __DEALLOC__(fert%nfert)
  __DEALLOC__(fert%tfert)
  __DEALLOC__(fert%mask)
  fert%curr_rec=-1
  
  ierr = nf_close(fert%ncid) ! ignore errors here
  fert%ncid = -1
end subroutine


! =============================================================================
! given interval [t0,t1), returns for each grid point fertilization that is to
! be applied during this interval, in kg N/(m2 year)
subroutine fert_data_get(fert, t0,t1,fert_amount)
  type(fert_data_type), intent(inout) :: fert  ! fertilization data
  type(time_type)     , intent(in)    :: t0,t1 ! time interval
  real                , intent(out)   :: fert_amount(:,:) 

  ! ---- local vars
  integer :: i1,i2 ! record indices
  integer :: ierr ! error code
  integer :: i,j
  real    :: w
  real    :: t0r, t1r ! real values of t0,t1, days since the beginning of
                      ! current time interval in the input data
  real, allocatable :: tfert_in(:,:) ! input buffer for tfert
  real, allocatable :: nfert_in(:,:) ! input buffer for nfert
  real, allocatable :: mask_in(:,:)  ! mask for interpolation
  type(horiz_interp_type) :: interp

  ! set initial value
  fert_amount(:,:) = 0
  
  ! get the time input data time interval for the current date
  call time_interp(t0, fert%time_in, w, i1, i2)
  ! make sure we have the proper record in the memory
  if (fert%curr_rec/=i1) then
     ! we need to read/interpolate  the record.
     allocate(tfert_in(fert%nlon_in,fert%nlat_in))
     allocate(nfert_in(fert%nlon_in,fert%nlat_in))
     allocate(mask_in (fert%nlon_in,fert%nlat_in))
     tfert_in=0
     nfert_in=0
     mask_in=0
     ! We have to construct interpolator every time, because the mask may change
     ! for every record

     ! read fertilization time
     ierr = nfu_get_rec(fert%ncid,"tfert",i1,tfert_in)
     if(ierr==NF_NOERR) then
        where (nfu_is_valid(tfert_in,fert%tfert_v))
          mask_in = 1
        elsewhere
          mask_in = 0
        end where
     else
        tfert_in = 0
        mask_in  = 1
     endif
     
     ! read fertilization amount
     __NF_ASRT__(nfu_get_rec(fert%ncid,"nfert",i1,nfert_in))
     ! mask is 1 only where both tfert and nfert are valid
     where (.not.nfu_is_valid(nfert_in,fert%nfert_v)) &
           mask_in = mask_in*0

     ! construct interpolator
     call horiz_interp_new(interp, fert%lonb_in, fert%latb_in, &
          lnd%glonb(lnd%is:lnd%ie+1,lnd%js:lnd%je+1), &
          lnd%glatb(lnd%is:lnd%ie+1,lnd%js:lnd%je+1), &
          interp_method='conservative',&
          mask_in=mask_in, is_latlon_in=.TRUE. )
     ! note that mask_out does not seem to be set if both input and output grids
     ! are lat-lon (horiz_interp_new_2d, line 546). So it's pretty useless in
     ! this case.
     
     ! interpolate input fields
     call horiz_interp(interp,nfert_in,fert%nfert)
     call horiz_interp(interp,tfert_in,fert%tfert)
     call horiz_interp(interp,mask_in,fert%mask)
     ! convert time to seconds (from days)
     fert%tfert = fert%tfert*86400
     where(fert%mask==0)
        fert%tfert = -9999.0
        fert%nfert = 0.0
     end where
     
     ! update the current record number
     fert%curr_rec = i1
     
     ! clean up memory
     deallocate(tfert_in,nfert_in,mask_in)
     call horiz_interp_del(interp)
  endif
  
  ! convert time intervals to real numbers: note that in FMS time is positive
  ! by definition [why, oh why???] so that t0-t1 == t1-t0, and an extra step is 
  ! necessary to get the right sign of the real time intervals...
  t0r = time_type_to_real(t0-fert%time_in(i1))
  if (t0<fert%time_in(i1)) t0r = -t0r
  t1r = time_type_to_real(t1-fert%time_in(i1))
  if (t1<fert%time_in(i1)) t1r = -t1r
  
  ! set up the fertilization values
  fert_amount = 0.0
  where(t0r<=fert%tfert.and.fert%tfert<t1r) &
      fert_amount = fert%nfert
end subroutine


! =============================================================================
subroutine nitrogen_fertilization(t0,t1,crop_fert,past_fert)
  type(time_type)     , intent(in)    :: t0,t1 ! time interval
  real                , intent(out)   :: crop_fert(:,:) ! crop fertilization input
  real                , intent(out)   :: past_fert(:,:) ! pasture fertilization input
  
  crop_fert = 0
  past_fert = 0
  if (do_crop_fertilization) &
     call fert_data_get(fert_crop,t0,t1,crop_fert)
  if (do_past_fertilization) &
     call fert_data_get(fert_past,t0,t1,past_fert)
end subroutine

end module nitrogen_sources_mod
