module nitrogen_data_mod

use fms_mod, only : &
     write_version_number, file_exist, open_namelist_file, check_nml_error, &
     close_file, stdlog, mpp_pe, mpp_root_pe, error_mesg, FATAL
use land_substances_mod, only : add_substance, inq_substance
use vegn_data_mod, only : NSPECIES, NCMPT, spdata

implicit none
private

! ==== public interfaces =====================================================
! ---- public constants
integer, public, parameter :: & ! types of nitrogen fixation calculations
     COMP_NFIX_PRESCRIBED  = 0, &
     COMP_NFIX_1           = 1, &
     COMP_NFIX_2           = 2, &
     COMP_NFIX_3           = 3

! ---- public data
public :: &
   do_nitrogen, do_nitrogen_plant_feedback, do_nitrogen_soil_feedback, &
   retrans_n, CtoN, n_fixrate_type, &
   a_lf, b_lf, f_lf_min, CtoN_seed, CtoN_litter_min, &
   lignin_lf, &
   soil_C_weight_content, kd_denit, Tr_denit, Q10_denit, St_denit, w_denit, &
   amm_buff, nitr_buff, DON_buff, &
   knup_plant, CtoN_SOM_slow, CtoN_SOM_res, &
   nu_decomp, eta_decomp, nitrif_rate, &
   fraction_dissolve, &
   nfixcost, k_phot, extinct_coeff, tau_fix, planning_horizon, &
   n_fixrate_option, nfixrate, ash_fraction, init_nfix_per_bliv, &
   digestive_waste_frac, n_deviation_cutoff, use_surface_runoff, &
   kmhalf, CtoN_target
integer, public :: isub_N = -1 ! substance index for nitrogen

! ---- public subroutine
public :: read_nitrogen_namelist

! ==== module constants ======================================================
character(len=*), private, parameter :: &
   version = '$Id: nitrogen_data.F90,v 1.1.2.9.2.2 2011/07/11 23:18:12 slm Exp $', &
   tagname = '$Name: nitro_20111003_slm $' ,&
   module_name = 'nitrogen_data_mod'
! ==== end of public interfaces ==============================================

integer :: idata,jdata ! iterators used in data initialization statements

! ---- namelist --------------------------------------------------------------
logical :: do_nitrogen                = .FALSE.
logical :: do_nitrogen_plant_feedback = .FALSE.
logical :: do_nitrogen_soil_feedback  = .FALSE.

character(len=32) :: n_fixrate_type = 'comp_nfix_1'

!       c4grass       c3grass    temp-decid      tropical     evergreen 
real :: retrans_n (0:NSPECIES-1) = &
     (/    0.40,         0.40,         0.40,         0.40,         0.30/)
real :: CtoN (0:NSPECIES-1,NCMPT); data ((CtoN(idata,jdata),idata=0,NSPECIES-1),jdata=1,NCMPT) &
      /     0.0,          0.0,          0.0,          0.0,          0.0, & ! reproduction - unused
          150.0,        150.0,        150.0,        150.0,        150.0, & ! sapwood
           35.0,         35.0,         35.0,         35.0,         50.0, & ! leaf
           50.0,         50.0,         50.0,         50.0,         60.0, & ! root
            0.0,          0.0,          0.0,          0.0,          0.0, & ! virtual leaf - unused
          500.0,        500.0,        500.0,        500.0,        500.0  / ! structural
real :: CtoN_seed = 30.0 ! seed C:N ratio. This is historical: it may be better to use 
              ! species values (CMPT_REPRO), but it's technically simpler right
              ! now to stick to the way it was implemented in LM3V. Request 
              ! Stefan's comment on that
real :: a_lf = 0.87 ! Fraction of litter-fall transferred to LF if lignin fraction is zero, unitless [Parton et al., 1993]
real :: b_lf = 0.0018 ! Change of LF fraction with the tissue lignin to N ratio, kgN/kgC [Parton et al., 1993]
real :: f_lf_min = 0.33 ! minimum fast litter fraction, unitless
real :: CtoN_litter_min = 150.0 ! minimum C:N ration in slow litter
real :: lignin_lf = 0.2 ! tissue lignin fraction, unitless [White et al., 2000]

real :: knup_plant(0:NSPECIES-1) = & ! N uptake rate constant of plants
     (/  3.1e5,         3.1e5,        1.8e4,        1.8e4,        1.8e4 /)
!real :: k_half_plant(N_NUTR)                          ! half saturation constant
real :: planning_horizon = 1.0  ! the time over which plant losses are buffered, years
!real :: phi(N_NUTR)                         ! parameter for down-regulation
!
real :: nfixcost = 9.2      ! carbon cost for nitrogen fixation; units=???
real :: nfixrate = -1       ! fixation rate for the n_fixrate_option=COMP_NFIX_PRESCRIBED, units=???
real :: k_phot   = 18.0     ! change rate of photosynthesis with N storage, unitless
real :: extinct_coeff = 0.5 ! light extinction coefficient for N fixers per unit LAI
                            ! why not use the actual extinction coefficient???
                            ! may be because it may not reflect the properties of 
                            ! fixers?
real :: tau_fix(0:NSPECIES-1) = & ! relaxation time for nitrogen fixation, years
     (/     1.0,          1.0,         20.0,         0.02,         20.0   /)
!
real :: soil_C_weight_content = 0.034 ! soil C content by weight, used for organic soil depth calculation
real :: amm_buff  = 10.0           ! buffering factor for ammonium
real :: nitr_buff = 1.0            ! buffering factor for nitrate
real :: DON_buff  = 10.0           ! stickiness factor, concentration in exported water is DON/don_buff
real :: kd_denit = 180 ! 
real :: Q10_denit = 2 ! 
real :: Tr_denit = 294.15 
real :: w_denit = 1 
real :: St_denit = 0 
!real :: p_buff
real :: nitrif_rate   = 180.0      ! base nitrification rate, 1/year
real :: CtoN_SOM_slow = 15.0       ! C:N ratio of slow soil organic matter compartment
real :: CtoN_SOM_res  = 10.0       ! C:N ratio of resistant soil organic matter compartment
!real :: mu_decomp                 ! N dependent decomposition rate of slow litter
!real :: k_half_n_immob            ! half saturation constant for stabilization
!real :: fraction_ndissolve(3)     ! fraction of N dissolved by litter decomposition and soil formation
real :: fraction_dissolve = 2.0e-2 ! fraction of decomposition of slow litter and soil pool that becomes soluble

real :: nu_decomp  = 0.0           ! parameter of N limitation in decomposition
real :: eta_decomp = 0.0           ! parameter of N influence on SOM production

real :: ash_fraction = 0.45
real :: init_nfix_per_bliv = 0.0   ! cold-start value 

real :: digestive_waste_frac = 0.0 ! fraction of harvested nitrogen that passes 
   ! through the digestive system of the animals very quickly and contributes 
   ! directly to litter (applied only to grazing)
real :: n_deviation_cutoff = 0.06  ! relative nitrogen status threshold for 
   ! merging vegetation tiles
logical :: use_surface_runoff=.FALSE. ! if true, surface runoff is included into drainage calculations

real :: kmhalf  = 1e-4   ! half saturation N concentration for immobilization into litter
real :: CtoN_target = 40.0 ! "target" C:N ration for immobilization into litter

namelist /land_nitrogen_nml/ &
   do_nitrogen, do_nitrogen_plant_feedback, do_nitrogen_soil_feedback, &
   retrans_n, CtoN, n_fixrate_type, nfixrate, CtoN_seed, &
   a_lf, b_lf, f_lf_min, CtoN_litter_min, &
   lignin_lf, &
   soil_C_weight_content, kd_denit, Tr_denit, Q10_denit, St_denit, w_denit, &
   amm_buff, nitr_buff, DON_buff, &
   knup_plant, CtoN_SOM_slow, CtoN_SOM_res, &
   nu_decomp, eta_decomp, nitrif_rate, &
   fraction_dissolve, &
   nfixcost, k_phot, extinct_coeff, tau_fix, planning_horizon, &
   ash_fraction, init_nfix_per_bliv, digestive_waste_frac, &
   n_deviation_cutoff, &
   use_surface_runoff, &
   kmhalf, CtoN_target

!---- end of namelist --------------------------------------------------------

! ==== module variables ======================================================
integer :: n_fixrate_option ! variable for efficient selection of the kind of 
                            ! nitrogen fixation calculations


contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! ============================================================================
subroutine read_nitrogen_namelist()
! ---- local variables
  integer :: unit, ierr, io
  integer :: i

  call write_version_number(version, tagname)
  if (file_exist('input.nml')) then
     unit = open_namelist_file()
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=land_nitrogen_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'land_nitrogen_nml')
     enddo
10   continue
     call close_file (unit)
  endif
  if (mpp_pe() == mpp_root_pe()) then
     unit = stdlog()
     write (unit, nml=land_nitrogen_nml)
     call close_file (unit)
  endif
  
  if (do_nitrogen) then
     call add_substance('N', 'nitrogen')
     call inq_substance('N', idx=isub_N )
  endif
  
  do i = 0,NSPECIES-1
     spdata(i)%retrans_n = retrans_n(i)
     spdata(i)%CtoN(:) = CtoN(i,:)
     spdata(i)%knup_plant = knup_plant(i)
     spdata(i)%tau_fix = tau_fix(i)
  enddo

  if (trim(n_fixrate_type)=='prescribed') then
     n_fixrate_option = COMP_NFIX_PRESCRIBED
  else if (trim(n_fixrate_type)=='comp_nfix_1') then
     n_fixrate_option = COMP_NFIX_1
  else if (trim(n_fixrate_type)=='comp_nfix_2') then
     n_fixrate_option = COMP_NFIX_2
  else if (trim(n_fixrate_type)=='comp_nfix_3') then
     n_fixrate_option = COMP_NFIX_3
  else
     call error_mesg('read_nutrogen_namelist','n_fixrate_type = "'//trim(n_fixrate_type)&
          //'" is invalid, use "prescribed", "comp_nfix_1", "comp_nfix_2", or "comp_nfix_3"',&
          FATAL)
  endif
  

end subroutine
end module nitrogen_data_mod
