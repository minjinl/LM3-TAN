module nitrogen_mod

#include "../shared/debug.inc"

use constants_mod, only: dens_h2o
use fms_mod, only : write_version_number, error_mesg, FATAL, WARNING
use time_manager_mod, only: time_type, time_type_to_real, get_date
use land_constants_mod, only : seconds_per_year
use land_data_mod,      only : land_state_type, lnd
use land_debug_mod,     only : is_watch_point, get_current_point
use land_tile_diag_mod, only : &
     register_tiled_diag_field, send_tile_data, diag_buff_type
use vegn_cohort_mod, only : vegn_cohort_type, A_function
use vegn_tile_mod, only : vegn_tile_type, n_store_optimum

use vegn_data_mod, only : CMPT_LEAF, CMPT_ROOT, CMPT_WOOD, CMPT_SAPWOOD, &
       SP_C3GRASS, SP_C4GRASS, SP_TROPICAL, LEAF_ON, LEAF_OFF, &
       spdata, SOMinput_frac_slow, SOMinput_frac_res, &
       k_litter_fast, k_litter_slow, k_SOM_slow, k_SOM_res

use nitrogen_data_mod, only: &
   COMP_NFIX_PRESCRIBED, COMP_NFIX_1, COMP_NFIX_2, COMP_NFIX_3, &
   isub_N,&
   a_lf, b_lf, f_lf_min, CtoN_litter_min, &
   lignin_lf, &
   soil_C_weight_content, kd_denit, Q10_denit, Tr_denit, w_denit, St_denit, &
   amm_buff, nitr_buff, DON_buff, &
   knup_plant, CtoN_SOM_slow, CtoN_SOM_res, &
   nu_decomp, eta_decomp, nitrif_rate, &
   fraction_dissolve, &
   nfixcost, k_phot, extinct_coeff, &
   n_fixrate_option, nfixrate, ash_fraction, &
   use_surface_runoff, &
   kmhalf, CtoN_target

implicit none
private

! ==== public interfaces =====================================================
public :: nitrogen_init
public :: getCtoNleaf, getCtoNroot, getCtoNwood, getCtoNsw
public :: litterfall, n_md, n_phot, n_fixation, n_excrete
public :: n_reallocate
public :: n_fire
public :: n_soil_depth
public :: n_soil_int
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), private, parameter :: &
   version = '$Id: nitrogen.F90,v 1.1.2.19.2.9.2.1 2011/10/04 00:03:24 slm Exp $', &
   tagname = '$Name: nleach_experimental_pcm_m2l_slm $' ,&
   module_name = 'nitrogen_mod', &
   diag_mod_name = 'vegn'


! ==== module variables ======================================================
real    :: dt_fast_yr ! fast (physical) time step, yr (year is defined as 365 days)

! diagnostic field IDs
integer :: &
     id_n_amm, id_n_nit, &
     id_netmin, id_nitrif, id_ammleach, id_nitrleach, id_nleach, id_donleach, &
     id_nuprate, id_ndep, id_ndep_nit, id_ndep_amm, id_ndep_don, id_nfert, id_nfm_nit, id_nfm_amm, id_nfm_don, &
     id_sorg_depth, id_sorg_T, id_sorg_theta, id_sorg_uptk, id_sorg_drain, &
     id_sorg_A, &
     id_immob, id_immob_SS, id_immob_SP, id_immob_LS, id_immob_LF, &
     id_decomp_C_SS, id_decomp_C_SP, id_decomp_C_LS, id_decomp_C_LF, & 
     id_nitrogen_in, id_nitrogen_out, id_n_conc, id_denit_rate

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! ============================================================================
subroutine nitrogen_init( id_lon, id_lat )
  integer        , intent(in) :: id_lon  ! ID of land longitude (X) axis  
  integer        , intent(in) :: id_lat  ! ID of land latitude (Y) axis
  
  call write_version_number(version, tagname)

  ! calculate fast time step, in years
  dt_fast_yr = time_type_to_real(lnd%dt_fast)/seconds_per_year
  
  ! initialize nitrogen diagnostic fields
  call nitrogen_diag_init(id_lon, id_lat, lnd%time)
  
end subroutine

! ============================================================================
subroutine nitrogen_diag_init(id_lon, id_lat, time)
  integer        , intent(in) :: id_lon  ! ID of land longitude (X) axis  
  integer        , intent(in) :: id_lat  ! ID of land latitude (Y) axis
  type(time_type), intent(in) :: time    ! initial time for diagnostic fields

  id_n_amm = register_tiled_diag_field ( diag_mod_name, 'n_amm', &
       (/id_lon, id_lat/), time,  'soil ammonium', 'kg N/m2', &
       missing_value=-999.0 )
  id_n_nit = register_tiled_diag_field ( diag_mod_name, 'n_nitr', &
       (/id_lon, id_lat/), time,  'soil nitrate', 'kg N/m2', &
       missing_value=-999.0 )
  id_n_conc = register_tiled_diag_field ( diag_mod_name, 'n_conc', &
       (/id_lon, id_lat/), time,  'available soil nitrogen concentration', 'kg N/m3', &
       missing_value=-999.0 )


  id_netmin = register_tiled_diag_field ( diag_mod_name, 'netmin',  &
       (/id_lon, id_lat/), time,  'net mineralization rate', 'kg N/(m2 year)', &
       missing_value=-999.0 )
  id_nitrif = register_tiled_diag_field ( diag_mod_name, 'nitrif',  &
       (/id_lon,id_lat/), time, 'nitrification rate', 'kg N/(m2 year)', &
       missing_value=-999.0 )
  id_nleach = register_tiled_diag_field ( diag_mod_name, 'nleach',  &
       (/id_lon,id_lat/), time, 'soil mineral nitrogen leaching rate', &
       'kg N/(m2 year)', missing_value=-999.0 )
  id_ammleach = register_tiled_diag_field ( diag_mod_name, 'amm_leach',  &
       (/id_lon,id_lat/), time, 'soil mineral ammonium leaching rate', &
       'kg N/(m2 year)', missing_value=-999.0 )
  id_nitrleach = register_tiled_diag_field ( diag_mod_name, 'nitr_leach',  &
       (/id_lon,id_lat/), time, 'soil mineral nitrate leaching rate', &
       'kg N/(m2 year)', missing_value=-999.0 )
  id_donleach = register_tiled_diag_field ( diag_mod_name, 'donleach',  &
       (/id_lon, id_lat/), time,  'dissolved organic nitrogen leaching rate', &
       'kg N/(m2 year)', missing_value=-999.0 )
  id_nuprate = register_tiled_diag_field ( diag_mod_name, 'nuprate',  &
       (/id_lon, id_lat/), time,  'plant nitrogen uptake', &
       'kg N/(m2 year)', missing_value=-999.0 )

  id_ndep = register_tiled_diag_field (diag_mod_name, 'ndep', (/id_lon,id_lat/), &
       time, 'total nitrogen deposition', 'kg N/(m2 year)', missing_value=-999.0)
  id_ndep_nit = register_tiled_diag_field (diag_mod_name, 'ndep_nit', (/id_lon,id_lat/), &
       time, 'nitrate deposition', 'kg N/(m2 year)', missing_value=-999.0)
  id_ndep_amm = register_tiled_diag_field (diag_mod_name, 'ndep_amm', (/id_lon,id_lat/), &
       time, 'ammonium deposition', 'kg N/(m2 year)', missing_value=-999.0)
  id_ndep_don = register_tiled_diag_field (diag_mod_name, 'ndep_don', (/id_lon,id_lat/), &
       time, 'DON input', 'kg N/(m2 year)', missing_value=-999.0)
  id_nfert = register_tiled_diag_field (diag_mod_name, 'nfert', (/id_lon,id_lat/), &
       time, 'nitrogen fertilization', 'kg N/(m2 year)', missing_value=-999.0)
  id_nfm_nit = register_tiled_diag_field (diag_mod_name, 'nfm_nit', (/id_lon,id_lat/), &
       time, 'nitrogen fertilization', 'kg N/(m2 year)', missing_value=-999.0)
  id_nfm_amm = register_tiled_diag_field (diag_mod_name, 'nfm_amm', (/id_lon,id_lat/), &
       time, 'nitrogen fertilization', 'kg N/(m2 year)', missing_value=-999.0)
  id_nfm_don = register_tiled_diag_field (diag_mod_name, 'nfm_don', (/id_lon,id_lat/), &
       time, 'nitrogen fertilization', 'kg N/(m2 year)', missing_value=-999.0)
  id_sorg_depth = register_tiled_diag_field (diag_mod_name, 'sorg_depth', (/id_lon,id_lat/), &
       time, 'depth of organic soil', 'm', missing_value=-999.0)
  id_sorg_T = register_tiled_diag_field (diag_mod_name, 'sorg_T', (/id_lon,id_lat/), &
       time, 'average temperature of organic soil', 'degK', missing_value=-999.0)
  id_sorg_theta = register_tiled_diag_field (diag_mod_name, 'sorg_theta', (/id_lon,id_lat/), &
       time, 'average moisture of organic soil', missing_value=-999.0)
  id_sorg_uptk = register_tiled_diag_field (diag_mod_name, 'sorg_uptk', (/id_lon,id_lat/), &
       time, 'uptake from organic soil', 'kg/(m2 s)', missing_value=-999.0)
  id_sorg_drain = register_tiled_diag_field (diag_mod_name, 'sorg_drain', (/id_lon,id_lat/), &
       time, 'drainage from organic soil', 'kg/(m2 s)', missing_value=-999.0)
  id_sorg_A = register_tiled_diag_field (diag_mod_name, 'sorg_A', (/id_lon,id_lat/), &
       time, 'decomposition factor due to organic soil temperature and moisture', &
       missing_value=-999.0)
  id_nitrogen_in = register_tiled_diag_field (diag_mod_name, 'nitrogen_in', (/id_lon,id_lat/), &
       time, 'nitrogen input accumulator', 'kg N/m2', missing_value=-999.0)
  id_nitrogen_out = register_tiled_diag_field (diag_mod_name, 'nitrogen_out', (/id_lon,id_lat/), &
       time, 'nitrogen output accumulator', 'kg N/m2', missing_value=-999.0)
   id_denit_rate = register_tiled_diag_field (diag_mod_name, 'denit_rate', (/id_lon,id_lat/), &
       time, 'denitrification rate', 'kgN/(m2 year)', missing_value=-999.0)
  id_immob =  register_tiled_diag_field (diag_mod_name, 'immob_N', (/id_lon,id_lat/), &
       time, 'rate of mineral N immobilization', 'kg N/(m2 year)', missing_value=-999.0)
  id_immob_SS =  register_tiled_diag_field (diag_mod_name, 'immob_N_SS', (/id_lon,id_lat/), &
       time, 'rate of mineral N immobilization to slow soil pool', 'kg N/(m2 year)', missing_value=-999.0)
  id_immob_SP =  register_tiled_diag_field (diag_mod_name, 'immob_N_SP', (/id_lon,id_lat/), &
       time, 'rate of mineral N immobilization to passive (recalcitrant) soil pool', 'kg N/(m2 year)', missing_value=-999.0)
  id_immob_LS =  register_tiled_diag_field (diag_mod_name, 'immob_N_LS', (/id_lon,id_lat/), &
       time, 'rate of mineral N immobilization to slow litter pool', 'kg N/(m2 year)', missing_value=-999.0)
  id_immob_LF =  register_tiled_diag_field (diag_mod_name, 'immob_N_LF', (/id_lon,id_lat/), &
       time, 'rate of mineral N immobilization to fast litter pool', 'kg N/(m2 year)', missing_value=-999.0)

  id_decomp_C_SS = register_tiled_diag_field (diag_mod_name, 'decomp_C_SS', (/id_lon,id_lat/), &
       time, 'rate of slow soil C decomposition', 'kg C/(m2 year)', missing_value=-999.0)
  id_decomp_C_SP = register_tiled_diag_field (diag_mod_name, 'decomp_C_SP', (/id_lon,id_lat/), &
       time, 'rate of passive (recalcitrant) soil C decomposition', 'kg C/(m2 year)', missing_value=-999.0)
  id_decomp_C_LS = register_tiled_diag_field (diag_mod_name, 'decomp_C_LS', (/id_lon,id_lat/), &
       time, 'rate of slow litter C decomposition', 'kg C/(m2 year)', missing_value=-999.0)
  id_decomp_C_LF = register_tiled_diag_field (diag_mod_name, 'decomp_C_LF', (/id_lon,id_lat/), &
       time, 'rate of fast litter C decomposition', 'kg C/(m2 year)', missing_value=-999.0)
end subroutine

! ============================================================================
function getCtoNleaf(cc) result(res)
  type(vegn_cohort_type), intent(in) :: cc
  real :: res ! return value

  res = spdata(cc%species)%CtoN(CMPT_LEAF)
  if (cc%bliving < 1.e-15) return

  if (cc%status == LEAF_ON .and. cc%nl > 0) then
    res = cc%bl/cc%nl
  else if (cc%status == LEAF_OFF .and. cc%nliving + cc%nstore > 0.) then
    res = cc%bliving/(cc%nliving + cc%nstore);
  else
    ! return the value from the data (namelist) 
  endif
end function

! ============================================================================
function getCtoNroot(cc) result(res)
  type(vegn_cohort_type), intent(in) :: cc
  real :: res ! return value
   
  res = spdata(cc%species)%CtoN(CMPT_ROOT)
  if (cc%bliving < 1.e-15) return

  if (cc%status == LEAF_ON .and. cc%nr > 0.) then
    res = cc%br/cc%nr;
  else if (cc%status == LEAF_OFF .and.  cc%nliving + cc%nstore > 0.) then
    res = cc%bliving/(cc%nliving + cc%nstore);
  else
    ! return the value from the data (namelist) 
  endif
end function

! ===========================================================================
function getCtoNwood(cc) result(res)
  type(vegn_cohort_type), intent(in) :: cc
  real :: res ! return value
   
  if (cc%nwood*cc%bwood > 0.) then
    res = cc%bwood/cc%nwood;
  else
    res = spdata(cc%species)%CtoN(CMPT_WOOD)
  endif
end function

! ===========================================================================
function getCtoNsw(cc) result(res)
  type(vegn_cohort_type), intent(in) :: cc
  real :: res

  res=spdata(cc%species)%CtoN(CMPT_SAPWOOD)
  if (cc%bliving < 1.e-15) return

  if (cc%status == LEAF_ON .and. cc%nsw > 0) then
    res=cc%bsw/(cc%nsw + MAX(cc%nstore,0.));
  else if (cc%status == LEAF_OFF .and. cc%nstore + cc%nliving > 0) then
    res=cc%bliving/(cc%nliving + cc%nstore);
  else
    ! return the value from the data (namelist) 
  endif
end function


! ===========================================================================
! calculates reduction of photosynthesis due to nitrogen limitation
real function n_phot(nsf)
  real, intent(in) :: nsf

  n_phot = 1.0 - exp(-MAX(MIN(nsf,1.),0.)*k_phot)
end function

! ============================================================================
! returns active soil depth for nitrogen calculations: uptake and leaching of
! nitrogen and conversion ov N amount to concentrations.
function n_soil_depth(vegn) result(res)
  type(vegn_tile_type), intent(in) :: vegn
  real :: res ! return value
  
  ! The soil depth is calculated assuming C weight content rc = 3.4% and average 
  ! soil density 1500 kg/m3 (main text, [23] on page 6)
  res = (vegn%fast_litter_C + vegn%slow_litter_C + vegn%slow_soil_C)/(soil_C_weight_content*1500) 
  res = MAX(res,0.05) ! to avoid singularity
  res = MIN(res,0.4)  ! truncate organic soil depth to avoid temperature-decomp
                      ! feedbacks in high latitudes
end function


! ============================================================================
! dynamics of mineral nitrogen pools
subroutine n_soil_int(vegn, diag, soil_depth, soilt, theta, &
     ndep_nit, ndep_amm, ndep_don, nfert, nfm_nit, nfm_amm, nfm_don, transp, drain, surf_runoff)
  type(vegn_tile_type), intent(inout) :: vegn
  type(diag_buff_type), intent(inout) :: diag
  real, intent(in) :: soil_depth ! depth of organic soil, m
  real, intent(in) :: soilt    ! soil temperature averaged over organic soil, degK
  real, intent(in) :: theta    ! soil moisture averaged over organic soil, unitless
  real, intent(in) :: ndep_nit, ndep_amm, ndep_don ! nitrate, ammonium, and DON deposition rates, kg N/(m2 yr)
  real, intent(in) :: nfert    ! nitrogen fertilization rate, kg N/yr
  real, intent(in) :: nfm_nit
  real, intent(in) :: nfm_amm
  real, intent(in) :: nfm_don
  real, intent(in) :: transp   ! uptake from organic soil, kg/(m2 s)
  real, intent(in) :: drain    ! drainage (water export other than uptake) from organic soil, kg/(m2 s)
  real, intent(in) :: surf_runoff ! surface runoff, kg/(m2 s)

  ! ---- constants
  real, parameter :: kphalf = 3.e-3  ! half saturation constant for plant uptake [kg/m3]
  real, parameter :: mobility = 1.0
  
  ! ---- local vars
  real :: amm_conc, nit_conc, n_conc ! concentrations of ammonia, nitrate, and total [kg N/m3]
  ! rate modifier factors and their derivatives
  real :: n_factor0, Dn_factorDn, n_factor ! influence on nitrogen on decomposition
  real :: q_factor0, Dq_factorDn, q_factor ! fraction of heavy material
  ! rates of the change of mineral nitrogen and their derivatives
  !       expl.estimate  deriv.W.r.t.amm     deriv.W.r.t.nit     updated value
  real :: amm_uptake0,   Damm_uptakeDamm,    Damm_uptakeDnit,    amm_uptake
  real :: nit_uptake0,   Dnit_uptakeDamm,    Dnit_uptakeDnit,    nit_uptake
  real :: amm_leach0,    Damm_leachDamm,     Damm_leachDnit,     amm_leach
  real :: nit_leach0,    Dnit_leachDamm,     Dnit_leachDnit,     nit_leach
  real :: amm_immob0,    Damm_immobDamm,     Damm_immobDnit,     amm_immob    ! immob. to soil pools
  real :: nit_immob0,    Dnit_immobDamm,     Dnit_immobDnit,     nit_immob    ! immob. to soil pools
  real :: amm_immob_LS0, Damm_immob_LS_Damm, Damm_immob_LS_Dnit, amm_immob_LS ! immob. to slow litter
  real :: nit_immob_LS0, Dnit_immob_LS_Damm, Dnit_immob_LS_Dnit, nit_immob_LS ! immob. to slow litter
  real :: amm_immob_LF0, Damm_immob_LF_Damm, Damm_immob_LF_Dnit, amm_immob_LF ! immob. to fast litter
  real :: nit_immob_LF0, Dnit_immob_LF_Damm, Dnit_immob_LF_Dnit, nit_immob_LF ! immob. to fast litter
  real :: nitrif0,       DnitrifDamm,        DnitrifDnit,        nitrif
  real :: mnrl0,         DmnrlDamm,          DmnrlDnit,          mnrl
  
  real, parameter :: T = 0.9 ! fraction of nstore/nstore_opt above which active N uptake decreases linearly    
  
  ! parts of linear system for implicit time stepping
  real :: RHS_amm, RHS_nit ! right-hand side
  real :: a11, a12
  real :: a21, a22
  real :: det ! determinant of the system
  real :: delta_amm, delta_nit ! tendencies of ammonium and nitrate
  real :: res_amm, res_nit ! reiduals, for debugging
  real :: delta_n_conc ! tendency of the nitrogen concentration
  real :: drainage ! total drainage used in the nitrogen tendency calculations, kg/(m2 s)

  real :: amm_in,  nit_in ! inputs of mineral N, kg N/(m2 year)
  real :: k_plant  ! uptake rate factor
  real :: sum_root
  real :: k_transp ! rate of transpiration uptake, 1/year
  real :: k_drain  ! rate of drainage, 1/year
  real :: immob_k  ! combination of factors for immobilization rate calculations
  real :: immob_k_LS ! combination of factors for calculations of immobilization rate into slow litter
  real :: immob_k_LF ! combination of factors for calculations of immobilization rate into fast litter
  real :: leach_k  ! combination of factors for DOM leaching calculations

  real :: A        ! A-factor = reduction of decomposition rate due 
                   ! to soil moisture and temperature, unitless
  ! various tendencies (LF=fast litter, LS=slow litter, SS=slow soil, and 
  ! SP=passive soil)
  real :: decomp_C_LF, decomp_C_LS, decomp_C_SS, decomp_C_SP
  real :: decomp_N_LF, decomp_N_LS, decomp_N_SS, decomp_N_SP
  real :: immob_C_SS, immob_C_SP
  real :: immob_N_SS, immob_N_SP, immob_N
  real :: leach_C_LF, leach_C_LS, leach_C_SS
  real :: leach_N_LF, leach_N_LS, leach_N_SS
  
  ! updated amounts of soluble nitrogen
  real :: new_N_LF_sol, new_N_LS_sol, new_N_SS_sol

   ! Used in the denitrification function
  real :: tempfactor_denit, waterfactor_denit
  real :: denit_rate
  real :: denit_rate0, Ddenit_rateDnit

  integer :: sp ! shorthand for species
  integer :: i,j,k,face
  integer :: r_uptake_factor ! factor that controls N uptake by roots 
  real :: tot_n0, tot_n1, delta_n ! for conservation check
  real, parameter :: tolerance = 1e-14
  character(512) :: message
  integer :: y,mo,d,h,m,s ! components of date

  !+ conservation checking part 1
  tot_n0 = vegn%fast_litter_N + vegn%slow_litter_N &
         + vegn%slow_soil_N + vegn%res_soil_N &
         + vegn%amm + vegn%nitr
  !- conservation checking part 1

  if(is_watch_point())then
     write(*,*)'###### n_soil_int: input data #####'
     __DEBUG2__(soilt,theta)
     __DEBUG4__(ndep_nit,ndep_amm,ndep_don,nfert)
     __DEBUG3__(nfm_nit,nfm_amm,nfm_don)
     __DEBUG3__(transp,drain,surf_runoff)
     __DEBUG2__(vegn%amm, vegn%nitr)
     __DEBUG2__(vegn%fast_litter_N, vegn%fast_litter_C)
     __DEBUG2__(vegn%slow_litter_N, vegn%slow_litter_C)
     __DEBUG2__(vegn%slow_soil_N, vegn%slow_soil_C)
     __DEBUG2__(vegn%res_soil_N, vegn%res_soil_C)
     write(*,*)'###### n_soil_int: end of input data #####'
  endif
  
  if (use_surface_runoff) then
     drainage = drain + surf_runoff
  else
     drainage = drain
  endif
  
  ! calculate available concentrations of ammonium and nitrogen
  amm_conc = max(vegn%amm,0.0)/(soil_depth*amm_buff)
  nit_conc = max(vegn%nitr,0.0)/(soil_depth*nitr_buff)
  n_conc   = amm_conc+nit_conc

  if(is_watch_point())then
     __DEBUG3__(n_conc, amm_conc, nit_conc)
  endif
  
  ! calculate plant uptake rate
  sum_root = 0
  do i = 1,vegn%n_cohorts
     sp = vegn%cohorts(i)%species
     sum_root = sum_root + &
        kphalf*vegn%cohorts(i)%Pr*vegn%cohorts(i)%bliving*spdata(sp)%knup_plant
  enddo
  if (vegn%cohorts(1)%nstore > vegn%cohorts(1)%nstore_opt) then
    r_uptake_factor = 0
  else if (vegn%cohorts(1)%nstore < vegn%cohorts(1)%nstore_opt*T) then
    r_uptake_factor = 1
  else
    r_uptake_factor = -1/(1-T)*vegn%cohorts(1)%nstore/vegn%cohorts(1)%nstore_opt + 1/(1-T)
  endif
  if (sum_root > 0.0) then
    k_transp = seconds_per_year*MAX(transp/(soil_depth*dens_h2o),0.0)
    k_plant = sum_root*r_uptake_factor*mobility/(kphalf+mobility*N_conc) + k_transp 
    amm_uptake0 = k_plant*amm_conc
    Damm_uptakeDamm = ((kphalf+mobility*nit_conc)/(kphalf+mobility*N_conc)**2*mobility*sum_root &
    *r_uptake_factor +  k_transp)/(soil_depth*amm_buff)
    Damm_uptakeDnit = -mobility**2*amm_conc/(kphalf+mobility*N_conc)**2*sum_root* r_uptake_factor &
    /(soil_depth*nitr_buff)
    nit_uptake0 = k_plant*nit_conc
    Dnit_uptakeDnit  = ((kphalf+mobility*amm_conc)/(kphalf+mobility*N_conc)**2*mobility &
    * sum_root*r_uptake_factor +  k_transp)/(soil_depth*nitr_buff)
    Dnit_uptakeDamm  = -mobility**2*nit_conc/(kphalf+mobility*N_conc)**2*sum_root &
    * r_uptake_factor/(soil_depth*amm_buff)
  else 
    k_plant = 0
    amm_uptake0 = 0.; nit_uptake0 = 0;
    Damm_uptakeDamm = 0; Damm_uptakeDnit = 0;
    Dnit_uptakeDamm = 0; Dnit_uptakeDnit = 0;
  endif
  
  ! leaching rate
  k_drain = seconds_per_year*MAX(drainage,0.0)/dens_h2o
  amm_leach0 = k_drain * amm_conc
  nit_leach0 = k_drain * nit_conc
  Damm_leachDamm = k_drain/(soil_depth*amm_buff) ; Damm_leachDnit = 0.0
  Dnit_leachDnit = k_drain/(soil_depth*nitr_buff); Dnit_leachDamm = 0.0
  
  ! decomposition factor due to soil temperature and moisture
  A = A_function(soilt,theta)

  ! nitrification rate
  nitrif0     = A*nitrif_rate*mobility*vegn%amm
  DnitrifDamm = A*nitrif_rate*mobility ; DnitrifDnit = 0.0
  
  ! denitrification rate
  tempfactor_denit = Q10_denit**((soilt-Tr_denit)/10.0)
  
  if (theta<St_denit) then
     waterfactor_denit = 0
  else 
  waterfactor_denit=min(1.0,max(0.0,(((theta-st_denit)/(1-st_denit))**w_denit)))
  endif

  denit_rate0     = kd_denit*waterfactor_denit*tempfactor_denit*mobility*vegn%nitr
  Ddenit_rateDnit = kd_denit*waterfactor_denit*tempfactor_denit*mobility ; 
  
  amm_in = ndep_amm+nfm_amm
  nit_in = ndep_nit+nfm_nit+nfert
  ! update mass balance
  vegn%nitrogen_in = vegn%nitrogen_in + (ndep_nit + ndep_amm + ndep_don + nfm_nit + nfm_amm + nfm_don)*dt_fast_yr
  
  ! update fast_litter_N with DON input
  vegn%fast_litter_N = vegn%fast_litter_N + (ndep_don + nfm_don)*dt_fast_yr
  
  ! mineralization - immobilization
  ! [X.1] decomposition factor due to nitrogen abundance and its derivative
  n_factor0 = 1+min(nu_decomp*n_conc, 0.35)
  if (nu_decomp*n_conc<0.35) then
     Dn_factorDn = nu_decomp
  else
     Dn_factorDn = 0
  endif
  ! [X.2] fraction of stabilized (heavy) material and its derivative
  q_factor0 = (SOMinput_frac_slow*n_conc/(CtoN_SOM_slow*(eta_decomp+n_conc)) &
           +  SOMinput_frac_res/CtoN_SOM_res)
  Dq_factorDn = SOMinput_frac_slow * eta_decomp &
                                  / (CtoN_SOM_slow*(eta_decomp+n_conc)**2) 
  if (n_conc>0) then
     ! immobilization to soil pools
     immob_k = A*k_litter_slow*vegn%slow_litter_C*n_factor0*q_factor0/n_conc
     amm_immob0 = immob_k*amm_conc
     Damm_immobDamm = A*k_litter_slow*vegn%slow_litter_C/(soil_depth*amm_buff*n_conc)&
     *  ( &
        nit_conc/n_conc * n_factor0 * q_factor0 + &
        amm_conc * Dn_factorDn * q_factor0 + &
        amm_conc * n_factor0 * Dq_factorDn &
        )
     Damm_immobDnit = A*k_litter_slow*vegn%slow_litter_C/(soil_depth*nitr_buff*n_conc)&
     *  ( &
       -amm_conc/n_conc * n_factor0 * q_factor0 +  &
        amm_conc * Dn_factorDn * q_factor0 +      &
        amm_conc * n_factor0 * Dq_factorDn        &
        ) 
     nit_immob0 = immob_k*nit_conc
     Dnit_immobDnit = A*k_litter_slow*vegn%slow_litter_C/(soil_depth*nitr_buff*n_conc)&
     *  ( &
        amm_conc/n_conc * n_factor0 * q_factor0 + &
        nit_conc * Dn_factorDn * q_factor0 + &
        nit_conc * n_factor0 * Dq_factorDn &
        )
     Dnit_immobDamm = A*k_litter_slow*vegn%slow_litter_C/(soil_depth*amm_buff*n_conc)&
     *  ( &
       -nit_conc/n_conc * n_factor0 * q_factor0 + &
        nit_conc * Dn_factorDn * q_factor0 +     &
        nit_conc * n_factor0 * Dq_factorDn       &
        ) 
     ! immobilization to slow litter:
     if (vegn%slow_litter_C > vegn%slow_litter_N * CtoN_target) then
        immob_k_LS = A*k_litter_slow*vegn%slow_litter_C                  &
                *(1.0/CtoN_target-vegn%slow_litter_N/vegn%slow_litter_C) &
                / (kmhalf+n_conc)
     else
        immob_k_LS = 0
     endif

     amm_immob_LS0 = immob_k_LS * n_factor0 * amm_conc 
     Damm_immob_LS_Damm = immob_k_LS / (soil_depth*amm_buff) &
        * (amm_conc*dn_factorDn + n_factor0*(kmhalf+nit_conc)/(kmhalf+n_conc))
     Damm_immob_LS_Dnit = immob_k_LS * amm_conc / (soil_depth*amm_buff) &
        * (Dn_factorDn - n_factor0/(kmhalf+n_conc))

     nit_immob_LS0 = immob_k_LS * nit_conc 
     Dnit_immob_LS_Dnit = immob_k_LS / (soil_depth*nitr_buff) &
        * (nit_conc*dn_factorDn + n_factor0*(kmhalf+amm_conc)/(kmhalf+n_conc))
     Dnit_immob_LS_Damm = immob_k_LS * nit_conc / (soil_depth*nitr_buff) &
        * (Dn_factorDn - n_factor0/(kmhalf+n_conc))

     ! immobilization to fast litter
     if (vegn%fast_litter_C > vegn%fast_litter_N*CtoN_target) then
        immob_k_LF = A*k_litter_fast*vegn%fast_litter_C &
                *(1.0/CtoN_target-vegn%fast_litter_N/vegn%fast_litter_C) &
                /(kmhalf+n_conc)
     else
        immob_k_LF = 0
     endif
     
     amm_immob_LF0 = immob_k_LF * amm_conc
     Damm_immob_LF_Damm = immob_k_LF/(soil_depth*amm_buff) &
         * (kmhalf+nit_conc) / (kmhalf+n_conc)
     Damm_immob_LF_Dnit = - immob_k_LF/(soil_depth*amm_buff) &
         * amm_conc/(kmhalf+n_conc)

     nit_immob_LF0 = immob_k_LF * nit_conc
     Dnit_immob_LF_Dnit = immob_k_LF/(soil_depth*nitr_buff) &
         * (kmhalf+amm_conc) / (kmhalf+n_conc)
     Dnit_immob_LF_Damm = - immob_k_LF/(soil_depth*nitr_buff) &
         * nit_conc/(kmhalf+n_conc)
  else
    amm_immob0    = 0; Damm_immobDamm     = 0; Damm_immobDnit     = 0;
    nit_immob0    = 0; Dnit_immobDamm     = 0; Dnit_immobDnit     = 0;
    amm_immob_LS0 = 0; Damm_immob_LS_Damm = 0; Damm_immob_LS_Dnit = 0;
    nit_immob_LS0 = 0; Dnit_immob_LS_Damm = 0; Dnit_immob_LS_Dnit = 0;
    amm_immob_LF0 = 0; Damm_immob_LF_Damm = 0; Damm_immob_LF_Dnit = 0;
    nit_immob_LF0 = 0; Dnit_immob_LF_Damm = 0; Dnit_immob_LF_Dnit = 0;
    Dq_factorDn = 0; Dn_factorDn = 0;
  endif
  
  ! rates of mineral nitrogen input due to decomposition
  mnrl0 = A * &
    ( k_litter_fast*vegn%fast_litter_N          &
    + k_litter_slow*n_factor0*vegn%slow_litter_N &
    + k_SOM_slow*vegn%slow_soil_N               &
    + k_SOM_res*vegn%res_soil_N                 &
    )
  DmnrlDamm = A*k_litter_slow*Dn_factorDn*vegn%slow_litter_N/(soil_depth*amm_buff)
  DmnrlDnit = A*k_litter_slow*Dn_factorDn*vegn%slow_litter_N/(soil_depth*nitr_buff)
  
  ! right-hand-side of the linear system
  RHS_amm = amm_in - amm_immob0 - amm_immob_LS0 - amm_immob_LF0 &
          - amm_uptake0 - amm_leach0 - nitrif0 &
          + mnrl0
  RHS_nit = nit_in - nit_immob0 - nit_immob_LS0 - nit_immob_LF0 &
          - nit_uptake0 - nit_leach0 + nitrif0-denit_rate0
  
  ! implicit time step matrix
  ! a11*delta_amm + a12*delta_nit  =  RHS_amm
  ! a21*delta_amm + a22*delta_nit  =  RHS_nit
  a11 = 1.0/dt_fast_yr &
      + Damm_immobDamm + Damm_immob_LS_Damm + Damm_immob_LF_Damm &
      + Damm_uptakeDamm + Damm_leachDamm + DnitrifDamm - DmnrlDamm 
  a12 = Damm_immobDnit + Damm_immob_LS_Dnit + Damm_immob_LF_Dnit &
      + Damm_uptakeDnit + Damm_leachDnit + DnitrifDnit - DmnrlDnit 
  a21 = Dnit_immobDamm + Dnit_immob_LS_Damm + Dnit_immob_LF_Damm &
      + Dnit_uptakeDamm + Dnit_leachDamm - DnitrifDamm
  a22 = 1.0/dt_fast_yr &
      + Dnit_immobDnit + Dnit_immob_LS_Dnit + Dnit_immob_LF_Dnit &
      + Dnit_uptakeDnit + Dnit_leachDnit - DnitrifDnit+Ddenit_rateDnit
  ! determinant of the system
  det = a11*a22-a12*a21
  if (det==0) &
     call error_mesg('n_soil_int',&
                'determinanat of the implicit time step system is zero',FATAL)
  delta_amm = (RHS_amm*a22 - RHS_nit*a12)/det
  delta_nit = (a11*RHS_nit - RHS_amm*a21)/det
  delta_n_conc = delta_amm/(soil_depth*amm_buff)+delta_nit/(soil_depth*nitr_buff)

  if (is_watch_point()) then
     res_amm = a11*delta_amm+a12*delta_nit-RHS_amm
     res_nit = a21*delta_amm+a22*delta_nit-RHS_nit
     __DEBUG3__(a11,a12,RHS_amm)
     __DEBUG3__(a21,a22,RHS_nit)
     __DEBUG2__(delta_amm,res_amm)
     __DEBUG2__(delta_nit,res_nit)
  endif

  ! calculate updated nitrogen fluxes
  amm_uptake = amm_uptake0 + delta_amm*Damm_uptakeDamm + delta_nit*Damm_uptakeDnit
  nit_uptake = nit_uptake0 + delta_amm*Dnit_uptakeDamm + delta_nit*Dnit_uptakeDnit
  amm_leach  = amm_leach0  + delta_amm*Damm_leachDamm  + delta_nit*Damm_leachDnit
  nit_leach  = nit_leach0  + delta_amm*Dnit_leachDamm  + delta_nit*Dnit_leachDnit
  amm_immob  = amm_immob0  + delta_amm*Damm_immobDamm  + delta_nit*Damm_immobDnit
  amm_immob_LS = amm_immob_LS0 + delta_amm*Damm_immob_LS_Damm + delta_nit*Damm_immob_LS_Dnit
  amm_immob_LF = amm_immob_LF0 + delta_amm*Damm_immob_LF_Damm + delta_nit*Damm_immob_LF_Dnit
  nit_immob  = nit_immob0  + delta_amm*Dnit_immobDamm  + delta_nit*Dnit_immobDnit
  nit_immob_LS = nit_immob_LS0 + delta_amm*Dnit_immob_LS_Damm + delta_nit*Dnit_immob_LS_Dnit
  nit_immob_LF = nit_immob_LF0 + delta_amm*Dnit_immob_LF_Damm + delta_nit*Dnit_immob_LF_Dnit
  nitrif     = nitrif0     + delta_amm*DnitrifDamm     + delta_nit*DnitrifDnit
  mnrl       = mnrl0       + delta_amm*DmnrlDamm       + delta_nit*DmnrlDnit
  denit_rate = denit_rate0  + delta_nit*Ddenit_rateDnit
  
  ! calculate tendencies: note that updated values of q_factor and n_factor 
  ! are used, and the product is also linearized in immob_C_SS
  ! updated value of n_factor and q_factor
  n_factor   = n_factor0 + Dn_factorDn*delta_n_conc
  q_factor   = q_factor0 + Dq_factorDn*delta_n_conc
  ! tendencies of litter carbon and nitrogen due to decomposition
  decomp_C_LF = A*k_litter_fast*vegn%fast_litter_C
  decomp_N_LF = A*k_litter_fast*vegn%fast_litter_N
  decomp_N_LS = A*k_litter_slow*vegn%slow_litter_N*n_factor
  decomp_C_LS = A*k_litter_slow*vegn%slow_litter_C*n_factor
  ! accumulate for soil acceleration
  vegn%K_ls_eq = vegn%K_ls_eq &
               + A*k_litter_slow*n_factor*dt_fast_yr
  ! tendencies of slow soil carbon and nitrogen due to immobilization 
  immob_N = A*k_litter_slow*vegn%slow_litter_C* &
     ( n_factor0*q_factor0 &
     + Dn_factorDn*delta_n_conc*q_factor0 &
     + n_factor0*Dq_factorDn*delta_n_conc &
     )
  ! tendencies of slow soil carbon and nitrogen due to decomposition
  decomp_C_SS = A*k_SOM_slow*vegn%slow_soil_C
  decomp_N_SS = A*k_SOM_slow*vegn%slow_soil_N
  ! tendencies of passive soil carbon and nitrogen due to immobilization 
  ! and decomposition
  immob_C_SP = decomp_C_LS*SOMinput_frac_res 
  immob_N_SP = immob_C_SP/CtoN_SOM_res
  immob_N_SS = immob_N - immob_N_SP
  immob_C_SS = immob_N_SS*CtoN_SOM_slow
  decomp_C_SP = A*k_SOM_res*vegn%res_soil_C
  decomp_N_SP = A*k_SOM_res*vegn%res_soil_N
  ! tendencies due to dissolved organic matter leaching
  leach_k = seconds_per_year*MAX(drainage,0.0) / &
      (dens_h2o*soil_depth*DON_buff)
  leach_N_LF = leach_k*vegn%fast_litter_fsol*vegn%fast_litter_N
  leach_C_LF = leach_k*vegn%fast_litter_fsol*vegn%fast_litter_C
  leach_N_LS = leach_k*vegn%slow_litter_fsol*vegn%slow_litter_N
  leach_C_LS = leach_k*vegn%slow_litter_fsol*vegn%slow_litter_C
  leach_N_SS = leach_k*vegn%slow_soil_fsol*vegn%slow_soil_N
  leach_C_SS = leach_N_SS*CtoN_SOM_slow
  ! accumulate for soil acceleration
  vegn%K_ls_eq = vegn%K_ls_eq &
               + leach_k*vegn%slow_litter_fsol*dt_fast_yr
  vegn%K_ss_eq = vegn%K_ss_eq &
               + leach_k*vegn%slow_soil_fsol*dt_fast_yr
  ! update amount of soluble nitrogen
  ! why production terms are calculated so strangely? Is my understanding correct?
  new_N_LF_sol = vegn%fast_litter_N * vegn%fast_litter_fsol + (&
                 ! production is zero for fast litter soluble nitrogen
                 - leach_N_LF - decomp_N_LF*vegn%fast_litter_fsol )*dt_fast_yr
  new_N_LS_sol = vegn%slow_litter_N * vegn%slow_litter_fsol + (&
                 + decomp_N_LS*fraction_dissolve*theta &  ! production
                 - leach_N_LS - decomp_N_LS*vegn%fast_litter_fsol )*dt_fast_yr
  new_N_SS_sol = vegn%slow_soil_N * vegn%slow_soil_fsol + (&
                 + immob_N_SS*fraction_dissolve*theta &   ! production
                 - leach_N_SS - decomp_N_SS*vegn%fast_litter_fsol )*dt_fast_yr
  ! update mineral nitrogen
  vegn%amm  = vegn%amm  + delta_amm
  vegn%nitr = vegn%nitr + delta_nit
  ! update litter pools
  vegn%fast_litter_C = vegn%fast_litter_C - (decomp_C_LF+leach_C_LF)*dt_fast_yr
  vegn%fast_litter_N = vegn%fast_litter_N - (decomp_N_LF+leach_N_LF)*dt_fast_yr &
                     + (amm_immob_LF+nit_immob_LF)*dt_fast_yr
  vegn%slow_litter_C = vegn%slow_litter_C - (decomp_C_LS+leach_C_LS)*dt_fast_yr
  vegn%slow_litter_N = vegn%slow_litter_N - (decomp_N_LS+leach_N_LS)*dt_fast_yr &
                     + (amm_immob_LS+nit_immob_LS)*dt_fast_yr
  ! update soil pools
  vegn%slow_soil_C = vegn%slow_soil_C + (immob_C_SS-decomp_C_SS-leach_C_SS)*dt_fast_yr
  vegn%slow_soil_N = vegn%slow_soil_N + (immob_N_SS-decomp_N_SS-leach_N_SS)*dt_fast_yr
  vegn%res_soil_C  = vegn%res_soil_C  + (immob_C_SP-decomp_C_SP)*dt_fast_yr
  vegn%res_soil_N  = vegn%res_soil_N  + (immob_N_SP-decomp_N_SP)*dt_fast_yr
  ! ---- update patch information for soil acceleration: should be before updates
  ! of soluble fractions since slow_litter_fsol is used here
  vegn%A_eq       = vegn%A_eq       + A*dt_fast_yr
  vegn%K_ss_eq    = vegn%K_ss_eq    + A*k_SOM_slow*dt_fast_yr
  vegn%ssc_in_eq  = vegn%ssc_in_eq  + immob_C_SS*dt_fast_yr
  vegn%resc_in_eq = vegn%resc_in_eq + immob_C_SP*dt_fast_yr
  ! update soluble fractions
  vegn%fast_litter_fsol = 0
  if (vegn%fast_litter_N > 0) then
     vegn%fast_litter_fsol = min(1.0,max(0.0,new_N_LF_sol/vegn%fast_litter_N))
  endif
  vegn%slow_litter_fsol = 0
  if (vegn%slow_litter_N > 0) then
     vegn%slow_litter_fsol = min(1.0,max(0.0,new_N_LS_sol/vegn%slow_litter_N))
  endif
  vegn%slow_soil_fsol = 0
  if (vegn%slow_soil_N > 0) then
     vegn%slow_soil_fsol = min(1.0,max(0.0,new_N_SS_sol/vegn%slow_soil_N))
  endif

  ! calculate heterotrophic respiration
  vegn%rh = decomp_C_LF + decomp_C_LS + decomp_C_SS + decomp_C_SP &
          - immob_C_SS - immob_C_SP &
          ! for now, include leaching in the respiration; later 
          ! move it through the rivers
          + leach_C_LF + leach_C_LS + leach_C_SS 

  ! update plant nitrogen: note that uptake goes into the first cohort, which
  ! obviously would not work right if we had more than one cohort.
  vegn%cohorts(1)%ngain = vegn%cohorts(1)%ngain &
          + (amm_uptake+nit_uptake)*dt_fast_yr

  !update nitrogen budget
  vegn%nitrogen_out = vegn%nitrogen_out + (leach_N_LF + leach_N_LS + leach_N_SS + &
       amm_leach + nit_leach+denit_rate)*dt_fast_yr

  vegn%donleach = (leach_N_LF+leach_N_LS+leach_N_SS)/seconds_per_year
  vegn%ammleach = amm_leach/seconds_per_year
  vegn%nitrleach = nit_leach/seconds_per_year

  ! ---- diagnostic section
  call send_tile_data(id_sorg_depth, soil_depth, diag)
  call send_tile_data(id_sorg_T, soilt, diag)
  call send_tile_data(id_sorg_theta, theta, diag)
  call send_tile_data(id_sorg_uptk, transp, diag)
  call send_tile_data(id_sorg_drain, drainage, diag)
  call send_tile_data(id_sorg_A, A, diag)
  if (.not. (abs(vegn%amm)<1e38)) then 
  call get_current_point(i,j,k)
  write (*,*) 'vegn%amm is bad', i,j,k,vegn%amm
  endif  
  call send_tile_data(id_n_amm, vegn%amm, diag)
  if (.not. (abs(vegn%nitr)<1e38)) then 
  call get_current_point(i,j,k)
  write (*,*) 'vegn%nitr is bad', i,j,k,vegn%nitr
  endif
  call send_tile_data(id_n_nit, vegn%nitr, diag)
  amm_conc = max(vegn%amm,0.0)/(soil_depth*amm_buff)
  nit_conc = max(vegn%nitr,0.0)/(soil_depth*nitr_buff)
  n_conc   = amm_conc+nit_conc
  call send_tile_data(id_n_conc, n_conc, diag)
  call send_tile_data(id_netmin, &
    mnrl-amm_immob-nit_immob-amm_immob_LS-nit_immob_LS-amm_immob_LF-nit_immob_LF, diag)
  call send_tile_data(id_immob, amm_immob+nit_immob+amm_immob_LS+nit_immob_LS+amm_immob_LF+nit_immob_LF, diag)
  call send_tile_data(id_immob_SS, immob_N_SS, diag)
  call send_tile_data(id_immob_SP, immob_N_SP, diag)
  call send_tile_data(id_immob_LS, amm_immob_LS+nit_immob_LS, diag)
  call send_tile_data(id_immob_LF, amm_immob_LF+nit_immob_LF, diag)
  call send_tile_data(id_decomp_C_SS, decomp_C_SS, diag)
  call send_tile_data(id_decomp_C_SP, decomp_C_SP, diag)
  call send_tile_data(id_decomp_C_LS, decomp_C_LS, diag)
  call send_tile_data(id_decomp_C_LF, decomp_C_LF, diag)
  call send_tile_data(id_nitrif, nitrif, diag)
  call send_tile_data(id_ammleach, amm_leach, diag)
  call send_tile_data(id_nitrleach, nit_leach, diag)
  if (.not. (abs(amm_leach)<1e38)) then 
  call get_current_point(i,j,k)
  write (*,*) 'amm_leach is bad', i,j,k,amm_leach
  endif
  if (.not. (abs(nit_leach)<1e38)) then 
  call get_current_point(i,j,k)
  write (*,*) 'nit_leach is bad', i,j,k,nit_leach
  endif
  call send_tile_data(id_nleach, amm_leach+nit_leach, diag)
  call send_tile_data(id_donleach, leach_N_LF+leach_N_LS+leach_N_SS, diag)
  call send_tile_data(id_nuprate, amm_uptake+nit_uptake, diag)
  call send_tile_data(id_ndep, ndep_nit+ndep_amm+ndep_don, diag)
  call send_tile_data(id_ndep_nit, ndep_nit, diag)
  call send_tile_data(id_ndep_amm, ndep_amm, diag)
  call send_tile_data(id_ndep_don, ndep_don, diag)
  call send_tile_data(id_nfert, nfert, diag)
  call send_tile_data(id_nfm_nit, nfm_nit, diag)
  call send_tile_data(id_nfm_amm, nfm_amm, diag)
  call send_tile_data(id_nfm_don, nfm_don, diag)
  call send_tile_data(id_nitrogen_in,vegn%nitrogen_in, diag)
  call send_tile_data(id_nitrogen_out,vegn%nitrogen_out, diag)
  if (.not. (abs(denit_rate)<1e38)) then 
  call get_current_point(i,j,k)
  write (*,*) 'denit_rate is bad', i,j,k,denit_rate
  endif
  call send_tile_data(id_denit_rate, denit_rate, diag)

  !+ conservation checking part 2
  tot_n1 = vegn%fast_litter_N + vegn%slow_litter_N &
         + vegn%slow_soil_N + vegn%res_soil_N &
         + vegn%amm + vegn%nitr
  delta_n=( ndep_nit + ndep_amm + ndep_don + nfert + nfm_nit + nfm_amm + nfm_don &
           -leach_N_LF-leach_N_LS-leach_N_SS &
           -amm_uptake-nit_uptake &
           -amm_leach-nit_leach &
          -denit_rate)*dt_fast_yr
  !- conservation checking part 2
  if(is_watch_point())then
     write(*,*)'###### n_soil_int: result #####'
     __DEBUG2__(vegn%amm, vegn%nitr)
     __DEBUG2__(vegn%fast_litter_N, vegn%fast_litter_C)
     __DEBUG2__(vegn%slow_litter_N, vegn%slow_litter_C)
     __DEBUG2__(vegn%slow_soil_N, vegn%slow_soil_C)
     __DEBUG2__(vegn%res_soil_N, vegn%res_soil_C)
     __DEBUG1__(vegn%rh)
     __DEBUG4__(decomp_C_LF,decomp_C_LS,decomp_C_SS,decomp_C_SP)
     __DEBUG2__(immob_C_SS,immob_C_SP)

     ! mineralization check: the difference below is an indicator of error
     ! in mineralization calculations
     __DEBUG4__(decomp_N_LF,decomp_N_LS,decomp_N_SS,decomp_N_SP)
     __DEBUG2__(decomp_N_LF+decomp_N_LS+decomp_N_SS+decomp_N_SP-mnrl,mnrl)

     ! immobilization check: the difference value below is an indicator of error 
     ! in immobilization calculations
     __DEBUG3__(immob_N_SS,immob_N_SP,immob_N_SS+immob_N_SP)
     __DEBUG3__(amm_immob,nit_immob,amm_immob+nit_immob)
     __DEBUG1__(amm_immob+nit_immob-immob_N_SS-immob_N_SP)

     __DEBUG3__(leach_C_LF,leach_C_LS,leach_C_SS)
     __DEBUG3__(leach_N_LF,leach_N_LS,leach_N_SS)
     __DEBUG2__(amm_uptake,nit_uptake)
     __DEBUG2__(amm_leach,nit_leach)
     __DEBUG1__(denit_rate)
     __DEBUG3__(tot_n0,tot_n1,delta_n)
     __DEBUG1__(tot_n1-tot_n0-delta_n)
     write(*,*)'###### n_soil_int: end of result #####'
  endif
  !+ conservation checking part 3
  !if(abs(tot_n1-tot_n0-delta_n)>tolerance) then
  !   call get_current_point(i,j,k,face)
  !   call get_date(lnd%time,y,mo,d,h,m,s)
  !   write(message,'(a,4(x,a,i4),4(x,a,g),x,a,i4.4,2("-",i2.2),x,i2.2,2(":",i2.2))')&
  !        'non-conservation of nitrogen', &
  !        'at i=',i,'j=',j,'k=',k,'face=',face, &
  !        'tot_n0',tot_n0,'tot_n1',tot_n1,'delta_n=',delta_n,&
  !        'difference=',abs(tot_n1-tot_n0-delta_n),&
  !        'time=',y,mo,d,h,m,s
  !   call error_mesg('n_soil_int',message,FATAL)
  !   
  !endif
  !- conservation checking part 3

end subroutine


! ============================================================================
! given carbon litter with specified C-to-N ratio, distributes it between slow 
! and fast litter pools, and updates litter nitrogen pools respectively
subroutine litterfall(vegn, litter_C, CtoN)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: litter_C ! amount of litter carbon
  real, intent(in) :: CtoN     ! C to N stoichiometric ratio in litter

  ! --- local vars
  real :: fast_C, slow_C ! contributions to respective litter carbon pools
  real :: fast_N, slow_N ! contributions to respective litter nitrogen pools
  real :: CtoN_slow ! C-to-N ration in slow litter
  real :: N_sol ! mass of soluble slow nitrogen

  !  divide carbon litter-fall into fast and slow portions
  fast_C = litter_C*MAX(a_lf - b_lf*lignin_lf*CtoN, f_lf_min)
  slow_C = litter_C - fast_C

  ! update carbon litter pools and slow litter pool lignin fraction
  vegn%fast_litter_C = vegn%fast_litter_C + fast_C
  vegn%slow_litter_C = vegn%slow_litter_C + slow_C
  
  ! carbon budget tracking
  vegn%flc_in    = vegn%flc_in    + fast_C
  vegn%slc_in    = vegn%slc_in    + slow_C
  vegn%slc_in_eq = vegn%slc_in_eq + slow_C
  
  ! update nitrogen litter pools
  CtoN_slow = MAX(CtoN, CtoN_litter_min) ! litter into slow has c:n ratio = CtoN_slow;
  fast_N = 0; slow_N = 0
  if (fast_C > 0 .and. slow_C > 0) then
    slow_N = slow_C/CtoN_slow
    fast_N = (fast_C+slow_C)/CtoN - slow_N
  else if (slow_C > 0.0) then
    slow_N = slow_C/CtoN
  else if (fast_C > 0.0) then
    fast_N = fast_C/CtoN
  endif
  
  ! update the litter pools of nitrogen, preserving the mass of nitrogen in
  ! the soluble fraction
  vegn%fast_litter_N = vegn%fast_litter_N + fast_N
  N_sol = vegn%slow_litter_N*vegn%slow_litter_fsol;
  vegn%slow_litter_N = vegn%slow_litter_N + slow_N
  if (slow_N > 0) then
    vegn%slow_litter_fsol = N_sol/vegn%slow_litter_N
  else
    vegn%slow_litter_fsol = 0.0
  endif
  
  ! nitrogen budget tracking
  vegn%fln_in    = vegn%fln_in    + fast_N;
  vegn%sln_in    = vegn%sln_in    + slow_N;
  vegn%sln_in_eq = vegn%sln_in_eq + slow_N;
   
end subroutine litterfall


! ============================================================================
! plant nitrogen turnover on fast time step
subroutine n_md (vegn, cc, diag, woodgain, md_leaf, CtoNleaf, md_root, CtoNroot, md_wood)
  type(vegn_tile_type), intent(inout) :: vegn
  type(diag_buff_type), intent(inout) :: diag
  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(in) :: woodgain
  real, intent(in) :: md_leaf
  real, intent(in) :: CtoNleaf
  real, intent(in) :: md_root
  real, intent(in) :: CtoNroot
  real, intent(in) :: md_wood

  ! local vars
  real :: fixcost   
  
  cc%ngain = cc%ngain - (md_leaf/CtoNleaf + md_root/CtoNroot + &
                         woodgain/spdata(cc%species)%CtoN(CMPT_WOOD));
  cc%nwoodgain = cc%nwoodgain &
               + woodgain/spdata(cc%species)%CtoN(CMPT_WOOD) &
               - md_wood/getCtoNwood(cc)
  
  ! budget
  vegn%vegN_out = vegn%vegN_out + &
    md_leaf/CtoNleaf + md_root/CtoNroot + md_wood/getCtoNwood(cc);
  
  select case (n_fixrate_option)
  case (COMP_NFIX_PRESCRIBED)
    ! 'prescribed'
    cc%nfixrate = nfixrate;
  case (COMP_NFIX_1) ! dynamic N fixation
    cc%nfixrate = cc%Pl*cc%bliving*cc%nfix_per_bliv;
  case (COMP_NFIX_2)
    ! Do nothing in this case, since cc%nfixrate is already calculated in
    ! n_fixation. It can't be computed here since it requires nstore_opt
  case (COMP_NFIX_3)
    cc%nfixrate = 1.8e-3*(1-exp(-cc%gpp*0.9));
  end select
  ! nitrogen fixation and carbon cost
  ! further N fixation rates can be obtained from input 
  ! fields (constant, but geographical pattern)
  ! or from namelist (globally constant)

  !  carbon costs associated with fixation
  fixcost = cc%nfixrate*nfixcost*dt_fast_yr;
  cc%carbon_gain = cc%carbon_gain - fixcost
  cc%carbon_loss = cc%carbon_loss + fixcost ! used in diagnostics only
  vegn%fast_litter_C = vegn%fast_litter_C + fixcost
  vegn%veg_out = vegn%veg_out + fixcost;
  vegn%flc_in  = vegn%flc_in + fixcost;
  cc%Ngain = cc%Ngain + cc%nfixrate*dt_fast_yr;
  vegn%vegn_in = vegn%vegn_in + cc%nfixrate*dt_fast_yr;
  vegn%nitrogen_in = vegn%nitrogen_in + cc%nfixrate*dt_fast_yr;
end subroutine n_md



! ===========================================================================
! This function is (ab)used to accommodate carbon and nitrogen in plants
! Processes that are here:
! a) defining nitrogen storage capacity 
! b) set thresholds for plant nitrogen uptake/shutdown of uptake (currently 
!    worth 1 year of turnover)
! c) defining down-regulation of photosynthesis due to N limitation
! d) root carbon exudation, in severe N limitation
! e) nitrogen fixation
! f) root nitrogen exudation if N is big
subroutine n_reallocate(vegn,cc)
  type(vegn_tile_type),   intent(inout) :: vegn
  type(vegn_cohort_type), intent(inout) :: cc

  real, parameter :: excrete_thresh = 1.7;
  ! ---- local vars
  real :: ndemand; ! nitrogen demand for conventional allocation
  real :: nsupply; ! current nitrogen in plant
  real :: excess;  ! excess nitrogen to be stored
  real :: excrete
  integer :: sp ! shorthand for current species
   
  if (is_watch_point()) then
     write(*,*)'############## n_reallocate input #################'
     __DEBUG2__(cc%species,cc%bliving)
     __DEBUG2__(cc%nliving,cc%nstore)
     __DEBUG3__(cc%Pl,cc%Pr,cc%Pr)
     __DEBUG1__(cc%Psw_alphasw)
     write(*,*)'########### end of n_reallocate input #############'
  endif

  sp = cc%species
  ! calculate current nitrogen in plant
  nsupply = cc%nliving + cc%nstore;
  
  if (cc%bliving > 0) then
    ! calculate nitrogen needed to satisfy prescribed tissue cc:N ratios 
    ndemand = cc%bliving * (                  &
         cc%Pl/spdata(sp)%CtoN(CMPT_LEAF)     &
       + cc%Psw/spdata(sp)%CtoN(CMPT_SAPWOOD) & 
       + cc%Pr/spdata(sp)%CtoN(CMPT_ROOT)     &
       )
    excess = nsupply - ndemand;
    ! storage is an excess of current nitrogen over demand
    cc%nstore  = excess;
    cc%nliving = ndemand;
    
    ! calculate optimum nitrogen storage
    cc%nstore_opt = n_store_optimum(cc);

    ! in grasses to avoid ridiculous nitrogen concentration:
    ! dump nitrogen into ground that is higher than a threshold
    if ( (sp == SP_C4GRASS.or.sp==SP_C3GRASS).and. &
          cc%nstore > excrete_thresh*cc%nstore_opt    ) then
      excrete = cc%nstore - excrete_thresh*cc%nstore_opt;
      call n_excrete(vegn,excrete);
      cc%nstore = excrete_thresh*cc%nstore_opt;
      cc%nsf = excrete_thresh;
      vegn%vegn_out = vegn%vegn_out + excrete;
    endif
    
    ! update plant nitrogen status
    cc%nsf = max(cc%nstore,0.0)/cc%nstore_opt;
  else ! bliving == 0; dump all nitrogen and reset nsf
    ! dump nliving
    excrete = max(cc%nliving,0.0);
    call n_excrete(vegn,excrete);
    vegn%vegn_out = vegn%vegn_out + excrete;
    cc%nliving    = cc%nliving    - excrete;
    ! and nstore
    excrete = max(cc%nstore,0.0);
    call n_excrete(vegn,excrete);
    vegn%vegn_out = vegn%vegn_out + excrete;
    cc%nstore     = cc%nstore     - excrete;

    cc%nsf = 0.;
  endif
  
  if (is_watch_point()) then
     __DEBUG3__(ndemand,nsupply,excess)
     __DEBUG1__(cc%nstore_opt)
     __DEBUG3__(cc%nstore,cc%nliving,cc%nsf)
     write(*,*)'########### end of n_reallocate #############'
  endif
end subroutine n_reallocate


! ===========================================================================
! dump some (typically residual) amount of nitrogen
subroutine n_excrete(vegn, ndump)
   type(vegn_tile_type), intent(inout) :: vegn
   real, intent(in) :: ndump
   
   character(256) :: message
   integer :: i,j,k,face

  if(ndump<0) then 
     call get_current_point(i,j,k,face)
     write(message,'(a,4(x,a,i4), x,a,g)')&
          'trying to dump negative amount of nitrogen', &
          'at i=',i,'j=',j,'k=',k,'face=',face, &
          'dumped amount=', ndump
     call error_mesg('n_excrete',message,FATAL)
  endif
  
  vegn%fast_litter_N = vegn%fast_litter_N + ndump;
  vegn%fln_in = vegn%fln_in + ndump;
end subroutine


! ===========================================================================
! given increments of nitrogen per tissue, depletes the tissue nitrogen pools
! and calculates concentration factor due to ash nitrogen enrichment
subroutine n_fire(vegn, cc, dnleaf, dnroot, dnsw, dnwood, smoke_fraction, &
	          conc_factor)
  type(vegn_tile_type)  , intent(inout) :: vegn
  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(in) :: dnleaf, dnroot, dnsw, dnwood, smoke_fraction
  real, intent(out):: conc_factor
  
  real :: dnsw_ !copy of dnsw, likely updated 
  real :: dnstore

  if (cc%status == LEAF_OFF) then
    cc%nlv = cc%nlv - (dnleaf + dnroot)
  else
    cc%nl = cc%nl - dnleaf
    cc%nr = cc%nr - dnroot
  endif

  ! calculate the loss fraction of store and update N loss from sap
  dnstore    = 0.    
  if (cc%nwood > 0.) dnstore = dnwood/cc%nwood*MAX(cc%nstore,0.)
  dnsw_      = dnsw - dnstore
  cc%nstore  = cc%nstore - dnstore
  cc%nsw     = cc%nsw - dnsw_
  cc%nwood   = cc%nwood - dnwood
  cc%nliving = cc%nliving - (dnleaf + dnroot + dnsw_)

  ! factor that changes C:N in litter fluxes, as nitrogen is enriched in ash
  conc_factor = (1-smoke_fraction)/                &    ! carbon over
                (1-smoke_fraction*(1-ash_fraction));    ! nitrogen

  vegn%csmoke_pool(isub_N) = vegn%csmoke_pool(isub_N) + &
     (dnleaf + dnroot + dnsw + dnwood)*smoke_fraction*(1-ash_fraction)

  vegn%vegn_out = vegn%vegn_out + dnleaf + dnroot + dnsw + dnwood;
end subroutine n_fire


! ===========================================================================
! given bliving, Pl, npp_prev_day, npp_prev_day_pot, etc...
! updates (prognostically) nfix_per_bliv
subroutine n_fixation(vegn, cc, nsf_old)
  type(vegn_tile_type), intent(inout) :: vegn
  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(in) :: nsf_old

  ! ---- local vars
  real :: nstore_opt;
  real :: fix_demand = 0.;
  
  real :: dt_fix_yr = 1.0/365.0; ! time step for fixation, years
  real :: deficit;
  real :: weight;
  real :: h1;
  real :: fixcost;
  real :: freespace;
  real :: CtoN_plant ! optimal (target) plant C:N ratio
!  real :: epsilon;

  if (cc%nsf > 0) then
    nstore_opt = max(cc%nstore,0.0)/cc%nsf;
  else
    nstore_opt = cc%bliving/35. ! assume some cton ratio
  endif

  select case (n_fixrate_option)
  case (COMP_NFIX_1)
    ! 1. TYPE, mechanistic n_fixation using N demand and LAI (light)

    ! calculate accumulated nitrogen deficit
    deficit = 0
    if (cc%bliving > 0) then
      CtoN_plant = (cc%bliving)/(cc%nliving + nstore_opt);
      deficit = (nsf_old-cc%nsf)*nstore_opt/dt_fix_yr &
	      + (cc%npp_previous_day_pot-cc%npp_previous_day)/CtoN_plant;
    endif
    
#if 1

    ! the need for n fixation scales with how photosynthesis is punished
    if (deficit + cc%nfix_per_bliv*cc%bliving*cc%Pl > 0) then
      weight = 1.2*k_phot*(1-n_phot(cc%nsf)) ! 20% higher demand of fixer
    else
      weight = 1.0;
    endif

    if (cc%nsf < 1.0 .and. cc%lai > 1e-12) then
      fix_demand = weight*(deficit/(cc%bliving*cc%Pl) + cc%nfix_per_bliv);
    else
      fix_demand = 0
    endif

    fix_demand = MAX(fix_demand,0.0)

    ! increase in fixation requires open space, except tropics
    if (cc%lai < 1.e-12 .or. cc%species==SP_TROPICAL) then
      freespace = 1.0
    else 
      freespace = exp(-extinct_coeff*cc%lai);
    endif
    
    ! constrain maximum fixation to entire NPP (which is a lot)
    if (cc%bliving > 0) &
      fix_demand = MIN(fix_demand,cc%npp_previous_day_pot/(nfixcost*cc%bliving*cc%Pl));
    
    cc%nfix_per_bliv = cc%nfix_per_bliv + &
      (fix_demand*freespace - cc%nfix_per_bliv)/spdata(cc%species)%tau_fix*dt_fix_yr

    ! avoid huge and unreasonable numbers in nitrogen fixation
    cc%nfix_per_bliv = MIN(cc%nfix_per_bliv ,1000.);
    cc%nfix_per_bliv = MAX(cc%nfix_per_bliv ,0.);

#else 
    ! is this branch needed?
    
    if (cc%lai > 1e-18 .and. cc%bliving>0) then

      if (nfixcost > 0.) then fixcost = nfixcost; else fixcost = 9.4; endif
      ! 1. calculate deficit that comes from previous day:
      if (deficit > 0.) fix_demand = fix_demand+cc%nfixrate; else fix_demand = 0.;
      ! 2. restore storage nitrogen over the planning horizon time scale:
      fix_demand = fix_demand + (1.-MIN(cc%nsf,1.))*nstore_opt/ndata->planning_horizon;
      ! 3. allow for 20 % more nitrogen richness in fixers
      fix_demand = fix_demand*1.2;
      ! 4. do not allow negative spiral of N limitation -> reducing leafmass -> reducing fixation
      ! epsilon    = MIN(cc%carbon_gain - cc%carbon_loss,0.)*cc%Pl;
      ! fix_demand = MAX(fix_demand,0.)*cc%bl/(cc%bl+epsilon); ! per unit leaf

      ! 5. fixation does pay off only if costs also contribute to increase in npp
      fix_demand = MIN(fix_demand,(cc%npp_previous_day_pot-cc%npp_previous_day)/fixcost);
      fix_demand = fix_demand/(cc%bliving);

      if (cc%species == 3) fpc = 0.; else fpc = 1.-exp(-extinct_coef*cc%lai);
      h1     = 1. - fpc*MIN(fix_demand,cc%nfix_per_bliv)/cc%nfix_per_bliv;
      kfixdt = h1/(spdata(cc%species)%tau_fix/dt_fix_yr)
      
      cc%nfix_per_bliv = cc%nfix_per_bliv + (fix_demand - cc%nfix_per_bliv)*(1.-exp(-kfixdt));

      ! avoid huge unreasonable numbers in nitrogen fixation
      cc%nfix_per_bliv  = MIN(cc%nfix_per_bliv ,1000.);
      cc%nfix_per_bliv  = MAX(cc%nfix_per_bliv ,0.);

    endif
#endif
   
    ! prevent that small plants perish because of N limitation
    if (cc%bliving < 1e-3 .and. cc%nsf < 1 .and. cc%bliving > 0) then
      cc%nfix_per_bliv = weight*MAX(deficit,0.)/(cc%bliving*cc%Pl)
      cc%nfixrate = 0.0;
    endif

  case (COMP_NFIX_2)
    ! 3rt TYPE: balance budget with instantaneous fixation but update 
    ! nfix_per_bliv with a rate constant, which helps to initialize actual 
    ! fixation.
    ! This formulation is thought for carbon only runs determines fixation rates 
    ! when N status gets very low and adds them up continuously in the variable 
    ! nfix_per_bliv.
    ! THEREFORE nfix_per_bliv HAS A DIFFERENT MEANING IN comp_nfix1 AND
    ! comp_nfix2: HERE nfix_per_bliv CAN THEN BE USED TO ESTIMATE FIXATION
    ! RATES W/O ADDING TOO MUCH SHOCK WHEN INTRODUCING N LIMITATION
    cc%nfixrate      = nstore_opt*max(0.2-cc%nsf,0.0)/dt_fix_yr
    if (is_watch_point()) then
      write(*,*)'############ comp_nfix_2 ################'
      __DEBUG2__(cc%nfix_per_bliv,fix_demand)
    endif
    if (cc%bliving > 0.) then
      cc%nfix_per_bliv = cc%nfix_per_bliv + fix_demand;
    endif
    if (is_watch_point()) then
      __DEBUG2__(nstore_opt,cc%nsf)
      __DEBUG2__(cc%nfix_per_bliv,fix_demand)
      write(*,*)'############ end of comp_nfix_2 ################'
    endif

  case (COMP_NFIX_3,COMP_NFIX_PRESCRIBED)
    ! directly into n_md module:
    ! cc%nfixrate = 1.8e-3*(1-exp(-cc%npp_previous_day*3.))
  end select
end subroutine n_fixation

end module
