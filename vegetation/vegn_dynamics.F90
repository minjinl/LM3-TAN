! ============================================================================
! updates carbon pools and rates on the fast time scale
! ============================================================================
module vegn_dynamics_mod

#include "../shared/debug.inc"

use fms_mod, only: write_version_number
use time_manager_mod, only: time_type

use land_constants_mod, only : seconds_per_year, mol_C
use land_tile_diag_mod, only : &
     register_tiled_diag_field, send_tile_data, diag_buff_type
use land_substances_mod, only : n_substances
use vegn_data_mod, only : spdata, &
     CMPT_VLEAF, CMPT_SAPWOOD, CMPT_ROOT, CMPT_WOOD, CMPT_LEAF, LEAF_ON, LEAF_OFF, &
     fsc_liv, fsc_wood, K1, K2, soil_carbon_depth_scale, C2B, agf_bs, &
     l_fract, a_lf, b_lf, f_lf_min, lignin_lf, k_litter_fast, k_litter_slow, &
     SOMinput_frac_slow, SOMinput_frac_res, k_SOM_slow, k_SOM_res, isub_C
use nitrogen_data_mod, only: do_nitrogen, do_nitrogen_plant_feedback, isub_N
use nitrogen_mod,  only: &
     getCtoNleaf, getCtoNroot, getCtoNwood, litterfall, n_phot, n_md, &
     n_soil_int, n_fixation, n_excrete, n_reallocate
use vegn_tile_mod, only: vegn_tile_type, update_biomass_pools, n_update_pool
use vegn_cohort_mod, only : vegn_cohort_type, A_function, &
     update_bio_living_fraction, update_species
     

use land_debug_mod, only : is_watch_point

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_dynamics_init

public :: vegn_carbon_int   ! fast time-scale integrator of vegetation carbon balance
public :: soil_carbon_int   ! fast time-scale integrator of soil carbon balance
public :: vegn_growth       ! slow time-scale re-distributor of accumulated carbon
public :: vegn_daily_npp    ! updates values of daily-average npp
public :: vegn_phenology    !
public :: vegn_biogeography !
public :: update_soil_pools !
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), private, parameter :: &
   version = '$Id: vegn_dynamics.F90,v 17.0.2.11 2010/08/23 20:10:04 slm Exp $', &
   tagname = '$Name: nitro_20111003_slm $' ,&
   module_name = 'vegn'
real, parameter :: GROWTH_RESP=0.333  ! fraction of npp lost as growth respiration


! ==== module data ===========================================================
real    :: dt_fast_yr ! fast (physical) time step, yr (year is defined as 365 days)

! diagnostic field IDs
integer :: id_npp, id_nep, id_gpp, id_npp_pot, id_gpp_pot
integer :: id_fast_soil_C, id_slow_soil_C, id_res_soil_C
integer :: id_fast_litter_C, id_slow_litter_C
integer :: id_rsoil, id_rsoil_fast
integer :: id_resp, id_resl, id_resr, id_resg, id_asoil
integer :: id_soilt, id_theta, id_litter
integer :: id_nfixrate


contains

! ============================================================================
subroutine vegn_dynamics_init(id_lon, id_lat, time, delta_time)
  integer        , intent(in) :: id_lon ! ID of land longitude (X) axis 
  integer        , intent(in) :: id_lat ! ID of land latitude (Y) axis
  type(time_type), intent(in) :: time       ! initial time for diagnostic fields
  real           , intent(in) :: delta_time ! fast time step, s

  call write_version_number(version, tagname)

  ! set up global variables
  dt_fast_yr = delta_time/seconds_per_year

  ! register diagnostic fields
  id_gpp = register_tiled_diag_field ( module_name, 'gpp',  &
       (/id_lon,id_lat/), time, 'gross primary production', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_npp = register_tiled_diag_field ( module_name, 'npp',  &
       (/id_lon,id_lat/), time, 'net primary production', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_nep = register_tiled_diag_field ( module_name, 'nep',  &
       (/id_lon,id_lat/), time, 'net ecosystem production', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_gpp_pot = register_tiled_diag_field ( module_name, 'gpp_pot',  &
       (/id_lon,id_lat/), time, 'non-nitrogen-limited gross primary production', &
       'kg C/(m2 year)', missing_value=-100.0 )
  id_npp_pot = register_tiled_diag_field ( module_name, 'npp_pot',  &
       (/id_lon,id_lat/), time, 'non-nitrogen-limited net primary production', &
       'kg C/(m2 year)', missing_value=-100.0 )
  id_litter = register_tiled_diag_field (module_name, 'litter', (/id_lon,id_lat/), &
       time, 'litter productivity', 'kg C/(m2 year)', missing_value=-100.0)
  id_fast_litter_C = register_tiled_diag_field ( module_name, 'flc',  &
       (/id_lon,id_lat/), time, 'fast litter carbon', 'kg C/m2', &
       missing_value=-100.0 )
  id_slow_litter_C = register_tiled_diag_field ( module_name, 'slc',  &
       (/id_lon,id_lat/), time, 'slow litter carbon', 'kg C/m2', &
       missing_value=-100.0 )
  id_fast_soil_C = register_tiled_diag_field ( module_name, 'fsc',  &
       (/id_lon,id_lat/), time, 'fast soil carbon', 'kg C/m2', &
       missing_value=-100.0 )
  id_slow_soil_C = register_tiled_diag_field ( module_name, 'ssc',  &
       (/id_lon,id_lat/), time, 'slow soil carbon', 'kg C/m2', &
       missing_value=-100.0 )
  id_res_soil_C = register_tiled_diag_field ( module_name, 'rsc',  &
       (/id_lon,id_lat/), time, 'resilient (passive) soil carbon', 'kg C/m2', &
       missing_value=-100.0 )
  id_resp = register_tiled_diag_field ( module_name, 'resp', (/id_lon,id_lat/), &
       time, 'respiration', 'kg C/(m2 year)', missing_value=-100.0 )
  id_resl = register_tiled_diag_field ( module_name, 'resl', (/id_lon,id_lat/), &
       time, 'leaf respiration', 'kg C/(m2 year)', missing_value=-100.0 )
  id_resr = register_tiled_diag_field ( module_name, 'resr', (/id_lon,id_lat/), &
       time, 'root respiration', 'kg C/(m2 year)', missing_value=-100.0 )
  id_resg = register_tiled_diag_field ( module_name, 'resg', (/id_lon,id_lat/), &
       time, 'growth respiration', 'kg C/(m2 year)', missing_value=-100.0 )
  id_rsoil = register_tiled_diag_field ( module_name, 'rsoil',  &
       (/id_lon,id_lat/), time, 'soil respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_rsoil_fast = register_tiled_diag_field ( module_name, 'rsoil_fast',  &
       (/id_lon,id_lat/), time, 'fast soil carbon respiration', 'kg C/(m2 year)', &
       missing_value=-100.0 )
  id_asoil = register_tiled_diag_field ( module_name, 'asoil',  &
       (/id_lon,id_lat/), time, 'aerobic activity modifier', &
       missing_value=-100.0 )
  id_soilt = register_tiled_diag_field ( module_name, 'tsoil_av',  &
       (/id_lon,id_lat/), time, 'average soil temperature for carbon decomposition', 'degK', &
       missing_value=-100.0 )
  id_theta = register_tiled_diag_field ( module_name, 'theta',  &
       (/id_lon,id_lat/), time, 'average soil wetness for carbon decomposition', 'm3/m3', &
       missing_value=-100.0 )

  id_nfixrate = register_tiled_diag_field ( module_name, 'nfixrate', (/id_lon,id_lat/), &
       time, 'plant nitrogen fixing rate', 'kg N/(m2 year)', missing_value=-999.0)

end subroutine vegn_dynamics_init


! ============================================================================
subroutine vegn_carbon_int(vegn, soilt, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: soilt ! average temperature of soil for soil carbon decomposition, deg K
  type(diag_buff_type), intent(inout) :: diag

  type(vegn_cohort_type), pointer :: cc
  real :: resp, resl, resr, resg ! respiration terms accumulated for all cohorts 
  real :: md_leaf, md_root, md_sw, md_wood; ! components of maintenance demand (carbon loss)
  real :: md_alive ! = md_leaf+md_root
  real :: gpp ! gross primary productivity per tile
  real :: gpp_pot, npp_pot ! non-nitrogen limited GPP and NPP per tile
  real :: cc_gpp_pot, cc_npp_pot ! non-nitrogen limited GPP and NPP per cohort
  real :: nfixrate ! tile-level nitrogen fixation rate
  integer :: sp ! shorthand for current cohort specie
  integer :: i
  real :: CtoNleaf, CtoNroot, CtoNwood ! C:N ratios for different compartments

  if(is_watch_point()) then
     write(*,*)'#### vegn_carbon_int ####'
     __DEBUG1__(soilt)
  endif

  !  update plant carbon
  vegn%npp = 0
  resp = 0 ; resl = 0 ; resr = 0 ; resg = 0 ; gpp = 0
  gpp_pot = 0 ; npp_pot = 0 ;
  nfixrate = 0 ;
  do i = 1, vegn%n_cohorts   
     cc => vegn%cohorts(i)
     sp = cc%species

     call eddy_npp(cc,soilt,cc_gpp_pot,cc_npp_pot);
     ! npp2 is for diagnostics and comparison
     cc%npp2 = cc%miami_npp;  ! treat miami npp as above+below npp
     
     cc%carbon_gain = cc%carbon_gain + cc%npp*dt_fast_yr;
     
     ! check if leaves/roots are present and need to be accounted in maintenance
     if(cc%status == LEAF_ON) then
        md_leaf  = cc%Pl * spdata(sp)%alpha(CMPT_LEAF) * cc%bliving * dt_fast_yr
        md_root  = cc%Pr * spdata(sp)%alpha(CMPT_ROOT) * cc%bliving * dt_fast_yr
        ! md_alive = md_root + md_leaf, but it is recalculated to preserve
        ! bitwise compatibility with original code
        md_alive = ( cc%Pl * spdata(sp)%alpha(CMPT_LEAF) + &
                     cc%Pr * spdata(sp)%alpha(CMPT_ROOT))* &
             cc%bliving*dt_fast_yr
     else
        md_leaf  = 0
        md_root  = 0
        md_alive = 0
     endif
     
     ! compute sapwood maintenance
     md_sw = cc%Psw_alphasw * cc%bliving * dt_fast_yr;
     
     ! compute branch and coarse wood losses for tree types
     md_wood =0;
     if (sp > 1) then
        md_wood = 0.6 *cc%bwood * spdata(sp)%alpha(CMPT_WOOD)*dt_fast_yr;
     endif
     
     ! cc%md = md_leaf + md_root + md_sw;
     ! cc%bwood_gain = cc%bwood_gain + md_sw - md_wood;
     ! but we recalculate md_sw here to preserve bitwise compatibility with the 
     ! original code
     cc%md = md_alive + cc%Psw_alphasw * cc%bliving * dt_fast_yr;
     cc%bwood_gain = cc%bwood_gain + cc%Psw_alphasw * cc%bliving * dt_fast_yr;
     cc%bwood_gain = cc%bwood_gain - md_wood;
     if (cc%bwood_gain < 0.0) cc%bwood_gain=0.0; ! potential non-conservation ?
     cc%carbon_gain = cc%carbon_gain - cc%md;
     cc%carbon_loss = cc%carbon_loss + cc%md; ! used in diagnostics only

     if (do_nitrogen) then
       ! distribute litter carbon between fast and slow litter pools
       CtoNleaf = getCtoNleaf(cc)/(1-spdata(sp)%retrans_n)
       CtoNroot = getCtoNroot(cc)
       CtoNwood = getCtoNwood(cc)
       call litterfall(vegn,md_leaf,CtoNleaf)
       call litterfall(vegn,md_root,CtoNroot)
       call litterfall(vegn,md_wood,CtoNwood)
       call n_md(vegn,cc,diag, md_sw,md_leaf,CtoNleaf,md_root,CtoNroot,md_wood)
     else
       ! distribute litter carbon between fast and slow soil pools: litter
       ! pools are bypassed in this case
       vegn%fast_soil_C = vegn%fast_soil_C +    fsc_liv *md_alive +    fsc_wood *md_wood;
       vegn%slow_soil_C = vegn%slow_soil_C + (1-fsc_liv)*md_alive + (1-fsc_wood)*md_wood;
     endif

     ! for budget tracking
!/*     cp->fsc_in+= data->fsc_liv*md_alive+data->fsc_wood*md_wood; */
!/*     cp->ssc_in+= (1.- data->fsc_liv)*md_alive+(1-data->fsc_wood)*md_wood; */
     vegn%fsc_in  = vegn%fsc_in + 1*md_alive+0*md_wood;
     vegn%ssc_in  = vegn%ssc_in + (1.- 1)*md_alive+(1-0)*md_wood;

     vegn%veg_in  = vegn%veg_in  + cc%npp*dt_fast_yr;
     vegn%veg_out = vegn%veg_out + md_alive+md_wood;

     if(is_watch_point()) then
        __DEBUG4__(cc%bl, cc%br, cc%bsw, cc%bwood)
        __DEBUG3__(cc%An_op, cc%An_cl, cc%lai)
        __DEBUG1__(cc%species)
        __DEBUG2__(cc%npp, cc%gpp)
        __DEBUG4__(cc%resp, cc%resl, cc%resr, cc%resg)
        __DEBUG2__(cc%carbon_gain, cc%carbon_loss)
        __DEBUG1__(cc%bwood_gain)
        __DEBUG2__(vegn%fast_soil_C,vegn%slow_soil_C)
     endif
     ! accumulate tile-level NPP and GPP
     vegn%npp = vegn%npp + cc%npp
     gpp = gpp + cc%gpp
     gpp_pot = gpp_pot + cc_gpp_pot
     npp_pot = npp_pot + cc_npp_pot
     ! accumulate respiration terms for tile-level reporting
     resp = resp + cc%resp ; resl = resl + cc%resl
     resr = resr + cc%resr ; resg = resg + cc%resg
     ! accumulate nitrogen output for tile-level reporting
     nfixrate = nfixrate + cc%nfixrate
  enddo

  vegn%age = vegn%age + dt_fast_yr;

  ! ---- diagnostic section
  call send_tile_data(id_gpp,gpp,diag)
  call send_tile_data(id_npp,vegn%npp,diag)
  call send_tile_data(id_gpp_pot,gpp_pot,diag)
  call send_tile_data(id_npp_pot,npp_pot,diag)
  call send_tile_data(id_resp, resp, diag)
  call send_tile_data(id_resl, resl, diag)
  call send_tile_data(id_resr, resr, diag)
  call send_tile_data(id_resg, resg, diag)
  call send_tile_data(id_soilt,soilt,diag)
  call send_tile_data(id_nfixrate,nfixrate,diag)
  
end subroutine vegn_carbon_int


! ============================================================================
! updates cohort biomass pools, LAI, SAI, and height using accumulated 
! carbon_gain and bwood_gain
subroutine vegn_growth (vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  real :: wood_demand ! N for heartwood
  real :: excrete ! waste nitrogen that goes into litter
  real :: nsf_old ! saved status of plant nitrogen
  integer :: i

  do i = 1, vegn%n_cohorts   
     cc => vegn%cohorts(i)
     if(is_watch_point()) then
        write(*,*)'############### vegn_growth ##################'
        __DEBUG2__(cc%bwood,cc%bliving)
        __DEBUG3__(cc%bwood_gain,cc%carbon_gain,cc%carbon_loss)
        __DEBUG2__(cc%nwood,cc%nliving)
        __DEBUG2__(cc%nwoodgain,cc%ngain)
     endif
     cc%bwood   = cc%bwood   + cc%bwood_gain
     cc%bliving = cc%bliving + cc%carbon_gain
     cc%nwood   = cc%nwood   + cc%nwoodgain;
     cc%nliving = cc%nliving + cc%ngain;
     if(cc%bliving < 0) then
       cc%bwood    = cc%bwood+cc%bliving
       cc%bliving  = 0
       if (cc%bwood < 0) &
          cc%bwood = 0 ! in principle, that's not conserving carbon

       if (do_nitrogen) then
          ! ---- part from n_growth
          wood_demand = cc%bwood/spdata(cc%species)%CtoN(CMPT_WOOD) - cc%nwood;
          cc%nwood = cc%nwood + wood_demand;
          excrete  = cc%nstore + cc%nliving - wood_demand;
          if (excrete > 0) then
             ! dump all nitrogen remaining in living tissues into the litter
             call n_excrete(vegn, excrete)
             cc%nstore  = 0 
             cc%nliving = 0
          else
             ! for the sake of conservation, store negative amount of nitrogen 
             ! in nstore
             cc%nstore  = excrete
             cc%nliving = 0
          endif
      
          vegn%vegn_out = vegn%vegn_out + excrete;
          ! ---- end of part from n_growth
       endif
     endif
     
     call update_biomass_pools(vegn,cc)
     if (do_nitrogen) then
        nsf_old = cc%nsf
        ! define nitro storage capacity and current plant nitro status
        call n_reallocate(vegn, cc)
        ! calculate fixation
        call n_fixation(vegn, cc, nsf_old)
        ! distribute the nitrogen among nl, nr, nlv, nsw.
        ! handle the situation when cc%nstore is less then zero: dump excess
        ! of bliving into fast_litter_C
        ! slm: note that it doesn't update biomasses, as it probably should. Ask
        ! Stefan about it.
        call n_update_pool(vegn,cc)
     endif
     
     cc%root_density = (cc%br + &
            (cc%bsw+cc%bwood+cc%blv)*(1-agf_bs))*C2B
     cc%Wl_max = spdata(cc%species)%cmc_lai*cc%lai
     cc%Ws_max = spdata(cc%species)%csc_lai*cc%lai
     
     ! update leaf age
     if (cc%status == LEAF_ON) then
        cc%leaf_age = cc%leaf_age + 1.0
        ! limit the maximum leaf age by the leaf time span (reciprocal of leaf 
        ! turnover rate alpha) for given species. alpha is in 1/year, factor of
        ! 365 converts the result to days.
        if (spdata(cc%species)%alpha(CMPT_LEAF) > 0) &
             cc%leaf_age = min(cc%leaf_age,365.0/spdata(cc%species)%alpha(CMPT_LEAF))
     endif
     if (is_watch_point()) then
       __DEBUG2__(cc%bliving,cc%bwood)
       __DEBUG4__(cc%Pl,cc%Pr,cc%Psw,cc%Psw_alphasw)
       __DEBUG4__(cc%bl,cc%blv,cc%br,cc%bsw)
       __DEBUG2__(cc%height,cc%lai)   
     endif
     ! reset carbon and nitrogen accumulation terms
     cc%carbon_gain = 0
     cc%carbon_loss = 0
     cc%bwood_gain  = 0
     cc%ngain       = 0
     cc%nwoodgain   = 0
  end do

end subroutine vegn_growth


! ============================================================================
subroutine soil_carbon_int(vegn, soil_depth, soilt, theta, &
                           ndep_nit, ndep_amm, ndep_don, nfert, nfm_nit, nfm_amm, nfm_don, uptake, &
                           drainage, surf_runoff, diag)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: soil_depth ! depth of soil for soil nitrogen calculations, m
  real, intent(in) :: soilt ! average active layer soil temperature, deg K 
  real, intent(in) :: theta ! average active layer soil moisture, unitless
  real, intent(in) :: ndep_nit, ndep_amm, ndep_don ! nitrate, ammonium, and DON deposition, kg N/(m2 year)
  real, intent(in) :: nfert, nfm_nit, nfm_amm, nfm_don ! nitrogen fertilization, kg N/(m2 year)
  real, intent(in) :: uptake ! vegetation uptake from active layer, kg/(m2 s) 
                             ! -- for nitrogen uptake
  real, intent(in) :: drainage ! water drainage from active layer, kg/(m2 s)
  real, intent(in) :: surf_runoff ! surface runoff, kg/(m2 s)
  type(diag_buff_type), intent(inout) :: diag
  
  ! --- local vars
  real :: fast_C_loss
  real :: slow_C_loss
  real :: A  ! decomposition rate reduction due to moisture and temperature
  
  ! **** NITROGEN feedback on carbon decomposition **** 
  ! decomposition of slow litter depends on mineral nitrogen content
  ! done in nutrient_int, and we need fluxes there
  if (do_nitrogen) then
    call n_soil_int(vegn, diag, soil_depth, soilt, theta, &
       ndep_nit, ndep_amm, ndep_don, nfert, nfm_nit, nfm_amm, nfm_don, uptake, drainage, surf_runoff)
  else
    A=A_function(soilt,theta);
        
    fast_C_loss = vegn%fast_soil_C*A*K1*dt_fast_yr;
    slow_C_loss = vegn%slow_soil_C*A*K2*dt_fast_yr;
    
    vegn%fast_soil_C = vegn%fast_soil_C - fast_C_loss;
    vegn%slow_soil_C = vegn%slow_soil_C - slow_C_loss;
  
    ! for budget check
    vegn%fsc_out = vegn%fsc_out + fast_C_loss;
    vegn%ssc_out = vegn%ssc_out + slow_C_loss;
  
    ! loss of C to atmosphere and leaching
    vegn%rh =   (fast_C_loss+slow_C_loss)/dt_fast_yr;
  
    ! accumulate decomposition rate reduction for the soil carbon restart output
    vegn%asoil_in = vegn%asoil_in + A
  
    ! ---- diagnostic section
    call send_tile_data(id_rsoil_fast, fast_C_loss/dt_fast_yr, diag)
    call send_tile_data(id_asoil, A, diag)
    call send_tile_data(id_theta, theta, diag)
  endif
  call send_tile_data(id_fast_litter_C, vegn%fast_litter_C, diag)
  call send_tile_data(id_slow_litter_C, vegn%slow_litter_C, diag)
  call send_tile_data(id_fast_soil_C, vegn%fast_soil_C, diag)
  call send_tile_data(id_slow_soil_C, vegn%slow_soil_C, diag)
  call send_tile_data(id_res_soil_C, vegn%res_soil_C, diag)
  call send_tile_data(id_rsoil, vegn%rh, diag)
end subroutine


! ============================================================================
subroutine eddy_npp(cc, tsoil, cc_gpp_pot, cc_npp_pot)
  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(in)  :: tsoil  
  real, intent(out) :: cc_gpp_pot ! non-nitrogen-limited GPP
  real, intent(out) :: cc_npp_pot ! non-nitrogen-limited NPP

  call plant_respiration(cc,tsoil);

  cc%gpp = (cc%An_op - cc%An_cl)*mol_C*cc%lai;
  cc%npp = cc%gpp - cc%resp;

!  if(cc%npp_previous_day > -0.00001/2500.0) then
  if(cc%npp_previous_day > 0) then
     cc%resg = GROWTH_RESP*cc%npp_previous_day;
     cc%npp  = cc%npp  - GROWTH_RESP*cc%npp_previous_day;
     cc%resp = cc%resp + GROWTH_RESP*cc%npp_previous_day;
  else
     cc%resg = 0;
  endif

  ! calculate non-nitrogen-limited values of GPP and NPP
  ! this code was in n_eddy_npp
  if (do_nitrogen_plant_feedback .and. cc%gpp > 0 .and. cc%nsf > 0) then
    cc_gpp_pot = cc%gpp/n_phot(cc%nsf);
    cc_npp_pot = cc%npp + (cc_gpp_pot-cc%gpp)*(1.0-GROWTH_RESP);     
  else
    cc_gpp_pot = cc%gpp
    cc_npp_pot = cc%npp
  endif
  
  ! accumulate NPP values over the day
  cc%npp_previous_day_tmp     = cc%npp_previous_day_tmp + cc%npp
  cc%npp_previous_day_pot_tmp = cc%npp_previous_day_pot_tmp + cc_npp_pot
  ! end of code from n_eddy_npp
end subroutine eddy_npp


! ============================================================================
subroutine plant_respiration(cc, tsoil)
  type(vegn_cohort_type), intent(inout) :: cc
  real, intent(in) :: tsoil
  
  real :: tf,tfs;
  real :: r_leaf, r_vleaf, r_stem, r_root
  
  integer :: sp ! shorthand for cohort species
  sp = cc%species

  tf = exp(3000.0*(1.0/288.16-1.0/cc%prog%Tv));
  tf = tf / ( &
            (1.0+exp(0.4*(5.0-cc%prog%Tv+273.16)))*&
            (1.0+exp(0.4*(cc%prog%Tv - 273.16-45.0)))&
            )

  tfs = exp(3000.0*(1.0/288.16-1.0/tsoil));
  tfs = tfs / ( &
              (1.0+exp(0.4*(5.0-tsoil+273.16)))* &
              (1.0+exp(0.4*(tsoil - 273.16-45.0)))&
              )

  r_leaf = -mol_C*cc%An_cl*cc%lai;
  r_vleaf = spdata(sp)%beta(CMPT_VLEAF)   * cc%blv*tf;
  r_stem  = spdata(sp)%beta(CMPT_SAPWOOD) * cc%bsw*tf;
  r_root  = spdata(sp)%beta(CMPT_ROOT)    * cc%br*tfs;
  
  cc%resp = r_leaf + r_vleaf + r_stem + r_root;
  cc%resl = r_leaf;
  cc%resr = r_root;
end subroutine plant_respiration


! ============================================================================
! calculates prev. day average NPP from accumulated values
subroutine vegn_daily_npp(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  integer :: n_fast_step;
  integer :: i
  type(vegn_cohort_type), pointer :: cc

  n_fast_step = 1.0/365.0/dt_fast_yr;
  do i = 1, vegn%n_cohorts   
     cc => vegn%cohorts(i)
     vegn%cohorts(i)%npp_previous_day=vegn%cohorts(i)%npp_previous_day_tmp/n_fast_step;
     vegn%cohorts(i)%npp_previous_day_tmp=0.0
     vegn%cohorts(i)%npp_previous_day_pot=vegn%cohorts(i)%npp_previous_day_pot_tmp/n_fast_step;
     vegn%cohorts(i)%npp_previous_day_pot_tmp=0.0
  enddo
end subroutine vegn_daily_npp


! =============================================================================
subroutine vegn_phenology(vegn, wilt)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: wilt ! ratio of wilting to saturated water content

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc
  real :: leaf_litter,root_litter;    
  real :: theta_crit; ! critical ratio of average soil water to sat. water
  real :: CtoNleaf, CtoNroot ! C:N ratios of the respective tissues
  integer :: i
  
  vegn%litter = 0

  do i = 1,vegn%n_cohorts   
     cc => vegn%cohorts(i)
     
     if(is_watch_point())then
        write(*,*)'####### vegn_phenology #######'
        __DEBUG4__(vegn%theta_av, wilt, spdata(cc%species)%cnst_crit_phen, spdata(cc%species)%fact_crit_phen)
        __DEBUG1__(cc%species)
        __DEBUG2__(vegn%tc_av,spdata(cc%species)%tc_crit)
     endif
     ! if drought-deciduous or cold-deciduous species
     ! temp=10 degrees C gives good growing season pattern        
     ! spp=0 is c3 grass,1 c3 grass,2 deciduous, 3 evergreen
     ! assumption is that no difference between drought and cold deciduous
     cc%status = LEAF_ON; ! set status to indicate no leaf drop
      
     if(cc%species < 4 )then! deciduous species
        ! actually either fact_crit_phen or cnst_crit_phen is zero, enforced
        ! by logic in the vegn_data.F90
        theta_crit = spdata(cc%species)%cnst_crit_phen &
              + wilt*spdata(cc%species)%fact_crit_phen
        theta_crit = max(0.0,min(1.0, theta_crit))
        if ( (vegn%theta_av < theta_crit) &
             .or.(vegn%tc_av < spdata(cc%species)%tc_crit) ) then
           cc%status = LEAF_OFF; ! set status to indicate leaf drop 
           cc%leaf_age = 0;
           
           leaf_litter = (1.0-l_fract)*cc%bl;
           root_litter = (1.0-l_fract)*cc%br;
           
           if (do_nitrogen) then
              CtoNleaf = getCtoNleaf(cc)/(1-spdata(cc%species)%retrans_n);
              CtoNroot = getCtonroot(cc);
	      call litterfall (vegn,leaf_litter,CtoNleaf);
	      call litterfall (vegn,root_litter,CtoNroot);
           else
              ! add to patch litter flux terms
              vegn%litter = vegn%litter + leaf_litter + root_litter;
              
              vegn%fast_soil_C = vegn%fast_soil_C +    fsc_liv *(leaf_litter+root_litter);
              vegn%slow_soil_C = vegn%slow_soil_C + (1-fsc_liv)*(leaf_litter+root_litter);
                
              ! vegn%fsc_in+=data->fsc_liv*(leaf_litter+root_litter);
              ! vegn%ssc_in+=(1.0-data->fsc_liv)*(leaf_litter+root_litter);
              vegn%fsc_in  = vegn%fsc_in  + leaf_litter+root_litter;
              vegn%veg_out = vegn%veg_out + leaf_litter+root_litter;
           endif
           
           cc%blv = cc%blv + l_fract*(cc%bl+cc%br);
           cc%bl  = 0.0;
           cc%br  = 0.0;
           cc%lai = 0.0;
           
           ! update state
           cc%bliving = cc%blv + cc%br + cc%bl + cc%bsw;
           cc%b = cc%bliving + cc%bwood ;

           if (do_nitrogen) then
! ---- code from n_phenology(vegn,cc,leaf_litter,CtoNleaf,root_litter,CtoNroot)
             cc%nlv = cc%nlv + cc%nl + cc%nr &
                    - leaf_litter/CtoNleaf - root_litter/CtoNroot;
             cc%nr = 0.;
             cc%nl = 0.;
             cc%nliving = cc%nliving - (leaf_litter/CtoNleaf + root_litter/CtoNroot);

             ! budget
             vegn%vegn_out = vegn%vegn_out + leaf_litter/CtoNleaf + root_litter/CtoNroot;
! ---- end of code from n_phenology
             call update_bio_living_fraction(cc);   
             call n_update_pool(vegn,cc)
             ! because n_update_pool sometimes resets bliving
             cc%bsw = cc%bliving*cc%Psw;
             cc%blv = cc%bliving*(cc%Pr + cc%Pl);
           else
             call update_bio_living_fraction(cc);
           endif
!slm:suspect-code:	! sfg: br and bl are set already here and not in growth
!slm:suspect-code:        else if (cc%status /= LEAF_ON .and. new_status == LEAF_ON ) then
!slm:suspect-code:	    cc%bl = cc%bliving*cc%Pl;
!slm:suspect-code:	    cc%br = cc%bliving*cc%Pr;
!slm:suspect-code:	    cc%blv = 0.;
!slm:suspect-code:	    cc%lai = cc%bl*spdata(cc%species)%specific_leaf_area;
!slm:suspect-code:	    cc%status = new_status;
        endif ! leaf drop criteria: theta_crit ans tc_crit
     endif ! species < 4 (deciduous)
  enddo
end subroutine vegn_phenology


! =============================================================================
subroutine vegn_biogeography(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  integer :: i

  do i = 1, vegn%n_cohorts   
     call update_species(vegn%cohorts(i), vegn%t_ann, vegn%t_cold, &
          vegn%p_ann*seconds_per_year, vegn%ncm, vegn%landuse)
  enddo
end subroutine

! =============================================================================
! The stuff below comes from she_update.c -- it looks like it belongs here, 
! since it is essentially a part of the carbon integration (update_patch_fast
! is only called immediately after carbon_int in lm3v)
! =============================================================================


! =============================================================================
subroutine update_soil_pools(vegn)
  type(vegn_tile_type), intent(inout) :: vegn
  
  ! ---- local vars
  real :: delta(n_substances);

  ! update fsc input rate so that intermediate fsc pool is never
  ! depleted below zero; on the other hand the pool can be only 
  ! depleted, never increased
  vegn%fsc_rate = MAX( 0.0, MIN(vegn%fsc_rate, vegn%fsc_pool/dt_fast_yr));
  delta = vegn%fsc_rate*dt_fast_yr;
  if (do_nitrogen) then
    vegn%fast_litter_C = vegn%fast_litter_C + delta(isub_C)
    vegn%fast_litter_N = vegn%fast_litter_N + delta(isub_N)
    vegn%fln_in = vegn%fln_in+delta(isub_N)
  else
    vegn%fast_soil_C   = vegn%fast_soil_C   + delta(isub_C);
  endif
  vegn%fsc_pool = vegn%fsc_pool - delta;

  ! update ssc input rate so that intermediate ssc pool is never
  ! depleted below zero; on the other hand the pool can be only 
  ! depleted, never increased
  vegn%ssc_rate = MAX(0.0, MIN(vegn%ssc_rate, vegn%ssc_pool/dt_fast_yr));
  delta = vegn%ssc_rate*dt_fast_yr;
  if (do_nitrogen) then
    vegn%slow_litter_C = vegn%slow_litter_C + delta(isub_C)
    vegn%slow_litter_N = vegn%slow_litter_N + delta(isub_N)
    vegn%sln_in = vegn%sln_in + delta(isub_N)
  else
    vegn%slow_soil_C   = vegn%slow_soil_C   + delta(isub_C);
  endif
  vegn%ssc_pool = vegn%ssc_pool - delta;
end subroutine update_soil_pools



end module vegn_dynamics_mod
