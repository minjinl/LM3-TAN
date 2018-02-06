! ============================================================================
! vegetation disturbances
! ============================================================================
module vegn_disturbance_mod

use land_constants_mod, only : seconds_per_year
use vegn_data_mod,   only : spdata, fsc_wood, fsc_liv, agf_bs, LEAF_OFF, isub_C
use vegn_tile_mod,   only : vegn_tile_type, update_biomass_pools
use vegn_cohort_mod, only : vegn_cohort_type, height_from_biomass, lai_from_biomass
use nitrogen_data_mod, only : do_nitrogen
use nitrogen_mod,    only : getCtoNleaf, getCtoNroot, getCtoNsw, getCtoNwood, &
        litterfall, n_fire

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_nat_mortality
public :: vegn_disturbance
public :: update_fuel
! =====end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: &
     version = '$Id: vegn_disturbance.F90,v 17.0.2.5 2010/08/12 21:20:18 slm Exp $', &
     tagname = '$Name: nitro_20111003_slm $', &
     module_name = 'vegn_disturbance_mod'

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine vegn_disturbance(vegn, dt)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: dt ! time since last disturbance calculations, s
  

  real, parameter :: BMIN = 1e-10; ! should be the same as in growth function
  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  real :: delta;
  real :: fraction_lost;
  real :: deltat
  integer :: i
  integer :: sp ! shorthand for cohort species
  real :: deltaleaf, deltaroot, deltasw, deltawood ! increments of the tissue biomasses
  real :: delta1,delta2, cton1, cton2 ! variables to remove entire biomass if bliving < BMIN
  real :: CtoNleaf, CtoNroot, CtoNsw, CtoNwood ! tissue C:N ratios
  real :: conc_factor ! N concentration factor due to ash enrichment of nitrogen
  real :: vlf

  deltat = dt/seconds_per_year ! convert time interval to years

  !  Disturbance Rates
  vegn%disturbance_rate(0) = 0.0;
  vegn%disturbance_rate(1) = 0.0;
  
  call calculate_patch_disturbance_rates(vegn)
   
  ! Fire disturbance implicitly, i.e.  not patch creating
  vegn%area_disturbed_by_fire = (1.0-exp(-vegn%disturbance_rate(1)*deltat));
  
  do i = 1,vegn%n_cohorts   
     cc => vegn%cohorts(i)
     sp = cc%species

     fraction_lost = 1.0-exp(-vegn%disturbance_rate(1)*deltat);	
      
     ! "dead" biomass : wood + sapwood
     delta = (cc%bwood+cc%bsw)*fraction_lost;
     if (do_nitrogen) then
       deltasw   = fraction_lost*cc%bsw
       deltawood = fraction_lost*cc%bwood
     else
       vegn%slow_soil_C = vegn%slow_soil_C + (1.0-spdata(sp)%smoke_fraction)*delta*(1-fsc_wood);
       vegn%fast_soil_C = vegn%fast_soil_C + (1.0-spdata(sp)%smoke_fraction)*delta*   fsc_wood;
       ! for budget tracking - temporarily not keeping wood and the rest separately,ens
       !      vegn%ssc_in+=delta*(1.0-spdata(sp)%smoke_fraction)*(1-fsc_wood); */
       !      vegn%fsc_in+=delta*(1.0-spdata(sp)%smoke_fraction)*fsc_wood; */
       vegn%ssc_in = vegn%ssc_in+(cc%bwood+cc%bsw)*fraction_lost *(1.0-spdata(sp)%smoke_fraction);
       !     vegn%fsc_in+=cc%bsw*fraction_lost *(1.0-spdata(sp)%smoke_fraction);
     endif
     cc%bwood = cc%bwood * (1-fraction_lost);
     cc%bsw   = cc%bsw   * (1-fraction_lost);
      
     vegn%csmoke_pool(isub_C) = vegn%csmoke_pool(isub_C) + spdata(sp)%smoke_fraction*delta;
      
     vegn%veg_out = vegn%veg_out+delta;
      
     !"alive" biomass: leaves, roots, and virtual pool
     delta = (cc%bl+cc%blv+cc%br)*fraction_lost;
     if (do_nitrogen) then
       vlf = cc%Pl/(cc%Pl+cc%Pr)
       deltaleaf = fraction_lost*(cc%bl +    vlf *cc%blv)
       deltaroot = fraction_lost*(cc%br + (1-vlf)*cc%blv)
     else
       vegn%fast_soil_C = vegn%fast_soil_C + (1.0-spdata(sp)%smoke_fraction)*delta*    fsc_liv ;
       vegn%slow_soil_C = vegn%slow_soil_C + (1.0-spdata(sp)%smoke_fraction)*delta*(1- fsc_liv);
       ! for budget tracking- temporarily keeping alive separate ens
       ! /*      vegn%fsc_in+=delta* fsc_liv; */
       ! /*      vegn%ssc_in+=delta* (1-fsc_liv); */
       vegn%fsc_in = vegn%fsc_in+delta*(1.0-spdata(sp)%smoke_fraction);
     endif
     cc%bl  = cc%bl  * (1-fraction_lost);
     cc%blv = cc%blv * (1-fraction_lost);
     cc%br  = cc%br  * (1-fraction_lost);
      
     vegn%csmoke_pool(isub_C) = vegn%csmoke_pool(isub_C) + spdata(sp)%smoke_fraction*delta;
      
     vegn%veg_out = vegn%veg_out+delta;
      
     !"living" biomass:leaves, roots and sapwood
     delta = cc%bliving*fraction_lost;
     cc%bliving = cc%bliving - delta;

     
     if (do_nitrogen) then
       CtoNleaf = getCtoNleaf(cc);
       CtoNroot = getCtoNroot(cc);
       CtoNsw   = getCtoNsw(cc);
       CtoNwood = getCtoNwood(cc);
       call n_fire(vegn, cc, deltaleaf/CtoNleaf, deltaroot/CtoNroot, &
	     deltasw/CtoNsw, deltawood/CtoNwood, spdata(sp)%smoke_fraction, conc_factor);
       call litterfall(vegn,(1.0-spdata(sp)%smoke_fraction)*deltaleaf,conc_factor*CtoNleaf);
       call litterfall(vegn,(1.0-spdata(sp)%smoke_fraction)*deltaroot,conc_factor*CtoNroot);
       call litterfall(vegn,(1.0-spdata(sp)%smoke_fraction)*deltasw  ,conc_factor*CtoNsw);
       call litterfall(vegn,(1.0-spdata(sp)%smoke_fraction)*deltawood,conc_factor*CtoNwood);
     endif

     if(cc%bliving < BMIN) then
        ! remove vegetation completely 	      
        if (do_nitrogen) then
           delta1 = cc%bliving
           delta2 = cc%bwood
           if (delta1 > 0.) then 
              cton1 = delta1/(cc%nliving + cc%nstore)
              call litterfall(vegn,delta1,cton1)
           endif
           if (delta2 > 0.) then 
              cton2 = delta2/(cc%nwood)
              call litterfall(vegn,delta2,cton2)
           endif
           vegn%vegn_out = vegn%vegn_out + cc%nliving + cc%nstore + cc%nwood
           cc%nliving=0.
           cc%nstore =0.
           cc%nwood  =0.
        else
           vegn%fast_soil_C = vegn%fast_soil_C + fsc_liv*cc%bliving+ fsc_wood*cc%bwood;
           vegn%slow_soil_C = vegn%slow_soil_C + (1.- fsc_liv)*cc%bliving+ (1-fsc_wood)*cc%bwood;
           vegn%fsc_in = vegn%fsc_in + cc%bwood+cc%bliving;
        endif
        vegn%veg_out = vegn%veg_out + cc%bwood+cc%bliving;
        
        cc%bliving = 0.;
        cc%bwood   = 0.;
     endif
     
     call update_biomass_pools(vegn,cc)
  enddo

  vegn%csmoke_rate = vegn%csmoke_pool; ! kg C/(m2 yr)
end subroutine vegn_disturbance

! ============================================================================
subroutine calculate_patch_disturbance_rates(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  real :: fire_prob;
  real :: fuel;

  fuel = vegn%fuel

#if SIMPLE_FIRE
  ! CALCULATE FIRE DISTURBANCE RATES
  vegn%disturbance_rate(1)=fire(vegn);
#else

  ! lambda is the number of drought months;
  fire_prob = vegn%lambda/(1.+vegn%lambda); 
  ! compute average fuel during fire months
  if (vegn%lambda > 0.00001 ) fuel = fuel/vegn%lambda;
  vegn%disturbance_rate(1) = fuel * fire_prob;

  ! put a threshold for very dry years for warm places
  if (vegn%t_ann > 273.16 .and. vegn%lambda > 3.)  vegn%disturbance_rate(1)=0.33;
#endif  
  
  if(vegn%disturbance_rate(1) > 0.33) vegn%disturbance_rate(1)=0.33;
  
  ! this is only true for the one cohort per patch case
  vegn%disturbance_rate(0) = spdata(vegn%cohorts(1)%species)%treefall_disturbance_rate;
  vegn%total_disturbance_rate = vegn%disturbance_rate(1)+vegn%disturbance_rate(0);
  
  vegn%fuel = fuel;
end subroutine calculate_patch_disturbance_rates


! ============================================================================
function fire(vegn) result(fireterm)
  real :: fireterm; ! return value
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  real :: precip_av ! average precipitation, mm/year

  fireterm = 0
  precip_av = vegn%p_ann * seconds_per_year;
!!$  vegn%ignition_rate = 0.00;             
  vegn%fuel = vegn%total_biomass;

  if(vegn%fuel>0.0) then
     if(precip_av < 400.+40.*(vegn%t_ann-273.16)) then
        fireterm = vegn%fuel*(400. + 40.*(vegn%t_ann-273.16) - precip_av);
     endif
  endif
end function


! ============================================================================
subroutine update_fuel(vegn, wilt)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: wilt ! ratio of wilting to saturated water content

  ! ---- local constants
  !  these three used to be in data 
  real, parameter :: fire_height_threashold = 100;
  real, parameter :: fp1 = 1.; ! disturbance rate per kgC/m2 of fuel
  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  real :: theta_crit; ! critical ratio of average soil water to sat. water
  real :: ignition_rate;
  real ::  babove;
  integer :: i

  do i = 1,vegn%n_cohorts   
     cc => vegn%cohorts(i)
     ! calculate theta_crit: actually either fact_crit_fire or cnst_crit_fire 
     ! is zero, enforced by logic in the vegn_data.F90
     theta_crit = spdata(cc%species)%cnst_crit_fire &
           + wilt*spdata(cc%species)%fact_crit_fire
     theta_crit = max(0.0,min(1.0, theta_crit))
     if((cc%height < fire_height_threashold) &
          .and.(vegn%theta_av < theta_crit)  &
          .and.(vegn%tsoil_av > 278.16)) then
        babove = cc%bl + agf_bs * (cc%bsw + cc%bwood + cc%blv);
        ! this is fuel available during the drought months only
        vegn%fuel = vegn%fuel + spdata(cc%species)%fuel_intensity*babove;	
     endif
  enddo

  ! note that the ignition rate calculation based on the value of theta_crit for 
  ! the last cohort -- currently it doesn't matter since we have just one cohort, 
  ! but something needs to be done about that in the future
  ignition_rate = 0.;
  if ( (vegn%theta_av < theta_crit) &
       .and. (vegn%tsoil_av>278.16)) ignition_rate = 1.;
  vegn%lambda = vegn%lambda + ignition_rate;

end subroutine update_fuel




! ============================================================================
subroutine vegn_nat_mortality(vegn, deltat)
  type(vegn_tile_type), intent(inout) :: vegn
  real, intent(in) :: deltat ! time since last mortality calculations, s
  
  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc    ! current cohort
  real :: delta; ! total biomass change due to mortality
  real :: deltasw, deltawood ! changes in respective biomasses
  real :: fraction_lost; ! fraction of biomass lost due to mortality
  real :: bdead; ! combined biomass pool
  real :: CtoNwood, CtoNsw ! C:N ratios for respective tissues
  real :: dnwood, dnsw, dnstore ! changes in respective nitrogen pools
  integer :: i
  
  vegn%disturbance_rate(0)        = 0.0; 
  vegn%area_disturbed_by_treefall = 0.0;
  
  do i = 1,vegn%n_cohorts
     cc => vegn%cohorts(i)
     ! Treat tree fall disturbance implicitly, i.e. not creating a new tile.
     ! note that this disturbance rate calculation only works for the one cohort per 
     ! tile case -- in case of multiple cohort disturbance rate perhaps needs to be 
     ! accumulated (or averaged? or something else?) over the cohorts.
     vegn%disturbance_rate(0) = spdata(cc%species)%treefall_disturbance_rate;
     vegn%area_disturbed_by_treefall = &
          1.0-exp(-vegn%disturbance_rate(0)*deltat/seconds_per_year);

     ! calculate combined biomass pools
     bdead  = cc%bsw + cc%bwood;
     ! ens need a daily PATCH_FREQ here, for now it is set to 48
     fraction_lost = 1.0-exp(-vegn%disturbance_rate(0)*deltat/seconds_per_year);     
      
     ! "dead" biomass : wood + sapwood
     delta = bdead*fraction_lost;

     if (do_nitrogen) then
       CtoNsw    = getCtoNsw(cc);
       CtoNwood  = getCtoNwood(cc);
       ! changes in sapwood and wood due to mortality
       deltasw   = cc%bsw*fraction_lost;
       deltawood = cc%bwood*fraction_lost;
       ! flux into litter pools
       call litterfall(vegn,deltasw,CtoNsw);
       call litterfall(vegn,deltawood,CtoNwood);

       dnstore    = fraction_lost*MAX(cc%nstore,0.);
       dnsw       = deltasw/CtoNsw - dnstore;
       dnwood     = deltawood/CtoNwood;

       cc%nstore  = cc%nstore - dnstore;
       cc%nsw     = cc%nsw - dnsw;
       cc%nwood   = cc%nwood - dnwood;
       cc%nliving = cc%nliving - dnsw;

       vegn%vegn_out = vegn%vegn_out + dnsw + dnwood + dnstore;
     else
       ! no nitrogen
       vegn%slow_soil_C = vegn%slow_soil_C + (1-fsc_wood)*delta;
       vegn%fast_soil_C = vegn%fast_soil_C +    fsc_wood *delta;
  
       ! for budget tracking -temporarily
       ! vegn%fsc_in+= cc%bsw*fraction_lost;
       vegn%ssc_in  = vegn%ssc_in  + (cc%bwood+cc%bsw)*fraction_lost;
     endif

     cc%bwood = cc%bwood * (1-fraction_lost);
     cc%bsw   = cc%bsw   * (1-fraction_lost);
     
     vegn%veg_out = vegn%veg_out + delta;

     ! note that fast "living" pools are not included into mortality because their 
     ! turnover is calculated separately

     cc%bliving = cc%bsw + cc%bl + cc%br + cc%blv;
     call update_biomass_pools(vegn,cc);
  enddo
     
end subroutine vegn_nat_mortality


end module vegn_disturbance_mod
