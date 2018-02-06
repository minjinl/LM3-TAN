module vegn_harvesting_mod

use fms_mod, only : write_version_number, string, error_mesg, FATAL, NOTE, &
     mpp_pe, write_version_number, file_exist, open_namelist_file, close_file, &
     check_nml_error, stdlog, mpp_root_pe
use mpp_io_mod, only : axistype, mpp_get_atts, mpp_get_axis_data, &
     mpp_open, mpp_close, MPP_RDONLY, MPP_WRONLY, MPP_ASCII
use land_substances_mod, only : n_substances
use vegn_data_mod, only : &
     N_LU_TYPES, LU_PAST, LU_CROP, LU_NTRL, LU_SCND, &
     HARV_POOL_PAST, HARV_POOL_CROP, HARV_POOL_CLEARED, HARV_POOL_WOOD_FAST, &
     HARV_POOL_WOOD_MED, HARV_POOL_WOOD_SLOW, CMPT_ROOT, &
     spdata, agf_bs, fsc_liv, fsc_wood, isub_C
use vegn_tile_mod, only : &
     vegn_tile_type, update_biomass_pools, vegn_tile_nitrogen
use vegn_cohort_mod, only : &
     vegn_cohort_type
use nitrogen_data_mod, only : do_nitrogen, isub_N, CtoN_litter_min, &
     digestive_waste_frac
use nitrogen_mod, only : getCtoNleaf, getCtoNroot, getCtoNsw
use land_data_mod, only : land_state_type, lnd
use land_debug_mod, only : check_conservation

implicit none
private

! ==== public interface ======================================================
public :: vegn_harvesting_init
public :: vegn_harvesting_end

public :: vegn_harvesting

public :: vegn_graze_pasture
public :: vegn_harvest_cropland
public :: vegn_cut_forest
! ==== end of public interface ===============================================

! ==== module constants =====================================================
character(len=*), parameter   :: &
     version = '$Id: vegn_harvesting.F90,v 17.0.2.12 2010/07/30 00:50:09 slm Exp $', &
     tagname = '$Name: nitro_20111003_slm $', &
     module_name = 'vegn_harvesting_mod'
real, parameter :: ONETHIRD = 1.0/3.0

! ==== module data ==========================================================

! ---- namelist variables ---------------------------------------------------
logical :: do_harvesting       = .TRUE.  ! if true, then harvesting of crops and pastures is done
real :: grazing_intensity      = 0.25    ! fraction of biomass removed each time by grazing
real :: grazing_residue        = 0.1     ! fraction of the grazed biomass transferred into soil pools
real :: frac_wood_wasted_harv  = 0.25    ! fraction of wood wasted while harvesting
real :: frac_wood_wasted_clear = 0.25    ! fraction of wood wasted while clearing land for pastures or crops
real :: frac_wood_fast         = ONETHIRD ! fraction of wood consumed fast
real :: frac_wood_med          = ONETHIRD ! fraction of wood consumed with medium speed
real :: frac_wood_slow         = ONETHIRD ! fraction of wood consumed slowly
real :: crop_seed_density      = 0.1     ! biomass of seeds left after crop harvesting, kg/m2
namelist/harvesting_nml/ do_harvesting, grazing_intensity, grazing_residue, &
     frac_wood_wasted_harv, frac_wood_wasted_clear, &
     frac_wood_fast, frac_wood_med, frac_wood_slow, &
     crop_seed_density

contains ! ###################################################################

! ============================================================================
subroutine vegn_harvesting_init
  integer :: unit, ierr, io

  call write_version_number(version, tagname)

  if (file_exist('input.nml')) then
     unit = open_namelist_file ( )
     ierr = 1;  
     do while (ierr /= 0)
        read (unit, nml=harvesting_nml, iostat=io, end=10)
        ierr = check_nml_error (io, 'harvesting_nml')
     enddo
10   continue
     call close_file (unit)
  endif
  
  if (mpp_pe() == mpp_root_pe()) then
     unit=stdlog()
     write(unit, nml=harvesting_nml)
  endif

  if (frac_wood_fast+frac_wood_med+frac_wood_slow/=1.0) then
     call error_mesg('vegn_harvesting_init', &
          'sum of frac_wood_fast, frac_wood_med, and frac_wood_slow must be 1.0',&
          FATAL)
  endif
end subroutine vegn_harvesting_init


! ============================================================================
subroutine vegn_harvesting_end
end subroutine vegn_harvesting_end


! ============================================================================
! harvest vegetation in a tile
subroutine vegn_harvesting(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  if (.not.do_harvesting) &
       return ! do nothing if no harvesting requested

  select case(vegn%landuse)
  case(LU_PAST)  ! pasture
     call vegn_graze_pasture    (vegn)
  case(LU_CROP)  ! crop
     call vegn_harvest_cropland (vegn)
  end select
end subroutine


! ============================================================================
subroutine vegn_graze_pasture(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars 
  real :: bdead0, balive0; ! initial combined biomass pools
  real :: bdead1, balive1; ! updated combined biomass pools
  real :: bgrazed ! grazed biomass
  real :: dNliv ! decrement of living nitrogen
  real :: dNlit, dNslow ! increments of total and slow litter nitrogen, respectively
  real :: CtoNleaf, CtoNroot ! C:N ratios
  type(vegn_cohort_type), pointer :: cc ! shorthand for the current cohort
  integer :: i
  real :: ntot0, ntot1, epsilon = 1.e-14

  ! nitrogen conservation check part 1
  ntot0 =  vegn_tile_nitrogen(vegn) + vegn%nitrogen_out - vegn%nitrogen_in

  balive0 = 0 ; bdead0 = 0 ;
  balive1 = 0 ; bdead1 = 0 ;

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1,vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! calculate total biomass pools for the patch
     balive0 = balive0 + cc%bl + cc%blv + cc%br
     bdead0  = bdead0  + cc%bwood + cc%bsw
     ! only potential leaves are consumed
     vegn%harv_pool(HARV_POOL_PAST,isub_C) = vegn%harv_pool(HARV_POOL_PAST,isub_C) + &
          cc%bliving*cc%Pl*grazing_intensity*(1-grazing_residue) ;
     cc%bliving = cc%bliving - cc%bliving*cc%Pl*grazing_intensity;

     if (do_nitrogen) then
       bgrazed = cc%bliving*cc%Pl*grazing_intensity
       vegn%fsc_pool(isub_C) = vegn%fsc_pool(isub_C) + bgrazed*grazing_residue*   fsc_liv
       vegn%ssc_pool(isub_C) = vegn%ssc_pool(isub_C) + bgrazed*grazing_residue*(1-fsc_liv)

       CtoNleaf = getCtoNleaf(cc); ! grazed portion has C:N ratio of leaves
       CtoNroot = getCtoNroot(cc); ! and wasted portion -- C:N ratio of roots
       vegn%harv_pool(HARV_POOL_PAST,isub_N) = vegn%harv_pool(HARV_POOL_PAST,isub_N) &
           + bgrazed/CtoNleaf*(1-grazing_residue)*(1-digestive_waste_frac);
    
       ! litter flux: C:N of waste is assumed to have root ratio
       dNliv = bgrazed*((1-grazing_residue)/CtoNleaf+grazing_residue/CtoNroot);
       dNlit = bgrazed*((1-grazing_residue)/CtoNleaf*digestive_waste_frac &
                        +grazing_residue/CtoNroot)
       dNslow = min(bgrazed*grazing_residue*(1-fsc_liv)/CtoN_litter_min,dNlit)
       cc%nliving = cc%nliving - dNliv;
       vegn%ssc_pool(isub_N) = vegn%ssc_pool(isub_N) + dNslow
       vegn%fsc_pool(isub_N) = vegn%fsc_pool(isub_N) - dNslow + dNlit
     endif
     ! redistribute leftover biomass among biomass pools
     call update_biomass_pools(vegn,cc);
 
     ! calculate new combined vegetation biomass pools
     balive1 = balive1 + cc%bl + cc%blv + cc%br
     bdead1  = bdead1  + cc%bwood + cc%bsw
  enddo

  if (.not.do_nitrogen) then
    ! update intermediate soil carbon pools
    vegn%fsc_pool(isub_C) = vegn%fsc_pool(isub_C) + &
         (   fsc_liv *(balive0-balive1)+    fsc_wood *(bdead0-bdead1))*grazing_residue;
    vegn%ssc_pool(isub_C) = vegn%ssc_pool(isub_C) + &
         ((1-fsc_liv)*(balive0-balive1)+ (1-fsc_wood)*(bdead0-bdead1))*grazing_residue;
  endif

  ! nitrogen conservation check part 2
  ntot1 =  vegn_tile_nitrogen(vegn) + vegn%nitrogen_out - vegn%nitrogen_in
  !call check_conservation('vegn_graze_pasture','nitrogen',ntot0,ntot1,epsilon,lnd%time)
  
end subroutine vegn_graze_pasture


! ================================================================================
subroutine vegn_harvest_cropland(vegn)
  type(vegn_tile_type), intent(inout) :: vegn

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc ! pointer to the current cohort
  real :: frac_harvested;    ! fraction of biomass harvested this time
  real :: bdead, balive, btotal; ! combined biomass pools
  integer :: i
  real :: wasted_living_C, wasted_wood_C
  real :: ntot0, ntot1, epsilon = 1.e-14

  ! nitrogen conservation check part 1
  ntot0 =  vegn_tile_nitrogen(vegn) + vegn%nitrogen_out - vegn%nitrogen_in
  
  balive = 0 ; bdead = 0
  ! calculate initial combined biomass pools for the patch
  do i = 1, vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! calculate total biomass pools for the patch
     balive = balive + cc%bl + cc%blv + cc%br
     bdead  = bdead  + cc%bwood + cc%bsw
  enddo
  btotal = balive+bdead;

  ! calculate harvested fraction: cut everything down to seed level
  frac_harvested = MIN(MAX((btotal-crop_seed_density)/btotal,0.0),1.0);

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1, vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! use for harvest only above ground living biomass and waste the correspondent below living and wood
     vegn%harv_pool(HARV_POOL_CROP,isub_C) = vegn%harv_pool(HARV_POOL_CROP,isub_C) + &
          cc%bliving*(cc%Pl + cc%Psw*agf_bs)*frac_harvested;
     if (.not.do_nitrogen) then
       vegn%fsc_pool(isub_C) = vegn%fsc_pool(isub_C) &
          + frac_harvested*(   fsc_liv *cc%bliving*cc%Pr +    fsc_wood *(cc%bwood + cc%bliving*cc%Psw*(1-agf_bs)));
       vegn%ssc_pool(isub_C) = vegn%ssc_pool(isub_C) &
          + frac_harvested*((1-fsc_liv)*cc%bliving*cc%Pr + (1-fsc_wood)*(cc%bwood + cc%bliving*cc%Psw*(1-agf_bs)));
     else
       vegn%harv_pool(HARV_POOL_CROP,isub_N) = vegn%harv_pool(HARV_POOL_CROP,isub_N) &
          + cc%bliving*frac_harvested*(cc%Pl/getCtoNleaf(cc) + cc%Psw*agf_bs/getCtoNsw(cc));
       wasted_living_C = frac_harvested*cc%bliving*(cc%Pr+cc%Psw*(1-agf_bs))
       wasted_wood_C   = frac_harvested*cc%bwood
       vegn%fsc_pool(isub_C) = vegn%fsc_pool(isub_C) &
          +    fsc_liv *wasted_living_C +    fsc_wood *wasted_wood_C
       vegn%ssc_pool(isub_C) = vegn%ssc_pool(isub_C) &
          + (1-fsc_liv)*wasted_living_C + (1-fsc_wood)*wasted_wood_C
       vegn%ssc_pool(isub_N) = vegn%ssc_pool(isub_N) &
          + wasted_living_C/CtoN_litter_min + cc%nwood*frac_harvested*(1-fsc_wood)
       vegn%fsc_pool(isub_N) = vegn%fsc_pool(isub_N) &
          + cc%bliving*frac_harvested*(cc%Pr/getCtoNroot(cc) + cc%Psw*(1-agf_bs)/getCtoNsw(cc)) &
          - wasted_living_C/CtoN_litter_min + cc%nwood*frac_harvested*fsc_wood
       cc%nliving = cc%nliving * (1-frac_harvested);
       if (cc%nstore>0) &
            cc%nstore = cc%nstore * (1-frac_harvested);
       cc%nwood   = cc%nwood   * (1-frac_harvested);
     endif

     cc%bliving = cc%bliving * (1-frac_harvested);
     cc%bwood   = cc%bwood   * (1-frac_harvested);
     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(vegn,cc);
  enddo

  ! nitrogen conservation check part 2
  ntot1 =  vegn_tile_nitrogen(vegn) + vegn%nitrogen_out - vegn%nitrogen_in
  !call check_conservation('vegn_harvest_cropland','nitrogen',ntot0,ntot1,epsilon,lnd%time)
  
end subroutine vegn_harvest_cropland


! ============================================================================
! for now cutting forest is the same as harvesting cropland --
! we basically cut down everything, leaving only seeds
subroutine vegn_cut_forest(vegn, new_landuse)
  type(vegn_tile_type), intent(inout) :: vegn
  integer, intent(in) :: new_landuse ! new land use type that gets assigned to
                                     ! the tile after the wood harvesting

  ! ---- local vars
  type(vegn_cohort_type), pointer :: cc ! pointer to the current cohort
  real :: frac_harvested;        ! fraction of biomass harvested this time
  real :: frac_wood_wasted       ! fraction of wood wasted during transition
  real :: wood_harvested(n_substances)  ! amount of harvested wood, kgC/m2
  real :: bdead, balive, btotal; ! combined biomass pools
  real :: delta(n_substances)
  integer :: i
  real :: ntot0, ntot1, epsilon = 1.e-14

  ! nitrogen conservation check part 1
  ntot0 =  vegn_tile_nitrogen(vegn) + vegn%nitrogen_out - vegn%nitrogen_in
  
  balive = 0 ; bdead = 0
  ! calculate initial combined biomass pools for the patch
  do i = 1, vegn%n_cohorts
     cc=>vegn%cohorts(i)
     ! calculate total biomass pools for the patch
     balive = balive + cc%bl + cc%blv + cc%br
     bdead  = bdead  + cc%bwood + cc%bsw
  enddo
  btotal = balive+bdead;

  ! calculate harvested fraction: cut everything down to seed level
  frac_harvested = MIN(MAX((btotal-crop_seed_density)/btotal,0.0),1.0);

  ! define fraction of wood wasted, based on the transition type
  if (new_landuse==LU_SCND) then
     frac_wood_wasted = frac_wood_wasted_harv
  else
     frac_wood_wasted = frac_wood_wasted_clear
  endif

  ! update biomass pools for each cohort according to harvested fraction
  do i = 1, vegn%n_cohorts
     cc => vegn%cohorts(i)

     ! calculate total amount of harvested wood, minus the wasted part
     wood_harvested(isub_C) = (cc%bwood+cc%bsw)*frac_harvested*(1-frac_wood_wasted)
     if (do_nitrogen) &
       wood_harvested(isub_N) = (cc%nwood+cc%bsw/getCtoNsw(cc)) &
                              * frac_harvested*(1-frac_wood_wasted)
     
     ! distribute harvested wood between pools
     if (new_landuse==LU_SCND) then
        ! this is harvesting, distribute between 3 different wood pools 
        vegn%harv_pool(HARV_POOL_WOOD_FAST,:) = vegn%harv_pool(HARV_POOL_WOOD_FAST,:) &
             + wood_harvested*frac_wood_fast
        vegn%harv_pool(HARV_POOL_WOOD_MED,:)  = vegn%harv_pool(HARV_POOL_WOOD_MED,:)  &
             + wood_harvested*frac_wood_med
        vegn%harv_pool(HARV_POOL_WOOD_SLOW,:) = vegn%harv_pool(HARV_POOL_WOOD_SLOW,:) &
             + wood_harvested*frac_wood_slow
     else
        ! this is land clearance (clearing land for pastures or crops): 
        ! everything goes into "cleared" pool
        vegn%harv_pool(HARV_POOL_CLEARED,:) = vegn%harv_pool(HARV_POOL_CLEARED,:) &
             + wood_harvested
     endif

     ! distribute wasted wood between fast and slow intermediate soil carbon 
     ! buffers according to fractions specified thorough the namelists
     if (.not.do_nitrogen) then
        delta(isub_C) = (cc%bwood+cc%bsw)*frac_harvested*frac_wood_wasted
     else
        delta(isub_C) = cc%bwood*frac_harvested*frac_wood_wasted
        delta(isub_N) = cc%nwood*frac_harvested*frac_wood_wasted
        ! sapwood is treated with living biomass in this case
     endif
     vegn%ssc_pool = vegn%ssc_pool + delta*(1-fsc_wood);
     vegn%fsc_pool = vegn%fsc_pool + delta*   fsc_wood ;

     ! distribute wasted living between fast and slow intermediate soil carbon 
     ! buffers according to fractions specified thorough the namelists
     if (.not.do_nitrogen) then
        delta(isub_C) = balive * frac_harvested;
     else
       delta(isub_C) = balive*frac_harvested + cc%bsw*frac_harvested*frac_wood_wasted
       delta(isub_N) = frac_harvested * &
                ( cc%bliving*(cc%Pl/getCtoNleaf(cc)+cc%Pr/getCtoNroot(cc))  &
                + cc%bsw/getCtoNsw(cc)*frac_wood_wasted )
       vegn%ssc_pool(isub_N) = vegn%ssc_pool(isub_N) &
                     + min(delta(isub_C)*(1-fsc_liv)/CtoN_litter_min,delta(isub_N))
       vegn%fsc_pool(isub_N) = vegn%fsc_pool(isub_N) &
                     - min(delta(isub_C)*(1-fsc_liv)/CtoN_litter_min,delta(isub_N)) &
                     + delta(isub_N)
     endif
     vegn%ssc_pool(isub_C) = vegn%ssc_pool(isub_C) + delta(isub_C)*(1-fsc_liv)
     vegn%fsc_pool(isub_C) = vegn%fsc_pool(isub_C) + delta(isub_C)*   fsc_liv
     
     ! note that in nitrogen code the wasted sapwood treated differently from
     ! the wasted wood : this is somewhat inconsistent with the treatment of
     ! the harvested wood and sapwood amounts.
     
     cc%bliving = cc%bliving * (1-frac_harvested);
     cc%bwood   = cc%bwood   * (1-frac_harvested);
     if (do_nitrogen) then
       cc%nliving = cc%nliving * (1-frac_harvested);
       if (cc%nstore>0) &
            cc%nstore = cc%nstore * (1-frac_harvested);
       cc%nwood   = cc%nwood   * (1-frac_harvested);
     endif
     ! redistribute leftover biomass between biomass pools
     call update_biomass_pools(vegn,cc);
  enddo

  ! nitrogen conservation check part 2
  ntot1 =  vegn_tile_nitrogen(vegn) + vegn%nitrogen_out - vegn%nitrogen_in
  !call check_conservation('vegn_cut_forest','nitrogen',ntot0,ntot1,epsilon,lnd%time)
  
end subroutine vegn_cut_forest

end module 
