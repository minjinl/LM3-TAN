module vegn_tile_mod

#include "../shared/debug.inc"

use fms_mod, only : &
     write_version_number, stdlog, error_mesg, FATAL, WARNING
use constants_mod, only : &
     tfreeze, hlf

use land_constants_mod, only : NBANDS
use land_io_mod, only : &
     init_cover_field
use land_tile_selectors_mod, only : &
     tile_selector_type, SEL_VEGN
use land_debug_mod, only : is_watch_point, get_current_point

use vegn_data_mod, only : &
     NSPECIES, MSPECIES, NCMPT, C2B, &
     read_vegn_data_namelist, spdata, &
     vegn_to_use,  input_cover_types, &
     mcv_min, mcv_lai, &
     vegn_index_constant, &
     agf_bs, BSEED, LU_NTRL, LU_SCND, N_HARV_POOLS, &
     LU_SEL_TAG, SP_SEL_TAG, NG_SEL_TAG, &
     SP_C3GRASS, SP_C4GRASS, SP_TROPICAL, &
     LEAF_OFF, LEAF_ON, &
     CMPT_LEAF, CMPT_ROOT, CMPT_SAPWOOD, CMPT_WOOD, &
     isub_C, &
     scnd_biomass_bins

use vegn_cohort_mod, only : vegn_cohort_type, vegn_phys_prog_type, &
     height_from_biomass, lai_from_biomass, update_bio_living_fraction, &
     cohort_uptake_profile, cohort_root_properties

use nitrogen_data_mod, only : do_nitrogen, do_nitrogen_plant_feedback, &
     planning_horizon, CtoN_seed, n_deviation_cutoff, isub_N

use land_substances_mod, only : n_substances

implicit none
private

! ==== public interfaces =====================================================
public :: vegn_tile_type

public :: new_vegn_tile, delete_vegn_tile
public :: vegn_tiles_can_be_merged, merge_vegn_tiles
public :: vegn_is_selected
public :: get_vegn_tile_tag
public :: vegn_tile_stock_pe
public :: vegn_tile_carbon ! returns total carbon per tile
public :: vegn_tile_nitrogen ! returns total nitrogen per tile
public :: vegn_tile_heat ! returns hate content of the vegetation
public :: vegn_tile_bliving ! returns living carbon per tile
public :: vegn_tile_nliving ! returns living nitrogen per tile

public :: read_vegn_data_namelist
public :: vegn_cover_cold_start

public :: vegn_uptake_profile
public :: vegn_root_properties
public :: vegn_data_rs_min
public :: vegn_seed_supply
public :: vegn_seed_demand

public :: vegn_tran_priority ! returns transition priority for land use 

public :: vegn_add_bliving
public :: update_derived_vegn_data  ! given state variables, calculate derived values
public :: update_biomass_pools, n_update_pool, n_store_optimum
! =====end of public interfaces ==============================================
interface new_vegn_tile
   module procedure vegn_tile_ctor
   module procedure vegn_tile_copy_ctor
end interface


! ==== module constants ======================================================
character(len=*), parameter   :: &
     version = '$Id: vegn_tile.F90,v 18.0.2.1.2.12.2.4.2.1 2011/10/04 00:03:24 slm Exp $', & 
     tagname = '$Name: nleach_experimental_pcm_m2l_slm $', &
     module_name = 'vegn_tile_mod'

! ==== types =================================================================
type :: vegn_tile_type
   integer :: tag ! kind of the tile
   integer :: landuse = LU_NTRL

   integer :: n_cohorts = 0
   type(vegn_cohort_type), pointer :: cohorts(:)=>NULL()

   real :: age=0 ! tile age

   ! fast litter pool
   real :: fast_litter_C    = 0.0 ! fast litter carbon pool, (kg C/m2)
   real :: fast_litter_N    = 0.0 ! fast litter nitrogen pool, (kg N/m2)
   real :: fast_litter_fsol = 0.0 ! fraction of soluble nitrogen in the fast litter pool
   ! slow litter pool
   real :: slow_litter_C    = 0.0 ! slow litter carbon pool, (kg C/m2)
   real :: slow_litter_N    = 0.0 ! slow litter nitrogen pool, (kg N/m2)
   real :: slow_litter_fsol = 0.5 ! fraction of soluble nitrogen in the slow litter pool
   ! fast soil pool
   real :: fast_soil_C      = 0.0 ! fast soil carbon pool, (kg C/m2)
   ! slow soil pool
   real :: slow_soil_C      = 0.0 ! slow soil carbon pool, (kg C/m2)
   real :: slow_soil_N      = 0.0 ! slow soil nitrogen pool, (kg N/m2)
   real :: slow_soil_fsol   = 0.5 ! fraction of soluble nitrogen in the slow soil pool
   ! passive (resilient?) soil pool
   real :: res_soil_C       = 0.0 ! resistant soil carbon pool, (kg C/m2)
   real :: res_soil_N       = 0.0 ! resistant soil nitrogen pool, (kg N/m2)
   ! mineral nitrogen pools
   real :: amm              = 0.0 ! ammonium, (kg N/m2)
   real :: nitr             = 0.0 ! nitrate, (kg N/m2)
   
   ! fields for smoothing out the contribution of the spike-type processes (e.g. 
   ! harvesting) to the soil carbon pools over some period of time
   real, pointer :: fsc_pool(:)=>NULL(), fsc_rate(:)=>NULL() ! for fast soil carbon
   real, pointer :: ssc_pool(:)=>NULL(), ssc_rate(:)=>NULL() ! for slow soil carbon

   real, pointer :: csmoke_pool(:)=>NULL() ! substance lost through fires, kg X/m2 
   real, pointer :: csmoke_rate(:)=>NULL() ! rate of release of the above to atmosphere, kg X/(m2 yr)

   real, pointer :: harv_pool(:,:)=>NULL() ! pools of harvested carbon, kg C/m2 (N_HARV_POOLS, SUBSTANCES)
   real, pointer :: harv_rate(:,:)=>NULL() ! rates of spending (release to the atmosphere), kg C/(m2 yr) (N_HARV_POOLS, SUBSTANCES)

   ! values for the diagnostic of carbon budget and soil carbon acceleration
   real :: asoil_in=0
   real :: ssc_in=0, ssc_out=0
   real :: resc_in=0, resc_out=0
   real :: fsc_in=0, fsc_out=0
   real :: veg_in=0, veg_out=0
   
   real :: flc_in=0, flc_out=0
   real :: slc_in=0, slc_out=0
   real :: sln_in=0, sln_out=0, sln_in0=0
   real :: ssn_in=0, ssn_out=0
   real :: resn_in=0, resn_out=0
   real :: fln_in=0, fln_out=0, fln_in0=0
   real :: vegN_in=0, vegN_out = 0
   real :: nitrogen_in=0, nitrogen_out=0
   ! NOTE: fln_in0 and sln_in0 are used to calculate the *rates* of nitrogen
   ! input to fast and slow litter pools. 
   ! Since fast_litter_N is updated in "litterfall", which is called in many
   ! places (and also in "n_excrete"), it was easier to store the the value of
   ! the cumulative fln_in at the start of time step in fln_in0, then at the end
   ! of the same time step calculate the rate as the (fln_in-fln_in0)/dt
   ! The same logic applies to slow_soil_N, sln_in and sln_in0

   real :: disturbance_rate(0:1) = 0 ! 1/year
   real :: lambda   ! cumulative drought months per year
   real :: fuel     ! fuel over dry months
   real :: litter   ! litter flux

   ! for equilibrium calculations
   real :: A_eq=0
   real :: K_ls_eq=0, K_ss_eq=0
   real :: slc_in_eq=0,  slc_out_eq
   real :: ssc_in_eq=0,  ssc_out_eq
   real :: resc_in_eq=0, resc_out_eq
   real :: sln_in_eq=0
   
   ! monthly accumulated/averaged values
   real :: theta_av ! relative soil_moisture availability not soil moisture
   real :: tsoil_av ! bulk soil temperature
   real :: tc_av    ! leaf temperature
   real :: precip_av! precipitation

   ! accumulation counters for long-term averages (monthly and annual). Having
   ! these counters in the tile is a bit stupid, since the values are the same for
   ! each tile, but it simplifies the current code, and they are going away when we
   ! switch to exponential averaging in any case.
   integer :: n_accum ! number of accumulated values for monthly averages
   integer :: nmn_acm ! number of accumulated values for annual averages
   ! annual-mean values
   real :: t_ann    ! annual mean T, degK
   real :: t_cold   ! average temperature of the coldest month, degK
   real :: p_ann    ! annual mean precipitation
   real :: ncm      ! number of cold months
   ! annual accumulated values
   real :: t_ann_acm  ! accumulated annual temperature for t_ann
   real :: t_cold_acm ! temperature of the coldest month in current year
   real :: p_ann_acm  ! accumulated annual precipitation for p_ann
   real :: ncm_acm    ! accumulated number of cold months


   ! it's probably possible to get rid of the fields below
   real :: npp=0 ! net primary productivity
   real :: nep=0 ! net ecosystem productivity
   real :: rh=0 ! soil carbon lost to the atmosphere
   real :: total_biomass !
   real :: area_disturbed_by_treefall
   real :: area_disturbed_by_fire
   real :: total_disturbance_rate
   real :: donleach  ! leaching of dissolved organic nitrogen, kg N/(m2 s)
   real :: ammleach  ! leaching of ammonium, kg N/(m2 s)
   real :: nitrleach ! leaching of nitrate, kg N/(m2 s)
end type vegn_tile_type

! ==== module data ===========================================================
real, public :: &
     cpw = 1952.0, & ! specific heat of water vapor at constant pressure
     clw = 4218.0, & ! specific heat of water (liquid)
     csw = 2106.0    ! specific heat of water (ice)

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! ============================================================================
function vegn_tile_ctor(tag) result(ptr)
  type(vegn_tile_type), pointer :: ptr ! return value
  integer, intent(in)  :: tag ! kind of tile

  allocate(ptr)
  ptr%tag = tag
  allocate(ptr%harv_pool(N_HARV_POOLS,n_substances)); ptr%harv_pool = 0
  allocate(ptr%harv_rate(N_HARV_POOLS,n_substances)); ptr%harv_rate = 0
  allocate(ptr%fsc_pool(n_substances)); ptr%fsc_pool = 0
  allocate(ptr%fsc_rate(n_substances)); ptr%fsc_rate = 0
  allocate(ptr%ssc_pool(n_substances)); ptr%ssc_pool = 0
  allocate(ptr%ssc_rate(n_substances)); ptr%ssc_rate = 0
  allocate(ptr%csmoke_pool(n_substances)); ptr%csmoke_pool = 0
  allocate(ptr%csmoke_rate(n_substances)); ptr%csmoke_rate = 0
  
end function vegn_tile_ctor

! ============================================================================
function vegn_tile_copy_ctor(vegn) result(ptr)
  type(vegn_tile_type), pointer :: ptr ! return value
  type(vegn_tile_type), intent(in) :: vegn ! return value

  allocate(ptr)
  ! copy all non-pointer members
  ptr=vegn
  ! copy pointer members (cohorts)
  allocate(ptr%cohorts(ptr%n_cohorts))
  ptr%cohorts(:) = vegn%cohorts(1:ptr%n_cohorts)
  ! copy harvesting pools and their spending rates
  allocate(ptr%harv_pool(N_HARV_POOLS,n_substances)); ptr%harv_pool = vegn%harv_pool
  allocate(ptr%harv_rate(N_HARV_POOLS,n_substances)); ptr%harv_rate = vegn%harv_rate
  ! copy intermediate pools and their spending rates
  allocate(ptr%fsc_pool(n_substances)); ptr%fsc_pool = vegn%fsc_pool
  allocate(ptr%fsc_rate(n_substances)); ptr%fsc_rate = vegn%fsc_rate
  allocate(ptr%ssc_pool(n_substances)); ptr%ssc_pool = vegn%ssc_pool
  allocate(ptr%ssc_rate(n_substances)); ptr%ssc_rate = vegn%ssc_rate
  allocate(ptr%csmoke_pool(n_substances)); ptr%csmoke_pool = vegn%csmoke_pool
  allocate(ptr%csmoke_rate(n_substances)); ptr%csmoke_rate = vegn%csmoke_rate
  
end function vegn_tile_copy_ctor

! ============================================================================
subroutine delete_vegn_tile(vegn)
  type(vegn_tile_type), pointer :: vegn

  deallocate(vegn%cohorts)
  deallocate(vegn%harv_pool)
  deallocate(vegn%harv_rate)
  deallocate(vegn%fsc_pool)
  deallocate(vegn%fsc_rate)
  deallocate(vegn%ssc_pool)
  deallocate(vegn%ssc_rate)
  deallocate(vegn%csmoke_pool)
  deallocate(vegn%csmoke_rate)  
  deallocate(vegn)
end subroutine delete_vegn_tile

! =============================================================================
function vegn_tiles_can_be_merged(vegn1,vegn2) result(response)
  logical :: response
  type(vegn_tile_type), intent(in) :: vegn1,vegn2

  real    :: b1, b2 
  integer :: i, i1, i2

  if (vegn1%landuse /= vegn2%landuse) then
     response = .false. ! different land use tiles can't be merged
  else if (vegn1%landuse == LU_SCND) then ! secondary vegetation tiles
     ! get tile wood biomasses
     b1 = get_vegn_tile_bwood(vegn1)
     b2 = get_vegn_tile_bwood(vegn2)
     ! find biomass bins where each the tiles belongs to
     i1 = 0 ; i2 = 0
     do i = 1, size(scnd_biomass_bins(:))
        if (b1>scnd_biomass_bins(i)) i1 = i
        if (b2>scnd_biomass_bins(i)) i2 = i
     enddo
     ! tiles can be merged only if biomasses belong to the same bin
     response = (i1 == i2)
     ! and then only if their nitrogen status is compatible
     response = response.and.n_tiles_can_be_merged(vegn1,vegn2)
  else
     response = .true. ! non-secondary tiles of the same land use type can always be merged
  endif
end function


! =============================================================================
! returns true is the nitrogen status in two tiles is close enough to be merged
function n_tiles_can_be_merged(vegn1,vegn2) result(response)
  logical :: response
  type(vegn_tile_type), intent(in) :: vegn1,vegn2

  ! ---- local vars
  real :: b1, b2, h1
  
  if (.not.do_nitrogen) then
     response = .TRUE.
     return
  endif
  
  b1 = n_status_for_merge(vegn1)
  b2 = n_status_for_merge(vegn2)
  
  ! tiles can be merged only if fixation is approximatively the same
  ! for extratropics
  if (b1+b2 > 0.0) then
     h1 = 2*abs(b1-b2)/abs(b1+b2) ! relative difference
     response = (h1 < n_deviation_cutoff .or. max(b1,b2) < 0.3e-4)
  else
     ! at least one of the species is tropical tree, allow merge
     response = .TRUE.
  endif
  
end function


! ============================================================================
! this is the nitrogen criterion for merge
! allow merge if 1) no n feedback -or- 2) tropical tree
! the decision of merge is based on difference of the value returned by this function
! thus return -1 for case 1) and 2) and nfix_per_bliv for other cases
real function n_status_for_merge(vegn)
  type(vegn_tile_type), intent(in) :: vegn
  
  if (vegn%cohorts(1)%species==SP_TROPICAL) then
    n_status_for_merge = -1.0
  else
    n_status_for_merge = vegn%cohorts(1)%nfix_per_bliv*vegn%cohorts(1)%bliving
  endif
end function


! ============================================================================
subroutine merge_vegn_tiles(t1,w1,t2,w2)
  type(vegn_tile_type), intent(in) :: t1
  type(vegn_tile_type), intent(inout) :: t2
  real, intent(in) :: w1, w2 ! relative weights
  
  ! ---- local vars
  real :: x1, x2 ! normalized relative weights
  real :: HEAT1, HEAT2 ! heat stored in respective canopies
  type(vegn_cohort_type), pointer :: c1, c2
  
  ! calculate normalized weights
  x1 = w1/(w1+w2)
  x2 = 1.0 - x1

  ! the following assumes that there is one, and only one, cohort per tile 
  c1 => t1%cohorts(1)
  c2 => t2%cohorts(1)
  ! define macro for merging cohort values
#define __MERGE__(field) c2%field = x1*c1%field + x2*c2%field
  HEAT1 = (clw*c1%prog%Wl + csw*c1%prog%Ws + c1%mcv_dry)*(c1%prog%Tv-tfreeze)
  HEAT2 = (clw*c2%prog%Wl + csw*c2%prog%Ws + c2%mcv_dry)*(c2%prog%Tv-tfreeze)
  __MERGE__(prog%Wl)
  __MERGE__(prog%Ws)

  __MERGE__(bl)      ! biomass of leaves, kg C/m2
  __MERGE__(blv)     ! biomass of virtual leaves (labile store), kg C/m2
  __MERGE__(br)      ! biomass of fine roots, kg C/m2
  __MERGE__(bsw)     ! biomass of sapwood, kg C/m2
  __MERGE__(bwood)   ! biomass of heartwood, kg C/m2
  __MERGE__(bliving) ! leaves, fine roots, and sapwood biomass, kg C/m2
  __MERGE__(carbon_gain) ! carbon gain during the day, kg C/m2
  __MERGE__(carbon_loss) ! carbon loss during the day, kg C/m2
  __MERGE__(bwood_gain) ! heartwood gain during the day, kg C/m2

  __MERGE__(nl)      ! nitrogen in leaves, kg N/m2
  __MERGE__(nlv)     ! nitrogen in virtual leaves (labile store), kg N/m2
  __MERGE__(nr)      ! nitrogen in fine roots, kg N/m2
  __MERGE__(nsw)     ! nitrogen in sapwood, kg N/m2
  __MERGE__(nwood)   ! nitrogen in heartwood, kg N/m2
  __MERGE__(nliving) ! leaves, fine roots, and sapwood nitrogen
  __MERGE__(nstore)  ! nitrogen in storage, kg N/m2
  __MERGE__(ngain)   ! nitrogen gain between vegn_growth calls, kg N/m2
  __MERGE__(nwoodgain) ! wood nitrogen gain between vegn_growth calls, kg N/m2

  __MERGE__(nsf)     ! measure of the fullness of extra N
  __MERGE__(nfix_per_bliv)

  ! should we do update_derived_vegn_data here? to get mcv_dry, etc
  call update_biomass_pools(t2,c2)

  ! calculate the resulting dry heat capacity
  c2%mcv_dry = max(mcv_min,mcv_lai*c2%lai)
  ! update canopy temperature -- just merge it based on area weights if the heat 
  ! capacities are zero, or merge it based on the heat content if the heat contents
  ! are non-zero
  if(HEAT1==0.and.HEAT2==0) then
     __MERGE__(prog%Tv)
  else
     c2%prog%Tv = (HEAT1*x1+HEAT2*x2) / &
          (clw*c2%prog%Wl + csw*c2%prog%Ws + c2%mcv_dry) + tfreeze
  endif

#undef  __MERGE__
! re-define macro for tile values
#define __MERGE__(field) t2%field = x1*t1%field + x2*t2%field

  __MERGE__(age);
  
  ! *_fsol variables are soluble fractions of nitrogen, so they have to be
  ! merged differently. First, we calculate the amount of soluble nitrogen
  ! in the resulting tile, and then normalize it by total amount of nitrogen
  t2%fast_litter_fsol = x1*t1%fast_litter_N*t1%fast_litter_fsol + x2*t2%fast_litter_N*t2%fast_litter_fsol
  t2%slow_litter_fsol = x1*t1%slow_litter_N*t1%slow_litter_fsol + x2*t2%slow_litter_N*t2%slow_litter_fsol
  t2%slow_soil_fsol   = x1*t1%slow_soil_N*t1%slow_soil_fsol   + x2*t2%slow_soil_N*t2%slow_soil_fsol 

  __MERGE__(fast_soil_C)! there's no nitrogen associated with this pool
  __MERGE__(slow_soil_C); __MERGE__(slow_soil_N)
  __MERGE__(res_soil_C);  __MERGE__(res_soil_N)
  
  __MERGE__(fast_litter_C); __MERGE__(slow_litter_C)
  __MERGE__(fast_litter_N); __MERGE__(slow_litter_N)

  ! convert soluble amounts to soluble fractions
  if (t2%fast_litter_N > 0 ) then
     t2%fast_litter_fsol = t2%fast_litter_fsol/t2%fast_litter_N
  else
     t2%fast_litter_fsol = 0
  endif
  if (t2%slow_litter_N > 0 ) then
     t2%slow_litter_fsol = t2%slow_litter_fsol/t2%slow_litter_N
  else
     t2%slow_litter_fsol = 0
  endif
  if (t2%slow_soil_N > 0 ) then
     t2%slow_soil_fsol = t2%slow_soil_fsol/t2%slow_soil_N
  else
     t2%slow_soil_fsol = 0
  endif
  
  __MERGE__(amm)
  __MERGE__(nitr)
  
  __MERGE__(fsc_pool); __MERGE__(fsc_rate)
  __MERGE__(ssc_pool); __MERGE__(ssc_rate)

  __MERGE__(csmoke_pool)
  __MERGE__(csmoke_rate)

  __MERGE__(harv_pool)
  __MERGE__(harv_rate)

  ! do we need to merge these?
  __MERGE__(asoil_in)
  __MERGE__(ssc_in); __MERGE__(ssc_out)
  __MERGE__(fsc_in); __MERGE__(fsc_out)
  __MERGE__(veg_in); __MERGE__(veg_out)
  
  ! or these?
  __MERGE__(disturbance_rate)
  __MERGE__(lambda)     ! cumulative drought months per year
  __MERGE__(fuel)       ! fuel over dry months
  __MERGE__(litter)     ! litter flux

  ! monthly accumulated/averaged values
  __MERGE__(theta_av)   ! relative soil_moisture availability not soil moisture
  __MERGE__(tsoil_av)   ! bulk soil temperature
  __MERGE__(tc_av)      ! leaf temperature
  __MERGE__(precip_av)  ! precipitation

  ! annual-mean values
  __MERGE__(t_ann)      ! annual mean T, degK
  __MERGE__(t_cold)     ! average temperature of the coldest month, degK
  __MERGE__(p_ann)      ! annual mean precipitation
  __MERGE__(ncm)        ! number of cold months

  ! annual accumulated values
  __MERGE__(t_ann_acm)  ! accumulated annual temperature for t_ann
  __MERGE__(t_cold_acm) ! temperature of the coldest month in current year
  __MERGE__(p_ann_acm)  ! accumulated annual precipitation for p_ann
  __MERGE__(ncm_acm)    ! accumulated number of cold months

#undef __MERGE__

end subroutine merge_vegn_tiles


! ============================================================================
! given a vegetation tile with the state variables set up, calculate derived
! parameters to get a consistent state
! NOTE: this subroutine does not call update_biomass_pools, although some 
! of the calculations are the same. The reason is because this function may 
! be used in the situation when the biomasses are not precisely consistent, for
! example when they come from the data override or from initial conditions.
subroutine update_derived_vegn_data(vegn)
  type(vegn_tile_type), intent(inout) :: vegn
  
  ! ---- local vars 
  type(vegn_cohort_type), pointer :: cc ! pointer to the current cohort
  integer :: i  ! cohort index
  integer :: sp ! shorthand for the vegetation species
  real    :: bs ! structural biomass: stem + structural roots, kg C/m2
  
  ! given that the cohort state variables are initialized, fill in
  ! the intermediate variables
  do i = 1,vegn%n_cohorts
    cc=>vegn%cohorts(i)
    
    sp = cc%species
    ! set the physiology type according to species
    cc%pt     = spdata(sp)%pt
    ! calculate total biomass, calculate height
    cc%b      = cc%bliving + cc%bwood
    cc%height = height_from_biomass(cc%b);
    ! update fractions of the living biomass
    call update_bio_living_fraction(cc)
    bs     = cc%bsw + cc%bwood;   

    if(sp<NSPECIES) then ! LM3V species
       ! calculate the leaf area index based on the biomass of leaves
       cc%lai = lai_from_biomass(cc%bl, sp)
       ! calculate the root density as the total biomass below ground, in
       ! biomass (not carbon!) units
       cc%root_density = (cc%br + (cc%bsw+cc%bwood+cc%blv)*(1-agf_bs))*C2B
    else
       cc%height        = spdata(sp)%dat_height
       cc%lai           = spdata(sp)%dat_lai
       cc%root_density  = spdata(sp)%dat_root_density
    endif
    cc%sai           = 0.035*cc%height
    cc%leaf_size     = spdata(sp)%leaf_size
    cc%root_zeta     = spdata(sp)%dat_root_zeta
    cc%rs_min        = spdata(sp)%dat_rs_min
    cc%leaf_refl     = spdata(sp)%leaf_refl
    cc%leaf_tran     = spdata(sp)%leaf_tran
    cc%leaf_emis     = spdata(sp)%leaf_emis
    cc%snow_crit     = spdata(sp)%dat_snow_crit
  
    ! set nitrogen state
    if (do_nitrogen) then
       cc%nsf = max(cc%nstore,0.0)/n_store_optimum(cc)
    endif
    ! putting this initialization within the cohort loop is probably incorrect 
    ! in case of multiple-cohort vegetation, however for a single cohort it works
    cc%Wl_max   = spdata(sp)%cmc_lai*cc%lai
    cc%Ws_max   = spdata(sp)%csc_lai*cc%lai
    cc%mcv_dry = max(mcv_min, mcv_lai*cc%lai)
  enddo
    
end subroutine update_derived_vegn_data

! ============================================================================
! calculate optimum nitrogen storage
! it has to be here (unfortunately) because it is called from 
! update_derived_vegn_data, and knows about vegn_cohort_type (and therefore
! can't be in nitrogen_data)
function n_store_optimum(cc) result(res)
  type(vegn_cohort_type), intent(in) :: cc
  real :: res ! return value
  
  integer :: sp ! shorthand for cc%species
  
  sp = cc%species

  res = cc%bliving*planning_horizon * ( &
       cc%Pl*spdata(sp)%alpha(CMPT_LEAF)*(1-spdata(sp)%retrans_n)/spdata(sp)%CtoN(CMPT_LEAF) &
     + cc%Pr*spdata(sp)%alpha(CMPT_ROOT)/spdata(sp)%CtoN(CMPT_ROOT) & 
     + cc%Psw_alphasw/spdata(sp)%CtoN(CMPT_WOOD) &
     )
end function

! ============================================================================
! redistribute living biomass pools in a given cohort, and update related 
! properties (height, lai, sai)
subroutine update_biomass_pools(vegn,cc)
  type(vegn_tile_type),   intent(inout) :: vegn
  type(vegn_cohort_type), intent(inout) :: cc

  cc%b      = cc%bliving + cc%bwood;
  cc%height = height_from_biomass(cc%b);
  call update_bio_living_fraction(cc);
  
  if (do_nitrogen) call n_update_pool(vegn, cc);

  cc%bsw = cc%Psw*cc%bliving;
  if(cc%status == LEAF_OFF) then
     cc%blv = cc%Pl*cc%bliving + cc%Pr*cc%bliving;
     cc%bl  = 0;
     cc%br  = 0;
  else
     cc%blv = 0;
     cc%bl  = cc%Pl*cc%bliving;
     cc%br  = cc%Pr*cc%bliving;
  endif
  cc%lai = lai_from_biomass(cc%bl,cc%species)
  cc%sai = 0.035*cc%height ! Federer and Lash,1978
end subroutine 

! ===========================================================================
! this function can and must be called after update_bio_living_fraction  
subroutine n_update_pool(vegn,cc)
  type(vegn_tile_type), intent(inout) :: vegn
  type(vegn_cohort_type), intent(inout) :: cc

  cc%nsw = cc%bliving*cc%Psw/spdata(cc%species)%CtoN(CMPT_SAPWOOD)

  if(cc%status == LEAF_OFF) then
    cc%nlv = cc%nliving - cc%nsw; ! that is nl + nr
    cc%nl = 0.0;
    cc%nr = 0.0;
  else
    cc%nl = cc%bliving*cc%Pl/spdata(cc%species)%CtoN(CMPT_LEAF) ! leaf have fixed C:N
    cc%nr = cc%nliving - cc%nsw - cc%nl;
    cc%nlv= 0.;
  endif
end subroutine n_update_pool

! ============================================================================
! returns the profiles of uptake used in the 'LINEAR' uptake option
subroutine vegn_uptake_profile(vegn, dz, uptake_frac_max, vegn_uptake_term)
  type(vegn_tile_type), intent(in)  :: vegn
  real,                 intent(in)  :: dz(:)
  real,                 intent(out) :: uptake_frac_max(:)
  real,                 intent(out) :: vegn_uptake_term(:)

  call cohort_uptake_profile(vegn%cohorts(1), dz, uptake_frac_max, vegn_uptake_term)
end subroutine


! ============================================================================
subroutine vegn_root_properties (vegn, dz, VRL, K_r, r_r)
  type(vegn_tile_type), intent(in)  :: vegn 
  real,                 intent(in)  :: dz(:)
  real, intent(out) :: &
       vrl(:), & ! volumetric fine root length, m/m3
       K_r,    & ! root membrane permeability per unit area, kg/(m3 s)
       r_r       ! radius of fine roots, m

  call cohort_root_properties(vegn%cohorts(1), dz, VRL, K_r, r_r)
end subroutine 


! ============================================================================
function vegn_data_rs_min ( vegn )
  real :: vegn_data_rs_min
  type(vegn_tile_type), intent(in)  :: vegn
  
  vegn_data_rs_min = vegn%cohorts(1)%rs_min
end function


! ============================================================================
function vegn_seed_supply ( vegn )
  real :: vegn_seed_supply
  type(vegn_tile_type), intent(in) :: vegn

  ! ---- local vars 
  real :: vegn_bliving
  integer :: i
  
  vegn_bliving = 0
  do i = 1,vegn%n_cohorts
     vegn_bliving = vegn_bliving + vegn%cohorts(i)%bliving
  enddo
  vegn_seed_supply = MAX (vegn_bliving-BSEED, 0.0)
  
end function 

! ============================================================================
function vegn_seed_demand ( vegn )
  real :: vegn_seed_demand
  type(vegn_tile_type), intent(in) :: vegn

  integer :: i

  vegn_seed_demand = 0
  do i = 1,vegn%n_cohorts
     if(vegn%cohorts(i)%bliving<BSEED.and.vegn%t_ann>253.16.and.vegn%p_ann>1E-6) then
        vegn_seed_demand = vegn_seed_demand + BSEED
     endif
  enddo
end function 

! ============================================================================
subroutine vegn_add_bliving ( vegn, delta )
  type(vegn_tile_type), intent(inout) :: vegn
  real :: delta ! increment of bliving
  
  character(512) :: message
  integer :: i,j,k,face

  vegn%cohorts(1)%bliving = vegn%cohorts(1)%bliving + delta
  if (vegn%cohorts(1)%bliving < 0) then
     call get_current_point(i,j,k,face)
     write(message,'(a,4(x,a,i4), 2(x,a,g))')&
        'resulting bliving is negative', &
        'at i=',i,'j=',j,'k=',k,'face=',face, &
        'delta=',delta,'result=',vegn%cohorts(1)%bliving
     call error_mesg('vegn_add_bliving', message, FATAL)
  endif
  
  if (do_nitrogen) then
     vegn%cohorts(1)%nliving = vegn%cohorts(1)%nliving + delta*CtoN_seed
     call get_current_point(i,j,k,face)
     if (vegn%cohorts(1)%nliving < 0) then
        write(message,'(a,4(x,a,i4), 2(x,a,g))')&
             'resulting nliving is negative', &
             'at i=',i,'j=',j,'k=',k,'face=',face, &
             'delta=',delta*CtoN_seed,'result=',vegn%cohorts(1)%nliving
        call error_mesg('vegn_add_bliving', message, WARNING)
     endif
  endif
  
  call update_biomass_pools(vegn,vegn%cohorts(1))
end subroutine 





! ============================================================================
! given a vegetation patch, destination kind of transition, and "transition 
! intensity" value, this function returns a fraction of tile that will parti-
! cipate in transition.
!
! this function must be contiguous, monotonic, its value must be within
! interval [0,1]
!
! this function is used to determine what part of each tile is to be converted
! to another land use kind; the equation is solved to get "transition intensity" 
! tau for which total area is equal to requested. Tau is, therefore, a dummy
! parameter, and only relative values of the priority functions for tiles 
! participating in transition have any meaning. For most transitions the priority 
! function is just equal to tau: therefore there is no preference, and all tiles
! contribute equally to converted area. For secondary vegetation harvesting, 
! however, priority also depends on wood biomass, and therefore tiles
! with high wood biomass are harvested first.
function vegn_tran_priority(vegn, dst_kind, tau) result(pri)
  real :: pri
  type(vegn_tile_type), intent(in) :: vegn
  integer             , intent(in) :: dst_kind
  real                , intent(in) :: tau

  real :: vegn_bwood
  integer :: i

  if (vegn%landuse==LU_SCND.and.dst_kind==LU_SCND) then ! secondary biomass harvesting
     vegn_bwood = 0
     do i = 1,vegn%n_cohorts
        vegn_bwood = vegn_bwood + vegn%cohorts(i)%bwood
     enddo
     pri = max(min(tau+vegn_bwood,1.0),0.0)
  else
     pri = max(min(tau,1.0),0.0)
  endif
end function 


! ============================================================================
function vegn_cover_cold_start(land_mask, glonb, glatb) result (vegn_frac)
! creates and initializes a field of fractional vegn coverage
! should be called for global grid; otherwise the part that fills
! missing data points may fail to find any good data, or do it in nproc-
! dependent way
  logical, intent(in) :: land_mask(:,:)    ! global land mask
  real,    intent(in) :: glonb(:,:), glatb(:,:)! boundaries of the global grid cells
  real,    pointer    :: vegn_frac (:,:,:) ! output-global map of vegn fractional coverage

  allocate( vegn_frac(size(land_mask,1),size(land_mask,2),MSPECIES))

  call init_cover_field(vegn_to_use, 'INPUT/cover_type.nc', 'cover','frac', &
       glonb, glatb, vegn_index_constant, input_cover_types, vegn_frac)
  
end function 

! =============================================================================
! returns true if tile fits the specified selector
function vegn_is_selected(vegn, sel)
  logical vegn_is_selected
  type(tile_selector_type),  intent(in) :: sel
  type(vegn_tile_type),      intent(in) :: vegn

  select case (sel%idata1)
  case (LU_SEL_TAG)
     vegn_is_selected = (sel%idata2 == vegn%landuse)
  case (SP_SEL_TAG)
     if (.not.associated(vegn%cohorts)) then
        vegn_is_selected = .FALSE.
     else
        vegn_is_selected = (sel%idata2 == vegn%cohorts(1)%species)
     endif
  case (NG_SEL_TAG)
     if (.not.associated(vegn%cohorts)) then
        vegn_is_selected = .FALSE.
     else
        vegn_is_selected = &
             ((vegn%cohorts(1)%species==SP_C4GRASS) .or.&
              (vegn%cohorts(1)%species==SP_C3GRASS)).and.&
             ((vegn%landuse==LU_NTRL).or. &
              (vegn%landuse==LU_SCND))
     endif
  case default
     vegn_is_selected = .FALSE.
  end select  
     
end function


! ============================================================================
! returns tag of the tile
function get_vegn_tile_tag(vegn) result(tag)
  integer :: tag
  type(vegn_tile_type), intent(in) :: vegn
  
  tag = vegn%tag
end function

! ============================================================================
! returns total wood biomass per tile 
function get_vegn_tile_bwood(vegn) result(bwood)
  real :: bwood
  type(vegn_tile_type), intent(in) :: vegn

  ! ---- local vars
  integer :: i

  bwood = 0
  do i = 1,vegn%n_cohorts
     bwood = bwood + vegn%cohorts(i)%bwood
  enddo
end function

! ============================================================================
! returns total living biomass per tile 
function vegn_tile_bliving(vegn) result(bliving)
  real :: bliving
  type(vegn_tile_type), intent(in) :: vegn

  ! ---- local vars
  integer :: i

  bliving = 0
  do i = 1,vegn%n_cohorts
     bliving = bliving + vegn%cohorts(i)%bliving
  enddo
end function

! ============================================================================
! returns total living nitrogen per tile 
function vegn_tile_nliving(vegn) result(nliving)
  real :: nliving
  type(vegn_tile_type), intent(in) :: vegn

  ! ---- local vars
  integer :: i

  nliving = 0
  do i = 1,vegn%n_cohorts
     nliving = nliving + vegn%cohorts(i)%nliving
  enddo
end function

! ============================================================================
subroutine vegn_tile_stock_pe (vegn, twd_liq, twd_sol  )
  type(vegn_tile_type),  intent(in)    :: vegn
  real,                  intent(out)   :: twd_liq, twd_sol
  integer n
  
  twd_liq = 0.
  twd_sol = 0.
  do n=1, vegn%n_cohorts
    twd_liq = twd_liq + vegn%cohorts(n)%prog%wl
    twd_sol = twd_sol + vegn%cohorts(n)%prog%ws
!      vegn_HEAT  = (mcv + clw*cohort%prog%Wl+ csw*cohort%prog%Ws)*(cohort%prog%Tv-tfreeze)

    enddo
end subroutine vegn_tile_stock_pe


! ============================================================================
! returns total carbon in the tile, kg C/m2
function vegn_tile_carbon(vegn) result(carbon) ; real carbon
  type(vegn_tile_type), intent(in)  :: vegn

  integer :: i

  carbon = 0
  do i = 1,vegn%n_cohorts
     carbon = carbon + &
          vegn%cohorts(i)%bl + vegn%cohorts(i)%blv + &
          vegn%cohorts(i)%br + vegn%cohorts(i)%bwood + &
          vegn%cohorts(i)%bsw + &
          vegn%cohorts(i)%carbon_gain + vegn%cohorts(i)%bwood_gain
          ! NOTE that carbon_loss is not included in total carbon calculations
          ! since it is (currently) just a diagnostics variable
  enddo
  carbon = carbon &
       + sum(vegn%harv_pool(:,isub_C)) &
       + vegn%fsc_pool(isub_C) + vegn%ssc_pool(isub_C) &
       + vegn%csmoke_pool(isub_C) &
       + vegn%fast_litter_C + vegn%slow_litter_C &
       + vegn%fast_soil_C + vegn%slow_soil_C + vegn%res_soil_C
end function


! ============================================================================
! returns total nitrogen in tile, kg N/m2
function vegn_tile_nitrogen(vegn) result(nitrogen) ; real nitrogen
  type(vegn_tile_type), intent(in)  :: vegn

  integer :: i
  nitrogen = 0
  do i = 1,vegn%n_cohorts
    nitrogen = nitrogen &
       + vegn%cohorts(i)%nliving &
       + vegn%cohorts(i)%nwood &
       + vegn%cohorts(i)%nstore &
       + vegn%cohorts(i)%ngain &
       + vegn%cohorts(i)%nwoodgain
  enddo
  if(isub_N > 0) nitrogen = nitrogen + &
       + sum(vegn%harv_pool(:,isub_N)) &
       + vegn%csmoke_pool(isub_N) &
       + vegn%fsc_pool(isub_N) + vegn%ssc_pool(isub_N) 
  nitrogen = nitrogen + &
       + vegn%fast_litter_N + vegn%slow_litter_N &
       + vegn%slow_soil_N + vegn%res_soil_N &
       + vegn%amm + vegn%nitr
end function

! ============================================================================
! returns heat content of the vegetation, J/m2
function vegn_tile_heat (vegn) result(heat) ; real heat
  type(vegn_tile_type), intent(in)  :: vegn

  integer :: i

  heat = 0
  do i = 1, vegn%n_cohorts
     heat = heat + &
          (clw*vegn%cohorts(i)%prog%Wl + &
             csw*vegn%cohorts(i)%prog%Ws + &
             vegn%cohorts(i)%mcv_dry)*(vegn%cohorts(i)%prog%Tv-tfreeze) - &
           hlf*vegn%cohorts(i)%prog%Ws
  enddo
end function

end module vegn_tile_mod
