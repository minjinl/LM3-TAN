module land_substances_mod

use fms_mod, only : error_mesg, FATAL

implicit none
private

! ==== public interfaces =====================================================
public :: land_substances_init ! initialize module
public :: land_substances_end ! clean-up

public :: n_substances  ! the number of registered substances
public :: add_substance ! adds a substance to the list
public :: inq_substance ! returns information about substance

interface inq_substance
   module procedure inq_subst_by_name
   module procedure inq_subst_by_idx
end interface
! ==== end of public interfaces ==============================================

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'land_substances_mod', &
     version     = '$Id: land_substances.F90,v 1.1.2.2 2010/04/13 15:14:22 slm Exp $', &
     tagname     = '$Name: nitro_20111003_slm $'

! ==== derived types =========================================================
type :: substance_type
  character(128) :: name
  character(128) :: long_name
end type

! ==== module variables ======================================================
logical :: module_is_initialized = .FALSE.
integer :: n_substances
type(substance_type), pointer :: subst(:)

contains


! ============================================================================
subroutine land_substances_init()
   allocate(subst(10)) ! allocate storage for substance info (arbitrary initial size)
   n_substances = 0    ! initialize substance counter
   module_is_initialized = .TRUE.
end subroutine


! ============================================================================
subroutine land_substances_end()
   deallocate(subst)
   n_substances = 0
   module_is_initialized = .FALSE.
end subroutine


! ============================================================================
subroutine add_substance(name, long_name)
   character(*), intent(in) :: name      ! name of the substance
   character(*), intent(in) :: long_name ! long name of the substance
   
   integer :: i
   
   type(substance_type), pointer :: new_subst(:)
   
   if (.not.module_is_initialized) call land_substances_init()
   
   ! check that the new substance name is unique
   do i = 1,n_substances
      if (trim(name)==trim(subst(i)%name)) &
         call error_mesg(module_name, &
                        'substance "'//trim(name)//'" already exists', FATAL)
   enddo
   
   ! allocate more space for substances, if necessary
   if (n_substances>=size(subst)) then
      allocate(new_subst(max(2*n_substances,1)))
      if (associated(subst)) then
         new_subst(1:n_substances) = subst(1:n_substances)
         deallocate(subst)
      endif
      subst=>new_subst
   endif
   
   n_substances = n_substances+1
   subst(n_substances)%name      = name
   subst(n_substances)%long_name = long_name
end subroutine


! ============================================================================
subroutine inq_subst_by_name(name, idx, long_name)
   character(*), intent(in) :: name
   integer, intent(out), optional :: idx
   character(*), intent(out), optional :: long_name
   
   integer :: i

   if (.not.module_is_initialized) call land_substances_init()

   do i = 1, n_substances
      if(subst(i)%name==trim(name)) then
         if (present(idx))       idx = i
         if (present(long_name)) long_name = subst(i)%long_name
         return
      endif
   enddo
   ! we get here only if the name is not found
   if (present(idx))       idx = -1
   if (present(long_name)) long_name = ' '
end subroutine


! ============================================================================
subroutine inq_subst_by_idx(idx, name, long_name)
   integer, intent(in) :: idx
   character(*), intent(out), optional :: name
   character(*), intent(out), optional :: long_name
   
   if (.not.module_is_initialized) call land_substances_init()

   if (idx>0.and.idx<=n_substances) then
      if (present(name)) name=subst(idx)%name
      if (present(long_name)) long_name=subst(idx)%long_name
   else
      if (present(name)) name=' '
      if (present(long_name)) long_name=' '
   endif
end subroutine

end module land_substances_mod