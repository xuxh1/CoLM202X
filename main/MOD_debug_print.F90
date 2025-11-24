#include "define.h"

module mod_debug_print
   use MOD_Precision, only: r8
   implicit none
   private

   !==============================
   ! Module persistent variables
   !==============================
   integer :: call_count = 0
   integer :: last_idate(3)
   logical :: first_call = .true.

   ! Proper initialization
   data last_idate / -999, -999, -999 /

   public :: print_var

   interface print_var
      module procedure print_var_real_scalar
      module procedure print_var_real_1d
      module procedure print_var_real_2d

      module procedure print_var_int_scalar
      module procedure print_var_int_1d
      module procedure print_var_int_2d

      module procedure print_var_logical_scalar
      module procedure print_var_logical_1d
      module procedure print_var_logical_2d

      module procedure print_var_char_scalar
      module procedure print_var_char_1d
      module procedure print_var_char_2d
   end interface print_var

contains
!=====================================================================
! Unified counter management
!=====================================================================
   subroutine print_var_entry(varname, idate)
      character(len=*), intent(in) :: varname
      integer, intent(in) :: idate(3)

      ! First call ever
      if (first_call) then
         call_count = 0
         last_idate = idate
         first_call = .false.
      endif

      ! If idate changed, reset counter
      if (any(idate /= last_idate)) then
         call_count = 0
         last_idate = idate
      endif

      ! Increase counter
      call_count = call_count + 1
   end subroutine print_var_entry


!=====================================================================
! Print control logic
!=====================================================================
   logical function should_print(varname, idate, j, layer_ok)
      character(len=*), intent(in)  :: varname
      integer, intent(in)           :: idate(3)
      integer, intent(in), optional :: j
      logical, intent(in), optional :: layer_ok

      should_print = .false.

      ! Force print if specific variable name
      ! if (trim(varname) == "ch4 errch4") then
      if (index(trim(varname), "SoilBiogeochemCompetition") > 0) then
         should_print = .true.
         return
      endif

      ! Only print at specific DOY
      ! if (idate(2) /= 92) return
      return

      ! Restrict by layer index
      if (present(j)) then
         if (j == 0 .or. j == 1 .or. j == 10) then
            should_print = .true.
         else
            should_print = .false.
         endif
      else
         should_print = .true.
      endif

      ! Force layer print
      if (present(layer_ok)) then
         if (layer_ok) should_print = .true.
      endif

   end function should_print


!=====================================================================
! Header and footer
!=====================================================================
   subroutine print_header(varname, idate, j, s)
      character(len=*), intent(in)           :: varname
      integer,          intent(in)           :: idate(3)
      integer,          intent(in), optional :: j, s

      integer :: year, doy, sod
      integer :: month, day, hour, minute, second
      integer :: mdays(12)

      year = idate(1)
      doy  = idate(2)
      sod  = idate(3)

      ! Determine days per month
      if (mod(year,400) == 0 .or. (mod(year,4) == 0 .and. mod(year,100) /= 0)) then
         mdays = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      else
         mdays = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      endif

      ! Convert DOY to month/day
      month = 1
      do while (doy > mdays(month))
         doy = doy - mdays(month)
         month = month + 1
      end do
      day = doy

      hour   = sod / 3600
      minute = mod(sod,3600) / 60
      second = mod(sod,60)

      print *, "===================================================================================="
      print '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2,"  ",I6.6)', year, month, day, hour, minute, second, call_count
      if (present(j)) print *, "layer =", j
      if (present(s)) print *, "ngas  =", s
      print *, trim(varname), " ="
   end subroutine print_header


   subroutine print_footer()
      print *, "===================================================================================="
   end subroutine print_footer


!=====================================================================
! REAL(r8) PRINTING
!=====================================================================
   subroutine print_var_real_scalar(var, varname, idate, j, layer_ok, s)
      real(r8), intent(in) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in) :: idate(3)
      integer, optional, intent(in) :: j, s
      logical, optional, intent(in) :: layer_ok

#ifdef SinglePoint
      if (.not. should_print(varname,idate,j,layer_ok)) return
      call print_var_entry(varname, idate)

      call print_header(varname,idate,j,s)
      write(*,'(A,ES24.16E3)') '                ', var
      call print_footer()
#endif
   end subroutine print_var_real_scalar


   subroutine print_var_real_1d(var, varname, idate, j, layer_ok, s)
      real(r8), intent(in) :: var(:)
      character(len=*), intent(in) :: varname
      integer, intent(in) :: idate(3)
      integer, optional, intent(in) :: j, s
      logical, optional, intent(in) :: layer_ok

      integer :: i, n
#ifdef SinglePoint
      if (.not. should_print(varname,idate,j,layer_ok)) return
      call print_var_entry(varname, idate)

      n = size(var)

      call print_header(varname,idate,j,s)
      if (n == 15) then
         do i = 1, n
            write(*,'(A,I0,A,ES14.6E3)') '  (', i-5, ') = ', var(i)
         end do
      else
         do i = 1, n
            write(*,'(A,I0,A,ES14.6E3)') '  (', i, ') = ', var(i)
         end do
      endif
      call print_footer()
#endif
   end subroutine print_var_real_1d


   subroutine print_var_real_2d(var, varname, idate, j, layer_ok, s)
      real(r8), intent(in) :: var(:,:)
      character(len=*), intent(in) :: varname
      integer, intent(in) :: idate(3)
      integer, optional, intent(in) :: j, s
      logical, optional, intent(in) :: layer_ok

      integer :: i, k, n1, n2
#ifdef SinglePoint
      if (.not. should_print(varname,idate,j,layer_ok)) return
      call print_var_entry(varname, idate)

      n1 = size(var,1)
      n2 = size(var,2)

      call print_header(varname,idate,j,s)
      do k = 1, n2
         write(*,'(A,I0)') '  row ', k
         do i = 1, n1
            write(*,'(A,I0,A,I0,A,ES14.6E3)') '    (', i, ',', k, ') = ', var(i,k)
         end do
      end do
      call print_footer()
#endif
   end subroutine print_var_real_2d


!=====================================================================
! INTEGER PRINTING
!=====================================================================
   subroutine print_var_int_scalar(var, varname, idate, j, layer_ok, s)
      integer, intent(in) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in) :: idate(3)
      integer, optional, intent(in) :: j, s
      logical, optional, intent(in) :: layer_ok

#ifdef SinglePoint
      if (.not. should_print(varname,idate,j,layer_ok)) return
      call print_var_entry(varname, idate)

      call print_header(varname,idate,j,s)
      write(*,'(A,I0)') '                ', var
      call print_footer()
#endif
   end subroutine print_var_int_scalar


   subroutine print_var_int_1d(var, varname, idate, j, layer_ok, s)
      integer, intent(in) :: var(:)
      character(len=*), intent(in) :: varname
      integer, intent(in) :: idate(3)
      integer, optional, intent(in) :: j, s
      logical, optional, intent(in) :: layer_ok

      integer :: i, n
#ifdef SinglePoint
      if (.not. should_print(varname,idate,j,layer_ok)) return
      call print_var_entry(varname, idate)

      n = size(var)

      call print_header(varname,idate,j,s)
      if (n == 15) then
         do i = 1, n
            write(*,'(A,I0,A,I0)') '  (', i-5, ') = ', var(i)
         end do
      else
         do i = 1, n
            write(*,'(A,I0,A,I0)') '  (', i, ') = ', var(i)
         end do
      endif

      call print_footer()
#endif
   end subroutine print_var_int_1d


   subroutine print_var_int_2d(var, varname, idate, j, layer_ok, s)
      integer, intent(in) :: var(:,:)
      character(len=*), intent(in) :: varname
      integer, intent(in) :: idate(3)
      integer, optional, intent(in) :: j, s
      logical, optional, intent(in) :: layer_ok

      integer :: i, k, n1, n2
#ifdef SinglePoint
      if (.not. should_print(varname,idate,j,layer_ok)) return
      call print_var_entry(varname, idate)

      n1 = size(var,1)
      n2 = size(var,2)

      call print_header(varname,idate,j,s)
      do k = 1, n2
         write(*,'(A,I0)') '  row ', k
         do i = 1, n1
            write(*,'(A,I0,A,I0,A,I0)') '    (', i, ',', k, ') = ', var(i,k)
         end do
      end do
      call print_footer()
#endif
   end subroutine print_var_int_2d


!=====================================================================
! LOGICAL PRINTING
!=====================================================================
   subroutine print_var_logical_scalar(var, varname, idate, j, layer_ok, s)
      logical, intent(in) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in) :: idate(3)
      integer, optional, intent(in) :: j, s
      logical, optional, intent(in) :: layer_ok

#ifdef SinglePoint
      if (.not. should_print(varname,idate,j,layer_ok)) return
      call print_var_entry(varname, idate)

      call print_header(varname,idate,j,s)
      write(*,'(A,L1)') '                ', var
      call print_footer()
#endif
   end subroutine print_var_logical_scalar


   subroutine print_var_logical_1d(var, varname, idate, j, layer_ok, s)
      logical, intent(in) :: var(:)
      character(len=*), intent(in) :: varname
      integer, intent(in) :: idate(3)
      integer, optional, intent(in) :: j, s
      logical, optional, intent(in) :: layer_ok

      integer :: i, n
#ifdef SinglePoint
      if (.not. should_print(varname,idate,j,layer_ok)) return
      call print_var_entry(varname, idate)

      n = size(var)

      call print_header(varname,idate,j,s)
      if (n == 15) then
         do i = 1, n
            write(*,'(A,I0,A,L1)') '  (', i-5, ') = ', var(i)
         end do
      else
         do i = 1, n
            write(*,'(A,I0,A,L1)') '  (', i, ') = ', var(i)
         end do
      endif

      call print_footer()
#endif
   end subroutine print_var_logical_1d


   subroutine print_var_logical_2d(var, varname, idate, j, layer_ok, s)
      logical, intent(in) :: var(:,:)
      character(len=*), intent(in) :: varname
      integer, intent(in) :: idate(3)
      integer, optional, intent(in) :: j, s
      logical, optional, intent(in) :: layer_ok

      integer :: i, k, n1, n2
#ifdef SinglePoint
      if (.not. should_print(varname,idate,j,layer_ok)) return
      call print_var_entry(varname, idate)

      n1 = size(var,1)
      n2 = size(var,2)

      call print_header(varname,idate,j,s)
      do k = 1, n2
         write(*,'(A,I0)') '  row ', k
         do i = 1, n1
            write(*,'(A,I0,A,I0,A,L1)') '    (', i, ',', k, ') = ', var(i,k)
         end do
      end do

      call print_footer()
#endif
   end subroutine print_var_logical_2d


!=====================================================================
! CHARACTER PRINTING
!=====================================================================
   subroutine print_var_char_scalar(var, varname, idate, j, layer_ok, s)
      character(len=*), intent(in) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in) :: idate(3)
      integer, optional, intent(in) :: j, s
      logical, optional, intent(in) :: layer_ok

#ifdef SinglePoint
      if (.not. should_print(varname,idate,j,layer_ok)) return
      call print_var_entry(varname, idate)

      call print_header(varname,idate,j,s)
      write(*,'(A,A)') '                ', trim(var)
      call print_footer()
#endif
   end subroutine print_var_char_scalar


   subroutine print_var_char_1d(var, varname, idate, j, layer_ok, s)
      character(len=*), intent(in) :: var(:)
      character(len=*), intent(in) :: varname
      integer, intent(in) :: idate(3)
      integer, optional, intent(in) :: j, s
      logical, optional, intent(in) :: layer_ok

      integer :: i, n
#ifdef SinglePoint
      if (.not. should_print(varname,idate,j,layer_ok)) return
      call print_var_entry(varname, idate)

      call print_header(varname,idate,j,s)
      n = size(var)
      do i = 1, n
         write(*,'(A,I0,A,A)') '  (', i, ') = ', trim(var(i))
      end do
      call print_footer()
#endif
   end subroutine print_var_char_1d


   subroutine print_var_char_2d(var, varname, idate, j, layer_ok, s)
      character(len=*), intent(in) :: var(:,:)
      character(len=*), intent(in) :: varname
      integer, intent(in) :: idate(3)
      integer, optional, intent(in) :: j, s
      logical, optional, intent(in) :: layer_ok

      integer :: i, k, n1, n2
#ifdef SinglePoint
      if (.not. should_print(varname,idate,j,layer_ok)) return
      call print_var_entry(varname, idate)

      n1 = size(var,1)
      n2 = size(var,2)

      call print_header(varname,idate,j,s)
      do k = 1, n2
         write(*,'(A,I0)') '  row ', k
         do i = 1, n1
            write(*,'(A,I0,A,I0,A,A)') '    (', i, ',', k, ') = ', trim(var(i,k))
         end do
      end do
      call print_footer()
#endif
   end subroutine print_var_char_2d


end module mod_debug_print
