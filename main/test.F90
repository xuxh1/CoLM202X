#include "define.h"

module MOD_calculation
	use MOD_Precision
	use MOD_Vars_Global, only : nl_soil
	implicit none
	save

	public  :: get_jwt_finundated

contains

   subroutine get_jwt_finundated (zwt,zi_soisno,jwt,sat,finundated)
      real(r8), intent(in) :: &
         zwt                            , &! the depth from ground (soil) surface to water table [m]
         zi_soisno(maxsnl:nl_soil)         ! interface level below a "z" level (m)

      integer, intent(out)  :: jwt                     ! index of the soil layer right above the water table (-)
      integer, intent(inout)  :: &
         sat                   , & ! 0 = unsatured, 1 = saturated
         finundated                ! fractional inundated area, = sat(0 or 1)

      integer  :: j

      ! compute jwt index
      ! The layer index of the first unsaturated layer,
      ! i.e., the layer right above the water table
      jwt = nl_soil
      ! allow jwt to equal zero when zwt is in top layer
      do j = 1, nl_soil
         if(zwt <= zi_soisno(j)) then
            jwt = j-1
            exit
         end if
      end do

      if (jwt==0) then
         sat = 1
      else
         sat = 0
      end if

      finundated = sat
      ! jwt = 0
