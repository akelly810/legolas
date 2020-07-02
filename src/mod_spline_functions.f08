! =============================================================================
!> @brief   Module that calculates the finite element basis functions.
!! @details The different basis functions and their derivatives
!!          used throughout the code are calculated here. All routines defined in
!!          this module simply return the basis functions for a specific point /p r
!!          in the interval <tt>(rj_lo, rj_hi)</tt>.
module mod_spline_functions
  use mod_global_variables, only: dp
  use mod_check_values, only: check_small_values
  implicit none

  private

  public :: quadratic_factors
  public :: quadratic_factors_deriv
  public :: cubic_factors
  public :: cubic_factors_deriv

contains


  !> @brief Calculates the quadratic basis functions.
  !! @param[in] r     current value for r in the interval
  !! @param[in] rj_lo left edge of the grid interval
  !! @param[in] rj_hi right edge of the grid interval
  !! @param[out] h_quadratic array containing the four quadratic basis functions in this interval
  subroutine quadratic_factors(r, rj_lo, rj_hi, h_quadratic)
    real(dp), intent(in)  ::  r, rj_lo, rj_hi
    real(dp), intent(out) ::  h_quadratic(4)

    h_quadratic(1) = 4.0d0 * (r - rj_lo) * (rj_hi - r) / (rj_hi - rj_lo)**2
    h_quadratic(2) = 0.0d0
    h_quadratic(3) = (2.0d0*r - rj_hi - rj_lo) * (r - rj_lo) / (rj_hi - rj_lo)**2
    h_quadratic(4) = (2.0d0*r - rj_hi - rj_lo) * (r - rj_hi) / (rj_hi - rj_lo)**2

    call check_small_values(h_quadratic)
  end subroutine quadratic_factors


  !> @brief Calculates the derivative of the quadratic basis functions.
  !! @param[in] r     current value for r in the interval
  !! @param[in] rj_lo left edge of the grid interval
  !! @param[in] rj_hi right edge of the grid interval
  !! @param[out] dh_quadratic_dr array containing the four derivatives of the
  !!                             quadratic basis functions in this interval
  subroutine quadratic_factors_deriv(r, rj_lo, rj_hi, dh_quadratic_dr)
    real(dp), intent(in)  ::  r, rj_lo, rj_hi
    real(dp), intent(out) ::  dh_quadratic_dr(4)

    dh_quadratic_dr(1) = 4.0d0 * (-2.0d0*r + rj_hi + rj_lo) / (rj_hi - rj_lo)**2
    dh_quadratic_dr(2) = 0.0d0
    dh_quadratic_dr(3) = (4.0d0*r - rj_hi - 3.0d0*rj_lo) / (rj_hi - rj_lo)**2
    dh_quadratic_dr(4) = (4.0d0*r - rj_lo - 3.0d0*rj_hi) / (rj_hi - rj_lo)**2

    call check_small_values(dh_quadratic_dr)
  end subroutine quadratic_factors_deriv


  !> @brief Calculates the cubic basis functions.
  !! @param[in] r     current value for r in the interval
  !! @param[in] rj_lo left edge of the grid interval
  !! @param[in] rj_hi right edge of the grid interval
  !! @param[out] h_cubic array containing the four cubic basis functions in this interval
  subroutine cubic_factors(r, rj_lo, rj_hi, h_cubic)
    real(dp), intent(in)  :: r, rj_lo, rj_hi
    real(dp), intent(out) :: h_cubic(4)

    h_cubic(1) =  3.0d0 * ( (r - rj_lo) / (rj_hi - rj_lo) )**2 &
                 -2.0d0 * ( (r - rj_lo) / (rj_hi - rj_lo) )**3
    h_cubic(2) =  3.0d0 * ( (rj_hi - r) / (rj_hi - rj_lo) )**2 &
                 -2.0d0 * ( (rj_hi - r) / (rj_hi - rj_lo) )**3
    h_cubic(3) = (r - rj_hi) * ( (r - rj_lo) / (rj_hi - rj_lo) )**2
    h_cubic(4) = (r - rj_lo) * ( (rj_hi - r) / (rj_hi - rj_lo) )**2

    call check_small_values(h_cubic)
  end subroutine cubic_factors


  !> @brief Calculates the derivative of the cubic basis functions.
  !! @param[in] r     current value for r in the interval
  !! @param[in] rj_lo left edge of the grid interval
  !! @param[in] rj_hi right edge of the grid interval
  !! @param[out] dh_cubic_dr array containing the four derivatives of the
  !!                         cubic basis functions in this interval
  subroutine cubic_factors_deriv(r, rj_lo, rj_hi, dh_cubic_dr)
    real(dp), intent(in)  :: r, rj_lo, rj_hi
    real(dp), intent(out) :: dh_cubic_dr(4)

    dh_cubic_dr(1) =  6.0d0 * (r - rj_lo) / (rj_hi - rj_lo)**2 &
                     -6.0d0 * (r - rj_lo)**2 / (rj_hi - rj_lo)**3
    dh_cubic_dr(2) = -6.0d0 * (rj_hi - r) / (rj_hi - rj_lo)**2 &
                     +6.0d0 * (rj_hi - r)**2 / (rj_hi - rj_lo)**3
    dh_cubic_dr(3) = ( 2.0d0*(r - rj_hi) * (r - rj_lo) + (r - rj_lo)**2 ) &
                     / (rj_hi - rj_lo)**2
    dh_cubic_dr(4) = ( 2.0d0*(r - rj_lo) * (r - rj_hi) + (r - rj_hi)**2 ) &
                     / (rj_hi - rj_lo)**2

    call check_small_values(dh_cubic_dr)
  end subroutine cubic_factors_deriv

end module mod_spline_functions
