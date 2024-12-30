! =============================================================================
!> Initial Condition
!!

module mod_initial_condition
  use mod_global_variables, only: dp
  use mod_physical_constants, only: dpi
  use mod_logging, only: logger, str
  use mod_settings, only: settings_t
  use mod_state_vector_component, only: sv_component_t

  implicit none

  ! Interface for state vector component profiles
  interface
    function profile_fcn(x) result(res)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: res(size(x))
    end function profile_fcn
  end interface

contains
  subroutine generate_ic(N, nodes, x0)
    !> The settings object is needed to get the basis functions
    ! type(settings_t), intent(in) :: settings
    !> Number of grid points
    integer, intent(in) :: N
    !> Array of grid points
    real(dp), intent(in) :: nodes(N)
    !> Assembled IC array
    real(dp), intent(out) :: x0(6*N)

    type(sv_component_t), allocatable :: components(:)
    integer :: i  ! loop index

    real(dp) :: rho_1(N)
    real(dp) :: rho_2(N)
    real(dp) :: v1_1(N)
    real(dp) :: v1_2(N)
    real(dp) :: T_1(N)
    real(dp) :: T_2(N)

    ! 1. Get state vector
    ! check cubic/quadratic
    ! get basis fcns

    ! 2. Compute expansion coefficients
    
    ! (1) rho
    call compute_quad_coeffs(rho, N, nodes, rho_1, rho_2)

    ! (2) v1
    call compute_cubic_coeffs(v1, dv1, N, nodes, v1_1, v1_2)
    
    ! (3) T
    call compute_quad_coeffs(T, N, nodes, T_1, T_2)


    ! 3. Assemble vector in FEM form

    ! Assembly of interleaved array
    do i = 1, N
      x0(6*(i-1) + 1) = rho_1(i)
      x0(6*(i-1) + 2) = rho_2(i)
      x0(6*(i-1) + 3) = v1_1(i)
      x0(6*(i-1) + 4) = v1_2(i)
      x0(6*(i-1) + 5) = T_1(i)
      x0(6*(i-1) + 6) = T_2(i)
    end do

  end subroutine generate_ic


  ! -----------------------------------------------------------------
  ! Reconstruct profiles on grid from FEM representation
  ! -----------------------------------------------------------------

  subroutine reconstruct_profiles(N, nodes, x, rho, v, T)
    integer, intent(in) :: N
    real(dp), intent(in) :: nodes(N)
    real(dp), intent(in) :: x(6*N)
    real(dp), intent(out) :: rho(N)
    real(dp), intent(out) :: v(N)
    real(dp), intent(out) :: T(N)

    ! do stuff

  end subroutine reconstruct_profiles

  ! -----------------------------------------------------------------
  ! Compute expansion coefficients
  ! -----------------------------------------------------------------
  
  ! rho, T
  subroutine compute_quad_coeffs(fcn, N, nodes, u1, u2)
    procedure(profile_fcn) :: fcn
    integer, intent(in)    :: N
    real(dp), intent(in)   :: nodes(N)
    real(dp), intent(out)  :: u1(N)
    real(dp), intent(out)  :: u2(N)

    integer :: i
    real(dp), allocatable :: midpoints(:)
    allocate(midpoints(N - 1))

    ! u1 = fcn(midpoints)
    midpoints = 0.5 * (nodes(2:) + nodes(:N-1))
    u1 = [(0.0d0, i = 1,N)]
    u1(2:) = fcn(midpoints)  ! 0 followed by N-1 midpoint values

    ! u2 = fcn(nodes)
    u2 = fcn(nodes)
  end subroutine compute_quad_coeffs


  ! v1
  subroutine compute_cubic_coeffs(fcn, dfcn, N, nodes, u1, u2)
    procedure(profile_fcn) :: fcn
    procedure(profile_fcn) :: dfcn
    integer, intent(in)    :: N
    real(dp), intent(in)   :: nodes(N)
    real(dp), intent(out)  :: u1(N)
    real(dp), intent(out)  :: u2(N)

    ! u1 = fcn(nodes)
    u1 = fcn(nodes)

    ! u2 = d_fcn(nodes)
    u2 = dfcn(nodes)
  end subroutine compute_cubic_coeffs


  ! -----------------------------------------------------------------
  ! Initial profiles
  ! -----------------------------------------------------------------
  function rho(x) result(res)
    real(dp), intent(in) :: x(:)
    real(dp) :: res(size(x))
    res = gaussian(x, 0.0d0, 1.0d0)
  end function rho

  function drho(x) result(res)
    real(dp), intent(in) :: x(:)
    real(dp) :: res(size(x))
    res = dgaussian(x, 0.0d0, 1.0d0)
  end function drho


  function T(x) result(res)
    real(dp), intent(in) :: x(:)
    real(dp) :: res(size(x))
    res = gaussian(x, 0.0d0, 1.0d0)
  end function T

  function dT(x) result(res)
    real(dp), intent(in) :: x(:)
    real(dp) :: res(size(x))
    res = dgaussian(x, 0.0d0, 1.0d0)
  end function dT


  function v1(x) result(res)
    real(dp), intent(in) :: x(:)
    real(dp) :: res(size(x))
    res = gaussian(x, 0.0d0, 1.0d0)
  end function v1

  function dv1(x) result(res)
    real(dp), intent(in) :: x(:)
    real(dp) :: res(size(x))
    res = dgaussian(x, 0.0d0, 1.0d0)
  end function dv1


  real(dp) elemental function gaussian(x, mean, std_dev)
    real(dp), intent(in) :: x
    real(dp), intent(in) :: mean
    real(dp), intent(in) :: std_dev

    gaussian = (1 / (std_dev * sqrt(2*dpi))) * exp(-0.5 * ((x - mean)/std_dev)**2)

  end function gaussian

  real(dp) elemental function dgaussian(x, mean, std_dev)
    real(dp), intent(in) :: x
    real(dp), intent(in) :: mean
    real(dp), intent(in) :: std_dev

    dgaussian = -((x - mean) / (std_dev**2)) * gaussian(x, mean, std_dev)

  end function dgaussian

end module mod_initial_condition