! =============================================================================
!> Initial Condition
!!

module mod_initial_condition
  use mod_global_variables, only: dp
  use mod_physical_constants, only: dpi
  use mod_logging, only: logger, str
  use mod_settings, only: settings_t
  use mod_state_vector_component, only: sv_component_t
  use mod_basis_functions, only: hcubic, hquad, basis_function
  use mod_grid, only: grid_t
  implicit none

  private

  public :: generate_ic
  public :: reconstruct_profiles

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

  !>
  !!
  !!
  subroutine reconstruct_profiles(N, nodes, N_fine, x_fine, x0_array, rho_out, v1_out, T_out)
    ! type(grid_t), intent(in)         :: grid
    integer, intent(in)              :: N
    real(dp), intent(in)             :: nodes(N)
    integer, intent(in)              :: N_fine
    real(dp), intent(in)             :: x_fine(N_fine)
    real(dp), intent(in), target     :: x0_array(6*N)
    ! reconstructed high-res variables for output / plotting
    real(dp), intent(out)            :: rho_out(N_fine)
    real(dp), intent(out)            :: v1_out(N_fine)
    real(dp), intent(out)            :: T_out(N_fine)

    integer :: i
    real(dp), pointer :: rho_1(:), rho_2(:)
    real(dp), pointer :: v1_1(:), v1_2(:)
    real(dp), pointer :: T_1(:), T_2(:)

    ! Extract components from interleaved array
    rho_1 => x0_array(1:6*N:6)
    rho_2 => x0_array(2:6*N:6)
    v1_1  => x0_array(3:6*N:6)
    v1_2  => x0_array(4:6*N:6)
    T_1   => x0_array(5:6*N:6)
    T_2   => x0_array(6:6*N:6)

    ! Zero the array
    rho_out = 0.0d0
    v1_out = 0.0d0
    T_out = 0.0d0

    ! Replace fine grid with ef_grid?

    do i = 1, N_fine
      rho_out(i) = approximate_u(x_fine(i), N, nodes, rho_1, rho_2, hquad)
      v1_out(i) = approximate_u(x_fine(i), N, nodes, v1_1, v1_2, hcubic)
      T_out(i) = approximate_u(x_fine(i), N, nodes, T_1, T_2, hquad)
    end do

  end subroutine reconstruct_profiles


  ! -----------------------------------------------------------------
  ! Reconstruction helper procedures
  ! -----------------------------------------------------------------

  !>
  !!
  pure real(dp) function approximate_u(x, N, nodes, u1, u2, basis_fcn) result(res)
    !>  
    real(dp), intent(in) :: x
    !> 
    integer, intent(in) :: N
    !> 
    real(dp), intent(in) :: nodes(N)
    !>
    real(dp), pointer, intent(in) :: u1(:)
    !>
    real(dp), pointer, intent(in) :: u2(:)
    !>
    procedure(basis_function) :: basis_fcn

    integer  :: i
    real(dp) :: xL, xR
    real(dp) :: h(4)

    i = find_element_index(N, x, nodes)
    ! If x is not in the domain spanned by the nodes, set to 0 and return
    if (i == 0) then
      res = 0.0d0
      return
    end if

    ! Coordinates of left and right nodes
    xL = nodes(i)
    xR = nodes(i+1)

    h = basis_fcn(x, xL, xR)

    res = u1(i + 1) * h(1) + &
          u1(i)     * h(2) + &
          u2(i + 1) * h(3) + &
          u2(i)     * h(4)  

  end function approximate_u



  !> Find which element interval [x_i, x_{i+1}] contains x.
  !! Return the index i, or if x is out of range return 0.
  !! This could probably be made quicker.
  pure integer function find_element_index(N, x, nodes) result(res)
    !>
    integer, intent(in)  :: N
    !>
    real(dp), intent(in) :: x
    !> 
    real(dp), intent(in) :: nodes(N)

    integer :: i
    res = 0

    do i = 1, N - 1
      if (nodes(i) <= x  .and. x <= nodes(i+1)) then
        res = i
        return
      else
        res = 0
      end if
    end do

  end function find_element_index


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

    gaussian = (1.0d0 / (std_dev * sqrt(2.0d0*dpi))) * exp(-0.5d0 * ((x - mean)/std_dev)**2)

  end function gaussian

  real(dp) elemental function dgaussian(x, mean, std_dev)
    real(dp), intent(in) :: x
    real(dp), intent(in) :: mean
    real(dp), intent(in) :: std_dev

    dgaussian = -((x - mean) / (std_dev**2)) * gaussian(x, mean, std_dev)

  end function dgaussian

end module mod_initial_condition