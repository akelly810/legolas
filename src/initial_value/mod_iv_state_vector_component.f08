module mod_iv_state_vector_component
  use mod_global_variables, only: dp, str_len_arr
  use mod_logging, only: logger
  use mod_iv_globals, only: profile_fcn
  use mod_state_vector_component, only: sv_component_t
  implicit none

  type, public :: iv_sv_component_t
    type(sv_component_t), pointer :: base => null()  ! pointer to existing sv_component_t instance
    logical :: is_bound = .false.

    real(dp), allocatable :: c1(:), c2(:)  ! coefficients

    ! Callback procedure pointers for profile function & derivative
    procedure(profile_fcn), nopass, pointer :: p_fcn  => null()
    procedure(profile_fcn), nopass, pointer :: p_dfcn => null()

    contains
      procedure, public :: bind_iv_component
      procedure, public :: compute_coeffs

      procedure, public :: compute_quad_coeffs
      procedure, public :: compute_cubic_coeffs

  end type iv_sv_component_t

  public :: new_iv_component

contains
  ! -----------------------------------------------------------------
  ! Constructor
  ! -----------------------------------------------------------------
  function new_iv_component() result(comp)
    type(iv_sv_component_t) :: comp
  end function new_iv_component


  subroutine bind_iv_component(self, base_instance, fcn, dfcn)
    class(iv_sv_component_t), intent(inout)  :: self
    type(sv_component_t), intent(in), target :: base_instance
    procedure(profile_fcn), pointer          :: fcn
    procedure(profile_fcn), pointer          :: dfcn

    if (.not. self%is_bound) then
      ! Make base point to the passed-in instance
      self%base => base_instance

      ! Store the function pointers
      self%p_fcn  => fcn
      self%p_dfcn => dfcn

      self%is_bound = .true.
    else
      call logger%error("iv_sv_component_t is already bound!")
    end if
  end subroutine bind_iv_component


  ! -----------------------------------------------------------------
  ! Compute coefficients
  ! -----------------------------------------------------------------
  subroutine compute_coeffs(self, N, nodes)
    class(iv_sv_component_t), intent(inout) :: self
    integer, intent(in)                     :: N
    real(dp), intent(in)                    :: nodes(N)

    ! TODO: remove.. this will probably be needed in the reconstruct code
    ! ! Use the parent's public methods to get the basis function pointer
    ! procedure(basis_function), pointer :: f
    ! call self%base%get_spline_function(1, f)

    ! Compute the coefficients associated to each basis fcn.
    select case(self%base%get_basis_function_name())
    case('QUADRATIC')
      call self%compute_quad_coeffs(N, nodes)
    case('CUBIC')
      call self%compute_cubic_coeffs(N, nodes)
    end select
  end subroutine compute_coeffs


  ! -----------------------------------------------------------------
  ! Compute expansion coefficients for each type
  ! -----------------------------------------------------------------
  subroutine compute_quad_coeffs(self, N, nodes)
    class(iv_sv_component_t), intent(inout) :: self
    integer, intent(in)    :: N
    real(dp), intent(in)   :: nodes(N)

    integer :: i
    real(dp), allocatable :: midpoints(:)
    allocate(midpoints(N - 1))

    ! c1 = fcn(midpoints)
    midpoints = 0.5 * (nodes(2:) + nodes(:N-1))
    self%c1 = [(0.0d0, i = 1,N)]
    self%c1(2:) = self%p_fcn(midpoints)  ! 0 followed by N-1 midpoint values

    ! c2 = fcn(nodes)
    self%c2 = self%p_fcn(nodes)

    deallocate(midpoints)
  end subroutine compute_quad_coeffs


  subroutine compute_cubic_coeffs(self, N, nodes)
    class(iv_sv_component_t), intent(inout) :: self
    integer, intent(in)    :: N
    real(dp), intent(in)   :: nodes(N)

    ! u1 = fcn(nodes)
    self%c1 = self%p_fcn(nodes)

    ! u2 = dfcn(nodes)
    self%c2 = self%p_dfcn(nodes)
  end subroutine compute_cubic_coeffs


end module mod_iv_state_vector_component