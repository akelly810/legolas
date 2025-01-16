module mod_iv_state_vector
  use mod_global_variables, only: dp
  use mod_logging, only: logger, str
  use mod_settings, only: settings_t
  use mod_state_vector, only: state_vector_t
  use mod_state_vector_names
  use mod_iv_globals, only: profile_fcn, iv_prof_fcn_ptr_t
  use mod_iv_state_vector_component, only: iv_sv_component_t, new_iv_component
  use mod_iv_initial_conditions, only: get_f_lists
  
  implicit none

  private

  ! Type to hold the component pointers
  type, private :: iv_comp_ptr_t
    class(iv_sv_component_t), pointer :: ptr
  end type iv_comp_ptr_t

  type, public :: iv_state_vector_t
    type(state_vector_t), pointer :: base => null()  ! pointer to existing state_vector_t instance
    class(iv_comp_ptr_t), allocatable :: components(:)
    logical :: is_initialised = .false.
    integer :: num_components
    integer :: stride

    real(dp), allocatable :: x0(:)
    complex(dp), allocatable :: x0_cplx(:)  ! FIXME: properly handle real vs complex mixing..

    contains
      procedure, public :: initialise_components
      procedure, public :: assemble_iv_array
      procedure, public :: disassemble_iv_array
      procedure, public :: reassemble_from_block

  end type iv_state_vector_t

  type(iv_sv_component_t), public, protected, target :: iv_rho1
  type(iv_sv_component_t), public, protected, target :: iv_v1
  type(iv_sv_component_t), public, protected, target :: iv_v2
  type(iv_sv_component_t), public, protected, target :: iv_v3
  type(iv_sv_component_t), public, protected, target :: iv_T1
  type(iv_sv_component_t), public, protected, target :: iv_a1
  type(iv_sv_component_t), public, protected, target :: iv_a2
  type(iv_sv_component_t), public, protected, target :: iv_a3

  public :: init_and_bind

  contains
  ! -----------------------------------------------------------------
  ! Constructor
  ! -----------------------------------------------------------------
  function init_and_bind(base_instance) result(iv_state_vector)
    type(state_vector_t), intent(in), target :: base_instance
    type(iv_state_vector_t) :: iv_state_vector

    ! Make base point to the passed-in instance.
    iv_state_vector%base => base_instance
  end function init_and_bind


  ! TODO: Test this
  subroutine initialise_components(self, physics_type)
    class(iv_state_vector_t), intent(inout) :: self
    character(len=*), intent(in) :: physics_type

    type(iv_prof_fcn_ptr_t), allocatable :: f_list(:), df_list(:)
    integer :: i

    if (self%is_initialised) then
      call logger%error("IV state vector is already initialised.")
    end if

    ! Initialise all of the components
    iv_rho1 = new_iv_component()
    iv_v1 = new_iv_component()
    iv_v2 = new_iv_component()
    iv_v3 = new_iv_component()
    iv_T1 = new_iv_component()
    iv_a1 = new_iv_component()
    iv_a2 = new_iv_component()
    iv_a3 = new_iv_component()

    ! Store pointers to the active components
    select case(physics_type)
    case("hd")
      self%num_components = 5
      allocate(self%components(self%num_components), &
               f_list(self%num_components), &
               df_list(self%num_components))

      self%components(1)%ptr => iv_rho1
      self%components(2)%ptr => iv_v1
      self%components(3)%ptr => iv_v2
      self%components(4)%ptr => iv_v3
      self%components(5)%ptr => iv_T1

      call get_f_lists(f_list, df_list)
    case("hd-1d")
      self%num_components = 3
      allocate(self%components(self%num_components), &
               f_list(self%num_components), &
               df_list(self%num_components))

      self%components(1)%ptr => iv_rho1
      self%components(2)%ptr => iv_v1
      self%components(3)%ptr => iv_T1

      call get_f_lists(f_list, df_list)
    case default
      self%num_components = 8
      allocate(self%components(self%num_components), &
               f_list(self%num_components), &
               df_list(self%num_components))

      self%components(1)%ptr => iv_rho1
      self%components(2)%ptr => iv_v1
      self%components(3)%ptr => iv_v2
      self%components(4)%ptr => iv_v3
      self%components(5)%ptr => iv_T1
      self%components(6)%ptr => iv_a1
      self%components(7)%ptr => iv_a2
      self%components(8)%ptr => iv_a3

      call get_f_lists(f_list, df_list)
    end select

    self%stride = 2 * self%num_components

    if (self%num_components /= size(self%base%components)) then
      call logger%error("Mismatch in number of state components used by IV module and main Legolas program.")
    end if

    do i = 1, self%num_components
      ! TODO: Add a 'name' attribute to the iv_components as well and check that iv_comp%name == sv_comp%name
      call self%components(i)%ptr%bind_iv_component(self%base%components(i)%ptr, &
                                                    fcn = f_list(i)%ptr, &
                                                    dfcn = df_list(i)%ptr)
    end do

    self%is_initialised = .true.

  end subroutine initialise_components


  ! -----------------------------------------------------------------
  ! Assemble IV array from profiles
  ! -----------------------------------------------------------------
  subroutine assemble_iv_array(self, N, nodes)
    !>
    class(iv_state_vector_t), intent(inout) :: self
    !> Number of grid points
    integer, intent(in) :: N
    !> Array of grid points
    real(dp), intent(in) :: nodes(N)

    integer :: i, j, idx

    ! FIXME: Think about this... move somewhere else?
    allocate(self%x0(self%stride * N))

    ! Compute coefficients for each component
    do i = 1, self%num_components
      call self%components(i)%ptr%compute_coeffs(N, nodes)
    end do

    ! Assembly of interleaved array into block format
    do i = 1, self%num_components  ! components
      do j = 1, 2                  ! c1 or c2
        idx = 2*(i - 1) + j        ! compute the offset in x0
        select case(j)
          case(1)  ! c1
            self%x0(idx : self%stride*N : self%stride) = self%components(i)%ptr%c1
          case(2)  ! c2
            self%x0(idx : self%stride*N : self%stride) = self%components(i)%ptr%c2
        end select
      end do
    end do

    ! FIXME: For now, just hold 2 separate arrays.. one for real and one for complex
    allocate(self%x0_cplx(self%stride * N))
    self%x0_cplx = cmplx(self%x0, 0.0d0, kind = dp)
  end subroutine assemble_iv_array


  !>
  !!
  subroutine reassemble_from_block(self, N_fine, x_fine, N, nodes)
    class(iv_state_vector_t), intent(inout) :: self
    !> Number of grid points for reassembly
    integer, intent(in) :: N_fine
    !> Array of grid points for reassembly
    real(dp), intent(in) :: x_fine(N_fine)
    !> Number of grid points in block structure
    integer, intent(in) :: N
    !> Array of grid points in block structure
    real(dp), intent(in) :: nodes(N)

    integer :: i

    ! Get each component to compute its profile
    do i = 1, self%num_components
      call self%components(i)%ptr%reconstruct_profile(N_fine, x_fine, N, nodes)
    end do

  end subroutine reassemble_from_block

  !> Gather each component's c1 and c2 coefficients from the
  !! block-format array x0. This is the inverse of assemble_iv_array.
  subroutine disassemble_iv_array(self, N)
    class(iv_state_vector_t), intent(inout) :: self
    !> # of grid points (nodes) for each component
    integer, intent(in) :: N
  
    integer :: i, j, idx
  
    ! Loop over all components and fill c1, c2 from x0
    do i = 1, self%num_components
      ! Make sure c1 and c2 are allocated
      if (.not. allocated(self%components(i)%ptr%c1)) then
        allocate(self%components(i)%ptr%c1(N), self%components(i)%ptr%c2(N))
      end if
  
      do j = 1, 2  ! c1 or c2
        idx = 2*(i - 1) + j  ! offset in the x0 array
        select case (j)
          case(1)  ! c1
            self%components(i)%ptr%c1 = self%x0(idx : self%stride*N : self%stride)
          case(2)  ! c2
            self%components(i)%ptr%c2 = self%x0(idx : self%stride*N : self%stride)
        end select
      end do
    end do
  
  end subroutine disassemble_iv_array
  

end module mod_iv_state_vector