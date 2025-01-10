module mod_iv_module
  use mod_global_variables, only: dp
  use mod_logging, only: logger, str
  use mod_matrix_structure, only: matrix_t
  use mod_settings, only: settings_t
  use mod_state_vector_component, only: sv_component_t
  use mod_grid, only: grid_t

  use mod_iv_state_vector, only: iv_state_vector_t, init_and_bind
  use mod_iv_solver, only: solve
  implicit none

  private

  type, public :: iv_module_t
    logical, private :: is_initialised
    type(settings_t), pointer, private :: settings
    type(grid_t), pointer, private :: grid
    type(iv_state_vector_t), private :: state_vec
  contains

  procedure, public :: initialise
  procedure, public :: solve_ivp

  end type iv_module_t

  public :: new_iv_module

contains

  function new_iv_module(settings, grid) result(iv_module)
    type(settings_t), target, intent(in) :: settings
    type(grid_t), target, intent(in) :: grid

    type(iv_module_t) :: iv_module

    iv_module%settings => settings
    iv_module%grid => grid
    iv_module%is_initialised = .false.
  end function new_iv_module


  subroutine initialise(self)
    class(iv_module_t), intent(inout) :: self

    if (self%is_initialised) return

    ! Initialise and setup the state vector
    self%state_vec = init_and_bind(self%settings%state_vector)
    call self%state_vec%initialise_components(self%settings%get_physics_type())

    ! Assemble the initial value array in block format
    call self%state_vec%assemble_iv_array(self%settings%grid%get_gridpts(), self%grid%base_grid)

    self%is_initialised = .true.
  end subroutine initialise


  subroutine solve_ivp(self, matrix_A, matrix_B)
    class(iv_module_t), intent(inout) :: self
    type(matrix_t) :: matrix_A
    type(matrix_t) :: matrix_B

    call solve(matrix_A, matrix_B, self%state_vec%x0_cplx, 0.01, 2.0)

  end subroutine solve_ivp

end module mod_iv_module