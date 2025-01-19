module mod_iv_module
  use mod_global_variables, only: dp
  use mod_logging, only: logger, str
  use mod_matrix_structure, only: matrix_t
  use mod_settings, only: settings_t
  use mod_state_vector_component, only: sv_component_t
  use mod_grid, only: grid_t
  use mod_iv_globals, only: linspace
  use mod_iv_state_vector, only: iv_state_vector_t, init_and_bind
  use mod_iv_solver, only: solve
  implicit none

  private

  type, public :: iv_module_t
    logical, private :: is_initialised
    type(settings_t), pointer, private :: settings
    type(grid_t), pointer, private :: grid
    type(iv_state_vector_t), private :: state_vec

    real, allocatable :: iv_grid(:) 

    complex(dp), allocatable :: snapshots(:,:)
  contains

  procedure, public :: initialise
  procedure, public :: solve_ivp

  procedure, public :: postprocess_snapshots

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

    ! Setup snapshots array
    allocate(self%snapshots(self%state_vec%stride * self%settings%grid%get_gridpts(), self%settings%iv%n_snapshots))

    self%is_initialised = .true.
  end subroutine initialise


  subroutine solve_ivp(self, matrix_A, matrix_B)
    class(iv_module_t), intent(inout) :: self
    type(matrix_t), intent(in) :: matrix_A
    type(matrix_t), intent(in) :: matrix_B

    ! Now call the solver
    call solve(matrix_A, matrix_B, self%state_vec%x0_cplx, self%settings, self%snapshots)
    
  end subroutine solve_ivp


  subroutine postprocess_snapshots(self)
    class(iv_module_t), intent(inout) :: self

    integer :: i_snap, i
    integer :: N_fine
    real(dp), allocatable :: iv_grid(:)

    ! Build fine grid for plotting
    N_fine = self%settings%iv%get_iv_gridpts()
    allocate(iv_grid(N_fine))

    iv_grid = linspace(self%settings%grid%get_grid_start(), self%settings%grid%get_grid_end(), N_fine)

    ! Loop over each snapshot
    do i_snap = 1, self%settings%iv%n_snapshots
      ! 1. Update the state vector
      self%state_vec%x0_cplx = self%snapshots(:, i_snap)
      self%state_vec%x0 = real(self%state_vec%x0_cplx)

      ! 2. Re-compute c1, c2 from x0_cplx
      call self%state_vec%disassemble_iv_array(self%settings%grid%get_gridpts())

      ! 3. Reconstruct
      call self%state_vec%reassemble_from_block( &
          N_fine, iv_grid, &
          self%settings%grid%get_gridpts(), &
          self%grid%base_grid )

      ! 4. Output the profile
      ! TODO: Integrate this with data io module
      open(unit=30, file="data/snapshot_"//trim(adjustl(str(i_snap)))//".txt", &
          status="unknown", action="write", form="formatted")

      do i = 1, size(self%state_vec%components(1)%ptr%profile)
        write(30, '(E12.5)') self%state_vec%components(1)%ptr%profile(i)
      end do

      close(30)
    end do

    deallocate(iv_grid)

  end subroutine postprocess_snapshots


end module mod_iv_module