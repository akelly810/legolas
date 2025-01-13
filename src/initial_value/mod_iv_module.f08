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

  procedure, public :: reconstruct_profiles
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
    allocate(self%snapshots(self%state_vec%stride, self%settings%iv%n_snapshots))

    self%is_initialised = .true.
  end subroutine initialise


  subroutine solve_ivp(self, matrix_A, matrix_B)
    class(iv_module_t), intent(inout) :: self
    type(matrix_t), intent(in) :: matrix_A
    type(matrix_t), intent(in) :: matrix_B

    ! Now call the solver
    call solve(matrix_A, matrix_B, self%state_vec%x0_cplx, self%settings%iv%step_size, self%settings%iv%t_end, &
                self%snapshots, self%settings%iv%snapshot_stride)
    
    ! TODO: Add mode for saving growth rates only.. plot 1st element of block-format vector x0

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

      ! 2. Re-compute c1, c2 from x0_cplx
      !    TODO: Make new procedure: state_vec%disassemble_iv_array(..)

      ! 3. Reconstruct
      call self%state_vec%reassemble_from_block( &
          N_fine, iv_grid, &
          self%settings%grid%get_gridpts(), &
          self%grid%base_grid )

      ! 4. Output the profile
      ! TODO: Integrate this with data io module
      open(unit=30, file="snapshot_"//trim(adjustl(str(i_snap)))//".txt", &
          status="unknown", action="write", form="formatted")

      do i = 1, size(self%state_vec%components(1)%ptr%profile)
        write(30, '(F8.5)') self%state_vec%components(1)%ptr%profile(i)
      end do

      close(30)
    end do

    deallocate(iv_grid)

  end subroutine postprocess_snapshots


  ! TODO: To be replaced by postprocess_snapshots()
  subroutine reconstruct_profiles(self)
    class(iv_module_t), intent(inout) :: self

    ! Contruct a fine linearly-spaced grid iv_grid for output
    integer :: N, i, j
    real(dp), allocatable :: iv_grid(:)
    real(dp) :: grid_start, grid_end, h

    N = self%settings%iv%get_iv_gridpts()
    grid_start = self%settings%grid%get_grid_start()
    grid_end = self%settings%grid%get_grid_end()
    
    allocate(iv_grid(N))
    h = (grid_end - grid_start) / (N - 1)
    do i = 1, N
      iv_grid(i) = grid_start + (i - 1) * h
    end do

    ! Handle profile reconstruction
    call self%state_vec%reassemble_from_block(N, iv_grid, self%settings%grid%get_gridpts(), self%grid%base_grid)

    ! Handle output
    ! TODO: Do this properly
    open( &
    unit=30, &
    file="output.txt", &
    access="stream", &
    status="unknown", &
    action="write", &
    form="formatted" &
    )

    do j = 1, size(self%state_vec%components(1)%ptr%profile)
      write(30, '(F6.4)') self%state_vec%components(1)%ptr%profile(j)
    end do

    close(30)

  end subroutine reconstruct_profiles

end module mod_iv_module