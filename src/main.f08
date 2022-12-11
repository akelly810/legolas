! =============================================================================
!> Main program for the Legolas finite element code.
!! Matrices, eigenvalues and left/right eigenvectors are defined here and passed
!! on to the different modules and submodules.
!!
!! <tt>Legolas</tt> is currently being developed by Niels Claes, Jordi De Jonghe
!! and Rony Keppens, at the Centre for mathematical Plasma-Astrophysics (CmPA),
!! KU Leuven, Belgium.
program legolas
  use mod_global_variables, only: dp, str_len
  use mod_matrix_structure, only: matrix_t
  use mod_matrix_manager, only: build_matrices
  use mod_solvers, only: solve_evp
  use mod_output, only: datfile_name, create_datfile
  use mod_logging, only: log_message, str
  use mod_console, only: print_console_info, print_whitespace
  use mod_inspections, only: handle_spurious_eigenvalues
  use mod_timing, only: timer_t, new_timer
  use mod_settings, only: settings_t, new_settings
  implicit none

  !> A matrix in eigenvalue problem wBX = AX
  type(matrix_t) :: matrix_A
  !> B matrix in eigenvalue problem wBX = AX
  type(matrix_t) :: matrix_B
  !> timer used by the whole program
  type(timer_t) :: timer
  !> dedicated settings type
  type(settings_t) :: settings
  !> array with eigenvalues
  complex(dp), allocatable  :: omega(:)
  !> matrix with right eigenvectors, column indices correspond to omega indices
  complex(dp), allocatable  :: eigenvecs_right(:, :)

  timer = new_timer()
  settings = new_settings()

  call timer%start_timer()
  call initialisation()
  timer%init_time = timer%end_timer()

  call print_console_info(settings)

  call timer%start_timer()
  call build_matrices(matrix_B, matrix_A, settings)
  timer%matrix_time = timer%end_timer()

  call log_message("solving eigenvalue problem...", level="info")
  call timer%start_timer()
  call solve_evp(matrix_A, matrix_B, settings, omega, eigenvecs_right)
  timer%evp_time = timer%end_timer()

  call timer%start_timer()
  call create_eigenfunctions()
  timer%eigenfunction_time = timer%end_timer()

  call timer%start_timer()
  call create_datfile(omega, matrix_A, matrix_B, eigenvecs_right, settings)
  timer%datfile_time = timer%end_timer()

  call cleanup()

  call print_timelog()

  if (settings%io%show_results) then
    call print_whitespace(1)
    call execute_command_line("python3 pylbo_wrapper.py -i " // trim(datfile_name))
  end if

contains

  !> Subroutine responsible for all initialisations.
  !! Allocates and initialises main and global variables, then the equilibrium state
  !! and eigenfunctions are initialised and the equilibrium is set.
  subroutine initialisation()
    use mod_global_variables, only: initialise_globals, &
     gamma, set_gamma, NaN, hall_mhd, x_start, x_end
    use mod_matrix_structure, only: new_matrix
    use mod_input, only: read_parfile, get_parfile
    use mod_equilibrium, only: initialise_equilibrium, set_equilibrium, hall_field
    use mod_console, only: print_logo

    character(len=5*str_len)  :: parfile
    integer   :: nb_evs

    real(dp) :: ratio
    ratio = NaN

    call initialise_globals()
    call get_parfile(parfile)
    call read_parfile(parfile, settings)
    call set_gamma(gamma)

    call print_logo()
    call log_message( &
      "the physics type is " // settings%get_physics_type(), level="info" &
    )
    call log_message( &
      "the state vector is " // str(settings%get_state_vector()), level="info" &
    )

    select case(settings%solvers%get_solver())
    case ("arnoldi")
      nb_evs = settings%solvers%number_of_eigenvalues
    case ("inverse-iteration")
      nb_evs = 1
    case default
      nb_evs = settings%dims%get_dim_matrix()
    end select
    call log_message("setting #eigenvalues to " // str(nb_evs), level="debug")
    allocate(omega(nb_evs))
    matrix_A = new_matrix(nb_rows=settings%dims%get_dim_matrix(), label="A")
    matrix_B = new_matrix(nb_rows=settings%dims%get_dim_matrix(), label="B")

    call initialise_equilibrium()
    call set_equilibrium()

    if (hall_mhd) then
      ratio = maxval(hall_field % hallfactor) / (x_end - x_start)
      if (ratio > 0.1d0) then
        call log_message("large ratio Hall scale / system scale: " // str(ratio), level="warning")
      end if
    end if

    ! Arnoldi solver needs this, since it always calculates an orthonormal basis
    if ( &
      settings%io%should_compute_eigenvectors() &
      .or. settings%solvers%get_solver() == "arnoldi" &
    ) then
      call log_message("allocating eigenvector arrays", level="debug")
      ! we need #rows = matrix dimension, #cols = #eigenvalues
      allocate(eigenvecs_right(settings%dims%get_dim_matrix(), nb_evs))
    else
      ! @note: this is needed to prevent segfaults, since it seems that in some
      ! cases for macOS the routine zgeev references the right eigenvectors even
      ! if they are not requested.
      call log_message("allocating eigenvector arrays as dummy", level="debug")
      allocate(eigenvecs_right(2, 2))
    end if
  end subroutine initialisation


  !> Initialises and calculates the eigenfunctions if requested.
  subroutine create_eigenfunctions()
    use mod_eigenfunctions, only: initialise_eigenfunctions, calculate_eigenfunctions

    if (settings%io%write_eigenfunctions) then
      call initialise_eigenfunctions(omega, settings)
      call calculate_eigenfunctions(eigenvecs_right, settings)
    end if
  end subroutine create_eigenfunctions


  !> Deallocates all main variables, then calls the cleanup
  !! routines of all relevant subroutines to do the same thing.
  subroutine cleanup()
    use mod_global_variables, only: radiative_cooling
    use mod_grid, only: grid_clean
    use mod_equilibrium, only: equilibrium_clean
    use mod_radiative_cooling, only: radiative_cooling_clean
    use mod_eigenfunctions, only: eigenfunctions_clean

    call matrix_A%delete_matrix()
    call matrix_B%delete_matrix()
    deallocate(omega)
    if (allocated(eigenvecs_right)) then
      deallocate(eigenvecs_right)
    end if

    call grid_clean()
    call equilibrium_clean()

    if (radiative_cooling) then
      call radiative_cooling_clean()
    end if
    call eigenfunctions_clean()
    call settings%delete()
  end subroutine cleanup


  subroutine print_timelog()
    use mod_logging, only: override_prefix_to_false
    real(dp) :: total_time

    call print_whitespace(1)
    call log_message("---------------------------------------------")
    override_prefix_to_false = .true.
    total_time = timer%get_total_time()

    call log_message("                << Time log >>")
    call log_message( &
      "Legolas finished in " // str(total_time) // " seconds", level="info" &
    )
    call log_message( &
      "   initialisation: " // str(timer%init_time) // " sec", level="info" &
    )
    call log_message( &
      "   matrix construction: " // str(timer%matrix_time) // " sec", &
      level="info" &
    )
    call log_message( &
      "   eigenvalue problem: " // str(timer%evp_time) // " sec", &
      level="info" &
    )
    call log_message( &
      "   eigenfunction assembly: " // str(timer%eigenfunction_time) // " sec", &
      level="info" &
    )
    call log_message( &
      "   datfile creation: " // str(timer%datfile_time) // " sec", level="info" &
    )
  end subroutine print_timelog

end program legolas
