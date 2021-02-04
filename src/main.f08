! =============================================================================
!> Main program for the Legolas finite element code.
!! Matrices, eigenvalues and left/right eigenvectors are defined here and passed
!! on to the different modules and submodules.
!!
!! <tt>Legolas</tt> is currently being developed by Niels Claes, Jordi De Jonghe
!! and Rony Keppens, at the Centre for mathematical Plasma-Astrophysics (CmPA),
!! KU Leuven, Belgium.
program legolas
  use mod_global_variables, only: dp, str_len, show_results, dry_run
  use mod_matrix_creation, only: create_matrices
  use mod_solvers, only: solve_evp
  use mod_output, only: datfile_name
  use mod_logging, only: log_message, print_console_info, print_whitespace
  use mod_inspections, only: handle_spurious_eigenvalues
  implicit none

  !> A matrix in eigenvalue problem wBX = AX
  complex(dp), allocatable  :: matrix_A(:, :)
  !> B matrix in eigenvalue problem wBX = AX
  real(dp), allocatable     :: matrix_B(:, :)
  !> array with eigenvalues
  complex(dp), allocatable  :: omega(:)
  !> matrix with right eigenvectors, column indices correspond to omega indices
  complex(dp), allocatable  :: eigenvecs_right(:, :)
  !> matrix with left eigenvectors, column indices correspond to omega indices
  complex(dp), allocatable  :: eigenvecs_left(:, :)

  call initialisation()
  call create_matrices(matrix_B, matrix_A)
  call print_console_info()

  if (.not. dry_run) then
    call log_message("solving eigenvalue problem...", level='info')
    call solve_evp(matrix_A, matrix_B, omega, eigenvecs_left, eigenvecs_right)
  else
    call log_message("running dry, overriding parfile and setting &
                      &eigenvalues to zero", level='info')
    omega = (0.0d0, 0.0d0)
  end if

  call handle_spurious_eigenvalues(omega)

  call finalise_results()
  call cleanup()

  if (show_results) then
    call print_whitespace(1)
    call execute_command_line("python3 pylbo_wrapper.py -i " // trim(datfile_name))
  end if

contains

  !> Subroutine responsible for all initialisations.
  !! Allocates and initialises main and global variables, then the equilibrium state
  !! and eigenfunctions are initialised and the equilibrium is set.
  subroutine initialisation()
    use mod_global_variables, only: initialise_globals, matrix_gridpts, hall_mhd
    use mod_input, only: read_parfile, get_parfile
    use mod_equilibrium, only: initialise_equilibrium, set_equilibrium, hall_field
    use mod_eigenfunctions, only: initialise_eigenfunctions
    use mod_logging, only: print_logo
    use mod_hallmhd, only: set_hallfactor

    character(len=str_len)  :: parfile

    call initialise_globals()
    call get_parfile(parfile)
    call read_parfile(parfile)

    call print_logo()

    allocate(matrix_A(matrix_gridpts, matrix_gridpts))
    allocate(matrix_B(matrix_gridpts, matrix_gridpts))
    allocate(omega(matrix_gridpts))
    allocate(eigenvecs_right(matrix_gridpts, matrix_gridpts))
    allocate(eigenvecs_left(matrix_gridpts, matrix_gridpts))

    call initialise_equilibrium()
    call initialise_eigenfunctions()
    call set_equilibrium()
    ! Set Hall factor if needed
    if (hall_mhd) then
      call set_hallfactor(hall_field)
    end if
  end subroutine initialisation


  !> Wraps up results and writes output.
  !! Makes a call to the eigenfunctions subroutine if specified in the parfile,
  !! then calls the output routines to write the datfile.
  subroutine finalise_results()
    use mod_global_variables, only: write_eigenfunctions
    use mod_output, only: create_datfile
    use mod_eigenfunctions, only: calculate_eigenfunctions

    if (write_eigenfunctions) then
      call calculate_eigenfunctions(eigenvecs_right)
    end if
    call create_datfile(omega, matrix_A, matrix_B)
  end subroutine finalise_results


  !> Deallocates all main variables, then calls the cleanup
  !! routines of all relevant subroutines to do the same thing.
  subroutine cleanup()
    use mod_global_variables, only: radiative_cooling
    use mod_grid, only: grid_clean
    use mod_equilibrium, only: equilibrium_clean
    use mod_radiative_cooling, only: radiative_cooling_clean
    use mod_eigenfunctions, only: eigenfunctions_clean

    deallocate(matrix_A)
    deallocate(matrix_B)
    deallocate(omega)
    deallocate(eigenvecs_left)
    deallocate(eigenvecs_right)

    call grid_clean()
    call equilibrium_clean()

    if (radiative_cooling) then
      call radiative_cooling_clean()
    end if
    call eigenfunctions_clean()
  end subroutine cleanup

end program legolas
