! =============================================================================
!> Initial value solver
!!
!!
module mod_iv_solvers
  use mod_logging, only: logger, str
  use mod_settings, only: settings_t
  use mod_matrix_structure, only: matrix_t
  use mod_global_variables, only: dp, NaN

  use mod_banded_matrix, only: banded_matrix_t, new_banded_matrix
  use mod_banded_operations, only: multiply
  use mod_linear_systems, only: solve_linear_system_complex_banded
  use mod_transform_matrix, only: array_to_banded
  implicit none

  real, parameter :: alpha = 0.5  ! trapezoidal method param
  integer         :: num_steps    ! number of time steps

  integer            :: i       ! loop index
  integer, parameter :: ku = 1  ! number of superdiagonals
  integer, parameter :: kl = 1  ! number of subdiagonals   

  ! Mx = rhs
  complex(8), allocatable :: rhs(:)         ! rhs is a vector
  complex(8), allocatable :: matrix_M(:,:)  ! dense matrix M
  type(banded_matrix_t) :: M                ! banded matrix

contains

  !>
  !
  subroutine solve_ivp(A, B, x, dt, t_end, hist)
    complex(8), dimension(:,:)              :: A, B            ! FEM matrices
    complex(8), dimension(:), intent(inout) :: x               ! initial condition, updated with final result
    real, intent(in)                        :: dt              ! timestep
    real, intent(in)                        :: t_end           ! end time
    complex(8), dimension(:,:), optional, intent(out) :: hist  ! save history for plotting (optional)

    ! Data setup
    num_steps = nint(t_end / dt)
    allocate(rhs, mold = x)
    allocate(matrix_M, mold = A)

    ! Trapezoidal (theta) method
    matrix_M = B - dt * alpha * A                       ! compute M matrix
    call array_to_banded(matrix_M, kl, ku, M)           ! convert dense to banded

    do i = 1, num_steps
      rhs = matmul(B + (1 - alpha) * dt * A, x)       ! compute rhs vector
      x = solve_linear_system_complex_banded(M, rhs)  ! solve banded system with zgbsv

      ! Save in history
      if (present(hist)) then
        hist(i,:) = x
      end if
    end do

    ! Unit tests require deallocation
    deallocate(rhs)
    deallocate(matrix_M)

  end subroutine solve_ivp


  !> Helper routine for logging matrices nicely to a file.
  !  Prints the real parts of the entries of a matrix of type COMPLEX(8)
  subroutine log_matrix(matrix, rows, cols, label)
    complex(8), intent(in) :: matrix(:,:)
    integer, intent(in) :: rows, cols
    character(len=*), intent(in), optional :: label
    integer :: i, j
    integer :: log_unit = 10  ! Arbitrary unit number for logging

    open(log_unit, file="matrix_log.txt", status="unknown")
    if (present(label)) then
        write(log_unit,*) "Matrix: ", trim(label)
    else
        write(log_unit,*) "Matrix:"
    end if

    ! Print the matrix row by row
    do i = 1, rows
      ! write(log_unit, "(A, I3, A)", advance="no") "Row ", i, ": "
      do j = 1, cols
        write(log_unit, "(F10.4)", advance="no") real(matrix(i, j))
      end do
      write(log_unit,*)  ! End the row
    end do

    close(log_unit)
end subroutine log_matrix

end module mod_iv_solvers