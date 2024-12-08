! =============================================================================
!> Initial value solver
!!
!!
module mod_iv_solvers
  use mod_logging, only: logger, str
  use mod_settings, only: settings_t
  use mod_matrix_structure, only: matrix_t
  use mod_global_variables, only: dp

  use mod_banded_matrix, only: banded_matrix_t, new_banded_matrix
  use mod_linear_systems, only: solve_linear_system_complex_banded
  use mod_transform_matrix, only: matrix_to_banded
  implicit none

  real, parameter :: alpha = 0.5  ! 0.0/0.5/1.0 for FW Euler / Trapezoidal method / BW Euler

contains

  !> Solve the initial value problem
  !! (AK) NOTE: Should move some logic out of here and take in a 'settings' object.
  subroutine solve_ivp(matrix_A, matrix_B, x, dt, t_end, hist)
    !> FEM matrix
    type(matrix_t)                          :: matrix_A
    !> FEM matrix
    type(matrix_t)                          :: matrix_B
    !> initial condition, gets updated with final result
    complex(dp), dimension(:), intent(inout) :: x
    !> time step size
    real, intent(in)                        :: dt
    !> simulation end time
    real, intent(in)                        :: t_end
    !> save history for plotting (optional)
    complex(dp), dimension(:,:), optional, intent(out) :: hist

    type(banded_matrix_t) :: A, B, M
    integer :: A_ku, A_kl  ! # upper diagonals, # lower diagonals
    integer :: B_ku, B_kl
    complex(dp) :: beta, gamma
    complex(dp), allocatable :: z(:)
    complex(dp), allocatable :: rhs(:)
    integer :: n            ! dimension of A, B, M, x
    integer :: num_steps    ! number of time steps
    integer :: i            ! loop index

    ! check input sanity
    if (.not. (matrix_A%matrix_dim == matrix_B%matrix_dim)) then
      call logger%error("A or B not square, or not compatible")
      return
    end if

    call matrix_A%get_nb_diagonals(ku=A_ku, kl=A_kl)
    call matrix_B%get_nb_diagonals(ku=B_ku, kl=B_kl)

    ! We will only work with banded matrices
    call matrix_to_banded(matrix_A, A_kl, A_ku, A)
    call matrix_to_banded(matrix_B, B_kl, B_ku, B)

    num_steps = nint(t_end / dt)
    allocate(rhs, mold = x)
    allocate(z, mold = x)
    M = new_banded_matrix(B%m, B%n, B%kl, B%ku)
    n = size(x)

    ! Trapezoidal (theta) method

    ! Compute M = B - dt * alpha * A
    ! 1. Copy B into M
    call zcopy(M%get_total_nb_elements(), B%AB, 1, M%AB, 1)

    ! 2. Scale A by (-dt*alpha) and add to M
    gamma = -dt * alpha
    M%AB = M%AB + gamma*A%AB

    do i = 1, num_steps      
      ! rhs = (B + beta * A)x = Bx + beta * Ax
      ! compute as 3 banded level 2 BLAS operations
      beta = (1 - alpha) * dt  ! scalar

      ! 1. compute rhs=Bx
      call zgbmv( &
        'N', &
        B%m, &
        B%n, &
        B%kl, &
        B%ku, &
        (1.0_dp, 0.0_dp), &
        B%AB, &
        size(B%AB, dim = 1), &
        x, &
        1, &
        (0.0_dp, 0.0_dp), &
        rhs, &
        1 &
      )

      ! 2. compute z=Ax
      call zgbmv( &
        'N', &
        A%m, &
        A%n, &
        A%kl, &
        A%ku, &
        (1.0_dp, 0.0_dp), &
        A%AB, &
        size(A%AB, dim = 1), &
        x, &
        1, &
        (0.0_dp, 0.0_dp), &
        z, &
        1 &
      )

      ! 3. compute rhs = rhs + beta*z
      call zaxpy(n, beta, z, 1, rhs, 1)

      ! Solve resulting banded system with zgbsv
      x = solve_linear_system_complex_banded(M, rhs)

      ! Save in history
      if (present(hist)) then
        hist(i,:) = x
      end if
    end do

    ! Unit tests require deallocation
    deallocate(rhs)
    deallocate(z)

  end subroutine solve_ivp


  !> Helper routine for logging matrices nicely to a file.
  !  Prints the real parts of the entries of a matrix of type COMPLEX(dp)
  subroutine log_matrix(matrix, rows, cols, label)
    complex(dp), intent(in) :: matrix(:,:)
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