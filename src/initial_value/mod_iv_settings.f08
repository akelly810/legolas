module mod_iv_settings
  use mod_global_variables, only: dp
  implicit none

  private

  type, public :: iv_settings_t
    integer, private :: rec_gridpts  ! gridpoints for reconstruction
    integer :: n_snapshots
    integer :: snapshot_stride      ! save every n-th snapshot

    ! Solver params
    real :: t_start
    real :: t_end
    integer :: n_steps
    real :: alpha
  contains
    procedure, public :: get_rec_gridpts
    procedure, public :: get_step_size
  end type iv_settings_t

  public :: new_iv_settings

contains

  function new_iv_settings() result(iv_settings)
    type(iv_settings_t) :: iv_settings

    ! Set defaults
    iv_settings%rec_gridpts = 100
    iv_settings%snapshot_stride = 10

    iv_settings%alpha = 0.52  ! 0.0/0.5/1.0 for FW Euler / Trapezoidal method / BW Euler

    iv_settings%t_start = 0.0
    iv_settings%t_end = 0.1
    iv_settings%n_steps = 500

    iv_settings%n_snapshots = floor(real(iv_settings%n_steps - 1) / real(iv_settings%snapshot_stride)) + 2

  end function new_iv_settings


  pure integer function get_rec_gridpts(self)
    class(iv_settings_t), intent(in) :: self
    get_rec_gridpts = self%rec_gridpts
  end function get_rec_gridpts


  pure real function get_step_size(self)
    class(iv_settings_t), intent(in) :: self
    get_step_size = (self%t_end - self%t_start) / real(self%n_steps)
  end function get_step_size

end module mod_iv_settings