module mod_iv_settings
  use mod_global_variables, only: dp
  implicit none

  private

  type, public :: iv_settings_t
    integer, private :: iv_gridpts  ! gridpoints for reconstruction
  contains
    procedure, public :: get_iv_gridpts
  end type iv_settings_t

  public :: new_iv_settings

contains

  function new_iv_settings() result(iv_settings)
    type(iv_settings_t) :: iv_settings

    ! Set defaults
    iv_settings%iv_gridpts = 100

  end function new_iv_settings


  pure integer function get_iv_gridpts(self)
    class(iv_settings_t), intent(in) :: self
    get_iv_gridpts = self%iv_gridpts
  end function get_iv_gridpts

end module mod_iv_settings