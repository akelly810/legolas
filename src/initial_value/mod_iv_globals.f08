module mod_iv_globals
  use mod_global_variables, only: dp  
  implicit none
  
  public

  abstract interface
    function profile_fcn(x) result(res)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: res(size(x))
    end function profile_fcn
  end interface

  ! Type to hold the profile function pointers
  type :: iv_prof_fcn_ptr_t
    procedure(profile_fcn), pointer, nopass :: ptr
  end type iv_prof_fcn_ptr_t


end module mod_iv_globals