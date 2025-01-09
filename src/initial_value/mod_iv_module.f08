module mod_iv_module
  use mod_global_variables, only: dp
  use mod_physical_constants, only: dpi
  use mod_logging, only: logger, str
  use mod_settings, only: settings_t
  use mod_state_vector_component, only: sv_component_t
  use mod_basis_functions, only: hcubic, hquad, basis_function
  use mod_grid, only: grid_t

  use mod_iv_state_vector, only: iv_state_vector_t, init_and_bind
  implicit none

  private


contains

  subroutine initialise_iv_module(settings)
    type(settings_t), intent(in) :: settings
    
    type(iv_state_vector_t) :: state_vec

    state_vec = init_and_bind(settings%state_vector)

  end subroutine
  

end module mod_iv_module