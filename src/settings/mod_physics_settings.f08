module mod_physics_settings
  use mod_global_variables, only: dp
  use mod_flow_settings, only: flow_settings_t, new_flow_settings
  use mod_cooling_settings, only: cooling_settings_t, new_cooling_settings
  use mod_gravity_settings, only: gravity_settings_t, new_gravity_settings
  use mod_resistivity_settings, only: resistivity_settings_t, new_resistivity_settings
  use mod_viscosity_settings, only: viscosity_settings_t, new_viscosity_settings
  use mod_conduction_settings, only: conduction_settings_t, new_conduction_settings
  use mod_hall_settings, only: hall_settings_t, new_hall_settings
  implicit none

  private

  type, public :: physics_t
    real(dp), private :: gamma
    logical :: is_incompressible
    type(flow_settings_t) :: flow
    type(cooling_settings_t) :: cooling
    type(gravity_settings_t) :: gravity
    type(resistivity_settings_t) :: resistivity
    type(viscosity_settings_t) :: viscosity
    type(conduction_settings_t) :: conduction
    type(hall_settings_t) :: hall

  contains

    procedure, public :: set_gamma
    procedure, public :: get_gamma
    procedure, public :: get_gamma_1
    procedure, public :: set_incompressible
    procedure, public :: set_defaults

    procedure, public :: enable_flow
    procedure, public :: enable_cooling
    procedure, public :: enable_gravity
    procedure, public :: enable_resistivity
    procedure, public :: enable_viscosity
    procedure, public :: enable_parallel_conduction
    procedure, public :: enable_perpendicular_conduction
    procedure, public :: enable_hall
  end type physics_t

  public :: new_physics_settings

contains

  pure function new_physics_settings() result(physics)
    type(physics_t) :: physics

    physics%is_incompressible = .false.
  end function new_physics_settings


  pure subroutine set_gamma(this, gamma)
    class(physics_t), intent(inout) :: this
    real(dp), intent(in) :: gamma
    this%gamma = gamma
  end subroutine set_gamma


  pure real(dp) function get_gamma(this)
    class(physics_t), intent(in) :: this
    get_gamma = this%gamma
  end function get_gamma


  pure real(dp) function get_gamma_1(this)
    class(physics_t), intent(in) :: this
    get_gamma_1 = this%get_gamma() - 1.0_dp
  end function get_gamma_1


  pure subroutine set_incompressible(this)
    class(physics_t), intent(inout) :: this
    this%is_incompressible = .true.
    call this%set_gamma(1.0e12_dp)
  end subroutine set_incompressible


  pure subroutine set_defaults(this)
    class(physics_t), intent(inout) :: this
    call this%set_gamma(5.0_dp / 3.0_dp)
    this%is_incompressible = .false.
  end subroutine set_defaults


  pure subroutine enable_flow(this)
    class(physics_t), intent(inout) :: this
    call this%flow%enable()
  end subroutine enable_flow


  pure subroutine enable_cooling(this, cooling_curve, interpolation_points)
    class(physics_t), intent(inout) :: this
    character(len=*), intent(in) :: cooling_curve
    integer, intent(in), optional :: interpolation_points

    if (present(interpolation_points)) then
      call this%cooling%set_interpolation_points(interpolation_points)
    end if
    call this%cooling%set_cooling_curve(cooling_curve)
    call this%cooling%enable()
  end subroutine enable_cooling


  pure subroutine enable_gravity(this)
    class(physics_t), intent(inout) :: this
    call this%gravity%enable()
  end subroutine enable_gravity


  pure subroutine enable_resistivity(this, fixed_eta_value)
    class(physics_t), intent(inout) :: this
    real(dp), intent(in), optional :: fixed_eta_value

    call this%resistivity%enable()
    if (present(fixed_eta_value)) then
      call this%resistivity%set_fixed_resistivity(fixed_eta_value)
    end if
  end subroutine enable_resistivity


  pure subroutine enable_viscosity(this, fixed_viscosity_value, viscous_heating)
    class(physics_t), intent(inout) :: this
    real(dp), intent(in), optional :: fixed_viscosity_value
    logical, intent(in), optional :: viscous_heating

    call this%viscosity%enable()
    if (present(fixed_viscosity_value)) then
      call this%viscosity%set_fixed_viscosity(fixed_viscosity_value)
    end if
    if (present(viscous_heating)) then
      call this%viscosity%enable_viscous_heating()
    end if
  end subroutine enable_viscosity


  pure subroutine enable_parallel_conduction(this, fixed_tc_para_value)
    class(physics_t), intent(inout) :: this
    real(dp), intent(in), optional :: fixed_tc_para_value

    call this%conduction%enable_para_conduction()
    if (present(fixed_tc_para_value)) then
      call this%conduction%set_fixed_tc_para(fixed_tc_para_value)
    end if
  end subroutine enable_parallel_conduction


  pure subroutine enable_perpendicular_conduction(this, fixed_tc_perp_value)
    class(physics_t), intent(inout) :: this
    real(dp), intent(in), optional :: fixed_tc_perp_value

    call this%conduction%enable_perp_conduction()
    if (present(fixed_tc_perp_value)) then
      call this%conduction%set_fixed_tc_perp(fixed_tc_perp_value)
    end if
  end subroutine enable_perpendicular_conduction


  pure subroutine enable_hall( &
    this, use_hall_substitution, electron_inertia, electron_fraction &
  )
    class(physics_t), intent(inout) :: this
    logical, intent(in), optional :: use_hall_substitution
    logical, intent(in), optional :: electron_inertia
    real(dp), intent(in), optional :: electron_fraction
    logical :: hall_substitution

    hall_substitution = .true.
    if (present(use_hall_substitution)) hall_substitution = use_hall_substitution
    call this%hall%enable(hall_substitution)

    if (present(electron_inertia)) call this%hall%enable_electron_inertia()
    if (present(electron_fraction)) then
      call this%hall%set_electron_fraction(electron_fraction)
    end if
  end subroutine enable_hall

end module mod_physics_settings
