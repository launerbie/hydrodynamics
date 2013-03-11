import numpy
from amuse.lab import *
N=100; Mtot = 1|units.MSun; Rvir=1|units.RSun; t_end=0.5|units.day; n_steps=10; write_hdf5=None
converter=nbody_system.nbody_to_si(Mtot, Rvir)
bodies = new_plummer_gas_model(N, convert_nbody=converter)
timerange = numpy.linspace(0, t_end.value_in(units.day), n_steps) | units.day

print bodies

def print_energies(hydro):
    Kh = hydro.kinetic_energy
    Th = hydro.thermal_energy
    Vh = hydro.potential_energy
    Kg = hydro.gas_particles.kinetic_energy()
    Tg = hydro.gas_particles.thermal_energy()
    Vg = hydro.gas_particles.potential_energy()
    K = hydro.particles.kinetic_energy()
    T = hydro.particles.thermal_energy()
    V = hydro.particles.potential_energy()
    total_hydro_energy = Kh+Th+Vh
    total_gas_particles_energy = Kg+Tg+Vg
    total_particles_energy = K+T+V
    print total_hydro_energy
    print total_gas_particles_energy
    print total_particles_energy

hydro = Fi(converter)
hydro.gas_particles.add_particles(bodies)

filename = 'plum.hdf5'
hydro_to_framework = hydro.gas_particles.new_channel_to(bodies)
write_set_to_file(bodies.savepoint(0.0 | t_end.unit),\
                         filename, "hdf5")
#print_energies(hydro)

for t in timerange:
    hydro.evolve_model(t)
    hydro_to_framework.copy()
    write_set_to_file(bodies.savepoint(t), filename, "hdf5")

#print_energies(hydro)

hydro.stop()
print bodies
bodies.get_timeline_of_attribute(bodies[-1].key,'z')





for line in dir(hydro):
    try:
        print line,"\t\t", type(eval('hydro.'+line))
    except AttributeError:
        print line,"\t\t", 'AttributeError'
    except Exception:
        print line,"\t\t", 'Exception'


for line in dir(hydro.gas_particles):
    try:
        print line,"\t\t", type(eval('hydro.gas_particles.'+line)), eval('hydro.gas_particles.'+line+'()')
    except AttributeError:
        print line,"\t\t", 'AttributeError'
    except Exception:
        print line,"\t\t", 'Exception'

for line in dir(hydro.particles):
    try:
        print line,"\t\t", type(eval('hydro.particles.'+line)), eval('hydro.particles.'+line+'()')
    except AttributeError:
        print line,"\t\t", 'AttributeError'
    except Exception:
        print line,"\t\t", 'Exception'

MODE_NORMAL 		<type 'str'>
MODE_PERIODIC_BOUNDARIES 		<type 'str'>
NBODY 		<type 'object'>
__class__ 		<type 'type'>
__del__ 		<type 'instancemethod'>
__delattr__ 		<type 'method-wrapper'>
__dict__ 		<type 'dict'>
__dir__ 		<type 'instancemethod'>
__doc__ 		<type 'str'>
__format__ 		<type 'builtin_function_or_method'>
__getattr__ 		<type 'instancemethod'>
__getattribute__ 		<type 'method-wrapper'>
__hash__ 		<type 'method-wrapper'>
__init__ 		<type 'instancemethod'>
__module__ 		<type 'str'>
__new__ 		<type 'builtin_function_or_method'>
__reduce__ 		<type 'builtin_function_or_method'>
__reduce_ex__ 		<type 'builtin_function_or_method'>
__repr__ 		<type 'method-wrapper'>
__setattr__ 		<type 'method-wrapper'>
__sizeof__ 		<type 'builtin_function_or_method'>
__str__ 		<type 'method-wrapper'>
__subclasshook__ 		<type 'builtin_function_or_method'>
__weakref__ 		<type 'weakref'>
_check_if_worker_is_up_to_date 		<type 'instancemethod'>
_create_new_grid 		<type 'instancemethod'>
_handlers 		<type 'list'>
_local_options 		<type 'dict'>
_start 		<type 'instancemethod'>
_stop 		<type 'instancemethod'>
_stop_worker 		<class 'amuse.support.methods.CodeMethodWrapper'>
all_literature_references_string 		<type 'instancemethod'>
before_get_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
before_set_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
center_of_mass_position 		center_of_mass_position 		AttributeError
center_of_mass_velocity 		center_of_mass_velocity 		AttributeError
channel 		<class 'amuse.rfi.channel.MpiChannel'>
channel_factory 		<type 'type'>
channel_type 		<type 'str'>
cleanup_code 		<type 'instancemethod'>
commit_parameters 		<type 'instancemethod'>
commit_particles 		<class 'amuse.support.methods.CodeMethodWrapper'>
define_converter 		<type 'instancemethod'>
define_errorcodes 		<type 'instancemethod'>
define_methods 		<type 'instancemethod'>
define_parameters 		<type 'instancemethod'>
define_particle_sets 		<type 'instancemethod'>
define_properties 		<type 'instancemethod'>
define_state 		<type 'instancemethod'>
delete_particle 		<class 'amuse.support.methods.CodeMethodWrapper'>
disable_stopping_condition 		<class 'amuse.support.methods.CodeMethodWrapper'>
dm_particles 		<class 'amuse.datamodel.particles.Particles'>
enable_stopping_condition 		<class 'amuse.support.methods.CodeMethodWrapper'>
ensure_stop_interface_at_exit 		<type 'instancemethod'>
errorcodes 		<type 'dict'>
evolve_model 		<class 'amuse.support.methods.CodeMethodWrapper'>
export2bibtex 		<type 'instancemethod'>
export2html 		<type 'instancemethod'>
gas_particles 		<class 'amuse.datamodel.particles.Particles'>
get_acc_tstp 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_acceleration 		get_acceleration 		Exception
get_adaptive_eps 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_alpha 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_balsara 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_begin_time 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_beta 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_bh_tol 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_center_of_mass_position 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_center_of_mass_velocity 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_colliding_particles 		<type 'instancemethod'>
get_consph 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_consthsm 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_cool_par 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_cosmo 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_courant 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_crionrate 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_data_directory 		<type 'instancemethod'>
get_density 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_dinternal_energy_dt 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_directsum 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_dtime 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_eps 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_eps2 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_eps_is_h 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_epsgas 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_epssph 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_feedback 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_fi_data_directory 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_firstsnap 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_fixthalo 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_freea 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_freeaexp 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_freetstp 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_freev 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_freevexp 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_gamma 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_gdgop 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_gdgtol 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_graineff 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_gravity_at_point 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_halofile 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_handler 		<type 'instancemethod'>
get_heat_par1 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_heat_par2 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_hupdatemethod 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_hydro_state_at_point 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_id_of_removed_sph_particle 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_index_of_first_particle 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_index_of_next_particle 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_indices_of_colliding_particles 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_internal_energy 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_isotherm 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_kinetic_energy 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_mass 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_masscrit 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_max_tbin 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_minppbin 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_name_of_current_state 		<type 'instancemethod'>
get_nn_tol 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_nsmooth 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_nsmtol 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_number_of_particles 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_number_of_sph_particles_removed 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_number_of_stopping_conditions_set 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_optdepth 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_output_directory 		<type 'instancemethod'>
get_pboxsize 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_periodic_boundaries_flag 		<type 'instancemethod'>
get_position 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_potential 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_potential_at_point 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_potential_energy 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_pressure 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_radiate 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_radius 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_removgas 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_rhomax 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_selfgrav 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_sfeff 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_sfmode 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_smoothing_length 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_smoothinput 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_sne_eff 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_sph_visc 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_sphinit 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_sqrttstp 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_star_tform 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_starform 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_state 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_state_sph 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_state_star 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_steplog 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_stepout 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_stopping_condition_info 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_stopping_condition_maximum_density_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_stopping_condition_maximum_internal_energy_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_stopping_condition_minimum_density_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_stopping_condition_minimum_internal_energy_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_stopping_condition_number_of_steps_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_stopping_condition_out_of_box_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_stopping_condition_particle_index 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_stopping_condition_timeout_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_targetnn 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_tbubble 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_tcollfac 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_thermal_energy 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_time 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_time_step 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_total_energy 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_total_mass 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_total_radius 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_tsnbeg 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_tstepcrit 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_tstpcr2 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_uentropy 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_unitl_in_kpc 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_unitm_in_msun 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_use_hydro 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_usequad 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_velocity 		<class 'amuse.support.methods.CodeMethodWrapper'>
get_verbosity 		<class 'amuse.support.methods.CodeMethodWrapper'>
has_stopping_condition 		<class 'amuse.support.methods.CodeMethodWrapper'>
hasoption 		<type 'instancemethod'>
initialize_code 		<type 'instancemethod'>
input_data_root_directory 		<type 'str'>
instances 		<type 'list'>
internal__get_message_polling_interval 		<class 'amuse.support.methods.CodeMethodWrapper'>
internal__set_message_polling_interval 		<class 'amuse.support.methods.CodeMethodWrapper'>
invoke_state_change 		<type 'instancemethod'>
is_stop_interfaces_registered 		<type 'bool'>
is_stopping_condition_enabled 		<class 'amuse.support.methods.CodeMethodWrapper'>
is_stopping_condition_set 		<class 'amuse.support.methods.CodeMethodWrapper'>
iter_options 		<type 'instancemethod'>
kinetic_energy 		<class 'amuse.units.quantities.ScalarQuantity'>
legacy_doc 		<type 'str'>
legacy_interface 		<class 'amuse.community.fi.interface.FiInterface'>
mode 		<type 'str'>
model_time 		<class 'amuse.units.quantities.ScalarQuantity'>
must_handle_state 		<type 'bool'>
must_start_worker 		<type 'bool'>
name_of_the_worker 		<type 'instancemethod'>
names_of_classes_with_references 		<type 'instancemethod'>
new_dm_particle 		<class 'amuse.support.methods.CodeMethodWrapper'>
new_particle 		<class 'amuse.support.methods.CodeMethodWrapper'>
new_sph_particle 		<class 'amuse.support.methods.CodeMethodWrapper'>
new_star_particle 		<class 'amuse.support.methods.CodeMethodWrapper'>
option_sections 		<type 'tuple'>
output_data_root_directory 		<type 'str'>
overridden 		<type 'instancemethod'>
parameters 		<class 'amuse.datamodel.parameters.ParametersWithDocs'>
particles 		<class 'amuse.datamodel.particles.ParticlesSuperset'>
polling_interval_in_milliseconds 		<type 'int'>
potential_energy 		<class 'amuse.units.quantities.ScalarQuantity'>
print_literature_references 		<type 'instancemethod'>
recommit_parameters 		<class 'amuse.support.methods.CodeMethodWrapper'>
recommit_particles 		<class 'amuse.support.methods.CodeMethodWrapper'>
redirect_file 		<type 'str'>
redirect_stderr_file 		<type 'str'>
redirect_stdout_file 		<type 'str'>
redirection 		<type 'str'>
redirection_filenames 		<type 'list'>
register_use 		<type 'instancemethod'>
reset 		<type 'instancemethod'>
set_acc_tstp 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_acceleration 		<type 'NoneType'>
set_adaptive_eps 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_alpha 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_balsara 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_begin_time 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_beta 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_bh_tol 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_consph 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_consthsm 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_cool_par 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_cosmo 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_courant 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_crionrate 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_directsum 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_dtime 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_eps 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_eps2 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_eps_is_h 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_epsgas 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_epssph 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_feedback 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_fi_data_directory 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_firstsnap 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_fixthalo 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_freea 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_freeaexp 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_freetstp 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_freev 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_freevexp 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_gamma 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_gdgop 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_gdgtol 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_graineff 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_halofile 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_heat_par1 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_heat_par2 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_hupdatemethod 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_internal_energy 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_isotherm 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_mass 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_masscrit 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_max_tbin 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_minppbin 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_nn_tol 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_nsmooth 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_nsmtol 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_optdepth 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_pboxsize 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_position 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_radiate 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_radius 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_removgas 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_rhomax 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_selfgrav 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_sfeff 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_sfmode 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_smoothing_length 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_smoothinput 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_sne_eff 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_sph_visc 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_sphinit 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_sqrttstp 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_star_tform 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_starform 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_state 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_state_sph 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_state_star 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_steplog 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_stepout 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_stopping_condition_maximum_density_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_stopping_condition_maximum_internal_energy_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_stopping_condition_minimum_density_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_stopping_condition_minimum_internal_energy_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_stopping_condition_number_of_steps_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_stopping_condition_out_of_box_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_stopping_condition_timeout_parameter 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_targetnn 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_tbubble 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_tcollfac 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_time_step 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_tsnbeg 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_tstepcrit 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_tstpcr2 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_uentropy 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_unitl_in_kpc 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_unitm_in_msun 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_use_hydro 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_usequad 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_velocity 		<class 'amuse.support.methods.CodeMethodWrapper'>
set_verbosity 		<class 'amuse.support.methods.CodeMethodWrapper'>
setup 		<type 'instancemethod'>
setup_particles 		<type 'instancemethod'>
star_particles 		<class 'amuse.datamodel.particles.Particles'>
state_machine 		<class 'amuse.support.state.StateMachine'>
stop 		<class 'amuse.support.methods.CodeMethodWrapper'>
stopping_conditions 		<class 'amuse.support.codes.stopping_conditions.StoppingConditions'>
synchronize_model 		<class 'amuse.support.methods.CodeMethodWrapper'>
thermal_energy 		<class 'amuse.units.quantities.ScalarQuantity'>
total_energy 		<class 'amuse.units.quantities.ScalarQuantity'>
total_mass 		total_mass 		AttributeError
total_radius 		total_radius 		AttributeError
unit_converter 		<class 'amuse.units.nbody_system.nbody_to_si'>
update_particle_set 		<type 'instancemethod'>
update_particles 		<type 'instancemethod'>
use_modules 		<type 'list'>

