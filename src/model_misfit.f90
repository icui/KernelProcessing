module misfit_subs
  use global, only : CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC, myrank, init_mpi
  implicit none

  integer, parameter :: nvars = 1
  character(len=500), dimension(nvars), parameter :: model_names = &
    (/character(len=500) :: "bulk_betah_kl_crust_mantle"/)
  character(len=500), dimension(nvars), parameter :: model_names2 = &
    (/character(len=500) :: "reg1/vsv"/)
  character(len=500), dimension(nvars), parameter :: sponge_names = &
    (/character(len=500) :: "reg1/spongestore"/)

  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, nvars) :: ref_model, &
                                                                          new_model, &
                                                                          sponge
  ! 6 parameter perturbation + 5 extra perturbation
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, nvars) :: perturb_model

  contains
  subroutine get_sys_args(ref_model_file, new_model_file, solver_file, sponge_file)
    use global, only : exit_mpi
    implicit none
    character(len=500), intent(in) :: ref_model_file, new_model_file, solver_file, sponge_file

    call getarg(1, ref_model_file)
    call getarg(2, new_model_file)
    call getarg(3, solver_file)
    call getarg(4, sponge_file)

    if(trim(ref_model_file) == '' .or. trim(new_model_file) == '') then
      call exit_mpi('Usage: xmodel_misfit ref_model_file new_model_file solver_data')
    endif

  end subroutine get_sys_args

end module misfit_subs

program main

  use mpi
  use adios_read_mod
  use global, only : CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC, myrank, init_mpi, &
    Parallel_ComputeL2normSquare
  use AdiosIO
  use misfit_subs
  implicit none

  character(len=500) :: ref_model_file, new_model_file, solver_file, sponge_file
  real(kind=CUSTOM_REAL) :: model_misfit
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: jacobian

  integer :: ier

  call init_mpi()

  call get_sys_args(ref_model_file, new_model_file, solver_file, sponge_file)
 
  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                                  "verbose=1", ier)


  if (trim(ref_model_file) == '_') then
    call read_bp_file_real(new_model_file, model_names, ref_model)
    perturb_model = ref_model
  else
    call read_bp_file_real(ref_model_file, model_names2, ref_model)
    call read_bp_file_real(new_model_file, model_names2, new_model)
    perturb_model = (ref_model - new_model)
  endif

  call read_bp_file_real(sponge_file, sponge_names, sponge)
  perturb_model = perturb_model * sponge

  call calculate_jacobian_matrix(solver_file, jacobian)
  call Parallel_ComputeL2normSquare(perturb_model, 1, jacobian, model_misfit)

  ! call write_bp_file(perturb_model, model_names, "KERNELS_GROUP", "sp.bp")

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

  if(myrank == 0) print *, "Model misfit:", model_misfit

end program main
