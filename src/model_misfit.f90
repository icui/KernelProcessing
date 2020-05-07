module perturb_subs
  use global, only : CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC, myrank, init_mpi
  implicit none

  integer, parameter :: nvars = 6
  character(len=500), dimension(nvars), parameter :: model_names = &
    (/character(len=500) :: "reg1/vpv", "reg1/vph", "reg1/vsv", &
                            "reg1/vsh", "reg1/eta", "reg1/rho"/)
  integer, parameter :: vpv_idx = 1, vph_idx = 2, vsv_idx = 3, &
                        vsh_idx = 4, eta_idx = 5, rho_idx = 6

  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, nvars) :: ref_model, &
                                                                          new_model
  ! 6 parameter perturbation + 5 extra perturbation
  ! don't change the order unless you know what you are doing
  character(len=500), dimension(6), parameter :: perturb_names = &
    (/character(len=500) :: "reg1/dvpvvpv", "reg1/dvphvph", "reg1/dvsvvsv", &
                            "reg1/dvshvsh", "reg1/detaeta", "reg1/drhorho"/)

  contains
  subroutine get_sys_args(ref_model_file, new_model_file)
    use global, only : myrank
    use global, only : exit_mpi
    implicit none
    character(len=500), intent(in) :: ref_model_file, new_model_file

    call getarg(1, ref_model_file)
    call getarg(2, new_model_file)

    if(trim(ref_model_file) == '' .or. trim(new_model_file) == '') then
      call exit_mpi('Usage: xmodel_perturbs ref_model_file new_model_file')
    endif

    if (myrank == 0) then
      print*, "ref model file: ", trim(ref_model_file)
      print*, "new model file: ", trim(new_model_file)
    endif
  end subroutine get_sys_args

end module perturb_subs

program main

  use mpi
  use adios_read_mod
  use global, only : CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC, myrank, init_mpi
  use AdiosIO
  use perturb_subs
  implicit none

  character(len=500) :: ref_model_file, new_model_file

  integer :: ier

  call init_mpi()

  if(myrank == 0) print*, "mpi done"

  call get_sys_args(ref_model_file, new_model_file)
 
  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                                  "verbose=1", ier)

  if(myrank == 0) print*, "Read ref model file"
  call read_bp_file_real(ref_model_file, model_names, ref_model)
  if(myrank == 0) print*, "Read new model file"
  call read_bp_file_real(new_model_file, model_names, new_model)
  if(myrank == 0) print*, "Done reading"

  perturb_model(:,:,:,:,1:6) = new_model - ref_model

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

  if(myrank == 0) print*, "Job finished"

end program main
