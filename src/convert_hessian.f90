! Ebru1: The choice of THRESHOLD value is somewhat subjective. It is not trivial to set it like the 20% of max value
! which may be OK for smaller scale studies but global scale needs a few trial&error to adjust this parameter for
! every iteration. Needs some more investigation..

! Ebru2: I find the preconditioner behave better after changing the order of smoothing and preconditioning in
! post-processing upon the suggestion by Ryan & Yanhua.
! However, I am still not convinced by Ryan's latest suggestion that preconditioner should be smoothed more than the
! gradients of other parameters that the preconditioner to be applied. I currently smooth the preconditioner and
! the other gradients in the same way.

module convert_hessian_sub
  use mpi
  use global, only : max_all_all_cr, min_all_all_cr, CUSTOM_REAL, exit_mpi, &
                     myrank, quantile_all_all_cr
  implicit none

  contains

  subroutine get_sys_args(input_hess, input_model, output_hess)
    character(len=*), intent(inout) :: input_hess, input_model

    call getarg(1, input_hess)
    call getarg(2, input_model)
    call getarg(3, output_hess)

    if(input_hess == '' .or. input_model == '' .or. output_hess == '') then
      call exit_mpi("Usage: xconvert_hessian input_hess input_model output_hess")
    endif

    if(myrank == 0) then
      write(*, *) "Input hessian: ", trim(input_hess)
      write(*, *) "Input model: ", trim(input_model)
      write(*, *) "Output hessian: ", trim(output_hess)
    endif

  end subroutine get_sys_args

end module convert_hessian_sub

program convert_hessian
  use mpi
  use adios_read_mod
  use AdiosIO
  use global, only : NGLLX, NGLLY, NGLLZ, NSPEC, myrank, CUSTOM_REAL
  use global, only : init_mpi
  use convert_hessian_sub

  implicit none

  character(len=500), parameter :: hess_names(3) = &
    (/character(len=500) :: "hess_kappa_kl_crust_mantle", "hess_mu_kl_crust_mantle", "hess_rho_kl_crust_mantle"/)

  character(len=500), parameter :: hess2_names(4) = &
    (/character(len=500) :: "hess_vs_kl_crust_mantle", "hess_vp_kl_crust_mantle", "hess_eta_kl_crust_mantle", &
                            "hess_rho_kl_crust_mantle"/)

  character(len=500), parameter :: model_names(3) = &
    (/character(len=500) :: "reg1/kappavstore", "reg1/muvstore", "reg1/muhstore"/)

  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, 3):: hess = 0.0, models = 0.0
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, 4):: hess2
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC):: kappa, mu
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC):: hess_kappa, hess_mu, hess_rho, hess_vp, hess_vs, hess_eta

  character(len=500) :: input_hess, input_model, output_hess
  integer:: ier

  call init_mpi()

  call get_sys_args(input_hess, input_model, output_hess)
  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, "verbose=1", ier)

  call read_bp_file_real(input_hess, hess_names, hess)
  call read_bp_file_real(input_model, model_names, models)

  kappa = models(:, :, :, :, 1)
  mu = sqrt((2.0*models(:, :, :, :, 2)**2 + models(:, :, :, :, 3)**2) / 3.0)

  hess_kappa = hess(:, :, :, :, 1)
  hess_mu = hess(:, :, :, :, 2)
  hess_rho = hess(:, :, :, :, 3)
  
  hess_vp = 4.0 * (1.0 + 4.0/3.0*(mu/kappa))**2 * hess_kappa
  hess_vs = 4.0 * hess_mu + 64.0 / 9.0 * (mu/kappa)**2 * hess_kappa
  hess_eta = hess_kappa + 4.0 / 9.0 * hess_mu

  hess2(:, :, :, :, 1) = hess_vs
  hess2(:, :, :, :, 2) = hess_vp
  hess2(:, :, :, :, 3) = hess_eta
  hess2(:, :, :, :, 4) = hess_rho

  call write_bp_file(hess2, hess2_names, "KERNEL_GOURPS", output_hess)

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program convert_hessian
