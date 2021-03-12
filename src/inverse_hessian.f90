! Ebru1: The choice of THRESHOLD value is somewhat subjective. It is not trivial to set it like the 20% of max value
! which may be OK for smaller scale studies but global scale needs a few trial&error to adjust this parameter for
! every iteration. Needs some more investigation..

! Ebru2: I find the preconditioner behave better after changing the order of smoothing and preconditioning in
! post-processing upon the suggestion by Ryan & Yanhua.
! However, I am still not convinced by Ryan's latest suggestion that preconditioner should be smoothed more than the
! gradients of other parameters that the preconditioner to be applied. I currently smooth the preconditioner and
! the other gradients in the same way.

module inverse_hessian_sub
  use mpi
  use global, only : max_all_all_cr, min_all_all_cr, CUSTOM_REAL, exit_mpi, &
                     myrank, quantile_all_all_cr
  implicit none

  contains

  subroutine get_sys_args(input_hess, input_model, output_hess, threshold_hess)
    character(len=*), intent(inout) :: input_kernel, input_hess, input_model, output_kernel
    real(kind=CUSTOM_REAL), intent(inout) :: threshold_hess

    character(len=20) :: threshold_str

    call getarg(1, input_hess)
    call getarg(2, input_model)
    call getarg(3, output_hess)
    call getarg(4, threshold_str)

    read(threshold_str, *) threshold_hess

    if(input_hess == '' .or. input_model == '' .or. output_hess == '') then
      call exit_mpi("Usage: xprecond_kernels input_hess input_model output_hess threshold_hess")
    endif

    if(myrank == 0) then
      write(*, *) "Input hessian: ", trim(input_hess)
      write(*, *) "Input model: ", trim(input_model)
      write(*, *) "Output hessian: ", trim(output_hess)
      write(*, *) "Threshold hessian: ", threshold_hess
    endif

  end subroutine get_sys_args

  subroutine prepare_hessian(hess, threshold, hess_inv)
    real(CUSTOM_REAL), dimension(:, :, :, :), intent(inout) :: hess, hess_inv
    real(CUSTOM_REAL), intent(in) :: threshold

    real(kind=CUSTOM_REAL):: maxh_all, minh_all, cutoff

    call max_all_all_cr(maxval(hess), maxh_all)
    call min_all_all_cr(minval(hess), minh_all)

    if ( maxh_all < 1.e-18 ) then
      call exit_mpi("hess max value < 1.e-18")
    endif

    if (myrank==0) then
      write(*, *) "Max and Min of hess: ", maxh_all, minh_all
    endif

    ! normalized hess
    hess = abs(hess) / maxh_all

    call quantile_all_all_cr(hess, threshold, cutoff)
    call max_all_all_cr(maxval(hess), maxh_all)
    call min_all_all_cr(minval(hess), minh_all)

    if (myrank==0) then
      write(*, *) 'min and max hess after norm', minh_all, maxh_all
      write(*, *) "Hessian Threshold quantile: ", threshold
      write(*, *) 'Hessian Threshold: ', cutoff
    endif

    where(hess > cutoff )
      hess_inv = 1.0_CUSTOM_REAL / hess
    elsewhere
      hess_inv = 1.0_CUSTOM_REAL / cutoff
    endwhere
  end subroutine prepare_hessian

end module inverse_hessian_sub

program precond_kernels
  use mpi
  use adios_read_mod
  use AdiosIO
  use global, only : NGLLX, NGLLY, NGLLZ, NSPEC, myrank, CUSTOM_REAL
  use global, only : init_mpi
  use inverse_hessian_sub

  implicit none

  character(len=500), parameter :: hess_names(3) = &
    (/character(len=500) :: "hess_kappa_kl_crust_mantle", "hess_mu_kl_crust_mantle", "hess_rho_kl_crust_mantle"/)

  character(len=500), parameter :: invhess_names(4) = &
    (/character(len=500) :: "invhess_vp_crust_mantle", "invhess_vs_crust_mantle", "invhess_eta_crust_mantle", "invhess_rho_crust_mantle"/)

  character(len=500), parameter :: model_names(3) = &
    (/character(len=500) :: "reg1/kappavstore", "reg1/muvstore", "reg1/muhstore"/)

  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, 3):: hess = 0.0, models = 0.0
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, 4):: invhess
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC):: hess_kappa, hess_mu, hess_rho, hess_vp, hess_vs, hess_eta, hess_inv

  character(len=500) :: input_hess, input_model, output_hess
  real(kind=CUSTOM_REAL) :: threshold_hess
  integer:: ier

  call init_mpi()

  call get_sys_args(input_hess, input_model, output_hess, threshold_hess)
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

  call prepare_hessian(hess_vs, threshold_hess, hess_inv)
  invhess(:, :, :, :, 1) = hess_inv

  call prepare_hessian(hess_vp, threshold_hess, hess_inv)
  invhess(:, :, :, :, 2) = hess_inv

  call prepare_hessian(hess_eta, threshold_hess, hess_inv)
  invhess(:, :, :, :, 3) = hess_inv

  call prepare_hessian(hess_rho, threshold_hess, hess_inv)
  invhess(:, :, :, :, 4) = hess_inv

  call write_bp_file(invhess, invhess_names, "KERNEL_GOURPS", output_hess)

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program precond_kernels
