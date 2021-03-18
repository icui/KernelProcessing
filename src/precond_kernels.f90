! Ebru1: The choice of THRESHOLD value is somewhat subjective. It is not trivial to set it like the 20% of max value
! which may be OK for smaller scale studies but global scale needs a few trial&error to adjust this parameter for
! every iteration. Needs some more investigation..

! Ebru2: I find the preconditioner behave better after changing the order of smoothing and preconditioning in
! post-processing upon the suggestion by Ryan & Yanhua.
! However, I am still not convinced by Ryan's latest suggestion that preconditioner should be smoothed more than the
! gradients of other parameters that the preconditioner to be applied. I currently smooth the preconditioner and
! the other gradients in the same way.

module precond_kernels_sub
  use mpi
  use global, only : max_all_all_cr, min_all_all_cr, CUSTOM_REAL, exit_mpi, &
                     myrank, quantile_all_all_cr
  implicit none

  contains

  subroutine get_sys_args(input_kernel, input_hess, output_kernel, threshold_hess)
    character(len=*), intent(inout) :: input_kernel, input_hess, output_kernel
    real(kind=CUSTOM_REAL), intent(inout) :: threshold_hess

    character(len=20) :: threshold_str

    call getarg(1, input_kernel)
    call getarg(2, input_hess)
    call getarg(3, output_kernel)
    call getarg(4, threshold_str)

    read(threshold_str, *) threshold_hess

    if(input_kernel == '' .or. input_hess == '' .or. output_kernel == '') then
      call exit_mpi("Usage: xprecond_kernels input_kernel input_hess input_model output_kernel threshold_hess")
    endif

    if(myrank == 0) then
      write(*, *) "Input kernel: ", trim(input_kernel)
      write(*, *) "Input hessian: ", trim(input_hess)
      write(*, *) "Output kernel: ", trim(output_kernel)
      write(*, *) "Threshold hessian: ", threshold_hess
    endif

  end subroutine get_sys_args

  subroutine prepare_hessian(hess, threshold, hess_inv)
    real(CUSTOM_REAL), dimension(:, :, :, :), intent(inout) :: hess, hess_inv
    real(CUSTOM_REAL), intent(in) :: threshold

    real(kind=CUSTOM_REAL):: maxh_all, minh_all, damp

    hess = abs(hess)
    call max_all_all_cr(maxval(hess), maxh_all)
    call min_all_all_cr(minval(hess), minh_all)

    if ( maxh_all < 1.e-18 ) then
      call exit_mpi("hess max value < 1.e-18")
    endif

    if (myrank==0) then
      write(*, *) "Max and Min of hess: ", maxh_all, minh_all
    endif

    ! normalized hess
    damp = maxh_all * threshold
    hess = (hess + damp) / (maxh_all + damp)
    
    call max_all_all_cr(maxval(hess), maxh_all)
    call min_all_all_cr(minval(hess), minh_all)

    if (myrank==0) then
      write(*, *) 'min and max hess after norm', minh_all, maxh_all
      write(*, *) "Hessian condition number: ", threshold
    endif

    hess_inv = 1.0_CUSTOM_REAL / hess
    
  end subroutine prepare_hessian

end module precond_kernels_sub

program precond_kernels
  use mpi
  use adios_read_mod
  use AdiosIO
  use global, only : NGLLX, NGLLY, NGLLZ, NSPEC, myrank, CUSTOM_REAL
  use global, only : init_mpi
  use precond_kernels_sub

  implicit none

  character(len=500), parameter :: kernel_names(5) = &
    (/character(len=500) :: "bulk_betah_kl_crust_mantle", "bulk_betav_kl_crust_mantle", &
                            "bulk_c_kl_crust_mantle", "eta_kl_crust_mantle", "rho_kl_crust_mantle"/)

  character(len=500), parameter :: hess_names(4) = &
    (/character(len=500) :: "hess_vs_kl_crust_mantle", "hess_vp_kl_crust_mantle", "hess_eta_kl_crust_mantle", &
                            "hess_rho_kl_crust_mantle"/)

  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, 3):: hess = 0.0
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, 5):: kernels = 0.0, kernels_precond = 0.0, hess_out = 0.0
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC):: hess_inv

  character(len=500) :: input_kernel, input_hess, output_kernel
  real(kind=CUSTOM_REAL) :: threshold_hess
  integer:: ier

  call init_mpi()

  call get_sys_args(input_kernel, input_hess, output_kernel, threshold_hess)
  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, "verbose=1", ier)

  call read_bp_file_real(input_kernel, kernel_names, kernels)
  call read_bp_file_real(input_hess, hess_names, hess)

  call prepare_hessian(hess(:, :, :, :, 1), threshold_hess, kernels_precond(:, :, :, :, 1))
  ! kernels_precond(:, :, :, :, 1) = kernels(:, :, :, :, 1) * hess_inv
  ! kernels_precond(:, :, :, :, 2) = kernels(:, :, :, :, 2) * hess_inv

  call prepare_hessian(hess(:, :, :, :, 2), threshold_hess, hess_inv)
  kernels_precond(:, :, :, :, 3) = kernels(:, :, :, :, 3) * hess_inv


  call prepare_hessian(hess(:, :, :, :, 3), threshold_hess, hess_inv)
  kernels_precond(:, :, :, :, 4) = kernels(:, :, :, :, 4) * hess_inv

  call prepare_hessian(hess(:, :, :, :, 4), threshold_hess, hess_inv)
  kernels_precond(:, :, :, :, 5) = kernels(:, :, :, :, 5) * hess_inv

  call write_bp_file(kernels_precond, kernel_names, "KERNEL_GOURPS", output_kernel)

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program precond_kernels
