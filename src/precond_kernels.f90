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

  subroutine get_sys_args(input_kernel, input_hess, output_kernel)
    character(len=*), intent(inout) :: input_kernel, input_hess, output_kernel

    call getarg(1, input_kernel)
    call getarg(2, input_hess)
    call getarg(4, output_kernel)

    if(input_kernel == '' .or. input_hess == '' .or. output_kernel == '') then
      call exit_mpi("Usage: xprecond_kernels input_kernel input_hess output_kernel")
    endif

    if(myrank == 0) then
      write(*, *) "Input kernel: ", trim(input_kernel)
      write(*, *) "Input hessian: ", trim(input_hess)
      write(*, *) "Output kernel: ", trim(output_kernel)
    endif

  end subroutine get_sys_args

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

  character(len=500), parameter :: invhess_names(4) = &
    (/character(len=500) :: "invhess_vs_kl_crust_mantle", "invhess_vp_kl_crust_mantle", "invhess_eta_kl_crust_mantle", &
                            "invhess_rho_kl_crust_mantle"/)

  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, 3):: invhess = 0.0
  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, 5):: kernels = 0.0, kernels_precond = 0.0

  character(len=500) :: input_kernel, input_hess, output_kernel
  integer:: ier

  call init_mpi()

  call get_sys_args(input_kernel, input_hess, output_kernel)
  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, "verbose=1", ier)

  call read_bp_file_real(input_kernel, kernel_names, kernels)
  call read_bp_file_real(input_hess, invhess_names, invhess)

  ! betah, betav / vs
  kernels_precond(:, :, :, :, 1) = kernels(:, :, :, :, 1) * invhess(:, :, :, :, 1)
  kernels_precond(:, :, :, :, 2) = kernels(:, :, :, :, 2) * invhess(:, :, :, :, 1)

  ! bulkc / vp
  kernels_precond(:, :, :, :, 3) = kernels(:, :, :, :, 3) * invhess(:, :, :, :, 2)
  
  ! eta
  kernels_precond(:, :, :, :, 4) = kernels(:, :, :, :, 4) * invhess(:, :, :, :, 3)

  ! rho
  kernels_precond(:, :, :, :, 5) = kernels(:, :, :, :, 5) * invhess(:, :, :, :, 4)

  call write_bp_file(kernels_precond, kernel_names, "KERNEL_GOURPS", output_kernel)

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program precond_kernels
