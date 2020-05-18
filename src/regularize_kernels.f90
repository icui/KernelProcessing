! Ebru1: The choice of THRESHOLD value is somewhat subjective. It is not trivial to set it like the 20% of max value
! which may be OK for smaller scale studies but global scale needs a few trial&error to adjust this parameter for
! every iteration. Needs some more investigation..

! Ebru2: I find the preconditioner behave better after changing the order of smoothing and preconditioning in
! post-processing upon the suggestion by Ryan & Yanhua.
! However, I am still not convinced by Ryan's latest suggestion that preconditioner should be smoothed more than the
! gradients of other parameters that the preconditioner to be applied. I currently smooth the preconditioner and
! the other gradients in the same way.

module regularize_kernels_sub
  use mpi
  use global, only : myrank, nprocs, NGLLX, NGLLY, NGLLZ, NSPEC, NGLOB, CUSTOM_REAL, exit_mpi, &
                     max_all_all_cr
  use AdiosIO

  implicit none

  ! ======================================================
  ! MODELS
  integer, parameter :: NMODELS = 6
  character(len=500), dimension(NMODELS), parameter :: model_names = &
    (/character(len=500) :: "reg1/vpv", "reg1/vph", "reg1/vsv", &
                            "reg1/vsh", "reg1/eta", "reg1/rho"/)
  integer, parameter :: vpv_idx=1, vph_idx=2, vsv_idx=3, vsh_idx=4, &
                        eta_idx=5, rho_idx=6
  character(len=500), dimension(NMODELS), parameter :: model_perturb_names = &
    (/character(len=150) :: "reg1/dvpvvpv","reg1/dvphvph","reg1/dvsvvsv", &
                            "reg1/dvshvsh","reg1/detaeta","reg1/drhorho"/)
  ! transverse isotropic model files
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NMODELS) :: models = 0.0

  ! ======================================================
  ! KERNELS

  integer, parameter :: NKERNELS = 6    !bulk_betah, bulk_betav, bulk_c, eta
  character(len=500), parameter :: kernel_names(NKERNELS) = &
    (/character(len=500) :: "hess_kl_crust_mantle", "bulk_betah_kl_crust_mantle", &
                            "bulk_betav_kl_crust_mantle", "bulk_c_kl_crust_mantle", &
                            "eta_kl_crust_mantle", "rho_kl_crust_mantle"/)
  integer, parameter :: hess_idx = 1, betah_kl_idx = 2, betav_kl_idx = 3, &
                        bulk_c_kl_idx = 4, eta_kl_idx = 5, rho_kl_idx = 6

  real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: kernels = 0.0, &
                                                                          kernels_damp = 0.0

  contains

  subroutine get_sys_args(input_file, input_model, output_file, step_fac)
    character(len=*), intent(inout) :: input_file, input_model, output_file
    real(kind=CUSTOM_REAL), intent(inout) :: step_fac

    character(len=20) :: step_fac_str

    call getarg(1, input_file)
    call getarg(2, input_model)
    call getarg(3, output_file)
    call getarg(4, step_fac_str)

    read(step_fac_str, *) step_fac

    if(input_file == '' .or. input_model == '' .or. output_file == '') then
      call exit_mpi("Usage: xregularize_kernels input_kernel input_model output_kernel")
    endif

    if(myrank == 0) then
      write(*, *) "Input kernel: ", trim(input_file)
      write(*, *) "Input model: ", trim(input_model)
      write(*, *) "Output kernel: ", trim(output_file)
      write(*, *) "Regularization factor: ", step_fac
    endif

  end subroutine get_sys_args

  subroutine regularize_kernel(step_fac)
    use global , only : FOUR_THIRDS
    ! DMP regularization:
    ! J = J0 + step_fac * ||m||^2
    ! K = K0 + step_fac * m
    ! H = H0 + step_fac

    real(kind=CUSTOM_REAL), intent(inout) :: step_fac
    real(kind=CUSTOM_REAL):: maxh_all, step_len

    integer :: i, j, k, ispec
    real(kind=CUSTOM_REAL) :: alphav, alphah, betav, betah, eta, rho, bulk_c
    real(kind=CUSTOM_REAL) :: betav_kl, betah_kl, bulk_c_kl, eta_kl, rho_kl

    call max_all_all_cr(maxval(abs(kernels(:, :, :, :, hess_idx))), maxh_all)

    step_len = maxh_all * step_fac

    if(myrank == 0) then
      write(*, *) "Regularization parameter:", step_len
    endif

    kernels_damp(:, :, :, :, hess_idx) = kernels(:, :, :, :, hess_idx) + step_len

    do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            ! initial model values
            alphav = models(i,j,k,ispec,vpv_idx)
            alphah = models(i,j,k,ispec,vph_idx)
            betav  = models(i,j,k,ispec,vsv_idx)
            betah  = models(i,j,k,ispec,vsh_idx)
            eta    = models(i,j,k,ispec,eta_idx)
            rho    = models(i,j,k,ispec,rho_idx)
            bulk_c = sqrt(alphav ** 2 - FOUR_THIRDS * betav ** 2)

            ! initial kernel values
            betav_kl  = kernels(i,j,k,ispec,betav_kl_idx)
            betah_kl  = kernels(i,j,k,ispec,betah_kl_idx)
            bulk_c_kl = kernels(i,j,k,ispec,bulk_c_kl_idx)
            eta_kl    = kernels(i,j,k,ispec,eta_kl_idx)
            rho_kl    = kernels(i,j,k,ispec,rho_kl_idx)

            ! regularized kernel values
            kernels_damp(i,j,k,ispec,betav_kl_idx) = betav_kl + step_len * betav
            kernels_damp(i,j,k,ispec,betah_kl_idx) = betah_kl + step_len * betah
            kernels_damp(i,j,k,ispec,bulk_c_kl_idx) = bulk_c_kl + step_len * bulk_c
            kernels_damp(i,j,k,ispec,eta_kl_idx) = eta_kl + step_len * eta
            kernels_damp(i,j,k,ispec,rho_kl_idx) = rho_kl + step_len * rho
          enddo
        enddo
      enddo
    enddo

  end subroutine regularize_kernel

end module regularize_kernels_sub

program regularize_kernels
  use mpi
  use adios_read_mod
  use AdiosIO
  use global, only : init_mpi
  use regularize_kernels_sub

  implicit none

  character(len=500) :: input_file, input_model, output_file
  real(kind=CUSTOM_REAL) :: step_fac
  integer:: ier

  call init_mpi()

  if(trim(kernel_names(hess_idx)) /= "hess_kl_crust_mantle") then
    call exit_mpi("hess_idx is wrong!")
  endif

  call get_sys_args(input_file, input_model, output_file, step_fac)

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  call read_bp_file_real(input_file, kernel_names, kernels)
  call read_bp_file_real(input_model, model_names, models)

  ! apply DMP to kernel and Hessian
  call regularize_kernel(step_fac)

  call write_bp_file(kernels_damp, kernel_names, "KERNEL_GOURPS", output_file)

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program regularize_kernels
