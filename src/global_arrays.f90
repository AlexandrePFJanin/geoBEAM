module global_arrays

  !=========================
  ! Dynamic memory allocation for global arrays
  ! Module added by A.JANIN 13.02.2026
  !=========================

  implicit none

  !=========================
  ! Dimensions runtime
  !=========================
  integer :: npts_glob
  integer :: ndis_glob

  !=========================
  ! Solver global parameter
  !=========================
  integer :: niter_solv   ! what is the current iteration number of the soler
  integer :: layer_solv   ! what is the current depth if the multi-layer on which the solver is

  !=========================
  ! Boundary conditions
  !=========================
  integer, allocatable :: ISPACE(:)
  real(4), allocatable :: SPACE(:)

  !=========================
  ! Influence matrix and derivates
  !=========================
  real(4), allocatable :: XMATRIX(:,:) ! The matrix: influence coeff + extra (last) column for b.c.s
  real(4), allocatable :: AMATRIX(:,:) ! Original coeff matrix only, never shrinked or reshaped
  integer :: NUM_Ds_SAVED ! keep track of NUM_Ds, index of column where the b.c.s are written
  logical :: debug
  real(4), allocatable :: DVEC(:,:)    ! Vectors of *displacement* on each element (along-strike,down-dip,normal-out): returned by SOLVE (so after one iteration only)
  real(4), allocatable :: DVEC2(:,:)   ! Same as DVEC but ALL relative displacements have been set to 0. (usefull in friction, to not count twice the effect of KODE>=10)
  real(4), allocatable :: DVEC2a(:,:)  ! Same as DVEC2 but for the first sub-layer of the iteration (for frictional elements only) when the solver is in STATE 2
  real(4), allocatable :: DVEC2b(:,:)  ! Same as DVEC2 but for the second sub-layer of the iteration (for dyndike elements only) when the solver is in STATE 2
  real(4), allocatable :: DVEC0(:,:)   ! Reference DVEC: vector DVEC at the previous iteration - used during the convergence search
  real(4), allocatable :: DVECI(:,:) ! Displacement built iterativelly when frictional elements are introduced

  !=========================
  ! Stresses on each dislocations
  !=========================
  real(4), allocatable :: SSTORED(:,:) ! not released (only when friction)
  real(4), allocatable :: SDRIVER(:,:) ! released for motion
  logical :: is_sdriver_updated

  !=========================
  ! Transformation matrices
  !=========================
  real(4), allocatable :: UG2P(:,:,:)
  real(4), allocatable :: SG2P(:,:,:)

  !=========================
  ! Element descriptors
  !=========================
  real(4), allocatable :: ZCE(:) ! z coordinates (depths) of the central point of elements
  real(4), allocatable :: XO(:),YO(:),ZO(:)
  real(4), allocatable :: C(:),S(:),DIP(:)
  real(4), allocatable :: CDIP(:),SDIP(:)
  real(4), allocatable :: BWX1(:),BWX2(:)
  integer, allocatable :: NBX1(:),NBX2(:)

  !=========================
  ! Station coordinates
  !=========================
  real(8), allocatable :: xu(:),yu(:),zu(:)

  !=========================
  ! Frictional status of the elements
  !=========================
  logical :: any_frictional    ! True if there is any frictional element
  integer(4), allocatable :: element_fstatus(:)       ! globally, what is the frictional status of the element (if the elements slides once in the iteration process, will be set as "sliding" for ever): that's the definitive fstatus.
  integer(4), allocatable :: element_iter_fstatus(:)  ! frictional status of each element only from one iteration to another (so, one element can goes from sliding to locked to sliding again with the iterative process)
  !  -10: unknown
  !   -2: not frictional
  !   -1: frictional but not determined
  !    0: locked
  !    1: sliding

  !=========================
  ! Dynamical dike elements
  !=========================
  logical :: any_dyndike  ! True if there is any dynamical dike element
  integer(4), allocatable :: element_dstatus(:)       ! dynamical dike status (0: not a dyndike element, 1: is a dyndike element)
  real(4), allocatable :: hidden_nstress(:)           ! True normal stress on the dynamical dike elements
  real(4), allocatable :: hidden_nstress_ref(:)       ! reference profile for the iterative scheme. At each increment, the displacement will be computed from (hidden_nstress - hidden_nstress_ref), a.k.a "hidden_nstress" at the last iteration to take into account only the changse
  integer :: tmp_niter                                ! Iter number of the last update of SDRIVER for the dyndike element
  real(4), allocatable :: tmp_SDRIVER(:)              ! SDRIVER element at the last call for dyndike element

  !=========================
  ! Final solution status
  !=========================
  integer(4), allocatable:: codeStatus(:)

  !=========================
  ! Locked (backed up) fields
  !=========================
  ! backup arrays used and allocated only when locking dyndike element (SOLVER state 2).
  ! no need of backuping XMATRIX as it's reset at each start of SOLVE
  logical :: is_F_locked           ! Flag keeping track of the locking status of the frictional elements
  logical :: is_D_locked           ! Flag keeping track of the locking statis of the dyndike elements
  logical, allocatable :: is_locked(:)                ! is this element locked (whatever if it is a frictional or a dyndike element)
  integer, allocatable :: locked_KODE(:)              ! backup KODE
  integer, allocatable :: locked_i_kode(:)            ! backup i_kode (mod: input_array)
  integer, allocatable :: locked_i_fcode(:)           ! backup i_fcode
  integer(4), allocatable :: locked_element_fstatus(:)   ! backup element_fstatus
  integer, allocatable :: locked_element_dstatus(:)   ! backup element_dstatus
  real(4), allocatable :: locked_BC(:,:)              ! backup BC

contains



  subroutine allocate_global_arrays(npts,ndis)
    implicit none
    integer, intent(in) :: npts, ndis
    npts_glob = npts
    ndis_glob = ndis
    ! Reserved space for b.c. (boundary condition) codes, and b.c.s
    allocate(ISPACE(ndis))
    allocate(SPACE(3*ndis))
    ! Influence coefficient matrix; last column contains the 
    ! b.c.s and then the solution (relative displacements)
    allocate(XMATRIX(3*ndis,3*ndis+1))
    XMATRIX = 0.0
    ! Original coeff matrix
    allocate(AMATRIX(3*ndis,3*ndis))
    ! Displacement vector on each element
    allocate(DVEC(ndis,3))
    allocate(DVEC2(ndis,3))
    allocate(DVEC2a(ndis,3))
    allocate(DVEC2b(ndis,3))
    allocate(DVEC0(ndis,3))
    allocate(DVECI(ndis,3))
    DVEC   = 0.0   ! Init the displacement on all elements to 0.0
    DVEC2  = 0.0   ! Init at 0
    DVEC2b = 0.0   ! Init at 0
    DVEC0  = 0.0   ! Init at 0, no displacement
    DVECI  = 0.0   ! Init at 0
    ! Stresses on each dislocation
    allocate(SSTORED(ndis,3))
    allocate(SDRIVER(ndis,3))
    SSTORED = 0.0  ! Initialize all elements to 0.0
    SDRIVER = 0.0  ! Initialize all elements to 0.0
    ! Displacement and stress transformation matrices from
    ! global to planar (in plane) coordinates
    allocate(UG2P(3,3,ndis)) ! displacement transformation matrix for each plane 
    allocate(SG2P(3,6,ndis)) ! stress transformation matrix for each plane 
    ! Element descriptor parameters; for each plane
    allocate(ZCE(ndis))
    allocate(XO(ndis),YO(ndis),ZO(ndis)) ! X,Y,Z reference point in global coordinates
    allocate(C(ndis),S(ndis),DIP(ndis))  ! cosine & sine of the strike (cw wrt N)
    allocate(CDIP(ndis),SDIP(ndis))   ! dip (wrt horizontal), cosine & sine of the dip
    allocate(BWX1(ndis),BWX2(ndis))   ! sub-element widths in the strike & dip directions
    allocate(NBX1(ndis),NBX2(ndis))   ! number of sub-elements in the strike & dip directions
    ! user coordinates
    allocate(xu(npts),yu(npts),zu(npts))
    ! frictional status of elements
    allocate(element_fstatus(ndis))
    allocate(element_iter_fstatus(ndis))
    element_fstatus = -10        ! init on "unknown"
    element_iter_fstatus = -10   ! init on "unknown"
    ! final solution status when existing
    allocate(codeStatus(2)) ! codeStatus(1) = converged?, 1=True ,0=No
                            ! codeStatus(2) = final Niter
    ! dynamical diking
    allocate(element_dstatus(ndis))
    allocate(hidden_nstress(ndis))
    allocate(hidden_nstress_ref(ndis))
    hidden_nstress_ref = 0. ! init at 0
    allocate(tmp_SDRIVER(ndis))
    tmp_SDRIVER = 0. ! init at 0
    ! locking
    allocate(is_locked(ndis))
    is_locked = .FALSE.
  end subroutine allocate_global_arrays


  subroutine deallocate_global_arrays()
    implicit none
    deallocate(ISPACE)
    deallocate(SPACE)
    deallocate(XMATRIX)
    deallocate(AMATRIX)
    deallocate(DVEC)
    deallocate(DVEC2)
    deallocate(DVEC2a)
    deallocate(DVEC2b)
    deallocate(DVEC0)
    deallocate(DVECI)
    deallocate(SSTORED)
    deallocate(SDRIVER)
    deallocate(UG2P)
    deallocate(SG2P)
    deallocate(ZCE)
    deallocate(XO,YO,ZO)
    deallocate(C,S,DIP)
    deallocate(CDIP,SDIP)
    deallocate(BWX1,BWX2)
    deallocate(NBX1,NBX2)
    deallocate(xu,yu,zu)
    deallocate(element_fstatus)
    deallocate(element_iter_fstatus)
    deallocate(codeStatus)
    deallocate(element_dstatus)
    deallocate(hidden_nstress)
    deallocate(hidden_nstress_ref)
    deallocate(tmp_SDRIVER)
    deallocate(is_locked)
  end subroutine deallocate_global_arrays



  subroutine reset_xmatrix()
    implicit none
    ! -------------------------------------------------------
    ! Reset XMATRIX from AMATRIX
    ! -------------------------------------------------------
    XMATRIX(:, 1:3*ndis_glob) = AMATRIX
    XMATRIX(:, 3*ndis_glob+1) = 0.
  end subroutine reset_xmatrix



  subroutine update_DVEC()
    use global_inputs
    implicit none
    integer :: K, IBC
    ! -------------------------------------------------------
    ! Update the vector DVEC: To be put just before each RETURN in SOLVE
    ! and after the restoration of XMATRIX last vector (to restore the
    ! reorganization induced by SHRINK)
    ! -------------------------------------------------------
    DO K = 1, ndis_glob
        IBC = (K-1)*3
        DVEC(K,1) = XMATRIX(IBC+1,NUM_Ds_SAVED)
        DVEC(K,2) = XMATRIX(IBC+2,NUM_Ds_SAVED)
        DVEC(K,3) = XMATRIX(IBC+3,NUM_Ds_SAVED)
    ENDDO
    ! debug mode
    if (debug) then
      write(*,*) ''
      write(*,*) ' [DEBUG] (updated DVEC)'
      write(*,*) '         KODE,   ID,   Ds (cm),   Dd (cm),   Dn (cm)'
      DO K = 1, ndis_glob
        IBC = (K-1)*3
        write(*,*) 'DVEC:',K,i_kode(K),DVEC(K,1),DVEC(K,2),DVEC(K,3)
      ENDDO
    endif
  end subroutine update_DVEC



  subroutine update_DVEC2(i_kode, unlock)
    implicit none
    integer :: K, IBC
    integer :: i_kode(:)
    logical :: skip, i_unlock
    logical, optional :: unlock
    ! -------------------------------------------------------
    ! Update the vector DVEC2
    ! Args:
    !     i_kode (array, dtype=int): array of KODE for each element
    !     unlock (logical, optional): Condition to consider or not
    !             the displacement of the locked elements. If unlock
    !             is set to True, then, take into account the displs
    !             of the locked elements. (Default: unlock=True).
    ! -------------------------------------------------------
    ! Treatment of the optional input
    if (present(unlock)) then
      i_unlock = unlock
    else
      i_unlock = .True.
    endif
    ! Debug mode header
    if (debug) then
      write(*,*)
      write(*,*) ' [DEBUG] (updated DVEC2)'
      write(*,*) '  -> Consider the locked elements? ', i_unlock
    endif
    ! Iteration on each element
    DO K = 1, ndis_glob
        IBC = (K-1)*3
        DVEC2(K,1) = DVEC(K,1)
        DVEC2(K,2) = DVEC(K,2)
        DVEC2(K,3) = DVEC(K,3)
        ! Determine if the element BC adjutement should be skipped or not
        skip = .False.
        IF (is_locked(K)) THEN
          IF (i_unlock) THEN
            skip = .True.
          ENDIF
        ENDIF
        ! Application of BCs
        IF (.NOT. skip) THEN
          IF (i_kode(K).EQ.10) THEN
            DVEC2(K,1) = 0.
            DVEC2(K,2) = 0.
            DVEC2(K,3) = 0.
          ELSEIF (i_kode(K).EQ.11.OR.i_kode(K).EQ.12) THEN
            DVEC2(K,3) = 0.
          ELSEIF (i_kode(K).EQ.13) THEN
            DVEC2(K,2) = 0.
            DVEC2(K,3) = 0.
          ELSEIF (i_kode(K).EQ.14) THEN
            DVEC2(K,1) = 0.
            DVEC2(K,3) = 0.
          ELSEIF (i_kode(K).EQ.15) THEN
            DVEC2(K,1) = 0.
            DVEC2(K,2) = 0.
          ENDIF
        ELSE
          CONTINUE  ! if the element is locked, copy all DVEC in DVEC2 assuming it's all good
        ENDIF
        ! debug mode
        if (debug) then
          write(*,*) 'DVEC2:',K,DVEC2(K,1),DVEC2(K,2),DVEC2(K,3)
        endif
    ENDDO
  end subroutine update_DVEC2



  subroutine update_DVEC2a(i_kode)
    implicit none
    integer :: K, IBC
    integer :: i_kode(:)
    ! -------------------------------------------------------
    ! Update the vector DVEC2a (for the sub-layer (a) of the
    ! solver when in STATE 2).
    ! -------------------------------------------------------
    if (debug) then
      write(*,*)
      write(*,*) ' [DEBUG] (updated DVEC2a)'
    endif
    DO K = 1, ndis_glob
        IBC = (K-1)*3
        DVEC2a(K,1) = DVEC(K,1)
        DVEC2a(K,2) = DVEC(K,2)
        DVEC2a(K,3) = DVEC(K,3)
        IF (i_kode(K).EQ.10) THEN
          DVEC2a(K,1) = 0.
          DVEC2a(K,2) = 0.
          DVEC2a(K,3) = 0.
        ELSEIF (i_kode(K).EQ.11.OR.i_kode(K).EQ.12) THEN
          DVEC2a(K,3) = 0.
        ELSEIF (i_kode(K).EQ.13) THEN
          DVEC2a(K,2) = 0.
          DVEC2a(K,3) = 0.
        ELSEIF (i_kode(K).EQ.14) THEN
          DVEC2a(K,1) = 0.
          DVEC2a(K,3) = 0.
        ELSEIF (i_kode(K).EQ.15) THEN
          DVEC2a(K,1) = 0.
          DVEC2a(K,2) = 0.
        ENDIF
        ! debug mode
        if (debug) then
          write(*,*) 'DVEC2a:',K,DVEC2a(K,1),DVEC2a(K,2),DVEC2a(K,3)
        endif
    ENDDO
  end subroutine update_DVEC2a
  
  
  
  subroutine update_DVEC2b(i_kode)
    implicit none
    integer :: K, IBC
    integer :: i_kode(:)
    ! -------------------------------------------------------
    ! Update the vector DVEC2b (for the sub-layer (b) of the
    ! solver when in STATE 2).
    ! -------------------------------------------------------
    if (debug) then
      write(*,*)
      write(*,*) ' [DEBUG] (updated DVEC2b)'
    endif
    DO K = 1, ndis_glob
        IBC = (K-1)*3
        DVEC2b(K,1) = DVEC(K,1)
        DVEC2b(K,2) = DVEC(K,2)
        DVEC2b(K,3) = DVEC(K,3)
        IF (i_kode(K).EQ.10) THEN
          DVEC2b(K,1) = 0.
          DVEC2b(K,2) = 0.
          DVEC2b(K,3) = 0.
        ELSEIF (i_kode(K).EQ.11.OR.i_kode(K).EQ.12) THEN
          DVEC2b(K,3) = 0.
        ELSEIF (i_kode(K).EQ.13) THEN
          DVEC2b(K,2) = 0.
          DVEC2b(K,3) = 0.
        ELSEIF (i_kode(K).EQ.14) THEN
          DVEC2b(K,1) = 0.
          DVEC2b(K,3) = 0.
        ELSEIF (i_kode(K).EQ.15) THEN
          DVEC2b(K,1) = 0.
          DVEC2b(K,2) = 0.
        ENDIF
        ! debug mode
        if (debug) then
          write(*,*) 'DVEC2b:',K,DVEC2b(K,1),DVEC2b(K,2),DVEC2b(K,3)
        endif
    ENDDO
  end subroutine update_DVEC2b



  subroutine copy_DVEC2_DVEC2a()
    implicit none
    ! -------------------------------------------------------
    ! Copy DVEC2 into DVEC2a
    ! -------------------------------------------------------
    DVEC2a(:,1) = DVEC2(:,1)
    DVEC2a(:,2) = DVEC2(:,2)
    DVEC2a(:,3) = DVEC2(:,3)
  end subroutine copy_DVEC2_DVEC2a



  subroutine copy_DVEC2_DVEC2b()
    implicit none
    ! -------------------------------------------------------
    ! Copy DVEC2 into DVEC2a
    ! -------------------------------------------------------
    DVEC2b(:,1) = DVEC2(:,1)
    DVEC2b(:,2) = DVEC2(:,2)
    DVEC2b(:,3) = DVEC2(:,3)
  end subroutine copy_DVEC2_DVEC2b



  subroutine update_XMATRIX_DISPL_DVECI()
    implicit none
    integer :: K, IBC
    ! -------------------------------------------------------
    ! Update the column of XMATRIX corresponding to the displacements
    ! with the elements of DVECI
    ! -------------------------------------------------------
    DO K = 1, ndis_glob
        IBC = (K-1)*3
        XMATRIX(IBC+1,NUM_Ds_SAVED) = DVECI(K,1)
        XMATRIX(IBC+2,NUM_Ds_SAVED) = DVECI(K,2)
        XMATRIX(IBC+3,NUM_Ds_SAVED) = DVECI(K,3)
    ENDDO
  end subroutine update_XMATRIX_DISPL_DVECI



  subroutine update_DVEC_from_DVEC2(i_kode)
    implicit none
    integer :: K, IBC
    integer :: i_kode(:)
    ! -------------------------------------------------------
    ! Update the vector DVEC from DVEC2: only needed when the
    ! solver is in STATE 2 with dyndike and frictional elements
    ! -------------------------------------------------------
    DO K = 1, ndis_glob
        IBC = (K-1)*3
        IF (i_kode(K).LT.10) THEN
          DVEC(K,1) = DVEC2(K,1)
          DVEC(K,2) = DVEC2(K,2)
          DVEC(K,3) = DVEC2(K,3)
        ELSEIF (i_kode(K).EQ.10) THEN
          continue ! transfer nothing because set to 0
        ELSEIF (i_kode(K).EQ.11.OR.i_kode(K).EQ.12) THEN
          DVEC(K,1) = DVEC2(K,1)
          DVEC(K,2) = DVEC2(K,2)
        ELSEIF (i_kode(K).EQ.13) THEN
          DVEC(K,1) = DVEC2(K,1)
        ELSEIF (i_kode(K).EQ.14) THEN
          DVEC(K,2) = DVEC2(K,2)
        ELSEIF (i_kode(K).EQ.15) THEN
          DVEC(K,3) = DVEC2(K,3)
        ENDIF
    ENDDO
  end subroutine update_DVEC_from_DVEC2



  subroutine sum_DVECI_DVEC2()
    implicit none
    ! -------------------------------------------------------
    ! Sum the vectors DVECI and DVEC2 (for solver in iterative states)
    ! -------------------------------------------------------
    DVECI(:,1) = DVECI(:,1) + DVEC2(:,1)
    DVECI(:,2) = DVECI(:,2) + DVEC2(:,2)
    DVECI(:,3) = DVECI(:,3) + DVEC2(:,3)
  end subroutine sum_DVECI_DVEC2



  subroutine sum_DVEC2a_DVEC2b()
    implicit none
    ! -------------------------------------------------------
    ! Sum the vectors DVEC2a and DVEC2b in DEVEC2 (for solver in STATE 2)
    ! -------------------------------------------------------
    DVEC2(:,1) = DVEC2a(:,1) + DVEC2b(:,1)
    DVEC2(:,2) = DVEC2a(:,2) + DVEC2b(:,2)
    DVEC2(:,3) = DVEC2a(:,3) + DVEC2b(:,3)
  end subroutine sum_DVEC2a_DVEC2b



  subroutine print_DVECI()
    implicit none
    integer :: K
    ! -------------------------------------------------------
    ! print DVECI in the terminal
    ! -------------------------------------------------------
    WRITE(*,*) '       ID,     Ds (cm),     Dd (cm),    Dn (cm)'
    DO K = 1, ndis_glob
        WRITE(*,*) K, DVECI(K,1), DVECI(K,2), DVECI(K,3)
    ENDDO
  end subroutine print_DVECI



  subroutine write_xmatrix(fname)
    use global_inputs
    implicit none
    character(7) :: fname
    integer :: K, L
    ! -------------------------------------------------------
    ! Dump XMATRIX to ASCII file for inspection
    ! -------------------------------------------------------
    open(unit=20, file=fname, status='unknown')
    WRITE(20,*) ' --- XMATRIX ---'
    WRITE(20,*) ''
    WRITE(20,*) 'SHAPE OF XMATRIX', SHAPE(XMATRIX)
    WRITE(20,*) 'LBOUND:', LBOUND(XMATRIX,1), LBOUND(XMATRIX,2)
    WRITE(20,*) 'UBOUND:', UBOUND(XMATRIX,1), UBOUND(XMATRIX,2)
    WRITE(20,*) ''
    WRITE(20,*) 'INDEX OF COLUMN CONTAINNIG B.C.S:'
    WRITE(20,*) NUM_Ds_SAVED
    WRITE(20,*) ''
    WRITE(20,*) 'BOUNDARY CONDITION VECTOR:'
    WRITE(20,*) 'ROW-ID,   PLANE-ID,    B.C-ID,  KODE,   VECTOR'
    DO K = LBOUND(XMATRIX,1), UBOUND(XMATRIX,1)
      WRITE(20, *) K, XMATRIX(K,NUM_Ds_SAVED)
    ENDDO
    WRITE(20,*) ''
    WRITE(20,*) 'XMATRIX:'
    DO K = LBOUND(XMATRIX,1), UBOUND(XMATRIX,1)
      DO L = LBOUND(XMATRIX,2), UBOUND(XMATRIX,2)
        WRITE(20,*) K, L, XMATRIX(K,L)
      ENDDO
    ENDDO
    close(20)
  end subroutine write_xmatrix



  subroutine write_amatrix(fname)
    implicit none
    character(7) :: fname
    integer :: K, L
    ! -------------------------------------------------------
    ! Dump AMATRIX to ASCII file for inspection
    ! -------------------------------------------------------
    open(unit=20, file=fname, status='unknown')
    WRITE(20,*) ' --- AMATRIX ---'
    WRITE(20,*) ''
    WRITE(20,*) 'SHAPE OF AMATRIX', SHAPE(AMATRIX)
    WRITE(20,*) 'LBOUND:', LBOUND(AMATRIX,1), LBOUND(AMATRIX,2)
    WRITE(20,*) 'UBOUND:', UBOUND(AMATRIX,1), UBOUND(AMATRIX,2)
    WRITE(20,*) ''
    WRITE(20,*) 'AMATRIX:'
    DO K = LBOUND(AMATRIX,1), UBOUND(AMATRIX,1)
      DO L = LBOUND(AMATRIX,2), UBOUND(AMATRIX,2)
        WRITE(20,*) K, L, AMATRIX(K,L)
      ENDDO
    ENDDO
    close(20)
  end subroutine write_amatrix



  subroutine allocate_locked_arrays(ndis)
    implicit none
    integer, intent(in) :: ndis
    is_F_locked = .FALSE.
    is_D_locked = .FALSE.
    allocate(locked_KODE(ndis))
    allocate(locked_i_kode(ndis))
    allocate(locked_i_fcode(ndis))
    allocate(locked_element_fstatus(ndis))
    allocate(locked_element_dstatus(ndis))
    allocate(locked_BC(3,ndis))
  end subroutine allocate_locked_arrays



  subroutine deallocate_locked_arrays()
    implicit none
    deallocate(locked_KODE)
    deallocate(locked_i_kode)
    deallocate(locked_i_fcode)
    deallocate(locked_element_fstatus)
    deallocate(locked_element_dstatus)
    deallocate(locked_BC)
  end subroutine deallocate_locked_arrays



  subroutine lock_F(ndis, KODE, BC)
    ! -----------------------------------
    ! Lock fricitonal elements so that no displacement will be computed for them
    ! To be called in DO_ALL
    ! -----------------------------------
    use global_inputs
    implicit none
    integer :: L, ndis
    integer, intent(inout) :: KODE(1,1,ndis)
    real(4), intent(inout) :: BC(3,1,1,ndis)
    ! Update flag
    is_F_locked = .TRUE.
    ! Backup fields
    locked_KODE(:) = KODE(1,1,:)
    locked_i_kode  = i_kode
    locked_i_fcode = i_fcode
    locked_element_fstatus = element_fstatus
    locked_BC(1,:) = BC(1,1,1,:)
    locked_BC(2,:) = BC(2,1,1,:)
    locked_BC(3,:) = BC(3,1,1,:)
    ! change BCs to locked KODE 10
    do L=1, ndis
      if (i_fcode(L) .GT. 0) then
        is_locked(L) = .TRUE.
        i_fcode(L)  = 0         ! frictionless
        element_fstatus(L) = -2 ! frictionless
        KODE(1,1,L) = 10        ! BC in displacement
        i_kode(L)   = 10        ! BC in displacement
        BC(1,1,1,L) = 0         ! no displacement (locked)
        BC(2,1,1,L) = 0         ! no displacement (locked)
        BC(3,1,1,L) = 0         ! no displacement (locked)
      endif
    enddo
  end subroutine lock_F
    


  subroutine unlock_F(ndis, KODE, BC)
    ! To be called in DO_ALL
    use global_inputs
    implicit none
    integer :: L, ndis
    integer, intent(inout) :: KODE(1,1,ndis)
    real(4), intent(inout) :: BC(3,1,1,ndis)
    ! Update flag
    is_F_locked = .FALSE.
    ! Restore fields
    do L=1, ndis
      if (locked_i_fcode(L) .GT. 0) then
        is_locked(L) = .FALSE.
        i_fcode(L)  = locked_i_fcode(L)
        element_fstatus(L) = locked_element_fstatus(L)
        KODE(1,1,L) = locked_KODE(L)
        i_kode(L)   = locked_i_kode(L)
        BC(1,1,1,L) = locked_BC(1,L)
        BC(2,1,1,L) = locked_BC(2,L)
        BC(3,1,1,L) = locked_BC(3,L)
      endif
    enddo
  end subroutine unlock_F



  subroutine lock_D(ndis, KODE, BC)
    ! -----------------------------------
    ! Lock dyndike elements so that no displacement will be computed for them
    ! To be called in DO_ALL
    ! -----------------------------------
    use global_inputs
    implicit none
    integer :: L, ndis
    integer, intent(inout) :: KODE(1,1,ndis)
    real(4), intent(inout) :: BC(3,1,1,ndis)
    ! Update flag
    is_D_locked = .TRUE.
    ! Backup fields
    locked_KODE(:) = KODE(1,1,:)
    locked_i_kode  = i_kode
    locked_element_dstatus = element_dstatus
    locked_BC(1,:) = BC(1,1,1,:)
    locked_BC(2,:) = BC(2,1,1,:)
    locked_BC(3,:) = BC(3,1,1,:)
    ! change BCs to locked KODE 10
    do L=1, ndis
      if (element_dstatus(L) .GT. 0) then
        is_locked(L) = .TRUE.
        element_dstatus(L) = 0  ! not dyndike
        KODE(1,1,L) = 10        ! BC in displacement
        i_kode(L)   = 10        ! BC in displacement
        BC(1,1,1,L) = 0         ! no displacement (locked)
        BC(2,1,1,L) = 0         ! no displacement (locked)
        BC(3,1,1,L) = 0         ! no displacement (locked)
      endif
    enddo
  end subroutine lock_D
    

  subroutine unlock_D(ndis, KODE, BC)
    ! To be called in DO_ALL
    use global_inputs
    implicit none
    integer :: L, ndis
    integer, intent(inout) :: KODE(1,1,ndis)
    real(4), intent(inout) :: BC(3,1,1,ndis)
    ! Update flag
    is_D_locked = .FALSE.
    ! Restore fields
    do L=1, ndis
      if (locked_element_dstatus(L) .GT. 0) then
        is_locked(L) = .FALSE.
        element_dstatus(L)  = locked_element_dstatus(L)
        KODE(1,1,L) = locked_KODE(L)
        i_kode(L)   = locked_i_kode(L)
        BC(1,1,1,L) = locked_BC(1,L)
        BC(2,1,1,L) = locked_BC(2,L)
        BC(3,1,1,L) = locked_BC(3,L)
      endif
    enddo
  end subroutine unlock_D





end module global_arrays
