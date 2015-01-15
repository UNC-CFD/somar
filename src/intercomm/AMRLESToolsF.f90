module AMRLESToolsF
use mpi

implicit none
private

public AMRLESToolsSetup, boundsCheck, recvDataFromAMR

save
integer :: amr2dns_leader, inter_comm
logical :: isSetup = .false.

! Which indices refer to x, y, z, and comp?
integer, parameter :: xi = 1
integer, parameter :: yi = 2
integer, parameter :: zi = 3
integer, parameter :: compi = 4

contains


! ------------------------------------------------------------------------------
! Gets the communication info from the user.
! ------------------------------------------------------------------------------
subroutine AMRLESToolsSetup (a_amr2dns_leader, a_inter_comm)
    integer, intent(in) :: a_amr2dns_leader, a_inter_comm

    amr2dns_leader = a_amr2dns_leader
    inter_comm = a_inter_comm
    isSetup = .true.
end subroutine AMRLESToolsSetup


! ------------------------------------------------------------------------------
! Checks the bounds of an array
! ------------------------------------------------------------------------------
integer function boundsCheck (dataFAB, iBoxLo, iBoxHi, iCompLo, iCompHi)
    ! Subroutine input/output parameters.
    double precision, dimension(:,:,:,:), intent(in) :: dataFAB
    integer, dimension(xi:zi), intent(in)            :: iBoxLo, iBoxHi
    integer, intent(in)                              :: iCompLo, iCompHi

    integer :: ierr

    ierr = 0
    if (size(dataFAB,xi) .lt. (iBoxHi(xi)-iBoxLo(xi)+1)) ierr = 1
    if (size(dataFAB,yi) .lt. (iBoxHi(yi)-iBoxLo(yi)+1)) ierr = 2
    if (size(dataFAB,zi) .lt. (iBoxHi(zi)-iBoxLo(zi)+1)) ierr = 4
    if (size(dataFAB,compi) .lt. (iCompHi-iCompLo+1))    ierr = 8

    boundsCheck = ierr
end function boundsCheck


! ------------------------------------------------------------------------------
! Gets a FAB from the AMR.
! ------------------------------------------------------------------------------
subroutine recvDataFromAMR (dataFAB, iBoxLo, iBoxHi, nComp, boxType, doAllocOpt)
    ! Subroutine input/output parameters.
    double precision, dimension(:,:,:,:), allocatable, intent(inout) :: dataFAB
    integer, dimension(xi:zi), intent(out)      :: iBoxLo, iBoxHi, boxType
    integer, intent(out)                        :: nComp
    logical, optional                           :: doAllocOpt

    ! Transfer variables
    integer, parameter                          :: metaDataCount = 10
    integer, dimension(metaDataCount)           :: metaData
    integer                                     :: sizeOfTransfer
    double precision, dimension(:), allocatable :: linearData

    ! Additional variables
    integer :: iCompLo, iCompHi, tmp, ierr
    logical :: doAlloc, needsAlloc


    ! Is this module set up for use?
    if (.not. isSetup) then
        print*, 'Please call setup before any other AMRLESToolsF subroutines.'
        call MAYDAYERROR()
    endif

    ! Recv metaData
    CALL MPI_BCAST(metaData, metaDataCount, MPI_INTEGER, &
                   amr2dns_leader, inter_comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
        print*, 'recvDataFromAMR: MPI_BCAST error when asking for metadata.'
        call MAYDAYERROR()
    endif

    ! Sort metaData
    iBoxLo(xi) = metaData(1)
    iBoxLo(yi) = metaData(2)
    iBoxLo(zi) = metaData(3)

    iBoxHi(xi) = metaData(4)
    iBoxHi(yi) = metaData(5)
    iBoxHi(zi) = metaData(6)

    nComp = metaData(7)

    boxType(xi) = metaData(8)
    boxType(yi) = metaData(9)
    boxType(zi) = metaData(10)

    ! Process metaData
    if (boxType(zi) .eq. -1) then
        ! Original data was 2D. In this case, we want to put the data into the
        ! xi and zi planes and leave yi as the streamwise component, to be
        ! handled by the user.
        tmp=iBoxLo(zi); iBoxLo(zi)=iBoxLo(yi); iBoxLo(yi)=tmp
        tmp=iBoxHi(zi); iBoxHi(zi)=iBoxHi(yi); iBoxHi(yi)=tmp
        tmp=boxType(zi); boxType(zi)=boxType(yi); boxType(yi)=tmp
    endif

    iCompLo = 1
    iCompHi = iCompLo+nComp-1


    ! By default, we will allocate space.
    if (.not. present(doAllocOpt)) then
        doAlloc = .true.
    else
        doAlloc = doAllocOpt
    endif

    ! Do we need to allocate space?
    needsAlloc = .true.
    if (allocated(dataFAB)) then
        ! Check bounds
        ierr = boundsCheck(dataFAB, iBoxLo, iBoxHi, iCompLo, iCompHi)
        if (ierr .eq. 0) then
            needsAlloc = .false.
        else
            deallocate(dataFAB)
        endif
    endif

    ! Allocate space, if needed.
    if (needsAlloc) then
        if (doAlloc) then
            allocate(dataFAB(iBoxLo(xi):iBoxHi(xi), &
                             iBoxLo(yi):iBoxHi(yi), &
                             iBoxLo(zi):iBoxHi(zi), &
                             iCompLo:iCompHi))
        else
            print*, 'recvDataFromAMR: dataFAB is not allocated and doAlloc is false.'
            call MAYDAYERROR()
        endif
    endif

    ! Double check bounds
    ierr = boundsCheck(dataFAB, iBoxLo, iBoxHi, iCompLo, iCompHi)
    if (ierr .ne. 0) then
        print*, 'recvDataFromAMR: dataFAB was not allocated properly.'
        call MAYDAYERROR()
    endif

    ! Recv data!
    sizeOfTransfer = size(dataFAB)
    ALLOCATE(linearData(sizeOfTransfer))

    CALL MPI_BCAST(linearData, sizeOfTransfer, MPI_DOUBLE_PRECISION, &
                   amr2dns_leader, inter_comm, ierr)

    if (ierr .ne. MPI_SUCCESS) then
        print*, 'recvDataFromAMR: MPI_BCAST error when asking for dataFAB.'
        call MAYDAYERROR()
    endif

    ! Change shape of data
    dataFAB(:,:,:,:) = reshape(linearData(:), shape(dataFAB))

end subroutine recvDataFromAMR

end module AMRLESToolsF
