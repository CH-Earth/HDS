module ascii_util
USE type_HDS
implicit none
private
public::read_csv
contains

 ! *********************************************************************************************************
 ! public subroutine read_csv: read csv file
 ! *********************************************************************************************************
 subroutine read_csv(infile,ixVec,varnames,xData,nLines,ierr,message)
 implicit none
 ! input/output
 character(*),intent(in)                    :: infile            ! filename
 integer(i4b),intent(in)                    :: ixVec(:)          ! indices in data structure
 character(*),intent(in)                    :: varnames(:)       ! desired variables
 type(hStruct), intent(inout), allocatable  :: xData(:)          ! data structure
 integer(i4b),intent(out)                   :: nLines            ! number of lines in the file
 integer(i4b),intent(out)                   :: ierr              ! error code
 character(*),intent(out)                   :: message           ! error message
 ! local variables
 integer(i4b)                               :: unt               ! file unit
 integer(i4b)                               :: iVar              ! variable index
 integer(i4b)                               :: nVar              ! number of variables
 integer(i4b)                               :: iLine             ! line index
 integer(i4b)                               :: nWords            ! number of variables in the header
 integer(i4b)                               :: length            ! length og a given variable name
 character(len=256)                         :: cHeader           ! header in the .csv file
 character(len=256), allocatable            :: headVarnames(:)   ! variable names in the header
 logical(lgt)      , allocatable            :: flagVarnames(:)   ! logical vector of matches for each variable
 integer(i4b)      , allocatable            :: ixMatch(:)        ! matches in the header (0,1)
 real(rkind)       , allocatable            :: lineData(:)       ! data for a given line
 integer(i4b)                               :: jxLoc(1)          ! indices of desired variables in the output
 integer(i4b), dimension(size(ixVec))       :: jxVec             ! indices of desired variables in the output
 character(len=256)                         :: cMessage          ! error message for downwind routine

 ! initialize errors
 ierr=0; message="read_csv/"

 ! get the number of variables
 nVar = size(ixVec)
 if(size(varnames) /= nVar)then; ierr=20; message=trim(message)//'mismatch length in indices and names'; return; endif

 ! get the number of lines in the file
 ! NOTE: Also get the header
 call get_lines(infile, header=cHeader, nLines=nLines, ierr=ierr, message=cMessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! convert the header to a vector of variable names
 nWords = count(transfer(cHeader, ' ', len_trim(cHeader)) == ',') + 1               ! number of words in the header (number of commas plus one)
 allocate(headVarnames(nWords), flagVarnames(nWords), lineData(nWords), stat=ierr)  ! define string array equal to the number of words
 if(ierr/=0)then; message=trim(message)//'problem allocating space'; return; endif  ! check allocation
 read(cHeader,*) (headVarnames(iVar), iVar=1,nWords)                                ! populate the string array

 ! get the index of the desired variables in the csv file
 ! NOTE: the "transfer" intrinsic is needed because maxloc does not support logical data types
 do iVar=1,nVar
  length       = len_trim(varnames(iVar))
  flagVarnames = varnames(iVar)(1:length) == headVarnames(:)(1:length)
  jxLoc        = maxloc(transfer(flagVarnames, ixMatch))
  jxVec(iVar)  = jxLoc(1)
 end do

 ! allocate space for the variables in the data structure
 allocate(xData(nVar), stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem allocating space for the variable dimension in the data structure'; return; endif
 do iVar=1,nVar
  allocate(xData(iVar)%dat(nLines), stat=ierr)
  if(ierr/=0)then; message=trim(message)//'problem allocating space for the data dimension in the data structure'; return; endif
 end do

 ! open file and read the header
 call file_open(trim(infile),unt,ierr,cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; end if
 read(unt,'(a)') cHeader

 ! read data
 do iLine=1,nLines
  read(unt,*) (lineData(iVar), iVar=1,nWords)
  forall(iVar=1:nVar) xData(iVar)%dat(iLine) = lineData(jxVec(iVar))
 end do

 ! close csv file
 close(unit=unt,iostat=ierr)
 if(ierr/=0)then;message=trim(message)//'problem closing file '//trim(infile); return; end if

 ! deallocate temporary variables
 deallocate(headVarnames, flagVarnames, lineData, stat=ierr)
 if(ierr/=0)then; message=trim(message)//'problem deallocating space'; return; endif

 end subroutine read_csv

 ! *********************************************************************************************************
 ! public subroutine file_open: open file
 ! *********************************************************************************************************
 subroutine file_open(infile,unt,ierr,message)
 implicit none
 ! declare dummy variables
 character(*),intent(in)              :: infile      ! filename
 integer(i4b),intent(out)             :: unt         ! file unit
 integer(i4b),intent(out)             :: ierr        ! error code
 character(*),intent(out)             :: message     ! error message
 ! declare local variables
 logical(lgt)                         :: xist        ! .TRUE. if the file exists
 logical(lgt)                         :: xopn        ! .TRUE. if the file is already open
 ! initialize errors
 ierr=0; message="file_open/"
 ! check if the file exists
 inquire(file=trim(infile),exist=xist) ! Check for existence of file
 if(.not.xist)then
   message=trim(message)//"FileNotFound[file='"//trim(infile)//"']"
   ierr=10; return
 end if
 ! check if the file is already open
 inquire(file=trim(infile),opened=xopn) ! Check if the file is open
 if(xopn)then
  message=trim(message)//"FileAlreadyOpen['"//trim(infile)//"']"
  ierr=20; return
 end if
 ! open file
 open(newunit=unt,file=trim(infile),status="old",action="read",iostat=ierr)
 if(ierr/=0)then
   message=trim(message)//"OpenError['"//trim(infile)//"']"
   ierr=20; return
 end if
 end subroutine file_open

 ! *********************************************************************************************************
 ! public subroutine get_lines: get the number of lines in a file
 ! *********************************************************************************************************
 subroutine get_lines(infile,header,nLines,ierr,message)
 implicit none
 ! declare dummy variables
 character(*),intent(in)              :: infile         ! filename
 character(*),intent(out),optional    :: header         ! file header
 integer(i4b),intent(out)             :: nLines         ! number of lines in the file
 integer(i4b),intent(out)             :: ierr           ! error code
 character(*),intent(out)             :: message        ! error message
! declare local variables
 integer(i4b)                         :: unt            ! file unit
 integer(i4b)                         :: iline          ! loop through lines in the file
 integer(i4b),parameter               :: maxLines=10000 ! maximum number of valid lines in a file
 character(len=2048)                  :: temp           ! character data or a given line
 integer(i4b)                         :: iend           ! index to indicate end of the file
 character(len=256)                   :: cmessage       ! error message for downwind routine
 ! initialize errors
 ierr=0; message='get_lines/'

 ! * open file
 call file_open(trim(infile),unt,ierr,cmessage)
 if(ierr/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! * read header
 if(present(header)) read(unt,'(a)') header

 ! * count non-comment lines in the file
 nLines=0  ! initialize the counter for the valid lines
 do iline=1,maxLines
  read(unt,'(a)',iostat=iend)temp; if(iend/=0)exit    ! read line of data
  if (temp(1:1)=='!' .or. temp == '')cycle            ! skip comment and empty lines
  nLines = nLines+1
  if (iline==maxLines)then; ierr=20; message=trim(message)//"exceedMaxLines"; return; end if
 end do  ! looping through the lines in the file (exit clause above will kick in)

 ! * close ascii file
 close(unit=unt,iostat=ierr)
 if(ierr/=0)then;message=trim(message)//'problem closing file '//trim(infile); return; end if

 end subroutine get_lines

end module ascii_util
