  module myncread
  use netcdf
  implicit none

  contains
  subroutine readnc_int(varname,var,ncID)
  implicit none
  integer,intent(out):: var !output an integer like nx
  integer            :: varID, i, status
  integer,intent(in) :: ncID
  character(len=*),intent(in):: varname
    status= nf90_inq_varID(ncID, varname, varID)
    status= nf90_get_var(ncID, varID, var)
  return
  end subroutine readnc_int
 
  subroutine readnc_real(varname,var,ncID)
  implicit none
  real,intent(out):: var !output an integer like nx
  integer            :: varID, i, status
  integer,intent(in) :: ncID
  character(len=*),intent(in):: varname
    status= nf90_inq_varID(ncID, varname, varID)
    status= nf90_get_var(ncID, varID, var)
  return
  end subroutine readnc_real
 
  subroutine readnc_int1d(varname,var,ncID)
  implicit none
  integer,allocatable,intent(out),dimension(:):: var !output an integer like nx
  integer            :: varID, i, dimIDs(1), status
  integer,intent(in) :: ncID
  integer:: vecLen
  character(len=*),intent(in):: varname
    status= nf90_inq_varID(ncID, varname, varID)
    status= nf90_inquire_variable(ncID, varID, dimids = dimIDs)
    status= nf90_inquire_dimension(ncID, dimIDs(1), len = vecLen)
    allocate(var(vecLen))
    status= nf90_get_var(ncID, varID, var)
  return
  end subroutine readnc_int1d

  subroutine readnc_real1d(varname,var,ncID)
  implicit none
  real,allocatable,intent(out),dimension(:):: var !output an integer like nx
  integer            :: varID, i, dimIDs(1), status
  integer,intent(in) :: ncID
  integer:: vecLen
  character(len=*),intent(in):: varname
    status= nf90_inq_varID(ncID, varname, varID)
    status= nf90_inquire_variable(ncID, varID, dimids = dimIDs)
    status= nf90_inquire_dimension(ncID, dimIDs(1), len = vecLen)
    allocate(var(vecLen))
    status= nf90_get_var(ncID, varID, var)
  return
  end subroutine readnc_real1d
  
 
  subroutine readnc_real2d(varname,var,ncID,nx,ny)
  implicit none
  integer,intent(in):: nx, ny
!  real,allocatable,dimension(:,:,:):: tempvar
  real,allocatable,intent(out),dimension(:,:):: var !output an integer like nx
  integer            :: varID, i,  status
  integer,intent(in) :: ncID
  character(len=*),intent(in):: varname
    status= nf90_inq_varID(ncID, varname, varID)
    !status= NF90_INQ_DIMID(ncID,'lon',dimid)
    !status= NF90_inquire_dimension(ncID,dimid,len=nx)
    !status= NF90_INQ_DIMID(ncID,'lat',dimid)
    !status= NF90_inquire_dimension(ncID,dimid,len=ny)
!    allocate(tempvar(1,leny,lenx))
      allocate(var(nx,ny))
    status= nf90_get_var(ncID, varID, var)
  return
  end subroutine readnc_real2d


  subroutine readnc_real3d(varname,var,ncID,nx,ny,nz)
  implicit none
  integer,intent(in):: nx, ny, nz
  real,allocatable,intent(out),dimension(:,:,:):: var !output an integer like nx
  integer            :: varID, i, status
  integer,intent(in) :: ncID
  integer:: vecLen(3)
  character(len=*),intent(in):: varname
    status= nf90_inq_varID(ncID, varname, varID)
    !status= NF90_INQ_DIMID(ncID,'lon',dimid)
    !status= NF90_inquire_dimension(ncID,dimid,len=lenx)
    !status= NF90_INQ_DIMID(ncID,'lat',dimid)
    !status= NF90_inquire_dimension(ncID,dimid,len=leny)
    !status= NF90_INQ_DIMID(ncID,'lev',dimid)
    !status= NF90_inquire_dimension(ncID,dimid,len=lenz)
    allocate(var(nx,ny,nz)) 
    status= nf90_get_var(ncID, varID, var)
!    var = tempvar(1,:,:,:)
  return
  end subroutine readnc_real3d
 
 
  subroutine handle_err(status)
  implicit none
  integer, intent ( in) :: status
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  return
  end subroutine handle_err
 
  end module myncread
