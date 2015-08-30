!to compile linking to netcdf library
!f90 -o myprogram myprogram.o -L/usr/local/netcdf/lib -lnetcdf
!status = nf90_open(path =
!"~/research/bogus/gcontrol.cam.0301-06.nc",ncid=ncid);
!      if(status /= nf90_NoErr) call handle_err(status)
!
!
!status= nf90_inq_varid(ncid,"Z3",Z3VarID) 
!      if(status /= nf90_NoErr) call handle_err(status)
!
!
!status = nf90_inq_dimids(ncid, "Z3", Z3DimIDs)
!      if (status /= nf90_noerr) call handle_err(status)
!
!
!status = nf90_inquire_variable(ncid, Z3VarID, ndims = numDims, natts =
!numAtts)
!
!status = nf90_inquire_variable(ncid, Z3VarID, ndims = numDims, natts =
!numAtts)
!      if(status /= nf90_NoErr) call handle_err(status)
!status = nf90_inquire_variable(ncid, Z3VarID, dimIDs =
!Z3DimIDs(:numDims))
!
!
!status = nf90_inquire_dimension(ncid, dimIDs(1), len = numLons)
!      if(status /= nf90_NoErr) call handle_err(status)
!status = nf90_inquire_dimension(ncid, dimIDs(2), len = numLats)
!      if(status /= nf90_NoErr) call handle_err(status)
!status = nf90_inquire_dimension(ncid, dimIDs(3), len = numTimes)
!      if(status /= nf90_NoErr) call handle_err(status)

