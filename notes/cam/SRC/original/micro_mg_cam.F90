module micro_mg_cam

!---------------------------------------------------------------------------------
!
!  CAM Interfaces for MG microphysics
!
!---------------------------------------------------------------------------------

use shr_kind_mod,   only: r8=>shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, pver, pverp
use physconst,      only: gravit, rair, tmelt, cpair, rh2o, rhoh2o, &
     latvap, latice, mwdry
use phys_control,   only: phys_getopts


use physics_types,  only: physics_state, physics_ptend, physics_ptend_init, &
     physics_state_copy, physics_ptend_copy, &
     physics_update, physics_state_dealloc, &
     physics_ptend_sum
use physics_buffer, only: physics_buffer_desc, pbuf_add_field, dyn_time_lvls, &
     pbuf_old_tim_idx, pbuf_get_index, dtype_r8, dtype_i4, &
     pbuf_get_field, pbuf_set_field
use constituents,   only: cnst_add, cnst_get_ind, &
     cnst_name, cnst_longname, sflxnam, apcnst, bpcnst, pcnst

use cldwat2m_macro, only: rhmini

use cam_history,    only: addfld, add_default, phys_decomp, outfld

use cam_logfile,    only: iulog
use abortutils,     only: endrun
use error_messages, only: handle_errmsg
use ref_pres,       only: top_lev=>trop_cloud_top_lev


implicit none
private
save

public :: &
     micro_mg_cam_readnl,          &
     micro_mg_cam_register,        &
     micro_mg_cam_init_cnst,       &
     micro_mg_cam_implements_cnst, &
     micro_mg_cam_init,            &
     micro_mg_cam_tend

integer :: micro_mg_version     = 1      ! Version number for MG.
integer :: micro_mg_sub_version = 0      ! Second part of version number.

logical :: microp_uniform = .false.

logical, public :: do_cldliq ! Prognose cldliq flag
logical, public :: do_cldice ! Prognose cldice flag

integer, parameter :: ncnst = 4       ! Number of constituents
character(len=8), parameter :: &      ! Constituent names
     cnst_names(ncnst) = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE'/)

integer :: &
     ixcldliq,      &! cloud liquid amount index
     ixcldice,      &! cloud ice amount index
     ixnumliq,      &! cloud liquid number index
     ixnumice        ! cloud ice water index

! Physics buffer indices for fields registered by this module
integer :: &
     cldo_idx,           &
     qme_idx,            &
     prain_idx,          &
     nevapr_idx,         &
     wsedl_idx,          &
     rei_idx,            &
     rel_idx,            &
     dei_idx,            &
     mu_idx,             &
     lambdac_idx,        &
     iciwpst_idx,        &
     iclwpst_idx,        &
     des_idx,            &
     icswp_idx,          &
     cldfsnow_idx,       &
     rate1_cw2pr_st_idx = -1, &
     ls_flxprc_idx,      &
     ls_flxsnw_idx,      &
     relvar_idx,         &
     cmeliq_idx,         &
     accre_enhan_idx

! Fields needed as inputs to COSP
integer :: &
     ls_mrprc_idx,    ls_mrsnw_idx,    &
     ls_reffrain_idx, ls_reffsnow_idx, &
     cv_reffliq_idx,  cv_reffice_idx

! Fields needed by Park macrophysics
integer :: &
     cc_t_idx,  cc_qv_idx, &
     cc_ql_idx, cc_qi_idx, &
     cc_nl_idx, cc_ni_idx, &
     cc_qlst_idx

! Used to replace aspects of MG microphysics
! (e.g. by CARMA)
integer :: tnd_qsnow_idx, tnd_nsnow_idx, re_ice_idx

! Index fields for precipitation efficiency.
integer :: acpr_idx, acgcme_idx, acnum_idx

! Physics buffer indices for fields registered by other modules
integer :: &
     ast_idx = -1,            &
     aist_idx = -1,           &
     alst_idx = -1,           &
     cld_idx = -1,            &
     concld_idx = -1

integer :: &
     naai_idx = -1,           &
     naai_hom_idx = -1,       &
     npccn_idx = -1,          &
     rndst_idx = -1,          &
     nacon_idx = -1,          &
     prec_str_idx = -1,       &
     snow_str_idx = -1,       &
     prec_pcw_idx = -1,       &
     snow_pcw_idx = -1,       &
     prec_sed_idx = -1,       &
     snow_sed_idx = -1

!===============================================================================
contains
!===============================================================================

subroutine micro_mg_cam_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Namelist variables
  logical :: micro_mg_do_cldice   = .true. ! do_cldice = .true., MG microphysics is prognosing cldice
  logical :: micro_mg_do_cldliq   = .true. ! do_cldliq = .true., MG microphysics is prognosing cldliq

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'micro_mg_cam_readnl'

  namelist /micro_mg_nl/ micro_mg_version, micro_mg_sub_version, &
       micro_mg_do_cldice, micro_mg_do_cldliq

  !-----------------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'micro_mg_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, micro_mg_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)

     ! set local variables
     do_cldice  = micro_mg_do_cldice
     do_cldliq  = micro_mg_do_cldliq

     ! Verify that version numbers are valid.
     select case (micro_mg_version)
     case (1)
        select case (micro_mg_sub_version)
        case(0)
           ! MG version 1.0
        case(5)
           ! MG version 1.5 - MG2 development
        case default
           call bad_version_endrun()
        end select
     case default
        call bad_version_endrun()
     end select

  end if

#ifdef SPMD
  ! Broadcast namelist variables
  call mpibcast(micro_mg_version,     1, mpiint, 0, mpicom)
  call mpibcast(micro_mg_sub_version, 1, mpiint, 0, mpicom)
  call mpibcast(do_cldice,            1, mpilog, 0, mpicom)
  call mpibcast(do_cldliq,            1, mpilog, 0, mpicom)
#endif

contains

  subroutine bad_version_endrun
    ! Endrun wrapper with a more useful error message.
    character(len=128) :: errstring
    write(errstring,*) "Invalid version number specified for MG microphysics: ", &
         micro_mg_version,".",micro_mg_sub_version
    call endrun(errstring)
  end subroutine bad_version_endrun

end subroutine micro_mg_cam_readnl

!================================================================================================

subroutine micro_mg_cam_register

  ! Register microphysics constituents and fields in the physics buffer.
  !-----------------------------------------------------------------------

  logical :: prog_modal_aero

  call phys_getopts( prog_modal_aero_out=prog_modal_aero )

  ! Register microphysics constituents and save indices.

  call cnst_add(cnst_names(1), mwdry, cpair, 0._r8, ixcldliq, &
       longname='Grid box averaged cloud liquid amount', is_convtran1=.true.)
  call cnst_add(cnst_names(2), mwdry, cpair, 0._r8, ixcldice, &
       longname='Grid box averaged cloud ice amount', is_convtran1=.true.)
  ! The next statements should have "is_convtran1=.true.", but this would change
  ! answers.
  call cnst_add(cnst_names(3), mwdry, cpair, 0._r8, ixnumliq, &
       longname='Grid box averaged cloud liquid number', is_convtran1=.false.)
  call cnst_add(cnst_names(4), mwdry, cpair, 0._r8, ixnumice, &
       longname='Grid box averaged cloud ice number', is_convtran1=.false.)

  ! Request physics buffer space for fields that persist across timesteps.

  call pbuf_add_field('CLDO','global',dtype_r8,(/pcols,pver,dyn_time_lvls/), cldo_idx)

  ! Physics buffer variables for convective cloud properties.

  call pbuf_add_field('QME',        'physpkg',dtype_r8,(/pcols,pver/), qme_idx)
  call pbuf_add_field('PRAIN',      'physpkg',dtype_r8,(/pcols,pver/), prain_idx)
  call pbuf_add_field('NEVAPR',     'physpkg',dtype_r8,(/pcols,pver/), nevapr_idx)

  call pbuf_add_field('WSEDL',      'physpkg',dtype_r8,(/pcols,pver/), wsedl_idx)

  call pbuf_add_field('REI',        'physpkg',dtype_r8,(/pcols,pver/), rei_idx)
  call pbuf_add_field('REL',        'physpkg',dtype_r8,(/pcols,pver/), rel_idx)

  ! Mitchell ice effective diameter for radiation
  call pbuf_add_field('DEI',        'physpkg',dtype_r8,(/pcols,pver/), dei_idx)
  ! Size distribution shape parameter for radiation
  call pbuf_add_field('MU',         'physpkg',dtype_r8,(/pcols,pver/), mu_idx)
  ! Size distribution shape parameter for radiation
  call pbuf_add_field('LAMBDAC',    'physpkg',dtype_r8,(/pcols,pver/), lambdac_idx)

  ! Stratiform only in cloud ice water path for radiation
  call pbuf_add_field('ICIWPST',    'physpkg',dtype_r8,(/pcols,pver/), iciwpst_idx)
  ! Stratiform in cloud liquid water path for radiation
  call pbuf_add_field('ICLWPST',    'physpkg',dtype_r8,(/pcols,pver/), iclwpst_idx)

  call pbuf_add_field('DES',        'physpkg',dtype_r8,(/pcols,pver/), des_idx)          ! Snow effective diameter for radiation
  call pbuf_add_field('ICSWP',      'physpkg',dtype_r8,(/pcols,pver/), icswp_idx)        ! In cloud snow water path for radiation
  ! Cloud fraction for liquid drops + snow
  call pbuf_add_field('CLDFSNOW ',  'physpkg',dtype_r8,(/pcols,pver,dyn_time_lvls/), cldfsnow_idx)

  if (prog_modal_aero) then
     call pbuf_add_field('RATE1_CW2PR_ST','physpkg',dtype_r8,(/pcols,pver/), rate1_cw2pr_st_idx)   ! rce 2010/05/01
  endif

  call pbuf_add_field('LS_FLXPRC',  'physpkg',dtype_r8,(/pcols,pverp/), ls_flxprc_idx)
  call pbuf_add_field('LS_FLXSNW',  'physpkg',dtype_r8,(/pcols,pverp/), ls_flxsnw_idx)


  ! Fields needed as inputs to COSP
  call pbuf_add_field('LS_MRPRC',   'physpkg',dtype_r8,(/pcols,pver/), ls_mrprc_idx)
  call pbuf_add_field('LS_MRSNW',   'physpkg',dtype_r8,(/pcols,pver/), ls_mrsnw_idx)
  call pbuf_add_field('LS_REFFRAIN','physpkg',dtype_r8,(/pcols,pver/), ls_reffrain_idx)
  call pbuf_add_field('LS_REFFSNOW','physpkg',dtype_r8,(/pcols,pver/), ls_reffsnow_idx)
  call pbuf_add_field('CV_REFFLIQ', 'physpkg',dtype_r8,(/pcols,pver/), cv_reffliq_idx)
  call pbuf_add_field('CV_REFFICE', 'physpkg',dtype_r8,(/pcols,pver/), cv_reffice_idx)

  ! CC_* Fields needed by Park macrophysics
  call pbuf_add_field('CC_T',     'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_t_idx)
  call pbuf_add_field('CC_qv',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_qv_idx)
  call pbuf_add_field('CC_ql',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_ql_idx)
  call pbuf_add_field('CC_qi',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_qi_idx)
  call pbuf_add_field('CC_nl',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_nl_idx)
  call pbuf_add_field('CC_ni',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_ni_idx)
  call pbuf_add_field('CC_qlst',  'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_qlst_idx)

  ! Additional pbuf for CARMA interface
  call pbuf_add_field('TND_QSNOW',  'physpkg',dtype_r8,(/pcols,pver/), tnd_qsnow_idx)
  call pbuf_add_field('TND_NSNOW',  'physpkg',dtype_r8,(/pcols,pver/), tnd_nsnow_idx)
  call pbuf_add_field('RE_ICE',     'physpkg',dtype_r8,(/pcols,pver/), re_ice_idx)

  ! Precipitation efficiency fields across timesteps.
  call pbuf_add_field('ACPRECL',    'global',dtype_r8,(/pcols/), acpr_idx)   ! accumulated precip
  call pbuf_add_field('ACGCME',     'global',dtype_r8,(/pcols/), acgcme_idx) ! accumulated condensation
  call pbuf_add_field('ACNUM',      'global',dtype_i4,(/pcols/), acnum_idx)  ! counter for accumulated # timesteps

  ! SGS variability
  call pbuf_add_field('RELVAR',     'global',dtype_r8,(/pcols,pver/), relvar_idx)
  call pbuf_add_field('ACCRE_ENHAN','global',dtype_r8,(/pcols,pver/), accre_enhan_idx)

end subroutine micro_mg_cam_register

!===============================================================================

function micro_mg_cam_implements_cnst(name)

  ! Return true if specified constituent is implemented by the
  ! microphysics package

  character(len=*), intent(in) :: name        ! constituent name
  logical :: micro_mg_cam_implements_cnst    ! return value

  ! Local workspace
  integer :: m
  !-----------------------------------------------------------------------

  micro_mg_cam_implements_cnst = any(name == cnst_names)

end function micro_mg_cam_implements_cnst

!===============================================================================

subroutine micro_mg_cam_init_cnst(name, q, gcid)

  ! Initialize the microphysics constituents, if they are
  ! not read from the initial file.

  character(len=*), intent(in)  :: name     ! constituent name
  real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
  integer,          intent(in)  :: gcid(:)  ! global column id
  !-----------------------------------------------------------------------

  if (micro_mg_cam_implements_cnst(name)) q = 0.0_r8

end subroutine micro_mg_cam_init_cnst

!===============================================================================

subroutine micro_mg_cam_init(pbuf2d)
  use time_manager, only: is_first_step
  use micro_mg1_0, only: micro_mg_init1_0 => micro_mg_init
  use micro_mg1_5, only: micro_mg_init1_5 => micro_mg_init

  !-----------------------------------------------------------------------
  !
  ! Initialization for MG microphysics
  !
  !-----------------------------------------------------------------------

  type(physics_buffer_desc), pointer :: pbuf2d(:,:)

  integer :: m, mm
  logical :: history_amwg         ! output the variables used by the AMWG diag package
  logical :: history_budget       ! Output tendencies and state variables for CAM4
  ! temperature, water vapor, cloud ice and cloud
  ! liquid budgets.
  integer :: budget_histfile      ! output history file number for budget fields

  character(128) :: errstring     ! return status (non-blank for error return)

  !-----------------------------------------------------------------------

  if (masterproc) then
     write(iulog,"(A,I2,A,I2)") "Initializing MG version ",micro_mg_version,".",micro_mg_sub_version
     if (.not. do_cldliq) &
          write(iulog,*) "MG prognostic cloud liquid has been turned off via namelist."
     if (.not. do_cldice) &
          write(iulog,*) "MG prognostic cloud ice has been turned off via namelist."
  end if

  select case (micro_mg_version)
  case (1)
     select case (micro_mg_sub_version)
     case (0)
        call micro_mg_init1_0( &
             r8, gravit, rair, rh2o, cpair, &
             rhoh2o, tmelt, latvap, latice, &
             rhmini, errstring)
     case (5)
        call micro_mg_init1_5( &
             r8, gravit, rair, rh2o, cpair, &
             tmelt, latvap, latice, rhmini, &
             microp_uniform, do_cldice, errstring)
     end select
  end select

  call handle_errmsg(errstring, subname="micro_mg_init")

  ! Register history variables
  do m = 1, ncnst
     call cnst_get_ind(cnst_names(m), mm)
     if ( any(mm == (/ ixcldliq, ixcldice /)) ) then
        ! mass mixing ratios
        call addfld(cnst_name(mm), 'kg/kg   ', pver, 'A', cnst_longname(mm)                   , phys_decomp)
        call addfld(sflxnam(mm),   'kg/m2/s ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp)
     else if ( any(mm == (/ ixnumliq, ixnumice /)) ) then
        ! number concentrations
        call addfld(cnst_name(mm), '1/kg    ', pver, 'A', cnst_longname(mm)                   , phys_decomp)
        call addfld(sflxnam(mm),   '1/m2/s  ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp)
     else
        call endrun( "micro_mg_cam_init: &
             &Could not call addfld for constituent with unknown units.")
     endif
  end do

  call addfld(apcnst(ixcldliq), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldliq))//' after physics'  , phys_decomp)
  call addfld(apcnst(ixcldice), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldice))//' after physics'  , phys_decomp)
  call addfld(bpcnst(ixcldliq), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldliq))//' before physics' , phys_decomp)
  call addfld(bpcnst(ixcldice), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldice))//' before physics' , phys_decomp)

  call addfld ('CME      ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap within the cloud'                      ,phys_decomp)
  call addfld ('PRODPREC ', 'kg/kg/s ', pver, 'A', 'Rate of conversion of condensate to precip'              ,phys_decomp)
  call addfld ('EVAPPREC ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling precip'                   ,phys_decomp)
  call addfld ('EVAPSNOW ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling snow'                     ,phys_decomp)
  call addfld ('HPROGCLD ', 'W/kg'    , pver, 'A', 'Heating from prognostic clouds'                          ,phys_decomp)
  call addfld ('FICE     ', 'fraction', pver, 'A', 'Fractional ice content within cloud'                     ,phys_decomp)
  call addfld ('ICWMRST  ', 'kg/kg   ', pver, 'A', 'Prognostic in-stratus water mixing ratio'                ,phys_decomp)
  call addfld ('ICIMRST  ', 'kg/kg   ', pver, 'A', 'Prognostic in-stratus ice mixing ratio'                  ,phys_decomp)

  ! MG microphysics diagnostics
  call addfld ('QCSEVAP  ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling cloud water'              ,phys_decomp)
  call addfld ('QISEVAP  ', 'kg/kg/s ', pver, 'A', 'Rate of sublimation of falling cloud ice'                ,phys_decomp)
  call addfld ('QVRES    ', 'kg/kg/s ', pver, 'A', 'Rate of residual condensation term'                      ,phys_decomp)
  call addfld ('CMEIOUT  ', 'kg/kg/s ', pver, 'A', 'Rate of deposition/sublimation of cloud ice'             ,phys_decomp)
  call addfld ('VTRMC    ', 'm/s     ', pver, 'A', 'Mass-weighted cloud water fallspeed'                     ,phys_decomp)
  call addfld ('VTRMI    ', 'm/s     ', pver, 'A', 'Mass-weighted cloud ice fallspeed'                       ,phys_decomp)
  call addfld ('QCSEDTEN ', 'kg/kg/s ', pver, 'A', 'Cloud water mixing ratio tendency from sedimentation'    ,phys_decomp)
  call addfld ('QISEDTEN ', 'kg/kg/s ', pver, 'A', 'Cloud ice mixing ratio tendency from sedimentation'      ,phys_decomp)
  call addfld ('PRAO     ', 'kg/kg/s ', pver, 'A', 'Accretion of cloud water by rain'                        ,phys_decomp)
  call addfld ('PRCO     ', 'kg/kg/s ', pver, 'A', 'Autoconversion of cloud water'                           ,phys_decomp)
  call addfld ('MNUCCCO  ', 'kg/kg/s ', pver, 'A', 'Immersion freezing of cloud water'                       ,phys_decomp)
  call addfld ('MNUCCTO  ', 'kg/kg/s ', pver, 'A', 'Contact freezing of cloud water'                         ,phys_decomp)
  call addfld ('MNUCCDO  ', 'kg/kg/s ', pver, 'A', 'Homogeneous and heterogeneous nucleation from vapor'     ,phys_decomp)
  call addfld ('MNUCCDOhet','kg/kg/s ', pver, 'A', 'Heterogeneous nucleation from vapor'                     ,phys_decomp)
  call addfld ('MSACWIO  ', 'kg/kg/s ', pver, 'A', 'Conversion of cloud water from rime-splintering'         ,phys_decomp)
  call addfld ('PSACWSO  ', 'kg/kg/s ', pver, 'A', 'Accretion of cloud water by snow'                        ,phys_decomp)
  call addfld ('BERGSO   ', 'kg/kg/s ', pver, 'A', 'Conversion of cloud water to snow from bergeron'         ,phys_decomp)
  call addfld ('BERGO    ', 'kg/kg/s ', pver, 'A', 'Conversion of cloud water to cloud ice from bergeron'    ,phys_decomp)
  call addfld ('MELTO    ', 'kg/kg/s ', pver, 'A', 'Melting of cloud ice'                                    ,phys_decomp)
  call addfld ('HOMOO    ', 'kg/kg/s ', pver, 'A', 'Homogeneous freezing of cloud water'                     ,phys_decomp)
  call addfld ('QCRESO   ', 'kg/kg/s ', pver, 'A', 'Residual condensation term for cloud water'              ,phys_decomp)
  call addfld ('PRCIO    ', 'kg/kg/s ', pver, 'A', 'Autoconversion of cloud ice'                             ,phys_decomp)
  call addfld ('PRAIO    ', 'kg/kg/s ', pver, 'A', 'Accretion of cloud ice by rain'                          ,phys_decomp)
  call addfld ('QIRESO   ', 'kg/kg/s ', pver, 'A', 'Residual deposition term for cloud ice'                  ,phys_decomp)
  call addfld ('MNUCCRO  ', 'kg/kg/s ', pver, 'A', 'Heterogeneous freezing of rain to snow'                  ,phys_decomp)
  call addfld ('PRACSO   ', 'kg/kg/s ', pver, 'A', 'Accretion of rain by snow'                               ,phys_decomp)
  call addfld ('MELTSDT  ', 'W/kg    ', pver, 'A', 'Latent heating rate due to melting of snow'              ,phys_decomp)
  call addfld ('FRZRDT   ', 'W/kg    ', pver, 'A', 'Latent heating rate due to homogeneous freezing of rain' ,phys_decomp)

  ! History variables for CAM5 microphysics
  call addfld ('MPDT     ', 'W/kg    ', pver, 'A', 'Heating tendency - Morrison microphysics'                ,phys_decomp)
  call addfld ('MPDQ     ', 'kg/kg/s ', pver, 'A', 'Q tendency - Morrison microphysics'                      ,phys_decomp)
  call addfld ('MPDLIQ   ', 'kg/kg/s ', pver, 'A', 'CLDLIQ tendency - Morrison microphysics'                 ,phys_decomp)
  call addfld ('MPDICE   ', 'kg/kg/s ', pver, 'A', 'CLDICE tendency - Morrison microphysics'                 ,phys_decomp)
  call addfld ('MPDW2V   ', 'kg/kg/s ', pver, 'A', 'Water <--> Vapor tendency - Morrison microphysics'       ,phys_decomp)
  call addfld ('MPDW2I   ', 'kg/kg/s ', pver, 'A', 'Water <--> Ice tendency - Morrison microphysics'         ,phys_decomp)
  call addfld ('MPDW2P   ', 'kg/kg/s ', pver, 'A', 'Water <--> Precip tendency - Morrison microphysics'      ,phys_decomp)
  call addfld ('MPDI2V   ', 'kg/kg/s ', pver, 'A', 'Ice <--> Vapor tendency - Morrison microphysics'         ,phys_decomp)
  call addfld ('MPDI2W   ', 'kg/kg/s ', pver, 'A', 'Ice <--> Water tendency - Morrison microphysics'         ,phys_decomp)
  call addfld ('MPDI2P   ', 'kg/kg/s ', pver, 'A', 'Ice <--> Precip tendency - Morrison microphysics'        ,phys_decomp)
  call addfld ('ICWNC    ', 'm-3     ', pver, 'A', 'Prognostic in-cloud water number conc'                   ,phys_decomp)
  call addfld ('ICINC    ', 'm-3     ', pver, 'A', 'Prognostic in-cloud ice number conc'                     ,phys_decomp)
  call addfld ('EFFLIQ_IND','Micron  ', pver, 'A', 'Prognostic droplet effective radius (indirect effect)'   ,phys_decomp)
  call addfld ('CDNUMC   ', '1/m2    ', 1,    'A', 'Vertically-integrated droplet concentration'             ,phys_decomp)
  call addfld ('MPICLWPI ', 'kg/m2   ', 1,    'A', 'Vertically-integrated &
       &in-cloud Initial Liquid WP (Before Micro)' ,phys_decomp)
  call addfld ('MPICIWPI ', 'kg/m2   ', 1,    'A', 'Vertically-integrated &
       &in-cloud Initial Ice WP (Before Micro)'    ,phys_decomp)

  ! Averaging for cloud particle number and size
  call addfld ('AWNC     ', 'm-3     ', pver, 'A', 'Average cloud water number conc'                         ,phys_decomp)
  call addfld ('AWNI     ', 'm-3     ', pver, 'A', 'Average cloud ice number conc'                           ,phys_decomp)
  call addfld ('AREL     ', 'Micron  ', pver, 'A', 'Average droplet effective radius'                        ,phys_decomp)
  call addfld ('AREI     ', 'Micron  ', pver, 'A', 'Average ice effective radius'                            ,phys_decomp)
  ! Frequency arrays for above
  call addfld ('FREQL    ', 'fraction', pver, 'A', 'Fractional occurrence of liquid'                          ,phys_decomp)
  call addfld ('FREQI    ', 'fraction', pver, 'A', 'Fractional occurrence of ice'                             ,phys_decomp)

  ! Average cloud top particle size and number (liq, ice) and frequency
  call addfld ('ACTREL   ', 'Micron  ', 1,    'A', 'Average Cloud Top droplet effective radius'              ,phys_decomp)
  call addfld ('ACTREI   ', 'Micron  ', 1,    'A', 'Average Cloud Top ice effective radius'                  ,phys_decomp)
  call addfld ('ACTNL    ', 'Micron  ', 1,    'A', 'Average Cloud Top droplet number'                        ,phys_decomp)
  call addfld ('ACTNI    ', 'Micron  ', 1,    'A', 'Average Cloud Top ice number'                            ,phys_decomp)

  call addfld ('FCTL     ', 'fraction', 1,    'A', 'Fractional occurrence of cloud top liquid'                ,phys_decomp)
  call addfld ('FCTI     ', 'fraction', 1,    'A', 'Fractional occurrence of cloud top ice'                   ,phys_decomp)

  call addfld ('LS_FLXPRC', 'kg/m2/s', pverp, 'A', 'ls stratiform gbm interface rain+snow flux', phys_decomp)
  call addfld ('LS_FLXSNW', 'kg/m2/s', pverp, 'A', 'ls stratiform gbm interface snow flux', phys_decomp)

  call addfld ('REL', 'micron', pver, 'A', 'MG REL stratiform cloud effective radius liquid', phys_decomp)
  call addfld ('REI', 'micron', pver, 'A', 'MG REI stratiform cloud effective radius ice', phys_decomp)
  call addfld ('LS_REFFRAIN', 'micron', pver, 'A', 'ls stratiform rain effective radius', phys_decomp)
  call addfld ('LS_REFFSNOW', 'micron', pver, 'A', 'ls stratiform snow effective radius', phys_decomp)
  call addfld ('CV_REFFLIQ', 'micron', pver, 'A', 'convective cloud liq effective radius', phys_decomp)
  call addfld ('CV_REFFICE', 'micron', pver, 'A', 'convective cloud ice effective radius', phys_decomp)

  ! diagnostic precip
  call addfld ('QRAIN   ','kg/kg   ',pver, 'A','Diagnostic grid-mean rain mixing ratio'         ,phys_decomp)
  call addfld ('QSNOW   ','kg/kg   ',pver, 'A','Diagnostic grid-mean snow mixing ratio'         ,phys_decomp)
  call addfld ('NRAIN   ','m-3     ',pver, 'A','Diagnostic grid-mean rain number conc'         ,phys_decomp)
  call addfld ('NSNOW   ','m-3     ',pver, 'A','Diagnostic grid-mean snow number conc'         ,phys_decomp)

  ! size of precip
  call addfld ('RERCLD   ','m      ',pver, 'A','Diagnostic effective radius of Liquid Cloud and Rain' ,phys_decomp)
  call addfld ('DSNOW   ','m       ',pver, 'A','Diagnostic grid-mean snow diameter'         ,phys_decomp)

  ! diagnostic radar reflectivity, cloud-averaged
  call addfld ('REFL  ','DBz  ',pver, 'A','94 GHz radar reflectivity'       ,phys_decomp)
  call addfld ('AREFL  ','DBz  ',pver, 'A','Average 94 GHz radar reflectivity'       ,phys_decomp)
  call addfld ('FREFL  ','fraction  ',pver, 'A','Fractional occurrence of radar reflectivity'       ,phys_decomp)

  call addfld ('CSRFL  ','DBz  ',pver, 'A','94 GHz radar reflectivity (CloudSat thresholds)'       ,phys_decomp)
  call addfld ('ACSRFL  ','DBz  ',pver, 'A','Average 94 GHz radar reflectivity (CloudSat thresholds)'       ,phys_decomp)
  call addfld ('FCSRFL  ','fraction  ',pver, 'A','Fractional occurrence of radar reflectivity (CloudSat thresholds)' &
       ,phys_decomp)

  call addfld ('AREFLZ ','mm^6/m^3 ',pver, 'A','Average 94 GHz radar reflectivity'       ,phys_decomp)

  ! Aerosol information
  call addfld ('NCAL    ','1/m3   ',pver, 'A','Number Concentation Activated for Liquid',phys_decomp)
  call addfld ('NCAI    ','1/m3   ',pver, 'A','Number Concentation Activated for Ice',phys_decomp)

  ! Average rain and snow mixing ratio (Q), number (N) and diameter (D), with frequency
  call addfld ('AQRAIN   ','kg/kg   ',pver, 'A','Average rain mixing ratio'         ,phys_decomp)
  call addfld ('AQSNOW   ','kg/kg   ',pver, 'A','Average snow mixing ratio'         ,phys_decomp)
  call addfld ('ANRAIN   ','m-3     ',pver, 'A','Average rain number conc'         ,phys_decomp)
  call addfld ('ANSNOW   ','m-3     ',pver, 'A','Average snow number conc'         ,phys_decomp)
  call addfld ('ADRAIN   ','Micron  ',pver, 'A','Average rain effective Diameter'         ,phys_decomp)
  call addfld ('ADSNOW   ','Micron  ',pver, 'A','Average snow effective Diameter'         ,phys_decomp)
  call addfld ('FREQR  ','fraction  ',pver, 'A','Fractional occurrence of rain'       ,phys_decomp)
  call addfld ('FREQS  ','fraction  ',pver, 'A','Fractional occurrence of snow'       ,phys_decomp)

  ! precipitation efficiency & other diagnostic fields
  call addfld('PE'    , '1',       1, 'A', 'Stratiform Precipitation Efficiency  (precip/cmeliq)',       phys_decomp )
  call addfld('APRL'  , 'm/s',     1, 'A', 'Average Stratiform Precip Rate over efficiency calculation', phys_decomp )
  call addfld('PEFRAC', '1',       1, 'A', 'Fraction of timesteps precip efficiency reported',           phys_decomp )
  call addfld('VPRCO' , 'kg/kg/s', 1, 'A', 'Vertical average of accretion rate',                         phys_decomp )
  call addfld('VPRAO' , 'kg/kg/s', 1, 'A', 'Vertical average of autoconversion rate',                    phys_decomp )
  call addfld('RACAU' , 'kg/kg/s', 1, 'A', 'Accretion/autoconversion ratio from vertical average',       phys_decomp )

  ! determine the add_default fields
  call phys_getopts(history_amwg_out           = history_amwg         , &
       history_budget_out         = history_budget       , &
       history_budget_histfile_num_out = budget_histfile)

  if (history_amwg) then
     call add_default ('FICE    ', 1, ' ')
     call add_default ('AQRAIN   ', 1, ' ')
     call add_default ('AQSNOW   ', 1, ' ')
     call add_default ('ANRAIN   ', 1, ' ')
     call add_default ('ANSNOW   ', 1, ' ')
     call add_default ('AREI     ', 1, ' ')
     call add_default ('AREL     ', 1, ' ')
     call add_default ('AWNC     ', 1, ' ')
     call add_default ('AWNI     ', 1, ' ')
     call add_default ('CDNUMC   ', 1, ' ')
     call add_default ('FREQR    ', 1, ' ')
     call add_default ('FREQS    ', 1, ' ')
     call add_default ('FREQL    ', 1, ' ')
     call add_default ('FREQI    ', 1, ' ')
     do m = 1, ncnst
        call cnst_get_ind(cnst_names(m), mm)
        call add_default(cnst_name(mm), 1, ' ')
        ! call add_default(sflxnam(mm),   1, ' ')
     end do
  end if

  if ( history_budget ) then
     call add_default ('EVAPSNOW ', budget_histfile, ' ')
     call add_default ('EVAPPREC ', budget_histfile, ' ')
     call add_default ('QVRES    ', budget_histfile, ' ')
     call add_default ('QISEVAP  ', budget_histfile, ' ')
     call add_default ('QCSEVAP  ', budget_histfile, ' ')
     call add_default ('QISEDTEN ', budget_histfile, ' ')
     call add_default ('QCSEDTEN ', budget_histfile, ' ')
     call add_default ('QIRESO   ', budget_histfile, ' ')
     call add_default ('QCRESO   ', budget_histfile, ' ')
     call add_default ('PSACWSO  ', budget_histfile, ' ')
     call add_default ('PRCO     ', budget_histfile, ' ')
     call add_default ('PRCIO    ', budget_histfile, ' ')
     call add_default ('PRAO     ', budget_histfile, ' ')
     call add_default ('PRAIO    ', budget_histfile, ' ')
     call add_default ('PRACSO   ', budget_histfile, ' ')
     call add_default ('MSACWIO  ', budget_histfile, ' ')
     call add_default ('MPDW2V   ', budget_histfile, ' ')
     call add_default ('MPDW2P   ', budget_histfile, ' ')
     call add_default ('MPDW2I   ', budget_histfile, ' ')
     call add_default ('MPDT     ', budget_histfile, ' ')
     call add_default ('MPDQ     ', budget_histfile, ' ')
     call add_default ('MPDLIQ   ', budget_histfile, ' ')
     call add_default ('MPDICE   ', budget_histfile, ' ')
     call add_default ('MPDI2W   ', budget_histfile, ' ')
     call add_default ('MPDI2V   ', budget_histfile, ' ')
     call add_default ('MPDI2P   ', budget_histfile, ' ')
     call add_default ('MNUCCTO  ', budget_histfile, ' ')
     call add_default ('MNUCCRO  ', budget_histfile, ' ')
     call add_default ('MNUCCCO  ', budget_histfile, ' ')
     call add_default ('MELTSDT  ', budget_histfile, ' ')
     call add_default ('MELTO    ', budget_histfile, ' ')
     call add_default ('HOMOO    ', budget_histfile, ' ')
     call add_default ('FRZRDT   ', budget_histfile, ' ')
     call add_default ('CMEIOUT  ', budget_histfile, ' ')
     call add_default ('BERGSO   ', budget_histfile, ' ')
     call add_default ('BERGO    ', budget_histfile, ' ')

     call add_default(cnst_name(ixcldliq), budget_histfile, ' ')
     call add_default(cnst_name(ixcldice), budget_histfile, ' ')
     call add_default(apcnst   (ixcldliq), budget_histfile, ' ')
     call add_default(apcnst   (ixcldice), budget_histfile, ' ')
     call add_default(bpcnst   (ixcldliq), budget_histfile, ' ')
     call add_default(bpcnst   (ixcldice), budget_histfile, ' ')

  end if

  ! physics buffer indices
  ast_idx      = pbuf_get_index('AST')
  aist_idx     = pbuf_get_index('AIST')
  alst_idx     = pbuf_get_index('ALST')
  cld_idx      = pbuf_get_index('CLD')
  concld_idx   = pbuf_get_index('CONCLD')

  naai_idx     = pbuf_get_index('NAAI')
  naai_hom_idx = pbuf_get_index('NAAI_HOM')
  npccn_idx    = pbuf_get_index('NPCCN')
  rndst_idx    = pbuf_get_index('RNDST')
  nacon_idx    = pbuf_get_index('NACON')

  prec_str_idx = pbuf_get_index('PREC_STR')
  snow_str_idx = pbuf_get_index('SNOW_STR')
  prec_sed_idx = pbuf_get_index('PREC_SED')
  snow_sed_idx = pbuf_get_index('SNOW_SED')
  prec_pcw_idx = pbuf_get_index('PREC_PCW')
  snow_pcw_idx = pbuf_get_index('SNOW_PCW')

  cmeliq_idx   = pbuf_get_index('CMELIQ')

  ! Initialize physics buffer fields for accumulating precip and condensation
  if (is_first_step()) then
     call pbuf_set_field(pbuf2d, cldo_idx,   0._r8)
     call pbuf_set_field(pbuf2d, cc_t_idx,   0._r8)
     call pbuf_set_field(pbuf2d, cc_qv_idx,  0._r8)
     call pbuf_set_field(pbuf2d, cc_ql_idx,  0._r8)
     call pbuf_set_field(pbuf2d, cc_qi_idx,  0._r8)
     call pbuf_set_field(pbuf2d, cc_nl_idx,  0._r8)
     call pbuf_set_field(pbuf2d, cc_ni_idx,  0._r8)
     call pbuf_set_field(pbuf2d, cc_qlst_idx,0._r8)
     call pbuf_set_field(pbuf2d, acpr_idx,   0._r8)
     call pbuf_set_field(pbuf2d, acgcme_idx, 0._r8)
     call pbuf_set_field(pbuf2d, acnum_idx,  0)
     call pbuf_set_field(pbuf2d, relvar_idx, 2._r8)
     call pbuf_set_field(pbuf2d, accre_enhan_idx, 1._r8)
  end if

end subroutine micro_mg_cam_init


!===============================================================================

subroutine micro_mg_cam_tend(state, ptend, dtime, pbuf)

  use micro_mg1_0, only: micro_mg_tend1_0 => micro_mg_tend
  use micro_mg1_5, only: micro_mg_tend1_5 => micro_mg_tend, &
       micro_mg_get_cols1_5 => micro_mg_get_cols

  type(physics_state),         intent(in)    :: state
  type(physics_ptend),         intent(out)   :: ptend
  real(r8),                    intent(in)    :: dtime
  type(physics_buffer_desc),   pointer       :: pbuf(:)

  ! Local variables
  logical :: microp_uniform = .false. ! True = configure microphysics for sub-columns
  ! False = use in regular mode w/o sub-columns
  integer :: lchnk, ncol

  integer :: i, k, itim_old, it

  real(r8), pointer :: naai(:,:)      ! ice nucleation number
  real(r8), pointer :: naai_hom(:,:)  ! ice nucleation number (homogeneous)
  real(r8), pointer :: npccn(:,:)     ! liquid activation number tendency
  real(r8), pointer :: rndst(:,:,:)
  real(r8), pointer :: nacon(:,:,:)

  real(r8), pointer :: prec_str(:)          ! [Total] Sfc flux of precip from stratiform [ m/s ]
  real(r8), pointer :: snow_str(:)          ! [Total] Sfc flux of snow from stratiform   [ m/s ]
  real(r8), pointer :: prec_sed(:)          ! Surface flux of total cloud water from sedimentation
  real(r8), pointer :: snow_sed(:)          ! Surface flux of cloud ice from sedimentation
  real(r8), pointer :: prec_pcw(:)          ! Sfc flux of precip from microphysics [ m/s ]
  real(r8), pointer :: snow_pcw(:)          ! Sfc flux of snow from microphysics [ m/s ]

  real(r8), pointer :: ast(:,:)          ! Relative humidity cloud fraction
  real(r8), pointer :: alst_mic(:,:)
  real(r8), pointer :: aist_mic(:,:)
  real(r8), pointer :: cldo(:,:)         ! Old cloud fraction
  real(r8), pointer :: nevapr(:,:)       ! Evaporation of total precipitation (rain + snow)
  real(r8), pointer :: relvar(:,:)       ! relative variance of cloud water
  real(r8), pointer :: accre_enhan(:,:)  ! optional accretion enhancement for experimentation
  real(r8), pointer :: prain(:,:)        ! Total precipitation (rain + snow)
  real(r8), pointer :: dei(:,:)          ! Ice effective diameter (meters) (AG: microns?)
  real(r8), pointer :: mu(:,:)           ! Size distribution shape parameter for radiation
  real(r8), pointer :: lambdac(:,:)      ! Size distribution slope parameter for radiation
  real(r8), pointer :: des(:,:)          ! Snow effective diameter (m)

  real(r8) :: rate1cld(pcols,pver) ! array to hold rate1ord_cw2pr_st from microphysics

  real(r8) :: tlat(pcols,pver)
  real(r8) :: qvlat(pcols,pver)
  real(r8) :: qcten(pcols,pver)
  real(r8) :: qiten(pcols,pver)
  real(r8) :: ncten(pcols,pver)
  real(r8) :: niten(pcols,pver)
  real(r8) :: prect(pcols)
  real(r8) :: preci(pcols)


  real(r8) :: evapsnow(pcols,pver)                    ! Local evaporation of snow
  real(r8) :: prodsnow(pcols,pver)                    ! Local production of snow
  real(r8) :: cmeice(pcols,pver)                      ! Rate of cond-evap of ice within the cloud
  real(r8) :: qsout(pcols,pver)                       ! Snow mixing ratio
  real(r8) :: rflx(pcols,pver+1)                   ! grid-box average rain flux (kg m^-2 s^-1)
  real(r8) :: sflx(pcols,pver+1)                   ! grid-box average snow flux (kg m^-2 s^-1)
  real(r8) :: qrout(pcols,pver)                     ! Rain mixing ratio
  real(r8) :: reff_rain(pcols,pver)                ! rain effective radius (um)
  real(r8) :: reff_snow(pcols,pver)                ! snow effective radius (um)
  real(r8) :: qcsevap(pcols,pver)                     ! Evaporation of falling cloud water
  real(r8) :: qisevap(pcols,pver)                     ! Sublimation of falling cloud ice
  real(r8) :: qvres(pcols,pver)                       ! Residual condensation term to remove excess saturation
  real(r8) :: cmeiout(pcols,pver)                     ! Deposition/sublimation rate of cloud ice
  real(r8) :: vtrmc(pcols,pver)                       ! Mass-weighted cloud water fallspeed
  real(r8) :: vtrmi(pcols,pver)                       ! Mass-weighted cloud ice fallspeed
  real(r8) :: qcsedten(pcols,pver)                    ! Cloud water mixing ratio tendency from sedimentation
  real(r8) :: qisedten(pcols,pver)                    ! Cloud ice mixing ratio tendency from sedimentation
  real(r8) :: prao(pcols,pver)
  real(r8) :: prco(pcols,pver)
  real(r8) :: mnuccco(pcols,pver)
  real(r8) :: mnuccto(pcols,pver)
  real(r8) :: msacwio(pcols,pver)
  real(r8) :: psacwso(pcols,pver)
  real(r8) :: bergso(pcols,pver)
  real(r8) :: bergo(pcols,pver)
  real(r8) :: melto(pcols,pver)
  real(r8) :: homoo(pcols,pver)
  real(r8) :: qcreso(pcols,pver)
  real(r8) :: prcio(pcols,pver)
  real(r8) :: praio(pcols,pver)
  real(r8) :: qireso(pcols,pver)
  real(r8) :: mnuccro(pcols,pver)
  real(r8) :: pracso (pcols,pver)
  real(r8) :: meltsdt(pcols,pver)
  real(r8) :: frzrdt (pcols,pver)
  real(r8) :: mnuccdo(pcols,pver)
  real(r8) :: nrout(pcols,pver)
  real(r8) :: nsout(pcols,pver)
  real(r8) :: refl(pcols,pver)   ! analytic radar reflectivity
  real(r8) :: arefl(pcols,pver)  !average reflectivity will zero points outside valid range
  real(r8) :: areflz(pcols,pver) !average reflectivity in z.
  real(r8) :: frefl(pcols,pver)
  real(r8) :: csrfl(pcols,pver)  !cloudsat reflectivity
  real(r8) :: acsrfl(pcols,pver) !cloudsat average
  real(r8) :: fcsrfl(pcols,pver)
  real(r8) :: rercld(pcols,pver) ! effective radius calculation for rain + cloud
  real(r8) :: ncai(pcols,pver)   ! output number conc of ice nuclei available (1/m3)
  real(r8) :: ncal(pcols,pver)   ! output number conc of CCN (1/m3)
  real(r8) :: qrout2(pcols,pver)
  real(r8) :: qsout2(pcols,pver)
  real(r8) :: nrout2(pcols,pver)
  real(r8) :: nsout2(pcols,pver)
  real(r8) :: drout2(pcols,pver) ! mean rain particle diameter (m)
  real(r8) :: dsout2(pcols,pver) ! mean snow particle diameter (m)
  real(r8) :: freqs(pcols,pver)
  real(r8) :: freqr(pcols,pver)
  real(r8) :: nfice(pcols,pver)

  real(r8) :: mnuccdohet(pcols,pver)

  ! physics buffer fields for COSP simulator
  real(r8), pointer :: mgflxprc(:,:)     ! MG grid-box mean flux_large_scale_cloud_rain+snow at interfaces (kg/m2/s)
  real(r8), pointer :: mgflxsnw(:,:)     ! MG grid-box mean flux_large_scale_cloud_snow at interfaces (kg/m2/s)
  real(r8), pointer :: mgmrprc(:,:)      ! MG grid-box mean mixingratio_large_scale_cloud_rain+snow at interfaces (kg/kg)
  real(r8), pointer :: mgmrsnw(:,:)      ! MG grid-box mean mixingratio_large_scale_cloud_snow at interfaces (kg/kg)
  real(r8), pointer :: mgreffrain(:,:)   ! MG diagnostic rain effective radius (um)
  real(r8), pointer :: mgreffsnow(:,:)   ! MG diagnostic snow effective radius (um)
  real(r8), pointer :: cvreffliq(:,:)    ! convective cloud liquid effective radius (um)
  real(r8), pointer :: cvreffice(:,:)    ! convective cloud ice effective radius (um)

  ! physics buffer fields used with CARMA
  real(r8), pointer, dimension(:,:) :: tnd_qsnow    ! external tendency on snow mass (kg/kg/s)
  real(r8), pointer, dimension(:,:) :: tnd_nsnow    ! external tendency on snow number(#/kg/s)
  real(r8), pointer, dimension(:,:) :: re_ice       ! ice effective radius (m)

  real(r8), pointer :: rate1ord_cw2pr_st(:,:) ! 1st order rate for direct conversion of
  ! strat. cloud water to precip (1/s)    ! rce 2010/05/01
  real(r8), pointer :: wsedl(:,:)        ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]

  ! For rrtm optics. specificed distribution.
  real(r8) :: mucon                                 ! Convective size distribution shape parameter
  real(r8) :: dcon                                  ! Convective size distribution effective radius (meters)
  real(r8) :: deicon                                ! Convective ice effective diameter (meters)


  real(r8), pointer :: CC_T(:,:)         ! Grid-mean microphysical tendency
  real(r8), pointer :: CC_qv(:,:)        ! Grid-mean microphysical tendency
  real(r8), pointer :: CC_ql(:,:)        ! Grid-mean microphysical tendency
  real(r8), pointer :: CC_qi(:,:)        ! Grid-mean microphysical tendency
  real(r8), pointer :: CC_nl(:,:)        ! Grid-mean microphysical tendency
  real(r8), pointer :: CC_ni(:,:)        ! Grid-mean microphysical tendency
  real(r8), pointer :: CC_qlst(:,:)      ! In-liquid stratus microphysical tendency

  real(r8), pointer :: qme(:,:)

  ! A local copy of state is used for diagnostic calculations
  type(physics_state) :: state_loc
  type(physics_ptend) :: ptend_loc

  real(r8) :: icecldf(pcols,pver)                     ! Ice cloud fraction
  real(r8) :: liqcldf(pcols,pver)                     ! Liquid cloud fraction (combined into cloud)

  real(r8), pointer :: rel(:,:)          ! Liquid effective drop radius (microns)
  real(r8), pointer :: rei(:,:)          ! Ice effective drop size (microns)
  real(r8) :: rel_fn(pcols,pver)         ! Ice effective drop size at fixed number (indirect effect) (microns)

  ! in-cloud water quantities adjusted for convective water
  real(r8) :: allcld_ice(pcols,pver)                 ! All-cloud cloud ice
  real(r8) :: allcld_liq(pcols,pver)                 ! All-cloud liquid

  real(r8), pointer :: cmeliq(:,:)

  real(r8), pointer :: cld(:,:)          ! Total cloud fraction
  real(r8), pointer :: concld(:,:)       ! Convective cloud fraction
  real(r8), pointer :: iciwpst(:,:)      ! Stratiform in-cloud ice water path for radiation
  real(r8), pointer :: iclwpst(:,:)      ! Stratiform in-cloud liquid water path for radiation
  real(r8), pointer :: cldfsnow(:,:)     ! Cloud fraction for liquid+snow
  real(r8), pointer :: icswp(:,:)        ! In-cloud snow water path

  real(r8) :: icimrst(pcols,pver)                     ! In stratus ice mixing ratio
  real(r8) :: icwmrst(pcols,pver)                     ! In stratus water mixing ratio
  real(r8) :: icinc(pcols,pver)                       ! In cloud ice number conc
  real(r8) :: icwnc(pcols,pver)                       ! In cloud water number conc

  real(r8) :: cdnumc(pcols)                           ! Vertically-integrated droplet concentration
  real(r8) :: iclwpi(pcols)                           ! Vertically-integrated in-cloud Liquid WP before microphysics
  real(r8) :: iciwpi(pcols)                           ! Vertically-integrated in-cloud Ice WP before microphysics

  ! Averaging arrays for effective radius and number....
  real(r8) :: efiout(pcols,pver)
  real(r8) :: efcout(pcols,pver)
  real(r8) :: ncout(pcols,pver)
  real(r8) :: niout(pcols,pver)
  real(r8) :: freqi(pcols,pver)
  real(r8) :: freql(pcols,pver)

  real(r8) :: icecldf_out(pcols,pver)                 ! Ice cloud fraction
  real(r8) :: liqcldf_out(pcols,pver)                 ! Liquid cloud fraction (combined into cloud)
  real(r8) :: icimrst_out(pcols,pver)                 ! In stratus ice mixing ratio
  real(r8) :: icwmrst_out(pcols,pver)                 ! In stratus water mixing ratio

  ! Average cloud top radius & number
  real(r8) :: ctrel(pcols)
  real(r8) :: ctrei(pcols)
  real(r8) :: ctnl(pcols)
  real(r8) :: ctni(pcols)
  real(r8) :: fcti(pcols)
  real(r8) :: fctl(pcols)

  real(r8) :: ftem(pcols,pver)

  ! Variables for precip efficiency calculation
  real(r8) :: minlwp        ! LWP threshold
  real(r8) :: rhow          ! density of water (1000 kg/m3)

  real(r8), pointer, dimension(:) :: acprecl ! accumulated precip across timesteps
  real(r8), pointer, dimension(:) :: acgcme  ! accumulated condensation across timesteps
  integer,  pointer, dimension(:) :: acnum   ! counter for # timesteps accumulated

  ! Variables for liquid water path and column condensation
  real(r8) :: tgliqwp(pcols)   ! column liquid
  real(r8) :: tgcmeliq(pcols)  ! column condensation rate (units)

  real(r8) :: pe(pcols)        ! precip efficiency for output
  real(r8) :: pefrac(pcols)    ! fraction of time precip efficiency is written out
  real(r8) :: tpr(pcols)       ! average accumulated precipitation rate in pe calculation

  ! variables for autoconversion and accretion vertical averages
  real(r8) :: vprco(pcols)     ! vertical average accretion
  real(r8) :: vprao(pcols)     ! vertical average autoconversion
  real(r8) :: racau(pcols)     ! ratio of vertical averages
  integer  :: cnt(pcols)       ! counters
  logical  :: lq(pcnst)

  real(r8) :: qc(pcols,pver)    ! cloud water mixing ratio (kg/kg)
  real(r8) :: qi(pcols,pver)    ! cloud ice mixing ratio (kg/kg)
  real(r8) :: nc(pcols,pver)    ! cloud water number conc (1/kg)
  real(r8) :: ni(pcols,pver)    ! cloud ice number conc (1/kg)

  integer :: nlev   ! number of levels where cloud physics is done
  integer :: mgncol ! size of mgcols
  integer, allocatable :: mgcols(:) ! Columns with microphysics performed

  character(128) :: errstring   ! return status (non-blank for error return)

  !-------------------------------------------------------------------------------

  ! Find the number of levels used in the microphysics.
  nlev  = pver - top_lev + 1

  lchnk = state%lchnk
  ncol  = state%ncol

  ! Physics buffer fields needed in micro_mg_tend
  call pbuf_get_field(pbuf, naai_idx, naai)
  call pbuf_get_field(pbuf, naai_hom_idx, naai_hom)
  call pbuf_get_field(pbuf, npccn_idx, npccn)
  call pbuf_get_field(pbuf, rndst_idx, rndst)
  call pbuf_get_field(pbuf, nacon_idx, nacon)

  call pbuf_get_field(pbuf, prec_str_idx, prec_str)
  call pbuf_get_field(pbuf, snow_str_idx, snow_str)
  call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw)
  call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw)
  call pbuf_get_field(pbuf, prec_sed_idx, prec_sed)
  call pbuf_get_field(pbuf, snow_sed_idx, snow_sed)

  itim_old = pbuf_old_tim_idx()
  call pbuf_get_field(pbuf, ast_idx,    ast,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
  call pbuf_get_field(pbuf, cldo_idx,   cldo,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
  call pbuf_get_field(pbuf, nevapr_idx, nevapr)
  call pbuf_get_field(pbuf, prain_idx, prain)
  call pbuf_get_field(pbuf, dei_idx, dei)
  call pbuf_get_field(pbuf, mu_idx, mu)
  call pbuf_get_field(pbuf, lambdac_idx, lambdac)
  call pbuf_get_field(pbuf, des_idx, des)
  call pbuf_get_field(pbuf, tnd_qsnow_idx, tnd_qsnow)
  call pbuf_get_field(pbuf, tnd_nsnow_idx, tnd_nsnow)
  call pbuf_get_field(pbuf, re_ice_idx, re_ice)
  call pbuf_get_field(pbuf, relvar_idx, relvar)
  call pbuf_get_field(pbuf, accre_enhan_idx, accre_enhan)

  ! Effective droplet radius
  call pbuf_get_field(pbuf, rel_idx,    rel    )
  call pbuf_get_field(pbuf, rei_idx,    rei    )

  ! Microphysics assumes 'liquid stratus frac = ice stratus frac
  !                      = max( liquid stratus frac, ice stratus frac )'.
  alst_mic => ast
  aist_mic => ast

  ! Output initial in-cloud LWP (before microphysics)

  iclwpi = 0._r8
  iciwpi = 0._r8

  do i = 1, ncol
     do k = top_lev, pver
        iclwpi(i) = iclwpi(i) + &
             min(state%q(i,k,ixcldliq) / max(0.0001_r8,ast(i,k)),0.005_r8) &
             * state%pdel(i,k) / gravit
        iciwpi(i) = iciwpi(i) + &
             min(state%q(i,k,ixcldice) / max(0.0001_r8,ast(i,k)),0.005_r8) &
             * state%pdel(i,k) / gravit
     end do
  end do

  call outfld('MPICLWPI', iclwpi, pcols, lchnk)
  call outfld('MPICIWPI', iciwpi, pcols, lchnk)

  ! Probably unnecessary. Leaving it just in case.
  cldo(:ncol,top_lev:pver)=ast(:ncol,top_lev:pver)

  ! Initialize local state from input.
  call physics_state_copy(state, state_loc)

  ! Initialize ptend for output.

  lq = .false.
  lq(1) = .true.
  lq(ixcldliq) = .true.
  lq(ixcldice) = .true.
  lq(ixnumliq) = .true.
  lq(ixnumice) = .true.

  ! the name 'cldwat' triggers special tests on cldliq
  ! and cldice in physics_update
  call physics_ptend_init(ptend, state%psetcols, "cldwat", &
       ls=.true., lq=lq)

  select case (micro_mg_version)
  case (1)
     select case (micro_mg_sub_version)
     case (0)

        qc = state_loc%q(:,:,ixcldliq)
        qi = state_loc%q(:,:,ixcldice)
        nc = state_loc%q(:,:,ixnumliq)
        ni = state_loc%q(:,:,ixnumice)

        call micro_mg_tend1_0( &
             microp_uniform, pcols, pver, ncol, top_lev, dtime, &
             state_loc%t, state_loc%q(:,:,1), qc, qi, nc,     &
             ni, state_loc%pmid, state_loc%pdel, ast, alst_mic,&
             relvar, accre_enhan,                             &
             aist_mic, rate1cld, naai, npccn,                 &
             rndst, nacon, tlat, qvlat, qcten,                &
             qiten, ncten, niten, rel, rel_fn,                &
             rei, prect, preci, nevapr, evapsnow,             &
             prain, prodsnow, cmeice, dei, mu,                &
             lambdac, qsout, des, rflx, sflx,                 &
             qrout, reff_rain, reff_snow, qcsevap, qisevap,   &
             qvres, cmeiout, vtrmc, vtrmi, qcsedten,          &
             qisedten, prao, prco, mnuccco, mnuccto,          &
             msacwio, psacwso, bergso, bergo, melto,          &
             homoo, qcreso, prcio, praio, qireso,             &
             mnuccro, pracso, meltsdt, frzrdt, mnuccdo,       &
             nrout, nsout, refl, arefl, areflz,               &
             frefl, csrfl, acsrfl, fcsrfl, rercld,            &
             ncai, ncal, qrout2, qsout2, nrout2,              &
             nsout2, drout2, dsout2, freqs, freqr,            &
             nfice, do_cldice, tnd_qsnow,                     &
             tnd_nsnow, re_ice, errstring)


     case (5)

        call micro_mg_get_cols1_5(ncol, nlev, top_lev, state%q(:,:,ixcldliq), &
             state%q(:,:,ixcldice), mgncol, mgcols)

        call micro_mg_tend1_5( &
             mgncol,   mgcols,   nlev,     top_lev,  dtime,              &
             state_loc%t,        state_loc%q(:,:,1),                     &
             state_loc%q(:,:,ixcldliq),    state_loc%q(:,:,ixcldice),    &
             state_loc%q(:,:,ixnumliq),    state_loc%q(:,:,ixnumice),    &
             relvar,             accre_enhan,                            &
             state_loc%pmid,     state_loc%pdel,     state_loc%pint,     &
             ast,                alst_mic,           aist_mic,           &
             rate1cld,           naai,     npccn,    rndst,    nacon,    &
             tlat,     qvlat,    qcten,    qiten,    ncten,    niten,    &
             rel,     rel_fn,  rei,               prect,    preci,    &
             nevapr,   evapsnow, prain,    prodsnow, cmeice,   dei,      &
             mu,       lambdac,  qsout,    des,      rflx,     sflx,     &
             qrout,              reff_rain,          reff_snow,          &
             qcsevap,  qisevap,  qvres,    cmeiout,  vtrmc,    vtrmi,    &
             qcsedten, qisedten, prao,     prco,     mnuccco,  mnuccto,  &
             msacwio,  psacwso,  bergso,   bergo,    melto,    homoo,    &
             qcreso,             prcio,    praio,    qireso,             &
             mnuccro,  pracso,   meltsdt,  frzrdt,   mnuccdo,            &
             nrout,    nsout,    refl,     arefl,    areflz,   frefl,    &
             csrfl,    acsrfl,   fcsrfl,             rercld,             &
             ncai,     ncal,     qrout2,   qsout2,   nrout2,   nsout2,   &
             drout2,   dsout2,   freqs,    freqr,    nfice,              &
             tnd_qsnow,          tnd_nsnow,          re_ice,             &
             errstring)

        call handle_errmsg(errstring, subname="micro_mg_tend1_5")
     end select
  end select

  call handle_errmsg(errstring, subname="micro_mg_tend")

  call physics_ptend_init(ptend_loc, state%psetcols, "micro_mg", &
       ls=.true., lq=lq)

  ! Set local tendency.
  ptend_loc%s(:ncol,top_lev:pver)          =  tlat(:ncol,top_lev:pver)
  ptend_loc%q(:ncol,top_lev:pver,1)        = qvlat(:ncol,top_lev:pver)
  ptend_loc%q(:ncol,top_lev:pver,ixcldliq) = qcten(:ncol,top_lev:pver)
  ptend_loc%q(:ncol,top_lev:pver,ixcldice) = qiten(:ncol,top_lev:pver)
  ptend_loc%q(:ncol,top_lev:pver,ixnumliq) = ncten(:ncol,top_lev:pver)
  ptend_loc%q(:ncol,top_lev:pver,ixnumice) = niten(:ncol,top_lev:pver)

  ! Sum into overall ptend
  call physics_ptend_sum(ptend_loc, ptend, ncol)

  ! Update local state
  call physics_update(state_loc, ptend_loc, dtime)

  ! Check to make sure that the microphysics code is respecting the flags that control
  ! whether MG should be prognosing cloud ice and cloud liquid or not.
  if (.not. do_cldice) then
     if (any(ptend%q(:ncol,top_lev:pver,ixcldice) /= 0.0_r8)) &
          call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
          " but micro_mg_tend has ice mass tendencies.")
     if (any(ptend%q(:ncol,top_lev:pver,ixnumice) /= 0.0_r8)) &
          call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
          " but micro_mg_tend has ice number tendencies.")
  end if
  if (.not. do_cldliq) then
     if (any(ptend%q(:ncol,top_lev:pver,ixcldliq) /= 0.0_r8)) &
          call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud liquid,"// &
          " but micro_mg_tend has liquid mass tendencies.")
     if (any(ptend%q(:ncol,top_lev:pver,ixnumliq) /= 0.0_r8)) &
          call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud liquid,"// &
          " but micro_mg_tend has liquid number tendencies.")
  end if

  call outfld('QRAIN',   qrout, pcols, lchnk)
  call outfld('QSNOW',   qsout, pcols, lchnk)
  call outfld('NRAIN',   nrout, pcols, lchnk)
  call outfld('NSNOW',   nsout, pcols, lchnk)
  call outfld('DSNOW',     des, pcols, lchnk)
  call outfld('REFL',     refl, pcols, lchnk)
  call outfld('AREFL',   arefl, pcols, lchnk)
  call outfld('AREFLZ', areflz, pcols, lchnk)
  call outfld('FREFL',   frefl, pcols, lchnk)
  call outfld('CSRFL',   csrfl, pcols, lchnk)
  call outfld('ACSRFL', acsrfl, pcols, lchnk)
  call outfld('FCSRFL', fcsrfl, pcols, lchnk)
  call outfld('RERCLD', rercld, pcols, lchnk)
  call outfld('NCAL',     ncal, pcols, lchnk)
  call outfld('NCAI',     ncai, pcols, lchnk)
  call outfld('AQRAIN', qrout2, pcols, lchnk)
  call outfld('AQSNOW', qsout2, pcols, lchnk)
  call outfld('ANRAIN', nrout2, pcols, lchnk)
  call outfld('ANSNOW', nsout2, pcols, lchnk)
  call outfld('ADRAIN', drout2, pcols, lchnk)
  call outfld('ADSNOW', dsout2, pcols, lchnk)
  call outfld('FREQR',   freqr, pcols, lchnk)
  call outfld('FREQS',   freqs, pcols, lchnk)
  call outfld('FICE',    nfice, pcols, lchnk)

  mnuccdohet = 0._r8
  do k=top_lev,pver
     do i=1,ncol
        if (naai(i,k) > 0._r8) then
           mnuccdohet(i,k) = mnuccdo(i,k) - (naai_hom(i,k)/naai(i,k))*mnuccdo(i,k)
        end if
     end do
  end do

  call pbuf_get_field(pbuf, ls_flxprc_idx, mgflxprc  )
  call pbuf_get_field(pbuf, ls_flxsnw_idx, mgflxsnw  )
  mgflxprc(:ncol,top_lev:pverp) = rflx(:ncol,top_lev:pverp) + sflx(:ncol,top_lev:pverp)
  mgflxsnw(:ncol,top_lev:pverp) = sflx(:ncol,top_lev:pverp)

  call pbuf_get_field(pbuf, ls_mrprc_idx, mgmrprc  )
  call pbuf_get_field(pbuf, ls_mrsnw_idx, mgmrsnw  )
  mgmrprc(:ncol,top_lev:pver) = qrout(:ncol,top_lev:pver) + qsout(:ncol,top_lev:pver)
  mgmrsnw(:ncol,top_lev:pver) = qsout(:ncol,top_lev:pver)

  call pbuf_get_field(pbuf, ls_reffrain_idx, mgreffrain  )
  call pbuf_get_field(pbuf, ls_reffsnow_idx, mgreffsnow  )
  mgreffrain(:ncol,top_lev:pver) = reff_rain(:ncol,top_lev:pver)
  mgreffsnow(:ncol,top_lev:pver) = reff_snow(:ncol,top_lev:pver)

  call pbuf_get_field(pbuf, cv_reffliq_idx, cvreffliq  )
  call pbuf_get_field(pbuf, cv_reffice_idx, cvreffice  )
  !! calculate effective radius of convective liquid and ice using dcon and deicon (not used by code, not useful for COSP)
  !! hard-coded as average of hard-coded values used for deep/shallow convective detrainment (near line 1502/1505)
  cvreffliq(:ncol,top_lev:pver) = 9.0_r8
  cvreffice(:ncol,top_lev:pver) = 37.0_r8

  call outfld( 'LS_REFFRAIN'  , mgreffrain,       pcols, lchnk )
  call outfld( 'LS_REFFSNOW'  , mgreffsnow,       pcols, lchnk )
  call outfld( 'CV_REFFLIQ'   , cvreffliq,       pcols, lchnk )
  call outfld( 'CV_REFFICE'   , cvreffice,       pcols, lchnk )

  ! Reassign rate1 if modal aerosols
  if (rate1_cw2pr_st_idx > 0) then
     call pbuf_get_field(pbuf, rate1_cw2pr_st_idx, rate1ord_cw2pr_st)
     rate1ord_cw2pr_st(:ncol,top_lev:pver) = rate1cld(:ncol,top_lev:pver)
  end if

  ! Sedimentation velocity for liquid stratus cloud droplet
  call pbuf_get_field(pbuf, wsedl_idx,      wsedl  )
  wsedl(:ncol,top_lev:pver) = vtrmc(:ncol,top_lev:pver)

  ! Assign default size distribution parameters for no-stratiform clouds (convection only)
  ! Also put into physics buffer for possible separate use by radiation
  dcon   = 25.e-6_r8
  mucon  = 5.3_r8
  deicon = 50._r8

  do k = top_lev, pver
     do i = 1, ncol
        ! Convert snow effective diameter to microns
        des(i,k) = des(i,k) * 1.e6_r8
        if ( ast(i,k) < 1.e-4_r8 ) then
           mu(i,k) = mucon
           lambdac(i,k) = (mucon + 1._r8)/dcon
           dei(i,k) = deicon
        end if
     end do
  end do

  ! Microphysical tendencies for use in the macrophysics at the next time step
  call pbuf_get_field(pbuf, cc_t_idx,    CC_t,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_qv_idx,   CC_qv,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_ql_idx,   CC_ql,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_qi_idx,   CC_qi,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_nl_idx,   CC_nl,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_ni_idx,   CC_ni,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cc_qlst_idx, CC_qlst, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

  call pbuf_get_field(pbuf, cmeliq_idx, cmeliq)

  CC_T(:ncol,top_lev:pver)    =  tlat(:ncol,top_lev:pver)/cpair
  CC_qv(:ncol,top_lev:pver)   = qvlat(:ncol,top_lev:pver)
  CC_ql(:ncol,top_lev:pver)   = qcten(:ncol,top_lev:pver)
  CC_qi(:ncol,top_lev:pver)   = qiten(:ncol,top_lev:pver)
  CC_nl(:ncol,top_lev:pver)   = ncten(:ncol,top_lev:pver)
  CC_ni(:ncol,top_lev:pver)   = niten(:ncol,top_lev:pver)
  CC_qlst(:ncol,top_lev:pver) = qcten(:ncol,top_lev:pver)/max(0.01_r8,alst_mic(:ncol,top_lev:pver))

  ! Net micro_mg_cam condensation rate
  call pbuf_get_field(pbuf, qme_idx, qme )

  qme(:ncol,top_lev:pver) = cmeliq(:ncol,top_lev:pver) + cmeiout(:ncol,top_lev:pver)

  ! For precip, accumulate only total precip in prec_pcw and snow_pcw variables.
  ! Other precip output variables are set to 0
  prec_pcw(:ncol) = prect(:ncol)
  snow_pcw(:ncol) = preci(:ncol)
  prec_sed(:ncol) = 0._r8
  snow_sed(:ncol) = 0._r8
  prec_str(:ncol) = prec_pcw(:ncol) + prec_sed(:ncol)
  snow_str(:ncol) = snow_pcw(:ncol) + snow_sed(:ncol)

  icecldf(:ncol,top_lev:pver) = ast(:ncol,top_lev:pver)
  liqcldf(:ncol,top_lev:pver) = ast(:ncol,top_lev:pver)


  ! ------------------------------------------------------------ !
  ! Compute in cloud ice and liquid mixing ratios                !
  ! Note that 'iclwp, iciwp' are used for radiation computation. !
  ! ------------------------------------------------------------ !

  call pbuf_get_field(pbuf, cld_idx,        cld,        start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, concld_idx,     concld,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
  call pbuf_get_field(pbuf, cldfsnow_idx,   cldfsnow,   start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

  call pbuf_get_field(pbuf, iciwpst_idx,    iciwpst  )
  call pbuf_get_field(pbuf, iclwpst_idx,    iclwpst  )
  call pbuf_get_field(pbuf, icswp_idx,      icswp  )

  icimrst_out = 0._r8
  icwmrst_out = 0._r8
  icinc = 0._r8
  icwnc = 0._r8
  iciwpst = 0._r8
  iclwpst = 0._r8
  icswp = 0._r8
  cldfsnow = 0._r8

  do k = top_lev, pver
     do i = 1, ncol
        ! Limits for in-cloud mixing ratios consistent with MG microphysics
        ! in-cloud mixing ratio maximum limit of 0.005 kg/kg
        icimrst(i,k)   = min( state_loc%q(i,k,ixcldice) / max(0.0001_r8,icecldf(i,k)),0.005_r8 )
        icwmrst(i,k)   = min( state_loc%q(i,k,ixcldliq) / max(0.0001_r8,liqcldf(i,k)),0.005_r8 )
        icinc(i,k)     = state_loc%q(i,k,ixnumice) / max(0.0001_r8,icecldf(i,k)) * &
             state_loc%pmid(i,k) / (287.15_r8*state_loc%t(i,k))
        icwnc(i,k)     = state_loc%q(i,k,ixnumliq) / max(0.0001_r8,liqcldf(i,k)) * &
             state_loc%pmid(i,k) / (287.15_r8*state_loc%t(i,k))
        ! Calculate micro_mg_cam cloud water paths in each layer
        ! Note: uses stratiform cloud fraction!
        iciwpst(i,k)   = min(state_loc%q(i,k,ixcldice)/max(0.0001_r8,ast(i,k)),0.005_r8) * state_loc%pdel(i,k) / gravit
        iclwpst(i,k)   = min(state_loc%q(i,k,ixcldliq)/max(0.0001_r8,ast(i,k)),0.005_r8) * state_loc%pdel(i,k) / gravit

        ! ------------------------------ !
        ! Adjust cloud fraction for snow !
        ! ------------------------------ !
        cldfsnow(i,k) = cld(i,k)
        ! If cloud and only ice ( no convective cloud or ice ), then set to 0.
        if( ( cldfsnow(i,k) .gt. 1.e-4_r8 ) .and. &
             ( concld(i,k)   .lt. 1.e-4_r8 ) .and. &
             ( state_loc%q(i,k,ixcldliq) .lt. 1.e-10_r8 ) ) then
           cldfsnow(i,k) = 0._r8
        end if
        ! If no cloud and snow, then set to 0.25
        if( ( cldfsnow(i,k) .lt. 1.e-4_r8 ) .and. ( qsout(i,k) .gt. 1.e-6_r8 ) ) then
           cldfsnow(i,k) = 0.25_r8
        end if
        ! Calculate in-cloud snow water path
        icswp(i,k) = qsout(i,k) / max( 0.0001_r8, cldfsnow(i,k) ) * state_loc%pdel(i,k) / gravit
     end do
  end do

  ! ------------------------------------- !
  ! Precipitation efficiency Calculation  !
  ! ------------------------------------- !

  ! precip efficiency calculation
  call pbuf_get_field(pbuf, acpr_idx,    acprecl)
  call pbuf_get_field(pbuf, acgcme_idx,   acgcme)
  call pbuf_get_field(pbuf, acnum_idx,   acnum)

  !-----------------------------------------------------------------------
  ! Liquid water path

  ! Compute liquid water paths, and column condensation
  tgliqwp(:ncol) = 0._r8
  tgcmeliq(:ncol) = 0._r8
  do k = top_lev, pver
     do i = 1, ncol
        tgliqwp(i)  = tgliqwp(i) + iclwpst(i,k)*cld(i,k)

        if (cmeliq(i,k) > 1.e-12_r8) then
           !convert cmeliq to right units:  kgh2o/kgair/s  *  kgair/m2  / kgh2o/m3  = m/s
           tgcmeliq(i) = tgcmeliq(i) + cmeliq(i,k) * (state_loc%pdel(i,k) / gravit) / rhoh2o
        end if
     end do
  end do

  ! note: 1e-6 kgho2/kgair/s * 1000. pa / (9.81 m/s2) / 1000 kgh2o/m3 = 1e-7 m/s
  ! this is 1ppmv of h2o in 10hpa
  ! alternatively: 0.1 mm/day * 1.e-4 m/mm * 1/86400 day/s = 1.e-9

  !-----------------------------------------------------------------------
  ! precipitation efficiency calculation  (accumulate cme and precip)

  minlwp = 0.01_r8        !minimum lwp threshold (kg/m3)

  ! zero out precip efficiency and total averaged precip
  pe(:ncol)     = 0._r8
  tpr(:ncol)    = 0._r8
  pefrac(:ncol) = 0._r8

  ! accumulate precip and condensation
  do i = 1, ncol

     acgcme(i)  = acgcme(i) + tgcmeliq(i)
     acprecl(i) = acprecl(i) + prec_str(i)
     acnum(i)   = acnum(i) + 1

     ! if LWP is zero, then 'end of cloud': calculate precip efficiency
     if (tgliqwp(i) < minlwp) then
        if (acprecl(i) > 5.e-8_r8) then
           tpr(i) = max(acprecl(i)/acnum(i), 1.e-15_r8)
           if (acgcme(i) > 1.e-10_r8) then
              pe(i) = min(max(acprecl(i)/acgcme(i), 1.e-15_r8), 1.e5_r8)
              pefrac(i) = 1._r8
           end if
        end if

        ! reset counters
        !        if (pe(i) /= 0._r8 .and. (pe(i) < 1.e-8_r8 .or. pe(i) > 1.e3_r8)) then
        !           write (iulog,*) 'PE:ANOMALY  pe, acprecl, acgcme, tpr, acnum ',pe(i),acprecl(i), acgcme(i), tpr(i), acnum(i)
        !        endif

        acprecl(i) = 0._r8
        acgcme(i)  = 0._r8
        acnum(i)   = 0
     end if               ! end LWP zero conditional

     ! if never find any rain....(after 10^3 timesteps...)
     if (acnum(i) > 1000) then
        acnum(i)   = 0
        acprecl(i) = 0._r8
        acgcme(i)  = 0._r8
     end if

  end do

  call outfld( 'PE'    , pe,     pcols, lchnk )
  call outfld( 'PEFRAC', pefrac, pcols, lchnk )
  call outfld( 'APRL'  , tpr,    pcols, lchnk )

  !-----------------------------------------------------------------------
  ! vertical average of non-zero accretion, autoconversion and ratio.
  ! vars: vprco(i),vprao(i),racau(i),cnt

  vprco = 0._r8
  cnt = 0
  do k = top_lev, pver
     vprco(:ncol) = vprco(:ncol) + prco(:ncol,k)
     where (prco(:ncol,k) /= 0._r8) cnt(:ncol) = cnt(:ncol) + 1
  end do

  where (cnt > 0) vprco = vprco/cnt

  vprao = 0._r8
  cnt = 0
  do k = top_lev, pver
     vprao(:ncol) = vprao(:ncol) + prao(:ncol,k)
     where (prao(:ncol,k) /= 0._r8) cnt(:ncol) = cnt(:ncol) + 1
  end do

  where (cnt > 0)
     vprao = vprao/cnt
     racau = vprco/vprao
  elsewhere
     racau = 0._r8
  end where

  racau = min(racau, 1.e10_r8)

  call outfld( 'VPRAO'    , vprao,     pcols, lchnk )
  call outfld( 'VPRCO'    , vprco,     pcols, lchnk )
  call outfld( 'RACAU'    , racau,     pcols, lchnk )

  ! --------------------- !
  ! History Output Fields !
  ! --------------------- !

  ! Column droplet concentration
  cdnumc(:ncol) = sum(state_loc%q(:ncol,top_lev:pver,ixnumliq) * &
       state_loc%pdel(:ncol,top_lev:pver)/gravit, dim=2)

  ! Averaging for new output fields
  efcout      = 0._r8
  efiout      = 0._r8
  ncout       = 0._r8
  niout       = 0._r8
  freql       = 0._r8
  freqi       = 0._r8
  liqcldf_out = 0._r8
  icecldf_out = 0._r8
  icwmrst_out = 0._r8
  icimrst_out = 0._r8

  do k = top_lev, pver
     do i = 1, ncol
        if ( liqcldf(i,k) > 0.01_r8 .and. icwmrst(i,k) > 5.e-5_r8 ) then
           efcout(i,k) = rel(i,k) * liqcldf(i,k)
           ncout(i,k)  = icwnc(i,k) * liqcldf(i,k)
           freql(i,k)  = liqcldf(i,k)
           liqcldf_out(i,k) = liqcldf(i,k)
           icwmrst_out(i,k) = icwmrst(i,k)
        end if
        if ( icecldf(i,k) > 0.01_r8 .and. icimrst(i,k) > 1.e-6_r8 ) then
           efiout(i,k) = rei(i,k) * icecldf(i,k)
           niout(i,k)  = icinc(i,k) * icecldf(i,k)
           freqi(i,k)  = icecldf(i,k)
           icecldf_out(i,k) = icecldf(i,k)
           icimrst_out(i,k) = icimrst(i,k)
        end if
     end do
  end do

  call outfld( 'AREL' , efcout,  pcols, lchnk )
  call outfld( 'AREI' , efiout,  pcols, lchnk )
  call outfld( 'AWNC' , ncout,   pcols, lchnk )
  call outfld( 'AWNI' , niout,   pcols, lchnk )
  call outfld( 'FREQL', freql,   pcols, lchnk )
  call outfld( 'FREQI', freqi,   pcols, lchnk )

  ! Cloud top effective radius and number.
  fcti  = 0._r8
  fctl  = 0._r8
  ctrel = 0._r8
  ctrei = 0._r8
  ctnl  = 0._r8
  ctni  = 0._r8
  do i = 1, ncol
     do k = top_lev, pver
        if ( liqcldf(i,k) > 0.01_r8 .and. icwmrst(i,k) > 1.e-7_r8 ) then
           ctrel(i) = rel(i,k) * liqcldf(i,k)
           ctnl(i)  = icwnc(i,k) * liqcldf(i,k)
           fctl(i)  = liqcldf(i,k)
           exit
        end if
        if ( icecldf(i,k) > 0.01_r8 .and. icimrst(i,k) > 1.e-7_r8 ) then
           ctrei(i) = rei(i,k) * icecldf(i,k)
           ctni(i)  = icinc(i,k) * icecldf(i,k)
           fcti(i)  = icecldf(i,k)
           exit
        end if
     end do
  end do

  call outfld( 'ACTREL'     , ctrel,     pcols, lchnk )
  call outfld( 'ACTREI'     , ctrei,     pcols, lchnk )
  call outfld( 'ACTNL'      , ctnl,      pcols, lchnk )
  call outfld( 'ACTNI'      , ctni,      pcols, lchnk )
  call outfld( 'FCTL'       , fctl,      pcols, lchnk )
  call outfld( 'FCTI'       , fcti,      pcols, lchnk )

  call outfld( 'MPDT'       , tlat,      pcols, lchnk )
  call outfld( 'MPDQ'       , qvlat,     pcols, lchnk )
  call outfld( 'MPDLIQ'     , qcten,     pcols, lchnk )
  call outfld( 'MPDICE'     , qiten,     pcols, lchnk )
  call outfld( 'ICINC'      , icinc,     pcols, lchnk )
  call outfld( 'ICWNC'      , icwnc,     pcols, lchnk )
  call outfld( 'EFFLIQ_IND' , rel_fn,    pcols, lchnk )
  call outfld( 'CDNUMC'     , cdnumc,    pcols, lchnk )

  call outfld('LS_FLXPRC', mgflxprc,    pcols, lchnk )
  call outfld('LS_FLXSNW', mgflxsnw,    pcols, lchnk )
  call outfld('REL', rel,    pcols, lchnk )
  call outfld('REI', rei,    pcols, lchnk )

  ! --------------------------------------------- !
  ! General outfield calls for microphysics       !
  ! --------------------------------------------- !

  call outfld( 'ICIMRST'  , icimrst_out, pcols, lchnk )
  call outfld( 'ICWMRST'  , icwmrst_out, pcols, lchnk )
  call outfld( 'CME'      , qme,         pcols, lchnk )
  call outfld( 'PRODPREC' , prain,       pcols, lchnk )
  call outfld( 'EVAPPREC' , nevapr,      pcols, lchnk )
  call outfld( 'EVAPSNOW' , evapsnow,    pcols, lchnk )
  call outfld( 'QCSEVAP'  , qcsevap,     pcols, lchnk )
  call outfld( 'QISEVAP'  , qisevap,     pcols, lchnk )
  call outfld( 'QVRES'    , qvres,       pcols, lchnk )
  call outfld( 'CMEIOUT'  , cmeiout,     pcols, lchnk )
  call outfld( 'VTRMC'    , vtrmc,       pcols, lchnk )
  call outfld( 'VTRMI'    , vtrmi,       pcols, lchnk )
  call outfld( 'QCSEDTEN' , qcsedten,    pcols, lchnk )
  call outfld( 'QISEDTEN' , qisedten,    pcols, lchnk )
  call outfld( 'PRAO'     , prao,        pcols, lchnk )
  call outfld( 'PRCO'     , prco,        pcols, lchnk )
  call outfld( 'MNUCCCO'  , mnuccco,     pcols, lchnk )
  call outfld( 'MNUCCTO'  , mnuccto,     pcols, lchnk )
  call outfld( 'MNUCCDO'  , mnuccdo,     pcols, lchnk )
  call outfld( 'MNUCCDOhet', mnuccdohet, pcols, lchnk )
  call outfld( 'MSACWIO'  , msacwio,     pcols, lchnk )
  call outfld( 'PSACWSO'  , psacwso,     pcols, lchnk )
  call outfld( 'BERGSO'   , bergso,      pcols, lchnk )
  call outfld( 'BERGO'    , bergo,       pcols, lchnk )
  call outfld( 'MELTO'    , melto,       pcols, lchnk )
  call outfld( 'HOMOO'    , homoo,       pcols, lchnk )
  call outfld( 'QCRESO'   , qcreso,      pcols, lchnk )
  call outfld( 'PRCIO'    , prcio,       pcols, lchnk )
  call outfld( 'PRAIO'    , praio,       pcols, lchnk )
  call outfld( 'QIRESO'   , qireso,      pcols, lchnk )
  call outfld( 'MNUCCRO'  , mnuccro,     pcols, lchnk )
  call outfld( 'PRACSO'   , pracso ,     pcols, lchnk )
  call outfld( 'MELTSDT'  , meltsdt,     pcols, lchnk )
  call outfld( 'FRZRDT'   , frzrdt ,     pcols, lchnk )

  ftem = 0._r8

  ftem(:ncol,top_lev:pver) =  qcreso(:ncol,top_lev:pver)
  call outfld( 'MPDW2V', ftem, pcols, lchnk )

  ftem(:ncol,top_lev:pver) =  melto(:ncol,top_lev:pver) - mnuccco(:ncol,top_lev:pver) - mnuccto(:ncol,top_lev:pver) - &
       bergo(:ncol,top_lev:pver) - homoo  (:ncol,top_lev:pver) - msacwio(:ncol,top_lev:pver)
  call outfld( 'MPDW2I', ftem, pcols, lchnk )

  ftem(:ncol,top_lev:pver) = -prao(:ncol,top_lev:pver) - prco(:ncol,top_lev:pver) - psacwso(:ncol,top_lev:pver) - &
       bergso(:ncol,top_lev:pver)
  call outfld( 'MPDW2P', ftem, pcols, lchnk )

  ftem(:ncol,top_lev:pver) =  cmeiout(:ncol,top_lev:pver) + qireso (:ncol,top_lev:pver)
  call outfld( 'MPDI2V', ftem, pcols, lchnk )

  ftem(:ncol,top_lev:pver) = -melto(:ncol,top_lev:pver) + mnuccco(:ncol,top_lev:pver) + mnuccto(:ncol,top_lev:pver) + &
       bergo(:ncol,top_lev:pver) + homoo  (:ncol,top_lev:pver) + msacwio(:ncol,top_lev:pver)
  call outfld( 'MPDI2W', ftem, pcols, lchnk )

  ftem(:ncol,top_lev:pver) = -prcio(:ncol,top_lev:pver) - praio  (:ncol,top_lev:pver)
  call outfld( 'MPDI2P', ftem, pcols, lchnk )

  ! ptend_loc is deallocated in physics_update above
  call physics_state_dealloc(state_loc)

end subroutine micro_mg_cam_tend

end module micro_mg_cam
