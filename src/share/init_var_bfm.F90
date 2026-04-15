!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! MODULE: init_var_bfm_local
!
! DESCRIPTION
!   Initialize state variables and model setting
!
! COPYING
!
!   Copyright (C) 2023 BFM System Team (bfm_st@cmcc.it)
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation.
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! INCLUDE
#include"cppdefs.h"


MODULE init_var_bfm_local
   USE global_mem,ONLY: RLEN,ZERO
   IMPLICIT NONE

   real(RLEN),parameter :: nc_ratio_default = 0.0126_RLEN    ! Redfield
   real(RLEN),parameter :: pc_ratio_default = 0.7862e-3_RLEN ! Redfield
   real(RLEN),parameter :: lc_ratio_default = 0.03_RLEN      ! standard diatom value
   real(RLEN),parameter :: sc_ratio_default = 0.0145_RLEN    ! Redfield
   real(RLEN),parameter :: hc_ratio_default = ZERO


   character(len=*),parameter :: bfm_init_fname      = 'BFM_General.nml'
   character(len=*),parameter :: bfm_init_fname_ice  = 'BFM_General.nml'
   character(len=*),parameter :: bfm_init_fname_ben  = 'BFM_General.nml'

   CONTAINS

   subroutine init_constituents( c,n,p,l,s,h,ppn,ppp,ppl,pps,pph,nc,pc,lc,sc,hc)
     use global_mem, only: RLEN,ZERO
     IMPLICIT NONE
     real(RLEN),dimension(:),intent(in)             :: c
     real(RLEN),intent(in),optional                 :: nc,pc,lc,sc,hc
     real(RLEN),dimension(:),intent(inout),optional :: n,p,l,s,h
     integer,intent(in),optional                    :: ppn,ppp,ppl,pps,pph
     real(RLEN)                                     :: nc_ratio,pc_ratio,lc_ratio,sc_ratio,hc_ratio
     
         nc_ratio = nc_ratio_default
         if (present(nc)) then
           if (nc>ZERO) nc_ratio = nc
         end if

         pc_ratio = pc_ratio_default
         if (present(pc)) then
           if (pc>ZERO) pc_ratio = pc
         end if

         lc_ratio = lc_ratio_default
         if (present(lc)) then
           if (lc>ZERO) lc_ratio = lc
         end if

         sc_ratio = sc_ratio_default
         if (present(sc)) then
           if (sc>ZERO) sc_ratio = sc
         end if

         hc_ratio = hc_ratio_default
         if (present(hc)) then
           if (hc>ZERO) hc_ratio = hc
         end if

         if (present(n) .and. present(ppn)) then
           if (ppn>0) then
             where (n==ZERO)
               n = nc_ratio*c
             end where
           end if
         end if
         if (present(p) .and. present(ppp)) then
           if (ppp>0) then
             where (p==ZERO)
               p = pc_ratio*c
             end where
           end if
         end if
         if (present(l) .and. present(ppl)) then
           if (ppl>0) then
             where (l==ZERO)
               l = lc_ratio*c
             end where
           end if
         end if
         if (present(s) .and. present(pps)) then
           if (pps>0) then
             where (s==ZERO)
               s = sc_ratio*c
             end where
           end if
         end if
         if (present(h) .and. present(pph)) then
           if (pph>0) then
             where (h==ZERO)
               h = hc_ratio*c
             end where
           end if
         end if
   
         return
   end subroutine init_constituents


!-----------------------------------------------------------------------
!
! ROUTINE: Initialise quotas of living and non-living organic components
!
! DESCRIPTION:
!   Call the initialization of internal ratios for the living
!   and non-living organic components
!
! INTERFACE
   subroutine ini_organic_quotas()
!
! USES
   use mem
   use mem_PelChem
   use mem_PelBac
   use mem_MesoZoo
   use mem_MicroZoo
   use mem_Phyto

   use constants, ONLY: HOURS_PER_DAY
   use mem_Param, ONLY: p_small,          &
                        CalcPelagicFlag,  &
                        CalcBenthicFlag,  &
                        CalcSeaiceFlag,   &
                        CalcPelChemistry, &
                        CalcTransportFlag
#if defined BENTHIC_BIO || defined BENTHIC_FULL
   use mem_BenBac
   use mem_BenOrganisms
#endif

#ifdef INCLUDE_SEAICE
   use mem_SeaiceAlgae
   use mem_SeaiceBac
   use mem_SeaiceZoo
#endif

   IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer              :: i
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   !---------------------------------------------
   ! Initialise other pelagic internal components
   ! with Redfield
   !---------------------------------------------

   do i = 1 , ( iiPelBacteria )
     call init_constituents( c=PelBacteria(i,iiC), &
       n=D3STATE(:,ppPelBacteria(i,iiN)), ppn=ppPelBacteria(i,iiN),  &
       p=D3STATE(:,ppPelBacteria(i,iiP)), ppp=ppPelBacteria(i,iiP),  &
       nc=p_qncPBA(i),  pc=p_qpcPBA(i) )
       qpcPBA(:,i)=p_qpcPBA(i) ; qncPBA(:,i)=p_qncPBA(i) ;
   end do

   do i = 1 , ( iiPhytoPlankton )
     call init_constituents( c=PhytoPlankton(i,iiC), &
       n=D3STATE(:,ppPhytoPlankton(i,iiN)), ppn=ppPhytoPlankton(i,iiN),  &
       p=D3STATE(:,ppPhytoPlankton(i,iiP)), ppp=ppPhytoPlankton(i,iiP),  &
       l=D3STATE(:,ppPhytoPlankton(i,iiL)), ppl=ppPhytoPlankton(i,iiL),  &
       s=D3STATE(:,ppPhytoPlankton(i,iiS)), pps=ppPhytoPlankton(i,iiS),  &
       nc=p_qncPPY(i),  pc=p_qpcPPY(i),  lc=p_qlcPPY(i),  sc=p_qscPPY(i) )
       qpcPPY(:,i)=p_qpcPPY(i) ; qncPPY(:,i)=p_qncPPY(i) ; qscPPY(:,i)=p_qscPPY(i) ; 
       qlcPPY(:,i)=p_qlcPPY(i) ; qccPPY(:,i)=p_qccPPY(i) ;
   end do

   do i = 1 , ( iiMesoZooPlankton )
     call init_constituents( c=MesoZooPlankton(i,iiC), &
       n=D3STATE(:,ppMesoZooPlankton(i,iiN)), ppn=ppMesoZooPlankton(i,iiN),  &
       p=D3STATE(:,ppMesoZooPlankton(i,iiP)), ppp=ppMesoZooPlankton(i,iiP),  &
       nc=p_qncMEZ(i),  pc=p_qpcMEZ(i) )
       qpcMEZ(:,i)=p_qpcMEZ(i) ; qncMEZ(:,i)=p_qncMEZ(i) ;
   end do

   do i = 1 , ( iiMicroZooPlankton )
     call init_constituents( c=MicroZooPlankton(i,iiC), &
       n=D3STATE(:,ppMicroZooPlankton(i,iiN)), ppn=ppMicroZooPlankton(i,iiN),  &
       p=D3STATE(:,ppMicroZooPlankton(i,iiP)), ppp=ppMicroZooPlankton(i,iiP),  &
       nc=p_qncMIZ(i),  pc=p_qpcMIZ(i) )
       qpcMIZ(:,i)=p_qpcMIZ(i) ; qncMIZ(:,i)=p_qncMIZ(i) ;
   end do

   do i = 1 , ( iiPelDetritus )
     call init_constituents( c=PelDetritus(i,iiC), &
       n=D3STATE(:,ppPelDetritus(i,iiN)), ppn=ppPelDetritus(i,iiN),  &
       p=D3STATE(:,ppPelDetritus(i,iiP)), ppp=ppPelDetritus(i,iiP),  &
       s=D3STATE(:,ppPelDetritus(i,iiS)), pps=ppPelDetritus(i,iiS) )
       qpcOMT(:,i)=pc_ratio_default ; qncOMT(:,i)=nc_ratio_default ; qscOMT(:,i)=sc_ratio_default ; 
      
   end do



   !---------------------------------------------
   ! Initialise other benthic internal components
   ! with Redfield
   !---------------------------------------------

   do i = 1 , ( iiBenDetritus )
     call init_constituents( c=BenDetritus(i,iiC), &
       n=D2STATE_BEN(:,ppBenDetritus(i,iiN)), ppn=ppBenDetritus(i,iiN),  &
       p=D2STATE_BEN(:,ppBenDetritus(i,iiP)), ppp=ppBenDetritus(i,iiP),  &
       s=D2STATE_BEN(:,ppBenDetritus(i,iiS)), pps=ppBenDetritus(i,iiS) )
       qpcBOM(:,i)=pc_ratio_default ; qncBOM(:,i)=nc_ratio_default ; qscBOM(:,i)=sc_ratio_default ; 
      
   end do




#ifdef INCLUDE_SEAICE
   !---------------------------------------------
   ! Initialise other seaice internal components
   ! with Redfield
   !---------------------------------------------



#endif
      return

end subroutine ini_organic_quotas

!-----------------------------------------------------------------------
!
! ROUTINE: Update quotas of Living and non-living organic components
!
! DESCRIPTION
!   Update the internal ratios to carbon content for the living
!   and non-living organic components
!
! INTERFACE
   subroutine upd_organic_quotas()
!
! USES
   use mem
   use mem_PelChem
   use mem_PelBac
   use mem_MesoZoo
   use mem_MicroZoo
   use mem_Phyto

   use constants, ONLY: HOURS_PER_DAY
   use mem_Param, ONLY: p_small,          &
                        CalcPelagicFlag,  &
                        CalcBenthicFlag,  &
                        CalcSeaiceFlag,   &
                        CalcPelChemistry, &
                        CalcTransportFlag
#ifdef INCLUDE_SEAICE
   use mem_SeaiceAlgae
   use mem_SeaiceBac
   use mem_SeaiceZoo
#endif

   IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer              :: i
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   !---------------------------------------------
   ! Update quotas for pelagic components
   !---------------------------------------------

   do i = 1 , iiPelBacteria
     if ( ppPelBacteria(i,iiN)  > 0) &
        qncPBA(:,i) = PelBacteria(i,iiN)/( p_small + PelBacteria(i,iiC))
     if ( ppPelBacteria(i,iiP)  > 0) &
        qpcPBA(:,i) = PelBacteria(i,iiP)/( p_small + PelBacteria(i,iiC))
   end do

   do i = 1 , iiPhytoPlankton
     if ( ppPhytoPlankton(i,iiN)  > 0) &
        qncPPY(:,i) = PhytoPlankton(i,iiN)/( p_small + PhytoPlankton(i,iiC))
     if ( ppPhytoPlankton(i,iiP)  > 0) &
        qpcPPY(:,i) = PhytoPlankton(i,iiP)/( p_small + PhytoPlankton(i,iiC))
     if ( ppPhytoPlankton(i,iiL)  > 0) &
        qlcPPY(:,i) = PhytoPlankton(i,iiL)/( p_small + PhytoPlankton(i,iiC))
     if ( ppPhytoPlankton(i,iiS)  > 0) &
        qscPPY(:,i) = PhytoPlankton(i,iiS)/( p_small + PhytoPlankton(i,iiC))
   end do

   do i = 1 , iiMesoZooPlankton
     if ( ppMesoZooPlankton(i,iiN)  > 0) &
        qncMEZ(:,i) = MesoZooPlankton(i,iiN)/( p_small + MesoZooPlankton(i,iiC))
     if ( ppMesoZooPlankton(i,iiP)  > 0) &
        qpcMEZ(:,i) = MesoZooPlankton(i,iiP)/( p_small + MesoZooPlankton(i,iiC))
   end do

   do i = 1 , iiMicroZooPlankton
     if ( ppMicroZooPlankton(i,iiN)  > 0) &
        qncMIZ(:,i) = MicroZooPlankton(i,iiN)/( p_small + MicroZooPlankton(i,iiC))
     if ( ppMicroZooPlankton(i,iiP)  > 0) &
        qpcMIZ(:,i) = MicroZooPlankton(i,iiP)/( p_small + MicroZooPlankton(i,iiC))
   end do

   do i = 1 , iiPelDetritus
     if ( ppPelDetritus(i,iiN)  > 0) &
        qncOMT(:,i) = PelDetritus(i,iiN)/( p_small + PelDetritus(i,iiC))
     if ( ppPelDetritus(i,iiP)  > 0) &
        qpcOMT(:,i) = PelDetritus(i,iiP)/( p_small + PelDetritus(i,iiC))
     if ( ppPelDetritus(i,iiS)  > 0) &
        qscOMT(:,i) = PelDetritus(i,iiS)/( p_small + PelDetritus(i,iiC))
   end do



   !---------------------------------------------
   ! Update quotas for benthic components
   !---------------------------------------------

   do i = 1 , iiBenDetritus
     if ( ppBenDetritus(i,iiN)  > 0) &
        qncBOM(:,i) = BenDetritus(i,iiN)/( p_small + BenDetritus(i,iiC))
     if ( ppBenDetritus(i,iiP)  > 0) &
        qpcBOM(:,i) = BenDetritus(i,iiP)/( p_small + BenDetritus(i,iiC))
     if ( ppBenDetritus(i,iiS)  > 0) &
        qscBOM(:,i) = BenDetritus(i,iiS)/( p_small + BenDetritus(i,iiC))
   end do




#ifdef INCLUDE_SEAICE
   !---------------------------------------------
   ! Update quotas for seaice components
   !---------------------------------------------



#endif
      return
end subroutine upd_organic_quotas


END MODULE init_var_bfm_local


!-----------------------------------------------------------------------
!
! ROUTINE: Initialise BFM variables
!
! DESCRIPTION:
!   Allocate BFM variables and give initial values of
!   parameters and state variables
!
! INTERFACE
   subroutine init_var_bfm(setup)
!
! USES
#ifndef NOT_STANDALONE
   use api_bfm
   use global_mem
#else
   use api_bfm, ONLY: GetLun
   use global_mem, ONLY: NMLUNIT
#endif
   use mem
   use mem_PelChem
   use mem_PelBac
   use mem_MesoZoo
   use mem_MicroZoo
   use mem_Phyto
   use init_var_bfm_local

   use constants, ONLY: HOURS_PER_DAY
   use mem_Param, ONLY: p_small,          &
                        CalcPelagicFlag,  &
                        CalcBenthicFlag,  &
                        CalcSeaiceFlag,   &
                        CalcPelChemistry, &
                        CalcTransportFlag

   use mem_Param, ONLY: AssignPelBenFluxesInBFMFlag
   use string_functions, ONLY: getseq_number,empty
   
#if defined BENTHIC_FULL && defined INCLUDE_BENPROFILES
   use mem_Param, ONLY: p_d_tot, p_sedlevels,p_sedsigma
#endif 
#ifdef INCLUDE_SEAICE
   use mem_SeaiceAlgae
   use mem_SeaiceBac
   use mem_SeaiceZoo
#endif 

   IMPLICIT NONE

  ! INPUT
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer,          intent(in)        :: setup

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   integer              :: icontrol,i,j,n,Flun
   integer,parameter    :: NSAVE=300  ! Maximum no variables which can be saved
   character(len=64),dimension(NSAVE):: var_save
   character(len=64),dimension(NSAVE):: ave_save
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Definition of Initial Pelagic (D3) state variables
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

   real(RLEN) :: O2o0, N1p0, N3n0, N4n0, O4n0, N5s0, N6r0, B1c0, B1n0, B1p0, &
    P1c0, P1n0, P1p0, P1l0, P1s0, P2c0, P2n0, P2p0, P2l0, P3c0, P3n0, P3p0, &
    P3l0, P4c0, P4n0, P4p0, P4l0, Z3c0, Z3n0, Z3p0, Z4c0, Z4n0, Z4p0, Z5c0, &
    Z5n0, Z5p0, Z6c0, Z6n0, Z6p0, R1c0, R1n0, R1p0, R2c0, R3c0, R6c0, R6n0, &
    R6p0, R6s0, O3c0, O3h0, O5c0


   namelist /bfm_init_nml/ O2o0, N1p0, N3n0, N4n0, O4n0, N5s0, N6r0, B1c0, &
    B1n0, B1p0, P1c0, P1n0, P1p0, P1l0, P1s0, P2c0, P2n0, P2p0, P2l0, P3c0, &
    P3n0, P3p0, P3l0, P4c0, P4n0, P4p0, P4l0, Z3c0, Z3n0, Z3p0, Z4c0, Z4n0, &
    Z4p0, Z5c0, Z5n0, Z5p0, Z6c0, Z6n0, Z6p0, R1c0, R1n0, R1p0, R2c0, R3c0, &
    R6c0, R6n0, R6p0, R6s0, O3c0, O3h0, O5c0


   namelist /bfm_init_nml/ surface_flux_method,       &
                           bottom_flux_method,        &
                           n_surface_fluxes, InitVar

   namelist /bfm_save_nml/ var_save, ave_save

   real(RLEN) :: Q1c0, Q1n0, Q1p0, Q11c0, Q11n0, Q11p0, Q6c0, Q6n0, Q6p0, &
    Q6s0, Q16c0, Q16n0, Q16p0, Q16s0


   namelist /bfm_init_nml_ben/ Q1c0, Q1n0, Q1p0, Q11c0, Q11n0, Q11p0, Q6c0, &
    Q6n0, Q6p0, Q6s0, Q16c0, Q16n0, Q16p0, Q16s0


#ifdef INCLUDE_SEAICE



#endif

   LEVEL2 'init_var_bfm'
   !---------------------------------------------
   ! Give zero initial values
   ! Overwritten by namelist parameters
   !---------------------------------------------
   surface_flux_method = -1
   bottom_flux_method = 0
   n_surface_fluxes = 1

   !---------------------------------------------
   ! Pelagic variables
   !---------------------------------------------

   O2o0 = _ZERO_
   N1p0 = _ZERO_
   N3n0 = _ZERO_
   N4n0 = _ZERO_
   O4n0 = _ZERO_
   N5s0 = _ZERO_
   N6r0 = _ZERO_
   B1c0 = _ZERO_
   B1n0 = _ZERO_
   B1p0 = _ZERO_
   P1c0 = _ZERO_
   P1n0 = _ZERO_
   P1p0 = _ZERO_
   P1l0 = _ZERO_
   P1s0 = _ZERO_
   P2c0 = _ZERO_
   P2n0 = _ZERO_
   P2p0 = _ZERO_
   P2l0 = _ZERO_
   P3c0 = _ZERO_
   P3n0 = _ZERO_
   P3p0 = _ZERO_
   P3l0 = _ZERO_
   P4c0 = _ZERO_
   P4n0 = _ZERO_
   P4p0 = _ZERO_
   P4l0 = _ZERO_
   Z3c0 = _ZERO_
   Z3n0 = _ZERO_
   Z3p0 = _ZERO_
   Z4c0 = _ZERO_
   Z4n0 = _ZERO_
   Z4p0 = _ZERO_
   Z5c0 = _ZERO_
   Z5n0 = _ZERO_
   Z5p0 = _ZERO_
   Z6c0 = _ZERO_
   Z6n0 = _ZERO_
   Z6p0 = _ZERO_
   R1c0 = _ZERO_
   R1n0 = _ZERO_
   R1p0 = _ZERO_
   R2c0 = _ZERO_
   R3c0 = _ZERO_
   R6c0 = _ZERO_
   R6n0 = _ZERO_
   R6p0 = _ZERO_
   R6s0 = _ZERO_
   O3c0 = _ZERO_
   O3h0 = _ZERO_
   O5c0 = _ZERO_


   !---------------------------------------------
   ! Benthic variables
   !---------------------------------------------

   Q1c0 = _ZERO_
   Q1n0 = _ZERO_
   Q1p0 = _ZERO_
   Q11c0 = _ZERO_
   Q11n0 = _ZERO_
   Q11p0 = _ZERO_
   Q6c0 = _ZERO_
   Q6n0 = _ZERO_
   Q6p0 = _ZERO_
   Q6s0 = _ZERO_
   Q16c0 = _ZERO_
   Q16n0 = _ZERO_
   Q16p0 = _ZERO_
   Q16s0 = _ZERO_



#ifdef INCLUDE_SEAICE
   !---------------------------------------------
   ! Seaice variables
   !---------------------------------------------



#endif

   !---------------------------------------------
   ! Initialize the structured array that 
   ! defines if a variable is initialized with 
   ! data. The namelist values override the
   ! assignment
   !---------------------------------------------
   InitVar = InputInfo(0, ZERO, "dummy.nc","dummy",ZERO,ZERO,ZERO,ZERO,.FALSE.,.FALSE.,.FALSE.,.FALSE.)

   !---------------------------------------------
   ! Open and read the namelist
   !---------------------------------------------
   icontrol=0
   NMLUNIT=GetLun()
   open(NMLUNIT,file=bfm_init_fname,action='read',status='old',err=100)
   var_save=""
   ave_save=""
   var_ave=.false.
   read(NMLUNIT,nml=bfm_save_nml,err=101)
   close(NMLUNIT)
   icontrol=1
100 if ( icontrol == 0 ) then
     LEVEL3 'I could not open ',trim(bfm_init_fname)
     LEVEL3 'The initial values of the BFM variables are set to ZERO'
     LEVEL3 'If thats not what you want you have to supply ',trim(bfm_init_fname)
     icontrol=1
  end if
101 if ( icontrol == 0 ) then
     FATAL 'I could not read bfm_save_nml'
     stop 'init_var_bfm'
   end if

   icontrol=0
   open(NMLUNIT,file=bfm_init_fname,action='read',status='old',err=102)
   read(NMLUNIT,nml=bfm_init_nml,err=103)
   close(NMLUNIT)
   icontrol=1
102 if ( icontrol == 0 ) then
     LEVEL3 'I could not open ',trim(bfm_init_fname)
     LEVEL3 'The initial values of the BFM variables are set to ZERO'
     LEVEL3 'If thats not what you want you have to supply ',trim(bfm_init_fname)
     icontrol=1
   end if
103 if ( icontrol == 0 ) then
     FATAL 'Misspelled variable in bfm_init_nml'
     stop 'init_var_bfm'
   end if


   icontrol=0
   open(NMLUNIT,file=bfm_init_fname_ben,action='read',status='old',err=104)
   read(NMLUNIT,nml=bfm_init_nml_ben,err=105)
   close(NMLUNIT)
   icontrol=1
104 if ( icontrol == 0 ) then
     LEVEL3 'I could not open ',trim(bfm_init_fname_ben)
     LEVEL3 'The initial values of the BFM variables are set to ZERO'
     LEVEL3 'If thats not what you want you have to supply ',trim(bfm_init_fname_ben)
     icontrol=1
   end if
105 if ( icontrol == 0 ) then
     FATAL 'Misspelled variable in bfm_init_nml_ben'
     stop 'init_var_bfm'
   end if


#ifdef INCLUDE_SEAICE
   icontrol=0
   open(NMLUNIT,file=bfm_init_fname_ice,action='read',status='old',err=106)
   read(NMLUNIT,nml=bfm_init_nml_ice,err=107)
   close(NMLUNIT)
   icontrol=1
106 if ( icontrol == 0 ) then
     LEVEL3 'I could not open ',trim(bfm_init_fname_ice)
     LEVEL3 'The initial values of the BFM variables are set to ZERO'
     LEVEL3 'If thats not what you want you have to supply ',trim(bfm_init_fname_ice)
     icontrol=1
   end if
107 if ( icontrol == 0 ) then
     FATAL 'Misspelled variable in bfm_init_nml_ice'
     stop 'init_var_bfm'
   end if

#endif

   !---------------------------------------------
   ! Check variable to be saved and
   ! set the corresponding flag value in var_ids
   !---------------------------------------------
   do i=1,NSAVE
      if (.NOT.empty(var_save(i))) then
            j=getseq_number(var_save(i),var_names,stEnd,.TRUE.)
            if ( j > 0 ) var_ids(j)=-1
      end if
      if ( .NOT.empty(var_save(i)) .AND. j==0 ) then
            STDERR 'Warning: variable ',trim(var_save(i)),' does not exist!'
      end if
   end do
   do i=1,NSAVE
      if (.NOT.empty(ave_save(i))) then
         j=getseq_number(ave_save(i),var_names,stEnd,.TRUE.)
         if ( .NOT.empty(ave_save(i)) .AND. j==0 ) then
            STDERR 'Warning: variable ',trim(ave_save(i)),' does not exist!'
         else if ( var_ids(j) <0 ) then
            STDERR 'Warning: Variable ',trim(ave_save(i)), &
               ' is already selected for output in var_save'
         else if ( j > 0 ) then
            var_ids(j) = -1
            var_ave(j) = .true.
            ave_ctl = .true.
         end if
      end if
   end do

   !---------------------------------------------
   ! Initialize BFM parameters
   !---------------------------------------------
   call Initialize

   !---------------------------------------------
   ! Initially set the number of sun hours
   ! equal to the number of hours in a day.
   !---------------------------------------------
   SUNQ = HOURS_PER_DAY

   !---------------------------------------------
   ! Initialise pelagic state variables
   ! also if using a benthic-only setup
   ! (for boundary conditions)
   !---------------------------------------------

     O2o = O2o0
     N1p = N1p0
     N3n = N3n0
     N4n = N4n0
     O4n = O4n0
     N5s = N5s0
     N6r = N6r0
     B1c = B1c0
     B1n = B1n0
     B1p = B1p0
     P1c = P1c0
     P1n = P1n0
     P1p = P1p0
     P1l = P1l0
     P1s = P1s0
     P2c = P2c0
     P2n = P2n0
     P2p = P2p0
     P2l = P2l0
     P3c = P3c0
     P3n = P3n0
     P3p = P3p0
     P3l = P3l0
     P4c = P4c0
     P4n = P4n0
     P4p = P4p0
     P4l = P4l0
     Z3c = Z3c0
     Z3n = Z3n0
     Z3p = Z3p0
     Z4c = Z4c0
     Z4n = Z4n0
     Z4p = Z4p0
     Z5c = Z5c0
     Z5n = Z5n0
     Z5p = Z5p0
     Z6c = Z6c0
     Z6n = Z6n0
     Z6p = Z6p0
     R1c = R1c0
     R1n = R1n0
     R1p = R1p0
     R2c = R2c0
     R3c = R3c0
     R6c = R6c0
     R6n = R6n0
     R6p = R6p0
     R6s = R6s0
     O3c = O3c0
     O3h = O3h0
     O5c = O5c0


   !---------------------------------------------
   ! Initialise benthic state variables
   ! also if using a ben-only setup
   ! (for boundary conditions)
   !---------------------------------------------

     Q1c = Q1c0
     Q1n = Q1n0
     Q1p = Q1p0
     Q11c = Q11c0
     Q11n = Q11n0
     Q11p = Q11p0
     Q6c = Q6c0
     Q6n = Q6n0
     Q6p = Q6p0
     Q6s = Q6s0
     Q16c = Q16c0
     Q16n = Q16n0
     Q16p = Q16p0
     Q16s = Q16s0


#ifdef INCLUDE_SEAICE
   !---------------------------------------------
   ! Initialise seaice state variables
   ! also if using a seaice-only setup
   ! (for boundary conditions)
   !---------------------------------------------



#endif

   !---------------------------------------------
   ! Check setup settings
   ! and finalize initialization
   !---------------------------------------------
   select case (setup)
      case (0)
         LEVEL2 "Fully coupled system, Pelagic, Benthic, Seaice"
      case (1) ! Pelagic only
         LEVEL2 "Pelagic-only setup (bio_setup=1), Switching off other systems & related fluxes"
         CalcBenthicFlag = .FALSE.
         CalcSeaiceFlag  = .FALSE.
         CalcBenDetritus = .FALSE.


      case (2) ! Benthic only
         LEVEL2 "Benthic-only setup (bio_setup=2), Switching off other systems & related fluxes"
         CalcPelagicFlag = .FALSE.
         CalcSeaiceFlag  = .FALSE.
         CalcPelBacteria = .FALSE.
         CalcPhytoPlankton = .FALSE.
         CalcMesoZooPlankton = .FALSE.
         CalcMicroZooPlankton = .FALSE.
         CalcPelDetritus = .FALSE.
         CalcInorganic = .FALSE.


      case (3) ! Pelagic-Benthic coupling
         LEVEL2 "Pelagic-Benthic setup (bio_setup=3), Switching off seaice system & related fluxes"
         CalcSeaiceFlag  = .FALSE.

      case (4) ! SeaIce only
         LEVEL2 "Seaice-only setup (bio_setup=4), Switching off other systems & related fluxes"
         CalcPelagicFlag = .FALSE.
         CalcBenthicFlag = .FALSE.
         CalcPelBacteria = .FALSE.
         CalcPhytoPlankton = .FALSE.
         CalcMesoZooPlankton = .FALSE.
         CalcMicroZooPlankton = .FALSE.
         CalcPelDetritus = .FALSE.
         CalcInorganic = .FALSE.

         CalcBenDetritus = .FALSE.

      case (5) ! Pelagic-SeaIce coupling
         LEVEL2 "Pelagic-Seaice setup (bio_setup=5), Switching off benthic system & related fluxes"
         CalcBenthicFlag = .FALSE.
         CalcBenDetritus = .FALSE.

   end select

   !---------------------------------------------
   ! Check benthic model
   !---------------------------------------------
   if (CalcBenthicFlag) then
#if ! defined BENTHIC_BIO && ! defined BENTHIC_FULL
        LEVEL3 "Benthic model is: simple nutrient return"
#endif
#if defined BENTHIC_BIO
        LEVEL3 "Benthic model is: benthos + intermediate nutrient return"
#endif
#if defined BENTHIC_FULL
        LEVEL3 "Benthic model is: benthos + Ruardij & Van Raaphorst"
#endif
   else
        LEVEL3 "Benthic model is: not used"
   endif

   !---------------------------------------------
   ! Check for transport flag
   !---------------------------------------------
#ifdef BFM_STANDALONE
   D3STATETYPE(:) = NOTRANSPORT
#else
   if (.NOT.CalcTransportFlag) D3STATETYPE(:) = NOTRANSPORT
#endif

   !----------------------------------------------------
   ! Zeroing of the switched off pelagic state variables
   !----------------------------------------------------

   do j = 1, iiPelBacteria
     if (.NOT.CalcPelBacteria(j)) then
       do i = 1,iiLastElement
         if ( ppPelBacteria(j,i) /= 0 ) then 
           D3STATE(:,ppPelBacteria(j,i)) = p_small
           D3STATETYPE(ppPelBacteria(j,i)) = OFF
         end if
       end do
     end if
   end do

   do j = 1, iiPhytoPlankton
     if (.NOT.CalcPhytoPlankton(j)) then
       do i = 1,iiLastElement
         if ( ppPhytoPlankton(j,i) /= 0 ) then 
           D3STATE(:,ppPhytoPlankton(j,i)) = p_small
           D3STATETYPE(ppPhytoPlankton(j,i)) = OFF
         end if
       end do
     end if
   end do

   do j = 1, iiMesoZooPlankton
     if (.NOT.CalcMesoZooPlankton(j)) then
       do i = 1,iiLastElement
         if ( ppMesoZooPlankton(j,i) /= 0 ) then 
           D3STATE(:,ppMesoZooPlankton(j,i)) = p_small
           D3STATETYPE(ppMesoZooPlankton(j,i)) = OFF
         end if
       end do
     end if
   end do

   do j = 1, iiMicroZooPlankton
     if (.NOT.CalcMicroZooPlankton(j)) then
       do i = 1,iiLastElement
         if ( ppMicroZooPlankton(j,i) /= 0 ) then 
           D3STATE(:,ppMicroZooPlankton(j,i)) = p_small
           D3STATETYPE(ppMicroZooPlankton(j,i)) = OFF
         end if
       end do
     end if
   end do

   do j = 1, iiPelDetritus
     if (.NOT.CalcPelDetritus(j)) then
       do i = 1,iiLastElement
         if ( ppPelDetritus(j,i) /= 0 ) then 
           D3STATE(:,ppPelDetritus(j,i)) = p_small
           D3STATETYPE(ppPelDetritus(j,i)) = OFF
         end if
       end do
     end if
   end do

   do j = 1, iiInorganic
     if (.NOT.CalcInorganic(j)) then
       do i = 1,iiLastElement
         if ( ppInorganic(j,i) /= 0 ) then 
           D3STATE(:,ppInorganic(j,i)) = p_small
           D3STATETYPE(ppInorganic(j,i)) = OFF
         end if
       end do
     end if
   end do



   !----------------------------------------------------
   ! Zeroing of the switched off benthic state variables
   !----------------------------------------------------

   do j = 1, iiBenDetritus
     if (.NOT.CalcBenDetritus(j)) then
       do i = 1,iiLastElement
         if ( ppBenDetritus(j,i) /= 0 ) then 
           D2STATE_BEN(:,ppBenDetritus(j,i)) = p_small
           D2STATETYPE_BEN(ppBenDetritus(j,i)) = OFF
         end if
       end do
     end if
   end do





#ifdef INCLUDE_SEAICE
   !----------------------------------------------------
   ! Zeroing of the switched off seaice state variables
   !----------------------------------------------------



#endif


   !---------------------------------------------
   ! Write defined variables to stdout
   !---------------------------------------------
#ifdef BFM_PARALLEL
   Flun = LOGUNIT
#else
   Flun = stderr
#endif

   if (setup == 0 .OR. setup == 1 .OR. setup == 3 .OR. setup == 5 ) then
      LEVEL3 'Pelagic variables:'
      write(Flun,155) 'ID','Var','Unit','Long Name','Flag'
      do n=stPelStateS,stPelStateE
        write(Flun,156) n,trim(var_names(n)),trim(var_units(n)) &
          ,trim(var_long(n)),D3STATETYPE(n-stPelStateS+1)
      end do

   endif

   if (setup == 0 .OR. ( setup >= 2 .AND. setup <= 3 ) ) then
      LEVEL3 'Benthic variables:'
      write(Flun,155) 'ID','Var','Unit','Long Name','Flag'
      do n=stBenStateS,stBenStateE
        write(Flun,156) n,trim(var_names(n)),trim(var_units(n)) &
          ,trim(var_long(n)),D2STATETYPE_BEN(n-stBenStateS+1)
      end do

#ifdef INCLUDE_BENPROFILES
      !---------------------------------------------
      ! initialize the vertical grid for benthic 
      ! nutrient profiles
      !---------------------------------------------
      LEVEL2 'Initialize the vertical grid for benthic profile diagnostics'
      LEVEL3 'Vertical sediment grid forced equal to model grid'
      p_sedlevels = NO_BOXES_Z_BEN
      call calc_sigma_depth(p_sedlevels,p_sedsigma,p_d_tot,seddepth)
      do n=1,p_sedlevels
         LEVEL3 n,seddepth(n)
      end do
#endif
   endif

#ifdef INCLUDE_SEAICE
   if (setup == 0 .OR. ( setup >= 4 .AND. setup <= 5 ) ) then
      LEVEL3 'Seaice variables:'
      write(Flun,155) 'ID','Var','Unit','Long Name','Flag'
      do n=stIceStateS,stIceStateE
        write(Flun,156) n,trim(var_names(n)),trim(var_units(n)) &
          ,trim(var_long(n)),D2STATETYPE_ICE(n-stIceStateS+1)
      end do

   endif
#endif

   return

155 FORMAT(10x,a4,1x,a5,1x,a12,1x,a40,1x,a10)
156 FORMAT(10x,i4,1x,a5,1x,a12,1x,a40,1x,i6)
   end subroutine init_var_bfm

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

