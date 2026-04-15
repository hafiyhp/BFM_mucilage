!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! MODULE: ModuleMem
!
! DESCRIPTION
!   Definition of Global Shared Memory
!   This module contains all the structural definitions of the BFM
!   and sets up the memory layout.
!   It is automatically generated from the prototype file 
!   BFM/proto/ModuleMem.proto by including the information from 
!   BFM/General/GlobalDefsBFM.model
!   Do not directly edit this code because changes will be lost at
!   any new compilation.
!
!! IMPORTANT NOTE:
!! Do not change the lines starting with two comment characters "!" 
!! These lines are used by the parser to generate the final module file
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
#include "DEBUG.h"
!
      module mem
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Modules can optionally use (import) other modules
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        USE BFM_ERROR_MSG, ONLY: BFM_ERROR
        use global_mem
        use constants
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Implicit typing is never allowed
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        implicit none
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Default all is private
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        private
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! State variables Info (pel)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

!    3d name                                                  description           unit
! ---------- ------------------------------------------------------------ ---------------
!        O2o                                                       Oxygen      mmol O2/m3
!        N1p                                                    Phosphate       mmol P/m3
!        N3n                                                      Nitrate       mmol N/m3
!        N4n                                                     Ammonium       mmol N/m3
!        O4n                                                 NitrogenSink       mmol N/m3
!        N5s                                                     Silicate      mmol Si/m3
!        N6r                                        Reduction Equivalents     mmol S--/m3
!        B1c                               Aerobic and Anaerobic Bacteria         mg C/m3
!        B1n                               Aerobic and Anaerobic Bacteria       mmol N/m3
!        B1p                               Aerobic and Anaerobic Bacteria       mmol P/m3
!        P1c                                                      Diatoms         mg C/m3
!        P1n                                                      Diatoms       mmol N/m3
!        P1p                                                      Diatoms       mmol P/m3
!        P1l                                                      Diatoms       mg Chl/m3
!        P1s                                                      Diatoms      mmol Si/m3
!        P2c                                                  Flagellates         mg C/m3
!        P2n                                                  Flagellates       mmol N/m3
!        P2p                                                  Flagellates       mmol P/m3
!        P2l                                                  Flagellates       mg Chl/m3
!        P3c                                            PicoPhytoplankton         mg C/m3
!        P3n                                            PicoPhytoplankton       mmol N/m3
!        P3p                                            PicoPhytoplankton       mmol P/m3
!        P3l                                            PicoPhytoplankton       mg Chl/m3
!        P4c                                          Large Phytoplankton         mg C/m3
!        P4n                                          Large Phytoplankton       mmol N/m3
!        P4p                                          Large Phytoplankton       mmol P/m3
!        P4l                                          Large Phytoplankton       mg Chl/m3
!        Z3c                                  Carnivorous Mesozooplankton         mg C/m3
!        Z3n                                  Carnivorous Mesozooplankton       mmol N/m3
!        Z3p                                  Carnivorous Mesozooplankton       mmol P/m3
!        Z4c                                   Omnivorous Mesozooplankton         mg C/m3
!        Z4n                                   Omnivorous Mesozooplankton       mmol N/m3
!        Z4p                                   Omnivorous Mesozooplankton       mmol P/m3
!        Z5c                                             Microzooplankton         mg C/m3
!        Z5n                                             Microzooplankton       mmol N/m3
!        Z5p                                             Microzooplankton       mmol P/m3
!        Z6c                         Heterotrophic Nanoflagellates (HNAN)         mg C/m3
!        Z6n                         Heterotrophic Nanoflagellates (HNAN)       mmol N/m3
!        Z6p                         Heterotrophic Nanoflagellates (HNAN)       mmol P/m3
!        R1c                              Labile Dissolved Organic Matter         mg C/m3
!        R1n                              Labile Dissolved Organic Matter       mmol N/m3
!        R1p                              Labile Dissolved Organic Matter       mmol P/m3
!        R2c                         Semi-labile Dissolved Organic Carbon         mg C/m3
!        R3c                     Semi-refractory Dissolved Organic Carbon         mg C/m3
!        R6c                                   Particulate Organic Matter         mg C/m3
!        R6n                                   Particulate Organic Matter       mmol N/m3
!        R6p                                   Particulate Organic Matter       mmol P/m3
!        R6s                                   Particulate Organic Matter      mmol Si/m3
!        O3c                                   Dissolved Inorganic Carbon         mg C/m3
!        O3h                                   Dissolved Inorganic Carbon      mmol eq/m3
!        O5c                                                      Calcite         mg C/m3


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! State variables Info (ice)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Definition of arrays which will hold all state variables and other
        ! global variables  used for exchange between submodels and/or output
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        
        real(RLEN),public,pointer,dimension(:,:) :: D3STATE

        real(RLEN),public,pointer,dimension(:,:)   :: D3SOURCE

        integer,public,pointer,dimension(:) :: D3STATETYPE

        real(RLEN),public,pointer,dimension(:,:) ::  PELSURFACE

        real(RLEN),public,pointer,dimension(:,:) ::  PELBOTTOM


        real(RLEN),public,pointer,dimension(:,:) :: D3DIAGNOS
        real(RLEN),public,pointer,dimension(:,:) :: D2DIAGNOS

        type ppFlux
          integer, allocatable, dimension(:) :: p
          integer, allocatable, dimension(:) :: dir
        end type ppFlux
        type(ppFlux), public,dimension(:,:),allocatable :: D3FLUX_MATRIX
        real(RLEN)  , public,dimension(:,:),allocatable :: D3FLUX_FUNC


#ifdef INCLUDE_SEAICE
        real(RLEN),public,pointer,dimension(:,:) :: D2STATE_ICE

        real(RLEN),public,pointer,dimension(:,:)   :: D2SOURCE_ICE

        integer,public,pointer,dimension(:) :: D2STATETYPE_ICE

        real(RLEN),public,pointer,dimension(:,:) :: D2DIAGNOS_ICE

        type(ppFlux), public,dimension(:,:),allocatable :: D2FLUX_MATRIX_ICE
        real(RLEN)  , public,dimension(:,:),allocatable :: D2FLUX_FUNC_ICE
#endif



        real(RLEN),public,pointer,dimension(:,:) :: D2STATE_BEN

        real(RLEN),public,pointer,dimension(:,:)   :: D2SOURCE_BEN

        integer,public,pointer,dimension(:) :: D2STATETYPE_BEN

        real(RLEN),public,pointer,dimension(:,:) :: D2DIAGNOS_BEN

        type(ppFlux), public,dimension(:,:),allocatable :: D2FLUX_MATRIX_BEN
        real(RLEN)  , public,dimension(:,:),allocatable :: D2FLUX_FUNC_BEN



        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! GLOBAL system CONSTANTS
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,parameter,public  ::iiPel   = 0
        integer,parameter,public  ::iiIce   = 700
        integer,parameter,public  ::iiBen   = 1000
        integer,parameter,public  ::iiReset = -1000

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! GLOBAL system CONSTANTS (pel)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,parameter,public  ::NO_D3_BOX_STATES   = 51

        integer,parameter,public  ::NO_D3_BOX_DIAGNOSS = 99

        integer,parameter,public  ::NO_D2_BOX_DIAGNOSS = 116

        integer,parameter,public  ::NO_D3_BOX_FLUX     = 17

        integer,public  ::NO_BOXES
        integer,public  ::NO_BOXES_X
        integer,public  ::NO_BOXES_Y
        integer,public  ::NO_BOXES_Z
        integer,public  ::NO_STATES
        integer,public  ::NO_BOXES_XY

#ifdef INCLUDE_SEAICE
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! GLOBAL system CONSTANTS (ice)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,parameter,public  ::NO_D2_BOX_STATES_ICE   = 0

        integer,parameter,public  ::NO_D2_BOX_DIAGNOSS_ICE = 0

        integer,parameter,public  ::NO_D2_BOX_FLUX_ICE     = 0


        integer,public  ::NO_BOXES_ICE
        integer,public  ::NO_BOXES_Z_ICE
        integer,public  ::NO_STATES_ICE
#endif

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! GLOBAL system CONSTANTS (ben)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,parameter,public  ::NO_D2_BOX_STATES_BEN   = 14

        integer,parameter,public  ::NO_D2_BOX_DIAGNOSS_BEN = 16

        integer,parameter,public  ::NO_D2_BOX_FLUX_BEN     = 1


        integer,public  ::NO_BOXES_BEN
        integer,public  ::NO_BOXES_Z_BEN
        integer,public  ::NO_STATES_BEN

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !! GLOBAL definition of Pelagic (D3/D2) state variables
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,parameter,public                :: ppNULL=0

        integer,parameter,public                :: ppO2o=1, ppN1p=2, ppN3n=3, &
    ppN4n=4, ppO4n=5, ppN5s=6, ppN6r=7, ppB1c=8, ppB1n=9, ppB1p=10, ppP1c=11, &
    ppP1n=12, ppP1p=13, ppP1l=14, ppP1s=15, ppP2c=16, ppP2n=17, ppP2p=18, &
    ppP2l=19, ppP2s=0, ppP3c=20, ppP3n=21, ppP3p=22, ppP3l=23, ppP3s=0, &
    ppP4c=24, ppP4n=25, ppP4p=26, ppP4l=27, ppP4s=0, ppZ3c=28, ppZ3n=29, &
    ppZ3p=30, ppZ4c=31, ppZ4n=32, ppZ4p=33, ppZ5c=34, ppZ5n=35, ppZ5p=36, &
    ppZ6c=37, ppZ6n=38, ppZ6p=39, ppR1c=40, ppR1n=41, ppR1p=42, ppR1s=0, &
    ppR2c=43, ppR2n=0, ppR2p=0, ppR2s=0, ppR3c=44, ppR3n=0, ppR3p=0, ppR3s=0, &
    ppR6c=45, ppR6n=46, ppR6p=47, ppR6s=48, ppO3c=49, ppO3h=50, ppO5c=51, &
    ppO5h=0

   real(RLEN),public,dimension(:),pointer  :: O2o, N1p, N3n, N4n, O4n, N5s, &
    N6r, B1c, B1n, B1p, P1c, P1n, P1p, P1l, P1s, P2c, P2n, P2p, P2l, P3c, P3n, &
    P3p, P3l, P4c, P4n, P4p, P4l, Z3c, Z3n, Z3p, Z4c, Z4n, Z4p, Z5c, Z5n, Z5p, &
    Z6c, Z6n, Z6p, R1c, R1n, R1p, R2c, R3c, R6c, R6n, R6p, R6s, O3c, O3h, O5c


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !! GLOBAL definition of Pelagic (D3) Flux variables
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,parameter,public                :: ppruPPYc=1, ppresPPYc=2, &
    ppresPBAc=3, ppresZOOc=4, pprcalCARc=5, ppruPPYn=6, ppruPPYp=7, ppruPPYs=8, &
    ppexPPYc=9, ppexPPYs=10, ppruZOOc=11, ppresPELo=12, ppremZOOn=13, &
    ppremZOOp=14, ppremPBAn=15, ppremPBAp=16, ppremN5s=17



#ifdef INCLUDE_SEAICE
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !! GLOBAL definition of SeaIce (D2) state variables
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=



#endif

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !! GLOBAL definition of SeaBen (D2) state variables
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,parameter,public                :: ppQ1c=1, ppQ1n=2, ppQ1p=3, &
    ppQ1s=0, ppQ11c=4, ppQ11n=5, ppQ11p=6, ppQ11s=0, ppQ6c=7, ppQ6n=8, ppQ6p=9, &
    ppQ6s=10, ppQ16c=11, ppQ16n=12, ppQ16p=13, ppQ16s=14

   real(RLEN),public,dimension(:),pointer  :: Q1c, Q1n, Q1p, Q11c, Q11n, &
    Q11p, Q6c, Q6n, Q6p, Q6s, Q16c, Q16n, Q16p, Q16s


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Constituent parameters:
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,parameter,public :: iiC=1, iiN=2, iiP=3, iiL=4, iiS=5, iiH=6, &
    iiLastElement=6


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Group parameters (pel):
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,parameter,public :: iiPelBacteria=1, iiB1=1
        integer,parameter,public :: iiPhytoPlankton=4, iiP1=1, iiP2=2, iiP3=3, iiP4=4
        integer,parameter,public :: iiMesoZooPlankton=2, iiZ3=1, iiZ4=2
        integer,parameter,public :: iiMicroZooPlankton=2, iiZ5=1, iiZ6=2
        integer,parameter,public :: iiPelDetritus=4, iiR1=1, iiR2=2, iiR3=3, iiR6=4
        integer,parameter,public :: iiInorganic=2, iiO3=1, iiO5=2


        logical,public :: CalcPelBacteria(iiPelBacteria) = .TRUE.
        logical,public :: CalcPhytoPlankton(iiPhytoPlankton) = .TRUE.
        logical,public :: CalcMesoZooPlankton(iiMesoZooPlankton) = .TRUE.
        logical,public :: CalcMicroZooPlankton(iiMicroZooPlankton) = .TRUE.
        logical,public :: CalcPelDetritus(iiPelDetritus) = .TRUE.
        logical,public :: CalcInorganic(iiInorganic) = .TRUE.


#ifdef INCLUDE_SEAICE
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Group parameters (ice):
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




#endif

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        ! Group parameters (ben):
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,parameter,public :: iiBenDetritus=4, iiQ1=1, iiQ11=2, iiQ6=3, iiQ16=4


        logical,public :: CalcBenDetritus(iiBenDetritus) = .TRUE.


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !  Global Variables (pel)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,public  :: BoxNumber
        integer,public  :: BoxNumberX
        integer,public  :: BoxNumberY
        integer,public  :: BoxNumberZ
        integer,public  :: BoxNumberXY

        real(RLEN),public                              :: LocalDelta ! 

        integer,public                                 :: InitializeModel ! 


#ifdef INCLUDE_SEAICE
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !  Global Variables (ice)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,public  :: BoxNumber_ice
        integer,public  :: BoxNumberX_ice
        integer,public  :: BoxNumberY_ice
        integer,public  :: BoxNumberXY_ice



#endif

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !  Global Variables (ben)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        integer,public  :: BoxNumber_ben
        integer,public  :: BoxNumberX_ben
        integer,public  :: BoxNumberY_ben
        integer,public  :: BoxNumberXY_ben




        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        !! GLOBAL definition of Pelagic (D3/D2) variables which can be outputted in netcdf
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-

!    3d name                                                  description           unit
! ---------- ------------------------------------------------------------ ---------------
!      Depth                                              Gridpoint Depth               m
!     Volume                                             Gridpoint Volume              m3
!       Area                                               Gridpoint Area              m2
!        ETW                                                  temperature               C
!        ESW                                                     Salinity               -
!       ERHO                                             Seawater Density           kg/m3
!        EIR                                                   Irradiance         uE/m2/s
!        ESS                                          Suspended Sediments            g/m3
!        EPR                                               Water Pressure            dbar
!        DIC                                   Dissolved Inorganic Carbon         umol/kg
!        CO2                                                      CO2(aq)         umol/kg
!       pCO2                                                 Oceanic pCO2            uatm
!       HCO3                                                  Bicarbonate         umol/kg
!        CO3                                                    Carbonate         umol/kg
!        ALK                                                   Alkalinity      umol eq/kg
!         pH                                                           pH               -
!      OCalc                                  Saturation state of Calcite               -
!      OArag                                Saturation state of Aragonite               -
!      ffCO2                                                 CO2 fugacity            uatm
!    totpelc                                        Total Mass in Pelagic             g C
!    totpeln                                        Total Mass in Pelagic             g N
!    totpelp                                        Total Mass in Pelagic             g P
!    totpels                                        Total Mass in Pelagic            g Si
!      cxoO2                                            Oxygen Saturation      mmol O2/m3
!     eO2mO2                                   Relative Oxygen saturation               -
!       Chla                                                Chlorophyll-a       mg Chl/m3
!    flPTN6r                        Pelagic Anaerobic Mineralization Rate    mmol O2/m3/d
!    flN3O4n                                 Pelagic Denitrification Rate     mmol N/m3/d
!    flN4N3n                                   Pelagic Nitrification Rate     mmol N/m3/d
!     sediR2                             Colloidal OMT sedimentation rate             m/d
!     sediR3                            Refractory OMT sedimentation rate             m/d
!     sediR6                                  Detritus sedimentation rate             m/d
!     sediO5                                   Calcite sedimentation rate             m/d
!       xEPS                                 Total Extinction Coefficient             1/m
!   ABIO_eps                               Abiotic Extinction Coefficient             1/m
!        BGE                                        Bac growth efficiency               -
!     BIOALK                                O3h biological rate of change         mmol/m3

! qpcPPY(iiP1)                                                      Diatoms     mmol P/mg C
! qpcPPY(iiP2)                                                  Flagellates     mmol P/mg C
! qpcPPY(iiP3)                                            PicoPhytoplankton     mmol P/mg C
! qpcPPY(iiP4)                                          Large Phytoplankton     mmol P/mg C
! qncPPY(iiP1)                                                      Diatoms     mmol N/mg C
! qncPPY(iiP2)                                                  Flagellates     mmol N/mg C
! qncPPY(iiP3)                                            PicoPhytoplankton     mmol N/mg C
! qncPPY(iiP4)                                          Large Phytoplankton     mmol N/mg C
! qscPPY(iiP1)                                                      Diatoms    mmol Si/mg C
! qscPPY(iiP2)                                                  Flagellates    mmol Si/mg C
! qscPPY(iiP3)                                            PicoPhytoplankton    mmol Si/mg C
! qscPPY(iiP4)                                          Large Phytoplankton    mmol Si/mg C
! qlcPPY(iiP1)                                                      Diatoms    mg Chl /mg C
! qlcPPY(iiP2)                                                  Flagellates    mg Chl /mg C
! qlcPPY(iiP3)                                            PicoPhytoplankton    mg Chl /mg C
! qlcPPY(iiP4)                                          Large Phytoplankton    mg Chl /mg C
! qccPPY(iiP1)                                                      Diatoms       mg C/mg C
! qccPPY(iiP2)                                                  Flagellates       mg C/mg C
! qccPPY(iiP3)                                            PicoPhytoplankton       mg C/mg C
! qccPPY(iiP4)                                          Large Phytoplankton       mg C/mg C
! qpcMEZ(iiZ3)                                  Carnivorous Mesozooplankton     mmol P/mg C
! qpcMEZ(iiZ4)                                   Omnivorous Mesozooplankton     mmol P/mg C
! qncMEZ(iiZ3)                                  Carnivorous Mesozooplankton     mmol N/mg C
! qncMEZ(iiZ4)                                   Omnivorous Mesozooplankton     mmol N/mg C
! qpcMIZ(iiZ5)                                             Microzooplankton     mmol P/mg C
! qpcMIZ(iiZ6)                         Heterotrophic Nanoflagellates (HNAN)     mmol P/mg C
! qncMIZ(iiZ5)                                             Microzooplankton     mmol N/mg C
! qncMIZ(iiZ6)                         Heterotrophic Nanoflagellates (HNAN)     mmol N/mg C
! qpcOMT(iiR1)                              Labile Dissolved Organic Matter     mmol N/mg C
! qpcOMT(iiR2)                         Semi-labile Dissolved Organic Carbon     mmol N/mg C
! qpcOMT(iiR3)                     Semi-refractory Dissolved Organic Carbon     mmol N/mg C
! qpcOMT(iiR6)                                   Particulate Organic Matter     mmol N/mg C
! qncOMT(iiR1)                              Labile Dissolved Organic Matter     mmol P/mg C
! qncOMT(iiR2)                         Semi-labile Dissolved Organic Carbon     mmol P/mg C
! qncOMT(iiR3)                     Semi-refractory Dissolved Organic Carbon     mmol P/mg C
! qncOMT(iiR6)                                   Particulate Organic Matter     mmol P/mg C
! qscOMT(iiR1)                              Labile Dissolved Organic Matter    mmol Si/mg C
! qscOMT(iiR2)                         Semi-labile Dissolved Organic Carbon    mmol Si/mg C
! qscOMT(iiR3)                     Semi-refractory Dissolved Organic Carbon    mmol Si/mg C
! qscOMT(iiR6)                                   Particulate Organic Matter    mmol Si/mg C
! qpcPBA(iiB1)                               Aerobic and Anaerobic Bacteria     mmol P/mg C
! qncPBA(iiB1)                               Aerobic and Anaerobic Bacteria     mmol N/mg C
! sediPPY(iiP1)                                                      Diatoms             m/d
! sediPPY(iiP2)                                                  Flagellates             m/d
! sediPPY(iiP3)                                            PicoPhytoplankton             m/d
! sediPPY(iiP4)                                          Large Phytoplankton             m/d
! sediMIZ(iiZ5)                                             Microzooplankton             m/d
! sediMIZ(iiZ6)                         Heterotrophic Nanoflagellates (HNAN)             m/d
! sediMEZ(iiZ3)                                  Carnivorous Mesozooplankton             m/d
! sediMEZ(iiZ4)                                   Omnivorous Mesozooplankton             m/d
! sunPPY(iiP1)                                                      Diatoms             1/d
! sunPPY(iiP2)                                                  Flagellates             1/d
! sunPPY(iiP3)                                            PicoPhytoplankton             1/d
! sunPPY(iiP4)                                          Large Phytoplankton             1/d
! eiPPY(iiP1)                                                      Diatoms               -
! eiPPY(iiP2)                                                  Flagellates               -
! eiPPY(iiP3)                                            PicoPhytoplankton               -
! eiPPY(iiP4)                                          Large Phytoplankton               -
! ELiPPY(iiP1)                                                      Diatoms            W/m2
! ELiPPY(iiP2)                                                  Flagellates            W/m2
! ELiPPY(iiP3)                                            PicoPhytoplankton            W/m2
! ELiPPY(iiP4)                                          Large Phytoplankton            W/m2


        integer,parameter,public :: ppDepth=1, ppVolume=2, ppArea=3, ppETW=4, &
    ppESW=5, ppERHO=6, ppEIR=7, ppESS=8, ppEPR=9, ppDIC=10, ppCO2=11, &
    pppCO2=12, ppHCO3=13, ppCO3=14, ppALK=15, pppH=16, ppOCalc=17, ppOArag=18, &
    ppffCO2=19, pptotpelc=20, pptotpeln=21, pptotpelp=22, pptotpels=23, &
    ppcxoO2=24, ppeO2mO2=25, ppChla=26, ppflPTN6r=27, ppflN3O4n=28, &
    ppflN4N3n=29, ppsediR2=30, ppsediR3=31, ppsediR6=32, ppsediO5=33, &
    ppxEPS=34, ppABIO_eps=35, ppBGE=36, ppBIOALK=37

        integer,public ::           ppqpcPPY(iiPhytoPlankton), &
    ppqncPPY(iiPhytoPlankton), ppqscPPY(iiPhytoPlankton), &
    ppqlcPPY(iiPhytoPlankton), ppqccPPY(iiPhytoPlankton), &
    ppqpcMEZ(iiMesoZooPlankton), ppqncMEZ(iiMesoZooPlankton), &
    ppqpcMIZ(iiMicroZooPlankton), ppqncMIZ(iiMicroZooPlankton), &
    ppqpcOMT(iiPelDetritus), ppqncOMT(iiPelDetritus), ppqscOMT(iiPelDetritus), &
    ppqpcPBA(iiPelBacteria), ppqncPBA(iiPelBacteria), &
    ppsediPPY(iiPhytoPlankton), ppsediMIZ(iiMicroZooPlankton), &
    ppsediMEZ(iiMesoZooPlankton), ppsunPPY(iiPhytoPlankton), &
    ppeiPPY(iiPhytoPlankton), ppELiPPY(iiPhytoPlankton)

   real(RLEN),public,dimension(:),pointer  :: Depth, Volume, Area, ETW, ESW, &
    ERHO, EIR, ESS, EPR, DIC, CO2, pCO2, HCO3, CO3, ALK, pH, OCalc, OArag, &
    ffCO2, totpelc, totpeln, totpelp, totpels, cxoO2, eO2mO2, Chla, flPTN6r, &
    flN3O4n, flN4N3n, sediR2, sediR3, sediR6, sediO5, xEPS, ABIO_eps, BGE, &
    BIOALK

        real(RLEN),public,dimension(:,:),pointer  :: qpcPPY, qncPPY, qscPPY, &
    qlcPPY, qccPPY, qpcMEZ, qncMEZ, qpcMIZ, qncMIZ, qpcOMT, qncOMT, qscOMT, &
    qpcPBA, qncPBA, sediPPY, sediMIZ, sediMEZ, sunPPY, eiPPY, ELiPPY


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

!    2d name                                                  description           unit
! ---------- ------------------------------------------------------------ ---------------
!     Area2d                                           2-D Gridpoint Area              m2
! ThereIsLight                                   Switch for day/night cycle               -
!       SUNQ                                           Daylength in hours               h
!      EWIND                                                   Wind speed             m/s
!       EICE                                             Sea-ice fraction               -
!      ETAUB                                                Bottom Stress           N m/s
!   EPCO2air                             Atmospheric CO2 Partial Pressure            uatm
! CO2airflux                                             Sea-air CO2 Flux       mmol/m2/d
!      dpco2                                      Air-Sea pCO2 difference            uatm
!       dpo2                                       Air-Sea pO2 difference             atm
!    totsysc                                                   total mass             g C
!    totsysn                                                   total mass             g N
!    totsysp                                                   total mass             g P
!    totsyss                                                   total mass            g Si



        integer,parameter,public :: ppArea2d=1, ppThereIsLight=2, ppSUNQ=3, &
    ppEWIND=4, ppEICE=5, ppETAUB=6, ppEPCO2air=7, ppCO2airflux=8, ppdpco2=9, &
    ppdpo2=10, pptotsysc=11, pptotsysn=12, pptotsysp=13, pptotsyss=14


   real(RLEN),public,dimension(:),pointer  :: Area2d, ThereIsLight, SUNQ, &
    EWIND, EICE, ETAUB, EPCO2air, CO2airflux, dpco2, dpo2, totsysc, totsysn, &
    totsysp, totsyss



#ifdef INCLUDE_SEAICE
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        !! GLOBAL definition of Seaice (D2) variables which can be outputted in netcdf
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-







#endif

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        !! GLOBAL definition of benthic (D2) variables which can be outputted in netcdf
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-
!    2d name                                                  description           unit
! ---------- ------------------------------------------------------------ ---------------
!    totbenc                                total mass present in benthos             g C
!    totbenn                                total mass present in benthos             g N
!    totbenp                                total mass present in benthos             g P
!    totbens                                total mass present in benthos            g Si

! qpcBOM(iiQ1)                       Surface Sediment Labile organic Matter     mmol N/mg C
! qpcBOM(iiQ11)                          Deep Sediment Labile organic Matter     mmol N/mg C
! qpcBOM(iiQ6)                  Surface Sediment Particulate organic Matter     mmol N/mg C
! qpcBOM(iiQ16)                     Deep Sediment Particulate organic Matter     mmol N/mg C
! qncBOM(iiQ1)                       Surface Sediment Labile organic Matter     mmol P/mg C
! qncBOM(iiQ11)                          Deep Sediment Labile organic Matter     mmol P/mg C
! qncBOM(iiQ6)                  Surface Sediment Particulate organic Matter     mmol P/mg C
! qncBOM(iiQ16)                     Deep Sediment Particulate organic Matter     mmol P/mg C
! qscBOM(iiQ1)                       Surface Sediment Labile organic Matter    mmol Si/mg C
! qscBOM(iiQ11)                          Deep Sediment Labile organic Matter    mmol Si/mg C
! qscBOM(iiQ6)                  Surface Sediment Particulate organic Matter    mmol Si/mg C
! qscBOM(iiQ16)                     Deep Sediment Particulate organic Matter    mmol Si/mg C


        integer,parameter,public :: pptotbenc=1, pptotbenn=2, pptotbenp=3, &
    pptotbens=4

        integer,public ::           ppqpcBOM(iiBenDetritus), &
    ppqncBOM(iiBenDetritus), ppqscBOM(iiBenDetritus)

   real(RLEN),public,dimension(:),pointer  :: totbenc, totbenn, totbenp, &
    totbens

        real(RLEN),public,dimension(:,:),pointer  :: qpcBOM, qncBOM, qscBOM


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !  boundary fluxes
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        real(RLEN),public,dimension(:),pointer  :: jsurO2o, jsurN1p, jsurN3n, &
    jsurN4n, jsurO4n, jsurN5s, jsurN6r, jsurB1c, jsurB1n, jsurB1p, jsurP1c, &
    jsurP1n, jsurP1p, jsurP1l, jsurP1s, jsurP2c, jsurP2n, jsurP2p, jsurP2l, &
    jsurP3c, jsurP3n, jsurP3p, jsurP3l, jsurP4c, jsurP4n, jsurP4p, jsurP4l, &
    jsurZ3c, jsurZ3n, jsurZ3p, jsurZ4c, jsurZ4n, jsurZ4p, jsurZ5c, jsurZ5n, &
    jsurZ5p, jsurZ6c, jsurZ6n, jsurZ6p, jsurR1c, jsurR1n, jsurR1p, jsurR2c, &
    jsurR3c, jsurR6c, jsurR6n, jsurR6p, jsurR6s, jsurO3c, jsurO3h, jsurO5c

        real(RLEN),public,dimension(:),pointer  :: jbotO2o, jbotN1p, jbotN3n, &
    jbotN4n, jbotO4n, jbotN5s, jbotN6r, jbotB1c, jbotB1n, jbotB1p, jbotP1c, &
    jbotP1n, jbotP1p, jbotP1l, jbotP1s, jbotP2c, jbotP2n, jbotP2p, jbotP2l, &
    jbotP3c, jbotP3n, jbotP3p, jbotP3l, jbotP4c, jbotP4n, jbotP4p, jbotP4l, &
    jbotZ3c, jbotZ3n, jbotZ3p, jbotZ4c, jbotZ4n, jbotZ4p, jbotZ5c, jbotZ5n, &
    jbotZ5p, jbotZ6c, jbotZ6n, jbotZ6p, jbotR1c, jbotR1n, jbotR1p, jbotR2c, &
    jbotR3c, jbotR6c, jbotR6n, jbotR6p, jbotR6s, jbotO3c, jbotO3h, jbotO5c


        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !  Other 3d-Global Variables 
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=




        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !  Other 2d-Global Variables (pel)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=       




#ifdef INCLUDE_SEAICE
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !  Other 2d-Global Variables (ice)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=       



#endif

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !  Other 2d-Global Variables (ben)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=       




        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !! SHARED GLOBAL FUNCTIONS (pel) (must be below contains)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        public flux, flux_vector, Source, Source_D3_vector, quota_flux

        public ppPelBacteria, PelBacteria, ppPhytoPlankton, PhytoPlankton, &
    ppMesoZooPlankton, MesoZooPlankton, ppMicroZooPlankton, MicroZooPlankton, &
    ppPelDetritus, PelDetritus, ppInorganic, Inorganic


#if defined INCLUDE_SEAICE
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !! SHARED GLOBAL FUNCTIONS (ice) (must be below contains)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        public Source_D2_vector_ice


#endif

        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        !! SHARED GLOBAL FUNCTIONS (ben) (must be below contains)
        !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        public Source_D2_vector_ben

        public ppBenDetritus, BenDetritus

        contains

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          !! Group Pelagic (D3) state functions
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


   function ppPelBacteria(n,constituent)
   
     IMPLICIT NONE
   
     integer ::ppPelBacteria
     integer, intent(IN)            :: n
     integer, intent(IN)            :: constituent
     
     integer,dimension(1 * iiLastElement ) :: pointers = (/ & 
             ppB1c ,ppB1n ,ppB1p ,0     ,0     ,0       &
          /)

     IF( n > 1 .OR. n == 0 ) THEN
       ppPelBacteria = 0
     ELSE IF( constituent > iiLastElement .OR. constituent == 0 ) THEN
       ppPelBacteria = 0
     ELSE
       ppPelBacteria = pointers( ( (n-1) * iiLastElement ) + constituent )
     ENDIF
   
   END function

   function ppPhytoPlankton(n,constituent)
   
     IMPLICIT NONE
   
     integer ::ppPhytoPlankton
     integer, intent(IN)            :: n
     integer, intent(IN)            :: constituent
     
     integer,dimension(4 * iiLastElement ) :: pointers = (/ & 
             ppP1c ,ppP1n ,ppP1p ,ppP1l ,ppP1s ,0     , &
          ppP2c ,ppP2n ,ppP2p ,ppP2l ,0     ,0     , &
          ppP3c ,ppP3n ,ppP3p ,ppP3l ,0     ,0     , &
          ppP4c ,ppP4n ,ppP4p ,ppP4l ,0     ,0       &
          /)

     IF( n > 4 .OR. n == 0 ) THEN
       ppPhytoPlankton = 0
     ELSE IF( constituent > iiLastElement .OR. constituent == 0 ) THEN
       ppPhytoPlankton = 0
     ELSE
       ppPhytoPlankton = pointers( ( (n-1) * iiLastElement ) + constituent )
     ENDIF
   
   END function

   function ppMesoZooPlankton(n,constituent)
   
     IMPLICIT NONE
   
     integer ::ppMesoZooPlankton
     integer, intent(IN)            :: n
     integer, intent(IN)            :: constituent
     
     integer,dimension(2 * iiLastElement ) :: pointers = (/ & 
             ppZ3c ,ppZ3n ,ppZ3p ,0     ,0     ,0     , &
          ppZ4c ,ppZ4n ,ppZ4p ,0     ,0     ,0       &
          /)

     IF( n > 2 .OR. n == 0 ) THEN
       ppMesoZooPlankton = 0
     ELSE IF( constituent > iiLastElement .OR. constituent == 0 ) THEN
       ppMesoZooPlankton = 0
     ELSE
       ppMesoZooPlankton = pointers( ( (n-1) * iiLastElement ) + constituent )
     ENDIF
   
   END function

   function ppMicroZooPlankton(n,constituent)
   
     IMPLICIT NONE
   
     integer ::ppMicroZooPlankton
     integer, intent(IN)            :: n
     integer, intent(IN)            :: constituent
     
     integer,dimension(2 * iiLastElement ) :: pointers = (/ & 
             ppZ5c ,ppZ5n ,ppZ5p ,0     ,0     ,0     , &
          ppZ6c ,ppZ6n ,ppZ6p ,0     ,0     ,0       &
          /)

     IF( n > 2 .OR. n == 0 ) THEN
       ppMicroZooPlankton = 0
     ELSE IF( constituent > iiLastElement .OR. constituent == 0 ) THEN
       ppMicroZooPlankton = 0
     ELSE
       ppMicroZooPlankton = pointers( ( (n-1) * iiLastElement ) + constituent )
     ENDIF
   
   END function

   function ppPelDetritus(n,constituent)
   
     IMPLICIT NONE
   
     integer ::ppPelDetritus
     integer, intent(IN)            :: n
     integer, intent(IN)            :: constituent
     
     integer,dimension(4 * iiLastElement ) :: pointers = (/ & 
             ppR1c ,ppR1n ,ppR1p ,0     ,0     ,0     , &
          ppR2c ,0     ,0     ,0     ,0     ,0     , &
          ppR3c ,0     ,0     ,0     ,0     ,0     , &
          ppR6c ,ppR6n ,ppR6p ,0     ,ppR6s ,0       &
          /)

     IF( n > 4 .OR. n == 0 ) THEN
       ppPelDetritus = 0
     ELSE IF( constituent > iiLastElement .OR. constituent == 0 ) THEN
       ppPelDetritus = 0
     ELSE
       ppPelDetritus = pointers( ( (n-1) * iiLastElement ) + constituent )
     ENDIF
   
   END function

   function ppInorganic(n,constituent)
   
     IMPLICIT NONE
   
     integer ::ppInorganic
     integer, intent(IN)            :: n
     integer, intent(IN)            :: constituent
     
     integer,dimension(2 * iiLastElement ) :: pointers = (/ & 
             ppO3c ,0     ,0     ,0     ,0     ,ppO3h , &
          ppO5c ,0     ,0     ,0     ,0     ,0       &
          /)

     IF( n > 2 .OR. n == 0 ) THEN
       ppInorganic = 0
     ELSE IF( constituent > iiLastElement .OR. constituent == 0 ) THEN
       ppInorganic = 0
     ELSE
       ppInorganic = pointers( ( (n-1) * iiLastElement ) + constituent )
     ENDIF
   
   END function

   function PelBacteria(n,constituent)
   
     IMPLICIT NONE
   
     real(RLEN),dimension(:),pointer ::PelBacteria
     integer, intent(IN)            :: n
     integer, intent(IN)            :: constituent
     
     if ( ppPelBacteria(n,constituent) > 0 ) then
       PelBacteria => D3STATE(:,ppPelBacteria(n,constituent))
     else
       PelBacteria => NULL()
     end if
   
   END function

   function PhytoPlankton(n,constituent)
   
     IMPLICIT NONE
   
     real(RLEN),dimension(:),pointer ::PhytoPlankton
     integer, intent(IN)            :: n
     integer, intent(IN)            :: constituent
     
     if ( ppPhytoPlankton(n,constituent) > 0 ) then
       PhytoPlankton => D3STATE(:,ppPhytoPlankton(n,constituent))
     else
       PhytoPlankton => NULL()
     end if
   
   END function

   function MesoZooPlankton(n,constituent)
   
     IMPLICIT NONE
   
     real(RLEN),dimension(:),pointer ::MesoZooPlankton
     integer, intent(IN)            :: n
     integer, intent(IN)            :: constituent
     
     if ( ppMesoZooPlankton(n,constituent) > 0 ) then
       MesoZooPlankton => D3STATE(:,ppMesoZooPlankton(n,constituent))
     else
       MesoZooPlankton => NULL()
     end if
   
   END function

   function MicroZooPlankton(n,constituent)
   
     IMPLICIT NONE
   
     real(RLEN),dimension(:),pointer ::MicroZooPlankton
     integer, intent(IN)            :: n
     integer, intent(IN)            :: constituent
     
     if ( ppMicroZooPlankton(n,constituent) > 0 ) then
       MicroZooPlankton => D3STATE(:,ppMicroZooPlankton(n,constituent))
     else
       MicroZooPlankton => NULL()
     end if
   
   END function

   function PelDetritus(n,constituent)
   
     IMPLICIT NONE
   
     real(RLEN),dimension(:),pointer ::PelDetritus
     integer, intent(IN)            :: n
     integer, intent(IN)            :: constituent
     
     if ( ppPelDetritus(n,constituent) > 0 ) then
       PelDetritus => D3STATE(:,ppPelDetritus(n,constituent))
     else
       PelDetritus => NULL()
     end if
   
   END function

   function Inorganic(n,constituent)
   
     IMPLICIT NONE
   
     real(RLEN),dimension(:),pointer ::Inorganic
     integer, intent(IN)            :: n
     integer, intent(IN)            :: constituent
     
     if ( ppInorganic(n,constituent) > 0 ) then
       Inorganic => D3STATE(:,ppInorganic(n,constituent))
     else
       Inorganic => NULL()
     end if
   
   END function


#if defined INCLUDE_SEAICE
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          !! Group Seaice (D2) state functions
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


#endif

          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
          !! Group Benthic (D2) state functions
          !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


   function ppBenDetritus(n,constituent)
   
     IMPLICIT NONE
   
     integer ::ppBenDetritus
     integer, intent(IN)            :: n
     integer, intent(IN)            :: constituent
     
     integer,dimension(4 * iiLastElement ) :: pointers = (/ & 
             ppQ1c ,ppQ1n ,ppQ1p ,0     ,0     ,0     , &
          ppQ11c ,ppQ11n ,ppQ11p ,0     ,0     ,0     , &
          ppQ6c ,ppQ6n ,ppQ6p ,0     ,ppQ6s ,0     , &
          ppQ16c ,ppQ16n ,ppQ16p ,0     ,ppQ16s ,0       &
          /)

     IF( n > 4 .OR. n == 0 ) THEN
       ppBenDetritus = 0
     ELSE IF( constituent > iiLastElement .OR. constituent == 0 ) THEN
       ppBenDetritus = 0
     ELSE
       ppBenDetritus = pointers( ( (n-1) * iiLastElement ) + constituent )
     ENDIF
   
   END function

   function BenDetritus(n,constituent)
   
     IMPLICIT NONE
   
     real(RLEN),dimension(:),pointer ::BenDetritus
     integer, intent(IN)            :: n
     integer, intent(IN)            :: constituent
     
     if ( ppBenDetritus(n,constituent) > 0 ) then
       BenDetritus => D2STATE_BEN(:,ppBenDetritus(n,constituent))
     else
       BenDetritus => NULL()
     end if
   
   END function



#include "FluxFunctions.h90"

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    ! end of contain section
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  end module mem

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
