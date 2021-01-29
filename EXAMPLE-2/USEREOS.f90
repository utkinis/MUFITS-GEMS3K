module USEREOS
    
    implicit none
    
    private
    
    integer,parameter::USEREOS_ncmax=2  ! Maximum value of nc (nc<=ncmax)
    integer,parameter::USEREOS_npmax=2  ! Maximum value of np (np<=npmax)
    
    integer,parameter::USEREOS_nc=2              ! - Number of components    
    integer,parameter::USEREOS_np=USEREOS_nc     ! - Maximum number of phases

    !<<<<<<<<<<<<<<<<<<< These data are loaded from the configuration file <<<<<<<<<<<<<<<<<<<<<<<<
    real(8)::Pref=0./0.                         ! - Reference pressure
    real(8)::Tref=0./0.                         ! - Reference temperature
    
    character(3)::CompNames(USEREOS_ncmax)=''   ! - The components names (3-byte characters) are also the phase names
    real(8)::MW(USEREOS_ncmax)=0./0.            ! - Molar weight of each component [g/mol]
    real(8)::DEN(USEREOS_ncmax)=0./0.           ! - Density at the reference pressure and temperature
    real(8)::HC(USEREOS_ncmax)=0./0.            ! - Molar heat capacity of each component  [kJ/mol]
    real(8)::Visc(USEREOS_ncmax)=0./0.          ! - Viscosity of each phase [cP]
    real(8)::Alpha(USEREOS_ncmax)=0./0.         ! - isothermal compressibility [1/bar]
    real(8)::Beta(USEREOS_ncmax)=0./0.          ! - heat expansion coefficient [1/K]
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    
    !<<<<<<<<<<< The shared library must export the following variables & routines <<<<<<<<<<<<<<<<
    public::USEREOS_ReadConfigurationFile   ! - reads the configuration file. The file name is supplied by MUFITS
    public::USEREOS_GetDimensions    ! - returns to MUFITS the number of components and the maximum number of phases
    public::USEREOS_GetGlobalParameters     ! - returns to MUFITS some general parameters of thermodynamic system
    public::USEREOS_PhaseEquilibrium             ! - calculates the phase equilibrium for a given pressure, temperature,
                                     !   and bulk composition
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
contains
    
    
    subroutine USEREOS_ReadConfigurationFile(filename,ierr,title) 
    !DEC$ ATTRIBUTES DLLEXPORT::USEREOS_ReadConfigurationFile
    
        character(*),intent(in)::filename        ! name of the configuration file specified by the USEREOS keyword.
        integer,intent(out)::ierr                ! if 0, then no errors
        character(80),intent(out)::title         ! Describes the shared library. This string is reported into the LOG-file
        
        !-----------------------------------------------------------------------------------------------------------

        integer,parameter::iport=3891
        integer::ic,jc
        logical::ex
        
        title='IMMISCIBLE TWO-PHASE FLOW'
        
        !<<<<<<<<<<<<<<<<<<<< reading the configuration file <<<<<<<<<<<<<<<<<<<<<<<<<<<<
        inquire(file=filename,exist=ex)
        if(not(ex)) then
            ierr=-1
            return
        endif            
        
        open(iport,file=filename,status='old')
        
        read(iport,*)  Pref,Tref;  Pref=Pref*1e5
        read(iport,*)  CompNames(1:USEREOS_nc)
        read(iport,*)  MW(1:USEREOS_nc);   MW=MW*1e-3
        read(iport,*)  DEN(1:USEREOS_nc)
        read(iport,*)  ALPHA(1:USEREOS_nc);  ALPHA=ALPHA/1e5
        read(iport,*)  BETA(1:USEREOS_nc)
        read(iport,*)  HC(1:USEREOS_nc);   HC=HC*1e3
        read(iport,*)  VISC(1:USEREOS_nc); VISC=VISC*1e-3        
        close(iport)
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        !<<<<<<<<<<<<<<<<<<<<<<<< checking data in the configuration file <<<<<<<<<<<<<<<
        if(ISNAN(Pref).or.ISNAN(Tref)) then
            ierr=2
            return
        endif            
        
        do ic=1,USEREOS_nc
            if(trim(CompNames(ic)).eq.'') then
                ierr=2
                return
            endif                
            if(MW(ic).le.0d0.or.ISNAN(MW(ic))) then
                ierr=2
                return
            endif     
            if(DEN(ic).le.0d0.or.ISNAN(DEN(ic))) then
                ierr=2
                return
            endif   
            if(ALPHA(ic).lt.0d0.or.ISNAN(ALPHA(ic))) then
                ierr=2
                return
            endif                 
            if(BETA(ic).lt.0d0.or.ISNAN(BETA(ic))) then
                ierr=2
                return
            endif                    
            if(HC(ic).le.0d0.or.ISNAN(HC(ic))) then
                ierr=2
                return
            endif        
            if(VISC(ic).le.0d0.or.ISNAN(VISC(ic))) then
                ierr=2
                return
            endif               
        enddo            
        
        do ic=2,USEREOS_nc
            do jc=1,ic-1
                if(trim(CompNames(ic)).eq.trim(CompNames(jc))) then
                    ierr=6
                    return
                endif                    
            enddo
        enddo       
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        ierr=0

    end subroutine USEREOS_ReadConfigurationFile
    
    !##############################################################################################################
    
    subroutine USEREOS_GetDimensions(NComponents,NPhaseMax,NAuxArray) 
    !DEC$ ATTRIBUTES DLLEXPORT::USEREOS_GetDimensions
    
        integer,intent(out)::NComponents                    ! - number of components of the thermodynamic system
        integer,intent(out)::NPhaseMax                      ! - maximum number of phases that the mixture can split into    
        integer,intent(out)::NAuxArray                      ! - number of auxiliary arrays
        !---------------------------------------------------------------------------------------------------------- 
        NComponents=USEREOS_nc
        NPhaseMax=USEREOS_np  
        NAuxArray=0
    
    end subroutine USEREOS_GetDimensions
    
    !##############################################################################################################
    
    subroutine USEREOS_GetGlobalParameters(nc,np,na,CMPnames,MolWeights,PHnames, &
                                                                 AuxArrNames,AuxArrDimens,opt) 
    !DEC$ ATTRIBUTES DLLEXPORT::USEREOS_GetGlobalParameters
        integer,intent(in)::nc,np,na
        character(3)::CMPnames(nc)                          ! - names of the components
        real(8)::MolWeights(nc)                             ! - molar weights
        character(3)::PHnames(np)                           ! - names of the phases
        character(8)::AuxArrNames(na)                       ! - names of auxiliary quantities. Each name must begin with '#'
        character(8)::AuxArrDimens(na)                      ! - dimensions of auxiliary quantities
        integer(1),intent(out)::opt(8)                      ! - Options (an option #i is on if opt(i)=1. Otherwise if opt(i)=0 the option is off.
                                                            !   if opt(1)=1, then the shared library must export the relative permeability
                                                            !   if opt(1)=0, then the relative permeability is specified in the MUFITS RUN-file.
                                                            !   if opt(2)=1, then the shared library must export the relative phae pressure 
                                                            !   i.e. the capillary pressure
                                                            !   if opt(2)=0, then the capillary pressure is specified in the MUFITS RUN-file.
        !----------------------------------------------------------------------------------------------------------

        CMPnames(1:USEREOS_nc)=CompNames(1:USEREOS_nc)
        MolWeights(1:USEREOS_nc)=MW(1:USEREOS_nc)
        
        PHnames(1:USEREOS_np)=CompNames(1:USEREOS_nc)
        
        opt=0
        
        !opt(1)=1  !- uncomment this line if the relative permeability will be exported from the shared library
        !opt(2)=1  !- uncomment this line if the capillary pressure curve will be exported from the shared library

    end subroutine USEREOS_GetGlobalParameters
    
    !##############################################################################################################
    
    subroutine USEREOS_PhaseEquilibrium(nc,np,na,Pres,Temp,z,NPhase,Props,PhaseID,AuxArr,Mode)
    !DEC$ ATTRIBUTES DLLEXPORT::USEREOS_PhaseEquilibrium
        integer,intent(in)::nc,np,na
        real(8),intent(in)::Pres    ! pressure  [Pa]
        real(8),intent(in)::Temp    ! temperature [K]
        real(8),intent(in)::z(nc)   ! bulk (molar at opt(1)=1, or mass at opt(1)=0) composition 
        real(8)::AuxArr(na)         ! auxiliary quantities 
                
        integer(1),intent(in)::Mode  ! if Mode=1, then the parameters below are used only for reporting fluid properties
                                     ! from the shared library to MUFITS
                                     ! if Mode=0, then the parameters below are used for 
                                     ! 1. importing from MUFITS the initial guess for the phase equilibrium
                                     ! 2. exporting from the shared library to MUFITS the phase equilibrium
        
                                    ! IMPORTANT NOTES: 
                                    ! 1. np must not be changed at Mode=0
                                    ! 2. the saturations are allowed to be negative at Mode=0
                                    ! 3. the order of phases must not change at Mode=0 (e.g. if phase #1 is gas on input, 
                                    !    then it must be gas on output)
        
        integer::NPhase
                                  
        real(8)::Props(6+nc,np)    ! Parameters of the phases
                                    ! 1 - density [kg/m3]
                                    ! 2 - enthalpy [J/kg]
                                    ! 3 - viscosity [Pa*s]
                                    ! 4 - saturation [-]
                                    ! 5 - relative permeability [-] (must be specified only if opt(1)=1)
                                    ! 6 - relative phase pressure [Pa] (must be specified only if opt(2)=1)
                                    ! 7:6+USEREOS_nc - phase composition [-] 
        
        integer(1)::PhaseID(np)  ! phase ID
        
        !--------------------------------------------------------------------------------------------------------------------
        
        integer::ip

        NPhase=2
        Props(7:8,1:2)=0d0
        
        !--- calculating density, enthalpy and viscosity of each phase ---
        do ip=1,np
            Props(1,ip)=DEN(ip)*(1d0+ALPHA(ip)*(Pres-Pref)-BETA(ip)*(Temp-Tref))
            Props(2,ip)=HC(ip)*(Temp-Tref)
            Props(3,ip)=VISC(ip)
            Props(6+ip,ip)=1d0
            PhaseID(ip)=ip
        enddo            
        
        !--- calculating saturation of each phase ---
        Props(4,1)=Props(1,2)*z(1)/(Props(1,2)*z(1)+Props(1,1)*z(2))
        Props(4,2)=1d0-Props(4,1)
        
        !--- calculating relative permeability and capillary pressure ---
        Props(5,1)=Props(4,1)**2
        Props(5,2)=Props(4,2)**2
        Props(6,1)=0d0
        Props(6,2)=Props(4,2)*1e5

    end subroutine USEREOS_PhaseEquilibrium
    
end module USEREOS    