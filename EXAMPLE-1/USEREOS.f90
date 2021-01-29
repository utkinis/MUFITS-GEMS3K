module USEREOS
    
    implicit none
    
    private
    
    integer,parameter::USEREOS_ncmax=16  ! Maximum value of nc (nc<=ncmax)
    integer,parameter::USEREOS_npmax=1   ! Maximum value of np (np<=npmax)
    
    real(8),parameter::Rgas=8.314299882  ! Universal gas constant 
    
    !<<<<<<<<<<<<<<<<<<< These data are loaded from the configuration file <<<<<<<<<<<<<<<<<<<<<<<<
    integer::USEREOS_nc=-1                      ! - Number of components
    character(3)::CompNames(USEREOS_ncmax)=''   ! - The components names (3-byte characters)
    real(8)::MW(USEREOS_ncmax)=0./0.            ! - Molar weight of each component [g/mol]
    real(8)::HC(USEREOS_ncmax)=0./0.            ! - Molar heat capacity of each component  [kJ/mol]
    real(8)::Visc(USEREOS_ncmax)=0./0.          ! - Viscosity of each gas [cP]
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    integer,parameter::USEREOS_np=1   ! Maximum number of phases
    character(3),parameter::PhaseNames(USEREOS_np)='GAS'
    
    !<<<<<<<<<<< The shared library must export the following variables & routines <<<<<<<<<<<<<<<<
    public::USEREOS_ReadConfigurationFile   ! - reads the configuration file. The file name is supplied by MUFITS
    public::USEREOS_GetDimensions           ! - returns to MUFITS the number of components and the maximum number of phases
    public::USEREOS_GetGlobalParameters     ! - returns to MUFITS some general parameters of thermodynamic system
    public::USEREOS_PhaseEquilibrium        ! - calculates the phase equilibrium for a given pressure, temperature,
                                            !   and bulk composition
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
contains
    
    subroutine USEREOS_ReadConfigurationFile(filename,ierr,title) 
    !DEC$ ATTRIBUTES DLLEXPORT::USEREOS_ReadConfigurationFile
    
        character(*),intent(in)::filename        ! path the configuration file specified by the USEREOS keyword.
        integer,intent(out)::ierr                ! if 0, then no errors
        character(80),intent(out)::title         ! Describes the shared library. This string is reported into the LOG-file
        
        !-----------------------------------------------------------------------------------------------------------

        integer,parameter::iport=3891
        integer::ic,jc
        logical::ex
        
        title='MIXTURE OF IDEAL GASES'

        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< reading the configuration file <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        inquire(file=filename,exist=ex)
        if(not(ex)) then
            ierr=-1
            return
        endif            
        open(iport,file=filename,status='old')
        
        read(iport,*)  USEREOS_nc
        if(USEREOS_nc.lt.1.or.USEREOS_nc.gt.USEREOS_ncmax) then
            ierr=1
            close(iport)
            return
        endif    
        
        read(iport,*)  CompNames(1:USEREOS_nc)
        read(iport,*)  MW(1:USEREOS_nc);   MW=MW*1e-3
        read(iport,*)  HC(1:USEREOS_nc);   
        read(iport,*)  VISC(1:USEREOS_nc); VISC=VISC*1e-3        
        close(iport)        
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< checking the data in the config file <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        do ic=1,USEREOS_nc
            if(trim(CompNames(ic)).eq.'') then
                ierr=2
                return
            endif                
            if(MW(ic).le.0d0.or.ISNAN(MW(ic))) then
                ierr=3
                return
            endif     
            if(HC(ic).le.0d0.or.ISNAN(HC(ic))) then
                ierr=4
                return
            endif        
            if(VISC(ic).le.0d0.or.ISNAN(VISC(ic))) then
                ierr=5
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
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
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
        NPhaseMax=USEREOS_npmax 
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
        
        PHnames(1)=PhaseNames(1)
        
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
        
        integer::ic
        real(8)::x(nc)  ! molar concentrations
        real(8)::M
        
        
        !--- caclulating molar composition ---
        M=0d0
        do ic=1,nc
            x(ic)=z(ic)/MW(ic)
            M=M+x(ic)
        enddo    
        x=x/M

        
        !--- calculating molar weight ---
        M=0d0
        do ic=1,nc
            M=M+MW(ic)*x(ic)
        enddo              
    
        !--- calculating density ---
        Props(1,1)=Rgas*Temp/Pres
        Props(1,1)=M/Props(1,1)    
        
        NPhase=1
        
        !--- calculating enthalpy and viscosity --
        Props(2:3,1)=0d0
        do ic=1,USEREOS_nc
            Props(2,1)=Props(2,1)+HC(ic)*x(ic)
            Props(3,1)=Props(3,1)+Visc(ic)*x(ic)
        enddo
        Props(2,1)=Props(2,1)*Temp/M
        Props(4,1)=1d0
        
        !--- calculating composition ---
        Props(7:6+USEREOS_nc,1)=z
        
        PhaseID(1)=1   ! is the GAS phase

    end subroutine USEREOS_PhaseEquilibrium
    
end module USEREOS    