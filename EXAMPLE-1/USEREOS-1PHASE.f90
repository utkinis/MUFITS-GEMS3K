module USEREOS
    
    implicit none
    
    private
    
    integer,parameter::ncmax=16  ! Maximum value of nc (nc<=ncmax)
    integer,parameter::npmax=1   ! Maximum value of np (np<=npmax)
    
    real(8),parameter::Rgas=8.314299882  ! Universal gas constant 
    
    !<<<<<<<<<<<<<<<<<<< These data are loaded from the configuration file <<<<<<<<<<<<<<<<<<<<<<<<
    integer::nc=-1                      ! - Number of components
    character(3)::CompNames(ncmax)=''   ! - The components names (3-byte characters)
    real(8)::MW(ncmax)=0./0.            ! - Molar weight of each component [g/mol]
    real(8)::HC(ncmax)=0./0.            ! - Molar heat capacity of each component  [kJ/mol]
    real(8)::Visc(ncmax)=0./0.          ! - Viscosity of each gas [cP]
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    character(3),parameter::PhaseNames(npmax)='GAS'
    
    integer::na=0   ! Size of the auxiliary array
    
    !<<<<<<<<<<< The shared library must export the following variables & routines <<<<<<<<<<<<<<<<
    public::ReadConfigurationFile           ! - reads the configuration file. The file name is supplied by MUFITS
    public::GetDimensions                   ! - returns to MUFITS the number of components and the maximum number of phases
    public::GetGlobalParameters             ! - returns to MUFITS some general parameters of the thermodynamic system
    public::PhaseEquilibrium                ! - calculates the phase equilibrium for a given pressure, temperature,
                                            !   and bulk composition
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
contains
    
    subroutine ReadConfigurationFile(filename,ierr,title) bind(c,name='ReadConfigurationFile')
    !DEC$ ATTRIBUTES DLLEXPORT ::ReadConfigurationFile
    
        character,intent(in)::filename(256)      ! path the configuration file specified by the USEREOS keyword.
        integer,intent(out)::ierr                ! if 0, then no errors
        character,intent(out)::title(80)         ! Describes the shared library. This string is reported into the LOG-file
        
        !-----------------------------------------------------------------------------------------------------------
        character(256)::filename_internal
        character(80)::title_internal

        integer,parameter::iport=3891
        integer::ic,jc,i
        logical::ex
        
        title_internal='MIXTURE OF IDEAL GASES'
        
        do i=1,80
            title(i)=title_internal(i:i)
        enddo            
        do i=1,256
            filename_internal(i:i)=filename(i)
        enddo             

        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< reading the configuration file <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        inquire(file=filename_internal,exist=ex)
        if(not(ex)) then
            ierr=-1
            return
        endif            
        open(iport,file=filename_internal,status='old')
        
        read(iport,*)  nc
        if(nc.lt.1.or.nc.gt.ncmax) then
            ierr=1
            close(iport)
            return
        endif    
        
        read(iport,*)  CompNames(1:nc)
        read(iport,*)  MW(1:nc);   MW=MW*1e-3
        read(iport,*)  HC(1:nc);   
        read(iport,*)  VISC(1:nc); VISC=VISC*1e-3        
        close(iport)        
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< checking the data in the config file <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        do ic=1,nc
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
        
        do ic=2,nc
            do jc=1,ic-1
                if(trim(CompNames(ic)).eq.trim(CompNames(jc))) then
                    ierr=6
                    return
                endif                    
            enddo
        enddo       
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        ierr=0
        

    end subroutine ReadConfigurationFile
    
    !##############################################################################################################
    
    subroutine GetDimensions(NComponents,NPhaseMax,NAuxArray) bind(c,name='GetDimensions')
    !DEC$ ATTRIBUTES DLLEXPORT ::GetDimensions
    
        integer,intent(out)::NComponents                    ! - number of components of the thermodynamic system
        integer,intent(out)::NPhaseMax                      ! - maximum number of phases that the mixture can split into    
        integer,intent(out)::NAuxArray                      ! - number of auxiliary arrays
        
        !---------------------------------------------------------------------------------------------------------- 
        NComponents=nc
        NPhaseMax=npmax 
        NAuxArray=na
    
    end subroutine GetDimensions
    
    !##############################################################################################################
    
    subroutine GetGlobalParameters(CMPnames,MolWeights,PHnames,AuxArrNames,AuxArrDimens,opt) bind(c,name='GetGlobalParameters')
    !DEC$ ATTRIBUTES DLLEXPORT ::GetGlobalParameters

        character::CMPnames(3,nc)                           ! - names of the components
        real(8)::MolWeights(nc)                             ! - molar weights
        character::PHnames(3,npmax)                         ! - names of the phases
        character::AuxArrNames(8,na)                        ! - names of auxiliary quantities. Each name must begin with '#'
        character::AuxArrDimens(8,na)                       ! - dimensions of auxiliary quantities
        integer(1),intent(out)::opt(8)                      ! - Options (an option #i is on if opt(i)=1. Otherwise if opt(i)=0 the option is off.
                                                            !   if opt(1)=1, then the shared library must export the relative permeability
                                                            !   if opt(1)=0, then the relative permeability is specified in the MUFITS RUN-file.
                                                            !   if opt(2)=1, then the shared library must export the relative phae pressure 
                                                            !   i.e. the capillary pressure
                                                            !   if opt(2)=0, then the capillary pressure is specified in the MUFITS RUN-file.
        !----------------------------------------------------------------------------------------------------------
        
        integer::i

        do i=1,3
        CMPnames(i,1:nc)=CompNames(1:nc)(i:i)
        enddo

        MolWeights(1:nc)=MW(1:nc)
        
        do i=1,3
        PHnames(i,1)=PhaseNames(1)(i:i)
        enddo
        
        opt=0
        
        !opt(1)=1  !- uncomment this line if the relative permeability will be exported from the shared library
        !opt(2)=1  !- uncomment this line if the capillary pressure curve will be exported from the shared library

    end subroutine GetGlobalParameters
    
    !##############################################################################################################
    
    subroutine PhaseEquilibrium(Pres,Temp,z,NPhase,Props,PhaseID,AuxArr,Mode) bind(c,name='PhaseEquilibrium')
    !DEC$ ATTRIBUTES DLLEXPORT ::PhaseEquilibrium
    
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
                                  
        real(8)::Props(6+nc,npmax)    ! Parameters of the phases
                                      ! 1 - density [kg/m3]
                                      ! 2 - enthalpy [J/kg]
                                      ! 3 - viscosity [Pa*s]
                                      ! 4 - saturation [-]
                                      ! 5 - relative permeability [-] (must be specified only if opt(1)=1)
                                      ! 6 - relative phase pressure [Pa] (must be specified only if opt(2)=1)
                                      ! 7:6+nc - phase composition [-] 
        
        integer(1)::PhaseID(npmax)  ! phase ID
        
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
        do ic=1,nc
            Props(2,1)=Props(2,1)+HC(ic)*x(ic)
            Props(3,1)=Props(3,1)+Visc(ic)*x(ic)
        enddo
        Props(2,1)=Props(2,1)*Temp/M
        Props(4,1)=1d0
        
        !--- calculating composition ---
        Props(7:6+nc,1)=z
        
        PhaseID(1)=1   ! is the GAS phase

    end subroutine PhaseEquilibrium
    
end module USEREOS    