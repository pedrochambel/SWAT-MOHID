!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : SWAT Model
! PROJECT       : SWAT
! MODULE        : Swat
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : November 2004
! REVISION      : Pedro Chambel Leitao
! DESCRIPTION   : Module to introduce new options to swat2000 model
!
!------------------------------------------------------------------------------


Module ModuleSwat

    use ModuleGlobalData
    use ModuleTime
    use ModuleTimeSerie
    use ModuleEnterData
	use ModuleHDF5

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartSwat
    private ::      AllocateInstance

    !Selector
    public  :: GetSwatPointer
    public  :: GetSwatInteger
    public  :: UnGetSwat
                     
    
    !Modifier
    public  :: ModifySwat

    !Destructor
    public  :: KillSwat                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjSwat 
    
    !Interfaces----------------------------------------------------------------
    private :: UnGetSwat3D_I
    private :: UnGetSwat3D_R8
    interface  UnGetSwat
        module procedure UnGetSwat3D_I
        module procedure UnGetSwat3D_R8
    end interface  UnGetSwat

    !Types---------------------------------------------------------------------
    private :: T_SwatOptions
    type       T_SwatOptions
        logical                                     :: SwatWithChanges
        logical                                     :: MeteoOutputOn
        logical                                     :: PlantOutputOn
        logical                                     :: OperationalOutputOn
        logical                                     :: NitrateOutputOn
        logical                                     :: PhosphorusOutputOn
		logical                                     :: WaterBalanceOn
        logical                                     :: SubBasinOn
        logical                                     :: ReachesOn
        logical                                     :: GlobalWaterOn
        logical                                     :: GlobalWaterCumulative
        logical                                     :: GlobalWaterM3
        logical                                     :: GlobalNutrientsOn
        logical                                     :: GlobalNutrientsCumulative
        logical                                     :: BasinAverage
        logical                                     :: UsleOutputOn
        logical                                     :: ChangeMGT
		logical										:: LandUse
		logical										:: NutrientPerHRU

    end type  T_SwatOptions

    type T_OutPut
        type (T_Time), dimension(:), pointer        :: OutTime
        integer                                     :: NextOutPut
        logical                                     :: Yes = .false.
    end type T_OutPut

    private :: T_Plant
    type       T_Plant
        real, dimension(:,:), pointer               :: Biomass                !bio_ms
        real, dimension(:,:), pointer               :: NitrogenFraction       !cnb
        real, dimension(:,:), pointer               :: NitrogenYield          !cnyld
        real, dimension(:,:), pointer               :: FractionOfHeatUnits    !phuacc
        real, dimension(:,:), pointer               :: LAI                    !alai
        real, dimension(:,:), pointer               :: FractionGrowth_StressP     !strsp
        real, dimension(:,:), pointer               :: FractionGrowth_StressN     !strsn
        real, dimension(:,:), pointer               :: FractionGrowth_StressTemp  !strstmp   
        real, dimension(:,:), pointer               :: FractionGrowth_StressWater !strsw 
    end type  T_Plant

    private :: T_RateFlux
    type       T_RateFlux
        real, dimension(:,:), pointer               :: ActiveOrgNToNitrate  !hmntl
        real, dimension(:,:), pointer               :: FreshOrgNToNitrate   !rmn2tl
        real, dimension(:,:), pointer               :: Denitrification      !wdntl
    end type  T_RateFlux

    private :: T_Nitrate
    type       T_Nitrate
        real, dimension(:,:), pointer               :: LateralFlow          !latno3
        real, dimension(:,:), pointer               :: Rain                 !no3pcp
        real, dimension(:,:), pointer               :: PlantUptake          !nplnt
        real, dimension(:,:), pointer               :: Percolation          !percn
        real, dimension(:,:), pointer               :: RunOff               !surqno3
    end type  T_Nitrate

    private :: T_OrganicN
    type       T_OrganicN
        real, dimension(:), pointer                 :: Dissolved            !
        real, dimension(:), pointer                 :: Particulated         !
    end type  T_OrganicN

    private :: T_NTransportFlux
    type       T_NTransportFlux
        type (T_Nitrate)                            :: Nitrate
        real, dimension(:,:), pointer               :: AutoFertilization    !auton
        real, dimension(:,:), pointer               :: Fertilization        !fertn
        real, dimension(:,:), pointer               :: Grazing              !grazn
        real, dimension(:,:), pointer               :: Volatilization       !rvol
        real, dimension(:,:), pointer               :: Nitrification        !rnit
		real, dimension(:,:), pointer               :: RunOffOrganicN       !sedorgn
    end type  T_NTransportFlux

    private :: T_PTransportFlux
    type       T_PTransportFlux
        real, dimension(:,:), pointer               ::  RunOffActiveMineral !sedminpa(:)   !active mineral phosphorus sorbed to sediment in surface runoff in HRU for day
        real, dimension(:,:), pointer               ::  RunOffStableMineral !sedminps(:)   !stable mineral phosphorus sorbed to sediment in surface runoff in HRU for day
        real, dimension(:,:), pointer               ::  RunOffOrganicMatter !sedorgp(:)    !organic phosphorus in surface runoff in HRU for the day
        real, dimension(:,:), pointer               ::  AutoFertilization   !autop         !amount of phosphorus applied in auto-fert application
        real, dimension(:,:), pointer               ::  Fertilization       !fertp         !total amount of phosphorus applied to soil in HRU on day
        real, dimension(:,:), pointer               ::  Grazing             !grazp         !total amount of nitrogen applied to soil during grazing in HRU on day
        real, dimension(:,:), pointer               ::  PlantUptake         !pplnt
    end type  T_PTransportFlux

    private :: T_Pools_Nitrogen
    type       T_Pools_Nitrogen
        real, dimension(:,:), pointer               :: Nitrate              !sol_no3
        real, dimension(:,:), pointer               :: FreshOrganicMatter   !sol_fon
        real, dimension(:,:), pointer               :: ActiveOrganicMatter  !sol_aorgn
        real, dimension(:,:), pointer               :: StableOrganicMatter  !sol_orgn
        real, dimension(:,:), pointer               :: SoilResidue          !sol_rsd
        real, dimension(:,:), pointer               :: Ammonium             !sol_nh3
    end type  T_Pools_Nitrogen

    private :: T_Pools_Phosphorus
    type       T_Pools_Phosphorus
        real, dimension(:,:), pointer               ::  ActiveMineral !sol_actp(:,:) !phosphorus stored in the active mineral phosphorus pool
        real, dimension(:,:), pointer               ::  FreshOrganicMatter !sol_fop(:,:)  !phosphorus stored in the fresh organic (residue) pool
        real, dimension(:,:), pointer               ::  OrganicMatter !sol_orgp(:,:) !phosphorus stored in the organic P pool
        real, dimension(:,:), pointer               ::  StableMineral !sol_stap(:,:) !phosphorus in the soil layer stored in the stable mineral phosphorus pool
    end type  T_Pools_Phosphorus

    private :: T_WaterFluxes
    type       T_WaterFluxes
        real, dimension(:), pointer                 :: GlobalToChannel      !wwy
        real, dimension(:), pointer                 :: GroundWaterToChannel !wgwq
        real, dimension(:), pointer                 :: RunOffToChannel      !wqd
        real, dimension(:), pointer                 :: LateralFlowToChannel !wlatq
        real, dimension(:), pointer                 :: LateralAndRunOffFlowToChannel !wlatq+wqd
        real, dimension(:), pointer                 :: Precipitation            !sub_subp
        real, dimension(:), pointer                 :: ActualEvapotranspiration !sub_etday
    end type  T_WaterFluxes

    private :: T_Meteo
    type       T_Meteo
        real, dimension(:,:), pointer               :: AverageTemperature           !tmpav (j)
        real, dimension(:,:), pointer               :: Precipitation                !subp (j)
        real, dimension(:,:), pointer               :: PotencialEvapoTranspiration  !pet_day
        real, dimension(:,:), pointer               :: CulturalTranspiration        !ep_max
        real, dimension(:,:), pointer               :: ActualEvapoTranspiration     !etday
        real, dimension(:,:), pointer               :: MinimumTemperature           !tmn(j)
        real, dimension(:,:), pointer               :: MaximumTemperature           !tmx(j)
        real, dimension(:,:), pointer               :: WindSpeed                    !u10(j)
        real, dimension(:,:), pointer               :: RelativeHumidity             !rhd(j)
        real, dimension(:,:), pointer               :: SolarRadiation               !hru_ra(j)
        real, dimension(:,:), pointer               :: Irrigation                   !aird(j)
        real, dimension(:,:), pointer               :: PrecipitationAsSnow          !snofall
        real, dimension(:,:), pointer               :: SnowMelt                     !snomlt
        real, dimension(:,:), pointer               :: SnowPack_mmH2O               !sno_hru
        real, dimension(:,:), pointer               :: SoilTemperature              !sol_tmp(2,j)
        real, dimension(:,:), pointer               :: SoilWater                    !sol_sw(j)
        real, dimension(:,:), pointer               :: AquiferWater                 !shallst(j)
    end type  T_Meteo

    private :: T_WatBal
    type       T_WatBal
		! pools of water (dstor)
        real, dimension(:,:), pointer               :: SnowPack_mmH2O						! sno_hru(j)            
        real, dimension(:,:), pointer               :: SnowPack_mmH2O_Old					! snoprev            
        real, dimension(:,:), pointer               :: SoilWater							! sol_sw(j)            
        real, dimension(:,:), pointer               :: SoilWaterOld							! swprev            
        real, dimension(:,:), pointer               :: AquiferWater							! shallst(j)            
        real, dimension(:,:), pointer               :: AquiferWaterOld						! shallstp            
        real, dimension(:,:), pointer               :: DeepAquifer							! deepst(j) 
        real, dimension(:,:), pointer               :: DeepAquiferOld						! deepstp            
        real, dimension(:,:), pointer               :: SurfaceWaterColumn					! surf_bs(1,j)           
        real, dimension(:,:), pointer               :: SurfaceWaterColumnOld				! bsprev            
        real, dimension(:,:), pointer               :: LateralFlowStorage					! bss(1,j)            
        real, dimension(:,:), pointer               :: LateralFlowStorageOld				! bssprev           
		! water fluxes (h2oloss)
        real, dimension(:,:), pointer               :: Precipitation						! subp(j)              
        real, dimension(:,:), pointer               :: SnowEvaporation						! snoev                  
        real, dimension(:,:), pointer               :: RunOffFlowToChannel					! qday                  
        real, dimension(:,:), pointer               :: LateralFlowToChannel					! latq(j)        
        real, dimension(:,:), pointer               :: ActualEvapoTranspiration				! etday          
        real, dimension(:,:), pointer               :: GrowndWaterFlowToChannel				! gw_q(j)            
        real, dimension(:,:), pointer               :: GrowndWaterFlowToAtmosphere			! revapday      
        real, dimension(:,:), pointer               :: PondWaterLoss						! twlpnd              
        real, dimension(:,:), pointer               :: WetlandWaterLoss						! twlwet                
        real, dimension(:,:), pointer               :: Irrigation							! aird(j)          
        real, dimension(:,:), pointer               :: RechargeFlowToDeepAndShalowAquifer	! rchrg(j)        
        real, dimension(:,:), pointer               :: DrainageTileFlowToChannel			! qtile         
        real, dimension(:,:), pointer               :: SoilFlowToRechargeFlow   	! sepbtm(j)  
		real, dimension(:,:), pointer               :: WaterBalanceError         
    end type  T_WatBal 


    private :: T_SubBasin
    type       T_SubBasin
        type (T_WaterFluxes)                        :: WaterFluxes
        type (T_OrganicN)                           :: OrganicN
        real, dimension(:), pointer                 :: OrganicP  !
        real, dimension(:), pointer                 :: MineralP     !
        real, dimension(:), pointer                 :: Sediments_concentration     !
        real, dimension(:), pointer                 :: Sediments_load     ! ton/ha
        real, dimension(:), pointer                 :: NitrateInRunOff      !wyno3
        integer, dimension(:  ),  pointer           :: WaterPoints1D
        real,dimension (:), pointer                 :: SubBasinFraction    !sub_fr - fraction of watershed in subbasin
    end type  T_SubBasin

    private :: T_SubBasin2
    type       T_SubBasin2
        type (T_WaterFluxes)                        :: WaterFluxes
        type (T_OrganicN)                           :: RunOffOrganicN
        real, dimension(:), pointer                 :: RunOffSolubleP        ! soluble P in surface runoff on day in subbasin
        real, dimension(:), pointer                 :: RunOffSedimentStableP ! stable mineral P attached to sediment removed in surface runoff
        real, dimension(:), pointer                 :: RunOffSedimentActiveP ! active mineral P attached to sediment removed in surface runoff
        real, dimension(:), pointer                 :: GrowndWaterSolubleP ! 
        real, dimension(:), pointer                 :: Sediments            ! ton/ha
        real, dimension(:), pointer                 :: NitrateInRunOff      ! sub_no3
        real, dimension(:), pointer                 :: NitrateInLateralFlow ! sub_latno3
        real, dimension(:), pointer                 :: NitrateInGrowndWaterFlow ! sub_gwno3
        integer, dimension(:  ),  pointer           :: WaterPoints1D
        real,dimension (:), pointer                 :: SubBasinFraction    !sub_fr - fraction of watershed in subbasin
    end type  T_SubBasin2

    private :: T_Reaches
    type       T_Reaches
        real, dimension(:), pointer                 :: FlowOut                  
        real, dimension(:), pointer                 :: algal_biomass           !    algae(:)    |mg alg/L   
        real, dimension(:), pointer                 :: ammonia                 !    ammonian(:) |mg N/L     
        real, dimension(:), pointer                 :: BOD                     !    bod(:)      |mg O2/L    
        real, dimension(:), pointer                 :: chlorophyll_a           !    chlora(:)   |mg chl-a/L 
        real, dimension(:), pointer                 :: dissolved_phosphorus    !    disolvp(:)  |mg P/L     
        real, dimension(:), pointer                 :: dissolved_oxygen        !    rch_dox(:)  |mg O2/L    
        real, dimension(:), pointer                 :: nitrate                 !    nitraten(:) |mg N/L     
        real, dimension(:), pointer                 :: nitrite                 !    nitriten(:) |mg N/L     
        real, dimension(:), pointer                 :: organic_nitrogen        !    organicn(:) |mg N/L     
        real, dimension(:), pointer                 :: organic_phosphorus      !    organicp(:) |mg P/L     
        real, dimension(:), pointer                 :: saturation_oxygen       !    soxy        |mg O2/L 
        real, dimension(:), pointer                 :: temperature             !    wtmp        |mg O2/L 
        real, dimension(:), pointer                 :: Sediments               !    sedyld      |metric tons 
    end type  T_Reaches

    private :: T_Usle
    type       T_Usle
        real, dimension(:,:), pointer               :: Musle_Sed            !sedyld - metric tons
        real, dimension(:,:), pointer               :: Usle_Sed             !usle   - metric tons/ha
        real, dimension(:,:), pointer               :: C_K_P_LS_rock        !cklsp  
        real, dimension(:,:), pointer               :: MuslePeakRunoffRate  !peakr  - m^3/s
        real, dimension(:,:), pointer               :: UsleRainfallIndex    !usle_ei  - 100(ft-tn in)/(acre-hr)
        real, dimension(:,:), pointer               :: K_LS_P_rock          !usle_mult  - none
    end type  T_Usle

    private :: T_GW
    type       T_GW
        real, dimension(:,:), pointer               :: NitrateConcentration ! no3gw(j)/gw_q(j)
!        real, dimension(:,:), pointer               :: AquiferWater							! shallst(j)            
!        real, dimension(:,:), pointer               :: AquiferWaterOld						! shallstp            
!        real, dimension(:,:), pointer               :: DeepAquifer							! deepst(j) 
!        real, dimension(:,:), pointer               :: DeepAquiferOld						! deepstp            
!        real, dimension(:,:), pointer               :: GrowndWaterFlowToChannel				! gw_q(j)            
!        real, dimension(:,:), pointer               :: GrowndWaterFlowToAtmosphere			! revapday      
!        real, dimension(:,:), pointer               :: Irrigation							! aird(j)          
!        real, dimension(:,:), pointer               :: RechargeFlowToDeepAndShalowAquifer	! rchrg(j)        
    end type  T_GW

    private :: T_HRU
    type       T_HRU
        type (T_GW)                                 :: GW
        type (T_Usle)                               :: Usle
        type (T_Meteo)                              :: Meteo
        type (T_Plant)                              :: Plant
        type (T_Pools_Nitrogen)                     :: PoolsN
        type (T_Pools_Phosphorus)                   :: PoolsP
        type (T_RateFlux)                           :: RateFlux
        type (T_NTransportFlux)                     :: NTransportFlux
        type (T_PTransportFlux)                     :: PTransportFlux
		type (T_WatBal)								:: WatBal
        real, dimension(:,:,:),pointer              :: HRU_output_values
        integer, dimension(:,:),  pointer           :: WaterPoints2D
        real                                        :: DenitrificationThreshold
        real                                        :: LateralFlowTravelTime  !lat_ttime [days]
 		integer										:: NumberOfLandUseOutputs
		character(len=4),dimension(:), pointer      :: LandUse
		real, dimension(:,:),pointer                :: LandUseNutrient
		real,   dimension(:),pointer                :: LandUseAreaHectare
        logical                                     :: AuxLogicalLandUse =.true.
    end type  T_HRU

    private :: T_SubBasin_Annual_Average
    type       T_SubBasin_Annual_Average
        type(T_SubBasin)                            :: SubBasin
        type(T_SubBasin2)                           :: SubBasin2
        type(T_HRU)                                 :: HRU
    end type  T_SubBasin_Annual_Average

    !Average all HRU in Sub-basins
    private :: T_SubBasin_Average
    type       T_SubBasin_Average
        type(T_SubBasin)                            :: SubBasin
        type(T_SubBasin2)                           :: SubBasin2
        type(T_HRU)                                 :: HRU
    end type  T_SubBasin_Average

    private :: T_Swat
    type       T_Swat
        integer                                     :: InstanceID
        type (T_Size3D)                             :: Size, WorkSize
        real(8), dimension(:, :, :),  pointer       :: Matrix
        type(T_Swat), pointer                       :: Next
        integer                                     :: ObjEnterData = 0
        integer                                     :: ObjMeteoTimeSerie = 0
        integer                                     :: ObjPlantTimeSerie = 0
        integer                                     :: ObjOperationalTimeSerie = 0
        integer                                     :: ObjUsleTimeSerie = 0
        integer                                     :: ObjNitrateTimeSerie = 0
        integer                                     :: ObjPhosphorusTimeSerie = 0
        integer                                     :: ObjWaterBalanceTimeSerie = 0
        integer                                     :: ObjSubBasinTimeSerie = 0
        integer                                     :: ObjReachesTimeSerie = 0
        integer                                     :: ObjGlobalWaterTimeSerie = 0
        integer                                     :: ObjGlobalNutrientTimeSerie = 0
        integer                                     :: ObjBasinAverage = 0
		integer                                     :: ObjLandUseTimeSerie        = 0
        integer                                     :: ObjTime
        integer                                     :: ObjHDF5_Reaches               = 0
        integer                                     :: ObjHDF5_SubBasin              = 0
        integer                                     :: BeginYear        !iyr
        integer                                     :: EndYear          !nbyr
        integer                                     :: BeginJulianDay   !idaf
        integer                                     :: EndJulianDay     !idal
        integer                                     :: MaxNumberHRU     !hrutot
        integer                                     :: MaxNumberSubBasin !mbb
        real                                        :: TotalBasinArea    !da_km
        real                                        :: NumberYearsSimulation
        real,dimension(:,:), pointer                :: HRUBasinFraction    !hru_fr - fraction of watershed in HRU
        real,dimension(:,:), pointer                :: dummy
        real,dimension(:,:), pointer                :: pcp     
        real,dimension(:,:), pointer                :: Flow    
        real,dimension(:,:), pointer                :: BaseFlow
        real,dimension(:,:), pointer                :: ET      
        real,dimension(:,:), pointer                :: Round   
        real,dimension(:,:), pointer                :: CountMonthDays
		integer,dimension(:,:), pointer             :: WaterStress05
		integer,dimension(:,:), pointer             :: WaterStress25
		integer,dimension(:,:), pointer             :: WaterStress50
        real,dimension(:),   pointer                :: RiverFlow

        real,dimension(:,:), pointer                :: N_Run_Off
        real,dimension(:,:), pointer                :: N_Percolated
        real,dimension(:,:), pointer                :: P_Run_Off	
        real,dimension(:,:), pointer                :: Sediment	

        real, dimension(:), pointer                 :: GlobalWaterFluxes
        real, dimension(:), pointer                 :: GlobalNutrients
        real, dimension(:), pointer                 :: Basin_Average !Average of all watershed
        real                                        :: Month
        real                                        :: Month_old
        type (T_Reaches)                            :: Reaches
        type (T_SubBasin)                           :: SubBasin
        type (T_SubBasin2)                          :: SubBasin2
        type (T_HRU)                                :: HRU
        type (T_SubBasin_Annual_Average)            :: SubBasin_Annual_Average
        type (T_SubBasin_Average)                   :: SubBasin_Average
        type (T_SwatOptions)                        :: SwatOptions
        type (T_OutPut)                             :: OutPut
		type (T_Time)                               :: CurrentTime
		type (T_Time)                               :: EndTime
		type (T_Time)                               :: BeginTime
    end type  T_Swat

    !Global Module Variables
    type (T_Swat), pointer                         :: FirstObjSwat
    type (T_Swat), pointer                         :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartSwat (ObjSwatID,                    &
                            iyr,                        &
                            nbyr,                       &
                            idaf,                       &
                            idal,                       &
                            MaxHruInSubBasin,           &
                            mbb,                        &
                            da_km,                      &
                            SwatWithChanges,            &
                            DenitrificationThreshold,   &
                            LateralFlowTravelTime,      &
                            ChangeMGT,                  &
                            STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjSwatID 
        integer, optional, intent(OUT)                  :: STAT     
        logical, intent(OUT)                            :: SwatWithChanges
        real                                            :: DenitrificationThreshold  
        real                                            :: LateralFlowTravelTime 
        logical                                         :: ChangeMGT

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_
        integer                                         :: STAT_CALL
        integer                                         :: iflag
        integer                                         :: iyr, nbyr, idaf, idal
        integer                                         :: MaxHruInSubBasin, mbb
        real                                            :: da_km

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_
 
        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mSWAT_)) then
            nullify (FirstObjSwat)
            call RegisterModule (mSWAT_) 
        endif

        call Ready(ObjSwatID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call StartupMohid('SWAT')

            call AllocateInstance

            call ReadKeywords

            SwatWithChanges = Me%SwatOptions%SwatWithChanges
            ChangeMGT       = Me%SwatOptions%ChangeMGT

            DenitrificationThreshold = Me%HRU%DenitrificationThreshold
            LateralFlowTravelTime    = Me%HRU%LateralFlowTravelTime

            !Returns ID
            ObjSwatID          = Me%InstanceID

            if (Me%SwatOptions%SwatWithChanges) then
                Me%BeginYear        = iyr
                Me%EndYear          = nbyr
                Me%BeginJulianDay   = idaf
                Me%EndJulianDay     = idal
                Me%MaxNumberHRU     = MaxHruInSubBasin
                Me%MaxNumberSubBasin = mbb
                Me%TotalBasinArea    = da_km

                call StartTime

				!Gets Output Time 
				call GetOutPutTime(Me%ObjEnterData,                                                 &
								   CurrentTime = Me%CurrentTime,                                  &
								   EndTime     = Me%EndTime,                                     &
								   keyword     = 'OUTPUT_TIME',                                  &
								   SearchType  = FromFile,                                       &
								   OutPutsTime = Me%OutPut%OutTime,                              &
								   OutPutsOn   = Me%OutPut%Yes,                                  &
								   STAT        = STAT_CALL)
				Me%OutPut%NextOutPut = 1  

                call AllocateVariables

                call Construct_Time_Serie

				if(Me%OutPut%Yes) call ConstructHDF5Output

				if(Me%SwatOptions%LandUse) call LandUseList

				open(unit=111,file='error.txt')
							write(111,999) "hru",",",&
							"sub",",",&
							"Year",",", &
							"Month",",",&
							"Day",",",&
							"pcp",",", &
							"Snow",",",&
							"RunOff",",",&
							"Lateral",",", &
							"ET"	,",", &			
							"GW"	,",", &		
							"GW_ET"	,",", &		
							"Pond"	,",", &	
							"Wetl"	,",", &		
							"Irri"	,",", &	
							"RchDpSh"	,",", &
							"Tile"	,",", &
							"SlRec"	,",", &
							"SW"	,",", &
							"SWOld" ,",", &
							"AqW" ,",", &						
							"AqWOld" ,",", &					
							"DpAq" ,",", &						
							"DpAqOld"	,",", &				
							"SWC"	,",", &			
							"SWCOld"	,",", &			
							"LatS"	,",", &			
							"LatSOld"	,",", &	
							"h2oloss"	,",", &		
							"dstor"
 999    format (a5,a1,a5,a1,a5,a1,a5,a1,a5,a1,&
			    a5,a1,a5,a1,a5,a1,a5,a1,a5,a1,&
			    a5,a1,a5,a1,a5,a1,a5,a1,a5,a1,&
			    a5,a1,a5,a1,a5,a1,a5,a1,a5,a1,&
			    a5,a1,a5,a1,a5,a1,a5,a1,a5,a1,&
			    a5,a1,a5,a1,a5,a1,a5,a1,a5,a1)

				open(unit=222,file='error2.txt')
							write(222,9999) "hru,",&
							"sub,",&
							"Year,", &
							"Month,",&
							"Hid_Year,",&
							"pcp,", &
							"Flow,", &			
							"BaseFlow,", &		
							"ET,"	, &		
							"Round,", &	
							"Count,", &	
							"DeepAqui,",	&
							"SW,",	&			
							"RivFlow,",	&			
							"Str05,"	,	&		
							"Str25,"	,	&		
							"Str50,"			

 9999    format (17a9)


				open(unit=333,file='N_P_Sed.txt')
							write(333,99999) "hru,",&
							"sub,",&
							"Year,", &
							"Month,",&
							"Hid_Year,",&
							"N_RunOff,", &
							"N_Percol,", &			
							"P_RunOff,", &	
							"Sediment"	

 99999    format (9a9)


            endif


            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleSwat - StartSwat - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartSwat
    !--------------------------------------------------------------------------

    subroutine ReadKeywords

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: iflag
        logical                                         :: dummy
        !------------------------------------------------------------------------

        call ConstructEnterData(EnterDataID = Me%ObjEnterData, &
                                FileName    = '..\res\OutPuts.dat', &
                                STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'StartSwat - ERR01'

        call GetData(Me%SwatOptions%SwatWithChanges,                                      &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='SWAT_WITH_CHANGES',                       &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR02'

        call GetData(Me%SwatOptions%MeteoOutputOn,                                          &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='METEO_OUT',                              &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR03'

        call GetData(Me%SwatOptions%PlantOutputOn,                                          &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='PLANTS_OUT',                              &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR03'

        call GetData(Me%SwatOptions%OperationalOutputOn,                             &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='OPERATIONAL',                   &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR03c'
	!	if (Me%SwatOptions%OperationalOutputOn) then
	!		if (.not.Me%SwatOptions%PlantOutputOn .or. .not.Me%SwatOptions%MeteoOutputOn) then
	!			stop 'METEO_OUT must be ON and PLANTS_OUT must be ON'
	!		endif
	!	end if


        call GetData(Me%SwatOptions%NitrateOutputOn,                                          &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='NITRATE_OUT',                              &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR04'

        call GetData(Me%SwatOptions%PhosphorusOutputOn,                         &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='PHOSPHORUS_OUT',                          &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR04b'

        call GetData(Me%SwatOptions%WaterBalanceOn,                             &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='HRU_WATER_BALANCE_OUT',                   &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR04b'

        call GetData(dummy,                                          &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='SUBBASIN_DISCHARGES_OUT',                              &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (iflag==1)stop 'Replace SUBBASIN_DISCHARGES_OUT by SUBBASIN_OUT'

        call GetData(Me%SwatOptions%SubBasinOn,                                          &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='SUBBASIN_OUT',                              &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR05'

        call GetData(Me%SwatOptions%ReachesOn,                                          &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='REACHES_OUT',                              &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR06'

        call GetData(Me%SwatOptions%GlobalWaterOn,                              &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='GLOBAL_WATER_OUT',                        &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR07'

        if (Me%SwatOptions%GlobalWaterOn)  then
            call GetData(Me%SwatOptions%GlobalWaterCumulative,                              &
                         Me%ObjEnterData,                                           &
                         flag           = iflag,                                    &
                         SearchType     = FromFile,                                 &
                         keyword        ='GLOBAL_WATER_CUMULATIVE',                        &
                         Default        =.FALSE.,                                   &
                         ClientModule   ='ModuleSwat',                              &
                         STAT           = STAT_CALL)             
            if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR07'

            call GetData(Me%SwatOptions%GlobalWaterM3,                              &
                         Me%ObjEnterData,                                           &
                         flag           = iflag,                                    &
                         SearchType     = FromFile,                                 &
                         keyword        ='GLOBAL_WATER_M3',                        &
                         Default        =.FALSE.,                                   &
                         ClientModule   ='ModuleSwat',                              &
                         STAT           = STAT_CALL)             
            if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR07'
        endif

        call GetData(Me%SwatOptions%GlobalNutrientsOn,                                          &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='GLOBAL_NUTRIENTS_OUT',                    &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR08'

        if (Me%SwatOptions%GlobalNutrientsOn)  then
            call GetData(Me%SwatOptions%GlobalNutrientsCumulative,              &
                         Me%ObjEnterData,                                       &
                         flag           = iflag,                                &
                         SearchType     = FromFile,                             &
                         keyword        ='GLOBAL_NUTRIENTS_CUMULATIVE',         &
                         Default        =.FALSE.,                               &
                         ClientModule   ='ModuleSwat',                          &
                         STAT           = STAT_CALL)             
            if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR07'
        endif

        call GetData(Me%SwatOptions%BasinAverage,                              &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='BASIN_AVERAGE',                           &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR07'

        if (Me%SwatOptions%BasinAverage) then
            if(.not.Me%SwatOptions%SubBasinOn) then
                stop 'StartSwat - ERR07b - please turn on SubBasins output: SUBBASIN_OUT: 1'
            endif
        endif

        call GetData(Me%SwatOptions%NutrientPerHRU,                        &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='NUTRIENT_PER_HRU',                   &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR07'

        if (Me%SwatOptions%NutrientPerHRU) then
            if((.not.Me%SwatOptions%PhosphorusOutputOn).and. (.not.Me%SwatOptions%NitrateOutputOn)) then
                stop 'StartSwat - ERR07c - please turn on NITRATE_OUT: 1 and PHOSPHORUS_OUT: 1'
            endif
        endif

        call GetData(Me%HRU%DenitrificationThreshold,                           &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='DENITRIFICATION_THRESHOLD',               &
                     Default        = 0.99,                                     &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR08'
        
        if (Me%HRU%DenitrificationThreshold /= 0.95) then
            write(*,*) 'SWAT value for DenitrificationThreshold is 0.95' 
            write(*,*) 'You are using a value of ', Me%HRU%DenitrificationThreshold 
            write(*,*) 'To change this value please enter keyword "DENITRIFICATION_THRESHOLD"'
            write(*,*) 'in file OutPuts.dat'
            write(*,*) '---------------######----------------'
        endif

        call GetData(Me%HRU%LateralFlowTravelTime,                              &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='LATERAL_FLOW_TRAVEL_TIME',                &
                     Default        = -99999.,                                  &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR08'

        if (Me%HRU%LateralFlowTravelTime > 0.0) then
            write(*,*) 'You are using a Lateral Flow Travel Time of ', Me%HRU%LateralFlowTravelTime 
            write(*,*) 'To use SWAT default values remove keyword "LATERAL_FLOW_TRAVEL_TIME" '
            write(*,*) 'in file OutPuts.dat'
            write(*,*) '---------------######----------------'
        
            Me%HRU%LateralFlowTravelTime = 1. - Exp(-1./Me%HRU%LateralFlowTravelTime) !equation obtained from subroutine hydroinit.f

        endif
             
        call GetData(Me%SwatOptions%UsleOutputOn,                               &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='USLE_OUT',                                &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR03'

        call GetData(Me%SwatOptions%ChangeMGT,                                  &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='CHANGE_MGT',                              &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR03'

        if (Me%SwatOptions%ChangeMGT) then
            write(*,*) 'You are replacing "kill/end of growing season" operation' 
            write(*,*) 'by "harvest and kill" operation' 
            write(*,*) 'This is only to be used with some specific managment files'
            write(*,*) 'Please check in mohid source safe, other tools\SWAT\Management files'
            write(*,*) 'To use SWAT default values remove keyword "CHANGE_MGT" '
            write(*,*) 'in file OutPuts.dat'
            write(*,*) '---------------######----------------'
        endif

        call GetData(Me%SwatOptions%LandUse,                                    &
                     Me%ObjEnterData,                                           &
                     flag           = iflag,                                    &
                     SearchType     = FromFile,                                 &
                     keyword        ='LAND_USE',                                &
                     Default        =.FALSE.,                                   &
                     ClientModule   ='ModuleSwat',                              &
                     STAT           = STAT_CALL)             
        if (STAT_CALL .NE. SUCCESS_)stop 'StartSwat - ERR33'



    end subroutine ReadKeywords
    !--------------------------------------------------------------------------
    subroutine StartTime

        real                                            :: DT = 86400.
        logical                                         :: VariableDT = .FALSE.
        integer                                         :: aux

        call SetDate (Me%BeginTime,Me%BeginYear,1,1,0,0,0 )

        Me%BeginTime = Me%BeginTime + (Me%BeginJulianDay-1)*86400

        aux = Me%BeginYear + Me%EndYear - 1

        call SetDate (Me%EndTime,aux,1,1,0,0,0 )

        Me%EndTime = Me%EndTime + (Me%EndJulianDay-1)*86400

        call StartComputeTime (TimeID     = Me%ObjTime,  &
                               BeginTime  = Me%BeginTime,&
                               EndTime    = Me%EndTime , &
                               DT         = DT,          &
                               VariableDT = VariableDT)

		call GetComputeCurrentTime(Me%ObjTime, Me%CurrentTime)

        Me%NumberYearsSimulation=(Me%EndTime - Me%BeginTime)/86400./365.25


    end subroutine StartTime
    !--------------------------------------------------------------------------

    subroutine AllocateVariables

        !Local-------------------------------------------------------------------
        integer                                         :: MaxNumberHRU
        integer                                         :: MaxNumberSubBasin
        !------------------------------------------------------------------------

        allocate(Me%HRUBasinFraction (1:Me%MaxNumberHRU,1:Me%MaxNumberSubBasin))

        allocate(Me%SubBasin%SubBasinFraction(1:Me%MaxNumberSubBasin))
        allocate(Me%SubBasin2%SubBasinFraction(1:Me%MaxNumberSubBasin))

        MaxNumberHRU      = Me%MaxNumberHRU

        allocate(Me%HRU%WaterPoints2D                (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
        allocate(Me%SubBasin%WaterPoints1D           (1:Me%MaxNumberSubBasin))
        allocate(Me%SubBasin2%WaterPoints1D           (1:Me%MaxNumberSubBasin))
        Me%Month = 0.0
        Me%Month_old = 0.0

        if (Me%SwatOptions%PlantOutputOn  .or.&
            Me%SwatOptions%NitrateOutputOn.or.&
            Me%SwatOptions%MeteoOutputOn  .or.&
            Me%SwatOptions%UsleOutputOn   .or.&
            Me%SwatOptions%WaterBalanceOn .or.&
            Me%SwatOptions%OperationalOutputOn .or.&
            Me%SwatOptions%PhosphorusOutputOn) then


            if (Me%SwatOptions%MeteoOutputOn.or.Me%OutPut%Yes.or.Me%SwatOptions%OperationalOutputOn) then
                allocate(Me%HRU%Meteo%AverageTemperature         (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%Precipitation              (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%PotencialEvapoTranspiration(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%CulturalTranspiration      (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%ActualEvapoTranspiration   (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%MinimumTemperature         (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%MaximumTemperature         (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%WindSpeed                  (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%RelativeHumidity           (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%SolarRadiation             (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%Irrigation                 (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%PrecipitationAsSnow        (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%SnowMelt                   (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%SnowPack_mmH2O             (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%SoilTemperature            (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%SoilWater                  (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Meteo%AquiferWater               (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
            endif

            if (Me%SwatOptions%PlantOutputOn.or.Me%OutPut%Yes.or.Me%SwatOptions%OperationalOutputOn) then
                allocate(Me%HRU%Plant%Biomass            (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Plant%NitrogenFraction   (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Plant%NitrogenYield      (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Plant%FractionOfHeatUnits(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Plant%LAI                (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Plant%FractionGrowth_StressP     (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Plant%FractionGrowth_StressN     (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Plant%FractionGrowth_StressTemp  (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Plant%FractionGrowth_StressWater (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
            endif

            if (Me%SwatOptions%UsleOutputOn) then
                allocate(Me%HRU%Usle%Musle_Sed            (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Usle%Usle_Sed             (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Usle%C_K_P_LS_rock        (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Usle%MuslePeakRunoffRate  (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Usle%UsleRainfallIndex    (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%Usle%K_LS_P_rock          (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))   
                allocate(Me%SubBasin_Annual_Average%HRU%Usle%Musle_Sed(1:1,1:Me%MaxNumberSubBasin))
                Me%SubBasin_Annual_Average%HRU%Usle%Musle_Sed = 0.0
            endif

            if (Me%SwatOptions%NitrateOutputOn) then
                allocate(Me%HRU%RateFlux%ActiveOrgNToNitrate        (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%RateFlux%FreshOrgNToNitrate         (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%RateFlux%Denitrification            (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))

                allocate(Me%HRU%NTransportFlux%Nitrate%LateralFlow  (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%NTransportFlux%Nitrate%Rain         (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%NTransportFlux%Nitrate%PlantUptake  (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%NTransportFlux%Nitrate%Percolation  (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%NTransportFlux%Nitrate%RunOff       (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%NTransportFlux%RunOffOrganicN       (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))

                allocate(Me%HRU%NTransportFlux%AutoFertilization    (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%NTransportFlux%Fertilization        (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%NTransportFlux%Grazing              (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%NTransportFlux%Volatilization       (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))

                allocate(Me%HRU%PoolsN%Nitrate                       (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%PoolsN%FreshOrganicMatter            (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%PoolsN%ActiveOrganicMatter           (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%PoolsN%StableOrganicMatter           (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%PoolsN%SoilResidue                   (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                
				allocate(Me%HRU%GW%NitrateConcentration              (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))   

                allocate(Me%HRU%PoolsN%Ammonium                      (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))

                allocate(Me%HRU%NTransportFlux%Nitrification         (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))

!                allocate(Me%SubBasin_Annual_Average%HRU%NTransportFlux%Nitrate%Percolation (1:1,1:Me%MaxNumberSubBasin))
!                Me%SubBasin_Annual_Average%HRU%NTransportFlux%Nitrate%Percolation = 0.0
             endif

            if (Me%SwatOptions%PhosphorusOutputOn) then
                allocate(Me%HRU%PoolsP%ActiveMineral              (1:MaxNumberHRU,1:Me%MaxNumberSubBasin)) 
                allocate(Me%HRU%PoolsP%FreshOrganicMatter         (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%PoolsP%OrganicMatter              (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%PoolsP%StableMineral              (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%PTransportFlux%AutoFertilization  (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%PTransportFlux%Fertilization      (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%PTransportFlux%Grazing            (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%PTransportFlux%RunOffActiveMineral(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%PTransportFlux%RunOffStableMineral(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%PTransportFlux%RunOffOrganicMatter(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%PTransportFlux%PlantUptake        (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
            endif


            if (Me%SwatOptions%WaterBalanceOn) then
		! pools of water (dstor)
                allocate(Me%HRU%WatBal%SnowPack_mmH2O							(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%SnowPack_mmH2O_Old						(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%SoilWater								(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%SoilWaterOld							(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%AquiferWater							(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%AquiferWaterOld						(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%DeepAquifer							(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%DeepAquiferOld							(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%SurfaceWaterColumn						(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%SurfaceWaterColumnOld					(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%LateralFlowStorage						(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%LateralFlowStorageOld					(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
		! water fluxes (h2oloss)
                allocate(Me%HRU%WatBal%Precipitation							(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%SnowEvaporation						(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%RunOffFlowToChannel					(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%LateralFlowToChannel					(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%ActualEvapoTranspiration				(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%GrowndWaterFlowToChannel				(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%GrowndWaterFlowToAtmosphere			(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%PondWaterLoss							(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%WetlandWaterLoss						(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%Irrigation								(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%RechargeFlowToDeepAndShalowAquifer		(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%DrainageTileFlowToChannel				(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%HRU%WatBal%SoilFlowToRechargeFlow			(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
				allocate(Me%HRU%WatBal%WaterBalanceError                    (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%pcp      (1:MaxNumberHRU,1:Me%MaxNumberSubBasin)) 
                allocate(Me%Flow     (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%BaseFlow (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%ET       (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%Round    (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%CountMonthDays    (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%WaterStress05     (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%WaterStress25     (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%WaterStress50     (1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%RiverFlow(1:Me%MaxNumberSubBasin))
				Me%pcp      = 0.0
				Me%Flow     = 0.0
				Me%BaseFlow = 0.0
				Me%ET       = 0.0
				Me%Round    = 0.0
				Me%CountMonthDays= 0.0
				Me%RiverFlow     = 0.0
				Me%WaterStress05= 0
				Me%WaterStress25 = 0
				Me%WaterStress50 = 0
			end if

            if (Me%SwatOptions%NutrientPerHRU) then
                allocate(Me%N_Run_Off   (1:MaxNumberHRU,1:Me%MaxNumberSubBasin)) 
                allocate(Me%N_Percolated(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%P_Run_Off	(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
                allocate(Me%Sediment	(1:MaxNumberHRU,1:Me%MaxNumberSubBasin))
				Me%N_Run_Off     = 0.0
				Me%N_Percolated  = 0.0
				Me%P_Run_Off	 = 0.0
				Me%Sediment  	 = 0.0
			end if



        endif

        if (Me%SwatOptions%SubBasinOn.or.Me%SwatOptions%ReachesOn) then
 
            MaxNumberSubBasin = Me%MaxNumberSubBasin
            allocate(Me%dummy                    (1,1:MaxNumberSubBasin))
            allocate(Me%Reaches%FlowOut          (1:Me%MaxNumberSubBasin))
                
            if (Me%SwatOptions%SubBasinOn) then
                allocate(Me%SubBasin%WaterFluxes%GlobalToChannel      (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin%WaterFluxes%GroundWaterToChannel (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin%WaterFluxes%RunOffToChannel      (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin%WaterFluxes%LateralFlowToChannel (1:Me%MaxNumberSubBasin))

                allocate(Me%SubBasin%NitrateInRunOff                   (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin%OrganicN%Dissolved               (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin%OrganicN%Particulated            (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin%OrganicP                         (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin%MineralP                         (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin%Sediments_concentration          (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin%Sediments_load                   (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin%WaterFluxes%LateralAndRunOffFlowToChannel (1:Me%MaxNumberSubBasin))
               ! allocate(Me%SubBasin_Annual_Average%HRU%Meteo%SoilWater(1:1,1:Me%MaxNumberSubBasin))
            endif

            if (Me%SwatOptions%SubBasinOn) then
                allocate(Me%SubBasin2%WaterFluxes%GlobalToChannel      (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin2%WaterFluxes%GroundWaterToChannel (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin2%WaterFluxes%RunOffToChannel      (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin2%WaterFluxes%LateralFlowToChannel (1:Me%MaxNumberSubBasin))

                allocate(Me%SubBasin2%NitrateInRunOff                  (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin2%NitrateInLateralFlow             (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin2%NitrateInGrowndWaterFlow         (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin2%RunOffOrganicN%Dissolved         (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin2%RunOffOrganicN%Particulated      (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin2%RunOffSolubleP                   (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin2%RunOffSedimentStableP            (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin2%RunOffSedimentActiveP            (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin2%GrowndWaterSolubleP              (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin2%Sediments                        (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin2%WaterFluxes%LateralAndRunOffFlowToChannel (1:Me%MaxNumberSubBasin))

                allocate(Me%SubBasin2%WaterFluxes%Precipitation           (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin2%WaterFluxes%ActualEvapotranspiration(1:Me%MaxNumberSubBasin))

                allocate(Me%SubBasin_Annual_Average%SubBasin2%NitrateInRunOff            (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin_Annual_Average%SubBasin2%RunOffOrganicN%Dissolved   (1:Me%MaxNumberSubBasin))  
                allocate(Me%SubBasin_Annual_Average%SubBasin2%RunOffOrganicN%Particulated(1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin_Annual_Average%SubBasin2%NitrateInLateralFlow       (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin_Annual_Average%SubBasin2%NitrateInGrowndWaterFlow   (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin_Annual_Average%SubBasin2%RunOffSolubleP             (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin_Annual_Average%SubBasin2%RunOffSedimentActiveP      (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin_Annual_Average%SubBasin2%RunOffSedimentStableP      (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin_Annual_Average%SubBasin2%GrowndWaterSolubleP        (1:Me%MaxNumberSubBasin))
                Me%SubBasin_Annual_Average%SubBasin2%NitrateInRunOff             = 0.0
                Me%SubBasin_Annual_Average%SubBasin2%RunOffOrganicN%Dissolved    = 0.0
                Me%SubBasin_Annual_Average%SubBasin2%RunOffOrganicN%Particulated = 0.0
                Me%SubBasin_Annual_Average%SubBasin2%NitrateInLateralFlow        = 0.0
                Me%SubBasin_Annual_Average%SubBasin2%NitrateInGrowndWaterFlow    = 0.0
                Me%SubBasin_Annual_Average%SubBasin2%RunOffSolubleP              = 0.0
                Me%SubBasin_Annual_Average%SubBasin2%RunOffSedimentActiveP       = 0.0
                Me%SubBasin_Annual_Average%SubBasin2%RunOffSedimentStableP       = 0.0
                Me%SubBasin_Annual_Average%SubBasin2%GrowndWaterSolubleP         = 0.0

                allocate(Me%SubBasin_Annual_Average%SubBasin2%WaterFluxes%GlobalToChannel         (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin_Annual_Average%SubBasin2%WaterFluxes%Precipitation           (1:Me%MaxNumberSubBasin))
                allocate(Me%SubBasin_Annual_Average%SubBasin2%WaterFluxes%ActualEvapotranspiration(1:Me%MaxNumberSubBasin))
                Me%SubBasin_Annual_Average%SubBasin2%WaterFluxes%GlobalToChannel       = 0.0
                Me%SubBasin_Annual_Average%SubBasin2%WaterFluxes%Precipitation                  = 0.0
                Me%SubBasin_Annual_Average%SubBasin2%WaterFluxes%ActualEvapotranspiration         = 0.0

                if (Me%OutPut%Yes) then
                    allocate(Me%SubBasin_Average%HRU%Meteo%SoilWater(1:1,1:Me%MaxNumberSubBasin))
                    allocate(Me%SubBasin_Average%HRU%Plant%LAI      (1:1,1:Me%MaxNumberSubBasin))
                endif
            endif

            if (Me%SwatOptions%ReachesOn) then
                allocate(Me%Reaches%algal_biomass                       (1:Me%MaxNumberSubBasin))
                allocate(Me%Reaches%ammonia                             (1:Me%MaxNumberSubBasin))
                allocate(Me%Reaches%BOD                                 (1:Me%MaxNumberSubBasin))
                allocate(Me%Reaches%chlorophyll_a                       (1:Me%MaxNumberSubBasin))
                allocate(Me%Reaches%dissolved_phosphorus                (1:Me%MaxNumberSubBasin))
                allocate(Me%Reaches%dissolved_oxygen                    (1:Me%MaxNumberSubBasin))
                allocate(Me%Reaches%nitrate                             (1:Me%MaxNumberSubBasin))
                allocate(Me%Reaches%nitrite                             (1:Me%MaxNumberSubBasin))
                allocate(Me%Reaches%organic_nitrogen                    (1:Me%MaxNumberSubBasin))
                allocate(Me%Reaches%organic_phosphorus                  (1:Me%MaxNumberSubBasin))
                allocate(Me%Reaches%saturation_oxygen                   (1:Me%MaxNumberSubBasin))
                allocate(Me%Reaches%temperature                         (1:Me%MaxNumberSubBasin))
                allocate(Me%Reaches%Sediments                           (1:Me%MaxNumberSubBasin))
            endif
        endif

        if (Me%SwatOptions%GlobalWaterOn) then
            allocate(Me%GlobalWaterFluxes(12))
            Me%GlobalWaterFluxes = 0.0
        endif

        if (Me%SwatOptions%GlobalNutrientsOn) then
            allocate(Me%GlobalNutrients(8))
            Me%GlobalNutrients = 0.0
        endif

        if (Me%SwatOptions%BasinAverage.and.Me%SwatOptions%SubBasinOn) then
            allocate(Me%Basin_Average(9))
            Me%Basin_Average = 0.0
        endif

    end subroutine AllocateVariables
    !--------------------------------------------------------------------------
    subroutine ConstructHDF5Output

        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        integer                                             :: HDF5_CREATE, STAT_CALL
        
        Me%OutPut%NextOutPut = 1  
        
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5_Reaches, "Reach.hdf5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleSwat - ERR10'

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5_SubBasin, "SubBasin.hdf5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleSwat - ERR11'
        
    end subroutine ConstructHDF5Output

    !--------------------------------------------------------------------------

    subroutine Construct_Time_Serie

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: Dummy
        character(len=StringLength)                         :: HRUTimeSerieLocationFile
        character(len=StringLength)                         :: SubBasinTimeSerieLocationFile
        character(len=StringLength)                         :: GlobalTimeSerieLocationFile
        character(len=StringLength)                         :: LandUseTimeSerieLocationFile

        !----------------------------------------------------------------------

        call GetData(Dummy,                                             &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModuleSwat',                         &
                     Default      = 'OutPuts.dat',                                 &
                     STAT         = STAT_CALL)
        if (iflag==1)                                                    &
            stop 'Replace TIME_SERIE_LOCATION with HRU_TIME_SERIE_LOCATION or SUBBASIN_TIME_SERIE_LOCATION'

        if (Me%SwatOptions%PlantOutputOn  .or.&
            Me%SwatOptions%NitrateOutputOn.or. &
            Me%SwatOptions%MeteoOutputOn  .or. &
            Me%SwatOptions%UsleOutputOn   .or.  &
            Me%SwatOptions%WaterBalanceOn .or.  &
            Me%SwatOptions%OperationalOutputOn .or.  &
            Me%SwatOptions%PhosphorusOutputOn) then
            call GetData(HRUTimeSerieLocationFile,                                             &
                         Me%ObjEnterData,iflag,                                             &
                         SearchType   = FromFile,                                           &
                         keyword      = 'HRU_TIME_SERIE_LOCATION',                              &
                         ClientModule = 'ModuleSwat',                         &
                         Default      = 'TS_HRU.dat',                                 &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                    &
                stop 'Construct_Plant_Time_Serie - ModuleSwat - ERR01'
        
            Me%HRU%WaterPoints2D      = 1

            if (Me%SwatOptions%MeteoOutputOn) &
                call Construct_Meteo_Time_Serie (HRUTimeSerieLocationFile)

            if (Me%SwatOptions%PlantOutputOn) &
                call Construct_Plant_Time_Serie (HRUTimeSerieLocationFile)
                !To obtain type of crop: ncr(nro(j),icr(j),j)

            if (Me%SwatOptions%UsleOutputOn) &
                call Construct_Usle_Time_Serie (HRUTimeSerieLocationFile)

            if (Me%SwatOptions%NitrateOutputOn) &
                call Construct_Nitrate_Time_Serie (HRUTimeSerieLocationFile)

            if (Me%SwatOptions%PhosphorusOutputOn) &
                call Construct_Phosphorus_Time_Serie (HRUTimeSerieLocationFile)

            if (Me%SwatOptions%WaterBalanceOn) &
                call Construct_WaterBalance_Time_Serie (HRUTimeSerieLocationFile)

            if (Me%SwatOptions%OperationalOutputOn) &
                call Construct_Operational_Time_Serie (HRUTimeSerieLocationFile)
        endif

        if (Me%SwatOptions%SubBasinOn.or.Me%SwatOptions%ReachesOn) then
            call GetData(SubBasinTimeSerieLocationFile,                                             &
                         Me%ObjEnterData,iflag,                                             &
                         SearchType   = FromFile,                                           &
                         keyword      = 'SUBBASIN_TIME_SERIE_LOCATION',                              &
                         ClientModule = 'ModuleSwat',                         &
                         Default      = 'TS_SubBasin.dat',                                 &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                    &
                stop 'Construct_Plant_Time_Serie - ModuleSwat - ERR01'
        
            Me%SubBasin%WaterPoints1D = 1
        
            if (Me%SwatOptions%SubBasinOn) &
                call Construct_SubBasin_Time_Serie(SubBasinTimeSerieLocationFile)

            if (Me%SwatOptions%ReachesOn) &
                call Construct_ReachesOn_Time_Serie(SubBasinTimeSerieLocationFile)

        endif

        if (Me%SwatOptions%GlobalWaterOn.or.Me%SwatOptions%GlobalNutrientsOn.or. Me%SwatOptions%BasinAverage) then

            call GetData(GlobalTimeSerieLocationFile,                                             &
                         Me%ObjEnterData,iflag,                                             &
                         SearchType   = FromFile,                                           &
                         keyword      = 'GLOBAL_TIME_SERIE_LOCATION',                              &
                         ClientModule = 'ModuleSwat',                         &
                         Default      = 'TS_Global.dat',                                 &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                    &
                stop 'Construct_Plant_Time_Serie - ModuleSwat - ERR01'
            
            if (Me%SwatOptions%GlobalWaterOn) &
            call Construct_GlobalWater_Time_Serie (GlobalTimeSerieLocationFile)
    
            if (Me%SwatOptions%GlobalNutrientsOn) &
            call Construct_GlobalNutrients_Time_Serie (GlobalTimeSerieLocationFile)

            if (Me%SwatOptions%BasinAverage.and.Me%SwatOptions%SubBasinOn) &
            call Construct_Average_Basin (GlobalTimeSerieLocationFile)

        endif            

        if (Me%SwatOptions%LandUse) then

            call GetData(LandUseTimeSerieLocationFile,                                             &
                         Me%ObjEnterData,iflag,                                             &
                         SearchType   = FromFile,                                           &
                         keyword      = 'LAND_USE_TIME_SERIE_LOCATION',                     &
                         ClientModule = 'ModuleSwat',                                       &
                         Default      = 'TS_LandUse.dat',                                    &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                    &
                stop 'Construct_Time_Serie - ModuleSwat - ERR01'
            
            call Construct_LandUse_Time_Serie (LandUseTimeSerieLocationFile)
    
        endif            

        
    end subroutine Construct_Time_Serie
    !--------------------------------------------------------------------------
    subroutine Construct_Meteo_Time_Serie (TimeSerieLocationFile)

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !----------------------------------------------------------------------

        !Allocates PropertyList
        allocate(PropertyList(1:17))


        PropertyList(1) = 'AverageTemp_ºC'         
        PropertyList(2) = 'Precipitation_mm'
        PropertyList(3) = 'PotencialET_mm'
        PropertyList(4) = 'CulturalT_mm'
        PropertyList(5) = 'ActualET_mm'
        PropertyList(6) = 'MinimumTemp_ºC'
        PropertyList(7) = 'MaximumTemp_ºC'
        PropertyList(8) = 'WindSpeed_m/s'         
        PropertyList(9) = 'RelativeHumidity'
        PropertyList(10) = 'SolarRadiation_MJ/m^2'   
        PropertyList(11) = 'Irrigation_mm'   
        PropertyList(12) = 'PrecipitationAsSnow_mmH2O'
        PropertyList(13) = 'SnowMelt_mmH2O'           
        PropertyList(14) = 'SnowPack_mmH2O'     
        PropertyList(15) = 'SoilTemperature_ºC'     
        PropertyList(16) = 'SoilWater_mmH2O'     
        PropertyList(17) = 'AquiferWater_mmH2O'
        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjMeteoTimeSerie,                                 &
                            Me%ObjTime,                                      &
                            TimeSerieLocationFile,                           &
                            PropertyList, "mto",                             &
                            WaterPoints2D = Me%HRU%WaterPoints2D,    &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Meteo_Time_Serie - ModuleSwat - ERR01'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Meteo_Time_Serie - ModuleSwat - ERR01'

        
    end subroutine Construct_Meteo_Time_Serie
    !--------------------------------------------------------------------------
    subroutine Construct_Plant_Time_Serie (TimeSerieLocationFile)

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !----------------------------------------------------------------------

        !Allocates PropertyList
        allocate(PropertyList(1:9))


        PropertyList(1) = 'Biomass_kg/ha'
        PropertyList(2) = 'NitrogenFraction'
        PropertyList(3) = 'NitrogenYield_kgN/kgYield'
        PropertyList(4) = 'FractionOfHeatUnits'
        PropertyList(5) = 'LAI'
        PropertyList(6) = 'StressP_FractGrowAchieved'
        PropertyList(7) = 'StressN_FractGrowAchieved'
        PropertyList(8) = 'StressTemp_FractGrowAchieved'
        PropertyList(9) = 'StressWater_FractGrowAchieved'

        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjPlantTimeSerie,                                 &
                            Me%ObjTime,                                      &
                            TimeSerieLocationFile,                           &
                            PropertyList, "plt",                             &
                            WaterPoints2D = Me%HRU%WaterPoints2D,    &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Plant_Time_Serie - ModuleSwat - ERR01'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Plant_Time_Serie - ModuleSwat - ERR01'

        
    end subroutine Construct_Plant_Time_Serie
      !--------------------------------------------------------------------------
    subroutine Construct_Operational_Time_Serie (TimeSerieLocationFile)

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !----------------------------------------------------------------------

        !Allocates PropertyList
        allocate(PropertyList(1:26))


        PropertyList(1) = 'AverageTemp_ºC'         
        PropertyList(2) = 'Precipitation_mm'
        PropertyList(3) = 'PotencialET_mm'
        PropertyList(4) = 'CulturalT_mm'
        PropertyList(5) = 'ActualET_mm'
        PropertyList(6) = 'MinimumTemp_ºC'
        PropertyList(7) = 'MaximumTemp_ºC'
        PropertyList(8) = 'WindSpeed_m/s'         
        PropertyList(9) = 'RelativeHumidity'
        PropertyList(10) = 'SolarRadiation_MJ/m^2'   
        PropertyList(11) = 'Irrigation_mm'   
        PropertyList(12) = 'PrecipitationAsSnow_mmH2O'
        PropertyList(13) = 'SnowMelt_mmH2O'           
        PropertyList(14) = 'SnowPack_mmH2O'     
        PropertyList(15) = 'SoilTemperature_ºC'     
        PropertyList(16) = 'SoilWater_mmH2O'     
        PropertyList(17) = 'AquiferWater_mmH2O'
        PropertyList(18) = 'Biomass_kg/ha'
        PropertyList(19) = 'NitrogenFraction'
        PropertyList(20) = 'NitrogenYield_kgN/kgYield'
        PropertyList(21) = 'FractionOfHeatUnits'
        PropertyList(22) = 'LAI'
        PropertyList(23) = 'StressP_FractGrowAchieved'
        PropertyList(24) = 'StressN_FractGrowAchieved'
        PropertyList(25) = 'StressTemp_FractGrowAchieved'
        PropertyList(26) = 'StressWater_FractGrowAchieved'
        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjOperationalTimeSerie,                                 &
                            Me%ObjTime,                                      &
                            TimeSerieLocationFile,                           &
                            PropertyList, "opt",                             &
                            WaterPoints2D = Me%HRU%WaterPoints2D,    &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Meteo_Time_Serie - ModuleSwat - ERR01'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Meteo_Time_Serie - ModuleSwat - ERR01'

        
    end subroutine Construct_Operational_Time_Serie

    !--------------------------------------------------------------------------
    subroutine Construct_Usle_Time_Serie (TimeSerieLocationFile)

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !----------------------------------------------------------------------

        !Allocates PropertyList
        allocate(PropertyList(1:6))


        PropertyList(1) = 'Musle_Sed_ton/ha'
        PropertyList(2) = 'Usle_Sed_ton/ha'
        PropertyList(3) = 'C_K_P_LS_rock'
        PropertyList(4) = 'MuslePeakRunoffRate_m3/s'
        PropertyList(5) = 'UsleRainfallIndex_100(ft-tn*in)/(acre-hr)'
        PropertyList(6) = 'K_LS_P_rock'


        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjUsleTimeSerie,                                 &
                            Me%ObjTime,                                      &
                            TimeSerieLocationFile,                           &
                            PropertyList, "usl",                             &
                            WaterPoints2D = Me%HRU%WaterPoints2D,    &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Plant_Time_Serie - ModuleSwat - ERR01'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Plant_Time_Serie - ModuleSwat - ERR01'

        
    end subroutine Construct_Usle_Time_Serie
    !--------------------------------------------------------------------------
    subroutine Construct_Nitrate_Time_Serie (TimeSerieLocationFile)

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !----------------------------------------------------------------------

        !Allocates PropertyList
        allocate(PropertyList(1:19))

        PropertyList(1)  = 'ActiveOrgNToNitrate_kgN/ha'     
        PropertyList(2)  = 'FreshOrgNToNitrate_kgN/ha'     
        PropertyList(3)  = 'Denitrification_kgN/ha'     
        PropertyList(4)  = 'Nitrate_LateralFlow_kgN/ha'
        PropertyList(5)  = 'Nitrate_Rain_kgN/ha'
        PropertyList(6)  = 'Nitrate_PlantUptake_kgN/ha'
        PropertyList(7)  = 'Nitrate_Percolation_kgN/ha'
        PropertyList(8)  = 'Nitrate_RunOff_kgN/ha'
        PropertyList(9)  = 'OrganicN_RunOff_kgN/ha'
        PropertyList(10) = 'AutoFertilization_kgN/ha'
        PropertyList(11) = 'Fertilization_kgN/ha'
        PropertyList(12) = 'Pools_Nitrate_kgN/ha' 
        PropertyList(13) = 'Pools_Fresh_OM_kgN/ha' 
        PropertyList(14) = 'Pools_Active_OM_kgN/ha' 
        PropertyList(15) = 'Pools_Stable_OM_kgN/ha' 
        PropertyList(16) = 'Grazing_kgN/ha'
        PropertyList(17) = 'Volatilization_kgN/ha'
        PropertyList(18) = 'GrowndWaterNitrateConcentration_mg/L'
        PropertyList(19) = 'SoilOrganicMatterResidue_kg/ha'
                      
        !Constructs TimeSerie'
        call StartTimeSerie(Me%ObjNitrateTimeSerie,                                 &
                            Me%ObjTime,                                      &
                            TimeSerieLocationFile,                           &
                            PropertyList, "nit",                             &
                            WaterPoints2D = Me%HRU%WaterPoints2D,                &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Nitrate_Time_Serie - ModuleSwat - ERR01'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Nitrate_Time_Serie - ModuleSwat - ERR01'

        
    end subroutine Construct_Nitrate_Time_Serie
    !--------------------------------------------------------------------------
    subroutine Construct_Phosphorus_Time_Serie (TimeSerieLocationFile)

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !----------------------------------------------------------------------

        !Allocates PropertyList
        allocate(PropertyList(1:11))

        PropertyList(1)  = 'ActiveMineral_kgP/ha'     
        PropertyList(2)  = 'FreshOrganicMatter_kgP/ha'     
        PropertyList(3)  = 'OrganicMatter_kgP/ha'     
        PropertyList(4)  = 'StableMineral_kgP/ha'
        PropertyList(5)  = 'AutoFertilization_kgP/ha'
        PropertyList(6)  = 'Fertilization_kgP/ha'
        PropertyList(7)  = 'Grazing_kgP/ha'
        PropertyList(8)  = 'RunOffActiveMineral_kgP/ha'
        PropertyList(9)  = 'RunOffStableMineral_kgP/ha'
        PropertyList(10) = 'RunOffOrganicMatter_kgP/ha'
        PropertyList(11) = 'PlantUptake_kgP/ha' 
                      
        !Constructs TimeSerie'
        call StartTimeSerie(Me%ObjPhosphorusTimeSerie,                                 &
                            Me%ObjTime,                                      &
                            TimeSerieLocationFile,                           &
                            PropertyList, "phs",                             &
                            WaterPoints2D = Me%HRU%WaterPoints2D,                &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Phosphorus_Time_Serie - ModuleSwat - ERR01'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Phosphorus_Time_Serie - ModuleSwat - ERR01'

        
    end subroutine Construct_Phosphorus_Time_Serie
    !--------------------------------------------------------------------------
    subroutine Construct_WaterBalance_Time_Serie (TimeSerieLocationFile)

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !----------------------------------------------------------------------

        !Allocates PropertyList
        allocate(PropertyList(1:26))

        PropertyList(1)  = 'SnowPack_mmH2O'
        PropertyList(2)  = 'SnowPack_mmH2O_Old'  
        PropertyList(3)  = 'SoilWater'
        PropertyList(4)  = 'SoilWaterOld'
        PropertyList(5)  = 'AquiferWater'
        PropertyList(6)  = 'AquiferWaterOld'
        PropertyList(7)  = 'DeepAquifer'
        PropertyList(8)  = 'DeepAquiferOld'
        PropertyList(9)  = 'SurfaceWaterColumn'
        PropertyList(10) = 'SurfaceWaterColumnOld'
        PropertyList(11) = 'LateralFlowStorage'
        PropertyList(12) = 'LateralFlowStorageOld'

        PropertyList(13) = 'Precipitation'
        PropertyList(14) = 'SnowEvaporation'
        PropertyList(15) = 'RunOffFlowToChannel'
        PropertyList(16) = 'LateralFlowToChannel'
        PropertyList(17) = 'ActualEvapoTranspiration'
        PropertyList(18) = 'GrowndWaterFlowToChannel'
        PropertyList(19) = 'GrowndWaterFlowToAtmosphere'
        PropertyList(20) = 'PondWaterLoss'
        PropertyList(21) = 'WetlandWaterLoss'
        PropertyList(22) = 'Irrigation'
        PropertyList(23) = 'RechargeDeepAndShalowAquifer'
        PropertyList(24) = 'DrainageTileFlowToChannel'
        PropertyList(25) = 'SoilFlowToRechargeFlow'
        PropertyList(26) = 'WaterBalanceError'
                      
        !Constructs TimeSerie'
        call StartTimeSerie(Me%ObjWaterBalanceTimeSerie,                     &
                            Me%ObjTime,                                      &
                            TimeSerieLocationFile,                           &
                            PropertyList, "wbl",                             &
                            WaterPoints2D = Me%HRU%WaterPoints2D,            &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Phosphorus_Time_Serie - ModuleSwat - ERR01'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Phosphorus_Time_Serie - ModuleSwat - ERR01'

        
    end subroutine Construct_WaterBalance_Time_Serie
    !--------------------------------------------------------------------------
    subroutine Construct_SubBasin_Time_Serie (TimeSerieLocationFile)

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !----------------------------------------------------------------------

        !Allocates PropertyList
        allocate(PropertyList(12))

        PropertyList(1)  = 'GlobalDischarge'     
        PropertyList(2)  = 'GroundWaterDischarge'     
        PropertyList(3)  = 'RunOffDischarge'     
        PropertyList(4)  = 'LateralFlowDischarge'
        PropertyList(5)  = 'ReachFlowOut'
        PropertyList(6)  = 'NitrateInRunOff'
        PropertyList(7)  = 'DONInRunOff'
        PropertyList(8)  = 'PONInRunOff'
        PropertyList(9)  = 'POPPInRunOff'
        PropertyList(10) = 'MineralPInRunOff'
        PropertyList(11) = 'Sediments_concentration'
        PropertyList(12) = 'LFAndRunOffDischarge'
                      
        !Constructs TimeSerie'
        call StartTimeSerie(Me%ObjSubBasinTimeSerie,                                 &
                            Me%ObjTime,                                      &
                            TimeSerieLocationFile,                           &
                            PropertyList, "dis",                             &
                            WaterPoints1D = Me%SubBasin%WaterPoints1D,    &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_SubBasin_Time_Serie - ModuleSwat - ERR01'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_SubBasin_Time_Serie - ModuleSwat - ERR01'

        
    end subroutine Construct_SubBasin_Time_Serie
    !--------------------------------------------------------------------------
    subroutine Construct_ReachesOn_Time_Serie (TimeSerieLocationFile)

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !----------------------------------------------------------------------

        !Allocates PropertyList
        allocate(PropertyList(14))

        PropertyList(1)  = 'FlowOut_m3/s'             
        PropertyList(2)  = 'algal_biomass_mg/l'           
        PropertyList(3)  = 'ammonia_mg/l'             
        PropertyList(4)  = 'BOD_mg/l'                 
        PropertyList(5)  = 'chlorophyll_a_mg/l'       
        PropertyList(6)  = 'dissolved_phosphorus_mg/l'
        PropertyList(7)  = 'dissolved_oxygen_mg/l'    
        PropertyList(8)  = 'nitrate_mg/l'             
        PropertyList(9)  = 'nitrite_mg/l'             
        PropertyList(10) = 'organic_nitrogen_mg/l'    
        PropertyList(11) = 'organic_phosphorus_mg/l'  
        PropertyList(12) = 'saturation_oxygen_mg/l'   
        PropertyList(13) = 'temperature_ºC'   
        PropertyList(14) = 'sediments_tons'   

        !Constructs TimeSerie'
        call StartTimeSerie(Me%ObjReachesTimeSerie,                                 &
                            Me%ObjTime,                                      &
                            TimeSerieLocationFile,                           &
                            PropertyList, "rch",                             &
                            WaterPoints1D = Me%SubBasin%WaterPoints1D,    &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_ReachesOn_Time_Serie - ModuleSwat - ERR01'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_ReachesOn_Time_Serie - ModuleSwat - ERR01'

        
    end subroutine Construct_ReachesOn_Time_Serie
    !--------------------------------------------------------------------------  
    subroutine Construct_GlobalWater_Time_Serie (TimeSerieLocationFile)

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList
        character(len=StringLength)                         :: GlobalWater

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !----------------------------------------------------------------------

        GlobalWater = 'GlobalWater'

        !Allocates PropertyList
        allocate(PropertyList(12))
        
        if(Me%SwatOptions%GlobalWaterM3) then
            PropertyList(1)  = 'precipitation_m3'             
            PropertyList(2)  = 'runoff_m3'           
            PropertyList(3)  = 'lateral_flow_m3'             
            PropertyList(4)  = 'percolation_m3'                 
            PropertyList(5)  = 'water_yield_m3'       
            PropertyList(6)  = 'act_evapotranspiration_m3'
            PropertyList(7)  = 'soil_stored_water_m3'
            PropertyList(8)  = 'snow_fall_freezing_rain_m3'    
            PropertyList(9)  = 'groundwater_m3'             
            PropertyList(10)  = 'deep_aquifer_recharge_m3'
            PropertyList(11) = 'pot_evapotranspiration_m3'             
            PropertyList(12) = 'drainage_tile_m3'  
        else
            PropertyList(1)  = 'precipitation_mm'             
            PropertyList(2)  = 'runoff_mm'           
            PropertyList(3)  = 'lateral_flow_mm'             
            PropertyList(4)  = 'percolation_mm'                 
            PropertyList(5)  = 'water_yield_mm'       
            PropertyList(6)  = 'act_evapotranspiration_mm'
            PropertyList(7)  = 'soil_stored_water_mm'
            PropertyList(8)  = 'snow_fall_freezing_rain_mm'    
            PropertyList(9)  = 'groundwater_mm'             
            PropertyList(10)  = 'deep_aquifer_recharge_mm'
            PropertyList(11) = 'pot_evapotranspiration_mm'             
            PropertyList(12) = 'drainage_tile_mm'    
        endif  

        !Constructs TimeSerie'
        call StartTimeSerie(Me%ObjGlobalWaterTimeSerie,                      &
                            Me%ObjTime,                                      &
                            TimeSerieLocationFile,                           &
                            PropertyList, "glb",                             &
!                            WaterPoints1D = Me%SubBasin%WaterPoints1D,    &
                            ResultFileName = GlobalWater,&
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_GlobalWater_Time_Serie - ModuleSwat - ERR01'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_GlobalWater_Time_Serie - ModuleSwat - ERR01'

        
    end subroutine Construct_GlobalWater_Time_Serie
    !--------------------------------------------------------------------------  

    subroutine Construct_Average_Basin (TimeSerieLocationFile)

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList
        character(len=StringLength)                         :: Average_Basin

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !----------------------------------------------------------------------

        Average_Basin = 'Average_Basin_kg_ha'

        !Allocates PropertyList
        allocate(PropertyList(9))
        
        PropertyList(1)  = 'RunOffSolubleP'
        PropertyList(2)  = 'RunOffSedimentStableP'
        PropertyList(3)  = 'RunOffSedimentActiveP'       
        PropertyList(4)  = 'GrowndWaterSolubleP'
        PropertyList(5)  = 'RunOffNitrate'
        PropertyList(6)  = 'RunOffOrganicN_Dissolved'
        PropertyList(7)  = 'RunOffOrganicN_Particulated'
        PropertyList(8)  = 'LateralFlowNitrate'
        PropertyList(9)  = 'GrowndWaterFlowNitrate'

        !Constructs TimeSerie'
        call StartTimeSerie(Me%ObjBasinAverage,                              &
                            Me%ObjTime,                                      &
                            TimeSerieLocationFile,                           &
                            PropertyList, "glb",                             &
!                            WaterPoints1D = Me%SubBasin%WaterPoints1D,    &
                            ResultFileName = Average_Basin,&
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Average_Basin - ModuleSwat - ERR01'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Average_Basin - ModuleSwat - ERR01'

        
    end subroutine Construct_Average_Basin
    !--------------------------------------------------------------------------  
    subroutine Construct_GlobalNutrients_Time_Serie (TimeSerieLocationFile)

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList
        character(len=StringLength)                         :: GlobalNutrient

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !----------------------------------------------------------------------

        GlobalNutrient = 'GlobalNutrient'

        !Allocates PropertyList
        allocate(PropertyList(8))

        PropertyList(1)  = 'sed_from_HRUs_ton_ha'
        PropertyList(2)  = 'sed_to_ponds_ton'
        PropertyList(3)  = 'sed_from_ponds_ton'  
        PropertyList(4)  = 'netChangeSedPonds_ton'
        PropertyList(5)  = 'sed_to_wetlands_ton'
        PropertyList(6)  = 'sed_from_wetlands_ton'
        PropertyList(7)  = 'netChangeSedWetlands_ton'
        PropertyList(8)  = 'N_added_agriculture_ton'

        !Constructs TimeSerie'
        call StartTimeSerie(Me%ObjGlobalNutrientTimeSerie,                      &
                            Me%ObjTime,                                      &
                            TimeSerieLocationFile,                           &
                            PropertyList, "glb",                             &
!                            WaterPoints1D = Me%SubBasin%WaterPoints1D,    &
                            ResultFileName = GlobalNutrient,&
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_GlobalNutrients_Time_Serie - ModuleSwat - ERR01'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_GlobalNutrients_Time_Serie - ModuleSwat - ERR01'

        
    end subroutine Construct_GlobalNutrients_Time_Serie

    !--------------------------------------------------------------------------  

    subroutine Construct_LandUse_Time_Serie (TimeSerieLocationFile)

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !----------------------------------------------------------------------

        !Allocates PropertyList
        allocate(PropertyList(41))

        PropertyList(1)  = 'Fert_N'		
        PropertyList(2)  = 'Fert_P'		
        PropertyList(3)  = 'AutoFert_N'		
        PropertyList(4)  = 'AutoFert_P'
        PropertyList(5)  = 'Grazing_N'
        PropertyList(6)  = 'Grazing_P'
        PropertyList(7)  = 'Contin_Fert_N'
        PropertyList(8)  = 'Contin_Fert_P'
        PropertyList(9)  = 'pcp_NO3_N'
        PropertyList(10) = 'Fixation_N'
        PropertyList(11) = 'FreshOM_NToNitrate'
        PropertyList(12) = 'ActiveOM_NToNitrate'
        PropertyList(13) = 'ActiveOM_NToStable'
        PropertyList(14) = 'FreshOM_PTo_'!'rmptl'
        PropertyList(15)  = "OrgainicP_To_LabilP" !'hmptl'
        PropertyList(16)  = 'LabilMineral_P_ActiveMiner'
        PropertyList(17)  = 'ActiveMineral_P_Stable_Mine'
        PropertyList(18)  = 'Denitrification'
        PropertyList(19)  = 'PlantUptake_N'
        PropertyList(20)  = 'PlantUptake_P'
        PropertyList(21)  = 'Organic_N_RunOff'
        PropertyList(22)  = 'Organic_P_RunOff'
        PropertyList(23)  = 'Mineral_P_RunOff' !sedminpa+sedminps
        PropertyList(24) =  'NO3_RunOff'		
        PropertyList(25) = 'NO3_LateralFlow'		
        PropertyList(26) = 'NO3_Percolation'		
        PropertyList(27) = 'NO3_GrowndWater'		
        PropertyList(28) = 'Soluble_P_RunOff'!'surqsolp(j)'		
        PropertyList(29) = 'Soluble_P_ToReach' !'minpgw(j)'		

        PropertyList(30) = 'Volatilization'	
        PropertyList(31) = 'FreshOM_NToNitrate'
        PropertyList(32) = 'Nitrate'    	
        PropertyList(33) = 'ActiveOM_N'   	
        PropertyList(34) = 'StableOM_N'  		
        PropertyList(35) = 'FreshOM_N'  	
        PropertyList(36) = 'Ammonium'  	
        PropertyList(37) = 'Nitrification'  	

        PropertyList(38) = 'ActiveMineral_P'
        PropertyList(39) = 'FreshOrganicMatter_P'
        PropertyList(40) = 'OrganicMatter_P' 
        PropertyList(41) = 'StableMineral_P'

        !Constructs TimeSerie'
        call StartTimeSerie(Me%ObjLandUseTimeSerie,                          &
                            Me%ObjTime,                                      &
                            TimeSerieLocationFile,                           &
                            PropertyList, "ldu",                             &
                            WaterPoints1D = Me%SubBasin%WaterPoints1D,       &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_LandUse_Time_Serie - ModuleSwat - ERR01'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_LandUse_Time_Serie - ModuleSwat - ERR01'

        
    end subroutine Construct_LandUse_Time_Serie
    !--------------------------------------------------------------------------  

	subroutine LandUseList

        !Local---------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        integer                                     :: ClientNumber
        logical                                     :: BlockFound
        character(len=StringLength)                 :: LandUseTimeSerieLocationFile
		integer										:: ObjEnterData
		integer										:: BlockNumber

        !Begin---------------------------------------------------------

        call GetData(LandUseTimeSerieLocationFile,                                             &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'LAND_USE_TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModuleSwat',                         &
                     Default      = 'TS_LandUse.dat',                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'LandUseList - ModuleSwat - ERR00'

        call ConstructEnterData(ObjEnterData, LandUseTimeSerieLocationFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'LandUseList - ModuleSwat - ERR00'

		Me%HRU%NumberOfLandUseOutputs = 0

do1 :   do
            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,              &
                                        '<BeginTimeSerie>', '<EndTimeSerie>',           &
                                        BlockFound, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'LandUseList - ModuleSwat - ERR00'

if1 :       if(STAT_CALL .EQ. SUCCESS_) then 
   
if2 :           if (BlockFound) then

					Me%HRU%NumberOfLandUseOutputs = Me%HRU%NumberOfLandUseOutputs + 1

                else

                    call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop 'LandUseList - ModuleSwat - ERR00'
                        
                    exit do1    !No more blocks

                end if if2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                if(STAT_CALL .ne. SUCCESS_)stop 'LandUseList - ModuleSwat - ERR00'
                    
            end if if1
        end do do1

		call KillEnterData(ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'LandUseList - ModuleSwat - ERR00'


        call ConstructEnterData(ObjEnterData, LandUseTimeSerieLocationFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'LandUseList - ModuleSwat - ERR00'



		allocate(Me%HRU%LandUse(Me%HRU%NumberOfLandUseOutputs))
		allocate(Me%HRU%LandUseAreaHectare(Me%HRU%NumberOfLandUseOutputs))
        Me%HRU%LandUseAreaHectare = 0
		allocate(Me%HRU%LandUseNutrient (80,Me%HRU%NumberOfLandUseOutputs))

		BlockNumber = 0

do2 :   do
            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,              &
                                        '<BeginTimeSerie>', '<EndTimeSerie>',           &
                                        BlockFound, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'LandUseList - ModuleSwat - ERR00'

if3 :       if(STAT_CALL .EQ. SUCCESS_) then 
   
if4 :           if (BlockFound) then

					BlockNumber = BlockNumber + 1

                    call GetData(Me%HRU%LandUse(BlockNumber),                              &
                                 ObjEnterData, iflag,                            &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'NAME',                             &
                                 default      = 'NOCR',                             &
                                 ClientModule = 'ModuleSWAT',                     &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'LandUseList - ModuleSwat - ERR00'


                else

                    call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop 'LandUseList - ModuleSwat - ERR00'
                        
                    exit do2    !No more blocks

                end if if4

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if3
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                if(STAT_CALL .ne. SUCCESS_)stop 'LandUseList - ModuleSwat - ERR00'
                    
            end if if3
        end do do2

		call KillEnterData(ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'LandUseList - ModuleSwat - ERR00'


	end subroutine LandUseList

    !--------------------------------------------------------------------------  

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Swat), pointer                         :: NewObjSwat
        type (T_Swat), pointer                         :: PreviousObjSwat


        !Allocates new instance
        allocate (NewObjSwat)
        nullify  (NewObjSwat%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjSwat)) then
            FirstObjSwat         => NewObjSwat
            Me                    => NewObjSwat
        else
            PreviousObjSwat      => FirstObjSwat
            Me                    => FirstObjSwat%Next
            do while (associated(Me))
                PreviousObjSwat  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjSwat
            PreviousObjSwat%Next => NewObjSwat
        endif

        Me%InstanceID = RegisterNewInstance (mSWAT_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine GetSwatWithChanges (ObjSwatID, SwatWithChanges, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSwatID
        logical                                         :: SwatWithChanges
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSwatID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            SwatWithChanges = Me%SwatOptions%SwatWithChanges

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetSwatWithChanges

    !--------------------------------------------------------------------------    
    
    !--------------------------------------------------------------------------
    subroutine GetSwatPointer (ObjSwatID, Matrix, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSwatID
        real(8), dimension(:, :, :),  pointer           :: Matrix
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSwatID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mSWAT_, Me%InstanceID)

            Matrix => Me%Matrix

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetSwatPointer
    
    !--------------------------------------------------------------------------
    
    subroutine GetSwatInteger (ObjSwatID, Int, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSwatID
        real                                            :: Int
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSwatID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Int = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetSwatInteger

    !--------------------------------------------------------------------------

    subroutine UnGetSwat3D_I(ObjSwatID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSwatID
        integer, dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSwatID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mSWAT_, Me%InstanceID, "UnGetSwat3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetSwat3D_I

    !--------------------------------------------------------------------------

    subroutine UnGetSwat3D_R8(ObjSwatID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSwatID
        real(8), dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSwatID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mSWAT_, Me%InstanceID,  "UnGetSwat3D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetSwat3D_R8

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifySwat(ObjSwatID,OutputHruMeteo,OutputHruQuality,OutputHruWatBal,OutputHruParameters, &
                          OutputSubBasin,OutputSubBasin2,OutputReaches,OutputGlobalBasin,OutputHRU_Cropname,&
						  HRU_output_values,iyr,ida,STAT)
    

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSwatID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL
        real, dimension(:,:,:),pointer              :: OutputHruMeteo
        real, dimension(:,:,:),pointer              :: OutputHruQuality
        real, dimension(:,:,:),pointer              :: OutputHruWatBal
        real, dimension(:,:,:),pointer              :: OutputHruParameters
        real, dimension(:,:),pointer                :: OutputSubBasin
        real, dimension(:,:),pointer                :: OutputSubBasin2
        real, dimension(:,:),pointer                :: OutputReaches
        real, dimension(:),pointer                  :: OutputGlobalBasin
        real, dimension(:,:,:),pointer              :: HRU_output_values
		character(len=4),  dimension(:,:),  pointer :: OutputHRU_Cropname 
        integer                                     :: j,i,l,k,n
        integer                                     :: iyr,ida
        real                                        :: Year, Month, Day, Hour, Minute, Second
        real                                        :: ConversionFactor, Aux, HRUFraction
		real                                        :: dstor, h2oloss,ToChannel
		real                                        :: DeepAqui,SW
		real                                        :: Month_old, Hid_year, year_old
		logical										:: Monthly_out

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSwatID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call GetComputeCurrentTime(Me%ObjTime, Me%CurrentTime)

			Me%Month_old = Me%Month
			year_old  = year

            call ExtractDate(Me%CurrentTime, Year, Me%Month, Day, Hour, Minute, Second)

			if(Me%Month > 10) then
				Hid_year = Year
			else
				Hid_year = Year - 1
			end if

			if ((Me%Month_old - Me%month<0) .or.(Me%Month_old - Me%month>10)) then
				Monthly_out = .true.
			else
				Monthly_out = .false.
			end if

		   Me%HRU%HRU_output_values => HRU_output_values

           do j=1,Me%MaxNumberHRU
                do i=1,Me%MaxNumberSubBasin
                    Me%HRUBasinFraction (j,i) = OutputHruParameters (1,j,i)
                    if(Me%HRUBasinFraction (j,i) > 1.0 .or. Me%HRUBasinFraction (j,i) <= 0) then
                        Me%HRU%WaterPoints2D(j,i) = 0
                    else
                        Me%HRU%WaterPoints2D(j,i) = 1
                    endif
                enddo
            enddo

            do i=1,Me%MaxNumberSubBasin
                Me%SubBasin%SubBasinFraction (i) = OutputSubBasin (13,i)
                if(Me%SubBasin%SubBasinFraction (i) > 1.0 .or. Me%SubBasin%SubBasinFraction (i) <= 0) then
                    Me%SubBasin%WaterPoints1D(i) = 0
                else
                    Me%SubBasin%WaterPoints1D(i) = 1
                endif
            end do

            do i=1,Me%MaxNumberSubBasin
                Me%SubBasin2%SubBasinFraction (i) = OutputSubBasin2 (13,i)
                if(Me%SubBasin2%SubBasinFraction (i) > 1.0 .or. Me%SubBasin2%SubBasinFraction (i) <= 0) then
                    Me%SubBasin2%WaterPoints1D(i) = 0
                else
                    Me%SubBasin2%WaterPoints1D(i) = 1
                endif
            end do

            if (Me%SwatOptions%ReachesOn) then
                do j=1,Me%MaxNumberSubBasin
                    if (Me%SubBasin%WaterPoints1D(j) == 1) then
                        Me%Reaches%FlowOut             (j) = OutputReaches (1,j)
                        Me%Reaches%algal_biomass       (j) = OutputReaches (2,j)
                        Me%Reaches%ammonia             (j) = OutputReaches (3,j)
                        Me%Reaches%BOD                 (j) = OutputReaches (4,j)
                        Me%Reaches%chlorophyll_a       (j) = OutputReaches (5,j)
                        Me%Reaches%dissolved_phosphorus(j) = OutputReaches (6,j)
                        Me%Reaches%dissolved_oxygen    (j) = OutputReaches (7,j)
                        Me%Reaches%nitrate             (j) = OutputReaches (8,j)
                        Me%Reaches%nitrite             (j) = OutputReaches (9,j)
                        Me%Reaches%organic_nitrogen    (j) = OutputReaches (10,j)
                        Me%Reaches%organic_phosphorus  (j) = OutputReaches (11,j)
                        Me%Reaches%saturation_oxygen   (j) = OutputReaches (12,j)
                        Me%Reaches%Temperature         (j) = OutputReaches (13,j)
                        Me%Reaches%Sediments           (j) = OutputReaches (14,j)
                    endif
                enddo
            endif

            do j=1,Me%MaxNumberHRU
                do i=1,Me%MaxNumberSubBasin
                    if(Me%HRU%WaterPoints2D(j,i) == 1) then
                        ![-] = [-] / [-] 
                        HRUFraction = Me%HRUBasinFraction (j,i)  / Me%SubBasin%SubBasinFraction (i)

                        if (Me%SwatOptions%MeteoOutputOn.or.Me%SwatOptions%OperationalOutputOn) then
                            Me%HRU%Meteo%AverageTemperature         (j,i) = OutputHruMeteo (1 ,j,i)    
                            Me%HRU%Meteo%Precipitation              (j,i) = OutputHruMeteo (2 ,j,i)    
                            Me%HRU%Meteo%PotencialEvapoTranspiration(j,i) = OutputHruMeteo (3 ,j,i)    
                            Me%HRU%Meteo%CulturalTranspiration      (j,i) = OutputHruMeteo (4 ,j,i)   
                            Me%HRU%Meteo%ActualEvapoTranspiration   (j,i) = OutputHruMeteo (5 ,j,i)   
                            Me%HRU%Meteo%MinimumTemperature         (j,i) = OutputHruMeteo (6 ,j,i)    
                            Me%HRU%Meteo%MaximumTemperature         (j,i) = OutputHruMeteo (7 ,j,i)    
                            Me%HRU%Meteo%WindSpeed                  (j,i) = OutputHruMeteo (8 ,j,i)    
                            Me%HRU%Meteo%RelativeHumidity           (j,i) = OutputHruMeteo (9 ,j,i)   
                            Me%HRU%Meteo%SolarRadiation             (j,i) = OutputHruMeteo (10,j,i)   
                            Me%HRU%Meteo%Irrigation                 (j,i) = OutputHruMeteo (11,j,i)   
                            Me%HRU%Meteo%PrecipitationAsSnow        (j,i) = OutputHruMeteo (12,j,i)
                            Me%HRU%Meteo%SnowMelt                   (j,i) = OutputHruMeteo (13,j,i)
                            Me%HRU%Meteo%SnowPack_mmH2O             (j,i) = OutputHruMeteo (14,j,i)
                            Me%HRU%Meteo%SoilTemperature            (j,i) = OutputHruMeteo (15,j,i)
                            Me%HRU%Meteo%SoilWater                  (j,i) = OutputHruMeteo (16,j,i)
                            Me%HRU%Meteo%AquiferWater               (j,i) = OutputHruMeteo (17,j,i)
                        endif 

                        if (Me%SwatOptions%PlantOutputOn.or.Me%SwatOptions%OperationalOutputOn) then
                            Me%HRU%Plant%Biomass                      (j,i) = OutputHruQuality (1 ,j,i)    
                            Me%HRU%Plant%NitrogenFraction             (j,i) = OutputHruQuality (2 ,j,i)    
                            Me%HRU%Plant%NitrogenYield                (j,i) = OutputHruQuality (3 ,j,i)    
                            Me%HRU%Plant%FractionOfHeatUnits          (j,i) = OutputHruQuality (4 ,j,i)   
                            Me%HRU%Plant%LAI                          (j,i) = OutputHruQuality (21,j,i)   
                            Me%HRU%Plant%FractionGrowth_StressP       (j,i) = OutputHruQuality (24,j,i)    
                            Me%HRU%Plant%FractionGrowth_StressN       (j,i) = OutputHruQuality (25,j,i)    
                            Me%HRU%Plant%FractionGrowth_StressTemp    (j,i) = OutputHruQuality (26,j,i)   
                            Me%HRU%Plant%FractionGrowth_StressWater   (j,i) = OutputHruQuality (27,j,i)   
                        endif 

                        if (Me%SwatOptions%UsleOutputOn) then
                            Me%HRU%Usle%Musle_Sed                      (j,i) = OutputHruQuality (22 ,j,i)    
                            Me%HRU%Usle%Usle_Sed                       (j,i) = OutputHruQuality (23 ,j,i)
                            Me%HRU%Usle%C_K_P_LS_rock                  (j,i) = OutputHruQuality (28 ,j,i)                      
                            Me%HRU%Usle%MuslePeakRunoffRate            (j,i) = OutputHruQuality (29 ,j,i)                          
                            Me%HRU%Usle%UsleRainfallIndex              (j,i) = OutputHruQuality (30 ,j,i)                          
                            Me%HRU%Usle%K_LS_P_rock                    (j,i) = OutputHruQuality (31 ,j,i)      
                            aux = Me%EndTime - Me%BeginTime
                            Me%SubBasin_Annual_Average%HRU%Usle%Musle_Sed    (1,i)   = Me%SubBasin_Annual_Average%HRU%Usle%Musle_Sed(1,i) + &
                                                                          Me%HRU%Usle%Musle_Sed(j,i) * HRUFraction / &
                                                                          Me%NumberYearsSimulation
                                            
                        endif 

                        if (Me%SwatOptions%NitrateOutputOn) then
                            Me%HRU%RateFlux%ActiveOrgNToNitrate       (j,i) = OutputHruQuality (5 ,j,i)    
                            Me%HRU%RateFlux%FreshOrgNToNitrate        (j,i) = OutputHruQuality (6 ,j,i)    
                            Me%HRU%RateFlux%Denitrification           (j,i) = OutputHruQuality (7 ,j,i)    

                            Me%HRU%NTransportFlux%Nitrate%LateralFlow (j,i) = OutputHruQuality (8 ,j,i)    
                            Me%HRU%NTransportFlux%Nitrate%Rain        (j,i) = OutputHruQuality (9 ,j,i)    
                            Me%HRU%NTransportFlux%Nitrate%PlantUptake (j,i) = OutputHruQuality (10,j,i)    
                            Me%HRU%NTransportFlux%Nitrate%Percolation (j,i) = OutputHruQuality (11,j,i)    
                            Me%HRU%NTransportFlux%Nitrate%RunOff      (j,i) = OutputHruQuality (12,j,i)    
                            Me%HRU%NTransportFlux%RunOffOrganicN      (j,i) = OutputHruQuality (47,j,i)    
 
                            Me%HRU%NTransportFlux%AutoFertilization   (j,i) = OutputHruQuality (13,j,i)    
                            Me%HRU%NTransportFlux%Fertilization       (j,i) = OutputHruQuality (14,j,i)    

                            Me%HRU%PoolsN%Nitrate                      (j,i) = OutputHruQuality (15,j,i)    
                            Me%HRU%PoolsN%FreshOrganicMatter           (j,i) = OutputHruQuality (16,j,i)    
                            Me%HRU%NTransportFlux%Grazing             (j,i) = OutputHruQuality (17,j,i)
                            Me%HRU%NTransportFlux%Volatilization      (j,i) = OutputHruQuality (18,j,i)
                            Me%HRU%PoolsN%ActiveOrganicMatter          (j,i) = OutputHruQuality (19,j,i)    
                            Me%HRU%PoolsN%StableOrganicMatter          (j,i) = OutputHruQuality (20,j,i)    
                            Me%HRU%GW%NitrateConcentration            (j,i) = OutputHruQuality (32,j,i)    
						    Me%HRU%PoolsN%SoilResidue                  (j,i) = OutputHruQuality (33,j,i)   
                            Me%HRU%PoolsN%Ammonium                     (j,i) = OutputHruQuality (34,j,i)    

                            Me%HRU%NTransportFlux%Nitrification       (j,i) = OutputHruQuality (35,j,i)    
                        
                            aux = Me%EndTime - Me%BeginTime
                            !kgN/ha/year = kgN/ha/year + kgN/ha/year
!                            Me%SubBasin_Annual_Average%HRU%NTransportFlux%Nitrate%Percolation(1,i)   = Me%SubBasin_Annual_Average%HRU%NTransportFlux%Nitrate%Percolation(1,i) + &
!                                                                          Me%HRU%NTransportFlux%Nitrate%Percolation(j,i) / &
!                                                                          Me%NumberYearsSimulation
                        endif


                        if (Me%SwatOptions%PhosphorusOutputOn) then
                            Me%HRU%PoolsP%ActiveMineral               (j,i) =  OutputHruQuality(36,j,i)  
                            Me%HRU%PoolsP%FreshOrganicMatter          (j,i) =  OutputHruQuality(37,j,i)  
                            Me%HRU%PoolsP%OrganicMatter               (j,i) =  OutputHruQuality(38,j,i)  
                            Me%HRU%PoolsP%StableMineral               (j,i) =  OutputHruQuality(39,j,i)  
                            Me%HRU%PTransportFlux%AutoFertilization   (j,i) =  OutputHruQuality(40,j,i)  !autop 
                            Me%HRU%PTransportFlux%Fertilization       (j,i) =  OutputHruQuality(41,j,i)  !fertp 
                            Me%HRU%PTransportFlux%Grazing             (j,i) =  OutputHruQuality(42,j,i)  !grazp 
                            Me%HRU%PTransportFlux%RunOffActiveMineral (j,i) =  OutputHruQuality(43,j,i)  !sedminpa(j)  
                            Me%HRU%PTransportFlux%RunOffStableMineral (j,i) =  OutputHruQuality(44,j,i)  !sedminps(j)  
                            Me%HRU%PTransportFlux%RunOffOrganicMatter (j,i) =  OutputHruQuality(45,j,i)  !sedorgp(j)   
                            Me%HRU%PTransportFlux%PlantUptake         (j,i) =  OutputHruQuality(46,j,i)  
                        endif
						
						if (Me%SwatOptions%WaterBalanceOn) then
		! pools of water (dstor)
                            Me%HRU%WatBal%SnowPack_mmH2O						(j,i) = OutputHruWatBal (1 ,j,i)
                            Me%HRU%WatBal%SnowPack_mmH2O_Old					(j,i) = OutputHruWatBal (2 ,j,i)
                            Me%HRU%WatBal%SoilWater								(j,i) = OutputHruWatBal (3 ,j,i)
                            Me%HRU%WatBal%SoilWaterOld							(j,i) = OutputHruWatBal (4 ,j,i)
                            Me%HRU%WatBal%AquiferWater							(j,i) = OutputHruWatBal (5 ,j,i)
                            Me%HRU%WatBal%AquiferWaterOld						(j,i) = OutputHruWatBal (6 ,j,i)
                            Me%HRU%WatBal%DeepAquifer							(j,i) = OutputHruWatBal (7 ,j,i)
                            Me%HRU%WatBal%DeepAquiferOld						(j,i) = OutputHruWatBal (8 ,j,i)
                            Me%HRU%WatBal%SurfaceWaterColumn					(j,i) = OutputHruWatBal (9 ,j,i)
                            Me%HRU%WatBal%SurfaceWaterColumnOld					(j,i) = OutputHruWatBal (10,j,i)
                            Me%HRU%WatBal%LateralFlowStorage					(j,i) = OutputHruWatBal (11,j,i)
                            Me%HRU%WatBal%LateralFlowStorageOld					(j,i) = OutputHruWatBal (12,j,i)
		! water fluxes (h2oloss)
                            Me%HRU%WatBal%Precipitation							(j,i) = OutputHruWatBal (13,j,i)
                            Me%HRU%WatBal%SnowEvaporation						(j,i) = OutputHruWatBal (14,j,i)
                            Me%HRU%WatBal%RunOffFlowToChannel					(j,i) = OutputHruWatBal (15,j,i)
                            Me%HRU%WatBal%LateralFlowToChannel					(j,i) = OutputHruWatBal (16,j,i)
                            Me%HRU%WatBal%ActualEvapoTranspiration				(j,i) = OutputHruWatBal (17,j,i)
                            Me%HRU%WatBal%GrowndWaterFlowToChannel				(j,i) = OutputHruWatBal (18,j,i)
                            Me%HRU%WatBal%GrowndWaterFlowToAtmosphere			(j,i) = OutputHruWatBal (19,j,i)
                            Me%HRU%WatBal%PondWaterLoss							(j,i) = OutputHruWatBal (20,j,i)
                            Me%HRU%WatBal%WetlandWaterLoss						(j,i) = OutputHruWatBal (21,j,i)
                            Me%HRU%WatBal%Irrigation							(j,i) = OutputHruWatBal (22,j,i)
                            Me%HRU%WatBal%RechargeFlowToDeepAndShalowAquifer	(j,i) = OutputHruWatBal (23,j,i)
                            Me%HRU%WatBal%DrainageTileFlowToChannel				(j,i) = OutputHruWatBal (24,j,i)
                            Me%HRU%WatBal%SoilFlowToRechargeFlow             	(j,i) = OutputHruWatBal (25,j,i)

		! pools of water (dstor)
							
							dstor = & 						
                            Me%HRU%WatBal%SnowPack_mmH2O						(j,i) - &
                            Me%HRU%WatBal%SnowPack_mmH2O_Old					(j,i) + &
                            Me%HRU%WatBal%SoilWater								(j,i) - &
                            Me%HRU%WatBal%SoilWaterOld							(j,i) + &
                            Me%HRU%WatBal%AquiferWater							(j,i) - &
                            Me%HRU%WatBal%AquiferWaterOld						(j,i) + &
                            Me%HRU%WatBal%DeepAquifer							(j,i) - &
                            Me%HRU%WatBal%DeepAquiferOld						(j,i) + &
                            Me%HRU%WatBal%SurfaceWaterColumn					(j,i) - &
                            Me%HRU%WatBal%SurfaceWaterColumnOld					(j,i) + &
                            Me%HRU%WatBal%LateralFlowStorage					(j,i) - &
                            Me%HRU%WatBal%LateralFlowStorageOld					(j,i) 
		! water fluxes (h2oloss)
                            h2oloss = & 
							Me%HRU%WatBal%Precipitation							(j,i) - &
                            Me%HRU%WatBal%SnowEvaporation						(j,i) - &
                            Me%HRU%WatBal%RunOffFlowToChannel					(j,i) - &
                            Me%HRU%WatBal%LateralFlowToChannel					(j,i) - &
                            Me%HRU%WatBal%ActualEvapoTranspiration				(j,i) - &
                            Me%HRU%WatBal%GrowndWaterFlowToChannel				(j,i) - &
                            Me%HRU%WatBal%GrowndWaterFlowToAtmosphere			(j,i) + &
                            Me%HRU%WatBal%PondWaterLoss							(j,i) + &
                            Me%HRU%WatBal%WetlandWaterLoss						(j,i) + &
                            Me%HRU%WatBal%Irrigation							(j,i) + &
                            Me%HRU%WatBal%RechargeFlowToDeepAndShalowAquifer	(j,i) - &
                            Me%HRU%WatBal%DrainageTileFlowToChannel				(j,i) - &
                            Me%HRU%WatBal%SoilFlowToRechargeFlow		        (j,i) 

							Me%HRU%WatBal%WaterBalanceError (j,i) =	dstor - h2oloss


!							write(111,777)	&
!							j														  ,",",&
!							i														  ,",",&
!							Year													  ,",",&
!							Month													  ,",",&
!							Day													      ,",",&
!							Me%HRU%WatBal%Precipitation							(j,i) ,",",&
!                            Me%HRU%WatBal%SnowEvaporation						(j,i) ,",",&
!                            Me%HRU%WatBal%RunOffFlowToChannel					(j,i) ,",",&
!                            Me%HRU%WatBal%LateralFlowToChannel					(j,i) ,",",&
!                            Me%HRU%WatBal%ActualEvapoTranspiration				(j,i) ,",",&
!                            Me%HRU%WatBal%GrowndWaterFlowToChannel				(j,i) ,",",&
!                            Me%HRU%WatBal%GrowndWaterFlowToAtmosphere			(j,i) ,",",&
!                            Me%HRU%WatBal%PondWaterLoss							(j,i) ,",",&
!                            Me%HRU%WatBal%WetlandWaterLoss						(j,i) ,",",&
!                            Me%HRU%WatBal%Irrigation							(j,i) ,",",&
!                            Me%HRU%WatBal%RechargeFlowToDeepAndShalowAquifer	(j,i) ,",",&
!                            Me%HRU%WatBal%DrainageTileFlowToChannel				(j,i) ,",",&
!                            Me%HRU%WatBal%SoilFlowToRechargeFlow		        (j,i) ,",",&
!                            Me%HRU%WatBal%SoilWater								(j,i) ,",",&				
!                            Me%HRU%WatBal%SoilWaterOld							(j,i) ,",",&			
!                            Me%HRU%WatBal%AquiferWater							(j,i) ,",",&			
!                            Me%HRU%WatBal%AquiferWaterOld						(j,i) ,",",&		
!                            Me%HRU%WatBal%DeepAquifer							(j,i) ,",",&			
!                            Me%HRU%WatBal%DeepAquiferOld						(j,i) ,",",&		
!                            Me%HRU%WatBal%SurfaceWaterColumn					(j,i) ,",",&	
!                            Me%HRU%WatBal%SurfaceWaterColumnOld					(j,i) ,",",&	
!                            Me%HRU%WatBal%LateralFlowStorage					(j,i) ,",",&	
!                            Me%HRU%WatBal%LateralFlowStorageOld					(j,i) ,",",&	
!							h2oloss													  ,",",&	
!							dstor
!
! 777    format (i5,a1,i5,a1,f8.2,a1,f8.2,a1,f8.2,&
!						 a1,f5.2,a1,f5.2,a1,f5.2,a1,f5.2,a1,f5.2,&
!						 a1,f5.2,a1,f5.2,a1,f5.2,a1,f5.2,a1,f5.2,&
!						 a1,f5.2,a1,f5.2,a1,f5.2,a1,f8.2,a1,f8.2,&
!						 a1,f5.2,a1,f5.2,a1,f8.2,a1,f8.2,a1,f5.2,&
!						 a1,f5.2,a1,f5.2,a1,f5.2,a1,f5.2,a1,f5.2)

							  DeepAqui =&
								Me%HRU%WatBal%DeepAquifer					(j,i)
							  SW =&
								Me%HRU%WatBal%SoilWater					(j,i)
					
							if (Monthly_out) then	
							
								year = year						

									write(222,888)	&
									j		,",",&
									i		,",",&
									year_old	,",",&
									Me%Month_Old	,",",&
									Hid_Year		,",",&
									Me%pcp (j,i)    ,",",&
									Me%Flow (j,i)   ,",",&
									Me%BaseFlow(j,i),",",&
									Me%ET  (j,i)    ,",",&
									Me%Round (j,i)  ,",",&
									Me%CountMonthDays (j,i)  ,",",&
									DeepAqui		,",",&
									SW              ,",",&
									Me%RiverFlow (i),",",&
									Me%WaterStress05 (j,i),",",&
									Me%WaterStress25 (j,i),",",&
									Me%WaterStress50 (j,i)

									Me%pcp      (j,i)= 0.0
									Me%Flow     (j,i)= 0.0
									Me%BaseFlow (j,i)= 0.0
									Me%ET       (j,i)= 0.0
									Me%Round    (j,i)= 0.0
									Me%CountMonthDays    (j,i)= 0.0
									Me%RiverFlow(i)= 0.0
									Me%WaterStress05 (j,i) = 0
									Me%WaterStress25 (j,i) = 0
									Me%WaterStress50 (j,i) = 0

							end if


							  Me%Flow (j,i) = Me%Flow (j,i) + & 
								Me%HRU%WatBal%RunOffFlowToChannel			(j,i) + &
								Me%HRU%WatBal%LateralFlowToChannel			(j,i) + &
								Me%HRU%WatBal%DrainageTileFlowToChannel		(j,i)
							  Me%ET (j,i) =   Me%ET (j,i) + &
								Me%HRU%WatBal%ActualEvapoTranspiration		(j,i) + &
								Me%HRU%WatBal%GrowndWaterFlowToAtmosphere	(j,i)
							  Me%pcp (j,i) =  Me%pcp (j,i) + &
								Me%HRU%WatBal%Precipitation					(j,i) - &
								Me%HRU%WatBal%SnowEvaporation				(j,i) 
							  Me%Round (j,i) = Me%Round (j,i) + &
								Me%HRU%WatBal%SoilFlowToRechargeFlow		      (j,i) - &  
								Me%HRU%WatBal%RechargeFlowToDeepAndShalowAquifer  (j,i) 	
							  Me%BaseFlow (j,i) = Me%BaseFlow (j,i) + &
								Me%HRU%WatBal%GrowndWaterFlowToChannel		(j,i) 
							  Me%CountMonthDays (j,i) = Me%CountMonthDays (j,i) + 1
							  ![hm3] = [hm3] + ([m3/s]/[day] * 86400 [s]/ 1 [day] ) * (1 [hm3] / 1e6 [m3])
							  Me%RiverFlow (i) = Me%RiverFlow (i)  + &
								Me%Reaches%FlowOut  (i) * 86400 * 1e-6

							  if(Me%HRU%Plant%FractionGrowth_StressWater   (j,i) < 0.05) then
									Me%WaterStress05 (j,i) = Me%WaterStress05 (j,i) + 1
							  end if

							  if(Me%HRU%Plant%FractionGrowth_StressWater   (j,i) < 0.25) then
									Me%WaterStress25 (j,i) = Me%WaterStress25 (j,i) + 1
							  end if

							  if(Me%HRU%Plant%FractionGrowth_StressWater   (j,i) < 0.5) then
									Me%WaterStress50 (j,i) = Me%WaterStress50 (j,i) + 1
							  end if


 888    format (i5,a1,i5,a1,f8.2,a1,f8.2,a1,f8.2,&
						 a1,f8.2,a1,f8.2,a1,f8.2,a1,f8.2,a1,f8.2,&
						 a1,f8.2, a1,f8.2, a1,f8.2, a1,e8.2,  &
						 a1,i2,a1,i2,a1,i2)

						end if

						if(Me%SwatOptions%NutrientPerHRU) then

							if (Monthly_out) then	
							
									write(333,8888)	&
									j					 ,",",&
									i					 ,",",&
									year_old			 ,",",&
									Me%Month_Old		 ,",",&
									Hid_Year			 ,",",&
									Me%N_Run_Off    (j,i),",",&
									Me%N_Percolated (j,i),",",&
									Me%P_Run_Off    (j,i),",",&
									Me%Sediment     (j,i)

									Me%N_Run_Off     (j,i)= 0.0
									Me%N_Percolated  (j,i)= 0.0
									Me%P_Run_Off	 (j,i)= 0.0
									Me%Sediment      (j,i)= 0.0

							end if

							  ![kg/ha] = [kg/ha] + [kg/ha]
							  Me%N_Run_Off (j,i) = Me%N_Run_Off (j,i) + & 
								Me%HRU%NTransportFlux%Nitrate%RunOff  (j,i) + &
								Me%HRU%NTransportFlux%RunOffOrganicN  (j,i)
							  Me%N_Percolated (j,i) =   Me%N_Percolated (j,i) + &
								Me%HRU%NTransportFlux%Nitrate%Percolation (j,i) + &
								Me%HRU%NTransportFlux%Nitrate%LateralFlow (j,i)
							  Me%P_Run_Off (j,i) =  Me%P_Run_Off (j,i) + &
								Me%HRU%PTransportFlux%RunOffActiveMineral (j,i) + &
								Me%HRU%PTransportFlux%RunOffStableMineral (j,i) + &
								Me%HRU%PTransportFlux%RunOffOrganicMatter (j,i)

							  ![ton/ha] = [ton/ha] + [ton/ha]
							  Me%Sediment (j,i) =  Me%Sediment (j,i) + &	
							    Me%HRU%Usle%Musle_Sed (j,i)

								!por isto para multiplas HRUs
								!* [km2] * [-] * [100ha/1km2] 
								!* Me%TotalBasinArea * Me%HRUBasinFraction (j,i) * 100 

 8888    format (i5,a1,i5,a1,f8.2,a1,f8.2,a1,f8.2,&
						 a1,f8.2,a1,f8.2,a1,f8.2,a1,f8.2)
						endif



                    endif    

                enddo
            enddo

            if (Me%SwatOptions%SubBasinOn) then
                do j=1,Me%MaxNumberSubBasin
                    if (Me%SubBasin%WaterPoints1D(j) == 1) then
                        Me%SubBasin%WaterFluxes%GlobalToChannel     (j) = OutputSubBasin (1,j)
                        Me%SubBasin%WaterFluxes%GroundWaterToChannel(j) = OutputSubBasin (2,j)
                        Me%SubBasin%WaterFluxes%RunOffToChannel     (j) = OutputSubBasin (3,j)
                        Me%SubBasin%WaterFluxes%LateralFlowToChannel(j) = OutputSubBasin (4,j)
                        !Me%Reaches%FlowOut                          (j) = OutputSubBasin (5,j)
                        Me%SubBasin%NitrateInRunOff                  (j) = OutputSubBasin (6,j)
                        Me%SubBasin%OrganicN%Dissolved              (j) = OutputSubBasin (7,j)
                        Me%SubBasin%OrganicN%Particulated           (j) = OutputSubBasin (8,j)
                        Me%SubBasin%OrganicP                        (j) = OutputSubBasin (9,j)
                        Me%SubBasin%MineralP                        (j) = OutputSubBasin (10,j)
                        Me%SubBasin%Sediments_concentration         (j) = OutputSubBasin (11,j)
                        Me%SubBasin%Sediments_load                  (j) = OutputSubBasin (12,j)
                        Me%SubBasin%WaterFluxes%LateralAndRunOffFlowToChannel (j) = &
                        Me%SubBasin%WaterFluxes%LateralFlowToChannel(j) + &
                        Me%SubBasin%WaterFluxes%RunOffToChannel(j)
					    !Me%SubBasin_Annual_Average%HRU%Meteo%SoilWater(1,j) = OutputHruMeteo (16,1,j)
                    end if
                enddo
            endif

            if (Me%SwatOptions%SubBasinOn) then
                if (Me%SwatOptions%BasinAverage) then
                    Me%Basin_Average = 0.0
                endif
                do j=1,Me%MaxNumberSubBasin
                    if (Me%SubBasin2%WaterPoints1D(j) == 1) then
                        Me%SubBasin2%WaterFluxes%GlobalToChannel     (j) = OutputSubBasin2 (1,j)
                        Me%SubBasin2%WaterFluxes%GroundWaterToChannel(j) = OutputSubBasin2 (2,j)
                        Me%SubBasin2%WaterFluxes%RunOffToChannel     (j) = OutputSubBasin2 (3,j)
                        Me%SubBasin2%WaterFluxes%LateralFlowToChannel(j) = OutputSubBasin2 (4,j)
                        !Me%Reaches%FlowOut                           (j) = OutputSubBasin (5,j)
                        Me%SubBasin2%NitrateInRunOff                 (j) = OutputSubBasin2 (6,j)
                        Me%SubBasin2%RunOffOrganicN%Dissolved        (j) = OutputSubBasin2 (7,j)
                        Me%SubBasin2%RunOffOrganicN%Particulated     (j) = OutputSubBasin2 (8,j)
                        Me%SubBasin2%RunOffSolubleP                  (j) = OutputSubBasin2 (9,j)
                        Me%SubBasin2%RunOffSedimentActiveP           (j) = OutputSubBasin2 (10,j)
                        Me%SubBasin2%Sediments                       (j) = OutputSubBasin2 (11,j)
                        Me%SubBasin2%WaterFluxes%LateralAndRunOffFlowToChannel (j) = &
                        Me%SubBasin2%WaterFluxes%LateralFlowToChannel(j) + &
                        Me%SubBasin2%WaterFluxes%RunOffToChannel(j)
                        Me%SubBasin2%NitrateInLateralFlow            (j) = OutputSubBasin2 (14,j)
                        Me%SubBasin2%NitrateInGrowndWaterFlow        (j) = OutputSubBasin2 (15,j)
                        Me%SubBasin2%RunOffSedimentStableP           (j) = OutputSubBasin2 (16,j)
                        Me%SubBasin2%GrowndWaterSolubleP             (j) = OutputSubBasin2 (17,j)
                        Me%SubBasin2%WaterFluxes%Precipitation                (j) = OutputSubBasin2 (18,j)
                        Me%SubBasin2%WaterFluxes%ActualEvapotranspiration     (j) = OutputSubBasin2 (19,j)
                        if (Me%OutPut%Yes) then
                            Me%SubBasin_Average%HRU%Meteo%SoilWater(1,j) = 0.0
                            Me%SubBasin_Average%HRU%Plant%LAI(1,j)       = 0.0
                            do i=1,Me%MaxNumberHRU
                                if(Me%HRU%WaterPoints2D(i,j) == 1) then
                                    HRUFraction = Me%HRUBasinFraction (i,j)  / Me%SubBasin%SubBasinFraction (j)
					                Me%SubBasin_Average%HRU%Meteo%SoilWater(1,j) = Me%SubBasin_Average%HRU%Meteo%SoilWater(1,j) + Me%HRU%Meteo%SoilWater (i,j) * HRUFraction
					                Me%SubBasin_Average%HRU%Plant%LAI(1,j)       = Me%SubBasin_Average%HRU%Plant%LAI(1,j) + Me%HRU%Plant%LAI(i,j)        * HRUFraction
                                end if
                            end do
                        endif
                        ![kgN/ha/year] = [kgN/ha/year] + [kgN/ha]/[year]
                        Me%SubBasin_Annual_Average%SubBasin2%NitrateInRunOff(j)   = Me%SubBasin_Annual_Average%SubBasin2%NitrateInRunOff(j)  + &
                                                                      Me%SubBasin2%NitrateInRunOff(j)  / &
                                                                      Me%NumberYearsSimulation
                        Me%SubBasin_Annual_Average%SubBasin2%RunOffOrganicN%Dissolved(j)   = Me%SubBasin_Annual_Average%SubBasin2%RunOffOrganicN%Dissolved(j)  + &
                                                                      Me%SubBasin2%RunOffOrganicN%Dissolved(j)  / &
                                                                      Me%NumberYearsSimulation
                        Me%SubBasin_Annual_Average%SubBasin2%RunOffOrganicN%Particulated(j)   = Me%SubBasin_Annual_Average%SubBasin2%RunOffOrganicN%Particulated(j)  + &
                                                                      Me%SubBasin2%RunOffOrganicN%Particulated(j)  / &
                                                                      Me%NumberYearsSimulation
                        Me%SubBasin_Annual_Average%SubBasin2%NitrateInLateralFlow(j)   = Me%SubBasin_Annual_Average%SubBasin2%NitrateInLateralFlow(j)  + &
                                                                      Me%SubBasin2%NitrateInLateralFlow(j)  / &
                                                                      Me%NumberYearsSimulation
                        Me%SubBasin_Annual_Average%SubBasin2%NitrateInGrowndWaterFlow(j)   = Me%SubBasin_Annual_Average%SubBasin2%NitrateInGrowndWaterFlow(j)  + &
                                                                      Me%SubBasin2%NitrateInGrowndWaterFlow(j)  / &
                                                                      Me%NumberYearsSimulation
                        Me%SubBasin_Annual_Average%SubBasin2%RunOffSolubleP(j)   = Me%SubBasin_Annual_Average%SubBasin2%RunOffSolubleP(j)  + &
                                                                      Me%SubBasin2%RunOffSolubleP(j)  / &
                                                                      Me%NumberYearsSimulation
                        Me%SubBasin_Annual_Average%SubBasin2%RunOffSedimentActiveP(j)   = Me%SubBasin_Annual_Average%SubBasin2%RunOffSedimentActiveP(j)  + &
                                                                      Me%SubBasin2%RunOffSedimentActiveP(j)  / &
                                                                      Me%NumberYearsSimulation
                        Me%SubBasin_Annual_Average%SubBasin2%RunOffSedimentStableP(j)   = Me%SubBasin_Annual_Average%SubBasin2%RunOffSedimentStableP(j)  + &
                                                                      Me%SubBasin2%RunOffSedimentStableP(j)  / &
                                                                      Me%NumberYearsSimulation
                        Me%SubBasin_Annual_Average%SubBasin2%GrowndWaterSolubleP(j)   = Me%SubBasin_Annual_Average%SubBasin2%GrowndWaterSolubleP(j)  + &
                                                                      Me%SubBasin2%GrowndWaterSolubleP(j)  / &
                                                                      Me%NumberYearsSimulation

                        Me%SubBasin_Annual_Average%SubBasin2%WaterFluxes%GlobalToChannel(j)   = Me%SubBasin_Annual_Average%SubBasin2%WaterFluxes%GlobalToChannel(j)  + &
                                                                      Me%SubBasin2%WaterFluxes%GlobalToChannel(j)  / &
                                                                      Me%NumberYearsSimulation

                        Me%SubBasin_Annual_Average%SubBasin2%WaterFluxes%Precipitation(j)   = Me%SubBasin_Annual_Average%SubBasin2%WaterFluxes%Precipitation(j)  + &
                                                                      Me%SubBasin2%WaterFluxes%Precipitation(j)  / &
                                                                      Me%NumberYearsSimulation

                        Me%SubBasin_Annual_Average%SubBasin2%WaterFluxes%ActualEvapotranspiration(j)   = Me%SubBasin_Annual_Average%SubBasin2%WaterFluxes%ActualEvapotranspiration(j)  + &
                                                                      Me%SubBasin2%WaterFluxes%ActualEvapotranspiration(j)  / &
                                                                      Me%NumberYearsSimulation


                        if (Me%SwatOptions%BasinAverage) then
                                ![kgN/ha] = [kgN/ha] + [kgN/ha] * [-]
                                Me%Basin_Average (1) = Me%Basin_Average (1) + Me%SubBasin2%RunOffSolubleP               (j) &
                                                                            * Me%SubBasin%SubBasinFraction (j) 
                                Me%Basin_Average (2) = Me%Basin_Average (2) + Me%SubBasin2%RunOffSedimentStableP        (j) &
                                                                            * Me%SubBasin%SubBasinFraction (j) 
                                Me%Basin_Average (3) = Me%Basin_Average (3) + Me%SubBasin2%RunOffSedimentActiveP        (j) &
                                                                            * Me%SubBasin%SubBasinFraction (j) 
                                Me%Basin_Average (4) = Me%Basin_Average (4) + Me%SubBasin2%GrowndWaterSolubleP          (j) &
                                                                            * Me%SubBasin%SubBasinFraction (j) 
                                Me%Basin_Average (5) = Me%Basin_Average (5) + Me%SubBasin2%NitrateInRunOff              (j) &
                                                                            * Me%SubBasin%SubBasinFraction (j) 
                                Me%Basin_Average (6) = Me%Basin_Average (6) + Me%SubBasin2%RunOffOrganicN%Dissolved     (j) &
                                                                            * Me%SubBasin%SubBasinFraction (j)
                                Me%Basin_Average (7) = Me%Basin_Average (7) + Me%SubBasin2%RunOffOrganicN%Particulated  (j) &
                                                                            * Me%SubBasin%SubBasinFraction (j) 
                                Me%Basin_Average (8) = Me%Basin_Average (8) + Me%SubBasin2%NitrateInLateralFlow         (j) &
                                                                            * Me%SubBasin%SubBasinFraction (j) 
                                Me%Basin_Average (9) = Me%Basin_Average (9) + Me%SubBasin2%NitrateInGrowndWaterFlow     (j) &
                                                                            * Me%SubBasin%SubBasinFraction (j) 
                        endif
 

                    end if
                enddo
            endif


            if (Me%SwatOptions%GlobalWaterOn) then

                if (Me%SwatOptions%GlobalWaterM3) then
                    ![m3/mm]            = [m/1e3mm] * [km2] * [1e6m2/km2]
                    ConversionFactor = 1e-3 * Me%TotalBasinArea * 1e6
                else
                    ![mm/mm]
                    ConversionFactor = 1.0
                endif
                
                if (Me%SwatOptions%GlobalWaterCumulative) then

                    Me%GlobalWaterFluxes (1) = Me%GlobalWaterFluxes (1) + OutputGlobalBasin(1)  * ConversionFactor !precipitation mm H2O'          
                    Me%GlobalWaterFluxes (2) = Me%GlobalWaterFluxes (2) + OutputGlobalBasin(3)  * ConversionFactor !runoff mm H2O'           
                    Me%GlobalWaterFluxes (3) = Me%GlobalWaterFluxes (3) + OutputGlobalBasin(4)  * ConversionFactor !lateral_flow mm H2O'           
                    Me%GlobalWaterFluxes (4) = Me%GlobalWaterFluxes (4) + OutputGlobalBasin(5)  * ConversionFactor !percolation mm H2O'            
                    Me%GlobalWaterFluxes (5) = Me%GlobalWaterFluxes (5) + OutputGlobalBasin(6)  * ConversionFactor !water_yield mm H2O'       
                    Me%GlobalWaterFluxes (6) = Me%GlobalWaterFluxes (6) + OutputGlobalBasin(7)  * ConversionFactor !act_evapotranspiration mm H2O'
                    Me%GlobalWaterFluxes (7) =                            OutputGlobalBasin(35) * ConversionFactor !soil_stored_water mm H2O' 
                    Me%GlobalWaterFluxes (8) = Me%GlobalWaterFluxes (8) + OutputGlobalBasin(39) * ConversionFactor !snow_fall_freezing_rain mm H2O'
                    Me%GlobalWaterFluxes (9) = Me%GlobalWaterFluxes (9) + OutputGlobalBasin(104)* ConversionFactor !groundwater mm H2O'
                    Me%GlobalWaterFluxes (10)= Me%GlobalWaterFluxes (10)+ OutputGlobalBasin(106)* ConversionFactor !deep_aquifer_recharge mm H2O   
                    Me%GlobalWaterFluxes (11)= Me%GlobalWaterFluxes (11)+ OutputGlobalBasin(108)* ConversionFactor !pot_evapotranspiration mm H2O' 
                    Me%GlobalWaterFluxes (12)= Me%GlobalWaterFluxes (12)+ OutputGlobalBasin(109)* ConversionFactor !drainage_tile mm H2O'    
                else
                    Me%GlobalWaterFluxes (1) = OutputGlobalBasin(1)  * ConversionFactor !precipitation mm H2O'          
                    Me%GlobalWaterFluxes (2) = OutputGlobalBasin(3)  * ConversionFactor !runoff mm H2O'           
                    Me%GlobalWaterFluxes (3) = OutputGlobalBasin(4)  * ConversionFactor !lateral_flow mm H2O'           
                    Me%GlobalWaterFluxes (4) = OutputGlobalBasin(5)  * ConversionFactor !percolation mm H2O'            
                    Me%GlobalWaterFluxes (5) = OutputGlobalBasin(6)  * ConversionFactor !water_yield mm H2O'       
                    Me%GlobalWaterFluxes (6) = OutputGlobalBasin(7)  * ConversionFactor !act_evapotranspiration mm H2O'
                    Me%GlobalWaterFluxes (7) = OutputGlobalBasin(35) * ConversionFactor !soil_stored_water mm H2O'
                    Me%GlobalWaterFluxes (8) = OutputGlobalBasin(39) * ConversionFactor !snow_fall_freezing_rain mm H2O'
                    Me%GlobalWaterFluxes (9) = OutputGlobalBasin(104)* ConversionFactor !groundwater mm H2O'
                    Me%GlobalWaterFluxes (10)= OutputGlobalBasin(106)* ConversionFactor !deep_aquifer_recharge mm H2O            
                    Me%GlobalWaterFluxes (11)= OutputGlobalBasin(108)* ConversionFactor !pot_evapotranspiration mm H2O' 
                    Me%GlobalWaterFluxes (12)= OutputGlobalBasin(109)* ConversionFactor !drainage_tile mm H2O'    
                endif

            endif

            if (Me%SwatOptions%GlobalNutrientsOn) then
                Me%GlobalNutrients(1) = OutputGlobalBasin(12)!sediment yield from HRUs TONS sed
                Me%GlobalNutrients(2) = OutputGlobalBasin(13)!sediment loading to ponds TONS sed
                Me%GlobalNutrients(3) = OutputGlobalBasin(14)!sediment loading from ponds TONS sed
                Me%GlobalNutrients(4) = OutputGlobalBasin(15)!net change in sediment level in ponds TONS sed
                Me%GlobalNutrients(5) = OutputGlobalBasin(16)!sediment loading to wetlands TONS sed
                Me%GlobalNutrients(6) = OutputGlobalBasin(17)!sediment loading to main channels from wetlands TONS sed
                Me%GlobalNutrients(7) = OutputGlobalBasin(18)!net change in sediment in wetlands TONS sed

                Me%GlobalNutrients(8) = 0.0

                do j=1,Me%MaxNumberHRU
                    do i=1,Me%MaxNumberSubBasin
                        if (Me%HRU%WaterPoints2D(j,i) == 1) then
                            !N added with agriculture [ton] = ([kg/ha] + [kg/ha] + [kg/ha]) 
                            !                            * [km2] * [-] * [100ha/1km2] * [1ton/1000kg]
                            Me%GlobalNutrients(8) = Me%GlobalNutrients(8) + (OutputHruQuality (13,j,i)    & !autofert
                                                                          +  OutputHruQuality (14,j,i)    & !fert
                                                                          +  OutputHruQuality (17,j,i) )  & !grazing
                                                                          * Me%TotalBasinArea * Me%HRUBasinFraction (j,i) &
                                                                          * 100 * 1e-3
                        endif
                    enddo
                enddo


            endif

            if(Me%SwatOptions%LandUse) then
                if(Me%HRU%AuxLogicalLandUse) then
                    Me%HRU%AuxLogicalLandUse = .false.
                    do j=1,Me%MaxNumberHRU
                        do i=1,Me%MaxNumberSubBasin
				            do k=1,Me%HRU%NumberOfLandUseOutputs
					            if (Me%HRU%LandUse(k) == OutputHRU_Cropname(j,i).and.Me%HRU%WaterPoints2D(j,i) == 1) then
                                    !                [ha] = [km2] * [-] * [100ha/1km2]
                                    Me%HRU%LandUseAreaHectare(k) = Me%HRU%LandUseAreaHectare(k) + &
                                                            Me%TotalBasinArea * Me%HRUBasinFraction (j,i) * 100
                                endif
                            end do
                        end do
                    end do
                endif

			    Me%HRU%LandUseNutrient = 0.0
                do j=1,Me%MaxNumberHRU
                    do i=1,Me%MaxNumberSubBasin
					    do k=1,Me%HRU%NumberOfLandUseOutputs
						    if (Me%HRU%LandUse(k) == OutputHRU_Cropname(j,i).and.Me%HRU%WaterPoints2D(j,i) == 1) then
                                ![-] = [km2] * [-] * [100ha/1km2] / [ha]
                                Aux = Me%TotalBasinArea * Me%HRUBasinFraction (j,i) * 100 / Me%HRU%LandUseAreaHectare(k)
							    do n = 1,68
                            
                                !         [kg/ha]               =  [kg/ha] + [kg/ha]  * [-]
								    Me%HRU%LandUseNutrient(n,k) = Me%HRU%LandUseNutrient(n,k) + Me%HRU%HRU_output_values (n,j,i) * Aux 
							    end do

                                Me%HRU%LandUseNutrient(69,k) = Me%HRU%LandUseNutrient(69,k) + Me%HRU%NTransportFlux%Volatilization (j,i) * Aux 
                                Me%HRU%LandUseNutrient(70,k) = Me%HRU%LandUseNutrient(70,k) + Me%HRU%RateFlux%FreshOrgNToNitrate   (j,i) * Aux 
                                Me%HRU%LandUseNutrient(71,k) = Me%HRU%LandUseNutrient(71,k) + Me%HRU%PoolsN%Nitrate                 (j,i) * Aux 
                                Me%HRU%LandUseNutrient(72,k) = Me%HRU%LandUseNutrient(72,k) + Me%HRU%PoolsN%ActiveOrganicMatter     (j,i) * Aux 
                                Me%HRU%LandUseNutrient(73,k) = Me%HRU%LandUseNutrient(73,k) + Me%HRU%PoolsN%StableOrganicMatter     (j,i) * Aux 
                                Me%HRU%LandUseNutrient(74,k) = Me%HRU%LandUseNutrient(74,k) + Me%HRU%PoolsN%FreshOrganicMatter      (j,i) * Aux 
                                Me%HRU%LandUseNutrient(75,k) = Me%HRU%LandUseNutrient(75,k) + Me%HRU%PoolsN%Ammonium                (j,i) * Aux 
                                Me%HRU%LandUseNutrient(76,k) = Me%HRU%LandUseNutrient(76,k) + Me%HRU%NTransportFlux%Nitrification  (j,i) * Aux 
                                Me%HRU%LandUseNutrient(77,k) = Me%HRU%LandUseNutrient(77,k) + Me%HRU%PoolsP%ActiveMineral           (j,i) * Aux 
                                Me%HRU%LandUseNutrient(78,k) = Me%HRU%LandUseNutrient(78,k) + Me%HRU%PoolsP%FreshOrganicMatter      (j,i) * Aux 
                                Me%HRU%LandUseNutrient(79,k) = Me%HRU%LandUseNutrient(79,k) + Me%HRU%PoolsP%OrganicMatter           (j,i) * Aux 
                                Me%HRU%LandUseNutrient(80,k) = Me%HRU%LandUseNutrient(80,k) + Me%HRU%PoolsP%StableMineral           (j,i) * Aux 
						    endif
					    enddo

                    enddo
                enddo
            endif

            call OutPut_TimeSeries


			if(Me%OutPut%Yes) call HDF5Output

            write(2,*) Year, Me%Month, Day, iyr,ida

            call ActualizeCurrentTime(Me%ObjTime,           &
                                      DT_Global = 43200.0,   &
                                      STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                              &
                stop 'ModifySwat - ModuleSwat - ERR01'

            call ActualizeCurrentTime(Me%ObjTime,           &
                                      DT_Global = 43200.0,   &
                                      STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                              &
                stop 'ModifySwat - ModuleSwat - ERR01'

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifySwat
    !--------------------------------------------------------------------------
    
    subroutine HDF5Output

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: j
        real, dimension(6)  , target                :: AuxTime
        real, dimension(:)  , pointer               :: TimePointer
        real, dimension(:)  , pointer               :: Dummy       

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        
        allocate (dummy(Me%MaxNumberSubBasin))

        if (Me%CurrentTime >= Me%OutPut%OutTime(Me%OutPut%NextOutPut)) then

            !Writes current time
            call ExtractDate   (Me%CurrentTime , AuxTime(1), AuxTime(2),         &
                                                AuxTime(3), AuxTime(4),         &
                                                AuxTime(5), AuxTime(6))
            TimePointer => AuxTime

			if (Me%SwatOptions%ReachesOn) Then
				call HDF5SetLimits  (Me%ObjHDF5_Reaches, 1, 6, STAT = STAT_CALL)
				if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR01'

				call HDF5WriteData  (Me%ObjHDF5_Reaches, "/Time", "Time",                   &
									 "YYYY/MM/DD HH:MM:SS",                         &
									 Array1D      = TimePointer,                    &
									 OutputNumber = Me%OutPut%NextOutPut,           &
									 STAT = STAT_CALL)
				if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR02'

				!Sets limits for next write operations
				call HDF5SetLimits   (Me%ObjHDF5_Reaches, 1, Me%MaxNumberSubBasin, STAT = STAT_CALL)
				if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR03'


				!Writes Flow values

				!Writes Flow X
				call HDF5WriteData   (Me%ObjHDF5_Reaches, "//Results/Flow",          &
									  "Flow", "m3/s",                              &
									  Array1D      = Me%Reaches%FlowOut,            &
									  OutputNumber = Me%OutPut%NextOutPut,          &
									  STAT = STAT_CALL)
				if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR09'
			end if

			if (Me%SwatOptions%SubBasinOn) then
				call HDF5SetLimits  (Me%ObjHDF5_SubBasin, 1, 6, STAT = STAT_CALL)
				if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR01'

				call HDF5WriteData  (Me%ObjHDF5_SubBasin, "/Time", "Time",                   &
									 "YYYY/MM/DD HH:MM:SS",                         &
									 Array1D      = TimePointer,                    &
									 OutputNumber = Me%OutPut%NextOutPut,           &
									 STAT = STAT_CALL)
				if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR02'

				!Sets limits for next write operations
				call HDF5SetLimits   (Me%ObjHDF5_SubBasin, 1, Me%MaxNumberSubBasin, STAT = STAT_CALL)
				if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR03'


				!Writes erosion
				call HDF5WriteData   (Me%ObjHDF5_SubBasin, "//Results/erosion",     &
									  "erosion", "ton/ha",                          &
									  Array1D      = Me%SubBasin%Sediments_load,    &
									  OutputNumber = Me%OutPut%NextOutPut,          &
									  STAT = STAT_CALL)
				if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR09'

                do j=1,Me%MaxNumberSubBasin
                    if (Me%SubBasin%WaterPoints1D(j) == 1) then
					    dummy (j) = Me%SubBasin_Average%HRU%Meteo%SoilWater(1,j)
                    end if
                enddo


				!Writes soil water 
				call HDF5WriteData   (Me%ObjHDF5_SubBasin, "//Results/SoilWater",     &
									  "SoilWater", "mm",                          &
									  Array1D      = dummy,    &
									  OutputNumber = Me%OutPut%NextOutPut,          &
									  STAT = STAT_CALL)
				if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR09'

                do j=1,Me%MaxNumberSubBasin
                    if (Me%SubBasin%WaterPoints1D(j) == 1) then
					    dummy (j) = Me%SubBasin_Average%HRU%Plant%LAI(1,j)
                    end if
                enddo


				!Writes soil water 
				call HDF5WriteData   (Me%ObjHDF5_SubBasin, "//Results/LAI",     &
									  "LAI", "m2/m2",                          &
									  Array1D      = dummy,    &
									  OutputNumber = Me%OutPut%NextOutPut,          &
									  STAT = STAT_CALL)
				if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR09'

			end if
          
            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5_Reaches, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR99'

            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5_SubBasin, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR99'

            Me%OutPut%NextOutPut = Me%OutPut%NextOutPut + 1

        endif
        
        deallocate (dummy)

    end subroutine HDF5Output
    !--------------------------------------------------------------------------

    subroutine OutPut_TimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        if (Me%SwatOptions%MeteoOutputOn) call MeteoTimeSeries

        if (Me%SwatOptions%PlantOutputOn) call PlantTimeSeries

        if (Me%SwatOptions%OperationalOutputOn) call OperationalTimeSeries

        if (Me%SwatOptions%UsleOutputOn) call UsleTimeSeries

        if (Me%SwatOptions%NitrateOutputOn) call NitrateTimeSeries

        if (Me%SwatOptions%PhosphorusOutputOn) call PhosphorusTimeSeries

        if (Me%SwatOptions%WaterBalanceOn) call WaterBalanceTimeSeries

        if (Me%SwatOptions%SubBasinOn) call SubBasinTimeSeries

        if (Me%SwatOptions%ReachesOn) call ReachesTimeSeries

		if (Me%SwatOptions%LandUse) call LandUseTimeSeries

        if (Me%SwatOptions%GlobalWaterOn) then
            call WriteTimeSerieLine (TimeSerieID = Me%ObjGlobalWaterTimeSerie,& 
                                     DataLine    = Me%GlobalWaterFluxes,      &
                                     STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                        &
                stop 'OutPut_TimeSeries - ModuleSwat - ERR01'
        endif

        if (Me%SwatOptions%GlobalNutrientsOn) then
            call WriteTimeSerieLine (TimeSerieID = Me%ObjGlobalNutrientTimeSerie,& 
                                     DataLine    = Me%GlobalNutrients,      &
                                     STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                        &
                stop 'OutPut_TimeSeries - ModuleSwat - ERR01'
        endif

        if (Me%SwatOptions%BasinAverage.and.Me%SwatOptions%SubBasinOn) then
            call WriteTimeSerieLine (TimeSerieID = Me%ObjBasinAverage,& 
                                     DataLine    = Me%Basin_Average,      &
                                     STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                        &
                stop 'OutPut_TimeSeries - ModuleSwat - ERR03'
        endif


    end subroutine OutPut_TimeSeries
    !--------------------------------------------------------------------------

    subroutine MeteoTimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%AverageTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%Precipitation, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR02'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%PotencialEvapoTranspiration, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR03'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%CulturalTranspiration      , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR04'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%ActualEvapoTranspiration, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR04'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%MinimumTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%MaximumTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR02'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%WindSpeed, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR03'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%RelativeHumidity, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR04'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%SolarRadiation, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR05'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%Irrigation, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR06'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%PrecipitationAsSnow, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR07'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%SnowMelt, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR08'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%SnowPack_mmH2O, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR09'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%SoilTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR09'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%SoilWater, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR10'
        call WriteTimeSerie(Me%ObjMeteoTimeSerie, &
            Data2D = Me%HRU%Meteo%AquiferWater, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR11'

    end subroutine MeteoTimeSeries

    !--------------------------------------------------------------------------

    subroutine PlantTimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------
        call WriteTimeSerie(Me%ObjPlantTimeSerie, &
            Data2D = Me%HRU%Plant%Biomass, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjPlantTimeSerie, &
            Data2D = Me%HRU%Plant%NitrogenFraction, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR02'
        call WriteTimeSerie(Me%ObjPlantTimeSerie, &
            Data2D = Me%HRU%Plant%NitrogenYield, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR03'
        call WriteTimeSerie(Me%ObjPlantTimeSerie, &
            Data2D = Me%HRU%Plant%FractionOfHeatUnits, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR04'
        call WriteTimeSerie(Me%ObjPlantTimeSerie, &
            Data2D = Me%HRU%Plant%LAI                , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR05'
        call WriteTimeSerie(Me%ObjPlantTimeSerie, &
            Data2D = Me%HRU%Plant%FractionGrowth_StressP, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR02'
        call WriteTimeSerie(Me%ObjPlantTimeSerie, &
            Data2D = Me%HRU%Plant%FractionGrowth_StressN, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR03'
        call WriteTimeSerie(Me%ObjPlantTimeSerie, &
            Data2D = Me%HRU%Plant%FractionGrowth_StressTemp, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR04'
        call WriteTimeSerie(Me%ObjPlantTimeSerie, &
            Data2D = Me%HRU%Plant%FractionGrowth_StressWater                , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR05'

    end subroutine PlantTimeSeries
    !--------------------------------------------------------------------------
    subroutine OperationalTimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%AverageTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%Precipitation, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR02'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%PotencialEvapoTranspiration, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR03'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%CulturalTranspiration      , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR04'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%ActualEvapoTranspiration, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR04'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%MinimumTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%MaximumTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR02'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%WindSpeed, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR03'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%RelativeHumidity, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR04'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%SolarRadiation, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR05'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%Irrigation, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR06'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%PrecipitationAsSnow, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR07'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%SnowMelt, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR08'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%SnowPack_mmH2O, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR09'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%SoilTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR09'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%SoilWater, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR10'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Meteo%AquiferWater, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MeteoTimeSeries - ModuleSwat - ERR11'

        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Plant%Biomass, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Plant%NitrogenFraction, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR02'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Plant%NitrogenYield, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR03'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Plant%FractionOfHeatUnits, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR04'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Plant%LAI                , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR05'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Plant%FractionGrowth_StressP, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR02'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Plant%FractionGrowth_StressN, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR03'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Plant%FractionGrowth_StressTemp, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR04'
        call WriteTimeSerie(Me%ObjOperationalTimeSerie, &
            Data2D = Me%HRU%Plant%FractionGrowth_StressWater                , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PlantTimeSeries - ModuleSwat - ERR05'

    end subroutine OperationalTimeSeries
    !--------------------------------------------------------------------------
    subroutine UsleTimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------
        call WriteTimeSerie(Me%ObjUsleTimeSerie, &
            Data2D = Me%HRU%Usle%Musle_Sed, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UsleTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjUsleTimeSerie, &
            Data2D = Me%HRU%Usle%Usle_Sed, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UsleTimeSeries - ModuleSwat - ERR02'
        call WriteTimeSerie(Me%ObjUsleTimeSerie, &
            Data2D = Me%HRU%Usle%C_K_P_LS_rock, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UsleTimeSeries - ModuleSwat - ERR03'
        call WriteTimeSerie(Me%ObjUsleTimeSerie, &
            Data2D = Me%HRU%Usle%MuslePeakRunoffRate, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UsleTimeSeries - ModuleSwat - ERR04'
        call WriteTimeSerie(Me%ObjUsleTimeSerie, &
            Data2D = Me%HRU%Usle%UsleRainfallIndex, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UsleTimeSeries - ModuleSwat - ERR04'
        call WriteTimeSerie(Me%ObjUsleTimeSerie, &
            Data2D = Me%HRU%Usle%K_LS_P_rock, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'UsleTimeSeries - ModuleSwat - ERR04'

    end subroutine UsleTimeSeries
    !--------------------------------------------------------------------------

    subroutine NitrateTimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%RateFlux%ActiveOrgNToNitrate        , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%RateFlux%FreshOrgNToNitrate        , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%RateFlux%Denitrification           , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%NTransportFlux%Nitrate%LateralFlow , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%NTransportFlux%Nitrate%Rain        , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%NTransportFlux%Nitrate%PlantUptake , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%NTransportFlux%Nitrate%Percolation , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%NTransportFlux%Nitrate%RunOff      , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%NTransportFlux%RunOffOrganicN      , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%NTransportFlux%AutoFertilization   , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%NTransportFlux%Fertilization       , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%PoolsN%Nitrate                      , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%PoolsN%FreshOrganicMatter           , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%PoolsN%ActiveOrganicMatter           , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%PoolsN%StableOrganicMatter           , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%NTransportFlux%Grazing             , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%NTransportFlux%Volatilization      , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%GW%NitrateConcentration            , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'
        call WriteTimeSerie(Me%ObjNitrateTimeSerie, &
            Data2D = Me%HRU%PoolsN%SoilResidue            , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'NitrateTimeSeries - ModuleSwat - ERR01'

    end subroutine NitrateTimeSeries

    !--------------------------------------------------------------------------

    subroutine PhosphorusTimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        call WriteTimeSerie(Me%ObjPhosphorusTimeSerie, &
            Data2D = Me%HRU%PoolsP%ActiveMineral                      , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PhosphorusTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjPhosphorusTimeSerie, &
           Data2D = Me%HRU%PoolsP%FreshOrganicMatter                 , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PhosphorusTimeSeries - ModuleSwat - ERR02'

        call WriteTimeSerie(Me%ObjPhosphorusTimeSerie, &
            Data2D = Me%HRU%PoolsP%OrganicMatter                      , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PhosphorusTimeSeries - ModuleSwat - ERR03'

        call WriteTimeSerie(Me%ObjPhosphorusTimeSerie, &
            Data2D = Me%HRU%PoolsP%StableMineral                      , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PhosphorusTimeSeries - ModuleSwat - ERR04'

        call WriteTimeSerie(Me%ObjPhosphorusTimeSerie, &
            Data2D = Me%HRU%PTransportFlux%AutoFertilization          , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PhosphorusTimeSeries - ModuleSwat - ERR05'

        call WriteTimeSerie(Me%ObjPhosphorusTimeSerie, &
            Data2D = Me%HRU%PTransportFlux%Fertilization              , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PhosphorusTimeSeries - ModuleSwat - ERR06'

        call WriteTimeSerie(Me%ObjPhosphorusTimeSerie, &
            Data2D = Me%HRU%PTransportFlux%Grazing                    , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PhosphorusTimeSeries - ModuleSwat - ERR07'

        call WriteTimeSerie(Me%ObjPhosphorusTimeSerie, &
            Data2D = Me%HRU%PTransportFlux%RunOffActiveMineral        , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PhosphorusTimeSeries - ModuleSwat - ERR08'

        call WriteTimeSerie(Me%ObjPhosphorusTimeSerie, &
            Data2D = Me%HRU%PTransportFlux%RunOffStableMineral        , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PhosphorusTimeSeries - ModuleSwat - ERR09'

        call WriteTimeSerie(Me%ObjPhosphorusTimeSerie, &
            Data2D = Me%HRU%PTransportFlux%RunOffOrganicMatter        , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PhosphorusTimeSeries - ModuleSwat - ERR10'

        call WriteTimeSerie(Me%ObjPhosphorusTimeSerie, &
            Data2D = Me%HRU%PTransportFlux%PlantUptake                , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PhosphorusTimeSeries - ModuleSwat - ERR11'

    end subroutine PhosphorusTimeSeries

    !--------------------------------------------------------------------------
    subroutine WaterBalanceTimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%SnowPack_mmH2O						, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%SnowPack_mmH2O_Old					, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%SoilWater							, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%SoilWaterOld						, STAT = STAT_CALL)	
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%AquiferWater						, STAT = STAT_CALL)	
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%AquiferWaterOld						, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%DeepAquifer							, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%DeepAquiferOld						, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%SurfaceWaterColumn					, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%SurfaceWaterColumnOld				, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%LateralFlowStorage					, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%LateralFlowStorageOld				, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'


        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%Precipitation						, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%SnowEvaporation						, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%RunOffFlowToChannel					, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%LateralFlowToChannel				, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%ActualEvapoTranspiration			, STAT = STAT_CALL)	
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%GrowndWaterFlowToChannel			, STAT = STAT_CALL)	
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%GrowndWaterFlowToAtmosphere			, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%PondWaterLoss						, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%WetlandWaterLoss					, STAT = STAT_CALL)	
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%Irrigation							, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%RechargeFlowToDeepAndShalowAquifer	, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%DrainageTileFlowToChannel			, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%SoilFlowToRechargeFlow			, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'

        call WriteTimeSerie(Me%ObjWaterBalanceTimeSerie, &
            Data2D = Me%HRU%WatBal%WaterBalanceError			, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterBalanceTimeSeries - ModuleSwat - ERR01'



    end subroutine WaterBalanceTimeSeries

    !--------------------------------------------------------------------------

    subroutine SubBasinTimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        call WriteTimeSerie (Me%ObjSubBasinTimeSerie, Data1D = Me%SubBasin%WaterFluxes%GlobalToChannel      )
        call WriteTimeSerie (Me%ObjSubBasinTimeSerie, Data1D = Me%SubBasin%WaterFluxes%GroundWaterToChannel )
        call WriteTimeSerie (Me%ObjSubBasinTimeSerie, Data1D = Me%SubBasin%WaterFluxes%RunOffToChannel      )
        call WriteTimeSerie (Me%ObjSubBasinTimeSerie, Data1D = Me%SubBasin%WaterFluxes%LateralFlowToChannel )
        call WriteTimeSerie (Me%ObjSubBasinTimeSerie, Data1D = Me%Reaches%FlowOut                           )
        call WriteTimeSerie (Me%ObjSubBasinTimeSerie, Data1D = Me%SubBasin%NitrateInRunOff                  )
        call WriteTimeSerie (Me%ObjSubBasinTimeSerie, Data1D = Me%SubBasin%OrganicN%Dissolved               )
        call WriteTimeSerie (Me%ObjSubBasinTimeSerie, Data1D = Me%SubBasin%OrganicN%Particulated            )
        call WriteTimeSerie (Me%ObjSubBasinTimeSerie, Data1D = Me%SubBasin%OrganicP                         )
        call WriteTimeSerie (Me%ObjSubBasinTimeSerie, Data1D = Me%SubBasin%MineralP                         )
        call WriteTimeSerie (Me%ObjSubBasinTimeSerie, Data1D = Me%SubBasin%Sediments_concentration          )
        call WriteTimeSerie (Me%ObjSubBasinTimeSerie, Data1D = Me%SubBasin%WaterFluxes%LateralAndRunOffFlowToChannel)

    end subroutine SubBasinTimeSeries
    !--------------------------------------------------------------------------
    subroutine ReachesTimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        call WriteTimeSerie (Me%ObjReachesTimeSerie, Data1D = Me%Reaches%FlowOut             )
        call WriteTimeSerie (Me%ObjReachesTimeSerie, Data1D = Me%Reaches%algal_biomass       )
        call WriteTimeSerie (Me%ObjReachesTimeSerie, Data1D = Me%Reaches%ammonia             )
        call WriteTimeSerie (Me%ObjReachesTimeSerie, Data1D = Me%Reaches%BOD                 )
        call WriteTimeSerie (Me%ObjReachesTimeSerie, Data1D = Me%Reaches%chlorophyll_a       )
        call WriteTimeSerie (Me%ObjReachesTimeSerie, Data1D = Me%Reaches%dissolved_phosphorus)
        call WriteTimeSerie (Me%ObjReachesTimeSerie, Data1D = Me%Reaches%dissolved_oxygen    )
        call WriteTimeSerie (Me%ObjReachesTimeSerie, Data1D = Me%Reaches%nitrate             )
        call WriteTimeSerie (Me%ObjReachesTimeSerie, Data1D = Me%Reaches%nitrite             )
        call WriteTimeSerie (Me%ObjReachesTimeSerie, Data1D = Me%Reaches%organic_nitrogen    )
        call WriteTimeSerie (Me%ObjReachesTimeSerie, Data1D = Me%Reaches%organic_phosphorus  )
        call WriteTimeSerie (Me%ObjReachesTimeSerie, Data1D = Me%Reaches%saturation_oxygen   )
        call WriteTimeSerie (Me%ObjReachesTimeSerie, Data1D = Me%Reaches%temperature         )
        call WriteTimeSerie (Me%ObjReachesTimeSerie, Data1D = Me%Reaches%Sediments           )
         
    end subroutine ReachesTimeSeries
    !--------------------------------------------------------------------------

    subroutine LandUseTimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
		real,dimension(:  ), pointer            :: dummy
		integer									:: i,k

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

		allocate(dummy(1:Me%HRU%NumberOfLandUseOutputs))

		do i=31,59

			do k=1,Me%HRU%NumberOfLandUseOutputs
				dummy(k) = Me%HRU%LandUseNutrient(i,k)
			enddo

			call WriteTimeSerie (Me%ObjLandUseTimeSerie, Data1D = dummy)

		enddo

		do i=69,80

			do k=1,Me%HRU%NumberOfLandUseOutputs
				dummy(k) = Me%HRU%LandUseNutrient(i,k)
			enddo

			call WriteTimeSerie (Me%ObjLandUseTimeSerie, Data1D = dummy)

		enddo

		deallocate(dummy)
         
    end subroutine LandUseTimeSeries
    !--------------------------------------------------------------------------

    subroutine WriteTimeSerieAux(Values,ObjSubBasinTimeSerie,SubBasin)
        !External--------------------------------------------------------------
        real,dimension(:,:), pointer            :: dummy
        real,dimension(:  ), pointer            :: Values
        integer                                 :: ObjSubBasinTimeSerie
        integer                                 :: SubBasin

        !Local-----------------------------------------------------------------
        integer                                 :: i
        integer                                 :: STAT_CALL

        !----------------------------------------------------------------------
        
        dummy => Me%dummy

        do i = 1, SubBasin
            dummy(1,i) = 0.0
            dummy(1,i) = Values       (i)
        enddo

        call WriteTimeSerie(ObjSubBasinTimeSerie,                    &
                            Data2D = dummy,                     &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                              &
            stop 'WriteTimeSerieAux - ModuleSwat - ERR01'


    end subroutine WriteTimeSerieAux
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillSwat(ObjSwatID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjSwatID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           
        integer                             :: STAT_CALL
        integer                             :: i
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSwatID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            if (Me%SwatOptions%SwatWithChanges) then

                if (Me%SwatOptions%MeteoOutputOn) then
                    call KillTimeSerie(Me%ObjMeteoTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSwat - ModuleSwat - ERR00'
                endif

                if (Me%SwatOptions%PlantOutputOn) then
                    call KillTimeSerie(Me%ObjPlantTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSwat - ModuleSwat - ERR01'
                endif

                if (Me%SwatOptions%OperationalOutputOn) then
                    call KillTimeSerie(Me%ObjOperationalTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSwat - ModuleSwat - ERR01'
                endif

                if (Me%SwatOptions%UsleOutputOn) then
                    open(unit=1,file='..\res\AverageAnnualMusle_Sed_ton_ha.dat')
                    do i=1,Me%MaxNumberSubBasin
                        write(1,*) i, Me%SubBasin_Annual_Average%HRU%Usle%Musle_Sed (1,i)
                    end do
                    close(1)
                    call KillTimeSerie(Me%ObjUsleTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSwat - ModuleSwat - ERR01'
                endif

                if (Me%SwatOptions%NitrateOutputOn) then
                    open(unit=1,file='..\res\AverageAnnualNitrateLoad_kgN_ha_year.dat')
                    write (1,*) 'Subbasin NitrogenInRunOff NitrogenLateralFlow  NitrogenGrowndWaterFlow'
                    do i=1,Me%MaxNumberSubBasin
                        write(1,*) i,  & !Me%SubBasin_Annual_Average%HRU%NTransportFlux%Nitrate%Percolation(1,i), &
                                      Me%SubBasin_Annual_Average%SubBasin2%NitrateInRunOff             (i)+ &     
                                      Me%SubBasin_Annual_Average%SubBasin2%RunOffOrganicN%Dissolved    (i)+ & 
                                      Me%SubBasin_Annual_Average%SubBasin2%RunOffOrganicN%Particulated (i), &
                                      Me%SubBasin_Annual_Average%SubBasin2%NitrateInLateralFlow        (i), &
                                      Me%SubBasin_Annual_Average%SubBasin2%NitrateInGrowndWaterFlow    (i)
                    end do
                    close(1)
                    call KillTimeSerie(Me%ObjNitrateTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSwat - ModuleSwat - ERR02'
                endif

                if (Me%SwatOptions%PhosphorusOutputOn) then
                    open(unit=1,file='..\res\AverageAnnualPhosphorusLoad_kgN_ha_year.dat')
                    write (1,*) 'Subbasin PhosphorusInRunOff PhosphorusGrowndWaterFlow'
                    do i=1,Me%MaxNumberSubBasin
                        write(1,*) i, Me%SubBasin_Annual_Average%SubBasin2%RunOffSolubleP        (i)+ &     
                                      Me%SubBasin_Annual_Average%SubBasin2%RunOffSedimentActiveP (i)+ & 
                                      Me%SubBasin_Annual_Average%SubBasin2%RunOffSedimentStableP (i), &
                                      Me%SubBasin_Annual_Average%SubBasin2%GrowndWaterSolubleP   (i)
                    end do
                    close(1)
                    call KillTimeSerie(Me%ObjPhosphorusTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSwat - ModuleSwat - ERR02b'
                endif

                if (Me%SwatOptions%SubBasinOn) then
                    open(unit=1,file='..\res\Flow_Pcp_ET.dat')
                    write (1,*) 'Subbasin Flow Precipitation ActualEvapotranspiration'
                    do i=1,Me%MaxNumberSubBasin
                        write(1,*) i, Me%SubBasin_Annual_Average%SubBasin2%WaterFluxes%GlobalToChannel          (i), &     
                                      Me%SubBasin_Annual_Average%SubBasin2%WaterFluxes%Precipitation            (i), & 
                                      Me%SubBasin_Annual_Average%SubBasin2%WaterFluxes%ActualEvapotranspiration (i)
                                      
                    end do
                    close(1)
				endif

                if (Me%SwatOptions%WaterBalanceOn) then
                    call KillTimeSerie(Me%ObjWaterBalanceTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSwat - ModuleSwat - ERR02c'
                endif

                if (Me%SwatOptions%SubBasinOn) then
                    call KillTimeSerie(Me%ObjSubBasinTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSwat - ModuleSwat - ERR03'
					if (Me%OutPut%Yes) then
						call KillHDF5 (Me%ObjHDF5_SubBasin, STAT = STAT_CALL)
						if (STAT_CALL /= SUCCESS_)                               &
							stop 'KillSwat - ModuleSwat - ERR99'
					endif
                endif

                if (Me%SwatOptions%ReachesOn) then
                    call KillTimeSerie(Me%ObjReachesTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSwat - ModuleSwat - ERR04'

					if (Me%OutPut%Yes) then
						call KillHDF5 (Me%ObjHDF5_Reaches, STAT = STAT_CALL)
						if (STAT_CALL /= SUCCESS_)                               &
							stop 'KillSwat - ModuleSwat - ERR99'
					endif
                endif

                if (Me%SwatOptions%GlobalWaterOn) then
                    call KillTimeSerie(Me%ObjGlobalWaterTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSwat - ModuleSwat - ERR05'
                endif

                if (Me%SwatOptions%GlobalNutrientsOn) then
                    call KillTimeSerie(Me%ObjGlobalNutrientTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSwat - ModuleSwat - ERR06'
                endif


                if (Me%SwatOptions%BasinAverage.and.Me%SwatOptions%SubBasinOn) then
                    call KillTimeSerie(Me%ObjBasinAverage, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSwat - ModuleSwat - ERR08'
                endif

                if (Me%SwatOptions%LandUse) then
                    call KillTimeSerie(Me%ObjLandUseTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSwat - ModuleSwat - ERR07'
                endif

            endif

            nUsers = DeassociateInstance(mSWAT_,  Me%InstanceID)

            if (nUsers == 0) then


                !Deallocates Instance
                call DeallocateInstance ()

                ObjSwatID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillSwat
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Swat), pointer          :: AuxObjSwat
        type (T_Swat), pointer          :: PreviousObjSwat

        !Updates pointers
        if (Me%InstanceID == FirstObjSwat%InstanceID) then
            FirstObjSwat => FirstObjSwat%Next
        else
            PreviousObjSwat => FirstObjSwat
            AuxObjSwat      => FirstObjSwat%Next
            do while (AuxObjSwat%InstanceID /= Me%InstanceID)
                PreviousObjSwat => AuxObjSwat
                AuxObjSwat      => AuxObjSwat%Next
            enddo

            !Now update linked list
            PreviousObjSwat%Next => AuxObjSwat%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjSwat_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSwat_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjSwat_ID > 0) then
            call LocateObjSwat (ObjSwat_ID)
            ready_ = VerifyReadLock (mSWAT_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjSwat (ObjSwatID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSwatID

        !Local-----------------------------------------------------------------

        Me => FirstObjSwat
        do while (associated (Me))
            if (Me%InstanceID == ObjSwatID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleSwat - LocateObjSwat - ERR01'

    end subroutine LocateObjSwat

    !--------------------------------------------------------------------------

end module ModuleSwat






