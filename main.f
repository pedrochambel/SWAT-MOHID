      !PCL include 'modparm.f'
      program main
!!    this is the main program that reads input, calls the main simulation
!!    model, and writes output.
!!    comment changes to test merging with trunk and c:\branch_test code
!!    two lines added to c:\branch_test code

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!         ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    date        |NA            |date simulation is performed where leftmost
!!                               |eight characters are set to a value of
!!                               |yyyymmdd, where yyyy is the year, mm is the 
!!                               |month and dd is the day
!!    isproj      |none          |special project code:
!!                               |1 test rewind (run simulation twice)
!!    time        |NA            |time simulation is performed where leftmost
!!                               |ten characters are set to a value of
!!                               |hhmmss.sss, where hh is the hour, mm is the 
!!                               |minutes and ss.sss is the seconds and
!!                               |milliseconds
!!    values(1)   |year          |year simulation is performed
!!    values(2)   |month         |month simulation is performed
!!    values(3)   |day           |day in month simulation is performed
!!    values(4)   |minutes       |time difference with respect to Coordinated
!!                               |Universal Time (ie Greenwich Mean Time)
!!    values(5)   |hour          |hour simulation is performed
!!    values(6)   |minutes       |minute simulation is performed
!!    values(7)   |seconds       |second simulation is performed
!!    values(8)   |milliseconds  |millisecond simulation is performed
!!    zone        |NA            |time difference with respect to Coordinated
!!                               |Universal Time (ie Greenwich Mean Time)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    prog        |NA            |program name and version
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    i           |none          |counter
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 


!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: date_and_time
!!    SWAT: getallo, allocate_parms, readfile, readfig
!!    SWAT: readbsn, std1, readwwq, readinpt, std2, storeinitial
!!    SWAT: openwth, headout, simulate, finalbal, writeaa, pestw 

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm
	use ModuleSwat !SwatWithChanges     !PCL 
      implicit none

      integer :: k, n !PCL
	
      MainObjSwatID = 0 !PCL      
      
      prog = "SWAT Apr 12 2013    VER 2012/Rev 591"

      write (*,1000)
 1000 format(1x,"               SWAT2012               ",/,             &
     &          "               Rev. 591               ",/,             &
     &          "      Soil & Water Assessment Tool    ",/,             &
     &          "               PC Version             ",/,             &
     &          " Program reading from file.cio . . . executing",/)

!! process input
		
      call getallo
      call allocate_parms
      call readfile
      call readbsn
      call readwwq
      if (fcstyr > 0 .and. fcstday > 0) call readfcst
      call readplant             !! read in the landuse/landcover database
      call readtill              !! read in the tillage database
      call readpest              !! read in the pesticide database
      call readfert              !! read in the fertilizer/nutrient database
      call readurban             !! read in the urban land types database
      call readseptwq            !! read in the septic types database     
      call readlup
      call readfig
      call readatmodep
      call readinpt
      call std1
      call std2
      call openwth
      call headout

      MaxHruInSubBasin = 0  !PCL
	do i=1,msub
	    if(MaxHruInSubBasin < hrutot(i)) then
              MaxHruInSubBasin = hrutot(i)
	    endif
	enddo

      call StartSwat (MainObjSwatID,  
     &                       iyr,     
     &                       nbyr,    
     &                       idaf,    
     &                       idal,
     &                       MaxHruInSubBasin,
     &                       msub,
     &                       da_km,
     &                       SwatWithChanges,
     &                       DenitrificationThreshold,
     &                       LateralFlowTravelTime) !PCL

	if (SwatWithChanges) then
		allocate(OutputHRU_Cropname     (MaxHruInSubBasin,msub)   )           
		allocate(HRU_output_values0     (mhruo)				      )           
		allocate(HRU_output_values      (mhruo,MaxHruInSubBasin,msub))				                 
		allocate(OutputHruMeteo         (17,MaxHruInSubBasin,msub))
		allocate(OutputHruWatBal        (25,MaxHruInSubBasin,msub))
          allocate(OutputHruQuality       (47,MaxHruInSubBasin,msub))
          allocate(OutputHruParameters    ( 1,MaxHruInSubBasin,msub))
	    allocate(OutputSubBasin         (13,msub)                 )
	    allocate(OutputSubBasin2        (19,msub)                 )
	    allocate(OutputReaches          (14,msub)                 )
	    allocate(OutputGlobalBasin      (mstdo )                  )
          if (LateralFlowTravelTime>0) then
              do i=1,mhru
                  lat_ttime (i) = LateralFlowTravelTime
	        enddo
	    endif

	endif
	
      !! convert integer to string for output.mgt file
      subnum = ""
      hruno = ""
      do i = 1, mhru
        write (subnum(i),fmt=' (i5.5)') hru_sub(i)
        write (hruno(i),fmt=' (i4.4)') hru_seq(i)  
      end do

      if (isproj == 2) then 
        hi_targ = 0.0
      end if

!! save initial values
      if (isproj == 1) then
        scenario = 2
        call storeinitial
      else if (fcstcycles > 1) then
        scenario =  fcstcycles
        call storeinitial
      else
        scenario = 1
      endif
        if (iclb /= 4) then
      do iscen = 1, scenario

     
        !! simulate watershed processes
        call simulate

        if (SwatWithChanges) call KillSwat (MainObjSwatID)
        
        !! perform summary calculations
        call finalbal
        call writeaa
        call pestw

        !!reinitialize for new scenario
        if (scenario > iscen) call rewind_init
      end do
         end if
      do i = 101, 109       !Claire 12/2/09: change 1, 9  to 101, 109.
        close (i)
      end do
      close(124)
      write (*,1001)
 1001 format (/," Execution successfully completed ")
	

        iscen=1
	stop
      end