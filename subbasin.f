      subroutine subbasin
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine controls the simulation of the land phase of the 
!!    hydrologic cycle

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name           |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    auto_wstr(:)   |none          |water stress factor which triggers auto
!!                                  |irrigation
!!    bio_e(:)       |(kg/ha)/      |biomass-energy ratio
!!                   |     (MJ/m**2)|The potential (unstressed) growth rate per
!!                                  |unit of intercepted photosynthetically
!!                                  |active radiation.
!!    canev          |mm H2O        |amount of water evaporated from canopy
!!                                  |storage
!!    ep_day         |mm H2O        |actual amount of transpiration that occurs
!!                                  |on day in HRU
!!    es_day         |mm H2O        |actual amount of evaporation (soil et) that
!!                                  |occurs on day in HRU
!!    gw_q(:)        |mm H2O        |groundwater contribution to streamflow from
!!                                  |HRU on current day
!!    hru_ra(:)      |MJ/m^2        |solar radiation for the day in HRU
!!    iida           |julian date   |day being simulated (current julian date)
!!    idplt(:)       |none          |land cover code from crop.dat
!!    igro(:)        |none          |land cover status code
!!                                  |0 no land cover currently growing
!!                                  |1 land cover growing
!!    inum1          |none          |subbasin number
!!    imp_trig(:)    |none          |release/impound action code:
!!                                  |0 begin impounding water
!!                                  |1 release impounded water
!!    irrsc(:)       |none          |irrigation source code:
!!                                  |1 divert water from reach
!!                                  |2 divert water from reservoir
!!                                  |3 divert water from shallow aquifer
!!                                  |4 divert water from deep aquifer
!!                                  |5 divert water from source outside
!!                                  |  watershed
!!    iurban(:)      |none          |urban simulation code:
!!                                  |0  no urban sections in HRU
!!                                  |1  urban sections in HRU, simulate using
!!                                  |   USGS regression equations
!!                                  |2  urban sections in HRU, simulate using
!!                                  |   build up/wash off algorithm
!!    latq(:)        |mm H2O        |total lateral flow in soil profile for the
!!                                  |day in HRU
!!    nafert(:)      |none          |sequence number of auto-fert application
!!                                  |within the year
!!    nair(:)        |none          |sequence number of auto-irrigation 
!!                                  |application within the year
!!    nfert(:)       |none          |sequence number of fertilizer application
!!                                  |within the year
!!    nirr(:)        |none          |sequence number of irrigation application
!!                                  |within the year
!!    nrelease(:)    |none          |sequence number of impound/release
!!                                  |operation within the year
!!    nro(:)         |none          |sequence number of year in rotation
!!    peakr          |m^3/s         |peak runoff rate
!!    pet_day        |mm H2O        |potential evapotranspiration on current
!!                                  |day in HRU
!!    phuacc(:)      |none          |fraction of plant heat units accumulated
!!    phubase(:)     |heat units    |base zero total heat units (used when no
!!                                  |land cover is growing)
!!                                  |pesticide application occurs
!!    pot_fr(:)      |km2/km2       |fraction of HRU area that drains into
!!                                  |pothole
!!    pot_vol(:)     |m**3 H2O      |current volume of water stored in the 
!!                                  |depression/impounded area
!!    precipday      |mm H2O        |precipitation for the day in HRU
!!    qday           |mm H2O        |surface runoff loading to main channel from
!!                                  |HRU for day
!!    qtile          |mm H2O        |drainage tile flow in soil layer for the 
!!                                  |day
!!    sci(:)         |none          |retention coefficient for CN method based
!!                                  |on plant ET
!!    sedyld(:)      |metric tons   |soil loss for day in HRU
!!    smx(:)         |none          |retention coefficient for CN method based
!!                                  |on soil moisture
!!    surfq(:)       |mm H2O        |surface runoff generated on day in HRU
!!    tmn(:)         |deg C         |minimum temperature for the day in HRU
!!    tmpav(:)       |deg C         |average temperature for the day in HRU
!!    tmx(:)         |deg C         |maximum temperature for the day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    albday      |none          |albedo, the fraction of the solar radiation
!!                               |reflected at the soil surface back into
!!                               |space
!!    etday       |mm H2O        |actual evapotranspiration occuring on day
!!                               |in HRU
!!    ihru        |none          |HRU number
!!    inflpcp     |mm H2O        |amount of precipitation that infiltrates
!!                               |into soil (enters soil)
!!    nafert(:)   |none          |sequence number of auto-fert application
!!                               |within the year
!!    nair(:)     |none          |sequence number of auto-irrigation 
!!                               |application within the year
!!    qdfr        |none          |fraction of water yield that is surface
!!                               |runoff
!!    qdr(:)      |mm H2O        |total amount of water entering main channel
!!                               |for day from HRU
!!    sci(:)      |none          |retention coefficient for CN method based
!!                               |on plant ET
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    d           |
!!    gma         |kPa/deg C     |psychrometric constant
!!    ho          |              |net radiation
!!    j           |none          |HRU number
!!    pet_alpha   |none          |alpha factor in Priestley-Taylor ET 
!!                               |equation
!!    tmpk        |deg K         |average temperature for the day in the HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Max
!!    SWAT: varinit, albedo, solt, surface, percmain, etpot, etact, fert
!!    SWAT: confert, graze, plantmod, nminrl, nitvol, pminrl, gwmod, apply, gwmod_deep
!!    SWAT: washp, decay, pestlch, enrsb, pesty, orgn, psed, nrain, nlch
!!    SWAT: solp, subwq, bacteria, urban, pothole, latsed, surfstor
!!    SWAT: substor, wetland, hrupond, irrsub, autoirr, watuse, watbal
!!    SWAT: sumv, virtual

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      integer :: j,sb,kk
      real :: tmpk, d, gma, ho, pet_alpha, aphu, phuop

	real :: VolumeToChannel, RunOffVolumeToChannel  !PCL
	real :: LayerDepth                              !PCL
	integer :: idum_, HRU_Properties       			!PCL

      ihru = 0
      ihru = hru1(inum1) 
      
      call sub_subbasin

      do iihru = 1, hrutot(inum1)

      j = 0
      j = ihru


      !!by zhang DSSAT tillage
      !!======================
      !!    deptil(:)   |mm  |depth of mixing caused by tillage operation
      !jj is hru number
      if (cswat == 2) then
          if (tillage_switch(ihru) .eq. 1) then
              if (tillage_days(ihru) .ge. 30) then
                    tillage_switch(ihru) = 0
                    tillage_days(ihru) = 0
              else
                    tillage_days(ihru) = tillage_days(ihru) + 1
              end if                
              !tillage_depth(ihru) = dtil
              !tillage_switch(ihru) = .TRUE. 
          end if
      end if
      !!by zhang DSSAT tillage  
      !!====================== 



      call varinit
      if (icr(j) <= 0) icr(j) = 1
      
      i_wtrhru = 0
      idplrot(icr(j),ihru) = idplt(j)
      if (idplt(j) /= 0) then
          if (cpnm(idplt(j)) == "WATR") then
              i_wtrhru = 1
          end if
      endif
      
	if (i_wtrhru == 1) then
         call water_hru
      else 

        !! Simulate land covers other than water

        !! update base zero total heat units
        if (tmpav(j) > 0. .and. phutot(hru_sub(j)) > 0.01) then
           phubase(j) = phubase(j) + tmpav(j) / phutot(hru_sub(j))
        end if
        
        call schedule_ops

        !! calculate albedo for day
        call albedo

        !! calculate soil temperature for soil layers
        call solt

!       if (ipot(j) /= j .and. imp_trig(nro(j),nrelease(j),j)==1)       &  Srini pothole
!
!     &        then             
          !! calculate surface runoff if HRU is not impounded or an 
          !! undrained depression--
          call surface

          !! add surface flow that was routed across the landscape on the previous day
       !!   qday = qday + surfq_ru(j)
       !!   surfq_ru(j) = 0.
          
          !! compute effective rainfall (amount that percs into soil)
          inflpcp = Max(0.,precipday - surfq(j))
!        end if
         
        !! perform management operations
        if (yr_skip(j) == 0) call operatn
          
        if (auto_wstr(j) > 1.e-6 .and. irrsc(j) > 2) call autoirr       !!NUBZ
        
        !! perform soil water routing
        call percmain

        !! compute evapotranspiration
        call etpot
!        if (pot_vol(j) < 1.e-6) call etact
        call etact

        !! compute water table depth using climate drivers
        call wattable

        !! new CN method
        if (icn == 1) then 
        sci(j) = sci(j) + pet_day*exp(-cncoef_sub(hru_sub(j))*sci(j)/   &
     &    smx(j)) - precipday + qday + qtile + latq(j) + sepbtm(j)
        else if (icn == 2) then 
        sci(j) = sci(j) + pet_day*exp(-cncoef_sub(hru_sub(j))*sci(j)/   &
     &    smx(j)) - precipday + qday + latq(j) + sepbtm(j) + qtile
        sci(j) = amin1(sci(j),smxco * smx(j))
        end if 
        
        !! apply fertilizer/manure in continuous fert operation
        if (icfrt(j) == 1) then
          ndcfrt(j) = ndcfrt(j) + 1
          call confert
        end if
        
        !! apply pesticide in continuous pest operation
        if (icpst(j) == 1) then 
          ndcpst(j) = ndcpst(j) + 1
          call conapply
        end if 
        
        !! remove biomass from grazing and apply manure
        if (igrz(j) == 1) then
          ndeat(j) = ndeat(j) + 1
          call graze
        end if
       
        !! compute crop growth
        call plantmod
        
        !! check for dormancy
        if (igro(j) == 1) call dormant
        !! compute actual ET for day in HRU
        etday = ep_day + es_day + canev

        !! write daily air and soil temperature file
        !! can be uncommmented if needed by user and also in readfile.f

!      write (120,12112) i,j,tmx(j),tmn(j),(sol_tmp(k,j),k=1,sol_nly(j))
!12112  format (2i4,12f8.2)

        !! compute nitrogen and phosphorus mineralization 

      if (cswat == 0) then
        call nminrl
	end if
	if (cswat == 1) then
		call carbon
	end if
	
	!! Add by zhang
	!!=================
	if (cswat == 2) then
	  call carbon_zhang2
	end if
	!! Add by zhang
	!!=================	

        call nitvol
        if (sol_P_model == 1) then
            call pminrl
        else
            call pminrl2
        end if

!!    compute biozone processes in septic HRUs
!!    if 1)current is septic hru and 2)  soil temperature is above zero
	  if (isep_opt(j)/=0.and.iyr>=isep_iyr(j)) then
	   if (sol_tmp(i_sep(j),j) > 0.) call biozone     
	  endif

        !! compute ground water contribution
        call gwmod
        call gwmod_deep

        !! compute pesticide washoff   
        if (precipday >= 2.54) call washp

        !! compute pesticide degradation
        call decay

        !! compute pesticide movement in soil
        call pestlch

        if (surfq(j) > 0. .and. peakr > 1.e-6) then
          if (precipday > 0.) then
            call enrsb(0)
            if (sedyld(j) > 0.) call pesty(0)

		  if (cswat == 0) then
			call orgn(0)
	    end if
	    if (cswat == 1) then
	    
		    call orgncswat(0)
		  end if
		  
		  !! Add by zhang
		  !! ====================
		  if (cswat == 2) then
		    call orgncswat2(0)
		  end if
		  !! Add by zhang
		  !! ====================

            call psed(0)
          end if
        end if

        !! add nitrate in rainfall to soil profile
        call nrain

        !! compute nitrate movement leaching
        call nlch

        !! compute phosphorus movement
        call solp

        !! compute chl-a, CBOD and dissolved oxygen loadings
        call subwq

        !! compute bacteria transport
        call bacteria

        !! compute loadings from urban areas
        if (urblu(j) > 0) then
	     if(ievent<3) then
	        call urban ! daily simulation
	     else
		     call urbanhr ! subdaily simulation J.Jeong 4/20/2009
	     endif
	  endif
	  
!! Srini Pothole
        !! compute undrained depression/impounded area (eg rice) processes
!        if (pot_fr(j) > 0.) then
!           if (ievent<3) then   
!          call pothole
!           else
!              call potholehr
!           endif
!        endif
        
        !! compute sediment loading in lateral flow and add to sedyld
        call latsed

        !! compute nutrient loading in groundwater flow
        call gwnutr
        call gw_no3

        !! lag nutrients and sediment in surface runoff
        call surfstor

        !! lag subsurface flow and nitrate in subsurface flow

        call substor

        !! add lateral flow that was routed across the landscape on the previous day
      !!  latq(j) = latq(j) + latq_ru(j)
      !!  latq_ru(j) = 0.
        
        !! compute reduction in pollutants due to edge-of-field filter strip
        if (vfsi(j) >0.)then
          call filter
          if (filterw(j) > 0.) call buffer
        end if
              if (vfsi(j) == 0. .and. filterw(j) > 0.) then 
                call filtw
                call buffer
              end if

	 !! compute reduction in pollutants due to in field grass waterway
         if (grwat_i(j) == 1) then
          call grass_wway
        end if

	 !! compute reduction in pollutants due to in fixed BMP eff
	   if (bmp_flag(j) == 1) then
          call bmpfixed
        end if


        !! compute water yield for HRU
        qdr(j) = qday + latq(j) + gw_q(j) + qtile + gw_qdeep(j)
        if (qdr(j) < 0.) qdr(j) = 0.
        if (qdr(j) > 0.) then
          qdfr = qday / qdr(j)
        else
          qdfr = 0.
        end if

        !! compute wetland processes
        call wetlan

        !! compute pond processes
        if (ievent<3) then
           call hrupond
        else
           call hrupondhr
        endif
        
!       Srini pothole        
        if (pot_fr(j) > 0.) call pothole
                
        xx = sed_con(j)+soln_con(j)+solp_con(j)+orgn_con(j)+orgp_con(j)
        if (xx > 1.e-6) then
          call urb_bmp
        end if

        if (SwatWithChanges) then
	    OutputHruMeteo  (11,iihru,inum1) = aird(j)
        endif

        !! consumptive water use (ponds, shallow aquifer, deep aquifer)
        call watuse

        !! perform water balance
        call watbal
        
        !! qdayout is surface runoff leaving the hru - after wetlands, ponds, and potholes
        qdayout(j) = qday

      endif

      !! perform output summarization
      call sumv

      !! summarize output for multiple HRUs per subbasin
      !! store reach loadings for new fig method
      call virtual
      aird(j) = 0.
      
            if (SwatWithChanges) then
              OutputHruMeteo  (1 ,iihru,inum1) = tmpav (j)
	        OutputHruMeteo  (2 ,iihru,inum1) = subp (j)
	        OutputHruMeteo  (3 ,iihru,inum1) = pet_day
	        OutputHruMeteo  (4 ,iihru,inum1) = ep_max
	        OutputHruMeteo  (5 ,iihru,inum1) = etday
	        OutputHruMeteo  (6 ,iihru,inum1) = tmn(j)
	        OutputHruMeteo  (7 ,iihru,inum1) = tmx(j)
	        OutputHruMeteo  (8 ,iihru,inum1) = u10(j)
	        OutputHruMeteo  (9 ,iihru,inum1) = rhd(j)
	        OutputHruMeteo  (10,iihru,inum1) = hru_ra(j)
	        !OutputHruMeteo  (11,iihru,inum1) = aird(j) !this must be done before subroutine Virtual
	        OutputHruMeteo  (12,iihru,inum1) = snofall
	        OutputHruMeteo  (13,iihru,inum1) = snomlt
	        OutputHruMeteo  (14,iihru,inum1) = sno_hru(j)
	        OutputHruMeteo  (15,iihru,inum1) = sol_tmp(2,j)

              LayerDepth = sol_z(1,j)
              do layer = 1, sol_nly(j)
	            OutputHruMeteo  (16,iihru,inum1) = 
     &            OutputHruMeteo  (16,iihru,inum1) + sol_st(layer,j) +
     &            sol_wp (layer,j) * LayerDepth
	            LayerDepth = sol_z(layer+1,j)-sol_z(layer,j)
	        enddo
              OutputHruMeteo  (17,iihru,inum1) = shallst(j)
              OutputHruQuality(1 ,iihru,inum1) = bio_ms (j)
              OutputHruQuality(2 ,iihru,inum1) = pltfr_n(j)
              OutputHruQuality(3 ,iihru,inum1) = 
     &		cnyld(idplt(j))
              OutputHruQuality(4 ,iihru,inum1) = phuacc(j)

              OutputHruQuality(5 ,iihru,inum1) = hmntl
              OutputHruQuality(6 ,iihru,inum1) = rmn2tl
              OutputHruQuality(7 ,iihru,inum1) = wdntl

              OutputHruQuality(8 ,iihru,inum1) = latno3(j)
              OutputHruQuality(9 ,iihru,inum1) = no3pcp
              OutputHruQuality(10,iihru,inum1) = nplnt(j)
              OutputHruQuality(11,iihru,inum1) = percn(j)
              OutputHruQuality(12,iihru,inum1) = surqno3(J)

              OutputHruQuality(13,iihru,inum1) = auton
              OutputHruQuality(14,iihru,inum1) = fertn
      
              do layer = 1, sol_nly(j)
                  OutputHruQuality(15,iihru,inum1) = 
     &            OutputHruQuality(15,iihru,inum1) + 
     &                                   sol_no3(layer,j)
	        enddo

              do layer = 1, sol_nly(j)
                  OutputHruQuality(16,iihru,inum1) = 
     &            OutputHruQuality(16,iihru,inum1) + 
     &                                   sol_fon(layer,j)
	        enddo

              OutputHruQuality(17,iihru,inum1) = grazn
              OutputHruQuality(18,iihru,inum1) = Get_rvol

              do layer = 1, sol_nly(j)
                  OutputHruQuality(19,iihru,inum1) = 
     &            OutputHruQuality(19,iihru,inum1) + 
     &                                   sol_aorgn(layer,j)
	        enddo

              do layer = 1, sol_nly(j)
                  OutputHruQuality(20,iihru,inum1) = 
     &            OutputHruQuality(20,iihru,inum1) + 
     &                                   sol_orgn(layer,j)
	        enddo

              OutputHruQuality(21,iihru,inum1) = laiday(j)

              !USLE
              OutputHruQuality(22,iihru,inum1) = sedyld(j) / 
     &                                                (hru_km(j) * 100)
              OutputHruQuality(23,iihru,inum1) = usle

              OutputHruQuality(24,iihru,inum1) = strsp(j)

              OutputHruQuality(25,iihru,inum1) = strsn(j)

              OutputHruQuality(26,iihru,inum1) = strstmp(j)

              OutputHruQuality(27,iihru,inum1) = strsw (j)

	        OutputHruQuality(28,iihru,inum1) = cklsp(j)

	        OutputHruQuality(29,iihru,inum1) = peakr

			OutputHruQuality(30,iihru,inum1) = usle_ei

	        OutputHruQuality(31,iihru,inum1) = usle_mult(j)


			if (gw_q(j)>0) then
			![kg]*[mg1e6/kg]/[ha]*[1e2*há/km2]*[km2])/([mm]*[m/mm1e+3]*[km2]*[1e6m/km2]*[dm3*1e3/m3])
			! cortando 1e6 do numerador com o 1e6 do denominador
			OutputHruQuality(32,iihru,inum1) = no3gw(j)*1e2*hru_km(j) /
     &                                          (gw_q(j) *
     &                                          1e-3*hru_km(j)*1e3)
			else
			OutputHruQuality(32,iihru,inum1) = 0.0
			endif

			OutputHruQuality(33,iihru,inum1) = sol_rsd(1,j)

              do layer = 1, sol_nly(j)
                  OutputHruQuality(34,iihru,inum1) = 
     &            OutputHruQuality(34,iihru,inum1) + 
     &                                   sol_nh3(layer,j)
	        enddo

              OutputHruQuality(35,iihru,inum1) = Get_rnit

              do layer = 1, sol_nly(j)
                  OutputHruQuality(36,iihru,inum1) = 
     &            OutputHruQuality(36,iihru,inum1) + 
     &                                   sol_actp(layer,j)
	        enddo
              do layer = 1, sol_nly(j)
                  OutputHruQuality(37,iihru,inum1) = 
     &            OutputHruQuality(37,iihru,inum1) + 
     &                                   sol_fop(layer,j)
	        enddo
              do layer = 1, sol_nly(j)
                  OutputHruQuality(38,iihru,inum1) = 
     &            OutputHruQuality(38,iihru,inum1) + 
     &                                   sol_orgp(layer,j)
	        enddo
              do layer = 1, sol_nly(j)
                  OutputHruQuality(39,iihru,inum1) = 
     &            OutputHruQuality(39,iihru,inum1) + 
     &                                   sol_stap(layer,j)
	        enddo

              OutputHruQuality(40,iihru,inum1) = autop 
              OutputHruQuality(41,iihru,inum1) = fertp 
              OutputHruQuality(42,iihru,inum1) = grazp 
              OutputHruQuality(43,iihru,inum1) = sedminpa(j)  
              OutputHruQuality(44,iihru,inum1) = sedminps(j)  
              OutputHruQuality(45,iihru,inum1) = sedorgp(j) 
              OutputHruQuality(46,iihru,inum1) = pplnt(j)  

			OutputHruQuality(47,iihru,inum1) = sedorgn (j)
	        
			OutputHruParameters (1,iihru,inum1) = hru_dafr (j)

			  idum_ = idplt(j)
			  if (idum_ > 0) then
				OutputHRU_cropname(iihru,inum1) = cpnm(idum_)
			  else
				OutputHRU_cropname(iihru,inum1) = "NOCR"
			  endif

			do HRU_Properties = 1, 68
			
			 HRU_output_values (HRU_Properties,iihru,inum1) = 
     &         HRU_output_values0(HRU_Properties)

			enddo
			      
			! pools of water (dstor)
			OutputHruWatBal  (1 ,iihru,inum1) = sno_hru(j) 
			OutputHruWatBal  (2 ,iihru,inum1) = snoprev 
			OutputHruWatBal  (3 ,iihru,inum1) = sol_sw(j) 
			OutputHruWatBal  (4 ,iihru,inum1) = swprev 
              OutputHruWatBal  (5 ,iihru,inum1) = shallst(j) 
              OutputHruWatBal  (6 ,iihru,inum1) = shallstp 
			OutputHruWatBal  (7 ,iihru,inum1) = deepst(j) 
			OutputHruWatBal  (8 ,iihru,inum1) = deepstp 
              OutputHruWatBal  (9 ,iihru,inum1) = surf_bs(1,j) 
              OutputHruWatBal  (10,iihru,inum1) = bsprev 
			OutputHruWatBal  (11,iihru,inum1) = bss(1,j) 
			OutputHruWatBal  (12,iihru,inum1) = bssprev

			! water fluxes (h2oloss)
			OutputHruWatBal  (13,iihru,inum1) = subp(j) 
			OutputHruWatBal  (14,iihru,inum1) = snoev 
			OutputHruWatBal  (15,iihru,inum1) = qday 
			OutputHruWatBal  (16,iihru,inum1) = latq(j) 
			OutputHruWatBal  (17,iihru,inum1) = etday 
			OutputHruWatBal  (18,iihru,inum1) = gw_q(j) 
			OutputHruWatBal  (19,iihru,inum1) = revapday 
			OutputHruWatBal  (20,iihru,inum1) = twlpnd 
			OutputHruWatBal  (21,iihru,inum1) = twlwet 
			OutputHruWatBal  (22,iihru,inum1) = aird(j) 
			OutputHruWatBal  (23,iihru,inum1) = rchrg(j) 
			OutputHruWatBal  (24,iihru,inum1) = qtile 
			OutputHruWatBal  (25,iihru,inum1) = sepbtm(j)

	endif
      
      ihru = ihru + 1
      end do

      !! route 2 landscape units
      if (ils2flag(inum1) > 0) then
      isub = inum1                        ! save the subbasin number
      
      !! calculate outputs from hillslope
      ihout1 = mhyd_bsn + (inum1 - 1) * 4 ! first outflow hyd number
      ihout = ihout1                      ! outflow hyd number
      inum1 = 1                           ! landscape unit number
      inum2 = isub                        ! subbasin number
      call routeunit                      ! hillslope unit
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      ihouts(ihout) = ihout
      
      !! calculate outputs from valley bottom
      inum1 = 2                           ! landscape unit number
      ihout = ihout + 1                   ! outflow hyd number
      sumdaru = 0.
      do j = 1, hrutot(isub)
        sumdaru = sumdaru + hru_km(j)
      end do 
      daru_km(inum2,inum1) = sumdaru
      call routeunit                      ! valley bottom unit
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      ihouts(ihout) = ihout
      
      !! route output from hillslope across valley bottom
      ihout = ihout + 1                   ! outflow hyd number
      inum1 = 2                           ! valley bottom landscape unit
      inum2 = ihout1                      ! inflow hyd=outlfow from hillslope
      inum3 = isub                        ! subbasin number
      rnum1 = 1.                          ! fraction overland flow
      iru_sub = 1                         ! route across landscape unit
      !! compute weighted K factor for sediment transport capacity
      sumk = 0.
      ovsl = 0.
      ovs = 0.
      do j = 1, hrutot(isub)
        sumk = sumk + usle_k(j) * hru_rufr(inum1,j)
        ovsl = ovsl + slsubbsn(j)
        ovs = ovs + hru_slp(j)
      end do 
      ovsl = ovsl / hrutot(isub)
      ovs = ovs / hrutot(isub)
      ru_k(isub,inum1) = sumk
      ru_ovsl(isub,inum1) = ovsl
      ru_ovs(isub,inum1) = ovs
      ru_ktc(isub,inum1) = 50.
      ru_a(isub,inum1) = daru_km(isub,1) / ru_ovsl(isub,inum1)
      call routels(iru_sub)               ! route across valley bottom
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      inum3s(ihout) = inum3
      ihouts(ihout) = ihout
      
      !! add routed with valley bottom loading
      inum1 = ihout                       ! hyd from routed 
      inum2 = ihout - 1                   ! hyd from loading
      ihout = ihout + 1                   ! outflow hyd number
      call addh                           ! add hyd's
      call sumhyd
      inum1s(ihout) = inum1
      inum2s(ihout) = inum2
      ihouts(ihout) = ihout
      
      !! save landscape routed output in place of subbasin output for routing
      varoute(isub,:) = varoute(ihout,:)
      end if
      
 1000 format(4i10,a10)
      if (SwatWithChanges) then
	  !Water from each subbasin that enters channel
	  ! [m3]   = [mm] * [m/1000mm] * [ha] * [-] * [10000m2/ha]
	  VolumeToChannel = sub_wyld(inum1) * da_ha * sub_fr(inum1) * 10. 
        ! [m3/s] = [m3] / [s]
        OutputSubBasin(1,inum1) =   
     &        VolumeToChannel / 86400.
        !GrowndWater from each subbasin that enters channel
	  ! [m3/s] = [mm] * [m/1000mm] * [ha] * [-] * [10000m2/ha] / [s]
        OutputSubBasin(2,inum1) =   
     &        sub_gwq(inum1) * da_ha * sub_fr(inum1) * 10. / 86400.
        !RunOff Water from each subbasin that enters channel
	  ! [m3]   = [mm] * [m/1000mm] * [ha] * [-] * [10000m2/ha]
        RunOffVolumeToChannel = sub_qd(inum1) * da_ha * sub_fr(inum1) * 
     &  10.
        ! [m3/s] = [m3] / [s]
        OutputSubBasin(3,inum1) =   
     &         RunOffVolumeToChannel / 86400.
        !Lateral Flow from each subbasin that enters channel
        OutputSubBasin(4,inum1) =   
     &        sub_latq(inum1) * da_ha * sub_fr(inum1) * 10. / 86400.
        if(RunOffVolumeToChannel>0) then
            !Nitrate concentration from RunOff
            ![mgN/dm3] = [1e6mgN/kgN] * [m3/1e3dm3] * [kgN/ha] * [ha] * [-] / [m3] 
            !ALL SOURCES OF NO3
            OutputSubBasin(6,inum1) =   
     &      1e+3 * 
     &      (sub_no3(inum1) + sub_latno3(inum1) + sub_gwno3(inum1)) 
     &      * da_ha * sub_fr(inum1) 
     &      / VolumeToChannel
	  else
            OutputSubBasin(6,inum1) = 0.0
	  endif
        if(RunOffVolumeToChannel>0) then
            !Nitrogen concentration from RunOff
            ![mgN/dm3] = [1e6mgN/kgN] * [m3/1e3dm3] * [kgN/ha] * [ha] * [-] / [m3] 
            OutputSubBasin(7,inum1) =   
     &      1e+3 * sub_yorgnDON(inum1) * da_ha * sub_fr(inum1) / 
     &      VolumeToChannel
	  else
            OutputSubBasin(7,inum1) = 0.0
	  endif
        if(RunOffVolumeToChannel>0) then
            !Nitrogen concentration from RunOff
            ![mgN/dm3] = [1e6mgN/kgN] * [m3/1e3dm3] * [ha] * [-] / [m3] 
            OutputSubBasin(8,inum1) =   
     &      1e+3 * sub_yorgnPON(inum1) * da_ha * sub_fr(inum1) / 
     &      VolumeToChannel
	  else
            OutputSubBasin(8,inum1) = 0.0
	  endif

        if(RunOffVolumeToChannel>0) then
            !Phosphorous organic concentration from RunOff
            ![mgP/dm3] = [1e6mgP/kgP] * [m3/1e3dm3] * ([kgP/ha] + [kgP/ha]) * [ha] * [-] / [m3] 
            OutputSubBasin(9,inum1) =   
     &      1e+3 * (sub_solp(inum1) + sub_sedps(inum1)) * da_ha *  
     &      sub_fr(inum1) / VolumeToChannel
            !Phosphorous mineral concentration from RunOff
            ![mgP/dm3] = [1e6mgP/kgP] * [m3/1e3dm3] * [kgP/ha] * [ha] * [-] / [m3] 
            !ALL SOURCES OF soluble P
            OutputSubBasin(10,inum1) =   
     &      1e+3 * (sub_sedpa(inum1) + sub_gwsolp(inum1)) 
     &      * da_ha * sub_fr(inum1) /
     &      VolumeToChannel
	  else
            OutputSubBasin(9, inum1) = 0.0
            OutputSubBasin(10,inum1) = 0.0
	  endif

        if(RunOffVolumeToChannel>0) then
            !Sediment concentration from RunOff
            ![mgSed/dm3] = [1e9mgSed/tonSed] * [m3/1e3dm3] * [tonSed]  / [m3] 
            OutputSubBasin(11,inum1) =   
     &      1e+6 * sedyld(j) /
     &      VolumeToChannel
	  else
            OutputSubBasin(11,inum1) = 0.0
	  endif

	  ![tonSed/ha]             = [tonSed]  / [ha]
        OutputSubBasin(12,inum1) = sedyld(j) / da_ha

        OutputSubBasin(13,inum1) = sub_fr(inum1)

            OutputSubBasin2(1,inum1) = sub_wyld(inum1)
            OutputSubBasin2(2,inum1) = sub_gwq(inum1) 
            OutputSubBasin2(3,inum1) = sub_qd(inum1)
            OutputSubBasin2(4,inum1) = sub_latq(inum1) 
            OutputSubBasin2(6,inum1) = sub_no3(inum1) 
            OutputSubBasin2(7,inum1) = sub_yorgnDON(inum1) 
            OutputSubBasin2(8,inum1) = sub_yorgnPON(inum1) 
            OutputSubBasin2(9,inum1) = sub_solp(inum1)   
            OutputSubBasin2(10,inum1) = sub_sedpa(inum1)
            OutputSubBasin2(11,inum1) = sedyld(j) 
            OutputSubBasin2(13,inum1) = sub_fr(inum1)
	      OutputSubBasin2(14,inum1) = sub_latno3(inum1)
	      OutputSubBasin2(15,inum1) = sub_gwno3(inum1)
	      OutputSubBasin2(16,inum1) = sub_sedps(inum1)
	      OutputSubBasin2(17,inum1) = sub_gwsolp(inum1)
	      OutputSubBasin2(18,inum1) = sub_subp(inum1)
	      OutputSubBasin2(19,inum1) = sub_etday(inum1)

	endif
      return
      end
