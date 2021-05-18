subroutine sample(&
	IA,& 
	IS,& 
	IE,&
	JA,& 
	JS,& 
	JE,&
	TA,& 
	TS,& 
	TE,&
	DX,&
	DY,&
	ATM_TEMP,&
	SFC_TEMP,&
	ATM_PRES,&
	SFC_PRES,&
	ATM_QV  ,&
	ATM_U   ,&
	ATM_V   ,&
	ATM_Z1  ,&
	PBL     ,&
	SFC_Z0M ,&
	SFC_Z0H ,&
	SFC_Z0E ,&
	SFC_DENS,&
	SFLX_MW ,&
	SFLX_MU ,&
	SFLX_MV )
	!bind(C)

  use scale
  use scale_precision
  use scale_const, only: &
    EPSvap => CONST_EPSvap
  use scale_bulkflux, only: &
    BULKFLUX_setup, &
    BULKFLUX
  use scale_atmos_saturation, only: &
    SATURATION_psat_all => ATMOS_SATURATION_psat_all

  implicit none

  integer,  intent(in) :: IA, IS, IE
  integer,  intent(in) :: JA, JS, JE
  integer,  intent(in) :: TA, TS, TE
  real(RP), intent(in) :: DX
  real(RP), intent(in) :: DY
  real(RP), intent(in) :: ATM_TEMP(IA,JA,TA)   ! temperature at the lowermost layer (cell center) [K]
  real(RP), intent(in) :: SFC_TEMP(IA,JA,TA)   ! temperature at the surface skin [K]
  real(RP), intent(in) :: ATM_PRES(IA,JA,TA)   ! pressure    at the lowermost layer (cell center) [Pa]
  real(RP), intent(in) :: SFC_PRES(IA,JA,TA)   ! pressure    at the surface atmosphere [Pa]
  real(RP), intent(in) :: ATM_QV  (IA,JA,TA)   ! qv          at the lowermost layer (cell center) [kg/kg]
  real(RP), intent(in) :: ATM_U   (IA,JA,TA)   ! velocity u  at the lowermost layer (cell center) [m/s]
  real(RP), intent(in) :: ATM_V   (IA,JA,TA)   ! velocity v  at the lowermost layer (cell center) [m/s]
  real(RP), intent(in) :: ATM_Z1  (IA,JA)      ! height of the lowermost grid from surface (cell center) [m]
  real(RP), intent(in) :: PBL     (IA,JA,TA)   ! depth of the PBL [m]
  real(RP), intent(in) :: SFC_Z0M (IA,JA)      ! surface roughness length (momentum) [m]
  real(RP), intent(in) :: SFC_Z0H (IA,JA)      ! surface roughness length (heat) [m]
  real(RP), intent(in) :: SFC_Z0E (IA,JA)      ! surface roughness length (vapor) [m]
  real(RP), intent(in) :: SFC_DENS(IA,JA,TA)   ! density     at the surface atmosphere [kg/m3]
  real(RP), intent(inout):: SFLX_MW (IA,JA,TA) ! surface flux for z-momentum    (area center)   [m/s*kg/m2/s]
  real(RP), intent(inout):: SFLX_MU (IA,JA,TA) ! surface flux for x-momentum    (area center)   [m/s*kg/m2/s]
  real(RP), intent(inout):: SFLX_MV (IA,JA,TA) ! surface flux for y-momentum    (area center)   [m/s*kg/m2/s]

  integer              :: i,j,t
  real(RP)             :: Ustar               ! friction velocity [m]
  real(RP)             :: Tstar               ! friction temperature [K]
  real(RP)             :: Qstar               ! friction mixing rate [kg/kg]
  real(RP)             :: Uabs                ! modified absolute velocity [m/s]
  real(RP)             :: Ra                  ! Aerodynamic resistance (=1/Ce) [1/s]
  real(RP)             :: SFC_QSAT            ! saturatad water vapor mixing ratio [kg/kg]
  real(RP)             :: SFC_QV              ! water vapor mixing ratio [kg/kg]
  real(RP)             :: SFC_PSAT (IA,JA,TA) ! saturatad water vapor pressure [Pa]
  real(RP)             :: FracU10             ! calculation parameter for U10 [-]
  real(RP)             :: FracT2              ! calculation parameter for T2 [-]
  real(RP)             :: FracQ2              ! calculation parameter for Q2 [-]
  real(RP) :: ATMOS_PHY_SF_BULK_beta = 1.0_RP   ! evaporation efficiency (0-1)

  !!
  call SCALE_init

  call BULKFLUX_setup( sqrt(DX**2+DY**2) )

  SFLX_MU = 0
  SFLX_MV = 0

  do t = TS, TE
    call SATURATION_psat_all( IA, IS, IE, JA, JS, JE, & ! [IN]
                                      SFC_TEMP(:,:,t), & ! [IN]
                                      SFC_PSAT(:,:,t)  ) ! [OUT]

    do j = JS, JE
    do i = IS, IE

      SFC_QSAT = EPSvap * SFC_PSAT(i,j,t) / ( SFC_PRES(i,j,t) - ( 1.0_RP-EPSvap ) * SFC_PSAT(i,j,t) )
      SFC_QV = ( 1.0_RP - ATMOS_PHY_SF_BULK_beta ) * ATM_QV(i,j,t) + ATMOS_PHY_SF_BULK_beta * SFC_QSAT

      call BULKFLUX( Ustar,           & ! [OUT]
                      Tstar,           & ! [OUT]
                      Qstar,           & ! [OUT]
                      Uabs,            & ! [OUT]
                      Ra,              & ! [OUT]
                      FracU10,         & ! [OUT]
                      FracT2,          & ! [OUT]
                      FracQ2,          & ! [OUT]
                      ATM_TEMP(i,j,t), & ! [IN]
                      SFC_TEMP(i,j,t), & ! [IN]
                      ATM_PRES(i,j,t), & ! [IN]
                      SFC_PRES(i,j,t), & ! [IN]
                      ATM_QV  (i,j,t), & ! [IN]
                      SFC_QV       ,   & ! [IN]
                      ATM_U   (i,j,t), & ! [IN]
                      ATM_V   (i,j,t), & ! [IN]
                      ATM_Z1  (i,j),   & ! [IN]
                      PBL     (i,j,t), & ! [IN]
                      SFC_Z0M (i,j),   & ! [IN]
                      SFC_Z0H (i,j),   & ! [IN]
                      SFC_Z0E (i,j)  )   ! [IN]

      SFLX_MU(i,j,t) = -SFC_DENS(i,j,t) * Ustar * Ustar / Uabs * ATM_U(i,j,t)
      SFLX_MV(i,j,t) = -SFC_DENS(i,j,t) * Ustar * Ustar / Uabs * ATM_V(i,j,t)

    end do
    end do
  end do

  call SCALE_finalize
  !stop

end subroutine sample
