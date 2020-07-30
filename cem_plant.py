import numpy as np
import math
import pandas as pd
from numpy.core import zeros
from volt import volt

def cem_plant(SMR, CC, chemical, Dry, echem, retro, burnH2, Sell_SCM, Sell_Iron, Sell_Alumina, Sell_Agrigate, cleanH,
              cleanE, DPkW, CO2_Tax, CH, mineCost, heatCost, maxSCM, eCost, kWhr_kg, r, CaO_Frac_Rock, Al2O3_Frac_Rock,
              SiO2_Frac_Rock, FeO_Frac_Rock, Fe2O3_Frac_Rock, MgO_Frac_Rock, CO2int, eCO2int, PPT_SCM, S_Cost, CF, TPY,
              SA_ratio, Eff, rev, PPT_F, PPT_Al, PPT_Agg, V, W, C_imp, Al_eff, Fe_eff, SCM_eff, OPC_eff, Agg_eff,
              CapEx_fac):
    #   ----------Plant Constraints----------
    ##print("DPkW: ", DPkW)
    ###print("CO2_Tax: ", CO2_Tax)
    TailShip = 1
    TPYcal = TPY * OPC_eff
    CO2int = ((not cleanH) * CO2int)
    eCO2int = ((not cleanE) * eCO2int)
    rock_fac = (1 - CaO_Frac_Rock - SiO2_Frac_Rock) / (Al2O3_Frac_Rock + FeO_Frac_Rock + Fe2O3_Frac_Rock + MgO_Frac_Rock)

    Al2O3_Frac_Rock = Al2O3_Frac_Rock * rock_fac
    FeO_Frac_Rock = FeO_Frac_Rock * rock_fac
    Fe2O3_Frac_Rock = Fe2O3_Frac_Rock * rock_fac
    MgO_Frac_Rock = MgO_Frac_Rock * rock_fac

    #   ----------Chemical Definitions----------
    #   ---OPC Definition (Type 5)---
    mFrac_C3S = 0.38  # mass fraction of C3S in clinker
    mFrac_C2S = 0.43  # mass fraction of C2S in clinker
    mFrac_C3A = 0.04  # mass fraction of C3A in clinker
    mFrac_C4AF = 0.09  # mass fraction of C4AF in clinker
    CaCO3_frac_dry = 0.56  # fraction of CaO in CaCO3

    #   ---Molar Masses---
    molMass_H2O = 18 / 1000000  # tonnes/mol H2O
    molMass_H2SO4 = 98 / 1000000  # tonnes/mol H2SO4
    molMass_SiO2 = 60 / 1000000  # tonnes/mol SiO2
    molMass_CaSO4 = 136 / 1000000  # tonnes/mol CaSO4
    molMass_F = 160 / 1000000  # tonnes/mol Fe2O3
    molMass_FO = 72 / 1000000  # tonnes/mol Fe2O3
    molMass_Al = 102 / 1000000  # tonnes/mol Al2O3
    molMass_Ca = 56 / 1000000  # tonnes/mol CaO
    molMass_Mg = 40 / 1000000  # tonnes/mol MgO
    molMass_CO2 = 44 / 1000000  # tonnes/mol MgO
    molMass_C3S = 228 / 1000000  # tonnes/mol C3S
    molMass_C2S = 172 / 1000000  # tonnes/mol C2S
    molMass_C3A = 270 / 1000000  # tonnes/mol C3A
    molMass_C4AF = 486 / 1000000  # tonnes/mol C4AF
    molMass_S = 32 / 1000000  # tonnes/mol Sulfur
    molMass_Al2SO43 = 342 / 1000000  # tonnes/mol Al2(SO4)3
    molMass_Fe2SO43 = 400 / 1000000  # tonnes/mol Fe2(SO4)3
    molMass_FeSO4 = 152 / 1000000  # tonnes/mol FeSO4
    molMass_MgSO4 = 120 / 1000000  # tonnes/mol MgSO4
    molMass_AlCl26H2O = 240 / 1000000  # tonnes/mol AlCl3.6H2O

    #   ----Calculate elemental oxide fraction in the cement----
    #   --percent CaO--
    
    CaO_Frac_Cem = molMass_Ca * 3 / molMass_C3S * mFrac_C3S + molMass_Ca * 2 / molMass_C2S * mFrac_C2S + molMass_Ca * 3 / molMass_C3A * mFrac_C3A + molMass_Ca * 4 / molMass_C4AF * mFrac_C4AF  
    
    

    # --percent SiO2--
    SiO2_Frac_Cem = molMass_SiO2 / molMass_C3S * mFrac_C3S + molMass_SiO2 / molMass_C2S * mFrac_C2S

    # --percent Al2O3--
    Al2O3_Frac_Cem = molMass_Al / molMass_C3A * mFrac_C3A + 102 / 486 * mFrac_C4AF

    # --percent Fe2O3--
    Fe2O3_Frac_Cem = molMass_F / molMass_C4AF * mFrac_C4AF

    molT_C3S = mFrac_C3S / molMass_C3S  # mol C3S /tonne clinker
    molT_C2S = mFrac_C2S / molMass_C2S  # mol C2S /tonne clinker
    molT_C3A = mFrac_C3A / molMass_C3A  # mol C3A /tonne clinker
    molT_C4AF = mFrac_C4AF / molMass_C4AF  # mol C4AF /tonne clinker

    mptCa_cem = CaO_Frac_Cem / molMass_Ca  # mol CaO per tonne OPC
    mptSi_cem = SiO2_Frac_Cem / molMass_SiO2  # mol SiO2 per tonne OPC
    mptAl_cem = Al2O3_Frac_Cem / molMass_Al  # mol Al2O3 per tonne OPC
    mptFe_cem = Fe2O3_Frac_Cem / molMass_F  # mol Fe2O3 per tonne OPC

    mptCa_Rock = CaO_Frac_Rock / molMass_Ca  # mols CaO in 1 T of mined rock
    mptAl_Rock = Al2O3_Frac_Rock / molMass_Al  # mols Al2O3 in 1 T of mined rock
    mptSi_Rock = SiO2_Frac_Rock / molMass_SiO2  # mols SiO2 in 1 T of mined rock
    mptFe_Rock = Fe2O3_Frac_Rock / molMass_F  # mols Fe2O3 in 1 T of mined rock
    mptFeO_Rock = FeO_Frac_Rock / molMass_FO  # mols Fe2O3 in 1 T of mined rock
    mptMg_Rock = MgO_Frac_Rock / molMass_Mg  # mols MgO in 1 T of mined rock

    # ----Masses----
    T_rock = mptCa_cem / mptCa_Rock  # total T rock needed

    # --conventional process mass of mining--
    mass_Dry = CaO_Frac_Cem / CaCO3_frac_dry + 1 - CaO_Frac_Cem
 
    # --conventional process CO2 emissions--
    mass_CO2 = mass_Dry - 1

    # --mass of SCM--
    # -Rock is turned into OPC + SCM except Fe and Mg fractions-
    mass_SCM = T_rock * (1 - (Fe2O3_Frac_Rock + FeO_Frac_Rock + MgO_Frac_Rock)) - (1 - Fe2O3_Frac_Cem)

    # --Rock can be oxidized into Fe2O3--
    mass_Fe = ((mptFe_Rock + mptFeO_Rock) * T_rock - mptFe_cem) * molMass_F
    mass_Al = ((mptAl_Rock) * T_rock - mptAl_cem) * molMass_Al
    mass_Agg = T_rock - 1

    # ----Acid and H2 needs----
    molH2SO4 = (3 * mptAl_Rock + mptMg_Rock + 3 * mptFe_Rock + mptFeO_Rock + mptCa_Rock) * T_rock  # mol sufuric acid per T OPC
    H_output = molH2SO4 * 2 / 1000 * TPYcal / 365  # kg H2/day
    HPT = H_output * 365 / TPYcal  # kg H2 per tonne clinker
    SAPT = molH2SO4 * molMass_H2SO4  # T H2SO4 per tonne clinker

    # ----------Thermodynamic calculations Chemical Process Leach via HCl----------

    # ---Heat of formation data---
    HCl = -92.3  # kJ/mol
    CaAl2Si2O8 = -4232.41  # kJ/mol
    SiO2 = -859.4  # kJ/mol
    CaCl2 = -795.0  # kJ/mol
    AlCl3 = -706.25  # kJ/mol
    H2OL = -285.8  # kJ/mol liquid water
    H2OG = -241.8  # kJ/mol water vapor
    Al2O3 = -1669.8  # kJ/mol
    CaO = -635.5  # kJ/mol
    H2SO4 = -811.3  # kJ/mol
    CaSO4 = -1437.2  # kJ/mol
    SO2 = -296.1  # kJ/mol
    C3S = -2929.2  # kJ/mol
    C2S = -2307.5  # kJ/mol
    C3A = -3588.6  # kJ/mol
    C4AF = -5092.89  # kJ/mol
    CaCO3 = -1207  # kJ/mol
    CO2 = -394  # kJ/mol
    Al2SO43 = -3440  # kJ/mol
    Al2Si2O7 = -3211  # kJ/mol
    CaSiO3 = -1630  # kJ/mol
    Fe2SO43 = -2583  # kJ/mol
    FeSO4 = -929  # kJ/mol
    MgSO4 = -1278.2  # kJ/mol
    Fe2O3 = -826.0  # kJ/mol
    MgO = -601.7  # kJ/mol
    MgO2H2 = -924.7
    MgCl2 = -641.8
    MgClOH = -783.25

    Mg2SiO4 = -2176.93  # kJ/mol
    Fe2SiO4 = -1478  # kJ/mol

    '''first leach with HCl, in water, this happens in water so a max temp of 105 C (boiling point of fuming HCl), excess
    HCl is either added in the begining or sequentially which will cause AlCl3 to precipitate out of solution
    with SiO2
    '''

    DH_AnDecomp = (4 * H2OL + 2 * AlCl3 + CaCl2 + 2 * SiO2 - CaAl2Si2O8 - 8 * HCl) * T_rock * mptAl_Rock  # kJ / T OPC
    DH_WoDecomp = (CaCl2 + SiO2 - CaSiO3 - 2 * HCl + H2OL) * (T_rock * (mptCa_Rock - mptAl_Rock))  # kJ / T OPC

    ''' Fe and mg is also leached out of solution, however thier heat of formation is not considered to a paucity of 
        data, this is therefore an underestimate of the heat produced
    '''
    #   Dry alum is produced from AlCl3 by pouring H2SO4 on the AlCl3
    DH_AS = (Al2SO43 - 2 * AlCl3 - 3 * H2SO4 + 6 * HCl) * T_rock * mptAl_Rock

    # --Gypsum is selectively precipitated, entropic heat of mixing is also not
    #   considered here, meaning this is a conservative estimate--
    DH_GypProd = (CaSO4 + 2 * HCl - H2SO4 - CaCl2) * mptCa_cem

    # --Gypsum is decomposed--
    DH_GypDecomp = (CaO + SO2 - CaSO4) * mptCa_cem

    ###########################################################
    # DH_GypDecomp = (CaO + H2SO4 - CaSO4 - H2OG) * mptCa_cem #
    ###########################################################

    ''' At this point Mg and Fe are left in solution as sulfate salts, they thenneed to be decomposed to regenerate the 
        sulfuric acid Iiron sulfates decompose around 500C, so this can be done with the contact process:
        https: // onlinelibrary.wiley.com / doi / pdf / 10.1002 / jctb.2720220505
    '''
    DH_Fe2O3Prod_cem = (Fe2O3 * 0.5 + SO2 - FeSO4) * mptFe_cem * 2
    DH_Fe2O3Prod = (Fe2O3 * 0.5 + SO2 - FeSO4) * mptFeO_Rock * T_rock - DH_Fe2O3Prod_cem

    ###################################################################
    #   DH_Fe2O3Prod_cem = ((Fe2O3 * 0.5 + H2SO4 - H2OG - FeSO4 - SiO2 / 2) * mptFeO_Rock * T_rock + ...
    #                      (Fe2SiO4 + 3 * (H2SO4 - H2OG) - Fe2SO43 - SiO2) * mptFe_Rock * T_rock) * mptFe_cem / (
    #                      mptFe_Rock + mptFeO_Rock);
    #
    #    DH_Fe2O3Prod = (Fe2O3 * 0.5 + H2SO4 - H2OG - FeSO4) * mptFeO_Rock * T_rock + ...
    #                   (Fe2O3 + 3 * (H2SO4 - H2OG) - Fe2SO43) * mptFe_Rock * T_rock -
    #                   DH_Fe2O3Prod_cem;
    ###################################################################
    DH_MgOProd = (MgO - MgSO4 + H2SO4 - H2OG) * T_rock * mptMg_Rock

    ##################################################################################
    DH_MgOProd = (Mg2SiO4 - 2 * MgSO4 - SiO2 + 2 * SO2) * T_rock * mptMg_Rock / 2
    ##################################################################################

    ''' The Iron oxide thermally decomposes at < 700C while the MgSO4 thermally deocmposes at around 1200C which means it
        is possible to separate the Fe via magnetif separation
    '''
    #   Kaolinite can be produced from alum for SCM or OPC
    KaolProd = (Al2Si2O7 + 3 * SO2 - Al2SO43 - SiO2) * (T_rock * mptAl_Rock - mptAl_cem) * (not Sell_Alumina)

    '''-------------------------------------------------------------------------------------------------------------#
        KaolProd = ((Al2O3 + 3 * (H2SO4 - H2OG) - Al2SO43) * (T_rock * mptAl_Rock - mptAl_cem) * (not Sell_Alumina))
       --------------------------------------------------------------------------------------------------------------#
    '''

    #   or Alumina can be produced
    DH_Al2O3_Prod = (Al2O3 + 3 * (H2SO4 - H2OG) - Al2SO43) * T_rock * mptAl_Rock * Sell_Alumina
    DH_Al2O3_Prod_cem = (Al2O3 + 3 * SO2 - Al2SO43) * mptAl_cem

    # ---Heat from Making Clinker---
    DH_C2S = (C2S - 2 * CaO - SiO2) * molT_C2S  # heat of formation C2S
    DH_C3S = (C3S - 3 * CaO - SiO2) * molT_C3S  # heat of formation C3S
    DH_C3A = (C3A - 3 * CaO - Al2O3) * molT_C3A  # heat of formation C3A
    DH_C4AF = (C4AF - 4 * CaO - Al2O3 - Fe2O3) * molT_C4AF  # heat of formation C4AF

    # ---Heat from thermally decomposing limestone---
    DH_Lime = (CaO + CO2 - CaCO3) * mptCa_cem  # heat of formation CaO from CaCO3

    # ---Acid Regeneration---
    #   Sulfuric acid is regenerated via the contact process
    S_WasteCost = (mptMg_Rock) * T_rock * molMass_S * S_Cost

    #   If SO2 was regenerated, how much energy is needed?
    S_RegenHeatCost = (DH_MgOProd / 1000000) * heatCost

    lossFrac = 0.1  # SO2 losses, that need to be replaced with sulfur

    # ---Contact process regenerates the H2SO4---
    DH_CP = (H2SO4 - SO2 - H2OL) * T_rock * (3 * mptAl_Rock + mptMg_Rock + 3 * mptFe_Rock + mptCa_Rock + mptFeO_Rock) * (not echem)

    # ------------------------
    DH_CP = 0
    # ------------------------

    #   If it is cheap, it may be good to waste the MgSO4 and buy sulfur
    if S_WasteCost < S_RegenHeatCost:
        S = mptMg_Rock * T_rock
    else:
        S = 0

    # ---Heat from burning sulfur including waste sulfur----
    SAmount = T_rock * (3 * mptAl_Rock + mptCa_Rock + 3 * mptFe_Rock + mptFeO_Rock) * lossFrac + S  # T Sulfur / T OPC
    DH_SB = SO2 * SAmount  # kJ / T OPC

    # -------------------------------------------------
    DH_SB = 0
    #DH_Fe2O3Prod = 0
    # -------------------------------------------------

    LH_Brimstone = 0

    if DH_Fe2O3Prod < 0:
        LH_Brimstone = DH_Fe2O3Prod / 1000000 + LH_Brimstone
        DH_Fe2O3Prod = 0

    if DH_MgOProd < 0:
        LH_Brimstone = DH_MgOProd / 1000000 + LH_Brimstone
        DH_MgOProd = 0

    DH_Brimstone = ( DH_GypDecomp + ( DH_C2S +  DH_C3S +  DH_C4AF +  DH_C3A) +  DH_Al2O3_Prod_cem +  DH_Fe2O3Prod_cem) / 1000000 + ( KaolProd +  DH_MgOProd +  DH_Fe2O3Prod) / 1000000
  
    # ----Sensible and Latent Heat requirements----
    # ---Define some constants---
    cp_rock = 0.002  # heat capacity of rock, GJ/T
    cp_H2O = 0.004  # heat capacity of water, GJ/T
    DT_H2O = 80  # delta T, heat to evaporate water degs C
    LV_H2O = 2.3  # latent heat of evaporating water GJ/T;
    DT_clinker = 850  # delta T to heat 600 to 1450 to make clinker
    DT_Kaolin = 580  # delta T to heat to froom 20 to 600 to make metakaolin delta T to thermally decompose iron to from 20 to 600 to make metakaolin
    DT_iron = 480
    GJ_H2 = 0.12  # GJ/kg H2
    SH_eff = 0.3  # 30% of sensible heat is lost each time
    DT_dryK = 780  # delta T for drying kaolinite
    KW = .13  # natural kaoinite is 13% water
    DT_cal = 880  # Delta T for the precalciner
    DT_clinker_dry = 550  # Delta T to bring to clinker temp for the dry process

    # ---Sensible heat required to make clinker in the brimstone process (GJ/T OPC)---
    SH_Cinker = (cp_rock * (mptAl_cem * molMass_Al + mptFe_cem * molMass_F + mptSi_cem * molMass_SiO2 + mptCa_cem * molMass_CaSO4) * DT_clinker) * SH_eff

    SH_CinkerPreheat = (cp_rock * (mptAl_cem * molMass_Al2SO43 + mptFe_cem * molMass_Fe2SO43 + mptSi_cem * molMass_SiO2 + mptCa_cem * molMass_CaSO4) * DT_Kaolin) * SH_eff

    SH_Iron = cp_rock * ((molMass_Fe2SO43 * mptFe_Rock + molMass_FeSO4 * mptFeO_Rock + molMass_MgSO4 * mptMg_Rock) * T_rock - mptFe_cem * molMass_Fe2SO43) * DT_iron * SH_eff

    SH_Al = cp_rock * (mptAl_Rock * T_rock - mptAl_cem) * (molMass_Al2SO43) * DT_Kaolin * SH_eff

    SH_kaolin = ((cp_rock * (mptAl_Rock * T_rock - mptAl_cem) * (molMass_SiO2 * 2) * DT_Kaolin) * SH_eff)

    SHBrimstone = (SH_Cinker + SH_kaolin * Sell_SCM + SH_Al * Sell_Alumina + SH_CinkerPreheat + SH_Iron)
  
    #   fracation of free water left in precipitate (GJ/T OPC)
    FreeWater = 0.1  # some people say 10% is possible, we assume 40%

    #   AlCl3.6H2O precipitates with 15% free water as determined experimentally
    AlWater = 0.1

    #   free water that needs to be boiled off of the precipitants plus water in solution
    Fe_Sol = 3  # ratio of sulfate salts to water (typical solubility)

    M_solution =   (mptCa_Rock * molMass_CaSO4 + mptSi_Rock * molMass_SiO2) * T_rock * FreeWater + (mptAl_Rock * molMass_AlCl26H2O * 2) * T_rock * AlWater + (mptFe_Rock * molMass_Fe2SO43 + mptFeO_Rock
                    * molMass_FeSO4 + mptMg_Rock * molMass_MgSO4) * T_rock / Fe_Sol 

    #   aluminum chloride forms a hexahydrate, so 12 waters need to be boiled off at these temps,
    #   CaSO4 will not form a hydrate, so it is not included
    AlumAd = 12  # 12 waters in a 2 AlCl3 aducts (where alum comes from)
    M_Water = (mptAl_Rock * AlumAd) * T_rock * molMass_H2O + M_solution

    ''' determine if water needs to be added based on amount of water in H2SO4
    if  M_Water <(SAPT/SA_ratio-SAPT)*echem M_Water = (SAPT/SA_ratio-SAPT)*echem
    elseif W < M_Water:
        M_Water = W
    '''
    # ----Sensible heat needed to heat water to boiling point----
    SH_water = M_Water * cp_H2O * DT_H2O  # latent heat to boil water
    LH_Brimstone = LH_Brimstone + M_Water * LV_H2O

    #   iron sulfate is added here because it decomposes at contact process temps
    #   https://onlinelibrary.wiley.com/doi/pdf/10.1002/jctb.2720220505

    QWater = (SH_water + LH_Brimstone + (DH_AnDecomp + DH_WoDecomp + DH_GypProd + DH_AS) / 1000000) * chemical + (DH_CP) / 1000000

    ''' Heat that is low temp enough to be provided by the contact process
        heat needed for brimstone process, water cannot be recycled so if contact
        process is used, it mught be outside the efficiency marker
    '''
    if QWater < 0:
        QWater = 0

    #   Add sensible, latent, and reaction heat together

    QBrimstone = DH_Brimstone + SHBrimstone + QWater - HPT * GJ_H2 * burnH2 * echem

    if QBrimstone < 0:
        QBrimstone = 0

    QWCem = (mptCa_cem * molMass_CaSO4 + mptSi_cem * molMass_SiO2) * FreeWater + ( mptAl_cem * molMass_AlCl26H2O * 2) * AlWater + (mptFe_cem * molMass_Fe2SO43) / Fe_Sol + (mptAl_cem * AlumAd) * molMass_H2O * (cp_H2O * DT_H2O + LV_H2O)

    QWAlumina = (((mptAl_Rock * T_rock - mptAl_cem) * molMass_AlCl26H2O * 2) * AlWater + ((mptAl_Rock * T_rock - mptAl_cem) * AlumAd) * molMass_H2O) * (cp_H2O * DT_H2O + LV_H2O)

    QWIron = (mptFe_Rock * molMass_Fe2SO43 + mptFeO_Rock * molMass_FeSO4 + mptMg_Rock * molMass_MgSO4) * T_rock / Fe_Sol - (mptFe_cem * molMass_Fe2SO43) / Fe_Sol

    QWSCM = SH_water + LH_Brimstone - QWCem - QWAlumina * Sell_Alumina - QWIron * Sell_Iron

    QBrimCem = (DH_GypDecomp + DH_C2S + DH_C3S + DH_C4AF + DH_C3A + DH_SB + DH_Al2O3_Prod_cem + DH_Fe2O3Prod_cem) / 1000000 + (SH_Cinker + SH_CinkerPreheat) * SH_eff + QWCem * echem

    QBrimSCM = KaolProd / 1000000 + SH_kaolin * SH_eff * (not Sell_Alumina) + QWSCM * echem
    QBrimAl = (DH_Al2O3_Prod / 1000000 + SH_Al * SH_eff + QWAlumina * echem) * Sell_Alumina
    QBrimFe = ((DH_Fe2O3Prod + SH_Iron * SH_eff) / 1000000 + QWIron) * Sell_Iron * echem

    #   return variables for graphing
    QBrim = np.array([DH_Brimstone, SHBrimstone, QWater])  # create an array
    QH2 = np.array([-HPT * GJ_H2 * burnH2 * echem, 0, 0])

    # ----Conventional Cement Production----
    SH_clay_norm = mass_SCM * (KW + 1) * DT_dryK * cp_rock * SH_eff
    LH_clay_norm = mass_SCM * KW * LV_H2O
    QClay_norm = SH_clay_norm + LH_clay_norm
    SH_cem_norm = (cp_rock * (DT_cal * mass_Dry + DT_clinker_dry) + (
            DH_C2S + DH_C3S + DH_C4AF + DH_C3A) / 1000000) * SH_eff
    QNormal = (DH_Lime) / 1000000 + SH_cem_norm + QClay_norm * CC * Dry
  
    en = QBrimstone / QNormal  # ratio of normal energy to Brimstone energy
 
    # ---Calcualte data to return---
    Alumina_HPT = 14
    Alumina_EPT = 150
    Iron_HPT = 0.012449659 / 0.74  # https://www.energy.gov/sites/prod/files/2013/11/f4/iron.pdf
    Al_conv = (Al2O3_Frac_Rock * T_rock - Al2O3_Frac_Cem) * Alumina_HPT * Sell_Alumina
    Fe_conv = ((Fe2O3_Frac_Rock + FeO_Frac_Rock) * T_rock - Fe2O3_Frac_Cem) * Iron_HPT * Sell_Iron

    SH_SMR = 0.0162  # GJ/kg H2
    LH_SMR = 0.0647  # GJ/kg H2
    RH_SMR = 0.0809  # GJ/kg H2

    # ---create input variables---
    input_a = ((DH_Lime) / 1000000 + RH_SMR * HPT * (Dry and SMR or echem and (not Dry))) + Al_conv + Fe_conv
    input_b = (SH_cem_norm + SH_clay_norm + SH_SMR * HPT * (Dry and SMR or echem and (not Dry)))
    input_c = (LH_clay_norm + LH_SMR * HPT * (Dry and SMR or echem and (not Dry)))
    QDry = np.array([input_a, input_b, input_c])

    '''
        CapEx of a sufuric acid plant scales acocording to the below relationship this is found in the following 
        reference. Our process needs heat exchangers and electricity generation similar to a modern contact process 
        plant, electricity generation, and electricity generation for 9.1889 GJ per tonne, therefore we used a contact
        process of equivelent size for the heat generation (which is larger than the acid needs, to model our electricity
        generation CapEx, contact process plant
    '''
    # ---Conservative---
    #   https://www.researchgate.net/publication/290390987_Technical_Cost_Comparison_of_Laterite_Treatment_Processes_Part_3
    #   include gas processing

    TSAPD = molH2SO4 * molMass_H2SO4 / 365 * TPYcal
    Gas_Process_Fac = 0.4
    CP_CapEx = 1.9229 * (TSAPD * molMass_S / molMass_H2SO4) ** 0.5629 * 1000000
    CP_CapEx = CP_CapEx + CP_CapEx * Gas_Process_Fac

    # ---Liberal---
    #   https://www.researchgate.net/publication/286199169_Design_of_a_Plant_to_Manufacture_Sulfuric_Acid_from_Sulfur
    #   CP_CapEx = 0.277524*TSAPD^0.5655*1.02^30*1000000

    #   Sulfuric Acid Agitation tank, hold 4.8 M T rock per day.
    #   https://kingriverresources.com.au/wp-content/uploads/austocks/krr/2019_03_21_KRR_1553127420.pdf
    RefTank_orMagSepCost = 21000000  # USD
    RefTankSize = 3.9452e+04  # T/day
    BrimTankSize = T_rock * TPY / 365 / 24

    #   one for HCl, one for H2SO4, one magnetic seporator (3 total)
    TankCost = BrimTankSize / RefTankSize * RefTank_orMagSepCost * 3 * (not Dry)

    #   http://www.sulphuric-acid.com/Sulphuric-Acid-on-the-Web/COM2010_Acid_SC1_Plant_Fundamentals_Louie.pdf
    if echem == 1:
        CP_CapEx = CP_CapEx * Gas_Process_Fac

    #   mass ratio of raw mixes Brimstone/convention
    if Dry == 0 or (Dry == 1 and CC == 1):
        massRatio = T_rock / mass_Dry
    else:
        massRatio = 1

    #   CapEx for the Dry Cement Process source: https://ieaghg.org/docs/General_Docs/Reports/2008-3.pdf
    #   We converted this to USD with a 1.29 EUR/USD exchange rate and a 2% annual inflation

    #   DryCapEx = 274660*TPY^0.5189

    ex = 1.29  # the exchange rate used to convert reference materials into USD

    RawCrush = (2.5 + 1.25 + 2.08) * ex * (massRatio - 1 * retro)
    Mill = (16.5 + 0.5 + 1.9 + 2 + 0.5 + 0.5 + 5.3) * ex * (massRatio - 1 * retro)
    Preheat = (5.4 + 0.8) * ex * (massRatio * OPC_eff - 1 * retro)
    PreCal = (0.5) * ex
    Kiln = (12) * ex * (massRatio * OPC_eff - 1 * retro)
    ClinkCool = (12 + 2 + 0.5 + 1.5 + 9 + 0.3) * ex * (massRatio * OPC_eff - 1 * retro)

    if Dry == 1:
        CoalPrep = (5) * ex
        PetCokePrep = (5) * ex

    else:
        CoalPrep = (5) * ex * en
        PetCokePrep = (5) * ex * en

    CemMill = (10 + 10) * (ex) * (not retro)
    PackAndLoad = (5 + 5 + 3) * (ex * (not retro))

    mass = mass_Dry

    if massRatio > 1.1:
        mass = T_rock

    EquipCapEx = RawCrush + Mill + Preheat + PreCal + Kiln + ClinkCool + CoalPrep + PetCokePrep + CemMill + PackAndLoad
 
    DesEng = (42) * ex
    Construc = (48) * ex
    Other = (8) * ex * (not retro)
    EPC = (17) * ex
    Installed = DesEng + Construc + Other + EPC + EquipCapEx
 

    TailCapEx = TPYcal * 8 / 1000000 * mass * TailShip
    MineCapEx = TPYcal * 0.5 / 1000000 * mass * TailShip
    Contig = (12) * ex
    Fees = (5) * ex
    OwnCost = (12) * ex * (not retro)

    CemCapEx = Installed + Contig + Fees + OwnCost + TailCapEx + MineCapEx
 
    CemCapEx = CemCapEx * 1000000 * TPYcal / 1000000  # linear cost relation with size

    #---Water electrolysis, from H2A model---
    CF_e = 1
    WE_V = 2
    WEDpkW = 53
    V_ratio = V / WE_V
    CD = volt(V, C_imp)
    CD_WE = 1.4
    CD_ratio = CD_WE / CD
    BaseHPT = 15000000
    InstallFac = 1.12
    Stacks = WEDpkW * (HPT * TPYcal * 20 / (365 * 24 * CF_e)) / 1000000
    GasManagement = 12.5 / CF_e
    WaterDeliver = 3 / CF_e
    Thermal = 5.5 / CF_e
    PowerElectronics = 22.5 * V_ratio / CF_e
    MechBoP = 5.5 / CF_e
    OtherBoP = 3
    lyser_CapEx = (Stacks + GasManagement + WaterDeliver + Thermal + PowerElectronics + MechBoP + OtherBoP) * InstallFac / BaseHPT * H_output * 365 * 1000000 * echem

    SMR_CapEx = 174611028 / 124638630
    SMR_CapEx = SMR_CapEx * HPT * TPYcal

    CapEx = (CemCapEx + lyser_CapEx * echem * (not Dry) + CP_CapEx * (not Dry) + TankCost *
             (not Dry) + SMR_CapEx * SMR) * CapEx_fac
 
    ePT = 100 * massRatio # electricity demand per tonne of cement for plant in kWhr / tonne clinker

    eH2 = H_output * 365 * 15 * (V / 0.56) / TPYcal *(not Dry) * echem # kWhr / T

    CO2 = mass_CO2 * Dry + (ePT + eH2 * echem) * eCO2int + CO2int * QBrimstone * \
          (not Dry) + CO2int * QNormal * Dry + HPT * SMR * 10 / 1000 * Dry

    # create input variables
 

    GHG = np.array([
        mass_CO2 + HPT * 10 / 1000 * ((SMR and Dry) or (echem and (not Dry))) * .5, # 1.
                 (ePT) * eCO2int, # 2.
                 (ePT + eH2 * echem) * eCO2int * (not Dry), #3.
                 CO2int * (QNormal + QClay_norm) + HPT * 10 / 1000 * ((SMR and Dry) or (echem and (not Dry))) * .5, #4.
                 CO2int * QBrimstone * (not Dry)]) #5.

    ePT_rock = ePT / T_rock

    GHG_cem = QBrimCem * CO2int + ePT_rock * eCO2int

    GHG_Al = (QBrimAl * CO2int + (Al2O3_Frac_Rock * T_rock - Al2O3_Frac_Cem) * ePT_rock * eCO2int) * Sell_Alumina

    GHG_Fe = (QBrimFe * CO2int + ((Fe2O3_Frac_Rock + FeO_Frac_Rock) * T_rock - Fe2O3_Frac_Cem) *
              ePT_rock * eCO2int) * Sell_Iron

    GHG_SCM = (QBrimSCM * CO2int + (T_rock - (((Fe2O3_Frac_Rock +  FeO_Frac_Rock) * T_rock - Fe2O3_Frac_Cem) +
            (Al2O3_Frac_Rock * T_rock - Al2O3_Frac_Cem) + 1)) * ePT_rock * eCO2int) * Sell_SCM
    GHG_H2 = eH2 *(not Dry) * echem * eCO2int

    GHG_Al_conv = ((Al2O3_Frac_Rock * T_rock - Al2O3_Frac_Cem) * Alumina_HPT * CO2int + Alumina_EPT * eCO2int) * Sell_Alumina
    GHG_Fe_conv = (((Fe2O3_Frac_Rock + FeO_Frac_Rock) * T_rock - Fe2O3_Frac_Cem) * Iron_HPT * CO2int) * Sell_Iron
    GHG_H2_conv = HPT * 10 / 1000 * echem
    GHG_Div = np.array([GHG_cem, GHG_Al, GHG_Fe, GHG_SCM, GHG_H2, GHG_Al_conv, GHG_Fe_conv, GHG_H2_conv]) * (not Dry)

    #   OpEx
    water_cost = 1 # $ / t US EPA
    water = molH2SO4 * molMass_H2O * (1 + 1 * echem) * water_cost * TPYcal

    #   if we make excess energy, convert it into electricity O & M rate as percent of CsapEx for cement industry(same
    #   source as cement CapEx)
    OMcem = 0.02
    #   O & M rate as percent of CsapEx for PEM industry(H2A source)
    OMpem = 0.03

    if Dry == 1 and CC == 1:
        #print("MINECOST == 15")
        mineCost = 15 # account for multiple mines in multiple locations

    SMR_OpEx = 1.15 * HPT * TPYcal # from H2A model

    TailCost = 0.39 # $ / T

    if Dry == 1:
        Tailings = TPYcal * 41 / 1000 * TailCost * TailShip
        
        Ship = mass_Dry * TPYcal * 1.5 * TailShip
        
        OpEx = OMcem * CemCapEx + mass_Dry * mineCost * TPYcal + (T_rock - 1) * TPYcal * mineCost * CC + QNormal * \
               TPYcal * heatCost + (ePT) * eCost * (TPYcal + (T_rock - 1) * CC) + CO2 * TPYcal * CO2_Tax + SMR_OpEx *\
               SMR + Tailings + Ship
    
        OpExMat = np.array([
            (OMcem * CemCapEx), #1.
            (mass_Dry * mineCost * TPYcal + (T_rock - 1) *TPYcal * mineCost * CC), #2.
            QNormal * TPYcal * heatCost, #3.
            (ePT) * eCost * (TPYcal + (T_rock - 1) * CC), #4.
            CO2 * TPYcal * CO2_Tax + SMR_OpEx * SMR, #5.
            0, #6.
            0, #7.
            Tailings, #8.
            Ship]) #9.
    else:
        Tailings = TPYcal * 41 / 1000 * TailCost * T_rock * TailShip
        
        Ship = T_rock * TPYcal * 1.5 * TailShip
        
        OpEx = OMcem * CemCapEx + OMpem * (lyser_CapEx + CP_CapEx) + T_rock * mineCost * TPYcal + QBrimstone * TPYcal\
               * heatCost + water + (ePT + eH2 * echem) * eCost * TPYcal + CO2 * TPYcal * CO2_Tax + S_Cost * SAmount +\
               Tailings + Ship
        
        OpExMat = np.array([
            (OMcem* CemCapEx + OMpem * (lyser_CapEx + CP_CapEx)),
            (T_rock * mineCost * TPYcal), QBrimstone * TPYcal * heatCost,
            ((ePT + eH2 * echem) * eCost * TPYcal),
            CO2 * TPYcal * CO2_Tax, water,
            S_Cost * SAmount,
            Tailings,
            Ship])

    Lifetime = 25 # years

    life_mat = (1 + r) ** np.arange(2, Lifetime + 2, 1) # (2:1:(Lifetime+1))  assume 1 year build time

    OpEx_tot = OpEx * 0
  
    OpEx_tot = OpEx_tot + OpEx * np.sum(1 / life_mat)
   
    OpExMat_tot = np.zeros(len(OpExMat), dtype=float) # CHECK CONVERSION: WHAT DATA TYPE DO YOU WANT THIS TO BE?
    
    OpExMat_tot = OpExMat_tot + OpExMat * np.sum(1 / life_mat)
    Cap_sum = np.sum((TPYcal)/ life_mat)
    Prod_rev = CapEx + OpEx_tot
    
    mass_SCMnew = 0
    if Dry == 1 and CC == 0 and SMR == 0:
        LCC = Prod_rev/ Cap_sum
        SCM = 0
        H2_value = 0
        CapExPT = CapEx / Cap_sum
        OpExPT = OpEx_tot / Cap_sum
    else:
        mass_SCMnew = mass_SCM

    if mass_SCMnew * SCM_eff > maxSCM:
        mass_SCMnew = maxSCM
    elif Dry == 1 and CC == 0:
        mass_SCMnew = 0


    SCM = (mass_SCMnew - mass_Al * Sell_Alumina) * PPT_SCM * Sell_SCM
    Fe = 0
    Al = 0
    Agg = 0

    if Dry == 0:
        Fe = mass_Fe * PPT_F * Sell_Iron * Fe_eff
        Al = mass_Al * PPT_Al * Sell_Alumina * Al_eff
        Agg = (mass_Agg - mass_SCMnew * Sell_SCM - mass_Fe * Sell_Iron - mass_Al * Sell_Alumina) * PPT_Agg * \
              Sell_Agrigate * Agg_eff

    H2_value = 0

    if (echem == 1 and Dry == 0 and burnH2 == 0) or (SMR == 1 and Dry == 1):
        H2_value = HPT * CH

    if rev < 0:
        LCC = Prod_rev / Cap_sum - H2_value - SCM - Fe - Al - Agg
    else:
        LCC = Prod_rev / Cap_sum - rev

    if Dry == 1 and CC == 0:
        SCM = 0

    #   $25 / T
    CapEx_Fe = mass_Fe * 240 * TPYcal # https: // www.investopedia.com / articles / investing / 030215 / how - iron - ore - market - works - supply - market - share.asp
    CapExMat = np.array([CemCapEx, lyser_CapEx, CP_CapEx, TankCost, SMR_CapEx, CapEx_Fe]) * CapEx_fac
    CapExPT = CapEx / Cap_sum
    OpExPT = OpEx_tot / Cap_sum
    OpExMat = OpExMat_tot / Cap_sum

    return LCC, SCM, H2_value, CapExPT, OpExPT, Fe, QDry, QBrim, QH2, CapExMat, GHG, OpExMat, Al, Agg, GHG_Div