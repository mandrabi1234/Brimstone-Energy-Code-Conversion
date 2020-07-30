import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
from cem_plant import cem_plant
from define_cem_sens_anal import Def_anal
from replace import replace
from volt import volt

# check for OpEx
# do a tornado plot

# water us, efficiency of Ca extraction, current density

# The following variables define the process

SMR = 0
CC = 0 # % if 1, calcined clay production is included in the dry process
chemical = 1 # if 1 the chemical separation process is used
Dry = 0 # if 1 the dry process (conventional) is used
echem = 0 # if 1 electrochemistry is used to regenerate H2SO4
cleanE = 0
cleanH = 0
retro = 0 #  if 1 retrofit an existing plant
burnH2 = 0 # if 1 H2 is burned to heat the kiln
Skarn = 0 # if 1, mining a skarn
Ave_basalt = 0
Sell_SCM = 0
Sell_Iron = 0
Sell_Alumina = 0
Sell_Aggregate = 0
tornado = 0

# the following variabbles are possible to plot against each other vs price, by turning a value to 1
# you are choosing to compare these variables. It is only possible to compare two variables

Revenue_switch = 0
DPkW_switch = 0 # % vary the amount of H2 made per day
CO2_Tax_switch = 0 # vary the current density of the lyser
CH_switch = 0 # vary the buying price of sulfur
mineCost_switch = 1 # vary the selling price of sulfuric acid
heatCost_switch = 0 # vary the efficiency of solar panels
maxSCM_switch = 0 # vary the electricity cost
eCost_switch = 1 # price to sell electricity
kWhr_kg_switch = 0 # vary the cost of storing energy
r_switch = 0 # vary the capacity factor of the plant

#Thank you very much for all that you have done

CaO_Frac_Rock_switch = 0 # fraction of CaO in the source rock
Al2O3_Frac_Rock_switch = 0 # fraction of Al2O3 in the source rock
SiO2_Frac_Rock_switch = 0 # fraction of SiO2 in the source rock
FeO_Frac_Rock_switch = 0 
Fe2O3_Frac_Rock_switch = 0 # fraction of Fe2O3 in the source rock
MgO_Frac_Rock_switch = 0 # fraction of MgO in the source rock

CO2int_switch = 0 # T CO2/GJ heat
eCO2int_switch = 0 # T CO2/kWhr
PPT_SCM_switch = 0 # $/T SCM
PPT_Al_switch = 0 # $/T SCM
PPT_Fe_switch = 0 # $/T SCM
PPT_Agg_switch = 0 # $/T SCM
S_Cost_switch = 0 # $/T SCM
CF_switch = 0 # Kiln Capacity Factor
TPY_switch = 0 # Kiln Output clinker in tonnes per year
SA_ratio_switch = 0 # concentration of H2SO4 from the electrolyzer
Eff_switch = 0 # heat to power efficiency for organic rankine cycle, if used
PPT_F_switch = 0

W_switch = 0
V_switch = 0
CD_switch = 0
Al_eff_switch = 0
Fe_eff_switch = 0
SCM_eff_switch = 0
OPC_eff_switch = 0
Agg_eff_switch = 0
CapEx_fac_switch = 0

#------------------------------------------------------------
#CHECK CONVERSION: GO OVER THIS CODE SNIPPET WITH CODY
log_plot = 0 #  if 1 the price of H2 is base 10 logged
switches = [DPkW_switch, CO2_Tax_switch, CH_switch, mineCost_switch,
 heatCost_switch, maxSCM_switch, eCost_switch, kWhr_kg_switch, 
 r_switch, CaO_Frac_Rock_switch, Al2O3_Frac_Rock_switch, 
 SiO2_Frac_Rock_switch, FeO_Frac_Rock_switch, 
 Fe2O3_Frac_Rock_switch, MgO_Frac_Rock_switch, CO2int_switch, 
 eCO2int_switch, PPT_SCM_switch, S_Cost_switch, CF_switch, 
 TPY_switch, SA_ratio_switch, Eff_switch, Revenue_switch, 
 PPT_F_switch, PPT_Al_switch, PPT_Agg_switch, V_switch, W_switch, 
 CD_switch, Al_eff_switch, Fe_eff_switch, SCM_eff_switch, 
 OPC_eff_switch, Agg_eff_switch, CapEx_fac_switch]#
#------------------------------------------------------------


x = 0 #  variable to graph

y = 0 # variable to graph

indexx = -1
indexy = -1

# call the define_sens_anal function which returns a cell array of values
# for each variable and constant, two variables are compared at a time
# while everything is kept constant
sens_var, constants, SENS = Def_anal(switches, Skarn, Ave_basalt, tornado)


# find the index of the variables to compare in the sens_var cell array
# last 2 values are the constants that would be used if the variables
# selected were not variables, and therefore can be ignored


print(len(sens_var)-2)
if tornado == 0:
  for i in range(len(sens_var) - 2):
    if ((indexx == -1) and (indexy == -1) and ((len(SENS[i])) > 1)):
      x = sens_var[i] # defines the variable
      indexx = i      
      sens_var[i] = -1 # replaces variable with -1            
    elif ((indexy == -1) and ((len(SENS[i])) > 1)):
      y = sens_var[i]
      indexy = i
      sens_var[i] = -1              
      break

  length_x = len(x)
  length_y = len(y)
  


  z = np.zeros((length_y, length_x))
  z_CO2 = np.zeros((length_y, length_x))
  z_en = np.zeros((length_y, length_x))
  zz = np.zeros((length_y, length_x)) # number of components in system
  w = np.zeros((length_x, 1))
  v = np.zeros((length_y, 1))
  q = np.zeros((length_x))



  # run the plant code in a loop, updating the value of the variables (x and y) 
  # to be compared every time
  
  print("IndexX", indexx)
  print("IndexY", indexy)
  print("--X--: ", x)
  print("--Y--", y)
  a = 0 # indexes
  for i in x: #THIS MIGHT BE WRONG (CHECK WITH A TEST RUN)
    print("i: ", i)
    b = 0 # % indexes
    print("b: ", b)
    for j in y: #THIS MIGHT BE WRONG (CHECK WITH A TEST RUN)
      inputs = (sens_var) # inputs = cell2mat(sens_var) # 

      inputs[indexx] = i # % update the proper input with the next value
      inputs[indexy]= j # % update the proper input with the next value

      print("INPUTS1: ", inputs[indexx])
      print("INPUTS2: ", inputs[indexy])
      [cost, _, _, _, _, _, _, _, _, _, _, _, _, _, _] = cem_plant(SMR, CC, chemical, Dry, echem, retro, burnH2, 
      Sell_SCM, Sell_Iron, Sell_Alumina, Sell_Aggregate, cleanH, cleanE, 
      inputs[0], inputs[1], inputs[2], inputs[3], inputs[4], inputs[5], 
      inputs[6], inputs[7], inputs[8], inputs[9], inputs[10], 
      inputs[11], inputs[12], inputs[13], inputs[14], inputs[15], 
      inputs[16], inputs[17], inputs[18], inputs[19], 
      inputs[20], inputs[21], inputs[22], inputs[23], inputs[24], 
      inputs[25], inputs[26], inputs[27], inputs[28], inputs[29], 
      inputs[30], inputs[31], inputs[32], inputs[33], inputs[34], 
      inputs[35])
#------------------------------
# CHECK CONVERSION: GO OVER THIS CODE SNIPPET WITH CODY
      z[b, a] = cost #z(b, a) = cost # the cost H2 array
      print("----Z-VALUES CHECKER----")
      print("b: ", b)
      print("a: ", a)
      print(z[b,a])
      print("---------------------")
      print(z)
      
      np.nan_to_num(z) # z(isnan(z)) = 0 (replace all NaN values with 0)
      # z_CO2(b, a) = CO2 # the CO2 production array
      # z_en(b, a) = Energy_needed; #the energy use array
      # zz(b, a,:) = LCH_mat# % grid of costs array
#------------------------------
      b = b + 1 
      #  zz(isnan(zz)) = 0#
        
    a = a + 1 
  
# find Base Case
# update the proper input with the next value

  inputs[indexx] = sens_var[len(sens_var) - 2] # update the proper input with the next value
        
  inputs[indexy] = sens_var[len(sens_var)-1]

  [cost, SCM_value, H2_value, CapExPT, OpExPT, Iron, QDry, QBrim, QH2, 
  CapExMat, GHG, OpExMat, Al, Agg, GHG_Div] = cem_plant(SMR, CC, chemical, Dry,
  echem, retro, burnH2, Sell_SCM, Sell_Iron, Sell_Alumina, Sell_Aggregate, cleanH,
  cleanE, constants[0], constants[1], constants[2], constants[3], constants[4], 
  constants[5], constants[6], constants[7], constants[8], 
  constants[9], constants[10], constants[11], constants[12], 
  constants[13], constants[14], constants[15], constants[16], 
  constants[17], constants[18], constants[19], constants[20], 
  constants[21], constants[22], constants[23], constants[24], 
  constants[25], constants[26], constants[27], constants[28], 
  constants[29], constants[30], constants[31], constants[32], 
  constants[33], constants[34], constants[35])

# find Dry Case
# update the proper input with the next value

  inputs[indexx] = sens_var[len(sens_var) - 2] # update the proper input with the next value
  inputs[indexy] = sens_var[len(sens_var) - 1]
  [Dry_base, _, _, DryCapExPT, _, _, DryQDry, _, _, DryCapExMat, DryGHG, 
  DryOpExMat, _, _, _] = cem_plant(0, 0, chemical, 1, echem, retro, 0, 
  Sell_SCM, Sell_Iron, Sell_Alumina, Sell_Aggregate, cleanH, cleanE, 
  constants[0], constants[1], constants[2], constants[3],
  constants[4], constants[5], constants[6], constants[7],
  constants[8], constants[9], constants[10], constants[11], 
  constants[12], constants[13], constants[14], constants[15], 
  constants[16], constants[17], constants[18], constants[19], constants[20],
  constants[21], constants[22], constants[23], constants[24], constants[25], 
  constants[26], constants[27], constants[28], constants[29], constants[30], constants[31], constants[32], constants[33], constants[34], constants[35])

  zlabel = 'Levelized Cost of Cement/T'
  # log the cost if necessary, unless there are negative values
  if log_plot == 1 and (np.sum(np.sum(abs(z))) <= np.sum(np.sum(z))):
    z = math.log(z)#
    zlabel = 'log cost per kg of hydrogen'#

#------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------DATA VISUALIZATIONS-----------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------

  #------Figure 1------
  fig1 = plt.figure()
  #left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
  #ax = fig1.add_axes([left, bottom, width, height])

  cp = plt.contourf(x, y, z)
  plt.colorbar(cp)

  
  #ax.set_title('Contour Plot')
  plt.xlabel('Mining Cost (USD/T)')
  plt.ylabel('Electricity Cost ($/kWhr)')
  plt.show(fig1)

  #------Figure 2------
  fig2 = plt.figure(figsize = (6,5))
  left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
  ax = fig2.add_axes([left, bottom, width, height])

  data = {"Dry USD/T OPC":Dry_base, "Break-Even Co-Product Value":((CapExPT) + OpExPT - Dry_base), 
  "USD/T OPC":cost, "Value SCM":SCM_value, "Value Iron":Iron, "Value Aluminum":Al, "Value Aggregate":Agg, 
  "Value H_2":H2_value, "CapEx per T":(CapExPT), "OpEx pet T":OpExPT} 

  key = list(data.keys())
  Y = list(data.values())

  bar = plt.bar(key, Y, width = 0.5)
  plt.xticks(rotation = 45)
  ylabel = ("$/T OPC")

  for index, value in enumerate(Y):
      plt.text(index -0.5, value, str(value))
  plt.show(fig2)

  #------Figure 3------
  fig3 = plt.figure()
  rc('font', weight = 'bold')
  DryCapEx = np.array([DryCapExMat[0], 0, 0, 0])
  BrimCapEx = np.array([CapExMat[0] *(not Dry), CapExMat[1] * echem *(not Dry), CapExMat[2] *(not Dry), CapExMat[3] *(not Dry)])


  bars1 = [DryCapEx[0], BrimCapEx[0]]
  bars2 = [DryCapEx[1], BrimCapEx[1]]
  bars3 = [DryCapEx[2], BrimCapEx[2]]
  bars4 = [DryCapEx[3], BrimCapEx[3]]
  #bars = np.add(bars1, bars2, bars3).tolist()
  r = [0, 1]

  names = ['Dry Process CapEx', 'Brimstone CapEx']

  barwidth = 0.5

  plt.bar(r, bars1, label = 'Cement + SCM Plant', color = '#b5ffb9', edgecolor = 'white', width = barwidth)
  plt.bar(r, bars2, label = 'Hyrogen', bottom = bars1, color = '#557f2d', edgecolor = 'white', width = barwidth)
  plt.bar(r, bars3, label = 'Acid Regeneration', bottom = bars2, color = '#f9bc86', edgecolor = 'white', width = barwidth)
  plt.bar(r, bars4, label = 'Acid Reactor', bottom = bars3, color = '#a3acff', edgecolor = 'white', width = barwidth)


  plt.xticks(r, names, fontweight = 'bold')
  plt.xticks(rotation = 45)
  plt.ylabel("CapEx, 1MT OPC Plant (USD)")
  plt.legend()
  plt.show(fig3)

  #------Figure 4------
  fig4 = plt.figure(figsize = (6,5))
  rc('font', weight = 'bold')

  Q = DryQDry
  out = QBrim*(not Dry)

  bars1 = [Q[0], QDry[0], out[0], QH2[0]]
  bars2 = [Q[1], QDry[1], out[1], QH2[1]]
  bars3 = [Q[2], QDry[2], out[2], QH2[2]]

  bars = np.add(bars1, bars2).tolist()
  r = [0, 1, 2, 3]

  names = ['Dry Process Heat', 'Equivalent OPC + SCM + Fe_2O_3 + Al_2_O_3 + H_2 Heat', 'Brimstone Process Heat', 'H_2 Combustion Heat']
  barWidth = 0.5

  plt.bar(r, bars1, color = '#b5ffb9', edgecolor = 'white', width = barWidth, label = 'Net Reaction Heat')
  plt.bar(r, bars2, bottom = bars1, color = '#f9bc86', edgecolor = 'white', width = barWidth, label = 'Net Sensible Heat')
  plt.bar(r, bars3, bottom = bars, color = '#a3acff', edgecolor = 'white', width = barWidth, label = 'Net Latent Heat')

  plt.xticks(r, names, fontweight = 'bold', rotation = 45)
  plt.ylabel("Heat Energy (GJ/T OPC)")
  plt.legend()
  plt.show(fig4)

  #------Figure 5------
  fig5 = plt.figure(figsize = (6,5))
  rc('font', weight = 'bold')

  DryCapEx = np.array([DryGHG[0], DryGHG[3], DryGHG[1]])
  BrimCapEx = np.array([0, GHG[4]*(not cleanH), GHG[2]*(not cleanE)])
  ECapEx = np.array([GHG[0], (GHG[3]+GHG_Div[5]+GHG_Div[7])*(not cleanH), GHG[1]*(not cleanE)])

  bars1 = [DryCapEx[0], BrimCapEx[0], ECapEx[0]]
  bars2 = [DryCapEx[1], BrimCapEx[1], ECapEx[1]]
  bars3 = [DryCapEx[2], BrimCapEx[2], ECapEx[2]]

  bars = np.add(bars1, bars2).tolist()

  r = [0, 1, 2]

  names = ['Dry Process', 'Brimstone Total', 'Equivalent OPC + SCM + Fe_2O_3 + Al_2_O_3 + H_2']
  barWidth = 0.5

  plt.bar(r, bars1, color = '#b5ffb9', edgecolor = 'white', width = barWidth, label = 'Process')
  plt.bar(r, bars2, bottom = bars1, color = '#f9bc86', edgecolor = 'white', width = barWidth, label = 'Heat')
  plt.bar(r, bars3, bottom = bars, color = '#a3acff', edgecolor = 'white', width = barWidth, label = 'Electricity')

  plt.xticks(r, names, fontweight = 'bold', rotation = 45)
  plt.ylabel("CO_2 Intensity (T CO_2/TOPC)")
  plt.legend()
  plt.show(fig5)

  #------Figure 6------
  # NOTE: Discuss graph output and labeling with Cody via Google Meets
  fig6 = plt.figure()
  rc('font', weight = 'bold')

  DryCapEx = DryOpExMat
  BrimCapEx = np.array([OpExMat[0] * (not Dry), OpExMat[1] * (not Dry),OpExMat[2] * (not Dry), OpExMat[3] * (not Dry), OpExMat[4] * (not Dry), OpExMat[5] * (not Dry), OpExMat[6] * (not Dry), OpExMat[7] * (not Dry), OpExMat[8] * (not Dry)])

  bars1 = [DryCapEx[0], BrimCapEx[0]]
  
  bars2 = [DryCapEx[1], BrimCapEx[1]]
  bars3 = [DryCapEx[2], BrimCapEx[2]]
  bars4 = [DryCapEx[3], BrimCapEx[3]]
  bars5 = [0,0]
  bars6 = [DryCapEx[5], BrimCapEx[5]]
  bars7 = [DryCapEx[6], BrimCapEx[6]]
  bars8 = [DryCapEx[7], BrimCapEx[7]]
  bars9 = [DryCapEx[8], BrimCapEx[8]]
  #bars = np.add(bars1, bars2, bars3, bars4, bars5, bars6, bars7, bars8).tolist()
  #print("bars: ", bars)

  print("bars1: ", bars1)
  print("bars2: ", bars2)
  print("bars3: ", bars3)
  print("bars4: ", bars4)
  print("bars5: ", bars5)
  print("bars6: ", bars6)
  print("bars7: ", bars7)
  print("bars8: ", bars8)
  print("bars9: ", bars9)

  r = [0, 1]
  names = ['Dry Process OpEX', 'Brimstone OpEx']
  barWidth = 0.65

  plt.bar(r, bars1, color='red', label = 'O&M', edgecolor='white', width=barWidth)
  plt.bar(r, bars2, bottom = bars1, color='green', label = 'Mining', edgecolor='white', width=barWidth)
  plt.bar(r, bars3, bottom = bars2, color='blue', label = 'Heat', edgecolor='white', width=barWidth)
  plt.bar(r, bars4, bottom = bars3, color='yellow', label = 'Electricity', edgecolor='white', width=barWidth)
  plt.bar(r, bars5, bottom = bars4, color='orange', label = 'CO_2 Tax', edgecolor='white', width=barWidth)
  plt.bar(r, bars6, bottom = bars5, color='purple', label = 'Water', edgecolor='white', width=barWidth)
  plt.bar(r, bars7, bottom = bars6, color='teal', label = 'Sulfur', edgecolor='white', width=barWidth )
  plt.bar(r, bars8, bottom = bars7, color='pink', edgecolor='white', width=barWidth)
  plt.bar(r, bars9, bottom = bars8, color='grey', edgecolor='white', width=barWidth)

  plt.xticks(r, names, fontweight = 'bold', rotation = 45)
  plt.ylabel = ("OpEx/T OPC, 1MTPY OPC Plant (USD)")
  plt.legend()
  plt.show(fig6)

  #------Figure 7------
  fig7 = plt.figure(figsize = (6,5))
  rc('font', weight = 'bold')

  data = {"Conv. Process Cement":np.sum(DryGHG), "Brimstone Cement":GHG_Div[0], 
  "Conv. Alumina":GHG_Div[5], "Brimstone Alimina":GHG_Div[1], "Conv. Iron Oxide":GHG_Div[6], "Brimstone Iron Oxide":GHG_Div[2], 
  "Conv. SCM":(np.sum(np.array([GHG[0], GHG[4] * (not cleanH), GHG[1] * (not cleanE)])) - np.sum(DryGHG)) * Sell_SCM, 
  "Brimstone_SCM":GHG[3], "Conv. H_2": GHG_Div[7], "Brimstone H_2":GHG_Div[4]}

  key = list(data.keys())
  Y = list(data.values())

  bar = plt.bar(key, Y, width = 0.5)
  plt.xticks(rotation = 45)
  ylabel = ("TCO_2/T OPC produced")

  for index, value in enumerate(Y):
      plt.text(index -0.5, value, str(value))

  plt.show(fig7)


#--Continue If-Else Statement (if tornado == 0)--
else:
    t_data = np.arange(1, len(constants)+1, 1)
    t_data = [[i] for i in t_data]
    for i in range(len(constants)):
      costs = np.array([0, 0])
      inputs = constants
      
      inputs[i] = (sens_var[i])
      [costs[0], _, _, _, _, _, _, _, _, _, _, _, _, _, _] = cem_plant(SMR, CC, chemical, Dry, echem, retro, burnH2, 
      Sell_SCM, Sell_Iron, Sell_Alumina, Sell_Aggregate, cleanH, 
      cleanE, inputs[0], inputs[1], inputs[2], inputs[3], inputs[4], inputs[5], 
      inputs[6], inputs[7], inputs[8], inputs[9], inputs[10], 
      inputs[11], inputs[12], inputs[13], inputs[14], inputs[15], 
      inputs[16], inputs[17], inputs[18], inputs[19], 
      inputs[20], inputs[21], inputs[22], inputs[23], inputs[24], 
      inputs[25], inputs[26], inputs[27], inputs[28], inputs[29], 
      inputs[30], inputs[31], inputs[32], inputs[33], inputs[34], inputs[35])
    
      inputs[i] = (sens_var[i])#
    
      [costs[1], _, _, _, _, _, _, _, _, _, _, _, _, _, _] = cem_plant(SMR, CC, chemical, Dry, echem, retro, burnH2, 
      Sell_SCM, Sell_Iron, Sell_Alumina, Sell_Aggregate, cleanH, 
      cleanE, inputs[0], inputs[1], inputs[2], inputs[3], inputs[4], inputs[5], 
      inputs[6], inputs[7], inputs[8], inputs[9], inputs[10], 
      inputs[11], inputs[12], inputs[13], inputs[14], inputs[15], 
      inputs[16], inputs[17], inputs[18], inputs[19], inputs[20], 
      inputs[21], inputs[22], inputs[23], inputs[24], 
      inputs[25], inputs[26], inputs[27], inputs[28], inputs[29], 
      inputs[30], inputs[31], inputs[32], inputs[33], inputs[34], 
      inputs[35])

      t_data[i] = [costs]#

# Names of the Y-axis ticks 
    names = np.array(['Electrolyzer Cost (USD/kW) 50->1000', 'CO_2 Tax (USD/T CO_2) -190->100', 
    'H_2 price (USD/kg)1->4', 'Mining Cost (USD/T) 3->10', 'Cost of Heat (USD/GJ) 1->10', 
    'SCM Sold Per T OPC (T SCM/T OPC) 0.4->5.3', 'Electricity Cost ($/kWhr) 0.01->0.1', 'kWhr per kg H_2 5->30', 
    'IRR 4%->30%', 'CO_2 Intensity of Heat (T CO_2/GJ) 0->0.2', 'CO_2 Intensity of Electricity (T CO_2/kWhr) 0->0.001', 
    'Price per T of SCM (USD/T) 0->100', 'Sulfur Cost (USD/T) 0->200', 'Plant Size (T OPC/yr) 1M->5M', 'Co-product Revenue per T OPC 0->500', 
    'Price per Tonne Fe_2O_3 0->100', 'Price per Tonne Al_2O_3 0->400', 'Price per Tonne Aggregate 0->20', 'Operating Cell Voltage (V) 0.2->1.2', 
    'Water Content (T/T OPC) 0.1-5', 'Improvement in Current Density (A/cm^2) 0.4-> -0.4', 
    'Al production efficiency 50%->100%', 'Fe production efficiency 50%->100%', 'SCM production efficiency 50%->100%', 
    'OPC production efficiency 50%->90%', 'Aggregate production efficiency 50%->100%', 'CapEx multiplier 0.5->2'])

    # ' Fraction CaO in Source Rock', ...
    # ' Fraction Al_2O_3 in Source Rock', ' Fraction SiO_2 in Source Rock', ...
    # ' Fraction FeO in Source Rock', 'Freaction Fe_2O_3 in Source Rock', ...
    # 'Fraction MgO, in Source Rock',
    # 'Capacity Factor',
    # 'Conventration of Produced Sulfuric Acid', ...
    # 'Heat to Power Efficiency',
    var1 = np.arange(1, 10, 1)
    var2 = np.arange(16, 20, 1)
    var3 = np.array([21])
    var4 = np.arange(24, len(t_data) + 1, 1)
    vars = np.hstack((var1, var2, var3, var4))

    
    a = 1#
    high = np.zeros(vars.shape)#
    low = np.zeros(vars.shape)#
    for i in vars: 
      n = t_data[i]#
      high[a] = np.maximum(n)#
      low[a] = np.minimum(n)#
      a = a + 1#
    [base, _, _, _, _, _, _, _, _, _, _, _, _, _, _] = cem_plant(SMR, CC, chemical, Dry, echem, retro, burnH2, 
  Sell_SCM, Sell_Iron, Sell_Alumina, Sell_Aggregate, cleanH, 
  cleanE, constants[0], constants[1], constants[2], constants[3], 
  constants[4], constants[5], constants[6], constants[7], constants[8], 
  constants[9], constants[10], constants[11], constants[12], 
  constants[13], constants[14], constants[15], constants[16], 
  constants[17], constants[18], constants[19], constants[20], 
  constants[21], constants[22], constants[23], constants[24], 
  constants[25], constants[26], constants[27], constants[28], 
  constants[29], constants[30], constants[31], constants[32], 
  constants[33], constants[34], constants[35])#

#--- FINAL FIGURE ---
# Note: unable to fully test out this visualization, but it should work as is.
    fig8 = plt.figure(figsize = (6,5))
    rc('font', weight = 'bold')
    high_sorts = np.sort(high)
    high_I = np.argsort(high)

    low_sorts = low[high_I]
    names_sort = names[high_I]
    ys = range(len(high_sorts))[::-1]

    for y, low_sort, high_sort in zip(ys, low_sorts, high_sorts):
          low_width = base-low_sort
          high_width = low_sort + high_sort - base
    plt.broken_barh([(low_sort, low_width), (base, high_width)], (y - 0.4, 0.8), 
    facecolors = ['white', 'white'], edgecolors = ['black', 'black'], linewidth =1)

    x = base + high_width/2
    if x <= base + 50:
        x = base + high_width + 50
    plt.text(x, y, str(high_sort), va = 'center', ha = 'center')
    plt.axvline(base, color = 'black')

    axes = plt.gca()
    axes.spines['left'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    axes.xaxis.set_ticks_position('top')

    plt.yticks(ys, names_sort)

    plt.xlim(base-1000, base+1000)
    plt.ylim(-1, len(names_sort))

    plt.xlabel("$USD/T OPC")
    plt.legend()
    plt.show(fig8)
