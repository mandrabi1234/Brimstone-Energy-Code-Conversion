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
#print("--6IX9INE--")
#print(sens_var)
#print(SENS)
#SENS[3] = sens_var[3]
#print(len(SENS[3]))
#sens_var = np.asarray(sens_var)
#print(constants)

# find the index of the variables to compare in the sens_var cell array
# last 2 values are the constants that would be used if the variables
# selected were not variables, and therefore can be ignored


#sens_var = np.array(sens_var)
#sens_var = sens_var.astype(int)
print(sens_var)
print(sens_var[6], len(SENS[6]))
if tornado == 0:
  for i in range(len(sens_var) - 2):
    print(sens_var[i])

    if ((indexx == -1) and (indexy == -1) and ((len(SENS[i])) > 1)):
      x = sens_var[i] # defines the variable
      indexx = 1      
      sens_var[i] = -1 # replaces variable with -1            
    elif ((indexy == -1) and ((len(SENS[i])) > 1)):
      y = sens_var[i]
      print("--Y--", y)
      indexy = i
      sens_var[i] = -1              
      break

  print('X',x)
  print('Y', y)
  length_x = len(x)
  length_y = len(y)
  
  #print("Length Check")
  #print(length_x)
  #print(length_y)


  z = np.zeros((length_y, length_x))
  #print("Z-length", z)
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
  a = 1 # indexes
  for i  in x: #THIS MIGHT BE WRONG (CHECK WITH A TEST RUN)
    b = 1 # % indexes
    for j in (y): #THIS MIGHT BE WRONG (CHECK WITH A TEST RUN)
      inputs = (sens_var) # inputs = cell2mat(sens_var) # 
      print("SENS_VAR")
      print(sens_var)
      print("")
      print("INPUTS")
      print(len(inputs))
      print(inputs)
      print(inputs[0])
      print(inputs[3])
      print("")
      inputs[indexx] = i # % update the proper input with the next value
      inputs[indexy]= j # % update the proper input with the next value

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
  print("CONSTANTS", constants)

  inputs[indexx] = sens_var[len(sens_var) - 2] # update the proper input with the next value
        
  inputs[indexy] = sens_var[len(sens_var)-1]
  print(inputs[indexy])
  print(inputs[indexx]) #
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
  if log_plot == 1 and (sum(sum(abs(z))) <= sum(sum(z))):
    z = math.log(z)#
    zlabel = 'log cost per kg of hydrogen'#
  
#------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------UNFINISHED CODE-----------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------

#cp = plt.contourf(x, y, z, 100)
#plt.colorbar(cp)
#plt.plot(cp)

'''
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(6,5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height]) 

start, stop, n_values = -8, 8, 800

x_vals = np.linspace(start, stop, n_values)
y_vals = np.linspace(start, stop, n_values)
X, Y = np.meshgrid(x_vals, y_vals)


Z = np.sqrt(X**2 + Y**2)

cp = plt.contourf(X, Y, Z)
plt.colorbar(cp)

ax.set_title('Contour Plot')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
plt.show()
'''
fig1 = plt.figure(figsize = (6,5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig1.add_axes([left, bottom, width, height])

cp = plt.contourf(x, y, z)
plt.colorbar(cp)

ax.set_title('Contour Plot')
ax.set_xlabel('x (cm')
ax.set_ylabel('y (cm)')
plt.show(fig1)
'''
        figure(1)
        % filled
        contour
        plot
        of
        co - variables
        to
        compare
        the
        costs
        of
        H2
        contourf(x, y, z)
        set(gca, 'Fontsize', 22)
        hold
        on
        c = colorbar()#
        labels = ...
        {'Electrolyzer Cost (USD/kW)', 'CO_2 Tax (USD/T CO_2)'...
        'H_2 price (USD/kg)', 'Mining Cost (USD/T)', ...
        'Cost of Heat (USD/GJ)', 'SCM Sold Per T OPC (T SCM/T OPC)', ...
        'Electricity Cost ($/kWhr)', 'kWhr per kg H_2', 'Discount Rate', ...
        ' Fraction CaO in Source Rock', ' Fraction Al_2O_3 in Source Rock', ...
        ' Fraction SiO_2 in Source Rock', ' Fraction FeO in Source Rock', ...
        'Freaction Fe_2O_3 in Source Rock', 'Fraction Mg_O, in Source Rock', ...
        'CO_2 Intensity of Heat (T CO_2/GJ)', ...
        'CO_2 Intensity of Electricity (T CO_2/kWhr)', ...
        'Price per T of SCM (USD/T)', 'Sulfur Cost (USD/T)', ...
        'Capacity Factor', 'Plant Size (T OPC/yr)', ...
        'Conventration of Produced Sulfuric Acid', ...
        'Heat to Power Efficiency', 'Co-product Revenue per T OPC'...
        'Price per Tonne Fe_2O_3', 'Price per Tonne Al_2_O_3', ...
        'Price per Tonne Aggregate', 'Operating Cell Voltage (V)', ...
        'Water Content (T/T OPC)', ...
        'Improvement in Current Density (A/cm^2)', ...
        'Al production efficiency', 'Fe production efficiency', ...
        'SCM production efficiency', 'OPC production efficiency', ...
        'Aggregate production efficiency'}#
        xlabel(labels
        {indexx})#
        ylabel(labels
        {indexy})#
        c.Label.String = zlabel#
'''
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print("             Data Visualization Check             ")
print("Dry_base", Dry_base)
print(CapExPT)
#print(sum(CapExPT))
print(((CapExPT) + OpExPT - Dry_base))
print(cost)
print(SCM_value)
print(Iron)
print(Al)
print(Agg)
print(H2_value)
#print(sum(CapExPT))
print(OpExPT)
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

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
'''
        % % % Value
        of
        All
        Products
        figure(2)
        key = {'Dry USD/T OPC', 'Break-Even Co-Product Value', 'USD/T OPC', ...
              'Value SCM', 'Value Iron', 'Value Aluminum', 'Value Aggrigate', ...
              'Value H_2', 'CapEx per T', 'OpEx pet T'}#
        Y = [Dry_base, (sum(CapExPT) + OpExPT - Dry_base), cost, SCM_value, Iron, ...
            Al, Agg, H2_value, sum(CapExPT), OpExPT]#
        bar(Y)#
        text(1: length(Y), Y, num2str(Y
        '),'
        vert
        ','
        bottom
        ','
        horiz
        ','
        center
        ')#
        box
        off
        hold
        on
        ylabel({'$/T OPC'})
        set(gca, 'XTickLabel', key, 'Fontsize', 22, 'XTickLabelRotation', 45)#
'''
'''
# y-axis bolded
rc ('font', weight = 'bold')

# Values of each group


plt.xticks()
plt.ylabel("CapEx, 1MT OPC Plant (USD)")
#plt.xlabel("")
plt.legend( loc = "upper right")
plt.show()

        % % % % % heat, add
        equivlent
        heat
        to
        graph!!!!!!!!!!!!!

        figure(3)
        key = {'Dry Process CapEx', 'Brimstone CapEx'}#
        stack = {'Cement + SCM Plant', 'Hydrogen', 'Acid Regeneration', 'Acid Reactor'}#
        DryCapEx = [DryCapExMat(1), 0, 0, 0]#
        BrimCapEx = [CapExMat(1) *not (Dry), CapExMat(2) * echem *not (Dry), ...
                    CapExMat(3) *not (Dry), CapExMat(4) *not (Dry)]#
        % ECapEx = [CapExMat(1) + CapExMat(6), CapExMat(5) * ((SMR & & Dry) | | (echem & &not (Dry))), 0, 0]#
        h = bar([DryCapEx#
        BrimCapEx], 'stacked')#
        hold
        on
        ylabel({'CapEx, 1MT OPC Plant (USD)'})
        set(gca, 'XTickLabel', key, 'Fontsize', 22, 'XTickLabelRotation', 45)#
        legend(h, stack)#

        figure(4)
        key = {'Dry Process Heat', ...
              'Equivalent OPC + SCM + Fe_2O_3 + Al_2_O_3 + H_2 Heat', ...
              'Brimstone Process Heat', 'H_2 Combustion Heat'}#
        stack = {'Net Reacton Heat', 'Net Sensible Heat', 'Net Latent Heat'}#
        Q = DryQDry#
        h = bar([Q#
        QDry#
        QBrim *
        not (Dry)#
        QH2], 'stacked')#
        hold
        on
        ylabel({'Heat Energy (GJ/T OPC)'})
        set(gca, 'XTickLabel', key, 'Fontsize', 22, 'XTickLabelRotation', 45)#
        legend(h, stack)#

        figure(5)
        key = {'Dry Process', 'Brimstone Total', 'Equivalent OPC + SCM + Fe_2O_3 + Al_2_O_3 + H_2'}#
        stack = {'Process', 'Heat', 'Electricity'}#
        DryCapEx = [DryGHG(1), DryGHG(4), DryGHG(2)]#
        BrimCapEx = [0, GHG(5) *not (cleanH), GHG(3) *not (cleanE)]#
        ECapEx = [GHG(1), (GHG(4) + GHG_Div(6) + GHG_Div(8)) *not (cleanH), GHG(2) *not (cleanE)]#
        h = bar([DryCapEx#
        BrimCapEx#
        ECapEx], 'stacked')#
        hold
        on
        ylabel({'CO_2 Intensity (T CO_2/TOPC)'})
        set(gca, 'XTickLabel', key, 'Fontsize', 22, 'XTickLabelRotation', 45)#
        legend(h, stack)#

        figure(6)
        key = {'Dry Process OpEx', 'Brimstone OpEx'}#
        stack = {'O&M', 'Mining', 'Heat', 'Electricity', 'CO_2 Tax', 'Water', 'Sulfur'}#
        DryCapEx = DryOpExMat#
        BrimCapEx = [OpExMat(1) *not (Dry), OpExMat(2) *not (Dry), ...
                    OpExMat(3) *not (Dry), OpExMat(4) *not (Dry), OpExMat(5) *not (Dry), ...
                    OpExMat(6) *not (Dry), OpExMat(7) *not (Dry), OpExMat(8) *not (Dry), OpExMat(9) *not (Dry)]#
        h = bar([DryCapEx#
        BrimCapEx], 'stacked')#
        hold
        on
        ylabel({'OpEx/T OPC, 1MTPY OPC Plant (USD)'})
        set(gca, 'XTickLabel', key, 'Fontsize', 22, 'XTickLabelRotation', 45)#
        legend(h, stack)#

        % % % % not correct
        for H2 case

fig2 = plt.figure(figsize = (6,5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig2.add_axes([left, bottom, width, height])

data = {"Conv. Process Cement":Dry_GHG, "Brimstone Cement":(GHG_Div(0)), 
"Conv. Alumina":GHG_Div(5), "Brimstone Alimina":GHG_Div(1), "Conv. Iron Oxide":GHG_Div(6), "Brimstone Iron Oxide":GHG_Div(2),
"Conv. SCM":(np.array([GHG(0), GHG(4) * (not cleanH), GHG(1) * (not cleanE)]), "Value H_2":H2_value, "CapEx per T":(CapExPT), "OpEx pet T":OpExPT} 

key = list(data.keys())
Y = list(data.values())
bar = plt.bar(key, Y, width = 0.5)
plt.xticks(rotation = 45)
ylabel = ("$/T OPC")
for index, value in enumerate(Y):
    plt.text(index -0.5, value, str(value))
plt.show(fig2)

        figure(7)
        key = {'Conv. Process Cement', 'Brimstone Cement', 'Conv. Alumina', ...
        'Brimstone Alimina', 'Conv. Iron Oxide', 'Brimstone Iron Oxide', ...
        'Conv. SCM', 'Brimstone SCM', 'Conv. H_2', 'Brimstone H_2'}#
        hb = bar([sum(DryGHG), GHG_Div(1), ...
        GHG_Div(6), GHG_Div(2), GHG_Div(7), GHG_Div(3), (sum([GHG(1), ...
        GHG(4) * not (cleanH), GHG(2) * not (cleanE)])- sum(DryGHG)) * Sell_SCM, ...
        GHG_Div(4), GHG_Div(8), GHG_Div(5)])#
        hold on
        ylabel({'TCO_2/T OPC produced'})
        set(gca, 'XTickLabel', key, 'Fontsize', 22, 'XTickLabelRotation', 45)#





else:
  t_data = cell(1, len(constants))#
  for i in range(len(constants):
    costs = [0, 0]#
    inputs = constants#
    inputs[i] = minimum(sens_var[i])
  [costs[1], _, _, _, _, _, _, _, _, _, _, _, _, _, _] = 
  cem_plant(SMR, CC, chemical, Dry, echem, retro, burnH2, 
  Sell_SCM, Sell_Iron, Sell_Alumina, Sell_Aggregate, cleanH, 
  cleanE, inputs[0], inputs[1], inputs[2], inputs[3], inputs[4], inputs[5], 
  inputs[6], inputs[7], inputs[8], inputs[9], inputs[10], 
  inputs[11], inputs[12], inputs[13], inputs[14], inputs[15], 
  inputs[16], inputs[17], inputs[18], inputs[19], 
  inputs[20], inputs[21], inputs[22], inputs[23], inputs[24], 
  inputs[25], inputs[26], inputs[27], inputs[28], inputs[29], 
  inputs[30], inputs[31], inputs[32], inputs[33], inputs[34], 
  inputs[35])

  inputs[i] = maximum(sens_var[i])#
  [costs[2], _, _, _, _, _, _, _, _, _, _, _, _, _, _] = 
  cem_plant(SMR, CC, chemical, Dry, echem, retro, burnH2, 
  Sell_SCM, Sell_Iron, Sell_Alumina, Sell_Aggregate, cleanH, 
  cleanE, inputs[0], inputs[1], inputs[2], inputs[3], inputs[4], inputs[5], 
  inputs[6], inputs[7], inputs[8], inputs[9], inputs[10], 
  inputs[11], inputs[12], inputs[13], inputs[14], inputs[15], 
  inputs[16], inputs[17], inputs[18], inputs[19], inputs[20], 
  inputs[21], inputs[22], inputs[23], inputs[24], 
  inputs[25], inputs[26], inputs[27], inputs[28], inputs[29], 
  inputs[30], inputs[31], inputs[32], inputs[33], inputs[34], 
  inputs[35])#
  t_data[i] = [costs]#


 # Names of the Y-axis ticks 
names = np.array['Electrolyzer Cost (USD/kW) 50->1000', 'CO_2 Tax (USD/T CO_2) -190->100', 
'H_2 price (USD/kg)1->4', 'Mining Cost (USD/T) 3->10', 'Cost of Heat (USD/GJ) 1->10', 
'SCM Sold Per T OPC (T SCM/T OPC) 0.4->5.3', 'Electricity Cost ($/kWhr) 0.01->0.1', 
'kWhr per kg H_2 5->30', 'IRR 4%->30%', 'CO_2 Intensity of Heat (T CO_2/GJ) 0->0.2', 
'CO_2 Intensity of Electricity (T CO_2/kWhr) 0->0.001', 'Price per T of SCM (USD/T) 0->100', 
'Sulfur Cost (USD/T) 0->200', 'Plant Size (T OPC/yr) 1M->5M', 'Co-product Revenue per T OPC 0->500', 
'Price per Tonne Fe_2O_3 0->100', 'Price per Tonne Al_2O_3 0->400', 'Price per Tonne Aggregate 0->20', 
'Operating Cell Voltage (V) 0.2->1.2', 'Water Content (T/T OPC) 0.1-5', 'Improvement in Current Density (A/cm^2) 0.4-> -0.4', 
'Al production efficiency 50%->100%', 'Fe production efficiency 50%->100%', 'SCM production efficiency 50%->100%', 
'OPC production efficiency 50%->90%', 'Aggregate production efficiency 50%->100%', 'CapEx multiplier 0.5->2']

# ' Fraction CaO in Source Rock', ...
# ' Fraction Al_2O_3 in Source Rock', ' Fraction SiO_2 in Source Rock', ...
# ' Fraction FeO in Source Rock', 'Freaction Fe_2O_3 in Source Rock', ...
# 'Fraction MgO, in Source Rock',
# 'Capacity Factor',
# 'Conventration of Produced Sulfuric Acid', ...
# 'Heat to Power Efficiency',

vars = [1:9, 16: 19, 21, 24: len(t_data)]#
a = 1#
high = np.zeros(shape(vars))#
low = np.zeros(shape(vars))#
for i in vars: 
  n = t_data[i]#
  high(a) = maximum(n)#
  low(a) = minimum(n)#
  a = a + 1#
[base, _, _, _, _, _, _, _, _, _, _, _, _, _, _] = 
cem_plant(SMR, CC, chemical, Dry, echem, retro, burnH2, 
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

  figure

  [high_sort, high_I] = sort(high, 'ascend')#
  low_sort = low(high_I)#
  names_sort = names(high_I)#
  h = barh(high_sort)#
  hold
  on
  barh(low_sort, 'r')
  bh = get(h, 'BaseLine')#
  set(bh, 'BaseValue', base)#
  set(gca, 'yticklabel', names_sort)
  set(gca, 'Ytick', [1: length(names)], 'YTickLabel', [1: length(names)])
  set(gca, 'yticklabel', names_sort)
  xlabel('$USD/T OPC')
  hold
  on
  end
'''