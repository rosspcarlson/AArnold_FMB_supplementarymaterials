%FBAmIbiomO2  260102


%% type I methanotroph analysis

%load model
%model_biomass = readCbModel('Models/MethanotrophT1Biomass.mat');
model_biomass = readCbModel('Models/MethanoT1Biomass42AGAM.mat');

% checks on model performance
verifyModel(model_biomass)
% verifyModel: The following problems have been encountered in the model structure

rxnno = 124;

%output matrix declaration
mIbiomOUT = zeros(rxnno + 15,6);

i = 1;


while i<7
%scale fluxes, methane uptake
%model_biomass = changeRxnBounds(model_biomass,'T1',10,'b'); % Methane uptake rate
%scale fluxes, growth rate
%model_biomass = changeRxnBounds(model_biomass,'biomass_typeI',0.05 + 0.05*i,'b'); %growth rate gradient
model_biomass = changeRxnBounds(model_biomass,'biomass_typeI',0.3,'b'); %growth rate fixed

%scale fluxes, nitrate as nitrogen source/electron acceptor 
%model_biomass = changeRxnBounds(model_biomass,'T25',0,'b'); % prevent nitrate from serving as electron acceptor
%model_biomass = changeRxnBounds(model_biomass,'T10',0,'b'); % prevent nitrate from serving as nitrogen source

%scale fluxes, O2 uptake
model_biomass = changeRxnBounds(model_biomass,'T4',2-0.025*i,'u'); % O2 uptake gradient
 
%scale fluxes, nongrowth associated maintenance energy
model_biomass = changeRxnBounds(model_biomass,'ATP_main',0.42,'l'); 

%control methane oxidation pathway
  
% %model_biomass = changeRxnBounds(model_biomass,'R1',0,'b'); % pMMO methane uptake rate
% model_biomass = changeRxnBounds(model_biomass,'R2',0,'b'); % sMMO methane uptake rate
% model_biomass = changeRxnBounds(model_biomass,'direct_coupling',0,'b'); % direct couple methane uptake rate

model_biomass = changeRxnBounds(model_biomass,'R1',0,'b'); % pMMO methane uptake rate
%model_biomass = changeRxnBounds(model_biomass,'R2',0,'b'); % sMMO methane uptake rate
model_biomass = changeRxnBounds(model_biomass,'direct_coupling',0,'b'); % direct couple methane uptake rate

% model_biomass = changeRxnBounds(model_biomass,'R1',0,'b'); % pMMO methane uptake rate
% model_biomass = changeRxnBounds(model_biomass,'R2',0,'b'); % sMMO methane uptake rate
% %model_biomass = changeRxnBounds(model_biomass,'direct_coupling',0,'b'); % direct couple methane uptake rate
      
% calculate max biomass yield on methane
%model_biomass = changeObjective(model_biomass,'biomass_typeI'); %growth
model_biomass = changeObjective(model_biomass,'T1'); %methane uptake
find(model_biomass.c) %confirms the objective function was set
%FBAsolution_maxBiomass = optimizeCbModel(model_biomass,'max') %for biomass
FBAsolution_maxBiomass = optimizeCbModel(model_biomass,'min') % for methane

% grab the relevant reactions and calculate yields
Cmol_acetate = 2;
Cmol_acetone = 3;
Cmol_biomass = 1;
Cmol_CO2 = 1;
Cmol_ethanol = 2;
Cmol_formaldehyde = 1;
Cmol_formate = 1;
Cmol_lactate = 3;
Cmol_methane = 1;
Cmol_methanol = 1;
Cmol_PHB = 4;
Cmol_pyruvate = 3;
Cmol_succinate = 4;

acetate = find(strcmp('T15',model_biomass.rxns));
acetone = find(strcmp('T17',model_biomass.rxns));
biomass = find(strcmp('biomass_typeI',model_biomass.rxns));
CO2 = find(strcmp('T26',model_biomass.rxns));
ethanol = find(strcmp('T18',model_biomass.rxns));
formaldehyde = find(strcmp('T14',model_biomass.rxns));
formate = find(strcmp('T6',model_biomass.rxns));
lactate = find(strcmp('T16',model_biomass.rxns));
methane = find(strcmp('T1',model_biomass.rxns));
methanol = find(strcmp('T5',model_biomass.rxns));
PHB_produced = find(strcmp('T23',model_biomass.rxns));
pyruvate = find(strcmp('T19',model_biomass.rxns));
succinate = find(strcmp('T20',model_biomass.rxns));

O2 = find(strcmp('T4',model_biomass.rxns));
GWP_20 = find(strcmp('GWP20',model_biomass.rxns));
GWP_100 = find(strcmp('GWP100',model_biomass.rxns));

maxYield_biomass_methane = (FBAsolution_maxBiomass.v(biomass)*Cmol_biomass) / FBAsolution_maxBiomass.v(methane)*Cmol_methane;
O2toMethane = (FBAsolution_maxBiomass.v(O2)) / FBAsolution_maxBiomass.v(methane)*Cmol_methane;
CO2Yield = (FBAsolution_maxBiomass.v(CO2)*Cmol_CO2) / FBAsolution_maxBiomass.v(methane)*Cmol_methane;
GWP100Yield = FBAsolution_maxBiomass.v(GWP_100) / FBAsolution_maxBiomass.v(methane)*Cmol_methane;
acetateYield = (FBAsolution_maxBiomass.v(acetate)*Cmol_acetate)/ (FBAsolution_maxBiomass.v(methane)*Cmol_methane);
acetoneYield = (FBAsolution_maxBiomass.v(acetone)*Cmol_acetone)/ (FBAsolution_maxBiomass.v(methane)*Cmol_methane);
ethanolYield = (FBAsolution_maxBiomass.v(ethanol)*Cmol_ethanol)/ (FBAsolution_maxBiomass.v(methane)*Cmol_methane);
formaldehydeYield = (FBAsolution_maxBiomass.v(formaldehyde)*Cmol_formaldehyde)/ (FBAsolution_maxBiomass.v(methane)*Cmol_methane);
formateYield = (FBAsolution_maxBiomass.v(formate)*Cmol_formate)/ (FBAsolution_maxBiomass.v(methane)*Cmol_methane);
lactateYield = (FBAsolution_maxBiomass.v(lactate)*Cmol_lactate)/ (FBAsolution_maxBiomass.v(methane)*Cmol_methane);
methanolYield = (FBAsolution_maxBiomass.v(methanol)*Cmol_methanol)/ (FBAsolution_maxBiomass.v(methane)*Cmol_methane);
pyruvateYield = (FBAsolution_maxBiomass.v(pyruvate)*Cmol_pyruvate)/ (FBAsolution_maxBiomass.v(methane)*Cmol_methane);
succinateYield = (FBAsolution_maxBiomass.v(succinate)*Cmol_succinate)/ (FBAsolution_maxBiomass.v(methane)*Cmol_methane);
fermentsYield = (acetateYield + acetoneYield + ethanolYield + formaldehydeYield + formateYield + lactateYield + methanolYield + pyruvateYield + succinateYield);
PHBYield = (FBAsolution_maxBiomass.v(PHB_produced)*Cmol_PHB)/ (FBAsolution_maxBiomass.v(methane)*Cmol_methane);


FBAflux = FBAsolution_maxBiomass.x;
mIbiomOUT(1:rxnno,i) = FBAflux ; %fluxes output
mIbiomOUT(rxnno+1,i) = maxYield_biomass_methane;
mIbiomOUT(rxnno+2,i) = O2toMethane;
mIbiomOUT(rxnno+3,i) = CO2Yield;  
mIbiomOUT(rxnno+4,i) = GWP100Yield;  
mIbiomOUT(rxnno+5,i) = acetateYield; 
mIbiomOUT(rxnno+6,i) = acetoneYield; 
mIbiomOUT(rxnno+7,i) = ethanolYield;  
mIbiomOUT(rxnno+8,i) = formaldehydeYield;
mIbiomOUT(rxnno+9,i) = formateYield;  
mIbiomOUT(rxnno+10,i) = lactateYield;  
mIbiomOUT(rxnno+11,i) = methanolYield;  
mIbiomOUT(rxnno+12,i) = pyruvateYield;  
mIbiomOUT(rxnno+13,i) = succinateYield;
mIbiomOUT(rxnno+14,i) = fermentsYield;
%mIbiomOUT(rxnno+15,i) = PHBYield;  % Type I do not make PHB

i = i+1 ;
end
                      
xlswrite(["260113_mI_biom42GAM_NO3_sMMO.xlsx"],mIbiomOUT);   
