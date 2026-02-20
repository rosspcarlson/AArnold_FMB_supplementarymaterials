%FBAmIIATPO2  260102


%% type II methanotroph analysis

% ATP yield on methane
%     
% load model. This one is for calculating ATP yields on methane, when
% 
model_ATP = readCbModel('Models/MethanotrophT2ATP.mat');

% checks on model performance
verifyModel(model_ATP)
% verifyModel: The following problems have been encountered in the model structure

% number of reactions in model for output definition
rxnno = 130;

%output
mIIATPOUT = zeros(rxnno + 15,10);

i = 1;

%loop for generating O2 gradient
while i<11 
%scale fluxes, methane uptake
model_ATP = changeRxnBounds(model_ATP,'T1',10,'b'); % Methane uptake rate

%scale fluxes, nitrate as nitrogen source/electron acceptor  
model_ATP = changeRxnBounds(model_ATP,'T25',0,'b'); % prevent nitrate from serving as electron acceptor
%scale fluxes, N2 as nitrogen source 
%model_ATP = changeRxnBounds(model_ATP,'T11',0,'b'); % prevent N2 from serving as nitrogen source

%scale fluxes, O2 uptake
model_ATP = changeRxnBounds(model_ATP,'T4',21-i ,'u'); % O2 uptake rate

%control methane oxidation pathway

%pMMO simulations
% %model_ATP = changeRxnBounds(model_ATP,'R1',0,'b'); % pMMO methane uptake rate
% model_ATP = changeRxnBounds(model_ATP,'R2',0,'b'); % sMMO methane uptake rate
% model_ATP = changeRxnBounds(model_ATP,'direct_coupling',0,'b'); % direct couple methane uptake rate

%sMMO simulations
% model_ATP = changeRxnBounds(model_ATP,'R1',0,'b'); % pMMO methane uptake rate
% %model_ATP = changeRxnBounds(model_ATP,'R2',0,'b'); % sMMO methane uptake rate
% model_ATP = changeRxnBounds(model_ATP,'direct_coupling',0,'b'); % direct couple methane uptake rate

%dcMMO simulations
model_ATP = changeRxnBounds(model_ATP,'R1',0,'b'); % pMMO methane uptake rate
model_ATP = changeRxnBounds(model_ATP,'R2',0,'b'); % sMMO methane uptake rate
%model_ATP = changeRxnBounds(model_biomass,'direct_coupling',0,'b'); % direct couple methane uptake rate

% calculate max ATP yield on methane
model_ATP = changeObjective(model_ATP,'ATP_main');
find(model_ATP.c) %confirms the objective function was set
FBAsolution_maxATP = optimizeCbModel(model_ATP,'max')

% relevant reactions for calculating yields
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

acetate = find(strcmp('T15',model_ATP.rxns));
acetone = find(strcmp('T17',model_ATP.rxns));
ATPmain = find(strcmp('ATP_main',model_ATP.rxns));
biomass = find(strcmp('biomass_typeI',model_ATP.rxns));
CO2 = find(strcmp('T26',model_ATP.rxns));
ethanol = find(strcmp('T18',model_ATP.rxns));
formaldehyde = find(strcmp('T14',model_ATP.rxns));
formate = find(strcmp('T6',model_ATP.rxns));
lactate = find(strcmp('T16',model_ATP.rxns));
methane = find(strcmp('T1',model_ATP.rxns));
methanol = find(strcmp('T5',model_ATP.rxns));
PHB_produced = find(strcmp('T23',model_ATP.rxns));
pyruvate = find(strcmp('T19',model_ATP.rxns));
succinate = find(strcmp('T20',model_ATP.rxns));

O2 = find(strcmp('T4',model_ATP.rxns));
GWP_20 = find(strcmp('GWP20',model_ATP.rxns));
GWP_100 = find(strcmp('GWP100',model_ATP.rxns));

maxYield_ATP_methane = (FBAsolution_maxATP.v(ATPmain) / FBAsolution_maxATP.v(methane)*Cmol_methane);
O2toMethane = (FBAsolution_maxATP.v(O2)) / FBAsolution_maxATP.v(methane)*Cmol_methane;
CO2Yield = (FBAsolution_maxATP.v(CO2)*Cmol_CO2) / FBAsolution_maxATP.v(methane)*Cmol_methane;
GWP100Yield = FBAsolution_maxATP.v(GWP_100) / FBAsolution_maxATP.v(methane)*Cmol_methane;
acetateYield = (FBAsolution_maxATP.v(acetate)*Cmol_acetate)/ (FBAsolution_maxATP.v(methane)*Cmol_methane);
acetoneYield = (FBAsolution_maxATP.v(acetone)*Cmol_acetone)/ (FBAsolution_maxATP.v(methane)*Cmol_methane);
ethanolYield = (FBAsolution_maxATP.v(ethanol)*Cmol_ethanol)/ (FBAsolution_maxATP.v(methane)*Cmol_methane);
formaldehydeYield = (FBAsolution_maxATP.v(formaldehyde)*Cmol_formaldehyde)/ (FBAsolution_maxATP.v(methane)*Cmol_methane);
formateYield = (FBAsolution_maxATP.v(formate)*Cmol_formate)/ (FBAsolution_maxATP.v(methane)*Cmol_methane);
lactateYield = (FBAsolution_maxATP.v(lactate)*Cmol_lactate)/ (FBAsolution_maxATP.v(methane)*Cmol_methane);
methanolYield = (FBAsolution_maxATP.v(methanol)*Cmol_methanol)/ (FBAsolution_maxATP.v(methane)*Cmol_methane);
pyruvateYield = (FBAsolution_maxATP.v(pyruvate)*Cmol_pyruvate)/ (FBAsolution_maxATP.v(methane)*Cmol_methane);
succinateYield = (FBAsolution_maxATP.v(succinate)*Cmol_succinate)/ (FBAsolution_maxATP.v(methane)*Cmol_methane);
fermentsYield = (acetateYield + acetoneYield + ethanolYield + formaldehydeYield + formateYield + lactateYield + methanolYield + pyruvateYield + succinateYield);
PHBYield = (FBAsolution_maxATP.v(PHB_produced)*Cmol_PHB)/ (FBAsolution_maxATP.v(methane)*Cmol_methane);

FBAflux = FBAsolution_maxATP.x;
mIIATPOUT(1:rxnno,i) = FBAflux ; %fluxes output
mIIATPOUT(rxnno+1,i) = maxYield_ATP_methane;
mIIATPOUT(rxnno+2,i) = O2toMethane;
mIIATPOUT(rxnno+3,i) = CO2Yield;  
mIIATPOUT(rxnno+4,i) = GWP100Yield;  
mIIATPOUT(rxnno+5,i) = acetateYield; 
mIIATPOUT(rxnno+6,i) = acetoneYield; 
mIIATPOUT(rxnno+7,i) = ethanolYield;  
mIIATPOUT(rxnno+8,i) = formaldehydeYield;
mIIATPOUT(rxnno+9,i) = formateYield;  
mIIATPOUT(rxnno+10,i) = lactateYield;  
mIIATPOUT(rxnno+11,i) = methanolYield;  
mIIATPOUT(rxnno+12,i) = pyruvateYield;  
mIIATPOUT(rxnno+13,i) = succinateYield;
mIIATPOUT(rxnno+14,i) = fermentsYield;
mIIATPOUT(rxnno+15,i) = PHBYield;   

i = i+1 ;
end

%name of Excel output file
                       
xlswrite(["260113_mII_ATP_dcMMO_N2.xlsx"],mIIATPOUT);
