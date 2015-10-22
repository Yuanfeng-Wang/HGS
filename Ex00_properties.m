%***********************************************************************************************************
%* HGS 1.3 
%* By Arnau Miro, Pau Manent, Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Example 00: Properties of single species and mixtures
clear; 

Ru=8.3144621/1000; % kJ/molK

% Using N2, verify that deltaH is aprox. equal to Cp*deltaT
[CpN2,~,~,~,~,~,hN2_1,~ ,~]=hgsprop({'N2'},[1],300,10);
[~   ,~,~,~,~,~,hN2_2,~ ,~]=hgsprop({'N2'},[1],301,10);

(hN2_2-hN2_1)
CpN2


% Using N2, verify that S2-S1 aprox. = cp.ln(T2/T1)-Ru*ln(P2/P1)

[CpN2,~,~,RN2,~,~,~,~,S1]=hgsprop({'N2'},[1],400,20);
[CpN2,~,~,RN2,~,~,~,~,S2]=hgsprop({'N2'},[1],500,10);
S2-S1

CpN2*log(500/400)-Ru*log(10/20)

% Using air (80% N2, 20% O2 in mol, or in volume) 
% verify @300K deltaH aprox. equal to Cp*deltaT
[CpA ,~,~,~,~,~,hA_1,~ ,~]=hgsprop({'N2' 'O2'},[0.8 0.2],300,10);
[~   ,~,~,~,~,~,hA_2,~ ,~]=hgsprop({'N2' 'O2'},[0.8 0.2],310,10);
(hA_2-hA_1)/10
CpA

% @300K and 1 bar verify:...

T=300 % K
P=1 % bar
m=10 % kg (random)

[Cp,Cv,MM  ,Rg ,gamma,a,H   ,G,S]=hgsprop({'N2' 'O2'},[0.8 0.2],T,P);

% ... that sound speed is correct
a
sqrt(gamma*Rg*1000*T) % (J/kgK)^(1/2)=m/s

% ... that Rg=Ru/MM 
1000*Ru/MM % kJ/kg K
Rg

% ... that Cp/Cv=gamma
Cp/Cv-gamma

% ... that g=h-T*s
v=1000*Rg*T/(P*1e5) % m^3/kg

% ... evaluate specific properties (per mass unit)
V=m*v
n=m*1000/MM % mol

h=H/m % kJ/kg
s=S/m % kJ/kgK
g=G/m % kJ/kg
u=h-P*1e5*v/1000

% ... check that g is correct
h-T*s
g

% verify molar cp of the mixture ..
CpO2=hgsprop({'O2'},[1],T,10) % kJ/molK
CpN2=hgsprop({'N2'},[1],T,10)

0.8*CpN2+0.2*CpO2-Cp % this has to be zero !

% compute cp and h per mass unit of N2

[~,~,MMN2,~,~,~,hN2,~,~]=hgsprop({'N2'},n,T,P);

CpN2kg=1000*CpN2/MMN2; % kJ/molK -> kJ/kgK
hN2kg=1000*hN2/MMN2; % kJ/mol -> kJ/kgK 

fprintf('N2: Cp=%f kJ/kgK hN2 =%f kJ/kg \n',CpN2kg, hN2kg );

TCO2=35+273  % K
pCO2=1.5     % bar

% Note that the program returns in units of kJ/molK. Here are converted to
% kJ/kgK
Cpkg=Cp/(MM/1000);

fprintf('cp=%f kJ/kgK \n',Cpkg);

% hgssingle computes H, G and S properties for a single element
hgssingle('O2','h',T,P)
