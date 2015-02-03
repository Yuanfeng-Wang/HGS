%***********************************************************************************************************
%* HGS 1.3
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Example 05b: Adiabatic H2 / O2 reaction using hgsTp
%
% Inlet: H2, O2 
% Outlet: H2O + (1/2)O2 at Tp

clear

species={'H2','O2','H2O','H','O','OH'};
Tr=350; % K 
P=10; % bar
nr=[2;1;0;0;0;0]; % mol

[Tp,np]=hgsTp(species,nr,298,1)
np/sum(np)