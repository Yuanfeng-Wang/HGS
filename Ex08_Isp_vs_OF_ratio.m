%***********************************************************************************************************
%* HGS 1.3 
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Example 08: ISP vs OF ratio

function Ex08_Isp_vs_OF_ratio

clear;


species={'H','H2','H2O','H2O2','HO2','O','O2','OH'};

% Inlet temperature as if the reactives were gas at 300K
Te=300 % K   reactives inlet temperature
Pc=50  % bar chamber pressure 
P2=0.1 % bar nozzle exit    
    
vrof=[]; % vector of OF ratios
visp=[]; % vector of specific impulses



for rof=2:0.25:8
    
    % evaluate mol of each specie at inlet for a given ROF ratio
    nO2=1;
    mO2=nO2*32;
    mH2=mO2/rof;
    nH2=mH2/2;
    
    ni_i=[0;... % H
    nH2;... % H2
    0;...   % H2O
    0;...   % H2O2
    0;...   % HO2
    0;...   % O
    nO2;... % O2
    0];     % OH

    ni_i=ni_i/sum(ni_i); % mole fractions  
    

    % Evaluate inlet properties with HGS assuming gas state
    [~,~,MM,~,~,~,H,~,~]=hgsprop(species,ni_i,Te,Pc);
    n=sum(ni_i); % mixture total number of mols (1)
    m=n*MM*1e-3 % mixture mass kg
    h1G=H/m;     % inlet mixture enthalpy in GAS state kJ/kgK 
                 % we evaluate it just for comparision

    % Inlet enthalpy as if reactives were satured liquid at 10 bar
    % O2 (NIST) hv(404.36 K)-hl(119.62 K)=14.3753 kJ/mol
    % H2 (NIST) hv(413.96K)-hl(31.39K)=10.9495 kJ/mol

    % Enthalpy of O2 liq at Tsat 10 bar (kJ/mol)
    hO2=hgssingle('O2','h',404.36,10)-14.3753; 
    % Enthalpy of H2 liq at Tsat 10 bar (kJ/mol)
    hH2=hgssingle('H2','h',413.96,10)-10.9495; 
    Hin=ni_i(2)*hH2+ni_i(7)*hO2;
    h1=Hin/m % inlet mixture enthalpy in LIQUID state kJ/kg 

    % We find temperature at nozzle inlet solving for Delta_H=0
    % hgsTp function can't be used as it assumes gas state
    
    Tc=fzero(@DeltaH,3000,optimset('Display','iter'));
    %Tc=myfzero(@DeltaH,300,4500,10,1e-4,1e-5,300);
    ni_calc=hgseq(species,ni_i,Tc,Pc);
    [Cp,Cv,MM,Rg,gamma,a,H,G,S]=hgsprop(species,ni_calc,Tc,Pc);
    m=sum(ni_calc)*MM*1e-3 % mixture mass kg (has to be as before!)
    s=S/m;    
        
    
    [ T2,n2 ] = hgsisentropic(species,ni_calc,Tc,Pc,P2,'shifting'  );
    
    [~,~,MM2,~,~,a2,H2,~,S2]=hgsprop(species,n2,T2,P2);
    m2=sum(n2)*MM2*1e-3; % mixture mass kg
    h2=H2/m2;
    
    %{ 
    %check..
    s
    s2=S2/m2 % kJ/kgK
    
    %}
    
    % here h1 is the chamber inlet enthalpy, that is equal to the
    % nozzle inlet enthalpy
    vt=sqrt(2*1000*(h1-h2)); % Enthalpy en J !

    Is=vt/9.81; % Is (optimal expansion, Pe=Pambient)    

    vrof(end+1)=rof;
    visp(end+1)=Is;
    
    plot(vrof,visp); xlabel('ratio O/F'); ylabel('Isp (s)');
    drawnow;
end



    function DeltaH=DeltaH(Tc)
        nc=hgseq(species,ni_i,Tc,Pc);
        [~,~,MMC,~,~,~,HC,~,~]=hgsprop(species,nc,Tc,Pc); 
        nC=sum(nc); % mixture total number of mols (1)
        mc=nC*MMC*1e-3; % mixture mass kg
        hc=HC/mc; % kJ/kgK
        DeltaH=hc-h1;
    end
    
end

