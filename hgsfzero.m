function [ xc,yc ] = hgsfzero( f,x1,x2, fchange,epsx,epsy,maxite,info )
%***********************************************************************************************************
%* HGS 1.3
%* By Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
% 
% Solver: optimized solver for hgs code.
% For any issues with the code see the documentation manual.
%
% Usage:
%       [ xc,yc ] = HGSFZERO( f,x1,x2, fchange,epsx,epsy,maxite,info )
%
% Inputs:
%
% Output:
%
% See also FSOLVE, FZERO
%
%   This code is part of the HGS TOOLBOX
%   OpenLLOP, UPC-ETSEIAT 2014-2015

y1=f(x1);
y2=f(x2);
if y1*y2 > 0 
    fprintf('myfzero.. problems... plotting function to be solved \n');
    xxv=linspace(x1,x2,10);
    for ii=1:length(xxv)
        yyv(ii)=f(xxv(ii));
    end
    figure
    plot(xxv,yyv);
    error('uhhh ? no sign change x1=%e y1=%e x2=%e y2=%e',x1,y1,x2,y2);
end

for i=1:maxite
    if x2-x1 > fchange % fchange should be about 1
        xc=(x1+x2)/2;
        fc=1;
    else
        xc=x1-y1*(x2-x1)/(y2-y1);
        fc=0;
    end
    
    yc=f(xc);
    if ~isempty(info)
        fprintf('myfzero i=%d fc=%d x1=%7.2e y1=%7.2e x2=%7.2e y2=%7.2e xc=%7.2e dx=%7.2e |yc|=%7.2e \n',...
            i,fc,x1,y1,x2,y2,xc,x2-x1,abs(yc) );
    end
    if abs(yc)<epsy && x2-x1<epsx
        break;
    end
    if yc*y1>0
        y1=yc;
        x1=xc;
        flag=1;
    else
        y2=yc;
        x2=xc;
        flag=2;
    end

end

end

