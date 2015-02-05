function [ xc,yc ] = hgssolve(fun,x0,solver,options)
%***********************************************************************************************************
%* HGS 1.3
%* By Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
% 
% Solver: solver selector for hgs code.
% For any issues with the code see the documentation manual.
%
%
% Inputs:
%
% Output:
%
% See also FSOLVE, FZERO
%
%   This code is part of the HGS TOOLBOX
%   OpenLLOP, UPC-ETSEIAT 2014-2015

% initial parameters if option is emtpy
if isempty(options) && strcmp(solver,'hgsfzero')
    options = struct('x2',5000,'fchange',2,'epsx',1e-1,'epsy',1e-4,'maxite',200,'info',[]);
elseif  isempty(options) % assuming fsolve or fzero
    options = optimset('Display','none');
end

% Get solver as function handle
solv = str2func(solver);

% Launch solver
[xc,yc,flag] = solv(fun,x0,options);

% If flag is negative then run error check code
if flag < 0, ErrorControl(fun,xc,x0,10), end

end

function ErrorControl(f,xc,x0,n)

% Get interval points
x1 = xc -x0/4;
x2 = xc + x0/4;

% Print Error message
fprintf('Problems... plotting function to be solved \n');

% Compute function values where the error has happened
xv=linspace(x1,x2,n); yv(1:n,1) = 0;
for ii=1:length(xv)
    yv(ii)=f(xv(ii));
end

% Plot function to be solved
 figure('Name','hgssolver error control',...
             'Color','w',...
             'NumberTitle','off');
xlabel('x values','Fontsize',16);
ylabel('function values','Fontsize',16);
grid on;

plot(xv,yv,'Linewidth',1.5);

% Raise error
error('uhhh ? no sign change x1=%e y1=%e x2=%e y2=%e',xv(1),yv(1),xv(n),yv(n));

end