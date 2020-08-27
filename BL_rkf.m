function [tp,xp,yp] = BL_rkf(dxdt,tspan,x0,yFcn,y0,varargin)

if nargin<5,error('at least 5 input arguments required'), end
if any(diff(tspan)<=0),error('tspan not ascending order'), end
if ~isempty(varargin)
    if any(strcmp(varargin,'RelTol')) == 1, k=find(strcmp(varargin,'RelTol')); RelTol = varargin{k+1}; end
    if any(strcmp(varargin,'hlim')) == 1, k=find(strcmp(varargin,'hlim')); hlim = varargin{k+1}; end
end        

if ~exist('hlim','var'), hlim = 1e-4; end % Standard maximum timestep
if ~exist('RelTol','var'), RelTol = 1e-3; end % Standard relative tolerance

n = length(tspan);
ti = tspan(1);tf = tspan(n);
tc = ti;          % Current time
x(1,:) = x0;
tp(1,1) = tc;       
xp(1,:) = x(1,:); 
yp(1,:) = y0; 
i = 1;
global tv0
tv0 = -1;
while tp(i) < tf
    hc = hlim;          % Reset current timestep to maximum allowable
    eps = 10*RelTol;    % Reset tolerance
    while eps>RelTol
        if tc+hc > tf, hc = tf-tc; end
        k1 = dxdt(tc,x(i,:));
        x1 = x(i,:) + k1*hc/4;
        k2 = dxdt(tc+hc/4,x1);
        x2 = x(i,:) + (3*k1/32+9/32*k2)*hc;
        k3 = dxdt(tc+3*hc/8,x2);
        x3 = x(i,:) + (1932*k1-7200*k2+7296*k3)/2197*hc;
        k4 = dxdt(tc+12/13*hc,x3);
        x4 = x(i,:) + (439/216*k1-8*k2+3680/513*k3-845/4104*k4)*hc;
        k5 = dxdt(tc+hc,x4);
        x5 = x(i,:) + (-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5)*hc;
        k6 = dxdt(tc+hc/2,x5);
        x6 = x(i,:) + (16/135*k1+6656/12825*k3+28561/56430*k4-9/50*k5+2/55*k6)*hc;
        eps = sqrt((1/360*k1-128/4275*k3-2197/75240*k4+1/50*k5+2/55*k6)*(1/360*k1-128/4275*k3-2197/75240*k4+1/50*k5+2/55*k6)')*hc; % Error between 5th and 4th order approximations
        hc = hc/2; % Step halving
    end
    hc = 2*hc; % Correcting
    x(i+1,:) = x6;
    tc = tc+hc;
    tp(i+1,1) = tc;       
    xp(i+1,:) = x(i+1,:);
    yp(i+1,:) = yFcn(tc,x(i+1,:));
    i=i+1;
end