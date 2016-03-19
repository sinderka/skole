function [u,v,fx,fy]=kovasznay_lapl
lambda  = -1;
u    = @(x,y) 1-exp(lambda*(x-0.5)).*cos(2*pi*y);
v    = @(x,y) (lambda/2/pi)*exp(lambda*(x-0.5)).*sin(2*pi*y);
p    = @(x,y) 0.5*(exp(2*lambda*(x-0.5)));
d2u  = @(x,y) (-lambda^2 + 4*pi^2)*exp(lambda*(x-0.5)).*cos(2*pi*y);
d2v  = @(x,y) (lambda^3/2/pi - 2*pi*lambda)*exp(lambda*(x-0.5)).*sin(2*pi*y);
px   = @(x,y) lambda*exp(2*lambda*(x-0.5));
py   = @(x,y) zeros(size(x)).*zeros(size(y));
fx   = @(x,y) -d2u(x,y);% + px(x,y);
fy   = @(x,y) -d2v(x,y);% + py(x,y);
end
