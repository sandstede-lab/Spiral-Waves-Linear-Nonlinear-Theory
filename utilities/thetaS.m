function [val,val_prime] = thetaS(x,s)
% Returns smoothed Heaviside function in same dimension as input x
% x: input vector
% s: smoothness parameter. Higher s corresponds to steaper transition

val = 0.5.*(1 + tanh(s.*x));
val_prime = 0.5 .* s .* (sech(s.*x).^2);



