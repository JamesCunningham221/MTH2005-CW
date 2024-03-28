function [ es ] = svp(T)

%Calculation of Saturation Vapour Pressure (es) 
%Numerous approximations available in the literature, e.g. Bolton, 1980

%This function expects to have temperature as an input, in Kelvins, e.g. if X is a vector input
%to the function and temperature is the second element, it would be defined
%as X(2). You will need to modify to fit into your own modelling framework.

%Coefficients
   A = 53.67957;
   B = 6743.769;
   C = 4.8451;

   es = 6.112 * exp( 17.65 * (T - 273.0)/(T - 29.5) );    %! NB: pressure in hPa 
   %Alternative Formulation 
   %es = exp(A - B/X(2) - C * log(X(2)));                       %! NB: pressure in hPa

   es = 100.0 * es;                      			%! Convert from hPa to Pa               


end

