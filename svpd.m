% ------------------------------------------------------------------------
% Compute saturation pressure over water surface
%
% Reference:
%
% Flatau,P.J., R.L. Walko, and W.R. Cotton, 1992: Polynomial fits to
%        saturation vapor pressure. J. Applied Meteorology, 31, 1507-1513.
%
% AOSC class note
%
% ------------------------------------------------------------------------

function [svpd] = svpd(temp)
global ps;

%  Input:
    % temp = temperature [K];
%
%  Output:
 %    vprs = saturation vapor pressure [Pa];
%
%  fitting coefficients from table 4 of Flatau et al. (1992)
%  valid for temperature range [-85,70]C
%
   a0 = 0.611583699e+3;
   a1 = 0.444606896e+2;
   a2 = 0.143177157e+1;
   a3 = 0.264224321e-1;
   a4 = 0.299291081e-3;
   a5 = 0.203154182e-5;
   a6 = 0.702620698e-8;
   a7 = 0.379534310e-11;
   a8 =-0.321582393e-13;
%
%  convert temperature from input kelvin to Celsius
%
   tc = temp - 273.15;
%
%  limit to the valid range
%
tc = max(-85.0, min(+70.0, tc));
%
%  below we still allow the calculation even beyond the valid bound
%  but print out a warning message
   if tc < -85.0
      disp('temperature < low limit -85 for valid output')
   elseif tc > +70.0
      disp('temperature > upper limit +70C for valid output')
   end
%
%  compute the saturation pressure by a polynomial fit
%
   svpd = a0+tc.*(a1+tc.*(a2+tc.*(a3+tc.*(a4+tc.*(a5+tc.*(a6+tc.*(a7+tc.*a8)))))))-ps;
%
%  all done
%
disp(svpd)
end