function [X_hat] = mod_slicer(X, mod_scheme)

N = length(X);
% modScheme
switch mod_scheme
    case 'QPSK'
        Xr = real(X);  Xi = imag(X);
        for i=1:N
            R(Xr<0) = -1;  I(Xi<0) = -1;
            R(0<=Xr) = 1;  I(0<=Xi) = 1;
        end
    case '16QAM'
        b = [-2 0 2];
        c = [-3 -1 1 3];
        Xr = real(X);  Xi = imag(X);
        for i=1:N
            R(Xr<b(1)) = c(1);  I(Xi<b(1)) = c(1);
            R(b(1)<=Xr&Xr<b(2)) = c(2);  I(b(1)<=Xi&Xi<b(2)) = c(2);
            R(b(2)<=Xr&Xr<b(3)) = c(3);  I(b(2)<=Xi&Xi<b(3)) = c(3);
            R(b(3)<=Xr) = c(4);  I(b(3)<=Xi) = c(4);
        end
end
X_hat = R;% + 1i*I;