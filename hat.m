function xhat = hat(x)
% Returns the hat map, mapping a 3-vector into a 3X3 skew-symmetric matrix
% equivalent to the cross-product operation

xhat = [  0   -x(3)  x(2)
         x(3)   0   -x(1)
        -x(2)  x(1)  0];
end
