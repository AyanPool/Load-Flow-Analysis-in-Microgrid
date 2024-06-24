function [Ybus] = Ybus_matrix(LF,LT,y_line0,nbus,nline)
% function to calculate Ybus matrix
Ybus = zeros(nbus);

for i = 1:nline
    
    Ybus(LF(i),LF(i)) = Ybus(LF(i),LF(i)) + y_line0(i) ;
    Ybus(LT(i),LT(i)) = Ybus(LT(i),LT(i)) + y_line0(i);
    Ybus(LF(i),LT(i)) = Ybus(LF(i),LT(i)) - y_line0(i);
    Ybus(LT(i),LF(i)) = Ybus(LF(i),LT(i));
    
end