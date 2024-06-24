function [Gbus] = Gbus_matrix(LF,LT,y_line0,nbus,nline)
% function to calculate Ybus matrix
Gbus = zeros(nbus);

for i = 1:nline
    
    Gbus(LF(i),LF(i)) = Gbus(LF(i),LF(i)) + y_line0(i) ;
    Gbus(LT(i),LT(i)) = Gbus(LT(i),LT(i)) + y_line0(i);
    Gbus(LF(i),LT(i)) = Gbus(LF(i),LT(i)) - y_line0(i);
    Gbus(LT(i),LF(i)) = Gbus(LF(i),LT(i));
    
end