function [J] = Jacobian_matrix(nDG,nPQ,nbus,Rv,V,Gbus,PL0,alpha)

J1 = zeros(nDG, nbus);
J2 = zeros(nDG, nDG);
J3 = zeros(nDG, nbus);
J4 = zeros(nDG, nDG);
J5 = zeros(nPQ, nbus);
J6 = zeros(nPQ, nDG);

% Calculation of J1 
for i = 1:nDG
    J1(i,i) = (Rv(i)/V(i).^2) - 1;
end

% Calculation of J2
for i = 1:nDG
    J2(i,i) = -Rv(i)/V(i);
end

% Calculation of J3 

for i = 1:nDG
    A = 0;
    for k=1:nbus
        if(k~=i)
            A = A + Gbus(i,k)*V(k);
        end
    end
    for j = 1:nbus
        if(j~=i)
            J3(i,j) = -Gbus(i,j)*V(i);            
        end
        
        J3(i,i) = -2*V(i)*Gbus(i,i) - A;
    end
end

% Calculation of J4
for i = 1:nDG
    J4(i,i) = 1;
end

% Calculation of J5

for i = 1:nPQ
    B = 0;
    for k=1:nbus
        if(k~=i+nDG)
            B = B + Gbus(i+nDG,k)*V(k);
        end
    end
    for j = 1:nbus
        if(j~=i+nDG)
            J5(i,j) = -Gbus(i+nDG,j)*V(j);
        end
        
        %J5(i,i+nDG) = -2*V(i+nDG)*Gbus(i+nDG,i+nDG) - B;
        J5(i,i+nDG) = -alpha*PL0(i+nDG)*V(i+nDG)^(alpha-1) -2*V(i+nDG)*Gbus(i+nDG,i+nDG) - B; % variable load
    end
end

% Jacobian Matrix J
J = [J1 J2; J3 J4; J5 J6];