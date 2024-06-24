clc;
clear all;
MG_system = 1;
[nDG,nPQ,Rv,R_line,PL0,r,LF,LT,nbus,nline] = system_data_DC(MG_system);
V0 = ones(nbus,1);
V = V0;
PG = ones(nDG,1);
z_line0 = r;
y_line0 = 1./z_line0;

alpha = 1.5;
PL = PL0 .* V.^alpha;

xold = [V0' PG'];
error = 0.1;
iter = 0;

[Gbus] = Gbus_matrix(LF,LT,y_line0,nbus,nline);
 
while(error > 1.0e-10)

    [J] = Jacobian_matrix(nDG,nPQ,nbus,Rv,V,Gbus,PL0,alpha);

    for i = 1:nDG
        f(i,1) = V0(i) - V(i) - Rv(i)*(PG(i)/V(i));
    end

    
    for i = 1:nDG
        A = 0;
        for j = 1:nbus
            if(j~=i)
                A = A + Gbus(i,j)*V(j);
            end
        end
        f(i+nDG,1) = PG(i) - PL(i) - Gbus(i,i)*V(i).^2 - V(i)*A ;
    end

   
    for i = 1:nPQ 
        C = 0;
        for j = 1:nbus
            if(j~=i+nDG)
                C = C + Gbus(i+nDG,j)*V(j);
            end
        end
        f(i+2*nDG,1) = - PL(i+nDG) - Gbus(i+nDG,i+nDG)*V(i+nDG).^2 - V(i+nDG)*C;
    end
    
     delX = -J\f;
     x = xold + delX';
     xold = x;
      
     error = max(abs(f)); 
     iter = iter+1;

     V = x(1:nbus);
    PG = x(nbus+1:nbus+nDG);

    err(iter) = error;
    VM(iter,:) = abs(V);
    PGM(iter,:) = PG;

end
    
%% Graphs

y = 1:iter;

figure(1);
plot(y, log10(err), 'LineWidth', 2); 
xlim = ([1,iter]);
xlabel('Iterations','FontSize',20); % Label for X-axis
ylabel('log_{10} error','FontSize',20); % Label for Y-axis
title('Variation of error with iterations'); % Title for the plot
grid on; % Display grid

figure(2);
plot(y,VM, 'LineWidth', 2)
xlabel('Iterations','FontSize',20); % Label for X-axis
ylabel('V_{BUS}','FontSize',20); % Label for Y-axis
title('Variation of Bus Voltage with iterations'); % Title for the plot
grid on; % Display grid
legend('V_1','V_2','V_3','V_4','V_5','V_6')

figure(3);
plot(y, PGM, 'LineWidth', 2); 
xlim = ([1,iter]);
xlabel('Iterations','FontSize',20); % Label for X-axis
ylabel('PG','FontSize',20); % Label for Y-axis
title('Variation of Active Power Generated with iterations'); % Title for the plot
grid on; % Display grid
legend('PG_1','PG_2','PG_3')