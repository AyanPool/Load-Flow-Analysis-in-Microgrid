% LOAD FLOW MG USING GAUSS SEIDEL
clc;
clear all;
clear global
MG_system = 1;% = 1 for 6 bus;  = 2 for 38 bus
[LF,LT,nbus,nDG,nline,r,x0,mp,nq,SL,PGmax,QGmax] = system_data_THESIS(MG_system);  
SL;
nPQ = nbus-nDG;

w0 = 1.0;
w = w0; % intital value of frequency assumed
V0 = ones(nbus,1);% specified bus voltages
V  = V0 ; % intital value of voltage; 
SG = zeros(nbus,1);
PG = zeros(nbus,1);
QG = zeros(nbus,1);

mpinv = 1./mp; 
nqinv = 1./nq;

C = sum(mpinv);

error = 0.1;
iter = 0;

%% Load Power
PL = real(SL);
QL = imag(SL);
PLOAD = sum(PL);

while (error > 1.0e-10)

    z_line0 = complex (r,x0*w);
    y_line0 = 1./z_line0;

    [Ybus] = Ybus_matrix(LF,LT,y_line0,nbus,nline);


%% calculation of PG ang QG using 12 and 13
    for k=1:nDG

        PG(k+nPQ) = mpinv(k)*(w0-w);
        QG(k+nPQ) = nqinv(k)*(V0(k)-abs(V(k)));

        if( PG(k+nPQ) > PGmax(k))
            PG(k+nPQ) = PGmax(k);
        end
        if(PG(k+nPQ) < 0)
            PG(k+nPQ) = 0;
        end
        if( QG(k+nPQ) > QGmax(k))
            QG(k+nPQ) = QGmax(k);
        end
        if(QG(k+nPQ) < 0)
            QG(k+nPQ) = 0;
        end
    
    end 
    SG = complex(PG,QG); 


%% calculate Pinj and Qinj
    Sinj = SG-SL;

%% Voltage calculation    
    for k=1:nbus
        A = 0;
        for n=1:nbus
            if n~=k
                A = A + Ybus(k,n)*V(n);
            end 
        end
        Vnew(k) = (1/Ybus(k,k)) * (conj(Sinj(k))/conj(V(k)) - A);  
        DV(k) = abs(Vnew(k) - V(k));
        V(k) = Vnew(k);
    end

%% calculation of line current
    I_line = (V(LF) - V(LT)).*y_line0;

%% Total Active Power Loss calculation
    Ploss = abs(I_line).^2.*r;
    PLOSS = sum(Ploss);

%% calculate w using eq 16
    w_new = (C*w0 - (PLOAD - PLOSS))/C;
    DV = abs(w_new - w);
    w = w_new;
    error = max(DV); % this is error

    iter = iter+1;
    err(iter) = error;
    omega(iter) =  w_new;

    VM(iter,:) = abs(V);
    PGM(iter,:) = PG;
    QGM(iter,:) = QG;
end    

%% Graphs

x = 1:iter;

figure(1);
plot(x, log10(err), 'LineWidth', 2); 
xlim = ([1,iter]);
xlabel('Iterations','FontSize',20); % Label for X-axis
ylabel('log_{10} error','FontSize',20); % Label for Y-axis
title('Variation of error with iterations'); % Title for the plot
grid on; % Display grid

figure(2);
plot(x,VM, 'LineWidth', 2)
xlabel('Iterations','FontSize',20); % Label for X-axis
ylabel('V_{BUS}','FontSize',20); % Label for Y-axis
title('Variation of Bus Voltage with iterations'); % Title for the plot
grid on; % Display grid
% for i = 1:nbus
%     legend_labels{i} = ['V_ ' num2str(i)]; % Generating legend labels dynamically
% end
% legend(legend_labels);

figure(3);
plot(x, PGM(:,nPQ+1:nbus), 'LineWidth', 2); 
xlim = ([1,iter]);
xlabel('Iterations','FontSize',20); % Label for X-axis
ylabel('PG','FontSize',20); % Label for Y-axis
title('Variation of Active Power Generated with iterations'); % Title for the plot
grid on; % Display grid
for i = 1:nDG
    legend_labels{i} = ['PG_ ' num2str(i+nPQ)]; % Generating legend labels dynamically
end
legend(legend_labels);

figure(4);
plot(x, QGM(:,nPQ+1:nbus), 'LineWidth', 2); 
xlim = ([1,iter]);
xlabel('Iterations','FontSize',20); % Label for X-axis
ylabel('QG','FontSize',20); % Label for Y-axis
title('Variation of Reactive Power Generated with iterations'); % Title for the plot
grid on; % Display grid
for i = 1:nDG
    legend_labels{i} = ['QG_ ' num2str(i+nPQ)]; % Generating legend labels dynamically
end
legend(legend_labels);

figure(5);
plot(x, omega, 'LineWidth', 2); 
xlim = ([1,iter]);
xlabel('Iterations','FontSize',20); % Label for X-axis
ylabel('w','FontSize',20); % Label for Y-axis
title('Variation of angular frequency with iterations'); % Title for the plot
grid on; % Display grid