%% AC MICRO GRID load flow program using Gauss Seidel 
% REFERENCE
%   A Simple and Accurate Approach to Solve the
%   Power Flow for Balanced Islanded Microgrids
%   2015 IEEE 15th International Conference on Environment and 
%   Electrical Engineering (EEEIC) 10-13 June 2015
%   DOI: 10.1109/EEEIC35013.2015
%--------------------------------------------------------------------------

clc;
clear ;
clear global
MG_system = 2;% = 1 for 6 bus;  = 2 for 38 bus
alpha = 2 ; beta = 2;kp = 1; kq = -1; % DATA FOR DEPENDENT LOADS
[LF,LT,nbus,BTYPE,nDG,nline,r,x0,mp,nq,PL0,QL0,PGmax,Qmax] = system_data_WITH_PV(MG_system);  


nPQ = nbus-nDG;
w0 = 1.0;
w = w0; % intital value of frequency assumed
V0 = 1.0*ones(nbus,1);% specified bus voltages
V0(nPQ+1:nbus) = 1.01;% Reference voltage for Generator buses
 

V  = V0 ; % intital value of voltage; 
VOLD = V0;
SG = zeros(nbus,1);
PG = zeros(nbus,1);
QG = zeros(nbus,1);
PL = PL0;
QL = QL0;
SL = complex(PL0,QL0);
Sinj = SG-SL;
z_line0 = complex (r,x0);
y_line0 = 1./z_line0;
mpinv = 1./mp; 
nqinv = 1./nq;
QGmax = zeros(nbus,1);
for i = nPQ+1 : nbus
    if (BTYPE(i) == 2)
        QGmax(i) = nqinv(i-nPQ)*0.01 ;% RECALCULATE THE QG LIMIT FOR PV BUSES
    else
        QGmax(i) = Qmax(i-nPQ);% QGMAX for PQ GEN BUSES
    end
end
error_max = 1e03;
mm = 0;

     %%

% 

%while ( mm < 2)
 while ( error_max > 1.0e-8)
    mm = mm+1;    % no of iterations
%        
       z_line = complex(r,x0*w);  % scale impedance  by w/w0       
       y_line = 1./z_line;   
       [Ybus] = Ybus_matrix(LF,LT,y_line,nbus,nline);
% %     

for i = 1:nbus
    FLAG = 0;
         xx = 1.0/Ybus(i,i) ;
         sum1 = 0.0;
         for j = 1:nbus
             if( j ~= i)
                 sum1 = sum1 + Ybus(i,j)*V(j);
             end
         end
         % checking for PV bus 
         if(BTYPE(i)== 2)
            
             VT = V;
             thi = angle(V(i));
             VT(i) = 1*complex(cos(thi),sin(thi));
             Ii  = (Ybus(i,:)* VT);
             Si = VT(i)*conj(Ii);
             Qi = imag(Si)  ;
            
             if(Qi > QGmax(i))
                 Qi = QGmax(i);
                 'MAX violated'
                 FLAG = 0;
             else 
                 if(Qi < -QGmax(i))
                     Qi = -QGmax(i);
                     'MIN violated'
                     FLAG = 0;
                 else
                     FLAG = 1 ;% THE Q IS WITHIN LIMITS AND IT IS A PV BUS
                 end      
             end
% set the injected power for PV bus for voltage calculations

            Sinj(i) = complex(real(SG(i)-SL(i)), Qi);    
         end
            V(i) = xx*(conj(Sinj(i)/V(i))-sum1) ;

         if( FLAG == 1&& BTYPE(i)== 2)
             thi = angle(V(i));
             V(i) = 1 * complex(cos(thi),sin(thi)) ;
         end

%% to ensure V(1) at angle 0.0 use the following statements
%  else subtract the angle of V(1) from all the voltages
%  
         if(i==1)  
             V(1) = abs(V(1))*complex(1,0); 
         end
         % update LOAD VALUES FOR V AND OMEGA DEPENDENT LOADS
         if(alpha ~= 0 ) 
               PL(i) = PL0(i)*(abs(V(i))/V0(i))^alpha*(1+kp*(w-w0));
               QL(i) = QL0(i)*(abs(V(i))/V0(i))^alpha*(1+kq*(w-w0));
               SL(i) = complex(PL(i),QL(i));
         end
end

%      update QG 
        for i = nPQ+1:nbus
             i;
            if (BTYPE(i) == 2)   
                QG(i) = Qi + QL(i)  ;
            else
                QG(i) = nqinv(i-nPQ)*(V0(i)-abs(V(i))) ; 
            end
        end

       Vm(:,mm) = abs(V) ;           
% %       
      PLT = sum (real(SL));
      QLT = sum (imag(SL));
%   % line loss calculations
      Iij = (V(LF)-V(LT)).* y_line;
      Sij =  V(LF).*conj(Iij);
      Sji = -V(LT).*conj(Iij);
      Sline_loss =  (abs(Iij).^2).*z_line;% just testing
      line_loss =   Sij+Sji;

      PLOSS = sum(real(line_loss));
      QLOSS = sum(imag(line_loss));
      PLS(mm) = PLOSS;
      QLS(mm) = QLOSS;
%  

%      %% calculated Generation

%      
      w = (sum(mpinv)*w0-(PLT+PLOSS))/sum(mpinv) ;
      omega(mm) = w;
      xaxs(mm) = mm-1;
     
     PG(nPQ+1:nbus)  = mpinv.*(w0-w);

     for i = nPQ +1 : nbus
         [PG QG];
         if(  PG (i) > PGmax(i-nPQ)) 
              PG (i)= PGmax(i-nPQ);
         end
         if(  QG (i) > QGmax(i))
              QG (i)= QGmax(i); 
         else
         if ( QG (i) < -QGmax(i))
              QG (i) = -QGmax(i);
         end
         end
         [PG QG];
     end
     
    SG = complex(PG,QG)   ;
    Sinj = SG-SL  ;
      PGP(:,mm) = PG(nPQ+1:nbus)   ;  
      QGP(:,mm) = QG(nPQ+1:nbus) ;
%    
%      
% calculation of max error
  if(mm >= 2)
  eVmax(mm)      = max (abs(Vm(:,mm)  -Vm(:,mm-1)));
  ePGmax (mm)    = max (abs(PGP(:,mm)-PGP(:,mm-1)));
  eQGmax (mm)    = max (abs(QGP(:,mm)-QGP(:,mm-1))); 
  error_max    = max(max(eVmax(mm),ePGmax(mm)),eQGmax(mm));
   max_err(mm) = error_max;
 end    
%     
end
% %  
figure(1) 
 plot(xaxs,  log10((max_err)),'LineWidth',2) ;
 ax = gca;
ax.FontSize = 13; 
xlabel('Iterations','FontSize',14);
ylabel('LOG_{10} error ','FontSize',14);
title('Variation of error with iterations','Color','k','FontSize',14)
 grid on

figure(2)
plot(xaxs,omega,'Color','b','LineWidth',2)
ax = gca;
ax.FontSize = 13; 
xlabel('Iterations','FontSize',14);
ylabel('w','FontSize',14);
title('Frquency variation with iterations','Color','k','FontSize',14)
grid on

figure(3)
plot(xaxs,Vm,'LineWidth',2)
ax = gca;
ax.FontSize = 13; 
xlabel('Iterations','FontSize',14);
ylabel('V_{BUS}','FontSize',14);
title('Variation of Bus Voltage with iterations','Color','k','FontSize',14)
if(MG_system == 1)
    legend('V_1','V_2','V_3','V_4','V_5','V_6')
end
grid on
figure(4)
grid on

plot(xaxs,PGP,'LineWidth',2)
ax = gca;
ax.FontSize = 13; 
xlabel('Iterations','FontSize',14);
ylabel('PG','FontSize',14);
title('Variation of Active Power Generated with iterations','Color','k','FontSize',14)
% if(MG_system == 1)
% legend('PG_1','PG_2','PG_3')
% else
% legend('PG_34','PG_35','PG_36','PG_37','PG_38')
% end
for i = 1:nDG
    legend_labels{i} = ['PG_ ' num2str(i+nPQ)]; % Generating legend labels dynamically
end
legend(legend_labels);
grid on

figure(5)
grid on
plot(xaxs,QGP,'LineWidth',2)
ax = gca;
ax.FontSize = 13; 
xlabel('Iterations','FontSize',14);
ylabel('QG','FontSize',14);
title('Variation of Reactive Power Generated with iterations','Color','k','FontSize',14)
% if(MG_system == 1)
% legend('QG_1','QG_2','QG_3')
% else
% legend('QG_34','QG_35','QG_36','QG_37','QG_38')
% end
for i = 1:nDG
    legend_labels{i} = ['QG_ ' num2str(i+nPQ)]; % Generating legend labels dynamically
end
legend(legend_labels);
grid on

figure(6)
grid on
plot(xaxs,PLS,'LineWidth',2)
hold on
plot(xaxs,QLS,'LineWidth',2)
ax = gca;
ax.FontSize = 13; 
xlabel('Iteration count','FontSize',14);
ylabel('P_loss/Q_loss (pu)','FontSize',14);
title('Loss variation with iterations','Color','k','FontSize',14)

legend('P_{LOSS}','Q_{LOSS}')
grid on


  IB = Ybus*V;
  SB = V.*conj(IB);
[abs(V)  angle(V)*180/pi;];
[abs(V)  angle(V)*180/pi-angle(V(1))*180/pi;];
VR = abs(V).*complex(cos(angle(V)-angle(V(1))),sin(angle(V)-angle(V(1))));