%% Teoria sterowania (TST) Projekt 2 Arkadiusz Piorkowski
clear all;

%% Dane 
T = 8e3;
PHI = [0.25 0.5;1 2];
PHI_xw = [0 1;0 1];
PHI_w = [sqrt(0.5) sqrt(0.5); -sqrt(0.5) sqrt(0.5)];
GAMMA = [1;4];
GAMMA_ni = [1;1];
C = [1 0];
n = 2;
nw = 2;
m = 1;
A = [PHI PHI_xw; zeros(nw,n) PHI_w];
B = [GAMMA;zeros(nw,m)];
C = [C ,zeros(1,nw)];
E = [zeros(n,1);GAMMA_ni];


%% Zaklocenie
% Wpisujemy wartosci Rv Re i Rve
Rv=1;%rozklad normalny (gaussowski) N(0,Rv)/czy to macierze kowariancji Q?
Re=1;%rozklad normalny (gaussowski) N(0,Re)/czy to macierze kowariancji R?
%Rve=0;%-0.072;

%lub obliczamy Rve z macierzy kowariancji przy zalozeniu ze randn to
%idealny rozklad normalny
v=sqrt(Rv)*randn(T,1);
e=sqrt(Re)*randn(T,1);
COV_VE=cov(v,e);
Rve=COV_VE(1,2);

%% Obserwator L typu deadbeat
Wo=obsv(A,C);%lub ctrb(A',C');
rank_Wo=rank(Wo);

Zo = [0; 0; 0; 0];
L_db = acker(A',C',Zo)';


%% Obserwator L minimalizujacy wariancje bledu estymacji stanu- filtr Kalmana
sys = ss(A,[B, ones(n+nw,1)],C,[],-1);
[Kest,L_kal,P_kal] = kalman(sys,Rv,Re,Rve);

% Sprawdzenie
P{1}=ones(n+nw,n+nw);
L{1}=zeros(n+nw,1);
for t=1:T-1
    L{t} = (A*P{t}*C'+Rve)/(C*P{t}*C'+Re);
    P{t+1} = A*P{t}*A'+Rv-L{t}*(C*P{t}*A'+Rve');
end

% Wyswietlenie P_kal i L_kal
%{
%reshape([P_kal{:}],4,4*T)
%reshape([L_kal{:}],4,T-1)
%P_kal{:},L_kal{:}
%}

%% Regulator K_x typu deadbeat

% Sterowalnosc - sprawdzenie
Wc=ctrb(PHI,GAMMA);%Wc=ctrb(A,B)
rank_Wc=rank(Wc);%=1
%zapewne bedzie trzeba wykorzystac pseudoodwrotnosc w celu zrzutowania na
%podprzestrzen


% Osobliwosc macierzy PHI (A) - nie mozna wykonac polecenia acker
%{
% Zc=[0;0];
% Kdb = acker(PHI,GAMMA,Zc);%nie mozna wykonac gdyz macierz PHI jest osobliwa
%(dla det(Iz-A+BK)=z^4 => (z*(z^2 - 2^(1/2)*z + 1)*(4*k1 + 16*k2 + 4*z - 9))/4=z^4)
% zatem przez porownanie wspolczynnikow det(Iz-PHI+GAMMA*K)=z^2
% k1*z - (9*z)/4 + 4*k2*z + z^2=z^2
% z^2 + (k1 - 9/4 + 4*k2)*z=z^2
% k1=9/4-4*k2 - dla dowolnego k2
%}
k2 = 0.5;
k1 = 9/4-4*k2;
Kx_db_check = [k1 k2];
%Wykonanie przez pseudoodwrotnosc:
%{
%Aby wykonac regulator deadbeat mozna rowaniez skorzystac z przyrownania
%rownania macierzowego -PHI+GAMMA*K do 0 : -PHI+GAMMA*K=0 w wyniku czego
%otrzymujemy rownanie na K w postaci K=pinv(GAMMA)*PHI - pseudoodwrotnosc
%ze wzgledu na brak mozliwosci odwrocenia macierzy GAMMA
Kx_db=pinv(GAMMA)*PHI; %-otrzymujemy ten sam wynik jak powyzej;
%}
Kx_db=pinv(GAMMA)*PHI;

%% Regulator K_x liniowo-kwadratowy

% Sterowanie optymalne Q (LQ)
Qx=2*eye(n);
QN=Qx;%koszt koncowy Jpi(x)=E sum([x' u'][Qx Qxu;Qxu' Qu][x; u]+xn'*Qn*xn);
Qu=0.1;
Qxu=zeros(n,m);
%[K,S] = dlqr(A,B,Qx,Qu,Qxu)
[Kx_lq, Sx_lq] = dlqr(PHI,GAMMA,Qx,Qu,Qxu);

% Sprawdzenie
%{
% S{T}=QN;
% for t=T-1:-1:1
%     K{t}=(B'*S{t+1}*B+Qu)\(B'*S{t+1}*A+Qxu');
%     S{t}=A'*S{t+1}*A+Qx-(A'*S{t+1}*B+Qxu)*K{t};
% end
%}
S_lq{T}=QN;
for t=T-1:-1:1
    K_lq{t}=(GAMMA'*S_lq{t+1}*GAMMA+Qu)\(GAMMA'*S_lq{t+1}*PHI+Qxu');
    S_lq{t}=PHI'*S_lq{t+1}*PHI+Qx-(PHI'*S_lq{t+1}*GAMMA+Qxu)*K_lq{t};
end

% Wyswietlenie K_lq
%{
% reshape([K_lq{:}],2,T-1);
% K_lq{:}
%}

%% Macierz K_w minimalizujaca wplyw zaklocen w(t) na stan x(t)

% Wyprowadzenie
%{
%podstawiajac do rownanaia:
%[x(t+1),w(t+1)]'=A*[x(t),w(t)]'+B*u(t)
%polityke sterowania w postaci: u(t)=-[Kx Kw]*[x(t),w(t)]'
%otrzymano równania:
%x(t+1)=(PHI-GAMMA*Kx)*x(t)+(PHI_xw-GAMMA*Kw)*w(t);
%w(t+1)=PHI_w*w(t);
%aby by³ minimalizowany wplyw zaklocen w(t) na stan x(t) nalezy przyrownac 
%wspolczynnik (PHI_xw-GAMMA*Kw) do 0
%PHI_xw-GAMMA*Kw=0 => GAMMA*Kw=PHI_xw => Kw = inv(GAMMA)*PHI_xw;
%odwrocenie jest niemozliwe, zatem trzeba wykorzystac pseudoodwrotnosc
%taki sam wynik otrzymuje sie stosujac dzielenie lewostronne Kw = GAMMA\PHI_xw;
%}
K_w = pinv(GAMMA)*PHI_xw;

%Wniosek dotyczacy regulatora deadbeat i K_w
%{
% Wykonujac zatem dzialanie K=pinv(A)*B otrzymujemy K=[Kx Kw] gdzie Kx
% oznacza regulator typu deadbeat:
% [x(t+1),w(t+1)]'=A*[x(t),w(t)]'+B*u(t)
% u(t)=-K*[x(t),w(t)]'
% y(t+1)=C[A*[x(t),w(t)]'-B*K*[x(t),w(t)]']
% y(t+1)=C*[A-B*K]*[x(t),w(t)]'
% A-BK=0;-regulator deadbeat w odniesieniu do x i w
% K=pinv(B)*A;
%}

%%  Petla procesu sterowania
%{
for i=1:T-1
%proces
[x(:,t+1);w(:,t+1)]=A*[x(:,t);w(:,t)]+B*u(t)+F*v(t);
y(:,t)=C*[x(:,t);w(:,t)]
%obserwator
[xest(:,t+1);west(:,t+1)] = A*[xest(:,t);west(:,t)]+B*u(t)+L*(y(:,t)-C*[xest(:,t);west(:,t)]);
%regulator  
u(t) = -[Kx K_w]*[xest(:,t+1);west(:,t+1)];
end
%ocena estymacji:
Error=[xest(:,:);west(:,:)]-[x(:,:);w(:,:)];
%}

% Inicjalizacja - wektory dla obserwatora i regulatora deadbeat
z_db_db(:,1)=[0;0;v(t);v(t)];
z_db_db_est = zeros(n+nw,T);
u_db_db = zeros(m,T);
y_db_db = zeros(1,T);
% Inicjalizacja - wektory dla obserwatora-filtr Kalmana i regulatora LQ
z_kal_lq=[0;0;v(t);v(t)];
z_kal_lq_est = zeros(n+nw,T);
u_kal_lq = zeros(m,T);
y_kal_lq = zeros(1,T);
% Inicjalizacja - wektory dla obserwatora deadbeat i regulatora LQ
z_db_lq=[0;0;v(t);v(t)];
z_db_lq_est = zeros(n+nw,T);
u_db_lq = zeros(m,T);
y_db_lq = zeros(1,T);
% Inicjalizacja - wektory dla obserwatora-filtr Kalmana i regulatora deadbeat
z_kal_db(:,1)=[0;0;v(t);v(t)];
z_kal_db_est = zeros(n+nw,T);
u_kal_db = zeros(m,T);
y_kal_db = zeros(1,T);

for t = 1:T-1
    
    % proces (obiekt sterowania) - obserwator deadbeat + regulator deadbeat
    z_db_db(:,t+1) = A*z_db_db(:,t)+B*u_db_db(t)+E*v(t);
    y_db_db(:,t) = C*z_db_db(:,t)+e(t);    
    % rownanie obserwatora deadbeat
    z_db_db_est(:,t+1) = A*z_db_db_est(:,t)+B*u_db_db(t)+L_db*(y_db_db(:,t)-C*z_db_db_est(:,t));
    % rownanie regulatora deadbeat  
    u_db_db(t+1) = -[Kx_db K_w]*z_db_db_est(:,t+1);    
    
    % proces (obiekt sterowania) - filtr Kalmana + regulator LQ
    z_kal_lq(:,t+1) = A*z_kal_lq(:,t)+B*u_kal_lq(t)+E*v(t);
    y_kal_lq(:,t) = C*z_kal_lq(:,t)+e(t);
    % rownanie obserwatora - filtr Kalmana
    z_kal_lq_est(:,t+1) = A*z_kal_lq_est(:,t)+B*u_kal_lq(t)+L_kal*(y_kal_lq(:,t)-C*z_kal_lq_est(:,t)); 
    % rownanie regulatora LQ
    u_kal_lq(t+1) = -[Kx_lq K_w]* z_kal_lq_est(:,t+1);  
    
    % proces (obiekt sterowania) - obserwator deadbeat + regulator LQ
    z_db_lq(:,t+1) = A*z_db_lq(:,t)+B*u_db_lq(t)+E*v(t);
    y_db_lq(:,t) = C*z_db_lq(:,t)+e(t);    
    % rownanie obserwatora deadbeat
    z_db_lq_est(:,t+1) = A*z_db_lq_est(:,t)+B*u_db_lq(t)+L_db*(y_db_lq(:,t)-C*z_db_lq_est(:,t));
    % rownanie regulatora LQ  
    u_db_lq(t+1) = -[Kx_lq K_w]*z_db_lq_est(:,t+1);    
    
    % proces (obiekt sterowania) - filtr Kalmana + regulator deadbeat
    z_kal_db(:,t+1) = A*z_kal_db(:,t)+B*u_kal_db(t)+E*v(t);
    y_kal_db(:,t) = C*z_kal_db(:,t)+e(t);
    % rownanie obserwatora - filtr Kalmana
    z_kal_db_est(:,t+1) = A*z_kal_db_est(:,t)+B*u_kal_db(t)+L_kal*(y_kal_db(:,t)-C*z_kal_db_est(:,t)); 
    % rownanie regulatora deadbeat
    u_kal_db(t+1) = -[Kx_db K_w]* z_kal_db_est(:,t+1);  
       
end

% Bledy estymacji
Error_db_db=z_db_db-z_db_db_est;
Error_kal_lq=z_kal_lq-z_kal_lq_est;
Error_db_lq=z_db_lq-z_db_lq_est;
Error_kal_db=z_kal_db-z_kal_db_est;
% Bledy kwadratowe(do oceny)
Error_db_db_sqr(1)=Error_db_db(1,:)*Error_db_db(1,:)';
Error_kal_lq_sqr(1)=Error_kal_lq(1,:)*Error_kal_lq(1,:)';
Error_db_lq_sqr(1)=Error_db_lq(1,:)*Error_db_lq(1,:)';
Error_kal_db_sqr(1)=Error_kal_db(1,:)*Error_kal_db(1,:)';
Error_db_db_sqr(2)=Error_db_db(2,:)*Error_db_db(2,:)';
Error_kal_lq_sqr(2)=Error_kal_lq(2,:)*Error_kal_lq(2,:)';
Error_db_lq_sqr(2)=Error_db_lq(2,:)*Error_db_lq(2,:)';
Error_kal_db_sqr(2)=Error_kal_db(2,:)*Error_kal_db(2,:)';

%% Wykresy - niepotrzebne


%% Zewnêtrzny sygna³ sterowania (feedforward)

% fala prostokatna
u_c=zeros(m,T);
okres=300;
amplituda=300;
licznik=0;
for i=1:T
    u_c(i)= amplituda*(2*mod(licznik,2)-1);   
    if mod(i,okres) == 0
        licznik=licznik+1;
    end
end

% Macierz K_c - wyjscie podaza za wejsciem-wyprowadzenie
%{
% wyprowadzone przy zalozeniu ze zastosowano regulator deadbeat
% x(t+1)=PHI*x(t)+GAMMA*u(t)
% y(t)=C*x(t)
% u(t)=u_fb(t)+u_ff(t)
% u(t)=-Kx*x(t)+Kc*u_c(t)
% Wprowadzajac do rownania stanu
% x(t+1)=PHI*x(t)+GAMMA*[-Kx*x(t)+Kc*u_c(t)]
% x(t+1)=(PHI-GAMMA*Kx)*x(t)+GAMMA*Kc*u_c(t);
% y(t+1)=C[(PHI-GAMMA*Kx)*x(t)+GAMMA*Kc*u_c(t)];
% y(t+1)=C(PHI-GAMMA*Kx)*x(t) + C*GAMMA*Kc*u_c(t);
% ze wzgledu iz Kx zostalo dobrane tak aby (PHI-GAMMA*Kx)=0, zatem:
% y(t+1)=C*GAMMA*Kc*u_c(t); 
% sygnal wyjsciowy ma podazac za sygnalem sterujacym, zatem C*GAMMA*Kc=1;
% Kc=inv(C*GAMMA)*1;
%}
K_c=(C(1:2)*GAMMA)\1;

% Petla procesu sterowania

% Inicjalizacja - wektory dla obserwatora i regulatora deadbeat
z_db_db_ff(:,1) = [0;0;v(t);v(t)];
z_db_db_ff_est = zeros(n+nw,T);
u_db_db_ff = zeros(m,T);
y_db_db_ff = zeros(1,T);
% Inicjalizacja - wektory dla obserwatora-filtr Kalmana i regulatora LQ
z_kal_lq_ff = [0;0;v(t);v(t)];
z_kal_lq_ff_est = zeros(n+nw,T);
u_kal_lq_ff = zeros(m,T);
y_kal_lq_ff = zeros(1,T);
% Inicjalizacja - wektory dla obserwatora deadbeat i regulatora LQ
z_db_lq_ff = [0;0;v(t);v(t)];
z_db_lq_ff_est = zeros(n+nw,T);
u_db_lq_ff = zeros(m,T);
y_db_lq_ff = zeros(1,T);
% Inicjalizacja - wektory dla obserwatora-filtr Kalmana i regulatora deadbeat
z_kal_db_ff(:,1) = [0;0;v(t);v(t)];
z_kal_db_ff_est = zeros(n+nw,T);
u_kal_db_ff = zeros(m,T);
y_kal_db_ff = zeros(1,T);

for t = 1:T-1
    
    % proces (obiekt sterowania) - obserwator deadbeat + regulator deadbeat
    z_db_db_ff(:,t+1) = A*z_db_db_ff(:,t)+B*u_db_db_ff(t)+E*v(t);
    y_db_db_ff(:,t) = C*z_db_db_ff(:,t)+e(t);    
    % rownanie obserwatora deadbeat
    z_db_db_ff_est(:,t+1) = A*z_db_db_ff_est(:,t)+B*u_db_db_ff(t)+L_db*(y_db_db_ff(:,t)-C*z_db_db_ff_est(:,t));
    % rownanie regulatora deadbeat  
    u_db_db_ff(t+1) = -[Kx_db K_w]*z_db_db_ff_est(:,t+1)+K_c*u_c(t);    
    
    % proces (obiekt sterowania) - filtr Kalmana + regulator LQ
    z_kal_lq_ff(:,t+1) = A*z_kal_lq_ff(:,t)+B*u_kal_lq_ff(t)+E*v(t);
    y_kal_lq_ff(:,t) = C*z_kal_lq_ff(:,t)+e(t);
    % rownanie obserwatora - filtr Kalmana
    z_kal_lq_ff_est(:,t+1) = A*z_kal_lq_ff_est(:,t)+B*u_kal_lq_ff(t)+L_kal*(y_kal_lq_ff(:,t)-C*z_kal_lq_ff_est(:,t)); 
    % rownanie regulatora LQ
    u_kal_lq_ff(t+1) = -[Kx_lq K_w]* z_kal_lq_ff_est(:,t+1)+K_c*u_c(t);  
    
    % proces (obiekt sterowania) - obserwator deadbeat + regulator LQ
    z_db_lq_ff(:,t+1) = A*z_db_lq_ff(:,t)+B*u_db_lq_ff(t)+E*v(t);
    y_db_lq_ff(:,t) = C*z_db_lq_ff(:,t)+e(t);    
    % rownanie obserwatora deadbeat
    z_db_lq_ff_est(:,t+1) = A*z_db_lq_ff_est(:,t)+B*u_db_lq_ff(t)+L_db*(y_db_lq_ff(:,t)-C*z_db_lq_ff_est(:,t));
    % rownanie regulatora LQ  
    u_db_lq_ff(t+1) = -[Kx_lq K_w]*z_db_lq_ff_est(:,t+1)+K_c*u_c(t);    
    
    % proces (obiekt sterowania) - filtr Kalmana + regulator deadbeat
    z_kal_db_ff(:,t+1) = A*z_kal_db_ff(:,t)+B*u_kal_db_ff(t)+E*v(t);
    y_kal_db_ff(:,t) = C*z_kal_db_ff(:,t)+e(t);
    % rownanie obserwatora - filtr Kalmana
    z_kal_db_ff_est(:,t+1) = A*z_kal_db_ff_est(:,t)+B*u_kal_db_ff(t)+L_kal*(y_kal_db_ff(:,t)-C*z_kal_db_ff_est(:,t)); 
    % rownanie regulatora deadbeat
    u_kal_db_ff(t+1) = -[Kx_db K_w]*z_kal_db_ff_est(:,t+1)+K_c*u_c(t);  
       
end

% Bledy estymacji
Error_db_db_ff = z_db_db_ff - z_db_db_ff_est;
Error_kal_lq_ff = z_kal_lq_ff - z_kal_lq_ff_est;
Error_db_lq_ff = z_db_lq_ff - z_db_lq_ff_est;
Error_kal_db_ff = z_kal_db_ff - z_kal_db_ff_est;
% Bledy kwadratowe(do oceny)
Error_db_db_ff_sqr(1) = Error_db_db_ff(1,:) * Error_db_db_ff(1,:)';
Error_kal_lq_ff_sqr(1) = Error_kal_lq_ff(1,:) * Error_kal_lq_ff(1,:)';
Error_db_lq_ff_sqr(1) = Error_db_lq_ff(1,:) * Error_db_lq_ff(1,:)';
Error_kal_db_ff_sqr(1) = Error_kal_db_ff(1,:) * Error_kal_db_ff(1,:)';
Error_db_db_ff_sqr(2) = Error_db_db_ff(2,:) * Error_db_db_ff(2,:)';
Error_kal_lq_ff_sqr(2) = Error_kal_lq_ff(2,:) * Error_kal_lq_ff(2,:)';
Error_db_lq_ff_sqr(2) = Error_db_lq_ff(2,:) * Error_db_lq_ff(2,:)';
Error_kal_db_ff_sqr(2) = Error_kal_db_ff(2,:) * Error_kal_db_ff(2,:)';
Error_db_db_ff_sqr(3) = (y_db_db_ff-u_c) * (y_db_db_ff-u_c)';
Error_kal_lq_ff_sqr(3) = (y_kal_lq_ff-u_c) * (y_kal_lq_ff-u_c)';
Error_db_lq_ff_sqr(3) = (y_db_lq_ff-u_c) * (y_db_lq_ff-u_c)';
Error_kal_db_ff_sqr(3) = (y_kal_db_ff-u_c) * (y_kal_db_ff-u_c)';

%% Wykresy

%wykresy dla przypadkow: kalman+LQ i deadbeat+deadbeat(dla
%innych kombinacji bledy sie powtarzaja)

fig1 = figure('Name','Zmienne stanu','units','normalized','outerposition',[0 0 0.8 0.8]);
subplot(2,1,1)
plot(z_db_db_ff_est(1,:))
hold on
plot(z_kal_lq_ff_est(1,:))
grid on
ylabel('x1_{db db},x1_{kal lq}')
xlabel('T')
legend('x1_{db db}','x1_{kal lq}')
title({'Przebieg zmiennych stanu estymowanych','x_1'})
subplot(2,1,2)
plot(z_db_db_ff_est(2,:))
hold on
plot(z_kal_lq_ff_est(2,:))
grid on
ylabel('x2_{db db},x2_{kal lq}')
xlabel('T')
legend('x2_{db db}','x2_{kal lq}')
title('x_2')


fig2 = figure('Name','Bledy estymacji stanu','units','normalized','outerposition',[0 0 0.8 0.8]);
subplot(2,1,1)
plot(Error_db_db_ff(1,:))
hold on
plot(Error_kal_lq_ff(1,:))
grid on
ylabel('Error x1_{db db}, Error x1_{kal lq}')
xlabel('T')
legend(strcat('Error x1_{db db}=',num2str(Error_db_db_ff_sqr(1),'%8.1e')),strcat('Error x1_{kal lq}=',num2str(Error_kal_lq_ff_sqr(1),'%8.1e')))
title({'Bledy estymacji stanu','x_1'})
subplot(2,1,2)
plot(Error_db_db_ff(2,:))
hold on
plot(Error_kal_lq_ff(2,:))
grid on
ylabel('Error x2_{db db}, Error x2_{kal lq}')
xlabel('T')
legend(strcat('Error x2_{db db}=',num2str(Error_db_db_ff_sqr(2),'%8.1e')),strcat('Error x2_{kal lq}=',num2str(Error_kal_lq_ff_sqr(2),'%8.1e')))
title('x_2')


fig3 = figure('Name','Wyjscie/wejscie procesu','units','normalized','outerposition',[0 0 0.8 0.8]);
subplot(2,1,1)
plot(y_db_db_ff)
hold on
plot(y_kal_lq_ff)
hold on
plot(u_c,'k')
grid on
ylabel('y_{db db},y_{kal lq},u_c')
xlabel('T')
legend(strcat('y_{db db}, Err=',num2str(Error_db_db_ff_sqr(3),'%8.2e')),strcat('y_{kal lq}, Err=',num2str(Error_kal_lq_ff_sqr(3),'%8.2e')),'u_c')
title({'Przebiegi sygnalu wyjsciowego i sterujacego','y,u_c'})
subplot(2,1,2)
plot(u_db_db_ff)
hold on
plot(u_kal_lq_ff)
grid on
ylabel('u_{db db},u_{kal lq}')
xlabel('T')
legend('u_{db db}','u_{kal lq}')
title('u_{db db}, u_{kal lq}')

%dodanie wykresow dla kombinacji filtr Kalmana + regulator deadbeat -
% ze wzgledu na najmniejszy blad odwzorowania
fig4 = figure('Name','Filtr Kalmana i regulator typu deadbeat','units','normalized','outerposition',[0 0 0.8 0.8]);
subplot(3,1,1)
plot(z_kal_db_ff_est(2,:))
hold on
plot(z_kal_db_ff_est(1,:))
grid on
ylabel('x2_{kal db},x1_{kal db}')
xlabel('T')
legend('x2_{kal db}','x1_{kal db}')
title({'Regulator deadbeat i filtr Kalmana','Przebiegi estymowanych zmiennych stanu'})

subplot(3,1,2)
plot(Error_kal_db_ff(1,:))
hold on
plot(Error_kal_db_ff(2,:))
grid on
ylabel('Error x1_{kal db}, Error x2_{kal db}')
xlabel('T')
legend(strcat('Error x1_{kal db}=',num2str(Error_kal_db_ff_sqr(1),'%8.1e')),strcat('Error x2_{kal db}=',num2str(Error_kal_db_ff_sqr(2),'%8.1e')))
title('Bledy estymacji stanu')

subplot(3,1,3)
plot(y_kal_db_ff)
hold on
plot(u_c,'k')
grid on
ylabel('y_{kal lq},u_c')
xlabel('T')
legend(strcat('y_{kal lq}, Err=',num2str(Error_kal_db_ff_sqr(3),'%8.2e')),'u_c')
title('Wyjscie procesu')

%% Wnioski

% Zastosowanie filtru Kalmana powoduje mniejsza wariancje bledu estymacji
% stanu w porownaniu do obserwatora typu deadbeat. Filtr Kalmana jest
% rozwiazaniem zadania optymalnej estymacji stanu. Jest to estymator, który
% minimalizuje wariancje estymacji stanu, korzystajac wylacznie z
% parametrow modelu, przy zalozeniu zaklocen gaussowskich dzialajacych na
% proces (ktore otrzymano poprzez zastosowanie funkcji randn).

% Nie zauwazono roznicy w jakosci estymacji stanu podczas porownania
% wynikow symulacji otrzymanych przy wykorzystaniu regulatora 
% liniowo-kwadratowego i typu deadbeat.

% Po wprowadzeniu do systemu zewnêtrznego sygna³u sterujacego (feedforward)
% postaci: u_ff(t)=K_c*u_c(t), zaprojektowano macierz K_c w taki sposób,
% aby wyjœcie ukladu pod¹¿a³o za sygna³em u_c(t). Wyznaczaj¹c macierz K_c
% pos³u¿ono siê za³o¿eniem wyprowadzonym dla regulatora deadbeat w postaci:
% -PHI+GAMMA*K=0 (w ogolnosci A-B*K=0) dzieki czemu w rownaniu wyjscia
% otrzymano: y(t+1)=C(PHI-GAMMA*Kx)*x(t) + C*GAMMA*Kc*u_c(t) (w ogolnosci:
% y(t+1)=C(A-B*Kx)*x(t) + C*B*Kc*u_c(t)). Zatem aby wyjscie podazalo za
% wejsciem nalezalo przyrowanc C*GAMMA*Kc do 1 w wyniku czego otrzymano
% Kc=inv(C*GAMMA)*1 (w ogolnosci Kc=inv(C*B)*1). 
% Niewiele lepsze wyniki (mniejsze bledy odwzorowania wejscia przez 
% wyjscie) otrzymano dla algorytmu regulacji typu deadbeat, co jest zgodne 
% z oczekiwaniami gdyz na jego podstawie wyzaczono macierz Kc.