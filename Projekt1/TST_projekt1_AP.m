%% Teoria sterowania (TST) Projekt 1 Arkadiusz Piorkowski
%{
Cel projektu:
Celem zadania jest zbadanie i zilustrowanie dynamiki ukladu liniowego 
x(t+1)=Ax(t), x(0)=x0 
A=[a_11 a_12 a_21 a_22] 
w przestrzeni fazowej.
Wymagania:
Badania i symulacje nalezy wykonac w srodowisku Matlab lub Octave,
wykorzystuj¹c dostêpne procedury numeryczne i graficzne. Stworzony skrypt
powinien umozliwiac:
• konstruowanie macierzy o zadanym widmie, 
• konstruowanie zbioru Z punktow poczatkowych rozlozonych na okregu, 
• obliczanie trajektorii ukladu (rozwiazan r-nia stanu) dla wybranych
punktow poczatkowych,
• obliczanie wektorow wlasnych v_i macierzy A, 
• ilustracje trajektorii w przestrzeni stanow (portret fazowy), 
• ilustracje obrazu AZ zbioru Z punktow polozonych na okregu jednostkowym,
• ilustracje wektorow lambda_iv_i, lambda_i <nalezy do> sigma(A), i=1,2,
• ilustracje pola wektorowego okreslajacego trajektorie ukladu.
Zadania badawcze 
Nalezy wykonac na nastepujace zadania: 
• podac interpretacje wartosci wlasnych i wektorow wlasnych, 
• zademonstrowac zaleznosc dynamiki ukladu od widma sigma(A)={lambda
<nalezy do> C: fi_A(lambda)=0 i wektorow wlasmych macierzy A.
%}

%% Konstruowanie macierzy o zadanym widmie
clear variables;
J = [ 1 0; 0 0.5]; %widmo rzeczywiste
P = [ 3 13; 11 1];
A = P*J/P; %macierz o zadanym widmie sigma(J)
%A = [-0.4 0.7; -0.7 -0.4]; %A=[a b;-b,a] => widmo a+-ib;

%% obliczanie wektorow wlasnych v_i macierzy A
%wartosci wlasne
[L] = eig(A);
%wektory wlasne(kolumnowo) i wartosci wlasne w postaci macierzy Jordana
[V, J] = eig(A);

%% Konstruowanie zbioru Z punktów poczatkowych rozlozonych na okregu
theta=0:pi/100:2*pi;
radius=1;
z=radius*exp(1i*theta);
Z=[real(z);imag(z)];

points=0:pi/4:2*pi;
x0=radius*exp(1i*points);
X0=[real(x0);imag(x0)];

fig1=figure(1);hold on; grid on;
plot(Z(1,:),Z(2,:),'b',X0(1,:),X0(2,:),'ko','MarkerFace','r');
axis equal;
hold on;
% title('Zbior Z punktow poczatkowych rozlozonych na okregu')
% xlabel('x_1');
% ylabel('x_2');
% print(fig1,'Obrazy/zbior_punktow_poczatkowych','-dpng','-r300');

%% obliczanie trajektorii ukladu (rozwiazan r-nia stanu) dla wybranych punktow poczatkowych

iter=1000; %liczba iteracji
x1=zeros(iter,length(points));
x2=zeros(iter,length(points));

for i=1:length(points)
    for j=1:iter
        if j==1
            x1(j,i)=X0(1,i);
            x2(j,i)=X0(2,i);
        else
            x1(j,i)=A(1,1)*x1(j-1,i) + A(1,2)*x2(j-1,i);
            x2(j,i)=A(2,1)*x1(j-1,i) + A(2,2)*x2(j-1,i);
        end
    end
end

%% ilustracja trajektorii w przestrzeni stanow (portret fazowy)
plot(real(x1),real(x2),'-o');
title(['Ilustracja trajektorii w przestrzeni stanow dla [\lambda_1 , \lambda_2]=',...
    strcat('[',num2str(L(1)),',',num2str(L(2)),']')])
xlabel('x_1');
ylabel('x_2');
%print(fig1, strcat('Obrazy/portret_fazowy_',num2str(10*J(1,1)),'_',num2str(10*J(2,2)),'.png'),'-dpng','-r300');
hold off

%% ilustracja obrazu AZ zbioru Z punktow polozonych na okregu jednostkowym

fig2=figure(2);hold on; grid on;
power=5;
for t=0:power
    AtZ=A^(t)*Z;
    plot(real(AtZ(1,:)),real(AtZ(2,:)));
    hold on;   
end
axis equal;

%% ilustracja wektorow lambda_i*v_i, lambda_i <nalezy do> sigma(A), i=1,2,

vect = zeros(length(L),length(L));

for i = 1:length(L)
    vect(:,i) = L(i) * V(:,i);
end

quiver(zeros(1,length(L)),zeros(1,length(L)),real(vect(1,:)),real(vect(2,:)),'AutoScale','off');

axis equal;
hold on

%% ilustracja pola wektorowego okreslajacego trajektorie ukladu

x1_lim=xlim;
x2_lim=ylim;

[X_1,X_2] = meshgrid(x1_lim(1):(x1_lim(2)-x1_lim(1))/15:x1_lim(2),x2_lim(1):(x2_lim(2)-x2_lim(1))/15:x2_lim(2));

Ident=eye(size(A,1),size(A,2));

for i=1:size(X_1,2)
    for j=1:size(X_2,1)
        dx=(A-Ident)*[X_1(i,j);X_2(i,j)];
        px1(i,j)=dx(1);
        px2(i,j)=dx(2);
    end
end
quiver(X_1,X_2,real(px1),real(px2),'Autoscale','on');

legend('Z','AZ','A^2Z','A^3Z','A^4Z','A^5Z','\lambda_i\nu_i','pole wekt.','Location','northeastoutside');
title({['Ilustracja wektorow \lambda_i\nu_i, obrazu A^tZ i pola wektorowego'],...
    ['dla [\lambda_1 , \lambda_2]=',strcat('[',num2str(L(1)),',',num2str(L(2)),']')]})
xlabel('x_1');
ylabel('x_2');
%print(fig2, strcat('Obrazy/pole_wektorowe_',num2str(10*J(1,1)),'_',num2str(10*J(2,2)),'.png'),'-dpng','-r300');
