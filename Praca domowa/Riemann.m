%Arkadiusz Piórkowski 259079- praca domowa TST 17Z
%przekszta³cenie steroidalne transformacji(wielomian/funckja wymierna)
%okregu zespolonego
%wymagany toolbox geom3d
clear
% wybor rodzaju transformacji: k==1 wielomian, k==2 funkcja wymierna
k=1;

%--------------------------------------------------------------------------
%okrêg zespolony o promieniu r
theta=0:pi/100:2*pi;
[m,n]=size(theta);
radius=sqrt(14);
circle=radius*exp(theta*1i);

%--------------------------------------------------------------------------
%wybór wprowadzenia wielomianu/funkcji wymiernej: p==1 pierwiastki/bieguny
%p==2 wspó³czynniki
p=1;
if k==1
    if p==1
        %pierwiastki wielomianu 3 stopnia
        root=[-2 + 1i, 1 - 1i, 2+3*1i];
        %wielomian zespolony
        polyn=poly(root); 
    elseif p==2
        %wspó³czynniki wielomianu 3 stopnia
        polyn=[1 1 1 1];
        %pierwiatki wielomianu
        root=roots(polyn);
    end
elseif k==2
    if p==1
        %pierwiastki funkcji wymiernej
        root=[-2 + 1i, 1 - 1i];
        %bieguny funkcji wymiernej
        pole=[2+3*1i];
        %licznik funkcji wymiernej
        nom=poly(root);
        %mianownik funkcji wymiernej
        den=poly(pole);
    elseif p==2
        %licznik funkcji wymiernej
        nom=[1 1 1];
        %mianownik funkcji wymiernej
        den=[1 -1i];
        %pierwiastki funkcji wymiernej
        root=roots(nom);
        %bieguny funkcji wymiernej
        pole=roots(den);
    end
end

%--------------------------------------------------------------------------
%obliczenie transformacji
if k==1
    P=polyn(1)*circle.^3+polyn(2)*circle.^2+polyn(3)*circle+polyn(4);
    zero=[0 0 0];
elseif k==2
    P=(nom(1)*circle.^2+nom(2)*circle+nom(3))./(den(1)*circle+den(2));
    zero=[0 0];
end

%--------------------------------------------------------------------------
%plaszczyzna Reimanna
sphere=[0 0 0.5 0.5]; %[Xc Yc Zc R];
%stworzenie lini przecinaj¹cych sfere
lines=zeros(n,6);
P_sphere=[sphere(1),sphere(2),sphere(3)+sphere(4)];
%punkty przeciêcia linii i sfery
Intersect_Points=zeros(n,3);
Point=zeros(2,3);
for d=1:n
    P_plane=[real(P(d)),imag(P(d)),0];
    lines(d,:)=createLine3d(P_plane,P_sphere);
    Point=intersectLineSphere(lines(d,:),sphere);
    Intersect_Points(d,:)=Point(1,:);
end


%--------------------------------------------------------------------------
%wykresy

fig=figure('Name','Transformacja okrêgu zespolonego i przekszta³cenie sterograficzne');  

%pierwszy wykres-okr¹g
subplot(2,2,1)

%rysowanie pierwiastków i biegunów funkcji przekszta³caj¹cej
plot(real(root),imag(root),'bo','MarkerSize',5);
hold on
if k==2
    plot(real(pole),imag(pole),'k*','MarkerSize',5);
    hold on
end

%rysowanie okrêgu
plot(real(circle),imag(circle),'r');
grid ON

%opis wykresu pierwszego
title('Okr¹g');
if radius>4
    xlim([-radius radius]);
    ylim([-radius radius]);
else
    xlim([-4 4]);
    ylim([-4 4]);
end
axis equal
xlabel('Re(z)');
ylabel('Im(z)');
if k==1
    legend('pierwiastki wielomianu','okr¹g','Location','SouthWest')
elseif k==2
    legend('pierwiastki f.wymiernej','bieguny f.wymiernej','okr¹g','Location','SouthWest')
end

%--------------------------------------------------------------------------
%drugi wykres-transformacja
subplot(2,2,2)

%rysowanie pierwiastków i biegunów funkcji przekszta³caj¹cej
plot(real(root),imag(root),'bo','MarkerSize',5);
hold on
if k==2
    plot(real(pole),imag(pole),'k*','MarkerSize',5);
    hold on
end

%rysowanie transformacji
plot(real(P),imag(P),'r')
grid ON
axis equal
x_lim=xlim;
y_lim=ylim;
%opis wykresu drugiego
if k==1
    title('Transformacja-wielomian');
elseif k==2
    title('Transformacja-f.wymierna');
end
xlabel('Re(z)');
ylabel('Im(z)');

%--------------------------------------------------------------------------
%trzeci wykres-przeciecie prostych z okregiem
subplot(2,2,3)
%rysowanie pierwiastków i biegunów funkcji przekszta³caj¹cej
plot3(real(root),imag(root),zero,'bo','MarkerSize',5);
hold on
if k==2
    plot3(real(pole),imag(pole),0,'k*','MarkerSize',5);
    hold on
end
Z=zeros(1,d);
%rysowanie transformacji
plot3(real(P),imag(P),Z,'r')
drawSphere(sphere, 'nPhi', 360, 'nTheta', 180);
l = light; 
hold on
for w=1:n
    drawLine3d(lines(w,:),'b');
end
axis([x_lim(1) x_lim(2) y_lim(1) y_lim(2) sphere(3)-sphere(4) sphere(3)+sphere(4)])
grid ON
view(3)
xlabel('Re(z)');
ylabel('Im(z)');
title('Proste przecinaj¹ce sfere');

%--------------------------------------------------------------------------
%czwarty wykres-p³aszczyzna Reimanna
subplot(2,2,4)
drawSphere(sphere, 'nPhi', 360, 'nTheta', 180);
drawSphere(sphere, 'linestyle', ':', 'facecolor', 'none');
axis([sphere(1)-sphere(4) sphere(1)+sphere(4) sphere(2)-sphere(4) sphere(2)+sphere(4) sphere(3)-sphere(4) sphere(3)+sphere(4)])
axis equal
l = light;
hold on
plot3(Intersect_Points(:,1),Intersect_Points(:,2),Intersect_Points(:,3),'LineWidth',3,'Color','r');

title('Transformacja-p³aszczyzna Reimanna');
xlabel('Re(z)');
ylabel('Im(z)');