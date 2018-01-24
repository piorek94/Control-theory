%Arkadiusz Piórkowski 259079- praca domowa TST 17Z
%animacja rysowania okrêgu(wielomian/wielomian)


clear
% wybor rodzaju transformacji: k==1 wielomian, k==2 funkcja wymierna
k=1;

%--------------------------------------------------------------------------
%okrêg zespolony o promieniu r
theta=0:pi/100:2*pi;
[m,n]=size(theta);
radius=2;
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
elseif k==2
    P=(nom(1)*circle.^2+nom(2)*circle+nom(3))./(den(1)*circle+den(2));
end

%--------------------------------------------------------------------------
%nowy film
% v=VideoWriter('point_by_point_polynominal.avi');
% open(v);

fig=figure('NumberTitle','Off','MenuBar','None','Name','Animacja transformacji P okrêgu zespolonego','position',[200 200 1100 500]);  

for d = 1:n
    %pierwszy wykres-okr¹g
    subplot(1,2,1)
    
    %rysowanie pierwiastków i biegunów funkcji przekszta³caj¹cej
    plot(real(root),imag(root),'bo','MarkerSize',5);
    hold on
    if k==2
        plot(real(pole),imag(pole),'k*','MarkerSize',5);
        hold on
    end
    
    %rysowanie okrêgu
    plot(real(circle(1:d)),imag(circle(1:d)),'r');
    hold on
    
    %rysowanie wektorów
    plot([real(circle(d)) real(root(1))],[imag(circle(d)) imag(root(1))],'b')
    hold on
    plot([real(circle(d)) real(root(2))],[imag(circle(d)) imag(root(2))],'b')
    hold on
    if k==1
        plot([real(circle(d)) real(root(3))],[imag(circle(d)) imag(root(3))],'b')
    elseif k==2
        plot([real(circle(d)) real(pole(1))],[imag(circle(d)) imag(pole(1))],'k')
    end
    hold off
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
    xlabel('Re(z)');
    ylabel('Im(z)');
    if k==1
       legend('pierwiastki wielomianu','okr¹g','Location','SouthWest') 
    elseif k==2
       legend('pierwiastki f.wymiernej','bieguny f.wymiernej','okr¹g','Location','SouthWest') 
    end
    
    %----------------------------------------------------------------------
    %drugi wykres-transformacja
    subplot(1,2,2)
    
    %rysowanie pierwiastków i biegunów funkcji przekszta³caj¹cej
    plot(real(root),imag(root),'bo','MarkerSize',5);
    hold on
    if k==2
        plot(real(pole),imag(pole),'k*','MarkerSize',5);
        hold on
    end
    
    %rysowanie transformacji
    plot(real(P(1:d)),imag(P(1:d)),'r')
    hold on
    
    %rysowanie wektorów
    plot([real(P(d)) real(root(1))],[imag(P(d)) imag(root(1))],'b')
    hold on
    plot([real(P(d)) real(root(2))],[imag(P(d)) imag(root(2))],'b')
    hold on
    if k==1
        plot([real(P(d)) real(root(3))],[imag(P(d)) imag(root(3))],'b')
    elseif k==2
        plot([real(P(d)) real(pole(1))],[imag(P(d)) imag(pole(1))],'k')
    end
    hold off
    grid ON
    axis equal
    
    %opis wykresu drugiego
    if k==1
        title('Transformacja-wielomian');
    elseif k==2
        title('Transformacja-f.wymierna');
    end
    xlabel('Re(z)');
    ylabel('Im(z)');        
    M(d)=getframe(fig);

    %zapis filmu
    %writeVideo(v,M(z));
end

%zamkniecie pliku
% close(v);
% 
% close all