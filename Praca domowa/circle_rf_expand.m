%Arkadiusz Piórkowski 259079- praca domowa TST 17Z
%animacja powiêkszania okrêgu(funckja wymierna)
clear

%okrêg zespolony o promieniu r
theta=0:pi/600:2*pi;
radius=0.5:0.02:4;
[p,q]=size(radius);
[m,n]=size(theta);

%pierwiastki i bieguny funckji wymiernej
roots_polyn=[-2 + 1i, 1 - 1i];
pole_polyn=[2+3*1i];
%mianownik i licznik funkcji wymiernej
polyn=poly(roots_polyn);
polyn_2=poly(pole_polyn);

%nowy film
v=VideoWriter('circle_rf_expand.avi');
open(v);

fig=figure('NumberTitle','Off','MenuBar','None','Name','Wplyw rozszerzania okrêgu zespolonoego na tranformacjê-f.wymierna','position',[200 200 1100 500]);

for k = 1:q
    circle=radius(k)*exp(theta*1i);
    P=(polyn(1)*circle.^2+polyn(2)*circle+polyn(3))./(polyn_2(1)*circle+polyn_2(2));   
        marker='r-';
        subplot(1,2,1)
        plot(real(roots_polyn),imag(roots_polyn),'bo','MarkerSize',5);
        hold on
        plot(real(pole_polyn),imag(pole_polyn),'go','MarkerSize',5);
        hold on
        plot(circle,marker,'MarkerSize',2)
        grid ON
        axis equal
        hold off
        legend('zera funckji wymiernej','bieguny funkcji wymiernej','okr¹g','Location','SouthWest')
        title('Okrag');
        xlabel('Re(z)');
        ylabel('Im(z)'); 
        
        subplot(1,2,2)
        plot(P,marker,'MarkerSize',2)
        grid ON
        axis equal
        hold off
        title('Transformacja-funkcja wymierna');
        xlabel('Re(z)');
        ylabel('Im(z)');        
        M(k)=getframe(fig);
        
        %zapis filmu
        writeVideo(v,M(k));
end

%zamkniecie pliku
close(v);

close all