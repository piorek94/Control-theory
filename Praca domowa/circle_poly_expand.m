%Arkadiusz Piórkowski 259079- praca domowa TST 17Z
%animacja powiêkszania okrêgu(wielomian)
clear

%okrêg zespolony o promieniu r
theta=0:pi/150:2*pi;
radius=0.5:0.02:5;
[p,q]=size(radius);
[m,n]=size(theta);

%pierwiastki wielomianu 3 stopnia
roots_polyn=[-2 + 1i, 1 - 1i, 2+3*1i];
%wielomian zespolony
polyn=poly(roots_polyn);

%nowy film
v=VideoWriter('circle_poly_expand.avi');
open(v);

fig=figure('NumberTitle','Off','MenuBar','None','Name','Wplyw rozszerzania okrêgu zespolonoego na tranformacjê-wielomian','position',[200 200 1100 500]);

for k = 1:q
    circle=radius(k)*exp(theta*1i);
    P=polyn(1)*circle.^3+polyn(2)*circle.^2+polyn(3)*circle+polyn(4);   
        marker='r-';
        subplot(1,2,1)
        plot(real(roots_polyn),imag(roots_polyn),'bo','MarkerSize',5);
        hold on
        plot(circle,marker,'MarkerSize',2)
        grid ON
        axis equal
        hold off
        legend('zera wielomianu','okr¹g','Location','SouthWest')
        title('Okrag');
        xlabel('Re(z)');
        ylabel('Im(z)'); 
        
        subplot(1,2,2)
        plot(P,marker,'MarkerSize',2)
        grid ON
        axis equal
        hold off
        title('Transformacja-wielomian');
        xlabel('Re(z)');
        ylabel('Im(z)');        
        M(k)=getframe(fig);
        
        %zapis filmu
        writeVideo(v,M(k));
end

%zamkniecie pliku
close(v);

close all