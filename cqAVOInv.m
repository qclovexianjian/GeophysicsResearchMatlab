function [ InvResult ] = cqAVOInv( cmp, dt, window, inc_ang, method, ifplot  )
% Uses different AVO inversion methods to do AVO inversion
%
% input:
% cmp = common-midpoint gather (must be aligned to the window center)
% dt = sampling rate in seconds
% window = [center, width]  AVO inversion will be done in this window and incident angle
%          will correspond to the window center line
% inc_ang = incident angles for each tace with reflection event at the window center
%           (in degree)
% method = 'Shuey'(default) Shuey's method AVO equation
%               rc(theta) = A + B*sin(theta).^2
%               Here we only use incident angle as theta, rather than averaged incident
%               and transmitted angle
% ifplot = true (default) to plot result
%
% Output
% InvResult = structure of the output
%   For method == 'Shuey'
%       A = Intercept or Normal incident reflection coefficient
%       B = Slope in front of sin(theta)
%       A_win = Portion of A in the window
%       B_win = Portion of B in the window

center = window(1); width = window(2);

t = 0:dt:(size(cmp,1)-1)*dt;
cmp_win = cmp( t<=(center+width/2) & t>=(center-width/2), : ); % only take that portion of the window

if strcmp(method,'Shuey')
    A = zeros(size(cmp,1),1);
    B = A;
    A_win = zeros(size(cmp_win,1),1);
    B_win = A_win; % Initialization of A and B in the window
    for n = 1:length(A_win)
        x = sin(inc_ang/180*pi).^2;
        y = cmp_win(n,:);
        p = polyfit(x,y,1);
        A_win(n) = p(2);
        B_win(n) = p(1);
    end
    A( t<=(center+width/2) & t>=(center-width/2) ) = A_win;
    B(  t<=(center+width/2) & t>=(center-width/2) ) = B_win;
    InvResult.A = A;
    InvResult.B = B;
    InvResult.A_win = A_win;
    InvResult.B_win = B_win;
    
    if ifplot
        figure; title(['Cross Plot: Shuey, Window Center: ',num2str(center),'(s)']);
        plot(A_win, B_win, 'rx','markersize',10);axis equal;
        xlabel('A'); ylabel('B');
        
        cqwva(cmp, dt, inc_ang);title('Input CMP Gather(Aligned)')
        ylim([center-width/2,center+width/2])
        xlabel('Incident Angle(Degree)'); ylabel('Time(s)');
        set(gcf,'pos',[2    42   791   319]);
        
        cqwva([A,B],dt)
        ylabel('Time(s)');set(gca,'xtick',[]);xlabel('A   |   B')
        ylim([center-width/2,center+width/2])
        set(gcf,'pos',[1373          42         250         319])
    end
end


end

