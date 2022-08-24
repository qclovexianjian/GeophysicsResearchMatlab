function [ tdelay, tori ] = cqwd2td( wd )
% estimate time delay in wd matrix by linear fitting

[~,tdelay] = max(wd);
tdelay = tdelay - 1;
tdelay = [0,tdelay];
tori = 0:size(wd,2);

[xData, yData] = prepareCurveData( tori, tdelay );

% Set up fittype and options.
ft = fittype( '0*(x<=a) + (k*x-k*a).*(x>a&x<(a+h))+(k*h)*(x>=(a+h));', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.600197393270656 0.172017444134905 0.0164175672011078];

% Fit model to data.
[fitresult] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'Linear_Delay' );
h = plot( fitresult, xData, yData );
legend( h, 'tdelay vs. tori', 'Linear_Delay', 'Location', 'NorthEast' );
% Label axes
xlabel tori
ylabel tdelay
grid on


end

