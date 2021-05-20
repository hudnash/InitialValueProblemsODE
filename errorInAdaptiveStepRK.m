clear all
clf
% dydt = @(t,y)(1-2*t*y)/(1+t^2); 
% RK4 and RK5 are embedded RK estimators of the y(t) function.
hold on
[tseries, yseries, abs_err] = adaptiveStep(@RK4,@RK5,.0001,.25,0,0,10); % Lower tolerance -> frequent steps
plot(tseries, yseries,'xr');
plot(tseries, abs_err,'-r');
[tseries, yseries, abs_err] = adaptiveStep(@RK4,@RK5,.01,.25,0,0,10); % Higher tolerance -> fewer steps
plot(tseries, yseries,'ob');
plot(tseries, abs_err,'--b');
legend('TOL = 0.0001: y(t)','TOL = 0.0001: ERROR','TOL = 0.01: y(t)','TOL = 0.01: ERROR');

% This version of the IVP solution function differs from the demonstration.
% It allows for the initial condition to be anywhere along the function,
% not exclusively isolated at the start of the range of t.

function [t, y, true_errors] = adaptiveStep(func1,func2,TOL,y0,t0,t_start,t_end) 
h = TOL; t(1) = t_start; % y4(1) = y0; y5(1) = y0; y(1) = y0; t(1) = t0; % These declarations rely on t0 <= t_start.
initialValue = 0;
if t_start >= t0
    y4(1) = y0; y5(1) = y0; y(1) = y0; t(1) = t0;
    initialValue = 1;
end
i = 1; % We cannot use t as the index variable because it is not an integer.
while t(i) < t_end - h % Term "-h" prevents loop from exceeding desired t range.
    if initialValue == 0 && t(i) >= t0
        y4(i) = y0; y5(i) = y0; y(i) = y0; t(i) = t0;
        initialValue = 1; % Initial condition has been applied.
                          % CAVEAT: Step size computation is ignored and
                          % restarted from this new interrupting point.
    end
    y4(i+1) = func1(y(i),t(i),h); y5(i+1) = func2(y(i),t(i),h);
    err = y5(i) - y4(i);
    if err ~= 0 % Counteract error case that RK4 == RK5.
        h = h*abs(TOL/err)^.2;
    end
    t(i+1) = t(i) + h;
    y(i+1) = func2(y(i),t(i),h);
    yAnalytical(i+1) = (0.25 + t(i+1))/(1 + t(i+1)^2);
    % calculate the true error
    true_errors(i+1) = abs(y(i+1)-yAnalytical(i+1));
    i=i+1;
end
end


