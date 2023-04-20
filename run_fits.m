function output = run_fits(promoter_data, time, frames_on, frames_base_low, frames_off, frames_tail, x0, plotting, extract_ON, base_for_OFF, bound_time)
% for 'on', will get two parameters:  i * f / 2 and t-on.  Will not get low  baseline: f * b/d/(d+f)
% for 'off', will get two parameters: i * f / 2 and t-off. Will not get high baseline: f * (b + i)/d/(d+f)
% i can be different between 'on' and 'off'
% for 'tail', will get -d

%we can either extract parameters for both on and off (extract_ON = 1)
%or just for off (extract_ON = 0)

%this separation is needed since some cells appear during the induction period
%For them, we don't know wha was the leakiness. Therefore, for extracting
%off parameters and tail, we forward average basal acivity of cells present
%at T = 0 (forwarded to function as base_for_OFF)
output = [];

t_on = time(frames_on);
t_off = time(frames_off) - time(frames_off(1));
t_tail = time(frames_tail) - time(frames_tail(1));

opts = optimset('MaxIter',10^6,'TolFun',10^(-6),'TolX',10^(-6));


%if extract_ON == 1 we are extracting ON and OFF kinetics
if(extract_ON == 1)
    %output is [i*f/2 t-on i*f/2 t-off -d] vector
    y_on = promoter_data(frames_on);
    base_low = mean(y_on(frames_base_low));
    fun = @(x) sseval(x,t_on,y_on,base_low);
    bestx_on = fminsearchbnd(fun,x0,[0 bound_time],[],opts);
    output = [output; bestx_on];
    
else %we extract only OFF kinetics
     %output is [NaN NaN i*f/2 t-off -d] vector
    output = [output; NaN; NaN]; 
    base_low = base_for_OFF;
end

y_off = promoter_data(frames_off);
base_high = mean(y_off(1));
fun = @(x) sseval(x,t_off,y_off,base_high);
bestx_off = fminsearchbnd(fun,x0,[-Inf bound_time],[],opts);
output = [output; bestx_off];

%fitting tail of fluorescent protein degradation as an exponential decay
y_tail = promoter_data(frames_tail);
f = fit(t_tail', (y_tail-base_low)', 'exp1');
output = [output; f.b];

if(plotting == 1)
    figure()
    plot(time, promoter_data, 'LineWidth', 9, 'Color', [0 0 0])
    hold all
    if(extract_ON == 1)
    %ploting fit for ON switch during t_on
    plot(t_on, fitfunsquare(bestx_on,t_on,base_low), 'LineWidth', 5)
    end
    %ploting fit for tau_off during frames_off
    plot(t_off + time(frames_off(1)), fitfunsquare(bestx_off,t_off,base_high), 'LineWidth', 5)
    %ploting fit for fluorescent protein degradation during frames_tail
    plot(t_tail + time(frames_tail(1)), base_low + f(t_tail), 'LineWidth', 5);
end

