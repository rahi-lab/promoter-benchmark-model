function [b, i, d, ton, toff, maxl_a] = extract_cell(expression, time, frames_base, frames_on, frames_off, frames_tail, x0, plotting, extract_ON, base_for_OFF, bound_time)
%this function extracts parameters of induction based on the expansion of
%model response into second order Taylor series

output = run_fits(expression, time, frames_on, frames_base, frames_off, frames_tail, x0, plotting, extract_ON, base_for_OFF, bound_time);

%output is [i*f/2 t-on i*f/2 t-off -d]
f = log(2)/16.54; %maturation half-time of yEVenus as measured by the cyclohexamide block
i = output(1)*2/f;
d = -output(5);
ton = output(2);
toff = output(4);

base_low = mean(expression(frames_base)); %basal level of expression

%base_low equals f*b/(d+f)/d (from Taylor series);
%to get b, we use the d extracted from tail
b = base_low*d*(d+f)/f;

maxl_a = max(expression);


