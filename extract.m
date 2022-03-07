function [B, I, D, T_ON, T_OFF, leakiness, maxl, lmaxt] =  extract(means, time, frames_base, frames_on, frames_off, frames_tail, x0)
%this function extracts parameters for induction kinetics of single cells

B = []; %array containing parameter b for single cell data
I = []; %array containing parameter i for single cell data
D = []; %array containing parameter d for single cell data
T_ON = []; %array containing parameter t-on for single cell data
T_OFF = []; %array containing parameter t-off for single cell data
leakiness = []; %array leakiness data for single cell data
maxl = []; %array containing max level of expression for single cell data
lmaxt = []; %array containing level of expression at t = 3.5 h for single cell data

N_rows = size(means,1);


for index = 1:N_rows
    
    %we first extract the ON induction with cells that are present from t =
    %1:max(frames_on)
    if(~isnan(sum(means(index, 1:max(frames_on)))))
        
        [b, i, d, t_on, t_off] = extract_cell(means(index,:), time, frames_base, frames_on, frames_off, frames_tail, x0, 0, 1, 0);
        
        B = [B b];
        I = [I i];
        T_ON = [T_ON t_on];
        
        leakiness = [leakiness mean(means(index, frames_base))];
        
    end
    
    %we then extract the parameters for OFF induction with cells 
    %that are present from t = frames_off(1):max(frames_tail)
    %in this case we don't have the value for the basal expression before
    %the induction for all cells so we use mean(B) for all; this does not make a difference for degradation rate
    
    if(~isnan(sum(means(index, frames_off(1):max(frames_tail)))))
        [b, i, d, t_on, t_off] = extract_cell(means(index,:), time, frames_base, frames_on, frames_off, frames_tail, x0, 0, 0, mean(B));
        D = [D d];
        T_OFF = [T_OFF t_off];
        maxl = [maxl max(means(index, frames_on(1):max(frames_tail)))];
    end
    
    %at the end we extract lmaxt for all cells present at t = frames_off(1)
    if(~isnan(means(index, frames_off(1))))
        lmaxt = [lmaxt means(index, frames_off(1))];
    end
    
end