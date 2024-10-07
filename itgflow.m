%% Code developed for the paper "Assessment of the Rate and Primary Pathway of Glymphatic Flow Through the Human Brain" submitted to Science Translational Medicine
%
%  Code-authors: Anders Wåhlin, Viktor Vigren Näslund, Anders Eklund
%
%  gFlow is the function for estimating k1 and k2 from which glymphatic flow
%  rate and ISF volume fraction can be calculated.

%  getEstCt is a support function that calculates brain tissue concentrations given
%  a matrix of CSF concentrations (Nrois X time), as well as a given k1 and
%  k2.
%
%  simulatedData.mat file contains simulated concentrations (q=55, ve=.25) for a
%  set of CSF (concCSF) and tissue ROIs (concTissue, in mMol), a vector
%  with time points for the interpolated curves (it, in hours), the time
%  step of the interpolation (dt, in hours), and the temporal position of the
%  original MRI measurments (points).
%
%% Main section
load simulatedData.mat

k0 = zeros(1,2)

fun1=@(k)gFlow(k,it,concCSF,concTissue,dt,points); %fit the model

[kest,fval,exitflag,output] = fminunc(fun1,k0);

q = kest(1)*tissueVolume %Glymphatic flow rate

ve = kest(1)/kest(2) %ISF volume fraction

%% Support functions

function estCt = getEstCt(k,it,ccsf,dt)

% Estimates brain tissue concentrations provided a CSF concentrations
% (Nrois X time), k1 and k2.
% dt is the temporal resolution in hrs

for i = 1:size(ccsf,1) %loop over ROIs
    estCt(i,:) = k(1)*convolution(ccsf(i,:)',exp(-k(2).*it'))*dt ;
end

end

%%
function sse = gFlow(k,it,ccsf,ct,dt,points)

% Calculates the difference between measured (ct) and estimated brain tissue
% concentrations provided CSF concentrations (ccsf), k1 and k2 (stored as a
% vector k).
% The vector points specifies time points where the difference is calculated, and dt is the temporal resolution in hrs

estCt = getEstCt(k,it,ccsf,dt);
sse = sum(sum(((ct(:,points)-estCt(:,points)).^2)));

end

%%
function c = convolution(a, b)
c = conv2(a(:), b(:), 'full');
c = c(1:size(a, 1));

if isrow(a)
    c = c.';
end
end
