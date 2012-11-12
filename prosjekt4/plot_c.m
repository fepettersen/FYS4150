infile = textread('isingresults_n0_temp_1.0000.txt');

E = infile(:,1);
M = infile(:,2);
flips = infile(:,3);
N = 1:length(E);

% plot(N,E)
% xlabel(' number of Monte Carlo cycles')
% ylabel(' average energy ')
% figure()
% plot(N,M)
% xlabel(' number of Monte Carlo cycles')
% ylabel(' average magnetization')
% figure()
% plot(N,flips,'ro')
% xlabel(' number of Monte Carlo cycles')
% ylabel(' number of accepted flips')

energy = -798.82;
eps = 1e-3;
counter = 0;

for i = 0.1*length(E):length(E)
   if (abs(E(i)-energy)<=eps)
       counter =counter +1
   end
end

counter