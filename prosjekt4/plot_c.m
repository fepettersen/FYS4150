infile1 = textread('isingresults_n0_temp_1.0000.txt');
infile2 = textread('isingresults_n0_temp_2.4000.txt');
E = infile1(:,1);
M = infile1(:,2);
flips1 = sum(infile1(:,3));
flips2 = sum(infile2(:,3));
N = 1:length(E);
flips = [flips1/length(flips1),flips2/length(flips2)];
temp = [1.0,2.4];
% plot(N,E)
% xlabel(' number of Monte Carlo cycles')
% ylabel(' average energy ')
% figure()
% plot(N,M)
% xlabel(' number of Monte Carlo cycles')
% ylabel(' average magnetization')
% figure()
plot(temp,flips','r')
xlabel(' temperature')
ylabel(' number of accepted flips')

% energy = -798.82;
% eps = 1e-3;
% counter = 0;
% 
% for i = 0.1*length(E):length(E)
%    if (abs(E(i)-energy)<=eps)
%        counter =counter +1
%    end
% end
% 
% counter