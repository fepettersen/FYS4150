infile = textread('other_results.txt');

E = infile(:,1);
M = infile(:,2);
flips = infile(:,3);
N = 1:length(E);
plot(N,E)
xlabel(' number of Monte Carlo cycles')
ylabel(' average energy ')
figure()
plot(N,M)
xlabel(' number of Monte Carlo cycles')
ylabel(' average magnetization')
figure()
plot(N,flips,'ro')
xlabel(' number of Monte Carlo cycles')
ylabel(' number of accepted flips')