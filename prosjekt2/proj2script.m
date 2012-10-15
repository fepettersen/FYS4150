n = linspace(50,750,15);
rho_max = linspace(4,15,12);

for i=1:length(n)
    for rho=1:length(rho_max)
        call = sprintf('./goggen %d %f',n(i),rho_max(rho));
        system(call);
    end
end
