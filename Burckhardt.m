function mu_ijx = Burckhardt(sijx, c1, c2, c3)
%BURCKHARDT Compute the tire force based on a burckhardt tire model

mu_ijx = c1*(1-exp(-c2*sijx)) - c3*sijx;
end

