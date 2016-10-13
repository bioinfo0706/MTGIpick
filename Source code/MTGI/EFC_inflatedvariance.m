function S = EFC_inflatedvariance(T,G1,G2,G2r,O,pth,th)

if nargin<6;
    th = 1e-5; %for numerical stability
end
if nargin<5;
    pth = 1e-6;
end

%% Enrichment
Zs = -norminv(pth);
aa = (Zs.^2).*G2r+(G2.^2);
bb = -(2.*O.*G2 + (Zs.^2).*G2r);
cc = (O.^2);
psp = (-bb + sqrt((bb.^2)-(4.*aa.*cc)))./(2.*aa);
psn = (-bb - sqrt((bb.^2)-(4.*aa.*cc)))./(2.*aa);
Pse = NaN*ones(size(psp));
Cmp = ((O-(psp.*G2))./sqrt(G2r.*psp.*(1-psp)));
Cmn = ((O-(psn.*G2))./sqrt(G2r.*psn.*(1-psn)));
Ip = abs(Cmp-Zs)<(O*th); 
In = abs(Cmn-Zs)<(O*th);
psn(In&Ip)=real(psn(In&Ip)); %one solution, remove imaginary part (numerical stability)
Pse(Ip) = psp(Ip);
Pse(In) = psn(In);
E = (G2.*Pse)./((G2.*G1)./T);
E(E<1|isnan(E))=1;

%% Depletion
Zs = norminv(pth);
aa = (Zs.^2).*G2r+(G2.^2);
bb = -(2.*O.*G2 + (Zs.^2).*G2r);
cc = (O.^2);
psp = (-bb + sqrt((bb.^2)-(4.*aa.*cc)))./(2.*aa);
psn = (-bb - sqrt((bb.^2)-(4.*aa.*cc)))./(2.*aa);
Psd = NaN*ones(size(psp));
Cmp = ((O-(psp.*G2))./sqrt(G2r.*psp.*(1-psp)));
Cmn = ((O-(psn.*G2))./sqrt(G2r.*psn.*(1-psn)));
Ip = abs(Cmp-Zs)<(O*th); 
In = abs(Cmn-Zs)<(O*th);
psn(In&Ip)=real(psn(In&Ip)); %one solution, remove imaginary part (numerical stability)
Psd(Ip) = psp(Ip);
Psd(In) = psn(In);
D = (G2.*Psd)./((G2.*G1)./T);
D(D>1|isnan(D))=1;

%% Score =
S = ones(size(E));
S(E>1) = E(E>1);
S(D<1) = D(D<1);
S = log2(S);
end