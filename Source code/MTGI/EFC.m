
function S = EFC(T,G1,G2,O,pth,th);

if nargin<6;
    th = 1e-5;
end
if nargin<5;
    pth = 1e-6;
end

%% Enrichment
Zs = -norminv(pth);
a = (Zs.^2).*G2+(G2.^2);
b = -(2.*O.*G2 + (Zs.^2).*G2);
c = (O.^2);
psp = (-b + sqrt((b.^2)-(4.*a.*c)))./(2.*a);
psn = (-b - sqrt((b.^2)-(4.*a.*c)))./(2.*a);
Pse = NaN*ones(size(psp));
Cmp = ((O-(psp.*G2))./sqrt(G2r.*psp.*(1-psp)));
Cmn = ((O-(psn.*G2))./sqrt(G2r.*psn.*(1-psn)));
Ip = abs(Cmp-Zs)<th;
In = abs(Cmn-Zs)<th;
psn(In&Ip)=real(psn(In&Ip)); %one solution, remove imaginary part (numerical stability)
Pse(Ip) = psp(Ip);
Pse(In) = psn(In);
E = (G2.*Pse)./((G2.*G1)./T);
E(E<1|isnan(E))=1;


%% Depletion
Zs = norminv(pth);
a = (Zs.^2).*G2+(G2.^2);
b = -(2.*O.*G2 + (Zs.^2).*G2);
c = (O.^2);
psp = (-b + sqrt((b.^2)-(4.*a.*c)))./(2.*a);
psn = (-b - sqrt((b.^2)-(4.*a.*c)))./(2.*a);
Psd = NaN*ones(size(psp));
Cmp = ((O-(psp.*G2))./sqrt(G2r.*psp.*(1-psp)));
Cmn = ((O-(psn.*G2))./sqrt(G2r.*psn.*(1-psn)));
Ip = abs(Cmp-Zs)<th;
In = abs(Cmn-Zs)<th;
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