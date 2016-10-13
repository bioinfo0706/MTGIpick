
function S = EFC_versus_iv(G2,G2r,O,Rr,pth,th,maxEFC,foldchange_or_additive)

if nargin<8;
    foldchange_or_additive = 1;
end
if nargin<7;
    maxEFC = 10;
end
if nargin<6;
    th = 1e-5;
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
Ip = abs(Cmp-Zs)<th;
In = abs(Cmn-Zs)<th;
psn(In&Ip)=real(psn(In&Ip)); %one solution, remove imaginary part (numerical stability)Pse(Ip) = psp(Ip);
Pse(In) = psn(In);
if foldchange_or_additive==1;
    E = (Pse.*G2)./Rr;
    E(E<1|isnan(E))=1;
else
    E = (Pse.*G2)-Rr;
    E = E./G2;
    E(E<0|isnan(E))=0;
end

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
Ip = abs(Cmp-Zs)<th;
In = abs(Cmn-Zs)<th;
psn(In&Ip)=real(psn(In&Ip)); %one solution, remove imaginary part (numerical stability)
Psd(Ip) = psp(Ip);
Psd(In) = psn(In);
if foldchange_or_additive==1;
    D = (Psd.*G2)./Rr;
    D(D>1|isnan(D))=1;
else
    D = (Psd.*G2)-Rr;
    D = D./G2;
    D(D>0|isnan(D))=0;
end


%% Score
if foldchange_or_additive==1;
    S = ones(size(E));
    S(E>1) = E(E>1);
    S(D<1) = D(D<1);
    S = log2(S);
    S(S>maxEFC) = maxEFC;
    S(S<-maxEFC) = -maxEFC;
else
    S = zeros(size(E));
    S(E>0) = E(E>0);
    S(D<0) = D(D<0);
end
end  