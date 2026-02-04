function [St_b, PSD_b] = PSD_binning_scr(St,PSD,St_scr,dsts)
%protection zone around the fundamental and first harmonic
%fundamental
%% 
dP = 0.03;
p1_str = find(St > St_scr-dP,1);
p1_end = find(St > St_scr+dP,1);
%first harmonic
p2_str = find(St > 2*St_scr-dP,1);
p2_end = find(St > 2*St_scr+dP,1);

%set bin width
%dst1 = 0.01;
%dst2 = 0.01; 
%dst3 = 0.015;
dst1 = dsts(1);
dst2 = dsts(2);
dst3 = dsts(3);
%break into section
St1 = St(1:p1_str);
St2 = St(p1_end:p2_str);
St3 = St(p2_end:end);
PSD1 = PSD(1:p1_str);
PSD2 = PSD(p1_end:p2_str);
PSD3 = PSD(p2_end:end);

for i = 1:3
    %% find the current case
    Stc = eval(sprintf("(St%d)",i));
    PSDc = eval(sprintf("(PSD%d)",i));
    dstc = eval(sprintf("(dst%d)",i));
    %% identify number of bins
    dSt = Stc(2)-Stc(1);
    Npts = floor(dstc/dSt);
    bin_str = 1:Npts:length(Stc)-Npts;
    Nbins = length(bin_str);
    PSDbc = zeros(1,Nbins);
    Stbc = zeros(1,Nbins);
for j = 1:Nbins %for each bin
    bin_end = bin_str(j) + Npts-1;
    St_bin = (Stc(bin_str(j))+Stc(bin_end))/2; %binned St
    PSD_bin = sum(PSDc(bin_str(j):bin_end))/Npts;%binned PSD
    PSDbc(j) = PSD_bin;
    Stbc(j) = St_bin;
end
    eval(sprintf("Stb%d = Stbc;",i))
    eval(sprintf("PSDb%d = PSDbc;",i))
end
%%
St_b = [Stb1,St(p1_str:p1_end)',Stb2,St(p2_str:p2_end)',Stb3];
PSD_b = [PSDb1,PSD(p1_str:p1_end)', PSDb2, PSD(p2_str:p2_end)',PSDb3];
%figure; semilogx(St_b, PSD_b)
end