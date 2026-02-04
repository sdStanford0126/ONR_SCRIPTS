%% script for reading FWH-output
%frequency nondimensionalization (direct) would be
% h/c_inf ==> need to convert to De/(Mj*c_inf)
% TODO: rework and redefine everything
% to keep inline with experimental results
clear; clc; close all
%% setup the dimensionalization conversion ratio
Pref = 20E-6; %refernce pressure, Pa

rho = 1.2; %kg/m^3
c   = 346.1; %m/s %sufficiently close (346.0613
Punit = rho*c^2;

NPR = 3;
gamma = 1.4;
Mj = jetMach_isen_ideal(NPR,gamma)%jet mach number acoustic 
h_s = 1;
L_s = 2*h_s;
De_s = 1.3*(L_s*h_s)^(0.625)/(L_s+h_s)^(0.25)
fact_St = De_s/Mj      

De = De_s * 0.475 /39.37; %m
h = De/De_s;
fact_St_p = De / (Mj*c)

% pick which side to do the data on

%% mesh sizes
meshsizes = ["252M_IPR_1.45"]

for ind_case = 1:length(meshsizes)

    meshsize = meshsizes(ind_case)

switch meshsize
    case "131M"
        indir= "./" + meshsize + "/";
        datafolder = indir + "131M_s9_fwh_out/"
        surface = "s9"
        bin_flag = true
    case "168M"
        indir= "./" + meshsize + "/";
        datafolder = indir + "/s11_fwh_400ATU_out/"
        surface = "s11"
        bin_flag=true
    case "217M"
        indir = "./" + meshsize + "/";
        datafolder = indir + "s11_fwh_150ATU_out/"
        surface = "s11"
        bin_flag = false
    case "252M"
        indir = "./" + meshsize + "/";
        datafolder = indir + "s11_250_550500/"
        surface = "s11"
        bin_flag = false
    case "252M_IPR_1.45"
        indir = "./" + meshsize + "/";
        datafolder = indir + "s11_fwh_out_250_478500/"
        surface = "s11"
        bin_flag = false
    otherwise
        fprintf("input: %s, not defined, stopping...", meshsize)
        return
end
%
% 
% meshsize="131M"
% indir="./" + meshsize + "/"
% datafolder = indir + "s11_fwh_out/"


sides = ["Major","Minor"]
for ind_side = 1:length(sides)
    side = sides(ind_side)

outdir = "./temp_img/"
%outdir = "/Users/steven/OneDrive/Stanford/ONR project/results/port_inject/168M/s11_img_200/"%"./"+surface+"_img/"

%datafolder = indir+sprintf("/%s_%s_fwh_out/",meshsize,surface)
%load expr data

OASPLfname=sprintf("OASPL_vs_ang_expr_NPR3p0_%s.mat",side)
%PSDfname=sprintf("PSD_angs_expr_NPR3p0_%s.mat",side)
PSDSPLfname=sprintf("PSDSPL_angs_expr_NPR3p0_%s.mat",side)
exprfolder="/Users/steven/Library/CloudStorage/OneDrive-Personal/Stanford/ONR project/results/port_expr_acoustics";
load(exprfolder+"/"+OASPLfname)
%load(exprfolder+"/"+PSDfname)
load(exprfolder+"/"+PSDSPLfname)

%setup image save names

%% setup
thetas = 40:150;
OASPL = zeros(length(thetas),1);

%%

fname_fmt = datafolder + "%s_Mic%03d_OASPL.dat"
for i = 1: length(thetas)
    fname = sprintf(fname_fmt,side,i)
    fid = fopen(fname);
    line = fgetl(fid);
    line = fgetl(fid);
    line = fgetl(fid);
    line = fgetl(fid);
    
    A = str2num(line);
    SPL_c = A(6);
    OASPL(i)=SPL_c;
    
end
 
%% plot OASPL data
OASPL_d = 10*log10(OASPL*Punit*Punit/(Pref)^2);
figure; grid on
plot(thetas,OASPL_d,LineWidth=2,Color='r',DisplayName=sprintf("%s sim",meshsize))
hold on;
plot(ffAngle,OASPLintgCorr,LineWidth=2,LineStyle="-.",Color='k',DisplayName="Expr.")
title("OASPL vs. $\theta$", FontSize=14,Interpreter="latex")
ylabel("$OASPL (dB)$",FontSize=14,Interpreter="latex")
xlabel("$angle(^\circ)$",FontSize=14,Interpreter="latex")
if side == "Major"
    ylim([106,118])
else
    ylim([104,118])
end
legend(Location="southwest",FontSize=14)
saveas(gcf,outdir+sprintf("OASPL_theta_%s_%s_%s.png",meshsize,surface,side))
hold off;
%% find angles on Kaurab prelim dataset
ind1 = find(thetas >= 50,1); %50 degrees
ind2 = find(thetas >= 90,1); %90 degrees
ind3 = find(thetas >= 150,1); %150 degrees

%lets get some additional angles here, especially as we transition from
%upstream angle to sideline angle

ind1e = find(ffAngle >= 50,1); %50 degrees
ind2e = find(ffAngle >= 90,1); %90 degrees
ind3e = find(ffAngle >= 150,1); %150 degrees
%% read PSD data at these 3 select angles
fname_fmt = datafolder+ side +"_Mic%03d_PSD.dat";
[f1,PSD1] = readPSDfile(fname_fmt,ind1);
[f2,PSD2] = readPSDfile(fname_fmt,ind2);
[f3,PSD3] = readPSDfile(fname_fmt,ind3);

%fname_fmt = datafolder+ side +"_Mic%03d_p_hat.dat";
%f1,phat1] = readPhatfile(fname_fmt,ind1);
%[f2,phat2] = readPhatfile(fname_fmt,ind2);
%[f3,phat3] = readPhatfile(fname_fmt,ind3);
%% scale the data according to units
f_bw = f1(2)-f1(1)
St_bw = f_bw*fact_St;
f_sim = 1/0.05;
f_bwp = f_bw * c/(0.475/39.37)
f_bwp1 = f_bw * fact_St / fact_St_p
f_simp = f_sim * c/h
%an accurate correction is would be with f_bwp = 50;
%PSD_unit = Punit^2/204.8e3; %using sampling frequency for normalization?

PSD_unit = Punit^2/(c/h * fact_St_p); 
St1 = f1*fact_St;
St2 = f2*fact_St;
St3 = f3*fact_St;



PSD_d1 = PSD1  * PSD_unit;
PSD_d2 = PSD2  * PSD_unit;
PSD_d3 = PSD3  * PSD_unit;

%PSD_d1 = (phat1.*conj(phat1)).*PSD_unit;
%PSD_d2 = (phat2.*conj(phat2)).*PSD_unit;
%PSD_d3 = (phat3.*conj(phat3)).*PSD_unit;

St_b1 = St1;
St_b2 = St2;
St_b3 = St3;

PSD_b1 = PSD_d1;
PSD_b2 = PSD_d2;
PSD_b3 = PSD_d3;
 Nbin = 1000;
 

 %[St_b1, PSD_b1] = PSD_binning(St1, PSD_d1, Nbin);
 %[St_b2, PSD_b2] = PSD_binning(St2, PSD_d2, Nbin);
 %[St_b3, PSD_b3] = PSD_binning(St3, PSD_d3, Nbin);


 %high freq bin averaging
 % [St_b1, PSD_b1] = PSDbinning_HI(St1, PSD_d1, 0.05,1);
 % [St_b2, PSD_b2] = PSDbinning_HI(St2, PSD_d2, 0.05,1);
 % [St_b3, PSD_b3] = PSDbinning_HI(St3, PSD_d3, 0.05,1);

 % %screech avoiding bin averaging

 if bin_flag
    St_scr = 0.356898
    dsts   = [0.01,0.01,0.015]
    [St_b1, PSD_b1] = PSD_binning_scr(St1,PSD_d1,St_scr,dsts);
    [St_b2, PSD_b2] = PSD_binning_scr(St2,PSD_d2,St_scr,dsts);
    [St_b3, PSD_b3] = PSD_binning_scr(St3,PSD_d3,St_scr,dsts);
 end
SPL1 = 10*log10(PSD_b1*St_bw/Pref^2);
SPL2 = 10*log10(PSD_b2*St_bw/Pref^2);
SPL3 = 10*log10(PSD_b3*St_bw/Pref^2);

%SPL1 = 10*log10(PSD_b1*50/Pref^2);
%SPL2 = 10*log10(PSD_b2*50/Pref^2);
%SPL3 = 10*log10(PSD_b3*50/Pref^2);

PSDex = 10.^(micSpectraCorr./10)./(50*fact_St_p).*Pref^2;

%%
thetas = 40:150;
figure;
avg_flag = false
for i = 1:length(ffAngle)
    th = ffAngle(i);
    fname_fmt = datafolder+ side +"_Mic%03d_PSD.dat";
    ind_PSD = find(abs(thetas-th)<0.5);
    assert(th == thetas(ind_PSD))
    [fc,PSDc] = readPSDfile(fname_fmt,ind_PSD,avg_flag);
    Stc = fc*fact_St;
    fmin = fc(2);
    Stcmin = Stc(2);
    PSD = PSDc;
    PSD_unit = Punit^2/(c/h*fact_St_p);
    PSDc = PSDc*PSD_unit;
    Stex = freqRange*fact_St_p;
    title_name = "PSD for " + side +  " $\theta =" + num2str(th) + "^\circ$";
    img_name = outdir+"PSD_" + num2str(th)  +"deg_"+surface+ "_" + side +"_"+meshsize+".png"
    figure; 
    loglog(Stc,PSDc,'r-',"DisplayName","sim",...
        LineWidth=1), hold on;
    loglog(Stex,PSDex(i,:),'-b',"DisplayName", "expr")
    xlim([1e-2,5])
    xlabel("$St$")
    ylabel("$PSD [Pa^2/St]$")
    title(title_name)
    legend()
    saveas(gcf,img_name);
end 


%% save data for surface comparison plotting
savename = sprintf("port_fwh_post_data_%s_%s_%s.mat",meshsize,surface,side)
indse = [ind1e,ind2e,ind3e];
save(savename,"St_b1","St_b2","St_b3","PSD_b1","PSD_b2","PSD_b3","thetas","OASPL_d","Stex","PSDex","OASPLintgCorr","ffAngle","side","meshsize","surface","side")
%%
figure;
semilogx(St_b1,SPL1,'-r',"DisplayName","$\theta = 50^\circ$",LineWidth=1)
hold on;
semilogx(St_b2,SPL2,'-b',"DisplayName","$\theta = 90^\circ$")
semilogx(St_b3,SPL3,'-k',"DisplayName","$\theta = 150^\circ$")
legend('Location','southwest',Interpreter='latex')
xlim([0.02,3])
%ylim([60,170])
xlabel("$St$",FontSize=14,Interpreter="latex")
ylabel("$SPL(dB)$",FontSize=14,Interpreter="latex")
title(sprintf("SPL at different angles, NPR = 3.0, %s mesh, baseline",meshsize),FontSize=14)
saveas(gcf,outdir+sprintf("SPL_diff_theta_%s_%s_%s.png",meshsize,surface,side))



%% direct SPL comparison with experimental results
%% 50 deg
figure;
loglog(St_b1,SPL1,'-r',"DisplayName","$\theta = 50^\circ$, Sim.",LineWidth=2)
hold on;
grid on;
gca.GridLineWidth = 2;
loglog(freqRange*fact_St_p,micSpectraCorr(ind1e,:),'-b',"DisplayName",...
    "$\theta = 50^\circ$, Expr no atten.",LineWidth=1)
legend('Location','southwest',Interpreter='latex',FontSize=14)
xlabel("$St$",FontSize=14,Interpreter="latex")
ylabel("$SPL(dB)$",FontSize=14,Interpreter="latex")
xlim([0.01,4.61])
title(sprintf("PSD at 50 deg, NPR = 3.0, %s mesh vs Expr, Baseline",meshsize),FontSize=14)
saveas(gcf,outdir+sprintf("PSD_%s_expr_50deg_%s_%s.png",meshsize,surface,side))
%% 90 deg
figure;
loglog(St_b2,SPL2,'-r',"DisplayName","$\theta = 90^\circ$, Sim.",LineWidth=2)
hold on;
grid on;
gca.GridLineWidth = 2;
loglog(freqRange*fact_St_p,micSpectraCorr(ind2e,:),'-b',"DisplayName",...
    "$\theta = 90^\circ$, Expr no atten.",LineWidth=1)

legend('Location','southwest',Interpreter='latex',FontSize=14)
xlim([0.01,4.61])
xlabel("$St$",FontSize=14,Interpreter="latex")
ylabel("$SPL(dB)$",FontSize=14,Interpreter="latex")
title(sprintf("PSD at 90 deg, NPR = 3.0, %s mesh vs Expr, Baseline",meshsize),FontSize=14)
saveas(gcf,outdir+sprintf("PSD_%s_expr_90deg_%s_%s.png",meshsize,surface,side))

%% 150 deg
figure;
loglog(St_b3,SPL3,'-r',"DisplayName","$\theta = 150^\circ$, Sim.",LineWidth=2)
hold on;
grid on;
gca.GridLineWidth = 2;
loglog(freqRange*fact_St_p,micSpectraCorr(ind3e,:),'-b',"DisplayName",...
    "$\theta = 150^\circ$, Expr no atten.",LineWidth=1)

legend('Location','southwest',Interpreter='latex',FontSize=14)
xlim([0.01,4.61])
xlabel("$St$",FontSize=14,Interpreter="latex")
ylabel("$SPL(dB)$",FontSize=14,Interpreter="latex")
title(sprintf("PSD at 150 deg, NPR = 3.0, %s mesh vs Expr, Baseline",meshsize),FontSize=14)
saveas(gcf,outdir+sprintf("PSD_%s_expr_150deg_%s_%s.png",meshsize,surface,side))
end
end