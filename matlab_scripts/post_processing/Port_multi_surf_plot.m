%% Jun 20 quick plot multiple surface

close all;clc;clear; 
setPlotPref(4.0,'latex',18);
%% 
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

%% 
f_sim = 1/0.05
St_sim = f_sim * fact_St
%%

indir = "./AIAA_Aviation_FW-H_post/FWH_data/"
outdir = "./AIAA_Aviation_FW-H_post/"
data_fmt = "port_fwh_post_data_%s_%s_%s.mat"

meshsizes = ["168M","217M","252M"]


sides = ["Major","Minor"]

for ind_side =  1:length(sides)
    side = sides(ind_side)
for i = 1:length(meshsizes)
    meshsize = meshsizes(i)
    surface = "s11"

    data_fname = indir + sprintf(data_fmt,meshsize,surface, side)
    load(data_fname)

    eval(sprintf("St1_%s = St_b1;",meshsize))
    eval(sprintf("St2_%s = St_b2;",meshsize))
    eval(sprintf("St3_%s = St_b3;",meshsize))


    eval(sprintf("PSD1_%s = PSD_b1;",meshsize))
    eval(sprintf("PSD2_%s = PSD_b2;",meshsize))
    eval(sprintf("PSD3_%s = PSD_b3;",meshsize))

    eval(sprintf("SPL1_%s = 10*log10(PSD_b1./Pref^2);",meshsize))
    eval(sprintf("SPL2_%s = 10*log10(PSD_b2./Pref^2);",meshsize))
    eval(sprintf("SPL3_%s = 10*log10(PSD_b3./Pref^2);",meshsize))
    
    eval(sprintf("thetas_%s = thetas;",meshsize))
    eval(sprintf("OASPL_%s = OASPL_d;",meshsize))
end
%St_ex_port = Stex;
%PSD_ex_port = PSDex;
%ffAngle_port = ffAngle;
%OASPL_ex_port = OASPLintgCorr; 
linstys = ["r-","m-","g-"]

 %% add baseline expr data here
OASPLfname=sprintf("OASPL_vs_ang_expr_NPR3p0_%s.mat",side)
%PSDfname=sprintf("PSD_angs_expr_NPR3p0_%s.mat",side)
PSDSPLfname=sprintf("PSDSPL_angs_expr_NPR3p0_%s.mat",side)
exprfolder="/Users/steven/Library/CloudStorage/OneDrive-Personal/Stanford/ONR project/results/baseline_expr_acoustics";
load(exprfolder+"/"+OASPLfname)
%load(exprfolder+"/"+PSDfname)
load(exprfolder+"/"+PSDSPLfname)

OASPL_ex_base = OASPLintgCorr;
ffAngle_base = ffAngle;
SPL_PSD_ex_base = micSpectraCorr;

PSD_ex_base = 10.^(SPL_PSD_ex_base./10)./(50*fact_St_p).*Pref^2;

St_ex_base = freqRange * fact_St_p;
ind1 = find(ffAngle == 50);
ind2 = find(ffAngle == 90);
ind3 = find(ffAngle == 150);

for i = 1:3
    i_str = num2str(i);
    eval(sprintf("PSD%s_base_ex = PSD_ex_base(ind%s,:);",i_str,i_str));
    eval(sprintf("SPL%s_base_ex = 10*log10(PSD%s_base_ex./Pref^2);",i_str,i_str));
end

%% add port expr data
OASPLfname=sprintf("OASPL_vs_ang_expr_NPR3p0_%s.mat",side)
%PSDfname=sprintf("PSD_angs_expr_NPR3p0_%s.mat",side)
PSDSPLfname=sprintf("PSDSPL_angs_expr_NPR3p0_%s.mat",side)
exprfolder="/Users/steven/Library/CloudStorage/OneDrive-Personal/Stanford/ONR project/results/port_expr_acoustics";
load(exprfolder+"/"+OASPLfname)
%load(exprfolder+"/"+PSDfname)
load(exprfolder+"/"+PSDSPLfname)

OASPL_ex_port = OASPLintgCorr;
ffAngle_port = ffAngle;
SPL_PSD_ex_port = micSpectraCorr;

PSD_ex_port = 10.^(SPL_PSD_ex_port./10)./(50*fact_St_p).*Pref^2;
St_ex_port = freqRange * fact_St_p;
ind1 = find(ffAngle == 50);
ind2 = find(ffAngle == 90);
ind3 = find(ffAngle == 150);

for i = 1:3
    i_str = num2str(i);
    eval(sprintf("PSD%s_port_ex = PSD_ex_port(ind%s,:);",i_str,i_str));
    eval(sprintf("SPL%s_port_ex = 10*log10(PSD%s_port_ex./Pref^2);",i_str,i_str));
end


%% add simulation baseline data here

BL_data_fmt = "baseline_fwh_post_data_144M_s9_%s.mat";
BL_data_Name = sprintf(indir + BL_data_fmt,side)
load(BL_data_Name)

St1_bl = St_b1;
St2_bl = St_b2;
St3_bl = St_b3;

PSD1_bl = PSD_b1;
PSD2_bl = PSD_b2;
PSD3_bl = PSD_b3;

OASPL_bl = OASPL_d;


SPL1_bl = 10*log10(PSD1_bl./Pref^2);

SPL2_bl = 10*log10(PSD2_bl./Pref^2);

SPL3_bl = 10*log10(PSD3_bl./Pref^2);

BL_data_fmt = "baseline_fwh_post_data_130M2_s9_%s.mat";
BL_data_Name = sprintf(indir + BL_data_fmt,side)
load(BL_data_Name)

St1_bl2 = St_b1;
St2_bl2 = St_b2;
St3_bl2 = St_b3;

PSD1_bl2 = PSD_b1;
PSD2_bl2 = PSD_b2;
PSD3_bl2 = PSD_b3;

OASPL_bl2 = OASPL_d;



%% OASPL
figure
for j = 1:length(meshsizes)
            meshsize = meshsizes(j)
        dSt_c = eval(sprintf("St%s_%s(2) - St%s_%s(1)"...
            ,num2str(i),meshsize,num2str(i),meshsize));
        df_s = dSt_c / fact_St;
        ATU = 1/df_s
        linsty = linstys(j);
        line_name = "Port, "+meshsize %+ ", T=" + num2str(ATU)
        eval(sprintf("plot(thetas_%s,OASPL_%s,""%s"",DisplayName=""%s"")", ...
            meshsize,meshsize,linsty,line_name))
        hold on
end
%plot(thetas,OASPL_bl2,'c-.',DisplayName="Baseline, 130M")
plot(thetas,OASPL_bl,'b:',DisplayName="Baseline, 144M")
plot(ffAngle_base,OASPL_ex_base,'k:',DisplayName="Baseline, Expr")
plot(ffAngle_port,OASPL_ex_port,'k-',DisplayName="Port, Expr")
legend(location="northwest")
xlabel("$\theta$")
ylabel("$OASPL, \, [dB]$")
xlim([40,150])
%box off
saveas(gcf,outdir + "/" + side + "_port_vs_base.png")
%% plotting PSD figures
titles = ["$PSD,\quad 50^\circ$",...
    "$PSD,\quad 90^\circ$","$PSD,\quad 150^\circ$"]
imgnames = [side + "_port_PSD_mesh_50.png", ...
    side + "_port_PSD_mesh_90.png",...
    side + "_port_PSD_mesh_150.png"]

for i = 1:3
    figure;
    for j = 1:length(meshsizes)
        meshsize = meshsizes(j)
        linsty = linstys(j);
        line_name = "Port, "+meshsize %+ ", T=" + num2str(ATU)
        eval(sprintf("loglog(St%d_%s,PSD%d_%s,""%s"",DisplayName=""%s"") ",...
            i,meshsize,i,meshsize,linsty,line_name))
        hold on
    end
    eval(sprintf("loglog(St%d_bl,PSD%d_bl, 'b:', DisplayName=""Baseline,144M"")",i,i));
    legend()
    legend(Location="southwest")
    title(titles(i)+" "+side,FontSize=16)
    xlabel("$St$")
    ylabel("$PSD,\quad [Pa^2/St]$")
    xlim([0.015,6])
    ylim([1e-4,5e4])
    %box off
    saveas(gcf,outdir+imgnames(i))
end




titles = ["$PSD,\quad 50^\circ$",...
    "$PSD,\quad 90^\circ$","$PSD,\quad 150^\circ$"]
outdir = "./AIAA_Aviation_FW-H_post/";
imgnames = [side + "_port_SPL_mesh_50.png", ...
    side + "_port_SPL_mesh_90.png",...
    side + "_port_SPL_mesh_150.png"];

for i = 1:3
    figure;
    for j = 1:length(meshsizes)
        
        meshsize = meshsizes(j)
        dSt_c = eval(sprintf("St%s_%s(2) - St%s_%s(1)"...
            ,num2str(i),meshsize,num2str(i),meshsize));
        df_s = dSt_c / fact_St;
        ATU = 1/df_s
        
        linsty = linstys(j);
        line_name = "Port, "+meshsize %+ ", T=" + num2str(ATU)
        eval(sprintf("semilogx(St%s_%s,SPL%s_%s,""%s"",DisplayName=""%s"")", ...
            num2str(i),meshsize,num2str(i),meshsize,linsty,line_name))
        hold on
    end
    %plot experimental result here
    eval(sprintf("semilogx(St%d_bl,SPL%d_bl, 'b:', DisplayName=""Baseline,144M"")",i,i));
    %eval(sprintf("semilogx(St_ex_base,SPL%d_base_ex,'k-',DisplayName=""Baseline,Expr"")",i));
    eval(sprintf("semilogx(St_ex_port,SPL%d_port_ex,'k-',DisplayName=""Port,Expr"")",i));
    %figure labeling
    title(titles(i)+" "+side,FontSize=16)
    legend(Location="southwest")
    xlabel("$St$")
    ylabel("$PSD,\quad [dB/St]$")
    xlim([0.015,4])
    ylim([60,150])
    %box off
    saveas(gcf,outdir+imgnames(i))
end

end