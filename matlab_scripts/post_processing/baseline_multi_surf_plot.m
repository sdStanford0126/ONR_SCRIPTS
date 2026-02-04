%% Jun 13 quick plot multiple surface

close all;clc;clear; setPlotPref(4,'latex',18);
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
indir = "./AIAA_Aviation_FW-H_post/FWH_data/"
data_fmt = indir + "/baseline_fwh_post_data_%s_%s_%s.mat"

%full set comparison
%meshsizes = ["47M", "68M","94M","130M","144M"]
%surfaces = ["s2","s9","s4","s9","s9"]
%linstys = ["g-","c-","b-","r-","m-"]


%side = "Major" 
%side = "Minor"

%subset comparison
%meshsizes = fliplr(["94M","130M","144M"])
%surfaces = fliplr(["s4","s9","s9"])
%linstys = fliplr(["b-.","g--","r-"])
%sides = ["Major","Minor"]
meshsizes = fliplr(["94M","130M","144M","151M","119M"])
surfaces = fliplr(["s4","s9","s9","s9","s9"])
linstys = fliplr(["b-.","g--","r-","c-x","m-o"])
sides = ["Major","Minor"]
for ind_side = 1:length(sides)
    side = sides(ind_side)
    assert(length(surfaces) == length(linstys) && length(linstys) == length(meshsizes))
for i = 1:length(meshsizes)
    meshsize = meshsizes(i);
    surface = surfaces(i);
    data_fname = sprintf(data_fmt,meshsize,surface, side)
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
St_ex_baseline = Stex;
PSDex_baseline = PSDex;
ffAngle_baseline = ffAngle;
OASPL_ex_baseline = OASPLintgCorr; 

ind1 = find(ffAngle == 50);
ind2 = find(ffAngle == 90);
ind3 = find(ffAngle == 150);

for i = 1:3
    i_str = num2str(i);
    eval(sprintf("PSD%s_ex = PSDex_baseline(ind%s,:);",i_str,i_str));
    eval(sprintf("SPL%s_ex = 10*log10(PSD%s_ex./Pref^2);",i_str,i_str));
end


%% add baseline expr data here (not needed here)
% OASPLfname=sprintf("OASPL_vs_ang_expr_NPR3p0_%s.mat",side)
% %PSDfname=sprintf("PSD_angs_expr_NPR3p0_%s.mat",side)
% PSDSPLfname=sprintf("PSDSPL_angs_expr_NPR3p0_%s.mat",side)
% exprfolder="/Users/steven/Library/CloudStorage/OneDrive-Personal/Stanford/ONR project/results/baseline_expr_acoustics";
% load(exprfolder+"/"+OASPLfname)
% %load(exprfolder+"/"+PSDfname)
% load(exprfolder+"/"+PSDSPLfname)
% 
% OASPL_ex_base = OASPLintgCorr;
% ffAngle_base = ffAngle;
% SPL_PSD_ex_base = micSpectraCorr;
% 
% PSD_ex_base = 10.^(SPL_PSD_ex_base./10)./(50*fact_St_p).*Pref^2;

%% plotting PSD figures
titles = ["$PSD,\quad 50^\circ$",...
    "$PSD,\quad 90^\circ$","$PSD,\quad 150^\circ$"]
outdir = "./";
imgnames = [side + "_baseline_PSD_mesh_50.png", ...
    side + "_baseline_PSD_mesh_90.png",...
    side + "_baseline_PSD_mesh_150.png"];

for i = 1:3
    figure;
    for j = 1:length(meshsizes)
        
        meshsize = meshsizes(j)
        dSt_c = eval(sprintf("St%s_%s(2) - St%s_%s(1)"...
            ,num2str(i),meshsize,num2str(i),meshsize));
        df_s = dSt_c / fact_St;
        ATU = 1/df_s
        
        linsty = linstys(j);
        line_name = meshsize %+ ", T=" + num2str(ATU)
        eval(sprintf("loglog(St%s_%s,PSD%s_%s,""%s"",DisplayName=""%s"")", ...
            num2str(i),meshsize,num2str(i),meshsize,linsty,line_name))
        hold on
    end
    %plot experimental result here
    eval(sprintf("loglog(Stex,PSD%s_ex,""k-"",DisplayName=""Expr."")",num2str(i)))
    %figure labeling
    title(titles(i)+" "+side,FontSize=16)
    legend(Location="southwest")
    xlabel("$St$")
    ylabel("$PSD,\quad [Pa^2/St]$")
    xlim([0.015,2])
    ylim([1e-4,5e4])
    saveas(gcf,outdir+imgnames(i))
end

titles = ["$PSD,\quad 50^\circ$",...
    "$PSD,\quad 90^\circ$","$PSD,\quad 150^\circ$"]
outdir = "./";
imgnames = [side + "_baseline_SPL_mesh_50.png", ...
    side + "_baseline_SPL_mesh_90.png",...
    side + "_baseline_SPL_mesh_150.png"];

for i = 1:3
    figure;
    for j = 1:length(meshsizes)
        
        meshsize = meshsizes(j)
        dSt_c = eval(sprintf("St%s_%s(2) - St%s_%s(1)"...
            ,num2str(i),meshsize,num2str(i),meshsize));
        df_s = dSt_c / fact_St;
        ATU = 1/df_s
        
        linsty = linstys(j);
        line_name = meshsize %+ ", T=" + num2str(ATU)
        eval(sprintf("semilogx(St%s_%s,SPL%s_%s,""%s"",DisplayName=""%s"")", ...
            num2str(i),meshsize,num2str(i),meshsize,linsty,line_name))
        hold on
    end
    %plot experimental result here
    eval(sprintf("semilogx(Stex,SPL%s_ex,""k-"",DisplayName=""Expr."")",num2str(i)))
    %figure labeling
    title(titles(i)+" "+side,FontSize=16)
    legend(Location="southwest")
    xlabel("$St$")
    ylabel("$PSD,\quad [dB/St]$")
    xlim([0.015,2])
    ylim([60,150])
    saveas(gcf,outdir+imgnames(i))
end



%% for OASPL

figure
for j = 1:length(meshsizes)
            meshsize = meshsizes(j)
        dSt_c = eval(sprintf("St%s_%s(2) - St%s_%s(1)"...
            ,num2str(i),meshsize,num2str(i),meshsize));
        df_s = dSt_c / fact_St;
        ATU = 1/df_s
        linsty = linstys(j);
        line_name = meshsize %+ ", T=" + num2str(ATU)
        eval(sprintf("plot(thetas_%s,OASPL_%s,""%s"",DisplayName=""%s"")", ...
            meshsize,meshsize,linsty,line_name))
        hold on
end
    title("OASPL "+side,FontSiz=16)
    plot(ffAngle_baseline,OASPL_ex_baseline,"k-",DisplayName="Expr.")
    legend(Location="northwest")
    xlabel("$\theta$")
    ylabel("$OASPL,\, [dB]$")
    saveas(gcf,outdir+side+"_OASPL.png")
end