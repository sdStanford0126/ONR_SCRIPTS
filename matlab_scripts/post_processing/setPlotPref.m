function setPlotPref(lw,inter,fontsize)
%setup default interpreter and linewidth
set(groot,'defaultLineLineWidth',lw) 
set(groot,'DefaultaxesLineWidth', 1.5) 
set(groot,'DefaultaxesFontSize', fontsize) 
set(groot,'DefaultTextFontSize', fontsize) 
set(groot, 'defaultAxesTickLabelInterpreter',inter);
set(groot,'defaultTextInterpreter',inter)
set(groot, 'defaultLegendInterpreter',inter);
end