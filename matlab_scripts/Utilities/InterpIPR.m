%% May 14th quick code on interpolated IPR data
%interpolating from two exisitng IPRs to find new IPR match experimental
%IPRds
IPR1 = 2.25
IPR2 = 1.45

MFR1 = 2.11
MFR2 = 1.03

MFR_tgt = 1.75

IPR_diff = IPR1-IPR2
MFR_diff = MFR1-MFR2

slope = IPR_diff/MFR_diff

IPR_tgt = IPR2 + (MFR_tgt-MFR2)*slope

gamma = 1.4;
P_tot = IPR_tgt * 1/gamma