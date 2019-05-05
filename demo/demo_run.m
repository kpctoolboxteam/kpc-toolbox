warning off

clear
fprintf('DEMO KPC-PH - medium variability / light-tail\n')
fprintf('Press any key to continue...\n')
pause
demo_kpcfit_ph_dec_pkt_1_udp

clear
fprintf('DEMO KPC - long-range dependence fit animation\n')
fprintf('Press any key to continue...\n')
pause
close all
demo_kpcfit_bcaug89_animate   

fprintf('DEMO KPC - long-range dependence fit\n')
fprintf('Press any key to continue...\n')
pause
close all
demo_kpcfit_bcaug89

fprintf('DEMO KPC - long-range dependence fit (3rd order)\n')
fprintf('Press any key to continue...\n')
pause
close all
demo_kpcfit_bcaug89_thirdord

clear
fprintf('DEMO KPC - long-range dependence fit, positive and negative acf - 8 states \n')
fprintf('Press any key to continue...\n')
pause
close all
demo_kpcfit_dec_pkt_1_udp8

clear
fprintf('DEMO KPC - long-range dependence fit, positive and negative acf - 16 states \n')
fprintf('Press any key to continue...\n')
pause
close all
demo_kpcfit_dec_pkt_1_udp16

