/***********************************************************************/
/*Replication Dataset for Gerber et al 2008 APSR Social Pressure Study */
/*Alan Gerber (alan.gerber@yale.edu)                                   */
/***********************************************************************/

****loading data****
clear
set mem 450m
use "Replication Dataset for Gerber et al 2008 APSR Social Pressure Study.dta", clear

****Table 1****
collapse (max) treatment hh_size (mean) g2002 g2000 p2004 p2002 p2000 sex yob, by(hh_id)
bysort treatment: sum hh_size g2002 g2000 p2004 p2002 p2000 sex yob

****MNL reported on p.37 of APSR article****
mlogit treatment hh_size g2002 g2000 p2004 p2002 p2000 sex yob, baseoutcome(0)
save "Replication Dataset for Gerber et al 2008 APSR Social Pressure Study-household level.dta"

****Table 2****
use "Replication Dataset for Gerber et al 2008 APSR Social Pressure Study.dta", clear
tabulate voted treatment

****Table 3****
gen control = treatment==0
gen hawthorne = treatment==1
gen civicduty = treatment==2
gen neighbors = treatment==3
gen self = treatment==4

****Table 3, model a****
regress voted civicduty hawthorne self neighbors, robust cluster(hh_id)

****Table 3, model b****
areg voted civicduty hawthorne self neighbors, absorb(cluster) robust cluster(hh_id)

****Table 3, model c****
areg voted civicduty hawthorne self neighbors g2002 g2000 p2004 p2004 p2002 p2000, absorb(cluster) robust cluster(hh_id)







