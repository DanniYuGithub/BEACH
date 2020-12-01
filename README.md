Biometric Exploratory Analysis Creation House (BEACH) is a shiny app that provides automation platform for users.

Before running BEACH, please make sure your computer is connected to internet and the following packages are installed.

dep.packages <- c("shiny", "DT", "haven", "xtable", "rtf", "plyr", "sas7bdat", "WriteXLS", "rJava");
na.packages <- dep.packages[!dep.packages %in% installed.packages()]
if (length(na.packages)>0) install.packages(na.packages);

if(!"sas7bdat.parso" %in% installed.packages()) devtools::install_github('BioStatMatt/sas7bdat.parso', force=TRUE)



Please set up your default internet browser as google chrome and then, in your R console, please run the following code to run BEACH locally.

library(shiny);
runGitHub("BEACH", "DanniYuGithub");


To install the package from R cran, please check the link https://cran.r-project.org/web/packages/BEACH/index.html
library(shiny); library(DT); library(BEACH); runBEACH()


Implementation: BEACH is written mainly in R/Shiny combined with some HTML & CSS code. Using this package, users can create customized analyses and make them available to end users who can perform interactive analyses and save analyses in the format of RTF, PDF or HTML. It allows developers to focus on R code for analysis, instead of dealing with html or shiny code. The available analysis modules in this package are DataQC, TrialDesign, SEM-BMK, visual_SDTM-LB, and LightON. 


A presentation slide about the foundmental idea of BEACH can be downloaded at https://www.pharmasug.org/proceedings/2018/AD/PharmaSUG-2018-AD05.pdf. 
