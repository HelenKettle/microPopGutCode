#compare 1.5.1 with 1.5
rm(list=ls())
graphics.off()

library(devtools)

#load new package----------------------------------------------
newpf='/home/helen/MicroPop/MicroPopPackageDec21/microPop/'
load_all(newpf)
print(packageVersion("microPop"))

build(newpf)

#instead of install.packages I found command line (linux prompt in dir with tar file) below worked better
#   R CMD INSTALL microPop_1.5.1.tar.gz



#install.packages('/home/helen/MicroPop/MicroPopPackageDec21/microPop',repos=NULL,type='source')
#install.packages('/home/helen/MicroPop/MicroPopPackageDec21/microPop_1.5.1.tar.gz',repos=NULL,type='source')

library(microPop)

browseVignettes(package="microPop") #choose html




