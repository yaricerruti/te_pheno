#!/usr/bin/Rscript

args = commandArgs(T)

packages = c("gkmSVM",
             "BSgenome.Hsapiens.UCSC.hg38.masked")

for(package in packages) {
  library(package, character.only = T)
}

############
## VARIABLES
############

input_consensus=args[1]
input_instances=args[2]
params_svmfnprfx=args[3]
output_weights_instances=args[4]


#######
#######
## MAIN
#######
#######


## testing on Dfam consensus
#start = Sys.time()
#gkmsvm_classify(seqfile = input_consensus,
#                svmfnprfx = params_svmfnprfx,
#                outfile = output_weights_consensus)
#end = Sys.time()
#print("gkmsvm_classify consensus:")
#print(start)
#print(end)
#print(end-start)
#gc()

## test on TE fasta sequences
start = Sys.time()
gkmsvm_classify(seqfile = input_instances,
                svmfnprfx = params_svmfnprfx,
                outfile = output_weights_instances)
end = Sys.time()
print("gkmsvm_classify instances:")
print(start)
print(end)
print(end-start)
gc()

##

