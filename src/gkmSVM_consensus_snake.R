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

input_bed_train=args[1]
input_consensus=args[2]
output_fasta_train=args[3]
output_neg_bed_train=args[4]
output_neg_fasta_train=args[5]
params_svmfnprfx=args[6]
output_kernel=args[7]
output_CVpredfn=args[8]
output_ROCfn=args[9]
output_PDFfn=args[10]
output_weights=args[11]

#######
#######
## MAIN
#######
#######

## generating null sequences from functional regions excluding TEs 
## ~ 2 Gb with HMM subset
## ~ 2 Gb with top 1000 subset
start = Sys.time()
genNullSeqs(inputBedFN = input_bed_train,
            nMaxTrials=10,
            xfold=1,
            genome=BSgenome.Hsapiens.UCSC.hg38.masked, ## masked from TEs
            outputPosFastaFN=output_fasta_train,
            outputBedFN=output_neg_bed_train,
            outputNegFastaFN=output_neg_fasta_train)

end = Sys.time()
print("genNullSeqs:")
print(start)
print(end)
print(end-start)
gc()

## generating kernel matrix
## maxed at 42.2 Gb, 296k fasta rows
## 2 minutes, 1 Gb with top 1000 (12k fasta rows)
start = Sys.time()
gkmsvm_kernel(posfile = output_fasta_train,
              negfile = output_neg_fasta_train,
              outfile = output_kernel)
end = Sys.time()
print("gkmsvm_kernel:")
print(start)
print(end)
print(end-start)
gc()

## training model
## ~ 1Gb remains occupied with top 1000
start = Sys.time()
gkmsvm_trainCV(kernelfn = output_kernel,
               posfn = output_fasta_train,
               negfn = output_neg_fasta_train,
               svmfnprfx = params_svmfnprfx,
               nCV=5,
               outputCVpredfn = output_CVpredfn,
               outputROCfn = output_ROCfn,
               outputPDFfn = output_PDFfn)
end = Sys.time()
print("gkmsvm_trainCV:")
print(start)
print(end)
print(end-start)
gc()

## testing model
## test on Dfam consensus
start = Sys.time()
gkmsvm_classify(seqfile = input_consensus,
                svmfnprfx = params_svmfnprfx,
                outfile = output_weights)
end = Sys.time()
print("gkmsvm_classify:")
print(start)
print(end)
print(end-start)
gc()

##
