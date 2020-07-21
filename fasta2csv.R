input <- readLines("/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_PPOL/bio_info_PPOL/DiscoSnpRad/2020/Galaxy77-[Fasta_with_DiscoSnpRAD_on_R1R2_clonefiltered].fasta")
output <- file("/Volumes/data/sync/poissons/genetique_poissons/RAD-seq_PPOL/bio_info_PPOL/DiscoSnpRad/2020/Galaxy77-[Fasta_with_DiscoSnpRAD_on_R1R2_clonefiltered].csv","w")

currentSeq <- 0
newLine <- 0

for(i in 1:length(input)) {
  if(strtrim(input[i], 1) == ">") {
    if(currentSeq == 0) {
      writeLines(paste(input[i],"\t"), output, sep="")
      currentSeq <- currentSeq + 1
    } else {
      writeLines(paste("\n",input[i],"\t", sep =""), output, sep="")
    }
  } else {
    writeLines(paste(input[i]), output, sep="")
  }
}

close(output)

