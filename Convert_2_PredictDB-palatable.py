from dosageparser import DosageParser

gtpath="/usr/users/akaphle/GTEx_VCF.files/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr####.gz"
samplepath="/usr/users/akaphle/GTEx_VCF.files/donor_ids.fam"
outfile="/usr/users/akaphle/PredictDB_pipeline/Klinikum_Predixscan/GTEx_Dosage/Genotype/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr####"

#Read Genotype
ds = DosageParser(gtpath, samplepath, 1, 1)
dosage = ds.dosage
snpinfo = ds.snpinfo
donorids = ds.sample_id
nsnps = ds.nsnps
nsample = ds.nsample

def write_dosages(outfile, snps, genotype,samplenames):
    # I'm ommiting the above line because genotype sorting is done outside the package...
    filtered_genotype = list()

    with open(outfile+".annot", 'w') as outstream2:
        with open(outfile+".gz", 'w') as outstream:
            headers = "Id "+" ".join(samplenames)+"\n"
            annotheader = "\t".join(["chr","pos","varID","refAllele","effectAllele","rsid"])
            outstream.write(headers)
            outstream2.write(annotheader+"\n")
            for i, snp in enumerate(snps):
                variant_id = "_".join([str(snp.chrm), str(snp.bp_location), snp.ref_allele, snp.alt_allele, "b37"])
                #annot_header = " ".join(["chr", "position", "VariantID", "RefAllele", "AlternativeAllele", "rsid", "rsid"])
                dosage_row = " ".join([variant_id] + list(map("{:.3f}".format, genotype[i])))
                outstream.write(dosage_row+"\n")

                annotline = "\t".join([str(snp.chrm), str(snp.bp_location), variant_id, snp.ref_allele, snp.alt_allele, snp.rsid])
                outstream2.write(annotline+"\n")

if __name__ == '__main__':
    # execute only if run as the entry point into the program
    write_dosages(outfile,snpinfo,dosage,donorids)
