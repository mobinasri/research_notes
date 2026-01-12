import os
import argparse
from cyvcf2 import VCF, Writer


def parse_args():
    p = argparse.ArgumentParser(
        description="For each sample, write variants where the sample is non-ref (not 0/0 and not missing) and "
                    "the minimum ALT allele frequency among the sample's ALT alleles is < AF threshold."
    )
    p.add_argument("--input", required=True, help="Input VCF path")
    p.add_argument("--output-dir", required=True, help="Output directory")
    p.add_argument("--samples", required=True,
                   help="Comma-separated list of sample names (e.g. HG01255,HG02280)")
    p.add_argument("--af-threshold", type=float, default=0.01,
                   help="AF threshold (default: 0.01). Keep variants with min ALT AF < threshold.")
    p.add_argument("--threads", type=int, default=4, help="cyvcf2 threads (default: 4)")
    return p.parse_args()


def main():
    args = parse_args()

    af_threshold = args.af_threshold
    samples = [s for s in args.samples.split(",") if s]
    input_vcf_path = args.input
    output_dir = args.output_dir
    threads = args.threads

    os.makedirs(output_dir, exist_ok=True)
    
    for sample in samples:
        input_reader = VCF(input_vcf_path,
                  samples = [sample],
                  threads=threads, 
                  gts012=True,
                  strict_gt=True)
    
        
        output_vcf_path = f"{output_dir}/{os.path.basename(input_vcf_path)}".replace(".vcf",f".{sample}.vcf")
        output_writer = Writer(output_vcf_path, input_reader)

        print(f"Writing {output_vcf_path}")

        for variant in input_reader:
            
            gt_type = variant.gt_types[0]
            # genotype 0 means hom-ref and 3 means at least one unknown allele
            if (gt_type == 3 or gt_type == 0):
                continue
            
            if isinstance(variant.INFO["AF"], tuple):
                af_list = list(variant.INFO["AF"])
            else:
                af_list = [variant.INFO["AF"]]
            
            # get allele indices
            alleles = set(variant.genotypes[0][:2])
            # remove ref allele is present
            alleles.discard(0)
    
            # allele index 1 is the first alternative allele and af_tuple does have frequency for ref allele 
            min_allele_freq = min([af_list[allele - 1] for allele in alleles])
            if min_allele_freq < af_threshold:
                output_writer.write_record(variant)
        
        input_reader.close()
        output_writer.close()

if __name__ == "__main__":
    main()
