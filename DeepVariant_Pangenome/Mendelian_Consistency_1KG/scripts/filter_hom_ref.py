#!/usr/bin/env python3

from cyvcf2 import VCF, Writer
import argparse
import time

def main():
    parser = argparse.ArgumentParser(
        description="Filter VCF to keep variants where at least one sample is non-hom-ref"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input VCF file (can be compressed)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output VCF file (can be compressed)"
    )
    parser.add_argument(
        "--samples",
        required=True,
        help="Comma-separated list of sample names"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads (default: 4)"
    )

    args = parser.parse_args()

    input_file = args.input
    output_file = args.output
    samples = args.samples.split(",")   # list of strings
    threads = args.threads

    start_time_secs = time.time()

    vcf = VCF(
        input_file,
        samples=samples,
        threads=threads,
        gts012=True,
        strict_gt=True
    )

    w = Writer(output_file, vcf)

    for variant in vcf:
        # keep variants where at least one sample is not hom-ref
        if (variant.num_called - variant.num_hom_ref) > 0:
            w.write_record(variant)

    vcf.close()
    w.close()

    print((time.time() - start_time_secs) / 60, "minutes")

if __name__ == "__main__":
    main()

