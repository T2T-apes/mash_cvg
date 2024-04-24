#!vim set syntax=python


_haps = ["chm13#1",
        "mGorGor1#M",
        "mGorGor1#P",
        "mPanPan1#M",
        "mPanPan1#P",
        "mPanTro3#1",
        "mPanTro3#2",
        "mPonAbe1#1",
        "mPonAbe1#2",
        "mPonPyg2#1",
        "mPonPyg2#2",
        "mSymSyn1#1",
        "mSymSyn1#2"]
haps = ["chm13#1"]

def get_inputs(wildcards):
    inputs = []
    for hap in haps:
        inputs.append("output/{hap}.filt_0_aln.pdf".format(hap=hap))
        inputs.append("output/{hap}.filt_10e6_aln.pdf".format(hap=hap))
        #inputs.append("output/{hap}/{hap}_cvg.filt_10e6_aln.parsed.bed.gz".format(hap=hap))
        #inputs.append("output/{hap}/{hap}_cvg.filt_10e6_aln.sorted.bed.gz".format(hap=hap))
        #inputs.append("output/{hap}/{hap}_cvg.filt_10e6_aln.sorted.bigwig".format(hap=hap))
        #inputs.append("output/{hap}/{hap}_contigs".format(hap=hap))
    return inputs


"""
/Users/petersudmant/Documents/science/sudmantlab/projects/primate_T2T/great_ape_T2T/wgatools/target/release/wgatools 

/Users/petersudmant/Documents/science/sudmantlab/projects/primate_T2T/great_ape_T2T/wgatools/target/release/wgatools pafcov  ~/tmp/ape_T2T/test_paf_chr20_10.aln.paf | ~/Documents/science/programs/UCSC/macOSX.x86_64/bedGraphPack stdin stdout
"""
#/Users/petersudmant/Documents/science/sudmantlab/projects/primate_T2T/great_ape_T2T/wgatools/target/release/wgatools

wgatools = "/Users/petersudmant/Documents/science/sudmantlab/projects/primate_T2T/great_ape_T2T/wgatools/target/release/wgatools"
bedGraphPack = "~/Documents/science/programs/UCSC/macOSX.x86_64/bedGraphPack"
bedGraphToBigWig = "~/Documents/science/programs/UCSC/macOSX.x86_64/bedGraphToBigWig "

rule all:
    input:
        get_inputs

rule contig_file:
    input:
        "output/{hap}/{hap}_filt_10e6_aln.paf"
    output:
        "output/{hap}/{hap}_contigs",
        "output/{hap}/{hap}_contigs_parsed"
    run:
        cmd  = """cat {input[0]} | awk '{{{{print $6,$7}}}}' | sort | uniq > {output[0]}"""
        #cmd = cmd.format(fn = input[0], 
        #                 fout = output[0])
        print(cmd)
        shell(cmd)
        cmd  = """cat {output[0]} | awk -F '#' '{{{{print $(NF)}}}}' > {output[1]}"""
        #cmd = cmd.format(fn = output[0], fout=output[1])
        print(cmd)
        shell(cmd)

rule make_wig:
    input:
        "output/{hap}/{hap}_cvg.filt_10e6_aln.sorted.bed",
        "output/{hap}/{hap}_contigs_parsed"
    output:
        "output/{hap}/{hap}_cvg.filt_10e6_aln.sorted.bigwig"
    run:
        shell("~/Documents/science/programs/UCSC/macOSX.x86_64/bedGraphToBigWig {input[0]} {input[1]} {output[0]}")
               
rule unzipped_sorted_bed:
    input:
        "output/{hap}/{hap}_cvg.filt_10e6_aln.sorted.bed.gz"
    output:
        "output/{hap}/{hap}_cvg.filt_10e6_aln.sorted.bed"
    run:
        shell("pigz -dc {input[0]} >{output[0]}")

rule sorted_bed:
    input:
        "output/{hap}/{hap}_cvg.filt_10e6_aln.bed.gz",
        "output/{hap}/{hap}_contigs_parsed"
    output:
        "output/{hap}/{hap}_cvg.filt_10e6_aln.sorted.bed.gz"
    run:
        cmd = ("pigz -dc {input[0]} | "
               "awk -F '#' -v OFS='\t' '{{{{print $3}}}}' | "
               "sort -k1,1 -k2,2n | "
               "pigz >{output[0]} ")
               #"~/Documents/science/programs/UCSC/macOSX.x86_64/bedGraphToBigWig stdin {input[1]} {output[0]}")
        shell(cmd)

rule cvg_bed:
    input:
        "output/{hap}/{hap}_filt_10e6_aln.paf"
    output:
        "output/{hap}/{hap}_cvg.filt_10e6_aln.bed.gz"
    run:
        cmd = "{wgatools} pafcov {input} | {bedGraphPack} stdin stdout | pigz >{output}"
        cmd = cmd.format(wgatools = wgatools, 
                         input = input[0], 
                         bedGraphPack=bedGraphPack,
                         output = output[0])
        shell(cmd)

rule combine:
    input:
        "output/{hap}/{hap}_filt_10e6_aln.paf"
    output:
        "output/{hap}.{filt}.pdf"
    run:
        cmd = "pdftk ./output/{hap}/pdfs/*{filt}*.pdf cat output ./output/{hap}.{filt}.pdf"
        cmd = cmd.format(hap=wildcards.hap,filt=wildcards.filt)
        shell(cmd)

rule make_plots:
    output:
        "output/{hap}/{hap}_filt_10e6_aln.paf"
    run:
        #fn_input = "/Users/petersudmant/tmp/ape_T2T/primates13.20231122_wfmash-v0.12.2-1-g545db3e-map/{hap}.map.paf"
        #fn_input = "/Users/petersudmant/tmp/ape_T2T/primates16.20231205_wfmash-v0.12.5/{hap}.map.paf"
        #fn_input = "/Users/petersudmant/tmp/ape_T2T/primates16.20231205_wfmash-v0.12.5/{hap}.aln.paf"
        fn_input = "/Users/petersudmant/tmp/ape_T2T/primates16.2024-03-25/{hap}.p70.aln.paf"
        fn_input = fn_input.format(hap=wildcards.hap)
        #outdir = "/Users/petersudmant/Documents/science/sudmantlab/projects/primate_T2T/great_ape_T2T/plot_wfmash_alns/outdir"
        outdir = "output/{hap}".format(hap=wildcards.hap)
        shell("mkdir -p {outdir}/pdfs".format(outdir=outdir))
        cmd = "Rscript plot_script.R {fn_input} {outdir} {hap}".format(fn_input = fn_input, outdir=outdir, hap=wildcards.hap)
        print(cmd)
        shell(cmd)
            



