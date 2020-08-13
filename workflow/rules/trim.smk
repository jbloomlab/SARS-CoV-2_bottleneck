### ======= Trim adaptors, poly A tails, and low quality reads ======= ###

rule trim_se:
    """
    Fast all-in-one processing of single `fastq` files. Automatic 
    adaptor trimming, low-qual base filtering, and reporting.

    Exclude reads < 70bp and cut poly X enabled (-x).
    """
    input: get_avaliable_fastqs
    output: 
        trimmed=join(config['trim_dir'], "{accession}", "{accession}.trimmed.fastq"),
        html=join(config['qc_dir'], "{accession}/fastp", "{accession}.fastp.html"),
        json=join(config['qc_dir'], "{accession}/fastp", "{accession}.fastp.json")
    conda: '../envs/trim.yml'
    log: "logs/fastp/{accession}.log"
    shell:
        """ 
        fastp \
            -i {input} \
            -l 70 \
            -x \
            --cut_tail_window_size 20 \
            -o {output.trimmed} \
            --html {output.html} \
            --json {output.json} &> {log}
        """

rule trim_pe:
    """
    Fast all-in-one processing of paired `fastq` files. Automatic 
    adaptor trimming, low-qual base filtering, and reporting.

    Exclude reads < 70bp and cut poly X enabled (-x).
    """
    input: get_avaliable_fastqs
    output:
        trimmed=[join(config['trim_dir'], "{accession}", "{accession}_1.trimmed.fastq"), 
                 join(config['trim_dir'], "{accession}", "{accession}_2.trimmed.fastq")],
        html=join(config['qc_dir'], "{accession}/fastp", "{accession}.fastp.html"),
        json=join(config['qc_dir'], "{accession}/fastp", "{accession}.fastp.json")
    conda: '../envs/trim.yml'
    log: "logs/fastp/{accession}.log"
    shell:
        """ 
        fastp \
            -i {input[0]} \
            -I {input[1]} \
            -l 70 \
            -x \
            --cut_tail_window_size 20 \
            -o {output.trimmed[0]} \
            -O {output.trimmed[1]} \
            --html {output.html} \
            --json {output.json} &> {log}
        """
