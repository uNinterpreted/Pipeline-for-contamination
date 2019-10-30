#!/usr/bin/perl -w

# This script intends to estimate three types of contamination in tumor-normal pair data.
# Designed for whole exome sequencing data.[Parameters should be adjusted].

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use FindBin qw($Bin); # $Bin = path/to/this/script;$Script = basename of this script
use lib $Bin; # locate to the relative path 
use Getopt::Long;
use SoftwareConfig;

our ($normal,$tumor,$outdir,$help,$pon);

GetOptions(
	'Normal|N=s' => \$normal,
	'Tumor|T=s' => \$tumor,
	'PanelOfNormals|PoN=s' => \$pon,
	'Outdir|O=s' => \$outdir,
	'help|h=s' => \$help,
);

my $usage="
contamination-estimated.pl:
	This pipeline can estimate the three types of contamination in tumor-normal pair data, including 
cross-species (by Centrifuge), cross-sample (by GATK CalculateContamination), tumor in normal (by DeTiN).

Usage: perl $0 -n R1,R2 -t R1,R2 -o /path/to/output 

Options:
	-N|--Normal               Input files (.fq) of normal samples, pair-end reads are sperated by comma [required].
	-T|--Tumor                Input files (.fq) of tumor samples, pair-end reads are sperated by comma [required].
	-PoN|--PanelOfNormals     The list of germline samples (recal.bam) to create somatic panel of normals, one sample per line.
				   If no list provided, the normal sample will be used instead (that is, no prefilter variants) [optional].
	-O|--Outdir               Path to output directory, SampleID is recommended [required].
	-h|--help                 Print this help information.
";

if($help){
	print $usage;
}

if( !defined $normal || !defined $tumor || !defined $outdir){
	print "ERROR: wrong input file or output dir.";
	print $usage;
	exit;
}


# read config file to obtain the absolute path of softwares and resource files
my $config="$Bin/software_files.config";

our $gnome=parse_config($config,"hg19_genome");
our $cosmic=parse_config($config,"hg19_cosmic");
our $dbsnp=parse_config($config,"hg19_dbsnp");
our $known_indels=parse_config($config,"hg19_indels");
our $germ_rsc=parse_config($config,"germline_resource");
our $small_comm=parse_config($config,"small_exac_common");
our $snp_inter=parse_config($config,"snp_interval");
our $target_bed=parse_config($config,"target_bed");
our $exac_sites=parse_config($config,"exac_sites");
# our $inter_list=  #consider to provide the interval list of targeted captured sequencing or whole exome sequencing


our $java_7=parse_config($config,"java_7");
our $java_8=parse_config($config,"java_8");

our $fastqc=parse_config($config,"fastqc");
# our $umi_tools=

our $bwa=parse_config($config,"bwa");
our $samtools=parse_config($config,"samtools");
our $picard=parse_config($config,"picard");

our $gatk4=parse_config($config,"gatk4");
our $gatk_protected=parse_config($config,"gatk_protected");

our $centrifuge=parse_config($config,"centrifuge");
our $mutect=parse_config($config,"MuTect");
our $detin=parse_config($config,"deTiN");


my $cmd="perl $0 ".join(" ",@ARGV);
our $program_start_time=`date "+%Y-%m-%d_%H:%M:%S"`;
our ($normal_R1,$normal_R2)=(split/,/,$normal);
our ($tumor_R1,$tumor_R2)=(split/,/,$tumor);
&checkdir($outdir);
$outdir=abs_path($outdir); # get the absolute path of all dirs and files
our $SampleName=basename($outdir); # use SampleID as the identifier
$normal_R1=abs_path($normal_R1);
$normal_R2=abs_path($normal_R2);
$tumor_R1=abs_path($tumor_R1);
$tumor_R2=abs_path($tumor_R2);
&run_log($cmd,"CMD");
chdir($outdir);

# Step1: quality control
# the process of trim the UMI and other adpaters is not included here.
&checkdir("1_quality");
chdir("1_quality");
$cmd="$fastqc -t 4 -o $outdir/1_quality $normal_R1 $normal_R2\n";
$cmd.="$fastqc -t 4 -o $outdir/1_quality $tumor_R1 $tumor_R2\n";

&monitor_jobs($cmd,"STEP1_quality_control");
chdir($outdir);


# Step2: Centrifuge estimate the cross-species contamination
&checkdir("2_cross-species");
chdir("2_cross-species");
my $index=dirname($centrifuge);
# index file to be update
$cmd="$centrifuge -p 4 --host-taxids 9606 -x $index/indices/p_compressed+h+v/p_compressed+h+v  -1 $normal_R1 -2 $normal_R2 -S normal.classification --report-file normal.report\n";
$cmd.="$centrifuge -p 4 --host-taxids 9606 -x $index/indices/p_compressed+h+v/p_compressed+h+v  -1 $tumor_R1 -2 $tumor_R2 -S tumor.classification --report-file tumor.report\n";

&monitor_jobs($cmd,"STEP2_Cross_species");
chdir($outdir);


# Step3: pre-process for variant discovery
&checkdir("3_pre-process");
chdir("3_pre-process");
undef $cmd;

foreach my $sample("tumor","normal"){
	if ( ! defined $cmd ){
		$cmd="$bwa mem -t 4 -M -R '\@RG\tID:$SampleName\tLB:library_name\tSM:${SampleName}_${sample}\tPL:ILLUMINA' $gnome ${sample}_R1 ${sample}_R2 |gzip > ${sample}.sam.gz\n";
	}else{
		$cmd.="$bwa mem -t 4 -M -R '\@RG\tID:$SampleName\tLB:library_name\tSM:${SampleName}_${sample}\tPL:ILLUMINA' $gnome ${sample}_R1 ${sample}_R2 |gzip > ${sample}.sam.gz\n";
	}
	$cmd.="$java_8 -Xmx8g -jar  $picard SortSam I=${sample}.sam.gz O=${sample}_sorted.bam CREATE_INDEX=true SO=coordinate\n";
	$cmd.="$java_8 -Xmx8g -jar $picard MarkDuplicates I=${sample}_sorted.bam O=${sample}_dedup.bam M=${sample}_dedup.metrics CREATE_INDEX=true\n";
	$cmd.="$gatk4 --java-options \"-Xmx8g\" BaseRecalibrator -I ${sample}_dedup.bam -R $gnome --known-sites $dbsnp --known-sites $known_indels -O ${sample}_recal_data.table\n";
	$cmd.="$gatk4 --java-options \"-Xmx8g\" ApplyBQSR -R $gnome -I ${sample}_dedup.bam --bqsr-recal-file ${sample}_recal_data.table -O ${sample}_recal.bam\n";
}

&monitor_jobs($cmd,"STEP3_Pre_process");
chdir($outdir);


# Step4: run GATK Mutect2 to call somatic variation and estimate cross-sample contamination
&checkdir("4_Mutect2");
chdir("4_Mutect2");

if ( defined $pon ){
	$pon=abs_path($pon);
	my $cmd;
	my $vcfs;
	my @list=`cat $pon`;
	chomp(@list);
	foreach my $bam(@list){
		my $flag=`$samtools view -H $bam |grep '\@RG' |head -1 |awk -F ' |:' '{for (i=1;i<=NF;i++) {if(\$i=="SM") print \$(i+1)}}'`;
		chomp($flag);
		if ( ! defined $cmd ){
			$cmd="$gatk4 --java-options \"-Xmx8g\" Mutect2 -R $gnome -I $bam -tumor $flag --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -O $flag.vcf.gz\n";
			$vcfs="-vcfs $flag.vcf.gz ";			
		}else{
			$cmd.="$gatk4 --java-options \"-Xmx8g\" Mutect2 -R $gnome -I $bam -tumor $flag --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -O $flag.vcf.gz\n";
			$vcfs.="-vcfs $flag.vcf.gz ";
		}
	}
	$cmd.="$gatk4 --java-options \"-Xmx8g\" CreateSomaticPanelOfNormals $vcfs -O ${SampleName}.PoN.vcf.gz\n";
}else{
	$cmd="$gatk4 --java-options \"-Xmx8g\" Mutect2 -R $gnome -I $outdir/3_pre-process/normal_recal.bam -tumor ${SampleName}_normal --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -O ${SampleName}_normal.vcf.gz\n";
	$cmd.="$gatk4 --java-options \"-Xmx8g\" CreateSomaticPanelOfNormals -vcfs ${SampleName}_normal.vcf.gz  -O ${SampleName}.PoN.vcf.gz\n";
}

&monitor_jobs($cmd,"STEP4_Create_PoN");
	
$cmd="$gatk4 --java-options \"-Xmx8g\" Mutect2 -R $gnome -I $outdir/3_pre-process/tumor_recal.bam -I $outdir/3_pre-process/normal_recal.bam -tumor ${SampleName}_tumor -normal ${SampleName}_normal -pon ${SampleName}.PoN.vcf.gz --germline-resource $germ_rsc --af-of-alleles-not-in-resource 0.0000025 --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -O somatic_${SampleName}.vcf.gz -bamout tumor_normal.bam\n";

&monitor_jobs($cmd,"STEP4_Run_Mutect2");

$cmd="$gatk4  GetPileupSummaries -I $outdir/3_pre-process/tumor_recal.bam -V $small_comm -O tumor_getpileupsummaries.table \n";
$cmd.="$gatk4 GetPileupSummaries -I $outdir/3_pre-process/normal_recal.bam -V $small_comm -O normal_getpileupsummaries.table \n";
$cmd.="$gatk4 CalculateContamination -I tumor_getpileupsummaries.table -matched normal_getpileupsummaries.table -O ${SampleName}_calculatecontamination.table \n";

&monitor_jobs($cmd,"STEP4_Calculate_Conta");

$cmd="$gatk4 FilterMutectCalls -V somatic_${SampleName}.vcf.gz --contamination-table ${SampleName}_calculatecontamination.table -O somatic_${SampleName}_oncefiltered.vcf.gz \n";
$cmd.="$gatk4 CollectSequencingArtifactMetrics -I $outdir/3_pre-process/tumor_recal.bam -O tumor_${SampleName}_artifact --FILE_EXTENSION \".txt\" -R $gnome \n";
$cmd.="$gatk4 FilterByOrientationBias -AM G/T -AM C/T -V somatic_${SampleName}_oncefiltered.vcf.gz -P tumor_${SampleName}_artifact.pre_adapter_detail_metrics.txt -O somatic_${SampleName}_twicefiltered.vcf.gz \n";

&monitor_jobs($cmd,"STEP4_Filter_Variants");
chdir($outdir);


# Step5: prepare input data for deTiN and run deTiN
&checkdir("5_deTiN");
chdir("5_deTiN");

$cmd="$java_7 -Xmx8g -jar $mutect --analysis_type MuTect --reference_sequence $gnome --cosmic $cosmic --dbsnp $dbsnp --input_file:normal $outdir/3_pre-process/normal_recal.bam --input_file:tumor $outdir/3_pre-process/tumor_recal.bam --out ${SampleName}.call_stats.out --coverage_file ${SampleName}.coverage.wig.txt --vcf ${SampleName}.vcf \n";
$cmd.="$java_8 -jar $gatk_protected GetBayesianHetCoverage --reference $gnome --snpIntervals $snp_inter --tumor $outdir/3_pre-process/tumor_recal.bam --tumorHets ${SampleName}_tumor.het --normal $outdir/3_pre-process/normal_recal.bam  --normalHets ${SampleName}_normal.het \n";

&monitor_jobs($cmd,"STEP5_Mutect1_HetCover");


if ( defined $pon ){
	$pon=abs_path($pon);
	my $cmd;
	my $povs;
	my @list=`cat $pon`;
	chomp(@list);
	foreach my $bam(@list){
		my $flag=`$samtools view -H $bam |grep '\@RG' |head -1 |awk -F '\t|:' '{for (i=1;i<=NF;i++) {if(\$i=="SM") print \$(i+1)}}'`;
		chomp($flag);
		if ( ! defined $cmd ){
			$cmd="$java_8 -jar $gatk_protected CalculateTargetCoverage --input $bam --transform PCOV --groupBy SAMPLE --targetInformationColumns FULL --keepduplicatereads true --output ${flag}_normal.pov.file --targets $target_bed \n";
			$povs="--input ${flag}_normal.pov.file ";
		}else{
			$cmd.="$java_8 -jar $gatk_protected CalculateTargetCoverage --input $bam --transform PCOV --groupBy SAMPLE --targetInformationColumns FULL --keepduplicatereads true --output ${flag}_normal.pov.file --targets $target_bed \n";
			$povs.="--input ${flag}_normal.pov.file ";
		}
	}
	$cmd.="$java_8 -jar $gatk_protected CombineReadCounts $povs -O combined-normal.tsv \n";
}else{
	$cmd="$java_8 -jar $gatk_protected CalculateTargetCoverage --input $outdir/3_pre-process/normal_recal.bam --transform PCOV --groupBy SAMPLE --targetInformationColumns FULL --keepduplicatereads true --output ${SampleName}_normal.pov.file --targets $target_bed \n";
	$cmd.="$java_8 -jar $gatk_protected CombineReadCounts --input ${SampleName}_normal.pov.file  -O combined-normal.tsv \n";
}
$cmd.="$java_8 -jar $gatk_protected CreatePanelOfNormals -I combined-normal.tsv -O normals.cnv.pon -noQC --minimumTargetFactorPercentileThreshold 5 --disableSpark \n";

&monitor_jobs($cmd,"STEP5_CNV_PoN");


$cmd="$java_8 -jar $gatk_protected CalculateTargetCoverage --input $outdir/3_pre-process/tumor_recal.bam --transform PCOV --groupBy SAMPLE --targetInformationColumns FULL --keepduplicatereads true --output ${SampleName}_tumor.pov.file --targets $target_bed \n";
$cmd.="$java_8 -jar $gatk_protected NormalizeSomaticReadCounts -I ${SampleName}_tumor.pov.file -PON normals.cnv.pon -PTN ${SampleName}_tumor.ptn.tsv -TN ${SampleName}_tumor.tn.tsv \n";
$cmd.="$java_8 -jar $gatk_protected PerformSegmentation -TN ${SampleName}_tumor.tn.tsv -O ${SampleName}_tumor.seg -LOG \n";
$cmd.="$java_8 -jar $gatk_protected CallSegment -TN ${SampleName}_tumor.tn.tsv -S ${SampleName}_tumor.seg -O ${SampleName}_tumor.called.seg \n";

&monitor_jobs($cmd,"STEP5_GATK_CNV");


$cmd="$java_8 -Xmx8g -jar $gatk_protected AllelicCNV --tumorHets ${SampleName}_tumor.het --tangentNormalized ${SampleName}_tumor.tn.tsv --segment ${SampleName}_tumor.called.seg --outputPrefix ${SampleName}_tumor \n ";
$cmd.="$java_8 -Xmx8g -jar $gatk_protected CallCNLoHAndSplits --tumorHets ${SampleName}_tumor.het --segments ${SampleName}_tumor-sim-final.seg --tangentNormalized ${SampleName}_tumor.tn.tsv --outputDir $outdir/5_deTiN/ACNV_results --rhoThreshold 0.2 --numIterations 10 --sparkMaster local[*] \n";

&monitor_jobs($cmd,"STEP5_ACNV");


&checkdir("deTiN_results");
$cmd="$detin --mutation_data_path ${SampleName}.call_stats.out --cn_data_path $outdir/5_deTiN/ACNV_results/${SampleName}_tumor-sim-final.acs.seg --tumor_het_data_path ${SampleName}_tumor.het --normal_het_data_path ${SampleName}_normal.het --exac_data_path $exac_sites --output_dir $outdir/5_deTiN/deTiN_results --output_name ${SampleName} \n";

&monitor_jobs($cmd,"STEP5_ACNV");
chdir($outdir);



sub monitor_jobs{
        my $cmd=shift;
        my $step=shift;
        &run_log($step,"START");
        my @unit_cmd=split /\n/,$cmd;
        foreach my $unit_cmd (@unit_cmd){
		next if /^\s|^$/;
		$unit_cmd.="1 >> $step.out 2 >> $step.err";
                my $check=system ($unit_cmd);
                &run_log($unit_cmd,"CMD");
                if($check!=0){
                        &run_log($step,"ERROR");
                        &run_log($unit_cmd,"ERROR");
                        die "$step error\n";
                }
        }
        &run_log($step,"DONE");
}


sub checkdir{
	my $dir=shift;
	if(!-d $dir){
		system("mkdir $dir");
	}else{
		die "Directory $dir existed.\n";
	} 
}

sub run_log{
	my $cmd=shift;
	my $record=shift;
	open OUT,">>$outdir/EstCon_$program_start_time.log";
	my $date=`date "+%Y-%m-%d %H:%M:%S"`;
	print OUT "[$record-$date] $cmd.\n";
}






