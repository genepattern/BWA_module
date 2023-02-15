use strict;
use warnings;
print "Loading files now!";
use Config;
use File::Basename;
use File::Path;
use File::Spec::Functions qw(catfile);
use IPC::Open3;
use IO::Select;
use IO::Handle;
use Getopt::Long qw(:config no_ignore_case);
use Archive::Extract;

# Prefer cmd-line archive tools for better large file support.
$Archive::Extract::PREFER_BIN=1;

my $bwa_exec = "bwa";
my $libdir = "";
my $index = "";
my $pair1 = "";
my $pair2 = "";
my $maxEditDist = "";
my $maxNumGap = "";
my $maxNumGapExt = "";
my $maxDelLen = "";
my $maxInDelLen = "";
my $seedLen = "";
my $maxSeedEditDist = "";
my $numThreads = "";
my $mismatchPen = "";
my $gapOpenPen = "";
my $gapExtenPen = "";
my $maxBestHits = "";
my $iterSearch = "";
my $trimRead = "";
my $illumina1_3 = "";
my $barcodeLen = "";
my $maxOccur = "";
my $maxInsertSize = "";
my $maxAlignments = "";
my $maxDCAlignments = "";
my $bamMapping = "";
my $outputPrefix = "";

my $result = GetOptions ("bwa:s" => \$bwa_exec,
                         "libdir:s" => \$libdir,
                         "w1:s" => \$index,
                         "p1:s" => \$pair1,
                         "p2:s" => \$pair2,
                         "n1:s" =>  \$maxEditDist,
                         "g:s" =>  \$maxNumGap,
                         "e:s" =>  \$maxNumGapExt,
                         "d:s" =>  \$maxDelLen,
                         "idl:s" =>  \$maxInDelLen,
                         "l:s" =>  \$seedLen,
                         "k:s" =>  \$maxSeedEditDist,
                         "t:s" =>  \$numThreads,
                         "m:s" =>  \$mismatchPen,
                         "o1:s" =>  \$gapOpenPen,
                         "E:s" =>  \$gapExtenPen,
                         "r:s" =>  \$maxBestHits,
                         "n2:s" =>  \$iterSearch,
                         "q:s" =>  \$trimRead,
                         "ilmn:s" =>  \$illumina1_3,
                         "B:s" =>  \$barcodeLen,
                         "u:s" =>  \$maxOccur,
                         "a:s" =>  \$maxInsertSize,
                         "x1:s" =>  \$maxAlignments,
                         "x2:s" =>  \$maxDCAlignments,
                         "bm:s"  =>  \$bamMapping,
                         "o2:s" =>  \$outputPrefix
                         );

my $custom_bwa_index_path = "";
my $read1_result = "";
my $read2_result = "";

if($pair1 eq "")
{
    wr_die("\nA read file must be specified");

}
if($outputPrefix eq "")
{
    wr_die("\nAn output prefix must be specified");
}
$outputPrefix =~ s/ /_/g; # Replace any spaces with '_' to avoid issues creating output files.

if($pair1 =~ /.*bam$/)
{
    wr_die("\nA bam mapping must be specified for a bam input file") if ($bamMapping eq "");
    wr_die("\nThe second reads.pair input file must be blank if the first is a bam file") if ($pair2 !~ /^\s*$/);
}

if($pair1 !~ /.*bam$/ && $bamMapping ne "")
{
    wr_die("\nA bam mapping should only be specified for a bam input file");
}

if($pair2 =~ /.*bam$/)
{
    wr_die("\nThe second reads.pair input file cannot be a BAM file.  For processing paired-end data, the two ends in a pair must be grouped together");
}


my $bwa_index = "";
if($index !~ /^\s*$/)
{
    my $bwa_index_path = $index;

    if($index =~ /.*zip$/)
    {
        $custom_bwa_index_path = "bwa_index/";

        my $ae = Archive::Extract->new( archive => $index );
        my $ok = $ae->extract( to => $custom_bwa_index_path ) or wr_die($ae->error);

        $bwa_index_path = $custom_bwa_index_path;
    }

    opendir(DIR, ${bwa_index_path}) or wr_die($!);

    my $bwa_index_name = "";
    while (my $file = readdir(DIR))
    {
        my $file_prefix = substr basename($file), 0, index(basename($file), '.');
        if($file !~ m/^\./)
        {
            $bwa_index_name = $file_prefix;
        }
    }
      
    closedir(DIR);
    #check if index name was found
    if($bwa_index_name eq "")
    {
        wr_die("\nUnable to detect the name of the BWA index.");
    }

    $bwa_index = catfile($bwa_index_path, $bwa_index_name);
}
else
{
    wr_die("\nA BWA index must be specified");
}

my @aln_cmd = ();
push (@aln_cmd, "bwa", "aln");

if($maxEditDist ne "")
{
    wr_die("\nmax.edit.distance must be an integer or a float") 
       unless $maxEditDist =~ /^[+-]?\d*(\.?\d+)?(e[+-]\d+)?$/;
    push (@aln_cmd, "-n", $maxEditDist);
}
if($maxNumGap ne "")
{
    wr_die("\nmax.num.gap must be an integer") unless $maxNumGap =~ /^\d+$/;
    push (@aln_cmd, "-o", $maxNumGap);
}
if($maxNumGapExt ne "")
{
    # Note: negative integer allowed
    wr_die("\nmax.gap.extension must be an integer") unless $maxNumGapExt =~ /^[+-]?\d+$/;
    push (@aln_cmd, "-e '$maxNumGapExt'");
}
if($maxDelLen ne "")
{
    wr_die("\nmax.deletion.length must be an integer") unless $maxDelLen =~ /^\d+$/;
    push (@aln_cmd, "-d", $maxDelLen);
}
if($maxInDelLen ne "")
{
    wr_die("\nmax.indel.length must be an integer") unless $maxInDelLen =~ /^\d+$/;
    push (@aln_cmd, "-i", $maxInDelLen);
}
if($seedLen ne "")
{
    wr_die("\nseed.length must be an integer") unless $seedLen =~ /^\d+$/;
    push (@aln_cmd, "-l", $seedLen);
}
if($maxSeedEditDist ne "")
{
    wr_die("\nmax.seed.edit.distance must be an integer") unless $maxSeedEditDist =~ /^\d+$/;
    push (@aln_cmd, "-k", $maxSeedEditDist);
}
if($numThreads ne "")
{
    wr_die("\nnum.threads must be an integer") unless $numThreads =~ /^\d+$/;
    push (@aln_cmd, "-t", $numThreads);
}
if($mismatchPen ne "")
{
    wr_die("\nmismatch.penalty must be an integer") unless $mismatchPen =~ /^\d+$/;
    push (@aln_cmd, "-M", $mismatchPen);
}
if($gapOpenPen ne "")
{
    wr_die("\ngap.open.penalty must be an integer") unless $gapOpenPen =~ /^\d+$/;
    push (@aln_cmd, "-O", $gapOpenPen);
}
if($gapExtenPen ne "")
{
    wr_die("\ngap.extension.penalty must be an integer") unless $gapExtenPen =~ /^\d+$/;
    push (@aln_cmd, "-E", $gapExtenPen);
}
if($maxBestHits ne "")
{
    wr_die("\nmax.best.hits must be an integer") unless $maxBestHits =~ /^\d+$/;
    push (@aln_cmd, "-R", $maxBestHits);
}
if($iterSearch eq "yes")  # Note: "yes" means "iterative search disabled"
{
    push (@aln_cmd, "-N");
}
if($trimRead ne "")
{
    wr_die("\ntrim.reads must be an integer") unless $trimRead =~ /^\d+$/;
    push (@aln_cmd, "-q", $trimRead);
}
if($illumina1_3 eq "yes")
{
    #format is Illumina 1.3+
    push (@aln_cmd, "-I");
}
if($barcodeLen ne "")
{
    wr_die("\nbarcode.length must be an integer") unless $barcodeLen =~ /^\d+$/;
    push (@aln_cmd, "-B", $barcodeLen);
}

my $bam1 = "";
my $bam2 = "";
if($bamMapping eq "single")
{
    #only use single-end reads in mapping if input file is in bam format
    $bam1 = "-b0";
}
elsif($bamMapping eq "first")
{
    #only use first read in a read pair in mapping if input file is in bam format
    $bam1 = "-b1";
}
elsif($bamMapping eq "second")
{
    #only use second read in a read pair in mapping if input file is in bam format
    $bam1 = "-b2";
}
elsif($bamMapping eq "paired")
{
    #use $pair1 as both first and second read pair in mapping if input file is in bam format
    $bam1 = "-b1";
    $bam2 = "-b2";
    $pair2 = $pair1;
}

my @sam_cmd = ();

#samse/sampe cmd
$read1_result = ${outputPrefix}."1.sai";
my $outputFile = "$outputPrefix.sam";
if($pair2 ne "") {
    push (@sam_cmd, "bwa", "sampe");
    push(@sam_cmd, "-f", $outputFile);
    if($maxOccur ne "")
    {
        wr_die("\nmax.occurrences must be an integer") unless $maxOccur =~ /^\d+$/;
        push (@sam_cmd, "-o", $maxOccur);
    }
    if($maxInsertSize ne "")
    {
        wr_die("\nmax.insert.size must be an integer") unless $maxInsertSize =~ /^\d+$/;
        push (@sam_cmd, "-a", $maxInsertSize);
    }

    if($maxAlignments ne "")
    {
        wr_die("\nmax.alignments must be an integer") unless $maxAlignments =~ /^\d+$/;
        push (@sam_cmd, "-n", $maxAlignments);
    }
    if($maxDCAlignments ne "")
    {
        wr_die("\nmax.dc.alignments must be an integer") unless $maxDCAlignments =~ /^\d+$/;
        push (@sam_cmd, "-N", $maxDCAlignments);
    }

    my @aln_cmd1 = @aln_cmd;
    push(@aln_cmd1, $bam1) if ($bam1 ne "");
    push(@aln_cmd1, "-f", $read1_result);
    runCmd(@aln_cmd1, $bwa_index, $pair1);

    my @aln_cmd2 = @aln_cmd;
    push(@aln_cmd2, $bam2) if ($bam2 ne "");
    $read2_result = ${outputPrefix}."2.sai";
    push(@aln_cmd2, "-f", $read2_result);
    runCmd(@aln_cmd2, $bwa_index, $pair2);

    push(@sam_cmd, $bwa_index, $read1_result, $read2_result, $pair1, $pair2);
}
else{
    push (@sam_cmd, "bwa", "samse");
    push(@sam_cmd, "-f", $outputFile);
    if($maxAlignments ne "")
    {
        wr_die("\nmax.alignments must be an integer") unless $maxAlignments =~ /^\d+$/;
        push (@sam_cmd, "-n", $maxAlignments);
    }
    if($maxOccur ne "")
    {
        print STDOUT "Ignoring 'max.occurrences' setting: does not apply for single-end reads.\n";
    }
    if($maxInsertSize ne "")
    {
        print STDOUT "Ignoring 'max.insert.size' setting: does not apply for single-end reads.\n";
    }
    if($maxDCAlignments ne "")
    {
        print STDOUT "Ignoring 'max.dc.alignments' setting: does not apply for single-end reads.\n";
    }

    push(@aln_cmd, $bam1) if ($bam1 ne "");
    push(@aln_cmd, "-f", $read1_result);
    runCmd(@aln_cmd, $bwa_index, $pair1);
    push(@sam_cmd, $bwa_index, $read1_result, $pair1);
}

runCmd(@sam_cmd);
verifySamOutput($outputFile);


sub verifySamOutput 
{
    my $samFile = shift;

    # Only doing the most basic check of the SAM output:
    my $haveReadLines = 0;  # Have we read any lines from the file?
    my $haveSeenHeader = 0;  # Have we seen any header lines?
    my $haveSeenContent = 0;  # Have we seen any content lines?
    open(SAMFILE, "<", $samFile) or die("Could not open sam file: $samFile");
    foreach my $line (<SAMFILE>)
    {
        if (! $haveReadLines) {$haveReadLines = 1; }
        
        if ($line =~ /^@/)
        {
            # Lines starting with '@' signify header lines.
            $haveSeenHeader = 1;
        }
        elsif ($line !~ /^\s*$/) 
        {
            # Everything else is content - we aren't going to check for valid content.
            $haveSeenContent = 1;
            last;
        }
    }
    close(SAMFILE) or warn "close failed $!";

    wr_die("An incomplete SAM file was created. Please check that a valid input file was provided.") 
        if (not ($haveReadLines && $haveSeenHeader && $haveSeenContent));
}

sub runCmd
{
    my (@cmd) = @_;
    
    #output command line to a separate file
    open (CMDLINEFILE, '>> cmdline.log');
    print CMDLINEFILE "Command: @cmd\n";

    print "About to execute: @cmd", "\n";
    
    ## The following is adapted from http://www.perlmonks.org/?node_id=150748
    my $Pin  = new IO::Handle;       $Pin->fdopen(10, "w");
    my $Pout = new IO::Handle;       $Pout->fdopen(11, "r");
    my $Perr = new IO::Handle;       $Perr->fdopen(12, "r");
    my $Proc = open3($Pin, $Pout, $Perr, @cmd);

    my $sel = IO::Select->new();
    $sel->add($Perr, $Pout);

    # Note that we are redirecting most STDERR output to STDOUT, only scanning for specific messages and
    # using these as a signal that an error truly occurred.
    while (my @ready = $sel->can_read)
    {
        foreach my $handle (@ready)
        {
            my ($count, $data);
            $count = sysread($handle, $data, 1024);
            if (not defined $count)
            {
                wr_die("Error reading BWA output");
            }
            if ($count == 0)
            {
                $sel->remove($handle);
                next;
            }
            if (fileno($handle) == fileno($Perr))
            {
                # process has printed something on standard error
                if ($data =~ m/Abort/i || $data =~ m/Not enough memory allocated!/i)
                {
                    #if child STDERR output looks like an error message, write to STDERR
                    print STDERR $data;
                }
                else
                {
                    # Otherwise just print it to the STDOUT of the parent process
                    print STDOUT $data;
                }
            }
            elsif (fileno($handle) == fileno($Pout))
            {
                # process has printed something on standard out
                print STDOUT $data;
            }
            else
            {
                # This should never happen...
                wr_die("Internal error; read from unknown filehandle");
            }
        }
    }
    
    close($Perr) or warn "close failed $!";
    close($Pout) or warn "close failed $!";

    waitpid($Proc, 0);
    if ($?) { wr_die("An error occurred when executing the command."); }
}

# More user-friendly errors than Perl's 'die' for use in the context of a wrapper script.
sub wr_die {
    print STDERR @_, "\n\n";
    exit(1);
}

END
{
    #delete custom BWA index files and .sai files in the current directory
    if ($custom_bwa_index_path ne "") {
        rmtree($custom_bwa_index_path);
    }
    if ($read1_result ne "") {
        unlink("$read1_result");
    }
    if ($read2_result ne "") {
        unlink("$read2_result");
    }
}
