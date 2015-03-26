use strict;
use warnings;

use Config;
use File::Basename;
use File::Path;
use File::Spec::Functions qw(catfile);
use IPC::Open3;
use IO::Select;
use IO::Handle;
use Getopt::Long;
use Archive::Extract;

# Prefer cmd-line archive tools for better large file support.
$Archive::Extract::PREFER_BIN=1;

my $bwa_exec = "";
my $libdir = "";
my $index = "";
my $read_file = "";
my $mate_read_file = "";
my $mScore = "";
my $mismatchPen = "";
my $gapOpenPen = "";
my $gapExtenPen = "";
my $numThreads = "";
my $bandWidth = "";
my $mThresh = "";
my $cThresh = "";
my $zHeur = "";
my $mInitSeed = "";
my $minSeed = "";
my $outputPrefix = "";

my $result = GetOptions ("bwa=s" => \$bwa_exec,
                         "libdir=s" => \$libdir,
                         "w1:s" => \$index,
                         "p=s"  => \$read_file,
                         "m:s"  => \$mate_read_file,
                         "a:s"  => \$mScore,
                         "b:s"  => \$mismatchPen,
                         "q:s"  => \$gapOpenPen,
                         "r:s"  => \$gapExtenPen,
                         "t1:s" => \$numThreads,
                         "w:s"  => \$bandWidth,
                         "t2:s" => \$mThresh,
                         "c:s"  => \$cThresh,
                         "z:s"  => \$zHeur,
                         "s:s"  => \$mInitSeed,
                         "n:s"  => \$minSeed,
                         "o:s" =>  \$outputPrefix
                         );

# Flag to signal that the input was invalid and should be removed on exit.
my $invalid_file = "false";

my $custom_bwa_index_path = "";

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

if($read_file eq "")
{
    print STDERR "\nA read file must be specified";
    exit(1);

}
if($outputPrefix eq "")
{
    print STDERR "\nAn output prefix must be specified";
    exit(1);

}

my @bwasw_cmd = ();
push (@bwasw_cmd, $bwa_exec, "bwasw");

if($mScore ne "")
{
    wr_die("\nmatch.score must be an integer") unless $mScore =~ /^\d+$/;
    push (@bwasw_cmd, "-a", $mScore);
}
if($mismatchPen ne "")
{
    wr_die("\nmismatch.penalty must be an integer") unless $mismatchPen =~ /^\d+$/;
    push (@bwasw_cmd, "-b", $mismatchPen);
}
if($gapOpenPen ne "")
{
    wr_die("\ngap.open.penalty must be an integer") unless $gapOpenPen =~ /^\d+$/;
    push (@bwasw_cmd, "-q", $gapOpenPen);
}
if($gapExtenPen ne "")
{
    wr_die("\ngap.extension.penalty must be an integer") unless $gapExtenPen =~ /^\d+$/;
    push (@bwasw_cmd, "-r", $gapExtenPen);
}
if($numThreads ne "")
{
    wr_die("\nnum.threads must be an integer") unless $numThreads =~ /^\d+$/;
    push (@bwasw_cmd, "-t", $numThreads);
}
if($bandWidth ne "")
{
    wr_die("\nband.width must be an integer") unless $bandWidth =~ /^\d+$/;
    push (@bwasw_cmd, "-w", $bandWidth);
}
if($mThresh ne "")
{
    wr_die("\nmin.score.threshold must be an integer") unless $mThresh =~ /^\d+$/;
    push (@bwasw_cmd, "-T", $mThresh);
}
if($cThresh ne "")
{
    wr_die("\nthreshold.coefficient must be an integer or a float") 
       unless $cThresh =~ /^[+-]?\d*(\.?\d+)?(e[+-]\d+)?$/;
    push (@bwasw_cmd, "-c", $cThresh);
}
if($zHeur ne "")
{
    wr_die("\nz.best.heuristics must be an integer") unless $zHeur =~ /^\d+$/;
    push (@bwasw_cmd, "-z", $zHeur);
}
if($mInitSeed ne "")
{
    wr_die("\nmax.sa.interval.size must be an integer") unless $mInitSeed =~ /^\d+$/;
    push (@bwasw_cmd, "-s", $mInitSeed);
}
if($minSeed ne "")
{
    wr_die("\nmin.num.seeds must be an integer") unless $minSeed =~ /^\d+$/;
    push (@bwasw_cmd, "-N", $minSeed);
}

my $outputFile = $outputPrefix.".sam";
$outputFile =~ s/ /_/g; # Replace any spaces with '_' to avoid issues creating output files.
push (@bwasw_cmd, "-f", $outputFile);
push (@bwasw_cmd, $bwa_index, "$read_file");

# Add the mate FASTQ if one was provided. 
$mate_read_file =~ s/^\s+//;
$mate_read_file =~ s/\s+$//;
if ($mate_read_file ne "") {
    push (@bwasw_cmd, "$mate_read_file");
}

#output command line to a separate file
open (CMDLINEFILE, '> cmdline.log');
print CMDLINEFILE "Command: @bwasw_cmd\n";

print "About to execute: @bwasw_cmd", "\n";

## The following is adapted from http://www.perlmonks.org/?node_id=150748
my $Pin  = new IO::Handle;       $Pin->fdopen(10, "w");
my $Pout = new IO::Handle;       $Pout->fdopen(11, "r");
my $Perr = new IO::Handle;       $Perr->fdopen(12, "r");
my $Proc = open3($Pin, $Pout, $Perr, @bwasw_cmd);

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
            #check for possible invalid input
            if ($data =~ m/read 0 sequences/i)
            {
                $invalid_file = "true";
                print STDERR $data;
                print STDERR "Please check that the input file is a valid FASTA and FASTQ file";
            }

            if ($data =~ m/Abort/i || $data =~ m/Not enough memory allocated!/i)
            {
                #if child STDERR output looks like an error message, write to STDERR
                print STDERR $data;
            }
            else
            {
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

# More user-friendly errors than Perl's 'die' for use in the context of a wrapper script.
sub wr_die {
    print STDERR @_, "\n\n";
    exit(1);
}

END
{
    #delete directory containing custom BWA index files from current directory
    if ($custom_bwa_index_path ne "") {
        rmtree($custom_bwa_index_path);
    }

    #check if invalid input file
    if($invalid_file eq "true")
    {
        if ($outputFile ne "") {
            unlink($outputFile);
        }
    }
}

