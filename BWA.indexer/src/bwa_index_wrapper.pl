use strict;
use warnings;

use Config;
use File::Basename;
use File::Path;
use IPC::Open3;
use IO::Select;
use IO::Handle;
use Getopt::Long;
use IO::Compress::Zip qw(:all);

my $bwa_exec = "";
my $libdir = "";
my $inputFile = "";
my $alg = "";
my $output = "";

my $result = GetOptions ("bwa:s" => \$bwa_exec,
                         "libdir:s" => \$libdir,
                         "f:s" => \$inputFile,
                         "a:s" => \$alg,
                         "o:s" => \$output
                         );

$output =~ s/ /_/g; # Replace any spaces with '_' to avoid issues creating output files.

my @cmd = ();
push (@cmd, "bwa", "index");
push (@cmd, "-a", $alg, "-p", $output, $inputFile);

#output command line to a separate file
open (CMDLINEFILE, '> cmdline.log');
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
    

# Zip the results to form an index bundle
zip [ glob("$output*") ] => "$output.zip"
    or die "Cannot create zip file: $ZipError" ;

# More user-friendly errors than Perl's 'die' for use in the context of a wrapper script.
sub wr_die {
    print STDERR @_, "\n\n";
    exit(1);
}
