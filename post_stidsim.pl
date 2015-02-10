#!/usr/bin/perl -w

use strict;

#GSL
use Math::GSL::RNG qw/:all/;
use Math::GSL::Randist qw/:all/;
use Math::GSL::CDF qw/:all/;
my $rng = Math::GSL::RNG->new();

#SQLite
use DBI;

#For parallel
my $n_threads=0; #Maximum by default

my @files;
my $curr_sp=-1;
my $sp_counter=0;
my $sequence_counter=1;
my $filehandwrite;
my $filehandread;
my $models;
my $settings;
my $partitions;
my $evolves;
my $backup=$/;
my $content;
my $locus;
my $n_digits=0;
my $trees;
my $itree;
#Truncated Gamma

sub rlrtrunc_exp
{
	my ($u,$alpha,$a,$b)=@_;
	my $p_a=gsl_cdf_exponential_P($a,1/$alpha);
	my $p_b=gsl_cdf_exponential_P($b,1/$alpha);
	
	
	return gsl_cdf_exponential_Pinv($p_a+($u*($p_b-$p_a)),$alpha);
}

sub rltrunc_exp
{
	my ($u,$alpha,$a)=@_;
	my $p_a=gsl_cdf_exponential_P($a,1/$alpha);
	my $p_b=1;
	
	
	return gsl_cdf_exponential_Pinv($p_a+($u*($p_b-$p_a)),$alpha);
}

#GTR Model R-Matrix
my ($a,$b,$c,$d,$e,$f);
#GTR Frequencies
my ($T,$C,$A,$G);
my $f_total=0;

my $length=0;

#Species_dependent parameters
my $shape_seqlength=0;
my $logscale_seqlength=0;
my $alpha=1;

if ($#ARGV != 1)
{
	die "Incorrect number of parameters, Usage: script.pl directory numberofcores\n";
}

(my $w_dir,$n_threads)=@ARGV;
chdir($w_dir) or die "Error changing the working dir to $w_dir\n";
$w_dir=~m/([^\/]*).?$/;
$w_dir=$1;

#Initializing SQlite
my $SQL_name="$w_dir.db";
my $SQL= DBI->connect("dbi:SQLite:dbname=$SQL_name","","",{ RaiseError => 1, AutoCommit => 0}) or die $DBI::errstr; #DB's Handler
$SQL->do("DROP TABLE IF EXISTS Sequences");
$SQL->do("PRAGMA foreign_keys=yes");
$SQL->do("CREATE TABLE Sequences(AID INTEGER PRIMARY KEY, SID INTEGER, n_ltree INTEGER, n_gtree INTEGER, length INTEGER, FOREIGN KEY(SID) REFERENCES Species_Trees(SID))");
$SQL->commit();						

opendir (my $dirs_handler, ".");

my @dirs = grep {-d "./$_" && ! /^\.{1,2}$/} readdir($dirs_handler);

foreach my  $dir (@dirs)
{
	$sp_counter=int($dir);
	print "\n\n\nTreating gene trees from the replicate $sp_counter\n";
	
	#Inside newdir
	chdir($dir) or die "Error changing the working dir\n";

	#Sampling Species_specific parameters and writing them into the SQLite DB.
	#$shape_seqlength=gsl_ran_flat($rng->raw(),5.7,7.3); # Average length drawn from ln(5.7)=300 to ln(7.3)=1500
	#$logscale_seqlength=gsl_ran_flat($rng->raw(),0.0,0.3);


	#$SQL->do("UPDATE Species_Trees SET Shape_seqlength=$shape_seqlength WHERE SID==$sp_counter");
	#$SQL->do("UPDATE Species_Trees SET Logscale_seqlength=$logscale_seqlength WHERE SID==$sp_counter");
	
	#INDELIBLE
	print "\t\nGenerating the INDELIBLE control.txt file\n";
	open($filehandwrite,">"."control.txt") or die "Error opening the file\n";
	
	@files=<g_trees*.trees>;
	
	$models="[TYPE] NUCLEOTIDE 1\n";
	$settings="[SETTINGS]\n[randomseed] 2478\n[fileperrep] FALSE\n";
	$trees='';
	$partitions='';
	$evolves='[EVOLVE] ';
	
	foreach my $file (@files)
	{
		open($filehandread,$file) or die "Error opening the file $file\n";
		$file=~m/g_trees(\d*)\.trees/;
		$n_digits=length($1);
		$locus=int($1);
		$/="";
		$itree=<$filehandread>;
		chomp($itree);
		close($filehandread);

		#Sampling exponential(1.2), truncated at .1
		$alpha=0;
		while ($alpha < 0.1) {
			$alpha=-log(gsl_rng_uniform_pos($rng->raw()))/1.2;
		}

		#Sampling Dirichlet (36 26 28 32) for frequencies (to avoid alignments without one of the bases)
		$A=gsl_ran_gamma($rng->raw(),36,1);
		$C=gsl_ran_gamma($rng->raw(),26,1);
		$G=gsl_ran_gamma($rng->raw(),28,1);
		$T=gsl_ran_gamma($rng->raw(),32,1);
		$f_total=$A+$C+$T+$G;
		$A/=$f_total;
		$C/=$f_total;
		$T/=$f_total;
		$G/=$f_total;
		
		#Sampling Dirichlet ( 16   3   5   5   6  15) for frequencies
		$a=gsl_ran_gamma($rng->raw(),16,1);
		$b=gsl_ran_gamma($rng->raw(),3,1);
		$c=gsl_ran_gamma($rng->raw(),5,1);
		$d=gsl_ran_gamma($rng->raw(),5,1);
		$e=gsl_ran_gamma($rng->raw(),6,1);
		$f=gsl_ran_gamma($rng->raw(),15,1);
		$f_total=$a+$b+$c+$d+$d+$e+$f;
		$a/=$f_total;
		$b/=$f_total;
		$c/=$f_total;
		$d/=$f_total;
		$e/=$f_total;
		$f/=$f_total;
		$a/=$f;
		$b/=$f;
		$c/=$f;
		$d/=$f;
		$e/=$f;
		$f/=$f;
		
		#Sampling Sequence length
		$length= int(gsl_ran_lognormal($rng->raw(),$shape_seqlength,$logscale_seqlength));

		$models.=sprintf("\[MODEL] GTR%.*d\n\t[submodel]  GTR %f %f %f %f %f\n\t[statefreq] %f %f %f %f\n\t[rates] 0 %f 0\n",$n_digits,$locus,$a,$b,$c,$d,$e,$T,$C,$A,$G,$alpha);
		$trees.=sprintf("\[TREE\] T%.*d %s\n",$n_digits,$locus,$itree);
		$partitions.=sprintf("\[PARTITIONS\] T%.*d \[T%.*d GTR%.*d %s\]\n",$n_digits,$locus,$n_digits,$locus,$n_digits,$locus,$length);
		$evolves.=sprintf("T%.*d 1 %.*d\n",$n_digits,$locus,$n_digits,$locus);
		
		#SQLite write
		$SQL->do("INSERT INTO Sequences VALUES ($sequence_counter,$sp_counter,$locus,1,$length)");
		#$SQL->do("UPDATE Species_Trees SET Alpha_site=$alpha WHERE SID==$sp_counter");
		
		$sequence_counter+=1;
	}
	
	print $filehandwrite $models,$settings,$trees,$partitions,$evolves;
	close($filehandwrite);
	
	print "\tFile created\n";
	
	chdir("..");
	
}

$SQL->commit();
$SQL->disconnect or die $DBI::errstr;

#print "Parallel sequence simulation\n";
#system("ls -d [0-9][0-9][0-9] | parallel -P $n_threads 'cd {} && indelible'");

